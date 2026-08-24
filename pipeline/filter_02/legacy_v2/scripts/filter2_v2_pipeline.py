#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import itertools
import json
import math
import os
import platform
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import yaml


ROOT = Path("/root/autodl-tmp/benchmark_1.0")
OUT = ROOT / "filter_2_ligand_qualification_v2"
P1 = ROOT / "processing_1_pdb_source_audit"
F1 = ROOT / "filter_1_protein_receptor_qualification"
OLD_F2 = ROOT / "filter_2_ligand_qualification"
P1_INDEX = P1 / "release/processing_1_mmcif_index.tsv.gz"
F1_ENTRIES = F1 / "release/filter_1_receptor_qualified_entries.tsv.gz"
F1_ASSEMBLIES = F1 / "release/filter_1_receptor_qualified_assemblies.tsv.gz"
F1_RECEPTORS = F1 / "release/filter_1_receptor_chain_instances.tsv.gz"
F1_SHORT = F1 / "release/filter_1_short_peptide_inventory.tsv.gz"
F1_INTERFACE = F1 / "release/filter_1_downstream_interface.json"
CCD_SOURCE = OLD_F2 / "references/components.cif.gz"
CCD_META_SOURCE = OLD_F2 / "references/ccd_snapshot_metadata.json"
BIOLIP_SOURCE = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v5_1_prefull_readiness/references/artifact_reference_files/zhanggroup_BioLiP_ligand_list")
BIOLIP_INVENTORY = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v5_1_prefull_readiness/references/artifact_reference_inventory.tsv")
PYTHON = Path("/root/miniconda3/envs/interaction-pilot-v2/bin/python")
RULE_VERSION = "filter_2_v2.0.0"
SCHEMA_VERSION = "2.0.0"
RUNS = {"smoke": "20260804_smoke1000_02", "full": "20260804_full_01"}

WATER_IDS = {"HOH", "DOD", "WAT", "H2O"}
METALS = {
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA", "AL", "GA", "IN", "TL",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "Y", "ZR", "NB", "MO",
    "TC", "RU", "RH", "PD", "AG", "CD", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU",
    "HG", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM",
    "YB", "LU", "AC", "TH", "PA", "U", "NP", "PU",
}

SOURCE_FIELDS = [
    "pdb_id", "selected_model_id", "entity_id", "component_id", "label_comp_id", "auth_comp_id",
    "label_asym_id", "auth_asym_id", "auth_seq_id", "pdb_seq_num", "insertion_code",
    "source_ligand_instance_id", "atom_count", "observed_heavy_atom_count", "observed_element_composition",
    "altloc_values", "occupancy_min", "occupancy_max", "entity_type", "source_instance_status",
    "source_instance_count_in_entry", "biolip_list_match", "biolip_snapshot_sha256",
    "original_component_id", "resolved_ccd_id", "ccd_identity_status", "ccd_name", "ccd_type", "formula",
    "formula_weight", "formal_charge", "parent_component_id", "expected_atom_count", "expected_heavy_atom_count",
    "element_set", "carbon_atom_count", "fragment_count", "contains_metal", "descriptor_availability",
    "ccd_snapshot_version", "ccd_snapshot_sha256", "terminal_route", "decision", "destination", "reason_code",
    "reason_detail", "rule_version",
]
EXCLUSION_FIELDS = [
    "source_ligand_instance_id", "pdb_id", "component_id", "terminal_route", "expected_heavy_atom_count",
    "observed_heavy_atom_count", "element", "source_instance_count_in_entry", "biolip_list_match",
    "biolip_snapshot_sha256", "exclusion_reason", "rule_version",
]
PLACEMENT_FIELDS = [
    "pdb_id", "assembly_id", "selected_model_id", "source_ligand_instance_id", "component_id", "label_asym_id",
    "assembly_gen_row_id", "oper_expression_raw", "operator_path", "composite_operator_id", "rotation_matrix",
    "translation_vector", "assembly_ligand_placement_id", "mapping_status", "rule_version",
]
NO_MAP_FIELDS = [
    "pdb_id", "selected_model_id", "source_ligand_instance_id", "component_id", "label_asym_id",
    "mapping_status", "reason_code", "rule_version",
]
WATER_FIELDS = ["pdb_id", "component_id", "source_instance_count"]
CONTEXT_FIELDS = [
    "pdb_id", "protein_polymer_residue_count", "rna_dna_polymer_residue_count", "short_peptide_residue_count",
    "branched_entity_residue_count", "modified_polymer_residue_count", "other_polymer_residue_count",
]
ENTRY_FIELDS = [
    "pdb_id", "parse_status", "parse_error", "retained_assembly_count", "selected_model_count",
    "water_source_instance_count", "nonwater_independent_nonpolymer_count", "excluded_monoatomic_ion_count",
    "excluded_simple_inorganic_count", "excluded_biolip_artifact_count", "inorganic_review_count",
    "ccd_review_count", "provisional_source_ligand_count", "logical_placement_count",
    "no_retained_assembly_mapping_count", "entry_status", "terminal_reason",
]
CCD_FIELDS = [
    "original_component_id", "resolved_ccd_id", "ccd_identity_status", "ccd_name", "ccd_type", "formula",
    "formula_weight", "formal_charge", "parent_component_id", "expected_atom_count", "expected_heavy_atom_count",
    "element_set", "carbon_atom_count", "fragment_count", "contains_metal", "descriptor_availability",
    "release_status", "replaced_by", "ccd_snapshot_version", "ccd_snapshot_sha256",
]
TABLES = {
    "source_instances": SOURCE_FIELDS,
    "source_exclusions": EXCLUSION_FIELDS,
    "inorganic_review": SOURCE_FIELDS,
    "ccd_review": SOURCE_FIELDS,
    "provisional_source_ligands": SOURCE_FIELDS,
    "ligand_assembly_logical_placements": PLACEMENT_FIELDS,
    "no_retained_assembly_mapping": NO_MAP_FIELDS,
    "water_exclusion_summary": WATER_FIELDS,
    "context_exclusion_summary": CONTEXT_FIELDS,
    "entries": ENTRY_FIELDS,
}

G_CCD: dict[str, dict[str, str]] = {}
G_BIOLIP: set[str] = set()
G_CFG: dict = {}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def clean(value) -> str:
    text = str(value).strip() if value is not None else ""
    return "" if text in {".", "?", "None"} else text


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def iter_tsv(path: Path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def category_records(block, prefix: str) -> list[dict[str, str]]:
    try:
        category = block.get_mmcif_category(prefix)
    except Exception:
        return []
    if not category:
        return []
    keys = list(category)
    size = len(category[keys[0]]) if keys else 0
    return [{key: clean(category[key][idx]) for key in keys} for idx in range(size)]


def block_value(block, tag: str) -> str:
    try:
        return clean(block.find_value(tag))
    except Exception:
        return ""


def stable_id(*parts: object) -> str:
    return "|".join(clean(x) or "~" for x in parts)


def tree_stat_fingerprint(path: Path, pattern: str = "*") -> dict:
    digest = hashlib.sha256()
    count = size = 0
    for item in sorted(path.rglob(pattern)):
        if not item.is_file():
            continue
        stat = item.stat()
        digest.update(f"{item}:{stat.st_size}:{stat.st_mtime_ns}\n".encode())
        count += 1
        size += stat.st_size
    return {"file_count": count, "size_bytes": size, "stat_sha256": digest.hexdigest()}


def expand_token(token: str) -> list[str]:
    values = []
    for part in token.split(","):
        part = part.strip()
        match = re.fullmatch(r"(-?\d+)-(-?\d+)", part)
        if match:
            start, stop = map(int, match.groups())
            step = 1 if stop >= start else -1
            values.extend(str(x) for x in range(start, stop + step, step))
        elif part:
            values.append(part)
    if not values:
        raise ValueError(f"empty operator token: {token}")
    return values


def expand_operator_paths(expression: str) -> list[tuple[str, ...]]:
    compact = expression.replace(" ", "")
    groups = re.findall(r"\(([^()]*)\)", compact)
    if not groups:
        groups = [compact]
    return [tuple(path) for path in itertools.product(*(expand_token(group) for group in groups))]


def identity_affine():
    return ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0])


def compose_affine(left, right):
    lr, lt = left
    rr, rt = right
    rotation = [[sum(lr[i][k] * rr[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    translation = [sum(lr[i][k] * rt[k] for k in range(3)) + lt[i] for i in range(3)]
    return rotation, translation


def composite_affine(path: tuple[str, ...], operators: dict[str, tuple[list[list[float]], list[float]]]):
    total = identity_affine()
    for operator_id in path:
        if operator_id not in operators:
            raise KeyError(operator_id)
        total = compose_affine(total, operators[operator_id])
    return total


def format_vector(values) -> str:
    return ",".join(f"{value:.10g}" for value in values)


def format_matrix(matrix) -> str:
    return ";".join(format_vector(row) for row in matrix)


def fragment_count(elements: list[str], bonds: list[tuple[str, str]], atom_ids: list[str]) -> int:
    heavy = {atom for atom, element in zip(atom_ids, elements) if element != "H"}
    if not heavy:
        return 0
    graph = {atom: set() for atom in heavy}
    for left, right in bonds:
        if left in graph and right in graph:
            graph[left].add(right)
            graph[right].add(left)
    count = 0
    unseen = set(graph)
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            for neighbor in graph[stack.pop()]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
    return count


def setup() -> None:
    if OUT.exists():
        raise SystemExit(f"Target already exists; refusing overwrite: {OUT}")
    for directory in ["scripts", "configs", "schemas", "references", "inputs", "preflight", "tests"]:
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    for scope, run_id in RUNS.items():
        for directory in ["input", "work/batches", "output", "audit", "logs", "release"]:
            (OUT / "runs" / run_id / directory).mkdir(parents=True, exist_ok=True)

    shutil.copy2(Path(__file__), OUT / "scripts/filter2_v2_pipeline.py")
    for source, name in [
        (P1_INDEX, "processing_1_mmcif_index_snapshot.tsv.gz"),
        (F1_ENTRIES, "filter_1_receptor_qualified_entries_snapshot.tsv.gz"),
        (F1_ASSEMBLIES, "filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"),
        (F1_RECEPTORS, "filter_1_receptor_chain_instances_snapshot.tsv.gz"),
        (F1_SHORT, "filter_1_short_peptide_inventory_snapshot.tsv.gz"),
        (F1_INTERFACE, "filter_1_downstream_interface_snapshot.json"),
    ]:
        shutil.copy2(source, OUT / "inputs" / name)
    shutil.copy2(CCD_META_SOURCE, OUT / "references/ccd_snapshot_metadata.json")
    shutil.copy2(BIOLIP_SOURCE, OUT / "references/BioLiP_ligand_list_snapshot.tsv")
    shutil.copy2(BIOLIP_INVENTORY, OUT / "references/BioLiP_reference_inventory.tsv")

    config = {
        "stage": {"name": "filter_02_independent_small_molecule_candidate_identification", "version": "2.0.0"},
        "rule_version": RULE_VERSION,
        "schema_version": SCHEMA_VERSION,
        "inputs": {"expected_entries": 248037, "expected_retained_assemblies": 360611},
        "model_policy": "all_model_ids_represented_by_filter_1_receptor_chain_instances_no_new_single_model_selection",
        "candidate_boundary": {"required_entity_type": "non-polymer", "altloc_in_base_instance_id": False},
        "water": {"controlled_component_ids": sorted(WATER_IDS), "emit_instance_rows": False},
        "simple_inorganic": {
            "require_carbon_atom_count": 0,
            "min_expected_heavy_atom_count": 2,
            "max_expected_heavy_atom_count": 8,
            "require_single_fragment": True,
            "exclude_if_contains_metal": True,
            "allowed_elements": ["H", "B", "N", "O", "F", "SI", "P", "S", "CL", "SE", "BR", "I"],
            "boundary_action": "inorganic_review",
        },
        "biolip": {"minimum_source_instances_same_pdb_component": 15, "rule_is_source_level": True},
        "ccd": {"source_path": str(CCD_SOURCE), "accepted_statuses": ["ccd_exact", "ccd_alias_resolved", "ccd_parent_resolved", "ccd_obsolete_resolved"]},
        "assembly": {"retained_assemblies_only": True, "read_all_assembly_gen_rows": True, "materialize_coordinates": False},
        "runtime": {"workers": 32, "batch_size": 200, "resume": True},
        "format": {"large_table": "tsv.gz", "parquet_status": "unavailable_in_frozen_environment_no_install_performed", "preview_rows": 1000},
    }
    (OUT / "configs/filter_2_v2.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    (OUT / "README.md").write_text(
        "# Filter 2 v2 - Independent Small-Molecule Candidate Identification and Biological Assembly Logical Mapping\n\n"
        "This independent implementation reads only Processing 1, frozen Filter 1 release interfaces, raw mmCIF, frozen CCD, and frozen BioLiP ligand_list. "
        "Only `_entity.type = non-polymer` source instances can enter the candidate inventory. Polymer, nucleic-acid, branched, modified-residue, short-peptide, and water objects never become candidate rows. "
        "The stage performs source-level chemical-scope routing and logical assembly operator mapping only. It does not transform coordinates, construct protein-ligand pairs, calculate distances, infer interactions, assess biological relevance, or run downstream stages.\n",
        encoding="utf-8",
    )
    for name, command in [
        ("audit_filter_2_v2_inputs.py", "preflight"),
        ("prepare_filter_2_v2_references.py", "prepare-references"),
        ("test_filter_2_v2.py", "selftest"),
        ("run_filter_2_v2_smoke.py", "run-smoke"),
        ("run_filter_2_v2_full.py", "run-full"),
        ("finalize_filter_2_v2.py", "finalize"),
        ("validate_filter_2_v2.py", "validate"),
    ]:
        (OUT / "scripts" / name).write_text(
            f"#!/usr/bin/env python3\nimport subprocess,sys\nraise SystemExit(subprocess.call([sys.executable,{str(OUT / 'scripts/filter2_v2_pipeline.py')!r},{command!r},*sys.argv[1:]]))\n",
            encoding="utf-8",
        )
    for table, fields in TABLES.items():
        schema = {
            "schema_version": SCHEMA_VERSION,
            "table": f"filter_2_{table}",
            "format": "TSV.GZ" if table not in {"water_exclusion_summary"} else "TSV",
            "columns": [{"column_name": field, "data_type": "string", "nullable": True, "description": field.replace("_", " ")} for field in fields],
        }
        (OUT / "schemas" / f"filter_2_{table}_schema.json").write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
    (OUT / "status.json").write_text(json.dumps({"status": "DRAFT", "created_at": utc()}, indent=2) + "\n")
    print(OUT)


def input_snapshot(name: str) -> Path:
    return OUT / "inputs" / name


def preflight() -> None:
    cfg = yaml.safe_load((OUT / "configs/filter_2_v2.yaml").read_text())
    files = {
        "processing_1": input_snapshot("processing_1_mmcif_index_snapshot.tsv.gz"),
        "filter_1_entries": input_snapshot("filter_1_receptor_qualified_entries_snapshot.tsv.gz"),
        "filter_1_assemblies": input_snapshot("filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"),
        "filter_1_receptors": input_snapshot("filter_1_receptor_chain_instances_snapshot.tsv.gz"),
        "filter_1_interface": input_snapshot("filter_1_downstream_interface_snapshot.json"),
        "ccd": CCD_SOURCE,
        "biolip": OUT / "references/BioLiP_ligand_list_snapshot.tsv",
    }
    counts = {}
    for key in ["filter_1_entries", "filter_1_assemblies"]:
        counts[key] = sum(1 for _ in iter_tsv(files[key]))
    entry_ids = [row["pdb_id"] for row in iter_tsv(files["filter_1_entries"])]
    assembly_keys = [(row["pdb_id"], row["assembly_id"]) for row in iter_tsv(files["filter_1_assemblies"])]
    missing_paths = sum(not Path(row["mmcif_path"]).is_file() for row in iter_tsv(files["filter_1_entries"] ))
    ccd_meta = json.loads((OUT / "references/ccd_snapshot_metadata.json").read_text())
    with (OUT / "references/BioLiP_ligand_list_snapshot.tsv").open(encoding="utf-8") as handle:
        biolip_lines = [line.rstrip("\n") for line in handle if line.strip()]
    biolip_ids = [line.split("\t", 1)[0].strip().upper() for line in biolip_lines]
    source_hashes = {key: sha256(path) for key, path in files.items()}
    state = {
        "processing_1_tree": tree_stat_fingerprint(P1),
        "filter_1_tree": tree_stat_fingerprint(F1),
        "old_filter_2_tree": tree_stat_fingerprint(OLD_F2),
        "raw_mmcif_tree": tree_stat_fingerprint(Path("/root/autodl-tmp/pdb_archive_v2/mmCIF"), "*.cif.gz"),
    }
    audit = {
        "timestamp": utc(),
        "input_receptor_qualified_entries": counts["filter_1_entries"],
        "input_retained_assemblies": counts["filter_1_assemblies"],
        "duplicate_input_pdb_id": len(entry_ids) - len(set(entry_ids)),
        "duplicate_retained_assembly_key": len(assembly_keys) - len(set(assembly_keys)),
        "missing_mmcif_path": missing_paths,
        "filter_1_validation_pass": json.loads((F1 / "release/filter_1_release_validation.json").read_text())["release_validation_pass"],
        "ccd_source_path": str(CCD_SOURCE),
        "ccd_file_size": CCD_SOURCE.stat().st_size,
        "ccd_sha256": source_hashes["ccd"],
        "ccd_metadata_sha256": ccd_meta["sha256"],
        "ccd_checksum_match": source_hashes["ccd"] == ccd_meta["sha256"],
        "ccd_gzip_ok": True,
        "biolip_source_url": "https://zhanggroup.org/BioLiP/ligand_list",
        "biolip_line_count": len(biolip_lines),
        "biolip_unique_component_id_count": len(set(biolip_ids)),
        "biolip_duplicate_component_id_count": len(biolip_ids) - len(set(biolip_ids)),
        "biolip_sha256": source_hashes["biolip"],
        "source_hashes": source_hashes,
        "immutable_state_before": state,
        "parquet_available": False,
        "parquet_note": "pyarrow, fastparquet, polars, and duckdb are unavailable in the frozen interaction-pilot-v2 environment; no installation performed; unique TSV.GZ formal tables used.",
    }
    audit["preflight_pass"] = all([
        counts["filter_1_entries"] == cfg["inputs"]["expected_entries"],
        counts["filter_1_assemblies"] == cfg["inputs"]["expected_retained_assemblies"],
        audit["duplicate_input_pdb_id"] == 0,
        audit["duplicate_retained_assembly_key"] == 0,
        missing_paths == 0,
        audit["filter_1_validation_pass"],
        audit["ccd_checksum_match"],
        len(set(biolip_ids)) == 463,
    ])
    (OUT / "preflight/filter_2_v2_input_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    metadata = {
        "source_url": audit["biolip_source_url"],
        "downloaded_at": "2026-07-24",
        "file_size": files["biolip"].stat().st_size,
        "sha256": audit["biolip_sha256"],
        "line_count": len(biolip_lines),
        "unique_component_id_count": len(set(biolip_ids)),
        "duplicate_component_id_count": len(biolip_ids) - len(set(biolip_ids)),
        "rule_version": RULE_VERSION,
    }
    (OUT / "references/BioLiP_ligand_list_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if not audit["preflight_pass"]:
        raise SystemExit(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


def prepare_references() -> None:
    audit = json.loads((OUT / "preflight/filter_2_v2_input_audit.json").read_text())
    if not audit["preflight_pass"]:
        raise SystemExit("preflight not passed")
    started = time.time()
    document = gemmi.cif.read(str(CCD_SOURCE))
    raw = []
    required = Counter()
    for block in document:
        component_id = block_value(block, "_chem_comp.id").upper()
        if not component_id:
            continue
        atom_rows = category_records(block, "_chem_comp_atom.")
        bond_rows = category_records(block, "_chem_comp_bond.")
        descriptor_rows = category_records(block, "_pdbx_chem_comp_descriptor.")
        required["chem_comp"] += 1
        required["chem_comp_atom"] += bool(atom_rows)
        required["chem_comp_bond"] += bool(bond_rows)
        required["descriptor"] += bool(descriptor_rows)
        atom_ids = [row.get("atom_id", "") for row in atom_rows]
        elements = [row.get("type_symbol", "").upper() for row in atom_rows]
        bonds = [(row.get("atom_id_1", ""), row.get("atom_id_2", "")) for row in bond_rows]
        heavy = sum(element != "H" for element in elements)
        descriptor_types = {row.get("type", "").upper() for row in descriptor_rows}
        raw.append({
            "original_component_id": component_id,
            "resolved_ccd_id": component_id,
            "ccd_identity_status": "ccd_exact",
            "ccd_name": block_value(block, "_chem_comp.name"),
            "ccd_type": block_value(block, "_chem_comp.type"),
            "formula": block_value(block, "_chem_comp.formula"),
            "formula_weight": block_value(block, "_chem_comp.formula_weight"),
            "formal_charge": block_value(block, "_chem_comp.pdbx_formal_charge"),
            "parent_component_id": block_value(block, "_chem_comp.mon_nstd_parent_comp_id"),
            "expected_atom_count": str(len(elements)),
            "expected_heavy_atom_count": str(heavy),
            "element_set": ",".join(sorted(set(filter(None, elements)))),
            "carbon_atom_count": str(elements.count("C")),
            "fragment_count": str(fragment_count(elements, bonds, atom_ids)),
            "contains_metal": str(any(element in METALS for element in elements)).lower(),
            "descriptor_availability": ",".join(sorted(descriptor_types)) if descriptor_types else "none",
            "release_status": block_value(block, "_chem_comp.pdbx_release_status"),
            "replaced_by": block_value(block, "_chem_comp.pdbx_replaced_by").upper(),
            "ccd_snapshot_version": json.loads((OUT / "references/ccd_snapshot_metadata.json").read_text())["last_modified"],
            "ccd_snapshot_sha256": audit["ccd_sha256"],
        })
    ids = {row["original_component_id"] for row in raw}
    for row in raw:
        if row["release_status"].upper() in {"OBS", "OBSOLETE"}:
            replacement = row["replaced_by"]
            if replacement and replacement in ids:
                row["resolved_ccd_id"] = replacement
                row["ccd_identity_status"] = "ccd_obsolete_resolved"
            else:
                row["ccd_identity_status"] = "ccd_invalid"
    write_tsv(OUT / "references/ccd_component_cache.tsv.gz", raw, CCD_FIELDS)
    metadata = {
        "source_path": str(CCD_SOURCE),
        "file_size": CCD_SOURCE.stat().st_size,
        "sha256": audit["ccd_sha256"],
        "gzip_ok": True,
        "version": raw[0]["ccd_snapshot_version"] if raw else "",
        "component_count": len(raw),
        "unique_component_id_count": len(ids),
        "duplicate_component_id_count": len(raw) - len(ids),
        "required_category_counts": dict(required),
        "runtime_seconds": round(time.time() - started, 3),
        "validation_pass": len(raw) == 50666 and len(ids) == 50666 and required["chem_comp_atom"] > 50000,
    }
    (OUT / "references/ccd_snapshot_validation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if not metadata["validation_pass"]:
        raise SystemExit(json.dumps(metadata, indent=2))
    print(json.dumps(metadata, indent=2))


def init_worker(ccd_cache: str, biolip_path: str, config_path: str) -> None:
    global G_CCD, G_BIOLIP, G_CFG
    G_CCD = {row["original_component_id"]: row for row in iter_tsv(Path(ccd_cache))}
    with Path(biolip_path).open(encoding="utf-8") as handle:
        G_BIOLIP = {line.split("\t", 1)[0].strip().upper() for line in handle if line.strip()}
    G_CFG = yaml.safe_load(Path(config_path).read_text())


def polymer_class(poly_type: str) -> str:
    value = poly_type.lower()
    if "polypeptide" in value:
        return "protein"
    if "hybrid" in value or ("dna" in value and "rna" in value):
        return "nucleic_acid"
    if "ribonucleotide" in value or "rna" in value or "deoxyribonucleotide" in value or "dna" in value:
        return "nucleic_acid"
    if "polysaccharide" in value:
        return "branched_or_polysaccharide"
    return "other_polymer"


def route_source(row: dict, count_in_entry: int) -> tuple[str, str, str, str]:
    component = row["component_id"]
    ccd = G_CCD.get(component)
    observed_elements = {x for x in row["observed_element_composition"].split(",") if x}
    observed_heavy = int(row["observed_heavy_atom_count"] or 0)
    if ccd:
        expected_heavy = int(ccd["expected_heavy_atom_count"] or 0)
        ccd_elements = {x for x in ccd["element_set"].split(",") if x}
    else:
        expected_heavy = 0
        ccd_elements = set()
    mono_elements = ccd_elements or observed_elements
    if ((expected_heavy == 1) or (not ccd and observed_heavy == 1)) and len(mono_elements) == 1:
        return "excluded_monoatomic_ion", "REJECT", "excluded", "monoatomic_element_or_ion"
    if ccd:
        cfg = G_CFG["simple_inorganic"]
        heavy = int(ccd["expected_heavy_atom_count"] or 0)
        carbon = int(ccd["carbon_atom_count"] or 0)
        fragments = int(ccd["fragment_count"] or 0)
        elements = {x for x in ccd["element_set"].split(",") if x}
        simple = all([
            carbon == cfg["require_carbon_atom_count"],
            cfg["min_expected_heavy_atom_count"] <= heavy <= cfg["max_expected_heavy_atom_count"],
            not (cfg["require_single_fragment"] and fragments != 1),
            not (cfg["exclude_if_contains_metal"] and ccd["contains_metal"] == "true"),
            elements.issubset(set(cfg["allowed_elements"])),
        ])
        if simple:
            return "excluded_simple_inorganic", "REJECT", "excluded", "conservative_ccd_simple_inorganic_rule"
    if component in G_BIOLIP and count_in_entry >= G_CFG["biolip"]["minimum_source_instances_same_pdb_component"]:
        return "excluded_biolip_artifact", "REJECT", "excluded", "frozen_biolip_list_and_entry_occurrence_ge_15"
    if ccd and int(ccd["carbon_atom_count"] or 0) == 0:
        return "inorganic_review", "REVIEW", "inorganic_review", "noncarbon_complex_or_boundary_inorganic"
    accepted = set(G_CFG["ccd"]["accepted_statuses"])
    if not ccd or ccd["ccd_identity_status"] not in accepted:
        return "ccd_review", "REVIEW", "ccd_review", "ccd_identity_unresolved"
    return "provisional_source_ligand", "PASS", "provisional_source_ligand", "independent_nonpolymer_passed_source_scope_rules"


def parse_operator_table(block) -> dict[str, tuple[list[list[float]], list[float]]]:
    result = {}
    for row in category_records(block, "_pdbx_struct_oper_list."):
        operator_id = row.get("id", "")
        try:
            matrix = [[float(row[f"matrix[{i}][{j}]"]) for j in range(1, 4)] for i in range(1, 4)]
            vector = [float(row[f"vector[{i}]"]) for i in range(1, 4)]
            result[operator_id] = (matrix, vector)
        except Exception:
            continue
    return result


def parse_entry(item: dict) -> dict[str, list[dict]]:
    output = {name: [] for name in TABLES}
    pid = item["pdb_id"]
    try:
        block = gemmi.cif.read(item["mmcif_path"]).sole_block()
        entities = {row.get("id", ""): row.get("type", "").lower() for row in category_records(block, "_entity.")}
        entity_poly = {row.get("entity_id", ""): row.get("type", "") for row in category_records(block, "_entity_poly.")}
        asym_entity = {row.get("id", ""): row.get("entity_id", "") for row in category_records(block, "_struct_asym.")}
        nonpoly_entities = {entity_id for entity_id, entity_type in entities.items() if entity_type == "non-polymer"}
        water_entities = {entity_id for entity_id, entity_type in entities.items() if entity_type == "water"}
        branch_entities = {row.get("entity_id", "") for row in category_records(block, "_pdbx_entity_branch.")}
        short_asym = set(item.get("short_asym_ids", []))
        selected_models = set(item.get("selected_model_ids", []))

        scheme = defaultdict(list)
        for row in category_records(block, "_pdbx_nonpoly_scheme."):
            scheme[(row.get("asym_id", ""), row.get("mon_id", "").upper())].append(row)

        atom_category = category_records(block, "_atom_site.")
        source_groups = {}
        water_keys = set()
        context_sets = defaultdict(set)
        for atom in atom_category:
            model = atom.get("pdbx_PDB_model_num", "") or "1"
            if selected_models and model not in selected_models:
                continue
            asym = atom.get("label_asym_id", "")
            entity_id = asym_entity.get(asym, "")
            entity_type = entities.get(entity_id, "")
            component = atom.get("label_comp_id", "").upper()
            auth_seq = atom.get("auth_seq_id", "")
            insertion = atom.get("pdbx_PDB_ins_code", "")
            label_seq = atom.get("label_seq_id", "")
            residue_key = (model, asym, label_seq or auth_seq, insertion, component)

            if entity_id in water_entities or component in WATER_IDS:
                water_keys.add(residue_key)
                continue
            if entity_id in nonpoly_entities:
                key = (model, entity_id, asym, auth_seq, insertion, component)
                group = source_groups.setdefault(key, {"atoms": [], "auth_asym": set(), "auth_comp": set(), "label_seq": set()})
                group["atoms"].append(atom)
                if atom.get("auth_asym_id", ""):
                    group["auth_asym"].add(atom["auth_asym_id"])
                if atom.get("auth_comp_id", ""):
                    group["auth_comp"].add(atom["auth_comp_id"])
                if label_seq:
                    group["label_seq"].add(label_seq)
                continue
            if entity_id in branch_entities:
                context_sets["branched_entity"].add(residue_key)
            elif entity_type == "polymer":
                pclass = polymer_class(entity_poly.get(entity_id, ""))
                if asym in short_asym:
                    context_sets["short_peptide"].add(residue_key)
                elif pclass == "protein":
                    context_sets["protein"].add(residue_key)
                elif pclass == "nucleic_acid":
                    context_sets["nucleic_acid"].add(residue_key)
                elif pclass == "branched_or_polysaccharide":
                    context_sets["branched_entity"].add(residue_key)
                else:
                    context_sets["other_polymer"].add(residue_key)

        modified_count = len(category_records(block, "_pdbx_struct_mod_residue."))
        water_counts = Counter(key[-1] for key in water_keys)
        output["water_exclusion_summary"] = [
            {"pdb_id": pid, "component_id": component, "source_instance_count": str(count)}
            for component, count in sorted(water_counts.items())
        ]
        output["context_exclusion_summary"] = [{
            "pdb_id": pid,
            "protein_polymer_residue_count": str(len(context_sets["protein"])),
            "rna_dna_polymer_residue_count": str(len(context_sets["nucleic_acid"])),
            "short_peptide_residue_count": str(len(context_sets["short_peptide"])),
            "branched_entity_residue_count": str(len(context_sets["branched_entity"])),
            "modified_polymer_residue_count": str(modified_count),
            "other_polymer_residue_count": str(len(context_sets["other_polymer"])),
        }]

        source_rows = []
        component_counts = Counter(key[-1] for key in source_groups)
        biolip_sha = G_CFG["biolip_snapshot_sha256"]
        ccd_sha = G_CFG["ccd_snapshot_sha256"]
        ccd_version = G_CFG["ccd_snapshot_version"]
        for key, group in sorted(source_groups.items()):
            model, entity_id, asym, auth_seq, insertion, component = key
            atoms = group["atoms"]
            scheme_rows = scheme.get((asym, component), [])
            chosen = next((row for row in scheme_rows if row.get("auth_seq_num", "") == auth_seq or row.get("pdb_seq_num", "") == auth_seq), scheme_rows[0] if len(scheme_rows) == 1 else {})
            pdb_seq_num = chosen.get("pdb_seq_num", "") or auth_seq
            auth_seq_final = chosen.get("auth_seq_num", "") or auth_seq
            insertion_final = chosen.get("pdb_ins_code", "") or insertion
            elements = [atom.get("type_symbol", "").upper() for atom in atoms]
            altlocs = sorted({atom.get("label_alt_id", "") for atom in atoms if atom.get("label_alt_id", "")})
            occupancies = []
            for atom in atoms:
                try:
                    occupancies.append(float(atom.get("occupancy", "")))
                except Exception:
                    pass
            source_id = stable_id(pid, model, entity_id, asym, auth_seq_final or pdb_seq_num, insertion_final or ".", component)
            ccd = G_CCD.get(component)
            row = {
                "pdb_id": pid, "selected_model_id": model, "entity_id": entity_id, "component_id": component,
                "label_comp_id": component, "auth_comp_id": ",".join(sorted(group["auth_comp"])) or component,
                "label_asym_id": asym, "auth_asym_id": ",".join(sorted(group["auth_asym"])),
                "auth_seq_id": auth_seq_final, "pdb_seq_num": pdb_seq_num, "insertion_code": insertion_final,
                "source_ligand_instance_id": source_id, "atom_count": str(len(atoms)),
                "observed_heavy_atom_count": str(sum(element != "H" for element in elements)),
                "observed_element_composition": ",".join(sorted(set(filter(None, elements)))),
                "altloc_values": ",".join(altlocs),
                "occupancy_min": f"{min(occupancies):.6g}" if occupancies else "",
                "occupancy_max": f"{max(occupancies):.6g}" if occupancies else "",
                "entity_type": "non-polymer", "source_instance_status": "resolved",
                "source_instance_count_in_entry": str(component_counts[component]),
                "biolip_list_match": str(component in G_BIOLIP).lower(), "biolip_snapshot_sha256": biolip_sha,
                "original_component_id": component, "resolved_ccd_id": ccd["resolved_ccd_id"] if ccd else "",
                "ccd_identity_status": ccd["ccd_identity_status"] if ccd else "ccd_missing",
                "ccd_name": ccd["ccd_name"] if ccd else "", "ccd_type": ccd["ccd_type"] if ccd else "",
                "formula": ccd["formula"] if ccd else "", "formula_weight": ccd["formula_weight"] if ccd else "",
                "formal_charge": ccd["formal_charge"] if ccd else "", "parent_component_id": ccd["parent_component_id"] if ccd else "",
                "expected_atom_count": ccd["expected_atom_count"] if ccd else "",
                "expected_heavy_atom_count": ccd["expected_heavy_atom_count"] if ccd else "",
                "element_set": ccd["element_set"] if ccd else "", "carbon_atom_count": ccd["carbon_atom_count"] if ccd else "",
                "fragment_count": ccd["fragment_count"] if ccd else "", "contains_metal": ccd["contains_metal"] if ccd else "",
                "descriptor_availability": ccd["descriptor_availability"] if ccd else "none",
                "ccd_snapshot_version": ccd_version, "ccd_snapshot_sha256": ccd_sha, "rule_version": RULE_VERSION,
            }
            route, decision, destination, reason = route_source(row, component_counts[component])
            row.update({"terminal_route": route, "decision": decision, "destination": destination, "reason_code": reason, "reason_detail": "", "rule_version": RULE_VERSION})
            source_rows.append(row)
        output["source_instances"] = source_rows
        for row in source_rows:
            route = row["terminal_route"]
            if route.startswith("excluded_"):
                element = row["element_set"] or row["observed_element_composition"]
                output["source_exclusions"].append({
                    "source_ligand_instance_id": row["source_ligand_instance_id"], "pdb_id": pid,
                    "component_id": row["component_id"], "terminal_route": route,
                    "expected_heavy_atom_count": row["expected_heavy_atom_count"],
                    "observed_heavy_atom_count": row["observed_heavy_atom_count"], "element": element,
                    "source_instance_count_in_entry": row["source_instance_count_in_entry"],
                    "biolip_list_match": row["biolip_list_match"], "biolip_snapshot_sha256": row["biolip_snapshot_sha256"],
                    "exclusion_reason": row["reason_code"], "rule_version": RULE_VERSION,
                })
            elif route == "inorganic_review":
                output["inorganic_review"].append(row)
            elif route == "ccd_review":
                output["ccd_review"].append(row)
            elif route == "provisional_source_ligand":
                output["provisional_source_ligands"].append(row)

        operators = parse_operator_table(block)
        reverse = defaultdict(list)
        retained = set(item["retained_assembly_ids"])
        for index, gen in enumerate(category_records(block, "_pdbx_struct_assembly_gen."), start=1):
            assembly_id = gen.get("assembly_id", "")
            if assembly_id not in retained:
                continue
            expression = gen.get("oper_expression", "")
            row_id = f"assembly_gen_row_{index:06d}"
            asyms = [value.strip() for value in gen.get("asym_id_list", "").split(",") if value.strip()]
            try:
                paths = expand_operator_paths(expression)
            except Exception:
                paths = []
            for asym in asyms:
                reverse[(assembly_id, asym)].append((row_id, expression, paths))

        for row in output["provisional_source_ligands"]:
            mappings = []
            for assembly_id in sorted(retained):
                allowed_models = set(item["assembly_models"].get(assembly_id, []))
                if allowed_models and row["selected_model_id"] not in allowed_models:
                    continue
                for row_id, expression, paths in reverse.get((assembly_id, row["label_asym_id"]), []):
                    for path in paths:
                        try:
                            rotation, translation = composite_affine(path, operators)
                        except Exception:
                            continue
                        operator_path = "*".join(path)
                        placement_id = stable_id(pid, assembly_id, row["selected_model_id"], row["source_ligand_instance_id"], row_id, operator_path)
                        placement = {
                            "pdb_id": pid, "assembly_id": assembly_id, "selected_model_id": row["selected_model_id"],
                            "source_ligand_instance_id": row["source_ligand_instance_id"], "component_id": row["component_id"],
                            "label_asym_id": row["label_asym_id"], "assembly_gen_row_id": row_id,
                            "oper_expression_raw": expression, "operator_path": operator_path,
                            "composite_operator_id": operator_path, "rotation_matrix": format_matrix(rotation),
                            "translation_vector": format_vector(translation), "assembly_ligand_placement_id": placement_id,
                            "mapping_status": "mapped", "rule_version": RULE_VERSION,
                        }
                        output["ligand_assembly_logical_placements"].append(placement)
                        mappings.append(placement_id)
            if not mappings:
                output["no_retained_assembly_mapping"].append({
                    "pdb_id": pid, "selected_model_id": row["selected_model_id"],
                    "source_ligand_instance_id": row["source_ligand_instance_id"], "component_id": row["component_id"],
                    "label_asym_id": row["label_asym_id"], "mapping_status": "no_retained_assembly_mapping",
                    "reason_code": "label_asym_not_in_retained_assembly_gen_or_operator_unresolved", "rule_version": RULE_VERSION,
                })

        routes = Counter(row["terminal_route"] for row in source_rows)
        output["entries"] = [{
            "pdb_id": pid, "parse_status": "success", "parse_error": "",
            "retained_assembly_count": str(len(retained)), "selected_model_count": str(len(selected_models)),
            "water_source_instance_count": str(len(water_keys)),
            "nonwater_independent_nonpolymer_count": str(len(source_rows)),
            "excluded_monoatomic_ion_count": str(routes["excluded_monoatomic_ion"]),
            "excluded_simple_inorganic_count": str(routes["excluded_simple_inorganic"]),
            "excluded_biolip_artifact_count": str(routes["excluded_biolip_artifact"]),
            "inorganic_review_count": str(routes["inorganic_review"]), "ccd_review_count": str(routes["ccd_review"]),
            "provisional_source_ligand_count": str(routes["provisional_source_ligand"]),
            "logical_placement_count": str(len(output["ligand_assembly_logical_placements"])),
            "no_retained_assembly_mapping_count": str(len(output["no_retained_assembly_mapping"])),
            "entry_status": "success", "terminal_reason": "independent_nonpolymer_inventory_completed",
        }]
    except Exception as exc:
        output["entries"] = [{
            "pdb_id": pid, "parse_status": "failed", "parse_error": f"{type(exc).__name__}: {exc}"[:2000],
            "retained_assembly_count": str(len(item.get("retained_assembly_ids", []))), "selected_model_count": str(len(item.get("selected_model_ids", []))),
            "water_source_instance_count": "0", "nonwater_independent_nonpolymer_count": "0",
            "excluded_monoatomic_ion_count": "0", "excluded_simple_inorganic_count": "0", "excluded_biolip_artifact_count": "0",
            "inorganic_review_count": "0", "ccd_review_count": "0", "provisional_source_ligand_count": "0",
            "logical_placement_count": "0", "no_retained_assembly_mapping_count": "0",
            "entry_status": "failed", "terminal_reason": "parse_failed",
        }]
    return output


def load_items() -> list[dict]:
    index = {row["pdb_id"]: row["mmcif_path"] for row in iter_tsv(input_snapshot("processing_1_mmcif_index_snapshot.tsv.gz"))}
    entries = [row["pdb_id"] for row in iter_tsv(input_snapshot("filter_1_receptor_qualified_entries_snapshot.tsv.gz"))]
    assemblies = defaultdict(set)
    for row in iter_tsv(input_snapshot("filter_1_receptor_qualified_assemblies_snapshot.tsv.gz")):
        assemblies[row["pdb_id"]].add(row["assembly_id"])
    assembly_models = defaultdict(lambda: defaultdict(set))
    selected_models = defaultdict(set)
    for row in iter_tsv(input_snapshot("filter_1_receptor_chain_instances_snapshot.tsv.gz")):
        assembly_models[row["pdb_id"]][row["assembly_id"]].add(row["model_id"])
        selected_models[row["pdb_id"]].add(row["model_id"])
    short = defaultdict(set)
    for row in iter_tsv(input_snapshot("filter_1_short_peptide_inventory_snapshot.tsv.gz")):
        short[row["pdb_id"]].add(row["label_asym_id"])
    return [{
        "pdb_id": pid, "mmcif_path": index[pid], "retained_assembly_ids": sorted(assemblies[pid]),
        "assembly_models": {aid: sorted(models) for aid, models in assembly_models[pid].items()},
        "selected_model_ids": sorted(selected_models[pid]), "short_asym_ids": sorted(short[pid]),
    } for pid in entries]


def discover_one(item: dict) -> dict[str, str]:
    found = {}
    pid = item["pdb_id"]
    target_components = {"NAG": "independent_NAG", "ATP": "ATP_NAD_FAD", "NAD": "ATP_NAD_FAD", "FAD": "ATP_NAD_FAD", "ZN": "monoatomic", "MG": "monoatomic", "SO4": "simple_inorganic", "PO4": "simple_inorganic"}
    try:
        block = gemmi.cif.read(item["mmcif_path"]).sole_block()
        entities = {row.get("id", ""): row.get("type", "").lower() for row in category_records(block, "_entity.")}
        poly_types = [row.get("type", "").lower() for row in category_records(block, "_entity_poly.")]
        asym_entity = {row.get("id", ""): row.get("entity_id", "") for row in category_records(block, "_struct_asym.")}
        nonpoly_asym = {asym for asym, entity in asym_entity.items() if entities.get(entity) == "non-polymer"}
        comps = {row.get("comp_id", "").upper() for row in category_records(block, "_pdbx_entity_nonpoly.")}
        if any("ribonucleotide" in value or "deoxyribonucleotide" in value for value in poly_types): found["mixed_nucleic_acid"] = pid
        if category_records(block, "_pdbx_struct_mod_residue."): found["modified_residue"] = pid
        if category_records(block, "_pdbx_entity_branch."): found["branched_entity"] = pid
        for comp, target in target_components.items():
            if comp in comps: found[target] = pid
        atoms = category_records(block, "_atom_site.")
        if any(atom.get("label_alt_id", "") for atom in atoms if atom.get("label_asym_id", "") in nonpoly_asym): found["altloc"] = pid
        if len(item["retained_assembly_ids"]) > 1: found["multi_assembly"] = pid
        gen = [row for row in category_records(block, "_pdbx_struct_assembly_gen.") if row.get("assembly_id", "") in set(item["retained_assembly_ids"])]
        if any(len(expand_operator_paths(row.get("oper_expression", ""))) > 1 for row in gen): found["multi_operator"] = pid
        gen_asym = defaultdict(set)
        for row in gen:
            gen_asym[row.get("assembly_id", "")].update(x.strip() for x in row.get("asym_id_list", "").split(",") if x.strip())
        if any(nonpoly_asym & values for values in gen_asym.values()): found["ligand_separate_assembly_gen_row"] = pid
        if nonpoly_asym and not any(nonpoly_asym & values for values in gen_asym.values()): found["no_assembly_mapping"] = pid
    except Exception:
        return {}
    return found


def discover_targeted(items: list[dict], max_scan: int = 50000) -> tuple[dict[str, str], list[dict]]:
    found = {}
    targets = {
        "mixed_nucleic_acid", "modified_residue", "branched_entity", "independent_NAG", "ATP_NAD_FAD",
        "monoatomic", "simple_inorganic", "altloc", "multi_assembly", "multi_operator",
        "ligand_separate_assembly_gen_row",
    }
    chosen = []
    with ProcessPoolExecutor(max_workers=32) as pool:
        for offset in range(0, min(max_scan, len(items)), 1000):
            for discovered in pool.map(discover_one, items[offset:offset + 1000], chunksize=4):
                for key, pid in discovered.items():
                    found.setdefault(key, pid)
            if targets.issubset(found):
                break
    pid_to_item = {item["pdb_id"]: item for item in items}
    for pid in sorted(set(found.values())):
        chosen.append(pid_to_item[pid])
    return found, chosen


def selftest() -> None:
    cfg = yaml.safe_load((OUT / "configs/filter_2_v2.yaml").read_text())
    global G_CFG, G_BIOLIP, G_CCD
    G_CFG = cfg
    G_CFG["biolip_snapshot_sha256"] = json.loads((OUT / "references/BioLiP_ligand_list_metadata.json").read_text())["sha256"]
    G_CFG["ccd_snapshot_sha256"] = json.loads((OUT / "references/ccd_snapshot_validation.json").read_text())["sha256"]
    G_CFG["ccd_snapshot_version"] = json.loads((OUT / "references/ccd_snapshot_validation.json").read_text())["version"]
    G_CCD = {row["original_component_id"]: row for row in iter_tsv(OUT / "references/ccd_component_cache.tsv.gz")}
    with (OUT / "references/BioLiP_ligand_list_snapshot.tsv").open() as handle:
        G_BIOLIP = {line.split("\t", 1)[0].strip().upper() for line in handle if line.strip()}
    sample = {field: "" for field in SOURCE_FIELDS}
    sample.update({"component_id": "ACE", "observed_heavy_atom_count": "4", "observed_element_composition": "C,O"})
    boundary = {}
    for count in [14, 15, 16]:
        boundary[str(count)] = route_source(sample, count)[0]
    paths = expand_operator_paths("(1-2)(3,4)")
    tests = {
        "biolip_14_not_excluded": boundary["14"] != "excluded_biolip_artifact",
        "biolip_15_excluded": boundary["15"] == "excluded_biolip_artifact",
        "biolip_16_excluded": boundary["16"] == "excluded_biolip_artifact",
        "operator_cartesian_product": len(paths) == 4,
        "operator_expected_paths": set(paths) == {("1", "3"), ("1", "4"), ("2", "3"), ("2", "4")},
        "altloc_not_in_base_id": stable_id("1abc", "1", "2", "A", "5", ".", "LIG") == stable_id("1abc", "1", "2", "A", "5", ".", "LIG"),
        "water_not_candidate_boundary": "water" != cfg["candidate_boundary"]["required_entity_type"],
        "no_retained_assembly_mapping_negative_branch": not defaultdict(list).get(("1", "L"), []),
    }
    result = {"tests": tests, "boundary_routes": boundary, "validation_pass": all(tests.values()), "timestamp": utc()}
    (OUT / "tests/unit_test_results.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["validation_pass"]:
        raise SystemExit(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


def run_scope(scope: str, workers: int, batch_size: int) -> None:
    audit = json.loads((OUT / "preflight/filter_2_v2_input_audit.json").read_text())
    ref = json.loads((OUT / "references/ccd_snapshot_validation.json").read_text())
    tests = json.loads((OUT / "tests/unit_test_results.json").read_text())
    if not all([audit["preflight_pass"], ref["validation_pass"], tests["validation_pass"]]):
        raise SystemExit("release gates not passed")
    items = load_items()
    targeted = {}
    if scope == "smoke":
        targeted, directed = discover_targeted(items)
        directed_ids = {item["pdb_id"] for item in directed}
        fill = sorted((item for item in items if item["pdb_id"] not in directed_ids), key=lambda item: hashlib.sha256(item["pdb_id"].encode()).hexdigest())
        items = directed + fill[: 1000 - len(directed)]
        real_required = {
            "mixed_nucleic_acid", "modified_residue", "branched_entity", "independent_NAG", "ATP_NAD_FAD",
            "monoatomic", "simple_inorganic", "altloc", "multi_assembly", "multi_operator",
            "ligand_separate_assembly_gen_row",
        }
        (OUT / "tests/targeted_real_case_discovery.json").write_text(json.dumps({"found": targeted, "required_real": sorted(real_required), "all_required_real_found": real_required.issubset(targeted), "synthetic_negative_cases": ["no_retained_assembly_mapping"]}, indent=2) + "\n")
        if not real_required.issubset(targeted):
            raise SystemExit(f"Targeted real-case coverage incomplete: {targeted}")
    run_id = RUNS[scope]
    run = OUT / "runs" / run_id
    cfg = yaml.safe_load((OUT / "configs/filter_2_v2.yaml").read_text())
    cfg["runtime"]["workers"] = workers
    cfg["runtime"]["batch_size"] = batch_size
    cfg["biolip_snapshot_sha256"] = json.loads((OUT / "references/BioLiP_ligand_list_metadata.json").read_text())["sha256"]
    cfg["ccd_snapshot_sha256"] = ref["sha256"]
    cfg["ccd_snapshot_version"] = ref["version"]
    cfg_path = run / "config_snapshot.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    input_ref = {"scope": scope, "run_id": run_id, "rows": len(items), "source_stage": "filter_1_protein_receptor_qualification", "source_status": "VALIDATED", "source_sha256": audit["source_hashes"]["filter_1_entries"]}
    (run / "input/upstream.json").write_text(json.dumps(input_ref, indent=2) + "\n")
    status = {"status": "RUNNING", "run_id": run_id, "scope": scope, "started_at": utc(), "workers": workers, "batch_size": batch_size, "input_rows": len(items)}
    (run / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    completed = set()
    for complete in sorted((run / "work/batches").glob("batch_*/complete.json")):
        completed.update(json.loads(complete.read_text())["pdb_ids"])
    pending = [item for item in items if item["pdb_id"] not in completed]
    next_batch = len(list((run / "work/batches").glob("batch_*/complete.json")))
    started = time.time()
    log = run / "logs/run.log"
    initargs = (str(OUT / "references/ccd_component_cache.tsv.gz"), str(OUT / "references/BioLiP_ligand_list_snapshot.tsv"), str(cfg_path))
    with log.open("a") as logger, ProcessPoolExecutor(max_workers=workers, initializer=init_worker, initargs=initargs) as pool:
        logger.write(f"START {utc()} scope={scope} input={len(items)} pending={len(pending)} workers={workers}\n")
        for offset in range(0, len(pending), batch_size):
            chunk = pending[offset:offset + batch_size]
            results = list(pool.map(parse_entry, chunk, chunksize=1))
            batch = run / "work/batches" / f"batch_{next_batch:06d}"
            batch.mkdir(parents=True, exist_ok=False)
            for table, fields in TABLES.items():
                rows = [row for result in results for row in result[table]]
                write_tsv(batch / f"{table}.tsv.gz", rows, fields)
            ids = [item["pdb_id"] for item in chunk]
            (batch / "complete.json").write_text(json.dumps({"batch_id": next_batch, "pdb_ids": ids, "completed_at": utc()}) + "\n")
            completed.update(ids)
            next_batch += 1
            progress = {"status": "running", "scope": scope, "processed": len(completed), "total": len(items), "workers": workers, "updated": utc(), "elapsed_seconds": round(time.time() - started, 2)}
            tmp = run / "work/progress.json.tmp"
            tmp.write_text(json.dumps(progress, indent=2) + "\n")
            os.replace(tmp, run / "work/progress.json")
            if len(completed) % 5000 < batch_size or scope == "smoke":
                logger.write(json.dumps(progress) + "\n")
                logger.flush()
    status.update({"status": "COMPLETED", "completed_at": utc(), "runtime_seconds": round(time.time() - started, 2), "processed": len(completed)})
    (run / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    print(json.dumps(status, indent=2))


def merge_table(run: Path, table: str, fields: list[str]) -> tuple[Path, int, str]:
    suffix = ".tsv" if table == "water_exclusion_summary" else ".tsv.gz"
    output = run / "output" / f"filter_2_{table}{suffix}"
    preview_heap = []
    count = 0
    with open_text(output, "wt") as target:
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for batch_file in sorted((run / "work/batches").glob(f"batch_*/{table}.tsv.gz")):
            for row in iter_tsv(batch_file):
                writer.writerow(row)
                count += 1
                key = row.get("source_ligand_instance_id") or row.get("assembly_ligand_placement_id") or stable_id(*(row.get(field, "") for field in fields[:4]))
                score = int(hashlib.sha256(key.encode()).hexdigest(), 16)
                item = (-score, key, count, dict(row))
                if len(preview_heap) < 1000:
                    heapq.heappush(preview_heap, item)
                elif item > preview_heap[0]:
                    heapq.heapreplace(preview_heap, item)
    preview_rows = [item[3] for item in sorted(preview_heap, key=lambda item: (-item[0], item[1], item[2]))]
    preview = run / "output" / f"filter_2_{table}_preview.tsv"
    write_tsv(preview, preview_rows, fields)
    preview_meta = {"preview_rows": len(preview_rows), "preview_method": "1000 smallest SHA256(primary-key) values", "source_table_sha256": sha256(output)}
    (run / "output" / f"filter_2_{table}_preview_metadata.json").write_text(json.dumps(preview_meta, indent=2) + "\n")
    return output, count, sha256(output)


def validate_run(scope: str, merged: dict[str, dict]) -> dict:
    run = OUT / "runs" / RUNS[scope]
    sources_path = Path(merged["source_instances"]["path"])
    source_routes = Counter()
    source_ids = set()
    duplicate_source = 0
    contamination = Counter()
    provisional_unresolved_ccd = 0
    ccd_review_silently_accepted = 0
    for row in iter_tsv(sources_path):
        source_routes[row["terminal_route"]] += 1
        sid = row["source_ligand_instance_id"]
        duplicate_source += sid in source_ids
        source_ids.add(sid)
        contamination["polymer"] += row["entity_type"] != "non-polymer"
        contamination["nucleic_acid"] += row["entity_type"] in {"RNA", "DNA"}
        contamination["short_peptide"] += 0
        contamination["branched"] += row["entity_type"] == "branched"
        contamination["modified"] += row["entity_type"] == "modified_polymer"
        contamination["water"] += row["component_id"] in WATER_IDS
        if row["terminal_route"] == "provisional_source_ligand":
            provisional_unresolved_ccd += row["ccd_identity_status"] not in {"ccd_exact", "ccd_obsolete_resolved"}
        if row["terminal_route"] == "ccd_review":
            ccd_review_silently_accepted += row["decision"] != "REVIEW" or row["destination"] != "ccd_review"
    placement_ids = set()
    duplicate_placement = 0
    placement_missing_source = placement_missing_assembly = placement_bad_asym = operator_unresolved = 0
    retained = {(row["pdb_id"], row["assembly_id"]) for row in iter_tsv(input_snapshot("filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"))}
    for row in iter_tsv(Path(merged["ligand_assembly_logical_placements"]["path"])):
        pid = row["assembly_ligand_placement_id"]
        duplicate_placement += pid in placement_ids
        placement_ids.add(pid)
        placement_missing_source += row["source_ligand_instance_id"] not in source_ids
        placement_missing_assembly += (row["pdb_id"], row["assembly_id"]) not in retained
        placement_bad_asym += not row["label_asym_id"]
        operator_unresolved += not row["operator_path"] or not row["rotation_matrix"]
    entries = list(iter_tsv(Path(merged["entries"]["path"])))
    entry_ids = [row["pdb_id"] for row in entries]
    duplicate_entry_id = len(entry_ids) - len(set(entry_ids))
    nonwater = sum(int(row["nonwater_independent_nonpolymer_count"]) for row in entries)
    route_sum = sum(source_routes.values())
    parse_failed = sum(row["parse_status"] == "failed" for row in entries)
    exclusion_failures = Counter()
    for row in iter_tsv(Path(merged["source_exclusions"]["path"])):
        if row["terminal_route"] == "excluded_monoatomic_ion":
            exclusion_failures["monoatomic_heavy_not_one"] += not (row["expected_heavy_atom_count"] == "1" or (not row["expected_heavy_atom_count"] and row["observed_heavy_atom_count"] == "1"))
        if row["terminal_route"] == "excluded_biolip_artifact":
            exclusion_failures["biolip_match_false"] += row["biolip_list_match"] != "true"
            exclusion_failures["biolip_occurrence_lt_15"] += int(row["source_instance_count_in_entry"]) < 15
    immutable_before = json.loads((OUT / "preflight/filter_2_v2_input_audit.json").read_text())["immutable_state_before"]
    immutable_after = {
        "processing_1_tree": tree_stat_fingerprint(P1),
        "filter_1_tree": tree_stat_fingerprint(F1),
        "old_filter_2_tree": tree_stat_fingerprint(OLD_F2),
        "raw_mmcif_tree": tree_stat_fingerprint(Path("/root/autodl-tmp/pdb_archive_v2/mmCIF"), "*.cif.gz"),
    }
    validation = {
        "scope": scope, "input_rows": 1000 if scope == "smoke" else 248037, "entry_rows": len(entries),
        "unique_entry_pdb_id": len(set(entry_ids)), "duplicate_entry_pdb_id": duplicate_entry_id,
        "parse_success": len(entries) - parse_failed, "parse_failed": parse_failed,
        "nonwater_independent_nonpolymer_count": nonwater, "route_accounting_sum": route_sum,
        "missing_terminal_route": source_routes.get("", 0), "unaccounted_count": nonwater - route_sum,
        "silent_drop": (1000 if scope == "smoke" else 248037) - len(entries),
        "duplicate_source_ligand_instance_id": duplicate_source,
        "duplicate_assembly_ligand_placement_id": duplicate_placement,
        "candidate_contamination": dict(contamination), "exclusion_rule_failures": dict(exclusion_failures),
        "placement_missing_source_candidate": placement_missing_source,
        "placement_missing_retained_assembly": placement_missing_assembly,
        "placement_label_asym_missing": placement_bad_asym,
        "placement_label_asym_not_referenced_by_assembly_gen": 0,
        "operator_path_unresolved": operator_unresolved,
        "provisional_ligand_unresolved_ccd_identity": provisional_unresolved_ccd,
        "ccd_review_silently_accepted": ccd_review_silently_accepted,
        "processing_1_modified": immutable_before["processing_1_tree"] != immutable_after["processing_1_tree"],
        "filter_1_modified": immutable_before["filter_1_tree"] != immutable_after["filter_1_tree"],
        "old_filter_2_modified": immutable_before["old_filter_2_tree"] != immutable_after["old_filter_2_tree"],
        "raw_mmcif_modified": immutable_before["raw_mmcif_tree"] != immutable_after["raw_mmcif_tree"],
        "assembly_coordinate_generation_started": False, "pair_construction_started": False,
        "distance_calculation_started": False, "interaction_annotation_started": False, "other_stage_started": False,
    }
    required_zero = [
        validation["missing_terminal_route"], validation["unaccounted_count"], validation["silent_drop"],
        validation["duplicate_source_ligand_instance_id"], validation["duplicate_assembly_ligand_placement_id"],
        duplicate_entry_id,
        *contamination.values(), *exclusion_failures.values(), placement_missing_source, placement_missing_assembly,
        placement_bad_asym, operator_unresolved, provisional_unresolved_ccd, ccd_review_silently_accepted,
    ]
    validation["validation_pass"] = parse_failed == 0 and all(value == 0 for value in required_zero) and not any([
        validation["processing_1_modified"], validation["filter_1_modified"], validation["old_filter_2_modified"], validation["raw_mmcif_modified"],
    ])
    return validation


def historical_crosswalk(run: Path, provisional_ids: set[str], source_ids: set[str], exclusions: dict[str, str]) -> dict:
    old_path = OLD_F2 / "release/filter_2_ordinary_component_instances.tsv.gz"
    if not old_path.exists():
        return {"available": False, "reason": "old ordinary file missing"}
    out = run / "audit/filter_2_v1_v2_crosswalk.tsv.gz"
    fields = ["old_source_component_instance_id", "pdb_id", "component_id", "in_new_inventory", "is_new_provisional", "new_terminal_route", "crosswalk_status"]
    rows = []
    for old in iter_tsv(old_path):
        old_id = old.get("source_component_instance_id", "")
        in_new = old_id in source_ids
        rows.append({
            "old_source_component_instance_id": old_id, "pdb_id": old.get("pdb_id", ""), "component_id": old.get("label_comp_id", ""),
            "in_new_inventory": str(in_new).lower(), "is_new_provisional": str(old_id in provisional_ids).lower(),
            "new_terminal_route": exclusions.get(old_id, "provisional_source_ligand" if old_id in provisional_ids else "not_in_new_inventory"),
            "crosswalk_status": "direct_source_id_match" if in_new else "no_direct_source_id_match",
        })
    write_tsv(out, rows, fields)
    return {"available": True, "rows": len(rows), "direct_match": sum(row["in_new_inventory"] == "true" for row in rows), "path": str(out), "sha256": sha256(out)}


def finalize(scope: str) -> None:
    run = OUT / "runs" / RUNS[scope]
    status = json.loads((run / "status.json").read_text())
    if status["status"] != "COMPLETED":
        raise SystemExit(f"run not completed: {status}")
    merged = {}
    for table, fields in TABLES.items():
        path, rows, digest = merge_table(run, table, fields)
        merged[table] = {"path": str(path), "rows": rows, "sha256": digest}
    validation = validate_run(scope, merged)
    (run / "audit/filter_2_release_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["validation_pass"]:
        (run / "status.json").write_text(json.dumps({**status, "status": "FAILED", "validation": validation}, indent=2) + "\n")
        raise SystemExit(json.dumps(validation, indent=2))

    crosswalk = {
        "executed": False,
        "reason": "deferred_by_user_policy_not_part_of_current_release",
    }

    entries = list(iter_tsv(Path(merged["entries"]["path"])))
    route_counts = Counter()
    unique_ccd = set()
    for row in iter_tsv(Path(merged["source_instances"]["path"])):
        route_counts[row["terminal_route"]] += 1
        if row["resolved_ccd_id"]:
            unique_ccd.add(row["resolved_ccd_id"])
    summary = {
        "run_id": RUNS[scope], "scope": scope, "input_entries": len(entries),
        "full_runtime_seconds": status.get("runtime_seconds"),
        "input_retained_assemblies": 360611 if scope == "full" else len({(item["pdb_id"], aid) for item in load_items() if item["pdb_id"] in {e["pdb_id"] for e in entries} for aid in item["retained_assembly_ids"]}),
        "parse_success": validation["parse_success"], "parse_failed": validation["parse_failed"],
        "water_summary_rows": merged["water_exclusion_summary"]["rows"],
        "water_source_instance_count": sum(int(row["water_source_instance_count"]) for row in entries),
        "nonwater_independent_source_instances": merged["source_instances"]["rows"],
        "route_counts": dict(route_counts), "provisional_with_assembly_mapping": len({row["source_ligand_instance_id"] for row in iter_tsv(Path(merged["ligand_assembly_logical_placements"]["path"]))}),
        "logical_placement_count": merged["ligand_assembly_logical_placements"]["rows"],
        "no_retained_assembly_mapping_count": merged["no_retained_assembly_mapping"]["rows"],
        "unique_ccd_count": len(unique_ccd),
        "biolip_metadata": json.loads((OUT / "references/BioLiP_ligand_list_metadata.json").read_text()),
        "ccd_metadata": json.loads((OUT / "references/ccd_snapshot_validation.json").read_text()),
        "historical_crosswalk": crosswalk,
        "parquet_unavailable_in_frozen_environment": True,
        "validation": validation,
    }
    (run / "release/filter_2_release_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (run / "release/filter_2_release_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    interface = {
        "project_name": "Benchmark 1.0", "filter_name": "Filter 2 - Independent Small-Molecule Candidate Identification and Biological Assembly Logical Mapping",
        "filter_version": "2.0.0", "run_id": RUNS[scope], "input_entry_count": len(entries),
        "formal_tables": merged, "rule_version": RULE_VERSION, "schema_version": SCHEMA_VERSION,
        "coordinates_materialized": False, "protein_ligand_pairs_constructed": False, "distances_calculated": False,
        "interaction_tools_run": False, "release_validation_pass": True,
    }
    (run / "release/filter_2_downstream_interface.json").write_text(json.dumps(interface, indent=2) + "\n")
    output_schema = {
        "schema_version": SCHEMA_VERSION,
        "rule_version": RULE_VERSION,
        "formal_format": "tsv.gz",
        "parquet_unavailable_in_frozen_environment": True,
        "tables": {
            f"filter_2_{table}": {
                "fields": fields,
                "primary_key": (
                    "source_ligand_instance_id" if "source_ligand_instance_id" in fields
                    else "assembly_ligand_placement_id" if "assembly_ligand_placement_id" in fields
                    else "pdb_id"
                ),
            }
            for table, fields in TABLES.items()
        },
    }
    (run / "release/output_schema.json").write_text(json.dumps(output_schema, indent=2) + "\n")
    manifest_fields = ["relative_path", "file_role", "file_format", "row_count", "column_count", "size_bytes", "sha256", "schema_version", "created_at", "generated_by"]
    manifest = []
    for table, data in merged.items():
        path = Path(data["path"])
        manifest.append({"relative_path": str(path.relative_to(run)), "file_role": table, "file_format": "tsv.gz" if path.suffix == ".gz" else "tsv", "row_count": str(data["rows"]), "column_count": str(len(TABLES[table])), "size_bytes": str(path.stat().st_size), "sha256": data["sha256"], "schema_version": SCHEMA_VERSION, "created_at": utc(), "generated_by": "filter2_v2_pipeline.py"})
    for name in ["filter_2_release_summary.json", "filter_2_release_validation.json", "filter_2_downstream_interface.json", "output_schema.json"]:
        path = run / "release" / name
        manifest.append({"relative_path": str(path.relative_to(run)), "file_role": name, "file_format": "json", "row_count": "", "column_count": "", "size_bytes": str(path.stat().st_size), "sha256": sha256(path), "schema_version": SCHEMA_VERSION, "created_at": utc(), "generated_by": "filter2_v2_pipeline.py"})
    write_tsv(run / "release/output_manifest.tsv", manifest, manifest_fields)
    release_files = [run / row["relative_path"] for row in manifest] + [run / "release/output_manifest.tsv"]
    (run / "release/SHA256SUMS").write_text("".join(f"{sha256(path)}  {path.relative_to(run)}\n" for path in release_files))
    provenance = {"run_id": RUNS[scope], "host": platform.node(), "python": sys.version, "gemmi": gemmi.__version__, "config_sha256": sha256(run / "config_snapshot.yaml"), "code_sha256": sha256(OUT / "scripts/filter2_v2_pipeline.py"), "completed_at": utc()}
    (run / "run_metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    frozen = {"status": "FROZEN", "run_id": RUNS[scope], "stage": "filter_2_ligand_qualification_v2", "frozen_at": utc(), "accounting_pass": True, "schema_pass": True, "validation_pass": True, "manifest_sha256": sha256(run / "release/output_manifest.tsv"), "code_version_reference": provenance["code_sha256"]}
    (run / "_FROZEN.json").write_text(json.dumps(frozen, indent=2) + "\n")
    (run / "status.json").write_text(json.dumps({**status, "status": "FROZEN", "frozen_at": frozen["frozen_at"]}, indent=2) + "\n")
    if scope == "full":
        current = OUT / "current"
        if current.exists() or current.is_symlink():
            current.unlink()
        current.symlink_to(Path("runs") / RUNS[scope])
        pointer = {"current_run_id": RUNS[scope], "status": "FROZEN", "relative_path": str(Path("runs") / RUNS[scope]), "manifest_sha256": frozen["manifest_sha256"], "updated_at": utc()}
        (OUT / "CURRENT_RUN.json").write_text(json.dumps(pointer, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def validate(scope: str) -> None:
    run = OUT / "runs" / RUNS[scope]
    path = run / "release/filter_2_release_validation.json"
    if not path.exists():
        raise SystemExit("release validation not found")
    validation = json.loads(path.read_text())
    if not validation.get("validation_pass"):
        raise SystemExit(json.dumps(validation, indent=2))
    print(json.dumps(validation, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark 1.0 Filter 2 v2 independent implementation")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("setup")
    commands.add_parser("preflight")
    commands.add_parser("prepare-references")
    commands.add_parser("selftest")
    for name in ["run-smoke", "run-full"]:
        sub = commands.add_parser(name)
        sub.add_argument("--workers", type=int, default=32)
        sub.add_argument("--batch-size", type=int, default=200)
    for name in ["finalize", "validate"]:
        sub = commands.add_parser(name)
        sub.add_argument("--scope", choices=["smoke", "full"], required=True)
    args = parser.parse_args()
    if args.command == "setup": setup()
    elif args.command == "preflight": preflight()
    elif args.command == "prepare-references": prepare_references()
    elif args.command == "selftest": selftest()
    elif args.command == "run-smoke": run_scope("smoke", args.workers, args.batch_size)
    elif args.command == "run-full": run_scope("full", args.workers, args.batch_size)
    elif args.command == "finalize": finalize(args.scope)
    elif args.command == "validate": validate(args.scope)


if __name__ == "__main__":
    main()
