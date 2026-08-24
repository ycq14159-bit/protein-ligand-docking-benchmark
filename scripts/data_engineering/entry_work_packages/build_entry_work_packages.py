#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import multiprocessing
import os
import platform
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import yaml

ROOT = Path("/root/autodl-tmp/benchmark_1.0/auxiliary_entry_work_packages")
BUILD_ID = "20260805_full_01"
BUILD = ROOT / "builds" / BUILD_ID
F2 = Path("/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v3/runs/20260804_full_01")
F1 = Path("/root/autodl-tmp/benchmark_1.0/filter_1_protein_receptor_qualification")
P1 = Path("/root/autodl-tmp/benchmark_1.0/processing_1_pdb_source_audit/release/processing_1_mmcif_index.tsv.gz")

P_PLACEMENTS = F2 / "output/ligand_assembly_logical_placements.tsv.gz"
P_PROVISIONAL = F2 / "output/provisional_source_ligands.tsv.gz"
P_ENTRIES = F2 / "output/entry_route_summary.tsv.gz"
P_F1_RECEPTORS = F1 / "release/filter_1_receptor_chain_instances.tsv.gz"
P_F1_SOURCE_CHAINS = F1 / "full/filter_1_source_chains.tsv.gz"
P_F1_ASSEMBLIES = F1 / "release/filter_1_receptor_qualified_assemblies.tsv.gz"

DATASETS = [
    "entry_ligand_placements", "entry_receptor_chain_instances", "entry_assembly_context",
    "entry_ligand_source_atoms", "entry_receptor_source_atoms",
]

G_CFG = {}
G_P1 = {}
G_ENTRY_IDS = []
G_PLACEMENTS = defaultdict(list)
G_SOURCES = defaultdict(dict)
G_RECEPTORS = defaultdict(list)
G_SOURCE_CHAIN = {}


def utc():
    return datetime.now(timezone.utc).isoformat()


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value in {".", "?", "None"} else value


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bucket_id(pdb_id, count=256):
    return int(hashlib.sha256(pdb_id.lower().encode()).hexdigest()[:8], 16) % count


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith(".gz") else Path(path).open(encoding="utf-8", newline="")


def iter_tsv(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def category_records(block, prefix):
    table = block.find_mmcif_category(prefix)
    if not table:
        return []
    tags = [tag[len(prefix):] for tag in table.tags]
    return [{tag: clean(value) for tag, value in zip(tags, row)} for row in table]


def expand_token(token):
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
        raise ValueError("empty operator token")
    return values


def expand_operator_paths(expression):
    compact = clean(expression).strip()
    if compact.startswith(";"):
        compact = compact[1:]
    if compact.endswith(";"):
        compact = compact[:-1]
    compact = re.sub(r"\s+", "", compact)
    groups = re.findall(r"\(([^()]*)\)", compact) or [compact]
    return [tuple(path) for path in itertools.product(*(expand_token(group) for group in groups))]


def identity_affine():
    return ([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0., 0., 0.])


def compose_affine(left, right):
    lr, lt = left; rr, rt = right
    rotation = [[sum(lr[i][k] * rr[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    translation = [sum(lr[i][k] * rt[k] for k in range(3)) + lt[i] for i in range(3)]
    return rotation, translation


def composite_affine(path, operators):
    total = identity_affine()
    for operator_id in path:
        total = compose_affine(total, operators[operator_id])
    return total


def parse_operator_records(block):
    result = {}
    for row in category_records(block, "_pdbx_struct_oper_list."):
        oid = row.get("id", "")
        try:
            matrix = [[float(row[f"matrix[{i}][{j}]"]) for j in range(1, 4)] for i in range(1, 4)]
            vector = [float(row[f"vector[{i}]"]) for i in range(1, 4)]
            if oid and all(math.isfinite(x) for line in matrix for x in line) and all(math.isfinite(x) for x in vector):
                result[oid] = (matrix, vector)
        except Exception:
            continue
    return result


def matrix_fields(rotation, translation):
    values = {f"r{i}{j}": float(rotation[i - 1][j - 1]) for i in range(1, 4) for j in range(1, 4)}
    values.update({f"t{i}": float(translation[i - 1]) for i in range(1, 4)})
    return values


def parse_matrix_string(value):
    return [[float(x) for x in row.split(",")] for row in value.split(";")]


def parse_vector_string(value):
    return [float(x) for x in value.split(",")]


def close_enough(a, b):
    cfg = G_CFG["numeric_tolerance"]
    return len(a) == len(b) and all(math.isclose(x, y, abs_tol=float(cfg["absolute"]), rel_tol=float(cfg["relative"])) for x, y in zip(a, b))


def as_int(value):
    try:
        return int(value)
    except Exception:
        return None


def as_float(value):
    try:
        answer = float(value)
        return answer if math.isfinite(answer) else None
    except Exception:
        return None


def normalized_insertion(value):
    value = clean(value)
    return "" if value.lower() in {"false", "none", "null"} else value


def parse_asym_id_list(value):
    text = clean(value).strip()
    if text.startswith(";"):
        text = text[1:]
    if text.endswith(";"):
        text = text[:-1]
    return [token.strip() for token in text.replace("\n", "").split(",") if token.strip()]


def schemas():
    s = pa.string(); i64 = pa.int64(); i16 = pa.int16(); f64 = pa.float64(); b = pa.bool_()
    matrix = [(f"r{x}{y}", f64) for x in range(1, 4) for y in range(1, 4)] + [(f"t{x}", f64) for x in range(1, 4)]
    return {
        "entry_ligand_placements": pa.schema([
            ("pdb_id", s), ("bucket_id", i16), ("assembly_id", s), ("model_id", s),
            ("filter_2_ligand_assembly_placement_id", s), ("filter_2_source_ligand_instance_id", s),
            ("component_id", s), ("entity_id", s), ("label_asym_id", s), ("auth_asym_id", s),
            ("auth_seq_id", s), ("insertion_code", s), ("assembly_gen_row_id", s),
            ("oper_expression_raw", s), ("operator_path", s), ("rotation_matrix", s),
            ("translation_vector", s), *matrix, ("ccd_identity_status", s),
            ("standard_total_atom_count", i64), ("expected_heavy_atom_count", i64),
            ("observed_heavy_atom_count", i64), ("formal_charge", i64), ("element_set", s),
        ]),
        "entry_receptor_chain_instances": pa.schema([
            ("pdb_id", s), ("bucket_id", i16), ("assembly_id", s), ("model_id", s),
            ("filter_1_chain_instance_id", s), ("filter_1_source_chain_key", s),
            ("entity_id", s), ("label_asym_id", s), ("auth_asym_id", s), ("polymer_class", s),
            ("assembly_gen_row_id", s), ("oper_expression_raw", s), ("operator_path", s),
            ("rotation_matrix", s), ("translation_vector", s), *matrix,
            ("declared_sequence_length", i64), ("observed_residue_count", i64),
            ("observed_fraction", f64), ("receptor_eligible", b), ("chain_role", s),
        ]),
        "entry_assembly_context": pa.schema([
            ("pdb_id", s), ("bucket_id", i16), ("assembly_id", s),
            ("assembly_details", s), ("assembly_method_details", s), ("oligomeric_details", s),
            ("oligomeric_count", i64), ("assembly_gen_row_id", s), ("oper_expression_raw", s),
            ("asym_id_list_raw", s), ("operator_path", s), ("composite_operator_id", s),
            ("rotation_matrix", s), ("translation_vector", s), *matrix,
        ]),
        "entry_ligand_source_atoms": pa.schema([
            ("pdb_id", s), ("bucket_id", i16), ("model_id", s),
            ("filter_2_source_ligand_instance_id", s), ("component_id", s), ("entity_id", s),
            ("label_asym_id", s), ("auth_asym_id", s), ("label_seq_id", s), ("auth_seq_id", s),
            ("insertion_code", s), ("source_atom_row_index", i64), ("atom_site_id", s),
            ("atom_pointer_status", s), ("group_PDB", s), ("label_entity_id", s),
            ("label_atom_id", s), ("auth_atom_id", s), ("type_symbol", s), ("alt_id", s),
            ("Cartn_x", f64), ("Cartn_y", f64), ("Cartn_z", f64), ("occupancy", f64),
            ("B_iso_or_equiv", f64), ("formal_charge_if_present", s),
        ]),
        "entry_receptor_source_atoms": pa.schema([
            ("pdb_id", s), ("bucket_id", i16), ("model_id", s), ("filter_1_source_chain_key", s),
            ("entity_id", s), ("label_asym_id", s), ("auth_asym_id", s),
            ("label_seq_id", s), ("auth_seq_id", s), ("insertion_code", s),
            ("label_comp_id", s), ("auth_comp_id", s), ("source_atom_row_index", i64),
            ("atom_site_id", s), ("atom_pointer_status", s), ("group_PDB", s),
            ("label_atom_id", s), ("auth_atom_id", s), ("type_symbol", s), ("alt_id", s),
            ("Cartn_x", f64), ("Cartn_y", f64), ("Cartn_z", f64), ("occupancy", f64),
            ("B_iso_or_equiv", f64),
        ]),
    }


SCHEMAS = schemas()


def write_parquet_atomic(path, rows, schema):
    if not rows:
        return 0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    table = pa.Table.from_pylist(rows, schema=schema)
    pq.write_table(table, tmp, compression=G_CFG["compression"], compression_level=int(G_CFG["compression_level"]),
                   version=G_CFG["parquet_version"], data_page_version=G_CFG["parquet_data_page_version"],
                   use_dictionary=True, write_statistics=True)
    tmp.replace(path)
    return table.num_rows


def setup():
    if ROOT.exists():
        raise SystemExit(f"target already exists: {ROOT}")
    for path in [ROOT / "configs", ROOT / "scripts", BUILD / "work/bucket_fragments", BUILD / "work/checkpoints",
                 BUILD / "tmp", BUILD / "logs", BUILD / "output", BUILD / "audit"]:
        path.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__), ROOT / "scripts/build_entry_work_packages.py")
    source_dir = Path(__file__).parent
    shutil.copy2(source_dir / "entry_work_packages_config.yaml", ROOT / "configs/build_v1.yaml")
    shutil.copy2(source_dir / "load_entry_work_package.py", ROOT / "scripts/load_entry_work_package.py")
    shutil.copy2(source_dir / "entry_work_packages_README.md", ROOT / "README.md")
    print(ROOT)


def input_paths():
    return {
        "filter2_frozen": F2 / "_FROZEN.json", "filter2_interface": F2 / "release/downstream_interface.json",
        "filter2_manifest": F2 / "release/output_manifest.tsv", "filter2_sha256sums": F2 / "release/SHA256SUMS",
        "filter2_placements": P_PLACEMENTS, "filter2_provisional": P_PROVISIONAL, "filter2_entries": P_ENTRIES,
        "filter1_summary": F1 / "release/filter_1_release_summary.json", "filter1_validation": F1 / "release/filter_1_release_validation.json",
        "filter1_receptors": P_F1_RECEPTORS, "filter1_source_chains": P_F1_SOURCE_CHAINS,
        "filter1_assemblies": P_F1_ASSEMBLIES, "processing1_index": P1,
    }


def preflight():
    global G_CFG
    G_CFG = yaml.safe_load((ROOT / "configs/build_v1.yaml").read_text())
    paths = input_paths()
    missing = [name for name, path in paths.items() if not path.exists()]
    if missing:
        raise SystemExit(f"missing inputs: {missing}")
    frozen = json.loads(paths["filter2_frozen"].read_text())
    interface = json.loads(paths["filter2_interface"].read_text())
    f1_validation = json.loads(paths["filter1_validation"].read_text())
    counts = {"placements": 0, "mapped_sources": set(), "active_pdb": set(), "active_keys": set()}
    for row in iter_tsv(P_PLACEMENTS):
        counts["placements"] += 1
        counts["mapped_sources"].add(row["source_ligand_instance_id"])
        counts["active_pdb"].add(row["pdb_id"])
        counts["active_keys"].add((row["pdb_id"], row["assembly_id"], row["selected_model_id"]))
    entry_count = sum(1 for _ in iter_tsv(P_ENTRIES))
    audit = {
        "timestamp": utc(), "filter2_status": frozen.get("status"),
        "filter2_validation_pass": frozen.get("validation_pass"),
        "filter1_validation_pass": f1_validation.get("release_validation_pass"),
        "input_entry_count": entry_count, "ligand_placement_count": counts["placements"],
        "mapped_source_ligand_count": len(counts["mapped_sources"]), "active_pdb_count": len(counts["active_pdb"]),
        "active_assembly_key_count": len(counts["active_keys"]),
        "expected_no_mapping_source_count": 1002,
        "input_sha256": {name: sha256(path) for name, path in paths.items()},
        "config_sha256": sha256(ROOT / "configs/build_v1.yaml"),
        "code_sha256": sha256(ROOT / "scripts/build_entry_work_packages.py"),
        "pyarrow_version": pa.__version__, "python_version": sys.version,
    }
    audit["preflight_pass"] = all([
        audit["filter2_status"] == "FROZEN", audit["filter2_validation_pass"] is True,
        audit["filter1_validation_pass"] is True, entry_count == 248037,
        counts["placements"] == 1151324, len(counts["mapped_sources"]) == 851966,
        interface.get("coordinates_materialized") is False, interface.get("pairs_constructed") is False,
    ])
    (BUILD / "audit/preflight.json").write_text(json.dumps(audit, indent=2) + "\n")
    (BUILD / "output/input_provenance.json").write_text(json.dumps(audit, indent=2) + "\n")
    if not audit["preflight_pass"]:
        raise SystemExit(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


def load_inputs():
    global G_CFG, G_P1, G_ENTRY_IDS, G_PLACEMENTS, G_SOURCES, G_RECEPTORS, G_SOURCE_CHAIN
    G_CFG = yaml.safe_load((ROOT / "configs/build_v1.yaml").read_text())
    G_P1 = {r["pdb_id"]: r for r in iter_tsv(P1)}
    G_ENTRY_IDS = [r["pdb_id"] for r in iter_tsv(P_ENTRIES)]
    mapped_ids = set()
    active_keys = set()
    for row in iter_tsv(P_PLACEMENTS):
        pid = row["pdb_id"]
        G_PLACEMENTS[pid].append(row)
        mapped_ids.add(row["source_ligand_instance_id"])
        active_keys.add((pid, row["assembly_id"], row["selected_model_id"]))
    for row in iter_tsv(P_PROVISIONAL):
        if row["source_ligand_instance_id"] in mapped_ids:
            G_SOURCES[row["pdb_id"]][row["source_ligand_instance_id"]] = row
    active_source_chain_keys = set()
    for row in iter_tsv(P_F1_RECEPTORS):
        key = (row["pdb_id"], row["assembly_id"], row["model_id"])
        if key in active_keys:
            G_RECEPTORS[row["pdb_id"]].append(row)
            active_source_chain_keys.add((row["pdb_id"], row["entity_id"], row["label_asym_id"]))
    for row in iter_tsv(P_F1_SOURCE_CHAINS):
        key = (row["pdb_id"], row["entity_id"], row["label_asym_id"])
        if key in active_source_chain_keys:
            G_SOURCE_CHAIN[key] = row
    return active_keys


def atom_common(pid, bid, index, atom, id_counts):
    atom_id = atom.get("id", "")
    pointer = "atom_site_id_unique" if atom_id and id_counts[atom_id] == 1 else "source_atom_row_index_fallback"
    return {
        "pdb_id": pid, "bucket_id": bid, "source_atom_row_index": index,
        "atom_site_id": atom_id, "atom_pointer_status": pointer, "group_PDB": atom.get("group_PDB", ""),
        "label_atom_id": atom.get("label_atom_id", ""), "auth_atom_id": atom.get("auth_atom_id", ""),
        "type_symbol": atom.get("type_symbol", "").upper(), "alt_id": atom.get("label_alt_id", ""),
        "Cartn_x": as_float(atom.get("Cartn_x")), "Cartn_y": as_float(atom.get("Cartn_y")),
        "Cartn_z": as_float(atom.get("Cartn_z")), "occupancy": as_float(atom.get("occupancy")),
        "B_iso_or_equiv": as_float(atom.get("B_iso_or_equiv")),
    }


def parse_entry(pid):
    bid = bucket_id(pid, int(G_CFG["bucket_count"]))
    output = {name: [] for name in DATASETS}
    exceptions = []
    placements = G_PLACEMENTS[pid]
    sources = G_SOURCES[pid]
    receptors = G_RECEPTORS[pid]
    path = G_P1[pid]["mmcif_path"]
    try:
        block = gemmi.cif.read(path).sole_block()
        operators = parse_operator_records(block)
        assembly_meta = {r.get("id", ""): r for r in category_records(block, "_pdbx_struct_assembly.")}
        active_aids = {r["assembly_id"] for r in placements}
        contexts = {}
        receptor_refs = defaultdict(list)
        for index, gen in enumerate(category_records(block, "_pdbx_struct_assembly_gen."), start=1):
            aid = gen.get("assembly_id", "")
            if aid not in active_aids:
                continue
            row_id = f"assembly_gen_row_{index:06d}"
            expr = gen.get("oper_expression", "")
            asyms = parse_asym_id_list(gen.get("asym_id_list", ""))
            try:
                paths = expand_operator_paths(expr)
            except Exception as exc:
                exceptions.append(("operator_reference_failure", f"{aid}:{row_id}:{exc}")); paths = []
            meta = assembly_meta.get(aid, {})
            for opath in paths:
                if any(op not in operators for op in opath):
                    exceptions.append(("operator_reference_failure", f"{aid}:{row_id}:{opath}")); continue
                rotation, translation = composite_affine(opath, operators)
                op_string = "*".join(opath)
                m = matrix_fields(rotation, translation)
                context = {
                    "pdb_id": pid, "bucket_id": bid, "assembly_id": aid,
                    "assembly_details": meta.get("details", ""), "assembly_method_details": meta.get("method_details", ""),
                    "oligomeric_details": meta.get("oligomeric_details", ""), "oligomeric_count": as_int(meta.get("oligomeric_count")),
                    "assembly_gen_row_id": row_id, "oper_expression_raw": expr,
                    "asym_id_list_raw": gen.get("asym_id_list", ""), "operator_path": op_string,
                    "composite_operator_id": op_string,
                    "rotation_matrix": ";".join(",".join(f"{x:.10g}" for x in line) for line in rotation),
                    "translation_vector": ",".join(f"{x:.10g}" for x in translation), **m,
                }
                output["entry_assembly_context"].append(context)
                contexts[(aid, row_id, op_string)] = context
                for asym in asyms:
                    receptor_refs[(aid, asym, "x".join(opath))].append(context)

        for placement in placements:
            sid = placement["source_ligand_instance_id"]
            source = sources.get(sid)
            if not source:
                exceptions.append(("official_id_missing", sid)); continue
            ckey = (placement["assembly_id"], placement["assembly_gen_row_id"], placement["operator_path"])
            context = contexts.get(ckey)
            if not context:
                exceptions.append(("operator_reference_failure", placement["assembly_ligand_placement_id"])); continue
            expected_r = [x for row in parse_matrix_string(placement["rotation_matrix"]) for x in row]
            actual_r = [context[f"r{i}{j}"] for i in range(1, 4) for j in range(1, 4)]
            if not close_enough(expected_r, actual_r) or not close_enough(parse_vector_string(placement["translation_vector"]), [context[f"t{i}"] for i in range(1, 4)]):
                exceptions.append(("matrix_mismatch", placement["assembly_ligand_placement_id"]))
            output["entry_ligand_placements"].append({
                "pdb_id": pid, "bucket_id": bid, "assembly_id": placement["assembly_id"],
                "model_id": placement["selected_model_id"],
                "filter_2_ligand_assembly_placement_id": placement["assembly_ligand_placement_id"],
                "filter_2_source_ligand_instance_id": sid, "component_id": placement["component_id"],
                "entity_id": source["entity_id"], "label_asym_id": source["label_asym_id"],
                "auth_asym_id": source["auth_asym_id"], "auth_seq_id": source["auth_seq_id"],
                "insertion_code": source["insertion_code"], "assembly_gen_row_id": placement["assembly_gen_row_id"],
                "oper_expression_raw": placement["oper_expression_raw"], "operator_path": placement["operator_path"],
                "rotation_matrix": placement["rotation_matrix"], "translation_vector": placement["translation_vector"],
                **{name: context[name] for name in [f"r{i}{j}" for i in range(1, 4) for j in range(1, 4)] + [f"t{i}" for i in range(1, 4)]},
                "ccd_identity_status": source["ccd_identity_status"],
                "standard_total_atom_count": as_int(source["standard_total_atom_count"]),
                "expected_heavy_atom_count": as_int(source["expected_heavy_atom_count"]),
                "observed_heavy_atom_count": as_int(source["observed_heavy_atom_count"]),
                "formal_charge": as_int(source["formal_charge"]), "element_set": source["element_set"],
            })

        for receptor in receptors:
            refs = receptor_refs.get((receptor["assembly_id"], receptor["label_asym_id"], receptor["operator_id"]), [])
            if not refs:
                exceptions.append(("operator_reference_failure", receptor["chain_instance_id"])); continue
            context = sorted(refs, key=lambda x: (x["assembly_gen_row_id"], x["operator_path"]))[0]
            source_key = "|".join([pid, receptor["model_id"], receptor["entity_id"], receptor["label_asym_id"]])
            source_chain = G_SOURCE_CHAIN.get((pid, receptor["entity_id"], receptor["label_asym_id"]), {})
            output["entry_receptor_chain_instances"].append({
                "pdb_id": pid, "bucket_id": bid, "assembly_id": receptor["assembly_id"], "model_id": receptor["model_id"],
                "filter_1_chain_instance_id": receptor["chain_instance_id"], "filter_1_source_chain_key": source_key,
                "entity_id": receptor["entity_id"], "label_asym_id": receptor["label_asym_id"], "auth_asym_id": receptor["auth_asym_id"],
                "polymer_class": receptor["polymer_class"], "assembly_gen_row_id": context["assembly_gen_row_id"],
                "oper_expression_raw": context["oper_expression_raw"], "operator_path": context["operator_path"],
                "rotation_matrix": context["rotation_matrix"], "translation_vector": context["translation_vector"],
                **{name: context[name] for name in [f"r{i}{j}" for i in range(1, 4) for j in range(1, 4)] + [f"t{i}" for i in range(1, 4)]},
                "declared_sequence_length": as_int(receptor["declared_sequence_length"]),
                "observed_residue_count": as_int(receptor["observed_residue_count"]),
                "observed_fraction": as_float(source_chain.get("observed_fraction")),
                "receptor_eligible": receptor["receptor_eligible"].lower() == "true", "chain_role": receptor["chain_role"],
            })

        entities = {r.get("id", ""): r.get("type", "").lower() for r in category_records(block, "_entity.")}
        asym_entity = {r.get("id", ""): r.get("entity_id", "") for r in category_records(block, "_struct_asym.")}
        scheme = defaultdict(list)
        for row in category_records(block, "_pdbx_nonpoly_scheme."):
            scheme[(row.get("asym_id", ""), row.get("mon_id", "").upper())].append(row)
        atoms = category_records(block, "_atom_site.")
        id_counts = Counter(a.get("id", "") for a in atoms if a.get("id", ""))
        coarse_sources = {(r["selected_model_id"], r["entity_id"], r["label_asym_id"], r["component_id"]) for r in sources.values()}
        candidate_lookup = {
            (r["selected_model_id"], r["entity_id"], r["label_asym_id"], r["auth_seq_id"],
             normalized_insertion(r["insertion_code"]), r["component_id"]): sid
            for sid, r in sources.items()
        }
        groups = defaultdict(list)
        receptor_atom_keys = {(r["model_id"], r["entity_id"], r["label_asym_id"]): "|".join([pid, r["model_id"], r["entity_id"], r["label_asym_id"]]) for r in receptors}
        for index, atom in enumerate(atoms, start=1):
            model = atom.get("pdbx_PDB_model_num", "") or "1"
            asym = atom.get("label_asym_id", "")
            entity = asym_entity.get(asym, atom.get("label_entity_id", ""))
            component = atom.get("label_comp_id", "").upper()
            coarse = (model, entity, asym, component)
            if coarse in coarse_sources:
                groups[(model, entity, asym, atom.get("auth_seq_id", ""), atom.get("pdbx_PDB_ins_code", ""), component)].append((index, atom))
            rkey = (model, entity, asym)
            if rkey in receptor_atom_keys:
                common = atom_common(pid, bid, index, atom, id_counts)
                if any(common[name] is None for name in ("Cartn_x", "Cartn_y", "Cartn_z")):
                    exceptions.append(("invalid_coordinate", f"receptor:{receptor_atom_keys[rkey]}:{index}"))
                common.update({
                    "model_id": model, "filter_1_source_chain_key": receptor_atom_keys[rkey], "entity_id": entity,
                    "label_asym_id": asym, "auth_asym_id": atom.get("auth_asym_id", ""),
                    "label_seq_id": atom.get("label_seq_id", ""), "auth_seq_id": atom.get("auth_seq_id", ""),
                    "insertion_code": atom.get("pdbx_PDB_ins_code", ""), "label_comp_id": component,
                    "auth_comp_id": atom.get("auth_comp_id", ""),
                })
                output["entry_receptor_source_atoms"].append(common)

        found_sources = set()
        for key, group in groups.items():
            model, entity, asym, auth_seq, insertion, component = key
            scheme_rows = scheme.get((asym, component), [])
            chosen = next((r for r in scheme_rows if r.get("auth_seq_num", "") == auth_seq or r.get("pdb_seq_num", "") == auth_seq), scheme_rows[0] if len(scheme_rows) == 1 else {})
            pdb_seq = chosen.get("pdb_seq_num", "") or auth_seq
            auth_seq_final = chosen.get("auth_seq_num", "") or auth_seq
            insertion_final = chosen.get("pdb_ins_code", "") or insertion
            sid = candidate_lookup.get((model, entity, asym, auth_seq_final or pdb_seq, normalized_insertion(insertion_final), component))
            if sid is None:
                continue
            found_sources.add(sid)
            formal_source = sources[sid]
            for index, atom in group:
                common = atom_common(pid, bid, index, atom, id_counts)
                if any(common[name] is None for name in ("Cartn_x", "Cartn_y", "Cartn_z")):
                    exceptions.append(("invalid_coordinate", f"ligand:{sid}:{index}"))
                common.update({
                    "model_id": model, "filter_2_source_ligand_instance_id": sid, "component_id": component,
                    "entity_id": entity, "label_asym_id": asym, "auth_asym_id": atom.get("auth_asym_id", ""),
                    "label_seq_id": atom.get("label_seq_id", ""), "auth_seq_id": auth_seq_final,
                    "insertion_code": formal_source["insertion_code"], "label_entity_id": atom.get("label_entity_id", entity),
                    "formal_charge_if_present": atom.get("pdbx_formal_charge", ""),
                })
                output["entry_ligand_source_atoms"].append(common)
        for sid in set(sources) - found_sources:
            exceptions.append(("ligand_instance_atom_match_failure", sid))

        receptor_source_expected = set(receptor_atom_keys.values())
        receptor_source_found = {r["filter_1_source_chain_key"] for r in output["entry_receptor_source_atoms"]}
        for source_key in receptor_source_expected - receptor_source_found:
            exceptions.append(("receptor_chain_atom_match_failure", source_key))

        for name, fields in {
            "entry_ligand_placements": ("pdb_id", "assembly_id", "model_id", "filter_2_ligand_assembly_placement_id"),
            "entry_receptor_chain_instances": ("pdb_id", "assembly_id", "model_id", "filter_1_chain_instance_id"),
            "entry_assembly_context": ("pdb_id", "assembly_id", "assembly_gen_row_id", "operator_path"),
            "entry_ligand_source_atoms": ("pdb_id", "filter_2_source_ligand_instance_id", "source_atom_row_index"),
            "entry_receptor_source_atoms": ("pdb_id", "filter_1_source_chain_key", "source_atom_row_index"),
        }.items():
            output[name].sort(key=lambda row: tuple(row[x] for x in fields))

        ligand_counts = Counter(r["filter_2_source_ligand_instance_id"] for r in output["entry_ligand_source_atoms"])
        receptor_counts = Counter(r["filter_1_source_chain_key"] for r in output["entry_receptor_source_atoms"])
        by_assembly_lig = Counter((r["assembly_id"], r["model_id"]) for r in output["entry_ligand_placements"])
        by_assembly_rec_atoms = Counter()
        transform = 0
        for r in output["entry_ligand_placements"]:
            transform += ligand_counts[r["filter_2_source_ligand_instance_id"]]
        for r in output["entry_receptor_chain_instances"]:
            count = receptor_counts[r["filter_1_source_chain_key"]]
            by_assembly_rec_atoms[(r["assembly_id"], r["model_id"])] += count
            transform += count
        pair_work = sum(by_assembly_lig[key] * by_assembly_rec_atoms[key] for key in by_assembly_lig)
        stats = {
            "pdb_id": pid, "bucket_id": bid, "mmcif_path": path, "has_active_work": True,
            "active_assembly_count": len({r["assembly_id"] for r in output["entry_ligand_placements"]}),
            "active_model_count": len({r["model_id"] for r in output["entry_ligand_placements"]}),
            "mapped_source_ligand_count": len(found_sources), "ligand_placement_count": len(output["entry_ligand_placements"]),
            "active_receptor_source_chain_count": len(receptor_source_found),
            "active_receptor_chain_instance_count": len(output["entry_receptor_chain_instances"]),
            "ligand_source_atom_count": len(output["entry_ligand_source_atoms"]),
            "receptor_source_atom_count": len(output["entry_receptor_source_atoms"]),
            "estimated_coordinate_transform_workload": transform,
            "estimated_pair_search_workload": pair_work, "exception_count": len(exceptions),
        }
        return output, stats, exceptions
    except Exception as exc:
        stats = {"pdb_id": pid, "bucket_id": bid, "mmcif_path": path, "has_active_work": True,
                 "active_assembly_count": 0, "active_model_count": 0, "mapped_source_ligand_count": 0,
                 "ligand_placement_count": 0, "active_receptor_source_chain_count": 0,
                 "active_receptor_chain_instance_count": 0, "ligand_source_atom_count": 0,
                 "receptor_source_atom_count": 0, "estimated_coordinate_transform_workload": 0,
                 "estimated_pair_search_workload": 0, "exception_count": 1}
        return output, stats, [("mmcif_parse_failure", f"{type(exc).__name__}:{exc}")]


def process_task(task):
    task_id, bid, pids = task
    marker = BUILD / "work/checkpoints" / f"task-{task_id:06d}.json"
    if marker.exists():
        return json.loads(marker.read_text())
    combined = {name: [] for name in DATASETS}
    stats = []
    exceptions = []
    for pid in pids:
        rows, stat, errors = parse_entry(pid)
        stats.append(stat)
        exceptions.extend({"pdb_id": pid, "error_type": kind, "error_message": detail} for kind, detail in errors)
        for name in DATASETS:
            combined[name].extend(rows[name])
    row_counts = {}
    for name in DATASETS:
        path = BUILD / "work/bucket_fragments" / name / f"bucket_id={bid:03d}" / f"part-{task_id:06d}.parquet"
        row_counts[name] = write_parquet_atomic(path, combined[name], SCHEMAS[name])
    stats_path = BUILD / "work/checkpoints" / f"task-{task_id:06d}-stats.json"
    stats_path.write_text(json.dumps({"stats": stats, "exceptions": exceptions}, separators=(",", ":")) + "\n")
    payload = {"task_id": task_id, "bucket_id": bid, "pdb_ids": pids, "row_counts": row_counts,
               "exception_count": len(exceptions), "completed_at": utc()}
    marker.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def init_queue(active_pdb):
    db = BUILD / "work/work_queue.sqlite"
    with sqlite3.connect(db) as conn:
        conn.execute("CREATE TABLE IF NOT EXISTS entry_status (pdb_id TEXT PRIMARY KEY,status TEXT,attempt_count INTEGER,worker_id TEXT,started_at TEXT,updated_at TEXT,finished_at TEXT,error_type TEXT,error_message TEXT)")
        conn.executemany("INSERT OR IGNORE INTO entry_status VALUES (?,?,?,?,?,?,?,?,?)",
                         [(pid, "pending" if pid in active_pdb else "skipped", 0, "", "", utc(), "", "", "") for pid in G_ENTRY_IDS])
    return db


def run_build():
    active_keys = load_inputs()
    active_pdb = sorted(G_PLACEMENTS)
    db = init_queue(set(active_pdb))
    grouped = defaultdict(list)
    for pid in active_pdb:
        grouped[bucket_id(pid, int(G_CFG["bucket_count"]))].append(pid)
    tasks = []
    task_id = 0
    size = int(G_CFG["rows_per_task"])
    for bid in range(int(G_CFG["bucket_count"])):
        for start in range(0, len(grouped[bid]), size):
            tasks.append((task_id, bid, grouped[bid][start:start + size])); task_id += 1
    pending = [task for task in tasks if not (BUILD / "work/checkpoints" / f"task-{task[0]:06d}.json").exists()]
    started = time.time()
    progress = {"status": "RUNNING", "active_pdb_count": len(active_pdb), "active_assembly_key_count": len(active_keys),
                "task_count": len(tasks), "completed_task_count": len(tasks) - len(pending), "started_at": utc()}
    (BUILD / "work/progress.json").write_text(json.dumps(progress, indent=2) + "\n")
    def record_result(result):
        progress["completed_task_count"] += 1
        progress["last_task_id"] = result["task_id"]
        progress["updated_at"] = utc()
        if progress["completed_task_count"] % 20 == 0:
            (BUILD / "work/progress.json").write_text(json.dumps(progress, indent=2) + "\n")
        with sqlite3.connect(db) as conn:
            status = "failed" if result["exception_count"] else "success"
            conn.executemany(
                "UPDATE entry_status SET status=?,attempt_count=attempt_count+1,updated_at=?,finished_at=? WHERE pdb_id=?",
                [(status, utc(), utc(), pid) for pid in result["pdb_ids"]],
            )

    # Large compressed entries and long-lived parser workers can retain substantial
    # memory. Isolate large tasks and recycle regular workers without changing any
    # record-level parsing or output semantics.
    large_limit = int(G_CFG.get("large_mmcif_bytes", 16 * 1024 * 1024))
    large_pending, regular_pending = [], []
    for task in pending:
        target = large_pending if any(
            Path(G_P1[pid]["mmcif_path"]).stat().st_size >= large_limit for pid in task[2]
        ) else regular_pending
        target.append(task)
    progress["large_pending_task_count"] = len(large_pending)
    progress["regular_pending_task_count"] = len(regular_pending)
    (BUILD / "work/progress.json").write_text(json.dumps(progress, indent=2) + "\n")

    context = multiprocessing.get_context("fork")
    for task in large_pending:
        with context.Pool(processes=1, maxtasksperchild=1) as pool:
            record_result(pool.apply(process_task, (task,)))

    with context.Pool(
        processes=int(G_CFG["workers"]),
        maxtasksperchild=int(G_CFG.get("max_tasks_per_worker", 4)),
    ) as pool:
        for result in pool.imap_unordered(process_task, regular_pending, chunksize=1):
            record_result(result)
    progress.update({"status": "COMPLETED", "completed_task_count": len(tasks), "completed_at": utc(), "runtime_seconds": round(time.time() - started, 2)})
    (BUILD / "work/progress.json").write_text(json.dumps(progress, indent=2) + "\n")
    print(json.dumps(progress, indent=2))


def manifest_schema():
    s = pa.string(); i64 = pa.int64(); i16 = pa.int16(); b = pa.bool_()
    return pa.schema([
        ("pdb_id", s), ("bucket_id", i16), ("mmcif_path", s), ("has_active_work", b),
        ("active_assembly_count", i64), ("active_model_count", i64), ("mapped_source_ligand_count", i64),
        ("ligand_placement_count", i64), ("active_receptor_source_chain_count", i64),
        ("active_receptor_chain_instance_count", i64), ("ligand_source_atom_count", i64),
        ("receptor_source_atom_count", i64), ("ligand_partition_path", s), ("receptor_partition_path", s),
        ("assembly_context_partition_path", s), ("ligand_atom_partition_path", s),
        ("receptor_atom_partition_path", s), ("estimated_coordinate_transform_workload", i64),
        ("estimated_pair_search_workload", i64), ("workload_class", s),
    ])


def hardlink_fragments():
    temp = BUILD / "tmp/output_finalizing"
    if temp.exists():
        shutil.rmtree(temp)
    temp.mkdir(parents=True)
    for name in DATASETS:
        source = BUILD / "work/bucket_fragments" / name
        target = temp / name
        for path in source.rglob("*.parquet"):
            destination = target / path.relative_to(source)
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.link(path, destination)
    return temp


def finalize():
    global G_CFG, G_P1, G_ENTRY_IDS
    G_CFG = yaml.safe_load((ROOT / "configs/build_v1.yaml").read_text())
    G_P1 = {r["pdb_id"]: r for r in iter_tsv(P1)}
    G_ENTRY_IDS = [r["pdb_id"] for r in iter_tsv(P_ENTRIES)]
    progress = json.loads((BUILD / "work/progress.json").read_text())
    if progress.get("status") != "COMPLETED":
        raise SystemExit("build is not complete")
    stats_by_pid = {}
    exceptions = []
    for path in sorted((BUILD / "work/checkpoints").glob("task-*-stats.json")):
        payload = json.loads(path.read_text())
        for row in payload["stats"]:
            stats_by_pid[row["pdb_id"]] = row
        exceptions.extend(payload["exceptions"])
    workloads = sorted(row["estimated_pair_search_workload"] for row in stats_by_pid.values())
    def percentile(p):
        return workloads[min(len(workloads) - 1, int((len(workloads) - 1) * p))] if workloads else 0
    cuts = {"p50": percentile(.50), "p90": percentile(.90), "p99": percentile(.99)}
    manifest = []
    for pid in sorted(G_ENTRY_IDS):
        bid = bucket_id(pid, int(G_CFG["bucket_count"]))
        row = stats_by_pid.get(pid, {"pdb_id": pid, "bucket_id": bid, "mmcif_path": G_P1[pid]["mmcif_path"], "has_active_work": False,
             "active_assembly_count": 0, "active_model_count": 0, "mapped_source_ligand_count": 0, "ligand_placement_count": 0,
             "active_receptor_source_chain_count": 0, "active_receptor_chain_instance_count": 0, "ligand_source_atom_count": 0,
             "receptor_source_atom_count": 0, "estimated_coordinate_transform_workload": 0, "estimated_pair_search_workload": 0})
        workload = row["estimated_pair_search_workload"]
        klass = "none" if not row["has_active_work"] else ("small" if workload <= cuts["p50"] else "medium" if workload <= cuts["p90"] else "large" if workload <= cuts["p99"] else "very_large")
        prefix = f"bucket_id={bid:03d}"
        manifest.append({**{k: row[k] for k in ["pdb_id", "bucket_id", "mmcif_path", "has_active_work", "active_assembly_count", "active_model_count", "mapped_source_ligand_count", "ligand_placement_count", "active_receptor_source_chain_count", "active_receptor_chain_instance_count", "ligand_source_atom_count", "receptor_source_atom_count"]},
                         "ligand_partition_path": f"entry_ligand_placements/{prefix}", "receptor_partition_path": f"entry_receptor_chain_instances/{prefix}",
                         "assembly_context_partition_path": f"entry_assembly_context/{prefix}", "ligand_atom_partition_path": f"entry_ligand_source_atoms/{prefix}",
                         "receptor_atom_partition_path": f"entry_receptor_source_atoms/{prefix}",
                         "estimated_coordinate_transform_workload": row["estimated_coordinate_transform_workload"],
                         "estimated_pair_search_workload": workload, "workload_class": klass})
    temp = hardlink_fragments()
    write_parquet_atomic(temp / "entry_work_manifest.parquet", manifest, manifest_schema())
    schema_payload = {
        "schema_version": G_CFG["schema_version"],
        "datasets": {
            **{name: [{"name": field.name, "type": str(field.type), "nullable": field.nullable} for field in schema]
               for name, schema in SCHEMAS.items()},
            "entry_work_manifest": [{"name": field.name, "type": str(field.type), "nullable": field.nullable}
                                    for field in manifest_schema()],
        },
    }
    (temp / "dataset_schemas.json").write_text(json.dumps(schema_payload, indent=2) + "\n")
    if exceptions:
        exception_schema = pa.schema([("pdb_id", pa.string()), ("error_type", pa.string()), ("error_message", pa.string())])
        write_parquet_atomic(temp / "entry_package_exceptions.parquet", sorted(exceptions, key=lambda x: (x["pdb_id"], x["error_type"])), exception_schema)
    output = BUILD / "output"
    backup = BUILD / "tmp/output_previous"
    if backup.exists(): shutil.rmtree(backup)
    if output.exists() and any(output.iterdir()): output.rename(backup)
    temp.rename(output)
    (BUILD / "audit/workload_boundaries.json").write_text(json.dumps(cuts, indent=2) + "\n")
    print(json.dumps({"manifest_rows": len(manifest), "active_rows": len(stats_by_pid), "exceptions": len(exceptions), "workload_boundaries": cuts}, indent=2))


def dataset_count(path):
    return ds.dataset(path, format="parquet").count_rows()


def validate_and_complete():
    global G_CFG
    G_CFG = yaml.safe_load((ROOT / "configs/build_v1.yaml").read_text())
    output = BUILD / "output"
    manifest_table = pq.read_table(output / "entry_work_manifest.parquet")
    manifest_rows = manifest_table.to_pylist()
    active = [r for r in manifest_rows if r["has_active_work"]]
    counts = {name: dataset_count(output / name) for name in DATASETS}
    placement_ids = ds.dataset(output / "entry_ligand_placements", format="parquet").to_table(columns=["filter_2_ligand_assembly_placement_id", "filter_2_source_ligand_instance_id", "pdb_id", "assembly_id", "model_id"])
    pids = placement_ids.column(0).to_pylist(); source_ids = placement_ids.column(1).to_pylist()
    receptor_ids = ds.dataset(output / "entry_receptor_chain_instances", format="parquet").to_table(columns=["filter_1_chain_instance_id", "polymer_class", "receptor_eligible"])
    rid = receptor_ids.column(0).to_pylist()
    active_keys = set(zip(placement_ids.column(2).to_pylist(), placement_ids.column(3).to_pylist(), placement_ids.column(4).to_pylist()))
    expected_receptor_ids = {
        row["chain_instance_id"] for row in iter_tsv(P_F1_RECEPTORS)
        if (row["pdb_id"], row["assembly_id"], row["model_id"]) in active_keys
    }
    actual_receptor_ids = set(rid)
    exception_path = output / "entry_package_exceptions.parquet"
    exception_count = pq.read_metadata(exception_path).num_rows if exception_path.exists() else 0
    validation = {
        "timestamp": utc(), "manifest_rows": len(manifest_rows), "active_pdb_count": len(active),
        "dataset_row_counts": counts, "ligand_placement_count_expected": 1151324,
        "ligand_placement_missing_or_extra": abs(counts["entry_ligand_placements"] - 1151324),
        "duplicate_ligand_placement_id": len(pids) - len(set(pids)),
        "mapped_source_ligand_count": len(set(source_ids)), "mapped_source_ligand_expected": 851966,
        "excluded_no_mapping_source_count": 1002,
        "duplicate_receptor_chain_instance_id": len(rid) - len(set(rid)),
        "missing_active_receptor_chain_instance": len(expected_receptor_ids - actual_receptor_ids),
        "extra_receptor_chain_instance": len(actual_receptor_ids - expected_receptor_ids),
        "nonprotein_receptor_count": sum(x.as_py() != "POLYPEPTIDE" for x in receptor_ids.column(1)),
        "ineligible_receptor_count": sum(x.as_py() is not True for x in receptor_ids.column(2)),
        "manifest_ligand_placement_sum": sum(r["ligand_placement_count"] for r in manifest_rows),
        "manifest_ligand_atom_sum": sum(r["ligand_source_atom_count"] for r in manifest_rows),
        "manifest_receptor_atom_sum": sum(r["receptor_source_atom_count"] for r in manifest_rows),
        "exception_count": exception_count, "exception_file_present": exception_path.exists(),
        "transformed_coordinates_generated": False, "protein_ligand_pairs_generated": False,
        "distance_or_contact_tables_generated": False, "ccd_atom_mapping_generated": False,
        "topology_restoration_generated": False,
    }
    validation["validation_pass"] = all([
        len(manifest_rows) == 248037, counts["entry_ligand_placements"] == 1151324,
        validation["duplicate_ligand_placement_id"] == 0, validation["mapped_source_ligand_count"] == 851966,
        validation["duplicate_receptor_chain_instance_id"] == 0,
        validation["missing_active_receptor_chain_instance"] == 0,
        validation["extra_receptor_chain_instance"] == 0,
        validation["nonprotein_receptor_count"] == 0, validation["ineligible_receptor_count"] == 0,
        validation["manifest_ligand_placement_sum"] == 1151324,
        validation["manifest_ligand_atom_sum"] == counts["entry_ligand_source_atoms"],
        validation["manifest_receptor_atom_sum"] == counts["entry_receptor_source_atoms"],
        exception_count == 0,
    ])
    (output / "validation_report.json").write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["validation_pass"]:
        raise SystemExit(json.dumps(validation, indent=2))
    parts = sorted(output.rglob("*.parquet"))
    total_size = sum(p.stat().st_size for p in parts)
    summary = {
        "build_id": BUILD_ID, "input_entries": 248037, "active_pdb_count": len(active),
        "active_assembly_key_count": len(active_keys),
        "ligand_placements": counts["entry_ligand_placements"], "mapped_source_ligands": len(set(source_ids)),
        "excluded_no_mapping_source_ligands": 1002, "active_receptor_chain_instances": counts["entry_receptor_chain_instances"],
        "active_receptor_source_chains": sum(r["active_receptor_source_chain_count"] for r in active),
        "ligand_source_atoms": counts["entry_ligand_source_atoms"], "receptor_source_atoms": counts["entry_receptor_source_atoms"],
        "assembly_context_rows": counts["entry_assembly_context"], "bucket_count": int(G_CFG["bucket_count"]),
        "parquet_part_files": len(parts), "total_size_bytes": total_size, "exception_count": 0,
        "validation_pass": True,
    }
    metadata = {
        "build_id": BUILD_ID, "created_at": utc(), "python": sys.version, "pyarrow": pa.__version__, "gemmi": gemmi.__version__,
        "platform": platform.platform(), "compression": G_CFG["compression"], "parquet_version": G_CFG["parquet_version"],
        "writer_options": G_CFG, "code_sha256": sha256(ROOT / "scripts/build_entry_work_packages.py"),
    }
    (output / "build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "build_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    manifest = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name not in {"output_manifest.tsv", "SHA256SUMS", "build_complete.json"}:
            manifest.append({"relative_path": str(path.relative_to(output)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    with (output / "output_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest)
    checksum_paths = [output / r["relative_path"] for r in manifest] + [output / "output_manifest.tsv"]
    (output / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(output)}\n" for p in checksum_paths))
    complete = {"status": "COMPLETE", "build_id": BUILD_ID, "completed_at": utc(), "validation_pass": True,
                "output_manifest_sha256": sha256(output / "output_manifest.tsv"), "scientific_freeze": False}
    (BUILD / "build_complete.json").write_text(json.dumps(complete, indent=2) + "\n")
    (ROOT / "build_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (ROOT / "CURRENT_BUILD.json").write_text(json.dumps({"current_build_id": BUILD_ID, "status": "COMPLETE", "scientific_freeze": False,
                                                          "relative_path": f"builds/{BUILD_ID}", "updated_at": utc()}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def selftest():
    assert bucket_id("1abc") == bucket_id("1ABC")
    assert expand_operator_paths("(1-2)(3,4)") == [("1", "3"), ("1", "4"), ("2", "3"), ("2", "4")]
    ops = {"1": identity_affine()}
    assert composite_affine(("1",), ops) == identity_affine()
    assert parse_asym_id_list(";A,B,C\n;") == ["A", "B", "C"]
    assert expand_operator_paths(";1-2\n;") == [("1",), ("2",)]
    assert len(SCHEMAS) == 5
    print(json.dumps({"selftest_pass": True, "tests": 6}))


def integration_test():
    load_inputs()
    pids = sorted(G_PLACEMENTS)[:17] + ["8at5", "9f59", "9v10"]
    totals = Counter()
    errors = []
    for pid in pids:
        rows, stats, entry_errors = parse_entry(pid)
        totals.update({name: len(value) for name, value in rows.items()})
        errors.extend({"pdb_id": pid, "error_type": kind, "detail": detail} for kind, detail in entry_errors)
        if stats["ligand_placement_count"] != len(G_PLACEMENTS[pid]):
            errors.append({"pdb_id": pid, "error_type": "placement_count_mismatch", "detail": str(stats)})
    report = {"tested_pdb_ids": pids, "dataset_row_counts": dict(totals), "errors": errors,
              "integration_test_pass": not errors}
    (BUILD / "audit/integration_test_20.json").write_text(json.dumps(report, indent=2) + "\n")
    if errors:
        raise SystemExit(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["setup", "preflight", "selftest", "integration-test", "build", "finalize", "validate"])
    args = parser.parse_args()
    {"setup": setup, "preflight": preflight, "selftest": selftest, "integration-test": integration_test, "build": run_build,
     "finalize": finalize, "validate": validate_and_complete}[args.command]()


if __name__ == "__main__":
    main()
