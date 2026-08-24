#!/usr/bin/env python3
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
import shutil
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import yaml

ROOT = Path("/root/autodl-tmp/benchmark_1.0")
V2 = ROOT / "filter_2_ligand_qualification_v2"
V2_RUN = V2 / "runs/20260804_full_01"
OUT = ROOT / "filter_2_ligand_qualification_v3"
RUN_ID = "20260804_full_01"
RUN = OUT / "runs" / RUN_ID
RULE_VERSION = "filter_2_v3.0.0_policy_revision"
SCHEMA_VERSION = "3.0.0"
METAL_POLICY_VERSION = "metal_elements_v1"
ACCEPTED_CCD = {"ccd_exact", "ccd_obsolete_resolved"}

SOURCE_INPUT = V2_RUN / "output/filter_2_source_instances.tsv.gz"
WATER_INPUT = V2_RUN / "output/filter_2_water_exclusion_summary.tsv"
ENTRY_INPUT = V2_RUN / "output/filter_2_entries.tsv.gz"
P1_INDEX = V2 / "inputs/processing_1_mmcif_index_snapshot.tsv.gz"
F1_ENTRIES = V2 / "inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz"
F1_ASSEMBLIES = V2 / "inputs/filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"
F1_RECEPTORS = V2 / "inputs/filter_1_receptor_chain_instances_snapshot.tsv.gz"
BIOLIP_INPUT = V2 / "references/BioLiP_ligand_list_snapshot.tsv"
BIOLIP_META = V2 / "references/BioLiP_ligand_list_metadata.json"
CCD_META = V2 / "references/ccd_snapshot_validation.json"

EXPECTED_HASHES = {
    "source_inventory": "618129f4afbfb902d1887defbe46f48d22e7bf970c2e4486486f6c347481fa55",
    "water_summary": "f5af11d656673a9a1cb03988fc024ed6c920315d7185b6d6404799741dbeb294",
    "entry_inventory": "5e087e107e63786bfea2a2c319bd4e3f5d3b188acb33b0efb71e431600778a48",
    "p1_index": "531b7c42219231dae986a866fde86b12a7548363130b2bb279e5afa0ae6317a1",
    "f1_entries": "11b714fabf38322b057a0b84bbfb9833a0c09e17fe994cde2e8005c8612dea59",
    "f1_assemblies": "9d7d37dd345529a4096e9e7460942e491cf205f9269556024fdfc3bd3970986e",
    "f1_receptors": "3fdb3865f1c148bafad5f944dbf449b0eb7e308274faf4b21a97a59fce0bdddb",
    "biolip": "e9b53b2b58657fa8fee3041560e841e85afcbba1f68879a0800799bfa0c7067b",
}

ROUTE_FIELDS = [
    "pdb_id", "selected_model_id", "entity_id", "component_id", "label_comp_id", "auth_comp_id",
    "label_asym_id", "auth_asym_id", "auth_seq_id", "pdb_seq_num", "insertion_code",
    "source_ligand_instance_id", "atom_count", "observed_heavy_atom_count", "observed_element_composition",
    "altloc_values", "occupancy_min", "occupancy_max", "entity_type", "source_instance_status",
    "source_instance_count_in_entry", "biolip_list_match", "biolip_snapshot_sha256",
    "original_component_id", "resolved_ccd_id", "ccd_identity_status", "ccd_name", "ccd_type", "formula",
    "formula_weight", "formal_charge", "parent_component_id", "expected_atom_count", "standard_total_atom_count",
    "expected_heavy_atom_count", "element_set", "carbon_atom_count", "fragment_count", "contains_metal",
    "descriptor_availability", "ccd_snapshot_version", "ccd_snapshot_sha256", "previous_terminal_route",
    "previous_reason_code", "terminal_route", "decision", "destination", "reason_code", "reason_detail",
    "matched_metal_elements", "metal_policy_version", "metal_elements_sha256", "rule_version",
]
PLACEMENT_FIELDS = [
    "pdb_id", "assembly_id", "selected_model_id", "source_ligand_instance_id", "component_id", "label_asym_id",
    "assembly_gen_row_id", "oper_expression_raw", "operator_path", "composite_operator_id", "rotation_matrix",
    "translation_vector", "assembly_ligand_placement_id", "mapping_status", "rule_version",
]
NO_MAP_FIELDS = [
    "pdb_id", "selected_model_id", "source_ligand_instance_id", "component_id", "label_asym_id",
    "mapping_status", "reason_code", "reason_detail", "compatible_assembly_count", "referenced_assembly_gen_rows",
    "oper_expression_parse_failure_count", "operator_id_missing_count", "operator_matrix_missing_count", "rule_version",
]
DELTA_FIELDS = ["source_ligand_instance_id", "pdb_id", "component_id", "old_route", "new_route", "changed", "rule_version"]
ENTRY_FIELDS = [
    "pdb_id", "input_source_instance_count", "excluded_atomic_component_count",
    "excluded_metal_containing_component_count", "excluded_biolip_artifact_count", "ccd_review_count",
    "provisional_source_ligand_count", "logical_placement_count", "no_mapping_count", "entry_status", "rule_version",
]
TABLES = {
    "source_terminal_routes": ROUTE_FIELDS,
    "source_exclusions": ROUTE_FIELDS,
    "excluded_atomic_components": ROUTE_FIELDS,
    "excluded_metal_containing_components": ROUTE_FIELDS,
    "excluded_biolip_artifacts": ROUTE_FIELDS,
    "ccd_review": ROUTE_FIELDS,
    "provisional_source_ligands": ROUTE_FIELDS,
    "ligand_assembly_logical_placements": PLACEMENT_FIELDS,
    "ligand_assembly_no_mapping": NO_MAP_FIELDS,
    "entry_route_summary": ENTRY_FIELDS,
    "route_delta_vs_previous": DELTA_FIELDS,
}

G_METALS = set()
G_BIOLIP = set()
G_METAL_SHA = ""


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if str(path).endswith(".gz") else Path(path).open(mode.replace("t", ""), encoding="utf-8", newline="")


def iter_tsv(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, rows, fields):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value in {".", "?", "None"} else value


def stable_id(*parts):
    return "|".join(clean(x) or "~" for x in parts)


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
        import re
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
    import re
    compact = expression.replace(" ", "")
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
        if operator_id not in operators:
            raise KeyError(operator_id)
        total = compose_affine(total, operators[operator_id])
    return total


def format_vector(values):
    return ",".join(f"{value:.10g}" for value in values)


def format_matrix(matrix):
    return ";".join(format_vector(row) for row in matrix)


def parse_operator_records(block):
    all_ids, valid = set(), {}
    for row in category_records(block, "_pdbx_struct_oper_list."):
        oid = row.get("id", "")
        if not oid:
            continue
        all_ids.add(oid)
        try:
            matrix = [[float(row[f"matrix[{i}][{j}]"]) for j in range(1, 4)] for i in range(1, 4)]
            vector = [float(row[f"vector[{i}]"]) for i in range(1, 4)]
            if all(math.isfinite(x) for line in matrix for x in line) and all(math.isfinite(x) for x in vector):
                valid[oid] = (matrix, vector)
        except Exception:
            pass
    return all_ids, valid


def int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def load_policy(path):
    payload = yaml.safe_load(Path(path).read_text())
    return {x.upper() for values in payload["groups"].values() for x in values}


def init_worker(metal_path, biolip_path):
    global G_METALS, G_BIOLIP, G_METAL_SHA
    G_METALS = load_policy(metal_path)
    G_METAL_SHA = sha256(metal_path)
    with Path(biolip_path).open(encoding="utf-8") as handle:
        G_BIOLIP = {line.split("\t", 1)[0].strip().upper() for line in handle if line.strip()}


def route_source(old):
    status = old.get("ccd_identity_status", "")
    total = int_or_none(old.get("expected_atom_count", ""))
    if status not in ACCEPTED_CCD:
        return "ccd_review", "REVIEW", "ccd_review", "ccd_identity_not_accepted", []
    if total is None or total < 1:
        return "ccd_review", "REVIEW", "ccd_review", "ccd_standard_atom_definition_missing", []
    if total == 1:
        return "excluded_atomic_component", "REJECT", "excluded_atomic_component", "single_atom_component", []
    elements = {x for x in old.get("element_set", "").split(",") if x}
    matched = sorted(elements & G_METALS)
    if matched:
        return "excluded_metal_containing_component", "REJECT", "metal_containing_special_subset", "ccd_component_contains_metal_element", matched
    membership = old.get("component_id", "").upper() in G_BIOLIP
    count = int(old.get("source_instance_count_in_entry", "") or 0)
    if membership and count >= 15:
        return "excluded_biolip_artifact", "REJECT", "excluded_biolip_artifact", "frozen_biolip_list_and_same_pdb_component_count_ge_15", []
    return "provisional_source_ligand", "PASS", "provisional_source_ligand", "accepted_multiatom_nonmetal_not_high_frequency_biolip", []


def transform_row(old):
    route, decision, destination, reason, matched = route_source(old)
    row = {field: old.get(field, "") for field in ROUTE_FIELDS}
    row.update({
        "standard_total_atom_count": old.get("expected_atom_count", ""),
        "previous_terminal_route": old.get("terminal_route", ""),
        "previous_reason_code": old.get("reason_code", ""),
        "terminal_route": route, "decision": decision, "destination": destination,
        "reason_code": reason, "reason_detail": "", "matched_metal_elements": ",".join(matched),
        "metal_policy_version": METAL_POLICY_VERSION, "metal_elements_sha256": G_METAL_SHA,
        "rule_version": RULE_VERSION,
    })
    return row


def build_assembly_evidence(block, retained):
    all_ids, operators = parse_operator_records(block)
    reverse = defaultdict(list)
    for index, gen in enumerate(category_records(block, "_pdbx_struct_assembly_gen."), start=1):
        aid = gen.get("assembly_id", "")
        if aid not in retained:
            continue
        row_id = f"assembly_gen_row_{index:06d}"
        expression = gen.get("oper_expression", "")
        asyms = [x.strip() for x in gen.get("asym_id_list", "").split(",") if x.strip()]
        try:
            paths, parse_error = expand_operator_paths(expression), ""
        except Exception as exc:
            paths, parse_error = [], f"{type(exc).__name__}:{exc}"
        for asym in asyms:
            reverse[(aid, asym)].append((row_id, expression, paths, parse_error))
    return all_ids, operators, reverse


def process_item(item):
    pid, old_rows = item["pdb_id"], item["source_rows"]
    result = {name: [] for name in TABLES}
    transformed = [transform_row(old) for old in old_rows]
    result["source_terminal_routes"] = transformed
    route_to_table = {
        "excluded_atomic_component": "excluded_atomic_components",
        "excluded_metal_containing_component": "excluded_metal_containing_components",
        "excluded_biolip_artifact": "excluded_biolip_artifacts",
        "ccd_review": "ccd_review", "provisional_source_ligand": "provisional_source_ligands",
    }
    for row in transformed:
        table = route_to_table[row["terminal_route"]]
        result[table].append(row)
        if row["terminal_route"].startswith("excluded_"):
            result["source_exclusions"].append(row)
        result["route_delta_vs_previous"].append({
            "source_ligand_instance_id": row["source_ligand_instance_id"], "pdb_id": pid,
            "component_id": row["component_id"], "old_route": row["previous_terminal_route"],
            "new_route": row["terminal_route"], "changed": str(row["previous_terminal_route"] != row["terminal_route"]).lower(),
            "rule_version": RULE_VERSION,
        })
    provisional = result["provisional_source_ligands"]
    if provisional:
        try:
            block = gemmi.cif.read(item["mmcif_path"]).sole_block()
            retained = set(item["retained_assembly_ids"])
            all_operator_ids, operators, reverse = build_assembly_evidence(block, retained)
            mapped_source_ids = set()
            for row in provisional:
                compatible = [aid for aid in sorted(retained) if not item["assembly_models"].get(aid) or row["selected_model_id"] in item["assembly_models"][aid]]
                counts = Counter()
                references = []
                for aid in compatible:
                    references.extend((aid, *ref) for ref in reverse.get((aid, row["label_asym_id"]), []))
                for aid, row_id, expression, paths, parse_error in references:
                    if parse_error:
                        counts["oper_expression_parse_failure"] += 1
                        continue
                    for path in paths:
                        missing = [op for op in path if op not in all_operator_ids]
                        invalid = [op for op in path if op in all_operator_ids and op not in operators]
                        if missing:
                            counts["operator_id_missing"] += 1
                            continue
                        if invalid:
                            counts["operator_matrix_missing"] += 1
                            continue
                        try:
                            rotation, translation = composite_affine(path, operators)
                        except Exception:
                            counts["other_mapping_failure"] += 1
                            continue
                        operator_path = "*".join(path)
                        placement_id = stable_id(pid, aid, row["selected_model_id"], row["source_ligand_instance_id"], row_id, operator_path)
                        result["ligand_assembly_logical_placements"].append({
                            "pdb_id": pid, "assembly_id": aid, "selected_model_id": row["selected_model_id"],
                            "source_ligand_instance_id": row["source_ligand_instance_id"], "component_id": row["component_id"],
                            "label_asym_id": row["label_asym_id"], "assembly_gen_row_id": row_id,
                            "oper_expression_raw": expression, "operator_path": operator_path,
                            "composite_operator_id": operator_path, "rotation_matrix": format_matrix(rotation),
                            "translation_vector": format_vector(translation), "assembly_ligand_placement_id": placement_id,
                            "mapping_status": "mapped", "rule_version": RULE_VERSION,
                        })
                        mapped_source_ids.add(row["source_ligand_instance_id"])
                mapped = row["source_ligand_instance_id"] in mapped_source_ids
                if not mapped:
                    if not compatible:
                        reason = "assembly_model_mismatch"
                    elif not references:
                        reason = "asym_not_referenced_by_retained_assembly"
                    elif counts["oper_expression_parse_failure"]:
                        reason = "oper_expression_parse_failure"
                    elif counts["operator_id_missing"]:
                        reason = "operator_id_missing"
                    elif counts["operator_matrix_missing"]:
                        reason = "operator_matrix_missing"
                    else:
                        reason = "other_mapping_failure"
                    result["ligand_assembly_no_mapping"].append({
                        "pdb_id": pid, "selected_model_id": row["selected_model_id"],
                        "source_ligand_instance_id": row["source_ligand_instance_id"], "component_id": row["component_id"],
                        "label_asym_id": row["label_asym_id"], "mapping_status": "no_mapping", "reason_code": reason,
                        "reason_detail": json.dumps(dict(counts), sort_keys=True), "compatible_assembly_count": str(len(compatible)),
                        "referenced_assembly_gen_rows": str(len(references)),
                        "oper_expression_parse_failure_count": str(counts["oper_expression_parse_failure"]),
                        "operator_id_missing_count": str(counts["operator_id_missing"]),
                        "operator_matrix_missing_count": str(counts["operator_matrix_missing"]), "rule_version": RULE_VERSION,
                    })
        except Exception as exc:
            for row in provisional:
                result["ligand_assembly_no_mapping"].append({
                    "pdb_id": pid, "selected_model_id": row["selected_model_id"], "source_ligand_instance_id": row["source_ligand_instance_id"],
                    "component_id": row["component_id"], "label_asym_id": row["label_asym_id"], "mapping_status": "no_mapping",
                    "reason_code": "other_mapping_failure", "reason_detail": f"{type(exc).__name__}:{exc}"[:1000],
                    "compatible_assembly_count": "0", "referenced_assembly_gen_rows": "0",
                    "oper_expression_parse_failure_count": "0", "operator_id_missing_count": "0",
                    "operator_matrix_missing_count": "0", "rule_version": RULE_VERSION,
                })
    counts = Counter(row["terminal_route"] for row in transformed)
    result["entry_route_summary"] = [{
        "pdb_id": pid, "input_source_instance_count": str(len(transformed)),
        "excluded_atomic_component_count": str(counts["excluded_atomic_component"]),
        "excluded_metal_containing_component_count": str(counts["excluded_metal_containing_component"]),
        "excluded_biolip_artifact_count": str(counts["excluded_biolip_artifact"]),
        "ccd_review_count": str(counts["ccd_review"]), "provisional_source_ligand_count": str(counts["provisional_source_ligand"]),
        "logical_placement_count": str(len(result["ligand_assembly_logical_placements"])),
        "no_mapping_count": str(len(result["ligand_assembly_no_mapping"])), "entry_status": "success", "rule_version": RULE_VERSION,
    }]
    return result


def load_items():
    index = {r["pdb_id"]: r["mmcif_path"] for r in iter_tsv(P1_INDEX)}
    entries = [r["pdb_id"] for r in iter_tsv(F1_ENTRIES)]
    assemblies = defaultdict(set)
    for r in iter_tsv(F1_ASSEMBLIES): assemblies[r["pdb_id"]].add(r["assembly_id"])
    models = defaultdict(lambda: defaultdict(set))
    for r in iter_tsv(F1_RECEPTORS): models[r["pdb_id"]][r["assembly_id"]].add(r["model_id"])
    return [{"pdb_id": pid, "mmcif_path": index[pid], "retained_assembly_ids": sorted(assemblies[pid]),
             "assembly_models": {aid: sorted(vals) for aid, vals in models[pid].items()}} for pid in entries]


def source_group_iterator():
    current, rows = None, []
    for row in iter_tsv(SOURCE_INPUT):
        pid = row["pdb_id"]
        if current is not None and pid != current:
            yield current, rows
            rows = []
        current = pid
        rows.append(row)
    if current is not None:
        yield current, rows


def setup():
    if OUT.exists():
        raise SystemExit(f"target already exists: {OUT}")
    for path in [OUT / "scripts", OUT / "config", OUT / "references", OUT / "schemas", OUT / "tests", RUN / "input", RUN / "work/batches", RUN / "output", RUN / "audit", RUN / "logs", RUN / "release"]:
        path.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__), OUT / "scripts/filter2_v3_pipeline.py")
    local_metal = Path(__file__).with_name("metal_elements_v1.yaml")
    shutil.copy2(local_metal, OUT / "config/metal_elements_v1.yaml")
    shutil.copy2(BIOLIP_INPUT, OUT / "references/BioLiP_ligand_list_snapshot.tsv")
    shutil.copy2(BIOLIP_META, OUT / "references/BioLiP_ligand_list_metadata.json")
    shutil.copy2(CCD_META, OUT / "references/CCD_snapshot_metadata.json")
    cfg = {
        "stage": {"name": "filter_2_ligand_qualification_v3", "version": "3.0.0"},
        "rule_version": RULE_VERSION, "schema_version": SCHEMA_VERSION,
        "accepted_ccd_statuses": sorted(ACCEPTED_CCD), "atomic_rule": "standard_total_atom_count == 1",
        "metal_policy": "config/metal_elements_v1.yaml",
        "biolip": {"same_pdb_component_source_instance_threshold": 15, "operator": ">="},
        "runtime": {"workers": 24, "batch_size": 200, "resume": True},
        "format": {"large_tables": "tsv.gz", "parquet_unavailable_in_frozen_environment": True, "preview_rows": 1000},
    }
    (OUT / "config/filter_2_v3.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    (RUN / "config_snapshot.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    for table, fields in TABLES.items():
        (OUT / "schemas" / f"{table}_schema.json").write_text(json.dumps({"schema_version": SCHEMA_VERSION, "table": table, "format": "TSV.GZ", "fields": fields}, indent=2) + "\n")
    print(OUT)


def preflight():
    if json.loads((V2_RUN / "_FROZEN.json").read_text()).get("status") != "FROZEN":
        raise SystemExit("v2 is not frozen")
    paths = {"source_inventory": SOURCE_INPUT, "water_summary": WATER_INPUT, "entry_inventory": ENTRY_INPUT,
             "p1_index": P1_INDEX, "f1_entries": F1_ENTRIES, "f1_assemblies": F1_ASSEMBLIES,
             "f1_receptors": F1_RECEPTORS, "biolip": BIOLIP_INPUT}
    actual = {k: sha256(v) for k, v in paths.items()}
    hash_pass = all(actual[k] == EXPECTED_HASHES[k] for k in EXPECTED_HASHES)
    entry_order = [r["pdb_id"] for r in iter_tsv(F1_ENTRIES)]
    order_index = {pid: i for i, pid in enumerate(entry_order)}
    source_count, unique_ids, previous_index, source_order_pass = 0, set(), -1, True
    for row in iter_tsv(SOURCE_INPUT):
        source_count += 1
        sid = row["source_ligand_instance_id"]
        unique_ids.add(sid)
        idx = order_index.get(row["pdb_id"], -1)
        if idx < 0: source_order_pass = False
        if idx < previous_index: source_order_pass = False
        previous_index = idx
    water_count = sum(int(r["source_instance_count"]) for r in iter_tsv(WATER_INPUT))
    metal_path = OUT / "config/metal_elements_v1.yaml"
    metals = load_policy(metal_path)
    ccd_meta = json.loads(CCD_META.read_text())
    ccd_path = Path(ccd_meta["source_path"])
    ccd_actual_sha = sha256(ccd_path)
    audit = {
        "timestamp": utc(), "old_frozen_run": str(V2_RUN), "old_run_status": "FROZEN",
        "input_entry_count": len(entry_order), "input_source_instance_count": source_count,
        "unique_source_instance_id_count": len(unique_ids), "duplicate_source_instance_id_count": source_count - len(unique_ids),
        "water_source_instance_count": water_count, "source_order_compatible_with_filter1_entries": source_order_pass,
        "input_paths": {k: str(v) for k, v in paths.items()}, "actual_sha256": actual, "expected_sha256": EXPECTED_HASHES,
        "input_sha256_pass": hash_pass, "metal_policy_version": METAL_POLICY_VERSION,
        "ccd_snapshot_path": str(ccd_path), "ccd_snapshot_sha256": ccd_actual_sha,
        "ccd_snapshot_expected_sha256": ccd_meta["sha256"], "ccd_snapshot_sha256_match": ccd_actual_sha == ccd_meta["sha256"],
        "metal_elements": sorted(metals), "metal_elements_count": len(metals), "metal_elements_sha256": sha256(metal_path),
        "parquet_unavailable_in_frozen_environment": True,
    }
    audit["preflight_pass"] = all([len(entry_order) == 248037, source_count == 2652746, len(unique_ids) == source_count,
                                    water_count == 529742, source_order_pass, hash_pass, audit["ccd_snapshot_sha256_match"]])
    (RUN / "audit/preflight.json").write_text(json.dumps(audit, indent=2) + "\n")
    (RUN / "input/input_provenance.json").write_text(json.dumps(audit, indent=2) + "\n")
    if not audit["preflight_pass"]: raise SystemExit(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


def selftest():
    init_worker(OUT / "config/metal_elements_v1.yaml", OUT / "references/BioLiP_ligand_list_snapshot.tsv")
    def row(component, total, elements, status="ccd_exact", count=1, member=False):
        return {"component_id": component, "expected_atom_count": str(total), "element_set": ",".join(elements),
                "ccd_identity_status": status, "source_instance_count_in_entry": str(count), "biolip_list_match": str(member).lower()}
    tests = {}
    for comp, element in {"NA":"NA","K":"K","MG":"MG","CA":"CA","ZN":"ZN","CL":"CL","BR":"BR","IOD":"I"}.items():
        tests[f"atomic_{comp}"] = route_source(row(comp, 1, [element]))[0] == "excluded_atomic_component"
    for comp, atoms, elements in [("NH4",5,["H","N"]),("NH3",4,["H","N"]),("OH",2,["H","O"]),("H2S",3,["H","S"]),("NH2",3,["H","N"])]:
        tests[f"not_atomic_{comp}"] = route_source(row(comp, atoms, elements))[0] != "excluded_atomic_component"
    for comp, elements in {"SF4":["FE","S"],"FES":["FE","S"],"F3S":["FE","S"],"VO4":["O","V"],"WO4":["O","W"],"MOO":["MO","O"],"IUM":["O","U"],"CPT":["CL","H","N","PT"],"BEF":["BE","F"],"AF3":["AL","F"],"OEX":["CA","MN","O"]}.items():
        tests[f"metal_{comp}"] = route_source(row(comp, 4, elements))[0] == "excluded_metal_containing_component"
    for comp, atoms, elements in [("SO4",5,["O","S"]),("PO4",5,["O","P"]),("NO3",4,["N","O"]),("CO3",4,["C","O"]),("OXY",2,["O"]),("NO",2,["N","O"]),("PEO",4,["H","O"]),("N2O",3,["N","O"]),("AZI",3,["N"])]:
        tests[f"low_frequency_{comp}"] = route_source(row(comp, atoms, elements, count=2))[0] == "provisional_source_ligand"
    tests["biolip_14"] = route_source(row("SO4",5,["O","S"],count=14))[0] == "provisional_source_ligand"
    tests["biolip_15"] = route_source(row("SO4",5,["O","S"],count=15))[0] == "excluded_biolip_artifact"
    tests["biolip_16"] = route_source(row("SO4",5,["O","S"],count=16))[0] == "excluded_biolip_artifact"
    tests["ccd_invalid_9eny_A1H57"] = route_source(row("A1H57",32,["C","H","N","O","S"],status="ccd_invalid"))[0] == "ccd_review"
    result = {"tests": tests, "test_count": len(tests), "validation_pass": all(tests.values()), "timestamp": utc()}
    (OUT / "tests/route_regression_tests.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["validation_pass"]: raise SystemExit(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


def run_full(workers, batch_size):
    if not json.loads((RUN / "audit/preflight.json").read_text())["preflight_pass"]: raise SystemExit("preflight failed")
    if not json.loads((OUT / "tests/route_regression_tests.json").read_text())["validation_pass"]: raise SystemExit("tests failed")
    cfg = yaml.safe_load((RUN / "config_snapshot.yaml").read_text()); cfg["runtime"].update({"workers": workers, "batch_size": batch_size})
    (RUN / "config_snapshot.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    items = load_items()
    completed = set()
    for complete in sorted((RUN / "work/batches").glob("batch_*/complete.json")):
        completed.update(json.loads(complete.read_text())["pdb_ids"])
    groups = iter(source_group_iterator())
    current = next(groups, None)
    order = {item["pdb_id"]: i for i, item in enumerate(items)}
    next_batch = len(list((RUN / "work/batches").glob("batch_*/complete.json")))
    started = time.time()
    status = {"status":"RUNNING","run_id":RUN_ID,"started_at":utc(),"workers":workers,"batch_size":batch_size,"input_entries":len(items)}
    (RUN / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    metal_path, biolip_path = OUT / "config/metal_elements_v1.yaml", OUT / "references/BioLiP_ligand_list_snapshot.tsv"
    pending_chunk = []
    with (RUN / "logs/run.log").open("a") as log, ProcessPoolExecutor(max_workers=workers, initializer=init_worker, initargs=(metal_path, biolip_path)) as pool:
        for item in items:
            pid = item["pdb_id"]
            rows = []
            if current and current[0] == pid:
                rows = current[1]; current = next(groups, None)
            elif current and order[current[0]] < order[pid]:
                raise RuntimeError(f"source order mismatch at {pid}: {current[0]}")
            if pid in completed: continue
            item = dict(item); item["source_rows"] = rows; pending_chunk.append(item)
            if len(pending_chunk) < batch_size: continue
            results = list(pool.map(process_item, pending_chunk, chunksize=1))
            batch = RUN / "work/batches" / f"batch_{next_batch:06d}"; batch.mkdir(parents=True)
            for table, fields in TABLES.items(): write_tsv(batch / f"{table}.tsv.gz", (row for result in results for row in result[table]), fields)
            ids = [x["pdb_id"] for x in pending_chunk]
            (batch / "complete.json").write_text(json.dumps({"batch_id":next_batch,"pdb_ids":ids,"completed_at":utc()}) + "\n")
            completed.update(ids); next_batch += 1; pending_chunk = []
            if len(completed) % 5000 < batch_size:
                progress = {"processed_entries":len(completed),"total_entries":len(items),"elapsed_seconds":round(time.time()-started,2),"updated_at":utc()}
                (RUN / "work/progress.json").write_text(json.dumps(progress, indent=2)+"\n"); log.write(json.dumps(progress)+"\n"); log.flush()
        if pending_chunk:
            results = list(pool.map(process_item, pending_chunk, chunksize=1))
            batch = RUN / "work/batches" / f"batch_{next_batch:06d}"; batch.mkdir(parents=True)
            for table, fields in TABLES.items(): write_tsv(batch / f"{table}.tsv.gz", (row for result in results for row in result[table]), fields)
            ids=[x["pdb_id"] for x in pending_chunk]; (batch/"complete.json").write_text(json.dumps({"batch_id":next_batch,"pdb_ids":ids,"completed_at":utc()})+"\n"); completed.update(ids)
    if current is not None: raise RuntimeError(f"unconsumed source group: {current[0]}")
    status.update({"status":"COMPLETED","completed_at":utc(),"runtime_seconds":round(time.time()-started,2),"processed_entries":len(completed)})
    (RUN / "status.json").write_text(json.dumps(status, indent=2)+"\n"); print(json.dumps(status, indent=2))


def merge_table(table, fields):
    out = RUN / "output" / f"{table}.tsv.gz"
    count = 0
    with gzip.open(out, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n"); writer.writeheader()
        for part in sorted((RUN / "work/batches").glob(f"batch_*/{table}.tsv.gz")):
            for row in iter_tsv(part): writer.writerow(row); count += 1
    return out, count


def deterministic_preview(path, fields, limit=1000):
    heap=[]
    keyfield = ("assembly_ligand_placement_id" if "assembly_ligand_placement_id" in fields
                else "source_ligand_instance_id" if "source_ligand_instance_id" in fields else "pdb_id")
    for sequence,row in enumerate(iter_tsv(path)):
        key=row.get(keyfield,""); score=int(hashlib.sha256(key.encode()).hexdigest(),16); item=(-score,key,sequence,row)
        if len(heap)<limit: heapq.heappush(heap,item)
        elif item>heap[0]: heapq.heapreplace(heap,item)
    rows=[x[3] for x in sorted(heap,key=lambda x:(-x[0],x[1],x[2]))]
    write_tsv(RUN/"output"/(path.name.replace(".tsv.gz","_preview.tsv")),rows,fields)


def parse_float_list(value):
    return [float(x) for x in value.split(",")]


def parse_matrix(value):
    return [parse_float_list(row) for row in value.split(";")]


def close_enough(a,b,abs_tol=1e-9,rel_tol=1e-9):
    # Operator vectors are serialized with 10 significant digits in the
    # placement table, so validation must account for that text round-trip.
    return len(a)==len(b) and all(
        math.isclose(x,y,abs_tol=abs_tol,rel_tol=rel_tol) for x,y in zip(a,b)
    )


def validate_placement_item(item):
    errors=Counter(); examples=[]
    try:
        block=gemmi.cif.read(item["mmcif_path"]).sole_block(); all_ids,operators=parse_operator_records(block)
        raw={}
        for idx,row in enumerate(category_records(block,"_pdbx_struct_assembly_gen."),start=1):
            raw[(row.get("assembly_id",""),f"assembly_gen_row_{idx:06d}")]=row
        for p in item["placements"]:
            key=(p["assembly_id"],p["assembly_gen_row_id"]); gen=raw.get(key)
            if not gen: errors["assembly_gen_row_missing"]+=1; continue
            asyms={x.strip() for x in gen.get("asym_id_list","").split(",") if x.strip()}
            if p["label_asym_id"] not in asyms: errors["label_asym_not_referenced"]+=1
            try: paths=expand_operator_paths(gen.get("oper_expression",""))
            except Exception: errors["oper_expression_parse_failure"]+=1; continue
            path=tuple(p["operator_path"].split("*"))
            if path not in paths: errors["operator_path_not_in_expression"]+=1; continue
            if any(op not in all_ids for op in path): errors["operator_id_missing"]+=1; continue
            if any(op not in operators for op in path): errors["operator_matrix_missing"]+=1; continue
            rotation,translation=composite_affine(path,operators)
            flat_expected=[x for row in rotation for x in row]; flat_actual=[x for row in parse_matrix(p["rotation_matrix"]) for x in row]
            if not close_enough(flat_expected,flat_actual): errors["rotation_matrix_mismatch"]+=1
            if not close_enough(translation,parse_float_list(p["translation_vector"])): errors["translation_vector_mismatch"]+=1
    except Exception as exc:
        errors["entry_validation_parse_failure"]+=1; examples.append(f"{item['pdb_id']}:{type(exc).__name__}:{exc}"[:500])
    return dict(errors),examples


def grouped_rows(path):
    current,rows=None,[]
    for row in iter_tsv(path):
        pid=row["pdb_id"]
        if current is not None and pid!=current: yield current,rows; rows=[]
        current=pid; rows.append(row)
    if current is not None: yield current,rows


def finalize_validate(workers):
    status=json.loads((RUN/"status.json").read_text())
    if status.get("status") not in {"COMPLETED", "VALIDATION_FAILED"}:
        raise SystemExit(f"full not completed or awaiting validation retry: {status}")
    merged={}
    for table,fields in TABLES.items():
        path,count=merge_table(table,fields); merged[table]={"path":str(path),"rows":count,"sha256":sha256(path)}; deterministic_preview(path,fields)
    init_worker(OUT/"config/metal_elements_v1.yaml",OUT/"references/BioLiP_ligand_list_snapshot.tsv")
    routes=Counter(); old_new=Counter(); old_new_pdb=defaultdict(set); old_new_ccd=defaultdict(set); comp_new=Counter(); comp_pdb=defaultdict(set)
    source_ids=set(); duplicate_source=0; route_errors=Counter(); provisional_ids=set()
    for row in iter_tsv(merged["source_terminal_routes"]["path"]):
        sid=row["source_ligand_instance_id"]; duplicate_source += sid in source_ids; source_ids.add(sid); routes[row["terminal_route"]]+=1
        expected=route_source(row)[0]
        if expected!=row["terminal_route"]: route_errors["independent_route_mismatch"]+=1
        key=(row["previous_terminal_route"],row["terminal_route"]); old_new[key]+=1; old_new_pdb[key].add(row["pdb_id"]); old_new_ccd[key].add(row["component_id"])
        comp_new[(row["component_id"],row["terminal_route"])]+=1; comp_pdb[(row["component_id"],row["terminal_route"])].add(row["pdb_id"])
        total=int_or_none(row["standard_total_atom_count"]); metals=set(row["element_set"].split(","))&G_METALS
        if row["terminal_route"]=="excluded_atomic_component" and not (row["ccd_identity_status"] in ACCEPTED_CCD and total==1): route_errors["atomic_rule_failure"]+=1
        if row["terminal_route"]=="excluded_metal_containing_component" and not (row["ccd_identity_status"] in ACCEPTED_CCD and total and total>=2 and metals): route_errors["metal_rule_failure"]+=1
        if row["terminal_route"]=="excluded_biolip_artifact" and not (row["ccd_identity_status"] in ACCEPTED_CCD and total and total>=2 and not metals and row["component_id"] in G_BIOLIP and int(row["source_instance_count_in_entry"])>=15): route_errors["biolip_rule_failure"]+=1
        if row["terminal_route"]=="provisional_source_ligand":
            provisional_ids.add(sid)
            if not (row["ccd_identity_status"] in ACCEPTED_CCD and total and total>=2 and not metals and not (row["component_id"] in G_BIOLIP and int(row["source_instance_count_in_entry"])>=15)): route_errors["provisional_rule_failure"]+=1
    delta_rows=[]
    for (old,new),count in sorted(old_new.items()): delta_rows.append({"old_route":old,"new_route":new,"instance_count":count,"unique_ccd_count":len(old_new_ccd[(old,new)]),"unique_pdb_count":len(old_new_pdb[(old,new)])})
    write_tsv(RUN/"release/route_delta_summary.tsv",delta_rows,["old_route","new_route","instance_count","unique_ccd_count","unique_pdb_count"])
    component_rows=[{"component_id":c,"new_route":r,"instance_count":n,"unique_pdb_count":len(comp_pdb[(c,r)])} for (c,r),n in sorted(comp_new.items())]
    write_tsv(RUN/"release/component_route_summary.tsv",component_rows,["component_id","new_route","instance_count","unique_pdb_count"])
    component_delta=Counter(); component_delta_pdb=defaultdict(set)
    for row in iter_tsv(merged["source_terminal_routes"]["path"]):
        key=(row["component_id"],row["previous_terminal_route"],row["terminal_route"]); component_delta[key]+=1; component_delta_pdb[key].add(row["pdb_id"])
    write_tsv(RUN/"release/component_route_delta_vs_20260804_full_01.tsv",({"component_id":c,"old_route":o,"new_route":n,"instance_count":v,"unique_pdb_count":len(component_delta_pdb[(c,o,n)])} for (c,o,n),v in sorted(component_delta.items())),["component_id","old_route","new_route","instance_count","unique_pdb_count"])
    placement_ids=set(); duplicate_placement=0; placement_source_fk=placement_assembly_fk=0; mapped_sources=set()
    retained={(r["pdb_id"],r["assembly_id"]) for r in iter_tsv(F1_ASSEMBLIES)}
    for row in iter_tsv(merged["ligand_assembly_logical_placements"]["path"]):
        pid=row["assembly_ligand_placement_id"]; duplicate_placement += pid in placement_ids; placement_ids.add(pid)
        placement_source_fk += row["source_ligand_instance_id"] not in provisional_ids
        placement_assembly_fk += (row["pdb_id"],row["assembly_id"]) not in retained
        mapped_sources.add(row["source_ligand_instance_id"])
    no_map_ids=set(); no_map_reasons=Counter(); duplicate_no_map=0; no_map_source_fk=0; unknown_no_map_reason=0
    allowed_no_map_reasons={"asym_not_referenced_by_retained_assembly","oper_expression_parse_failure","operator_id_missing","operator_matrix_missing","assembly_model_mismatch","other_mapping_failure"}
    for row in iter_tsv(merged["ligand_assembly_no_mapping"]["path"]):
        sid=row["source_ligand_instance_id"]; duplicate_no_map += sid in no_map_ids; no_map_ids.add(sid)
        no_map_source_fk += sid not in provisional_ids; no_map_reasons[row["reason_code"]]+=1
        unknown_no_map_reason += row["reason_code"] not in allowed_no_map_reasons
    placement_accounting_missing=len(provisional_ids-(mapped_sources|no_map_ids)); placement_accounting_overlap=len(mapped_sources&no_map_ids)
    items=load_items(); groups=iter(grouped_rows(Path(merged["ligand_assembly_logical_placements"]["path"]))); current=next(groups,None); tasks=[]
    for item in items:
        rows=[]
        if current and current[0]==item["pdb_id"]: rows=current[1]; current=next(groups,None)
        if rows:
            x=dict(item); x["placements"]=rows; tasks.append(x)
    placement_errors=Counter(); error_examples=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for errors,examples in pool.map(validate_placement_item,tasks,chunksize=20): placement_errors.update(errors); error_examples.extend(examples[:3])
    allowed_routes={"excluded_atomic_component","excluded_metal_containing_component","excluded_biolip_artifact","ccd_review","provisional_source_ligand"}
    forbidden_old=routes.get("excluded_simple_inorganic",0)+routes.get("inorganic_review",0)
    validation={
        "input_entries":248037,"input_source_instances":2652746,"water_source_instances":529742,
        "source_inventory_sha256":sha256(SOURCE_INPUT),"source_inventory_sha256_match":sha256(SOURCE_INPUT)==EXPECTED_HASHES["source_inventory"],
        "terminal_route_counts":dict(routes),"terminal_accounting_sum":sum(routes.values()),"silent_drop":2652746-sum(routes.values()),
        "unknown_terminal_route_count":sum(v for k,v in routes.items() if k not in allowed_routes),"forbidden_old_route_count":forbidden_old,
        "duplicate_source_instance_id":duplicate_source,"independent_route_errors":dict(route_errors),
        "duplicate_placement_id":duplicate_placement,"placement_source_fk_failure":placement_source_fk,"placement_retained_assembly_fk_failure":placement_assembly_fk,
        "placement_accounting_missing":placement_accounting_missing,"placement_accounting_overlap":placement_accounting_overlap,
        "duplicate_no_mapping_source_id":duplicate_no_map,"no_mapping_source_fk_failure":no_map_source_fk,
        "unknown_no_mapping_reason_count":unknown_no_map_reason,
        "independent_raw_assembly_placement_errors":dict(placement_errors),"placement_validation_error_examples":error_examples,
        "no_mapping_reason_counts":dict(no_map_reasons),"regression_tests_pass":json.loads((OUT/"tests/route_regression_tests.json").read_text())["validation_pass"],
        "assembly_coordinates_generated":False,"protein_ligand_pairs_generated":False,"distance_calculated":False,"interaction_tools_executed":False,
    }
    zeros=[validation["silent_drop"],validation["unknown_terminal_route_count"],forbidden_old,duplicate_source,duplicate_placement,
           placement_source_fk,placement_assembly_fk,placement_accounting_missing,placement_accounting_overlap,
           duplicate_no_map,no_map_source_fk,unknown_no_map_reason,*route_errors.values(),*placement_errors.values()]
    validation["validation_pass"]=all(x==0 for x in zeros) and validation["terminal_accounting_sum"]==2652746 and validation["source_inventory_sha256_match"] and validation["regression_tests_pass"]
    (RUN/"release/validation_report.json").write_text(json.dumps(validation,indent=2)+"\n")
    route_delta_json={"transitions":delta_rows,"specified_components":{}}
    specified={"SO4","PO4","NO3","CO3","OXY","NO","PEO","NH4","SF4","FES","VO4","WO4","CPT"}
    for (c,o,n),v in component_delta.items():
        if c in specified: route_delta_json["specified_components"].setdefault(c,[]).append({"old_route":o,"new_route":n,"instance_count":v,"unique_pdb_count":len(component_delta_pdb[(c,o,n)])})
    (RUN/"release/route_delta_summary.json").write_text(json.dumps(route_delta_json,indent=2)+"\n")
    if not validation["validation_pass"]:
        (RUN/"status.json").write_text(json.dumps({**status,"status":"VALIDATION_FAILED","validation":validation},indent=2)+"\n")
        raise SystemExit(json.dumps(validation,indent=2))
    summary={"run_id":RUN_ID,"rule_version":RULE_VERSION,"runtime_seconds":status.get("runtime_seconds"),"input_source_instances":2652746,
             "water_source_instances":529742,"terminal_route_counts":dict(routes),"provisional_sources_with_mapping":len(mapped_sources),
             "logical_placement_count":merged["ligand_assembly_logical_placements"]["rows"],"no_mapping_count":merged["ligand_assembly_no_mapping"]["rows"],
             "metal_policy_version":METAL_POLICY_VERSION,"metal_elements_sha256":sha256(OUT/"config/metal_elements_v1.yaml"),
             "biolip_sha256":sha256(OUT/"references/BioLiP_ligand_list_snapshot.tsv"),"ccd_metadata":json.loads((OUT/"references/CCD_snapshot_metadata.json").read_text()),
             "parquet_unavailable_in_frozen_environment":True,"validation_pass":True}
    (RUN/"release/release_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    interface={"stage":"filter_2_ligand_qualification_v3","run_id":RUN_ID,"status":"FROZEN","formal_tables":merged,
               "provisional_source_ligands":merged["provisional_source_ligands"],"ligand_assembly_logical_placements":merged["ligand_assembly_logical_placements"],
               "coordinates_materialized":False,"pairs_constructed":False,"validation_pass":True}
    (RUN/"release/downstream_interface.json").write_text(json.dumps(interface,indent=2)+"\n")
    run_metadata={"run_id":RUN_ID,"host":platform.node(),"python":sys.version,"gemmi":gemmi.__version__,"code_sha256":sha256(OUT/"scripts/filter2_v3_pipeline.py"),"completed_at":utc()}
    (RUN/"run_metadata.json").write_text(json.dumps(run_metadata,indent=2)+"\n")
    extras=[RUN/"release/route_delta_summary.tsv",RUN/"release/component_route_summary.tsv",RUN/"release/component_route_delta_vs_20260804_full_01.tsv",RUN/"release/route_delta_summary.json",RUN/"release/validation_report.json",RUN/"release/release_summary.json",RUN/"release/downstream_interface.json",RUN/"config_snapshot.yaml",OUT/"config/metal_elements_v1.yaml",OUT/"references/BioLiP_ligand_list_snapshot.tsv",OUT/"references/BioLiP_ligand_list_metadata.json",OUT/"references/CCD_snapshot_metadata.json",RUN/"input/input_provenance.json",RUN/"run_metadata.json"]
    extras.extend(sorted((RUN/"output").glob("*_preview.tsv")))
    extras.extend(sorted((OUT/"schemas").glob("*_schema.json")))
    manifest=[]
    for table,data in merged.items(): manifest.append({"relative_path":str(Path(data["path"]).relative_to(OUT)),"file_role":table,"row_count":data["rows"],"size_bytes":Path(data["path"]).stat().st_size,"sha256":data["sha256"]})
    for path in extras: manifest.append({"relative_path":str(path.relative_to(OUT)),"file_role":path.name,"row_count":"","size_bytes":path.stat().st_size,"sha256":sha256(path)})
    write_tsv(RUN/"release/output_manifest.tsv",manifest,["relative_path","file_role","row_count","size_bytes","sha256"])
    hash_paths=[OUT/row["relative_path"] for row in manifest]+[RUN/"release/output_manifest.tsv"]
    (RUN/"release/SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(OUT)}\n" for p in hash_paths))
    frozen={"status":"FROZEN","run_id":RUN_ID,"frozen_at":utc(),"validation_pass":True,"accounting_pass":True,"manifest_sha256":sha256(RUN/"release/output_manifest.tsv"),"code_sha256":sha256(OUT/"scripts/filter2_v3_pipeline.py")}
    (RUN/"_FROZEN.json").write_text(json.dumps(frozen,indent=2)+"\n")
    (RUN/"status.json").write_text(json.dumps({**status,"status":"FROZEN","frozen_at":frozen["frozen_at"],"validation":validation},indent=2)+"\n")
    (OUT/"CURRENT_RUN.json").write_text(json.dumps({"current_run_id":RUN_ID,"status":"FROZEN","relative_path":f"runs/{RUN_ID}","manifest_sha256":frozen["manifest_sha256"],"updated_at":utc()},indent=2)+"\n")
    current=OUT/"current"
    if current.exists() or current.is_symlink(): current.unlink()
    current.symlink_to(Path("runs")/RUN_ID)
    print(json.dumps(summary,indent=2))


def main():
    parser=argparse.ArgumentParser(); sub=parser.add_subparsers(dest="command",required=True)
    sub.add_parser("setup"); sub.add_parser("preflight"); sub.add_parser("selftest")
    run=sub.add_parser("run-full"); run.add_argument("--workers",type=int,default=24); run.add_argument("--batch-size",type=int,default=200)
    fin=sub.add_parser("finalize-validate"); fin.add_argument("--workers",type=int,default=24)
    args=parser.parse_args()
    if args.command=="setup": setup()
    elif args.command=="preflight": preflight()
    elif args.command=="selftest": selftest()
    elif args.command=="run-full": run_full(args.workers,args.batch_size)
    elif args.command=="finalize-validate": finalize_validate(args.workers)


if __name__=="__main__": main()
