#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import math
import os
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN = ROOT / "runs/20260812_full_01"
P3 = Path("/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/runs/20260811_full_01/output")
P2 = Path("/root/autodl-tmp/benchmark_1.0/processing_2_assembly_ready_structure_preparation/runs/20260810_full_01/output")
VAL = RUN / "work/normalized_validation/by_bucket"
MMCIF = Path("/root/autodl-tmp/pdb_archive_v2/mmCIF")
OUT = RUN / "work/quality_batches"

EXPECTED_HEAVY = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6,
    "GLN": 9, "GLU": 9, "GLY": 4, "HIS": 10, "ILE": 8,
    "LEU": 8, "LYS": 9, "MET": 8, "MSE": 8, "PHE": 11,
    "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12,
    "VAL": 7, "SEC": 6, "PYL": 9,
}

ENTRY_THRESHOLDS = {"resolution": 3.0, "r_work": 0.40, "r_free": 0.45, "r_gap": 0.05}
LIGAND_THRESHOLDS = {"rscc": 0.80, "rsr": 0.30, "occupancy": 0.80}
POCKET_THRESHOLDS = {"rscc": 0.80, "rsr": 0.30, "occupancy": 0.80}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", ".", "?", "none", "false", "nan"}:
        return ""
    return text


def finite(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def bool_value(value) -> bool:
    return bool(value) if not pd.isna(value) else False


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def bucket_paths(root: Path, bucket: int) -> list[str]:
    directory = root / f"bucket_id={bucket:03d}"
    return [str(path) for path in sorted(directory.glob("*.parquet"))]


def read_bucket(root: Path, bucket: int, columns: list[str] | None = None) -> pd.DataFrame:
    paths = bucket_paths(root, bucket)
    if not paths:
        return pd.DataFrame(columns=columns or [])
    table = ds.dataset(paths, format="parquet").to_table(columns=columns)
    return table.to_pandas(split_blocks=True, self_destruct=True)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), temporary, compression="zstd")
    os.replace(temporary, path)


def method_class(methods: list[str]) -> str:
    values = [clean(value).lower() for value in methods if clean(value)]
    if len(values) != 1:
        return "hybrid_or_multiple_methods" if values else "unknown"
    value = values[0]
    if "x-ray" in value or "xray" in value:
        return "xray"
    if "electron microscopy" in value or "cryo" in value:
        return "cryo_em"
    if "nmr" in value:
        return "nmr"
    return "other"


def category_rows(block, category_name: str) -> list[dict]:
    try:
        category = block.get_mmcif_category(category_name)
    except Exception:
        return []
    if not category:
        return []
    keys = list(category)
    size = max((len(category[key]) for key in keys), default=0)
    return [{key: clean(category[key][index]) for key in keys} for index in range(size)]


def parse_structure_metadata(pdb_id: str) -> tuple[dict, list[dict]]:
    path = MMCIF / pdb_id[1:3] / f"{pdb_id}.cif.gz"
    base = {
        "pdb_id": pdb_id,
        "mmcif_path": str(path),
        "mmcif_parse_status": "PARSE_FAILED",
        "mmcif_parse_error": "",
        "experimental_method_raw": "",
        "experimental_method_class": "unknown",
        "model_count": None,
    }
    unobserved = []
    if not path.exists():
        base["mmcif_parse_status"] = "MISSING"
        return base, unobserved
    try:
        block = gemmi.cif.read_file(str(path)).sole_block()
        methods = category_rows(block, "_exptl.")
        method_values = [row.get("method", "") for row in methods]
        base["experimental_method_raw"] = ";".join(value for value in method_values if value)
        base["experimental_method_class"] = method_class(method_values)
        atom = block.get_mmcif_category("_atom_site.")
        models = {clean(value) for value in atom.get("pdbx_PDB_model_num", []) if clean(value)}
        base["model_count"] = len(models) if models else 1
        for source, category_name in (
            ("unobserved_residue", "_pdbx_unobs_or_zero_occ_residues."),
            ("unobserved_atom", "_pdbx_unobs_or_zero_occ_atoms."),
        ):
            for row in category_rows(block, category_name):
                unobserved.append({
                    "pdb_id": pdb_id,
                    "unobserved_type": source,
                    "model_id": row.get("PDB_model_num", ""),
                    "polymer_flag": row.get("polymer_flag", ""),
                    "occupancy_flag": row.get("occupancy_flag", ""),
                    "auth_asym_id": row.get("auth_asym_id", ""),
                    "label_asym_id": row.get("label_asym_id", ""),
                    "auth_seq_id": row.get("auth_seq_id", ""),
                    "label_seq_id": row.get("label_seq_id", ""),
                    "insertion_code": row.get("PDB_ins_code", ""),
                    "component_id": row.get("auth_comp_id", row.get("label_comp_id", "")),
                    "atom_id": row.get("auth_atom_id", row.get("label_atom_id", "")),
                })
        base["mmcif_parse_status"] = "PARSE_SUCCESS"
    except Exception as exc:
        base["mmcif_parse_error"] = f"{type(exc).__name__}:{exc}"[:2000]
    return base, unobserved


def residue_id_parts(value: str) -> tuple[str, str, str, str]:
    parts = str(value).split("|")
    parts += [""] * (4 - len(parts))
    return clean(parts[0]), clean(parts[1]), clean(parts[2]), clean(parts[3]).upper()


def validation_candidates(frame: pd.DataFrame) -> tuple[dict, dict]:
    auth = defaultdict(list)
    label = defaultdict(list)
    for row in frame.to_dict("records"):
        model = clean(row.get("model_id")) or "1"
        component = clean(row.get("component_id")).upper()
        auth_key = (
            clean(row.get("pdb_id")).lower(), model, clean(row.get("auth_asym_id")),
            clean(row.get("auth_seq_id")), clean(row.get("insertion_code")), component,
        )
        label_key = (
            clean(row.get("pdb_id")).lower(), model, clean(row.get("label_asym_id")),
            clean(row.get("label_seq_id")), component,
        )
        auth[auth_key].append(row)
        label[label_key].append(row)
    return auth, label


def choose_alt(rows: list[dict], selected_alt: str) -> tuple[dict | None, str]:
    if not rows:
        return None, "NOT_FOUND"
    selected_alt = clean(selected_alt)
    exact = [row for row in rows if clean(row.get("alt_id")) == selected_alt]
    if len(exact) == 1:
        return exact[0], "MAPPED_EXACT_ALTLOC"
    if len(exact) > 1:
        return None, "AMBIGUOUS_EXACT_ALTLOC"
    blanks = [row for row in rows if not clean(row.get("alt_id"))]
    if len(rows) == 1:
        return rows[0], "MAPPED_SINGLE_ROW_ALTLOC_FALLBACK"
    if not selected_alt and len(blanks) == 1:
        return blanks[0], "MAPPED_BLANK_ALTLOC"
    return None, "AMBIGUOUS_ALTLOC"


def map_validation(
    auth_index: dict,
    label_index: dict,
    *,
    pdb_id: str,
    model_id: str,
    auth_asym_id: str,
    label_asym_id: str,
    auth_seq_id: str,
    label_seq_id: str,
    insertion_code: str,
    component_id: str,
    selected_alt: str,
) -> tuple[dict | None, str, str]:
    pdb_id = clean(pdb_id).lower()
    model_id = clean(model_id) or "1"
    component_id = clean(component_id).upper()
    auth_key = (pdb_id, model_id, clean(auth_asym_id), clean(auth_seq_id), clean(insertion_code), component_id)
    row, status = choose_alt(auth_index.get(auth_key, []), selected_alt)
    if row is not None:
        return row, status, "AUTH_KEY"
    if status.startswith("AMBIGUOUS"):
        return None, status, "AUTH_KEY"
    label_key = (pdb_id, model_id, clean(label_asym_id), clean(label_seq_id), component_id)
    row, status = choose_alt(label_index.get(label_key, []), selected_alt)
    return row, status, "LABEL_KEY" if row is not None or status.startswith("AMBIGUOUS") else "NONE"


def ligand_context(ligand_atoms: pd.DataFrame, topology: pd.DataFrame) -> pd.DataFrame:
    if ligand_atoms.empty:
        return pd.DataFrame()
    key = "filter_2_ligand_assembly_placement_id"
    columns = [
        key, "filter_2_source_ligand_instance_id", "pdb_id", "assembly_id", "model_id",
        "component_id", "entity_id", "label_asym_id", "auth_asym_id", "label_seq_id",
        "auth_seq_id", "insertion_code", "operator_path", "selected_altloc_original",
        "expected_heavy_atom_count", "observed_heavy_atom_count",
    ]
    first = ligand_atoms[columns].drop_duplicates(key).copy()
    occupancy = ligand_atoms.groupby(key, sort=False)["occupancy"].mean().rename("coordinate_mean_occupancy")
    first = first.merge(occupancy, on=key, how="left")
    selected_alt = (
        ligand_atoms.groupby(key, sort=False)["selected_altloc_original"]
        .agg(lambda values: Counter(clean(value) for value in values if clean(value)).most_common(1)[0][0]
             if any(clean(value) for value in values) else "")
        .rename("_selected_altloc")
    )
    first = first.merge(selected_alt, on=key, how="left")
    first["selected_altloc_original"] = first["_selected_altloc"]
    first = first.drop(columns=["_selected_altloc"])
    topo_columns = [
        "source_ligand_instance_id", "mapping_status", "missing_heavy_atom_count",
        "topology_status", "rdkit_parse_success", "rdkit_sanitize_success", "rdkit_error",
    ]
    first = first.merge(topology[topo_columns].drop_duplicates("source_ligand_instance_id"),
                        left_on="filter_2_source_ligand_instance_id",
                        right_on="source_ligand_instance_id", how="left")
    return first


def build_receptor_residue_context(receptor_atoms: pd.DataFrame, needed: set[tuple[str, str]]) -> pd.DataFrame:
    if receptor_atoms.empty or not needed:
        return pd.DataFrame()
    label_seq, auth_seq, ins, comp = zip(*[residue_id_parts(value) for _, value in needed])
    needed_frame = pd.DataFrame({
        "filter_1_chain_instance_id": [item[0] for item in needed],
        "protein_residue_id": [item[1] for item in needed],
        "target_label_seq_id": label_seq,
        "target_auth_seq_id": auth_seq,
        "target_insertion_code": ins,
        "target_component_id": comp,
    })
    atoms = receptor_atoms.copy()
    atoms["_label"] = atoms["label_seq_id"].map(clean)
    atoms["_auth"] = atoms["auth_seq_id"].map(clean)
    atoms["_ins"] = atoms["insertion_code"].map(clean)
    atoms["_comp"] = atoms["label_comp_id"].map(clean).str.upper()
    merged = needed_frame.merge(
        atoms,
        left_on=["filter_1_chain_instance_id", "target_label_seq_id", "target_auth_seq_id", "target_insertion_code", "target_component_id"],
        right_on=["filter_1_chain_instance_id", "_label", "_auth", "_ins", "_comp"],
        how="left",
    )
    group_columns = ["filter_1_chain_instance_id", "protein_residue_id"]
    rows = []
    for (chain_id, residue_id), group in merged.groupby(group_columns, sort=False, dropna=False):
        valid = group[group["pdb_id"].notna()]
        base = valid.iloc[0] if not valid.empty else group.iloc[0]
        component = clean(base.get("target_component_id")).upper()
        observed = valid.loc[valid["type_symbol"].fillna("").str.upper() != "H", "label_atom_id"].nunique()
        # CCD monomer tables include polymer leaving/terminal atoms.  Using their
        # raw atom count for an internal residue creates a systematic +1 error.
        expected = EXPECTED_HEAVY.get(component)
        selected_values = [clean(value) for value in valid.get("selected_altloc_original", []) if clean(value)]
        selected_alt = Counter(selected_values).most_common(1)[0][0] if selected_values else ""
        rows.append({
            "chain_instance_id": chain_id,
            "protein_residue_id": residue_id,
            "pdb_id": clean(base.get("pdb_id")).lower(),
            "model_id": clean(base.get("model_id")) or "1",
            "auth_asym_id": clean(base.get("auth_asym_id")),
            "label_asym_id": clean(base.get("label_asym_id")),
            "auth_seq_id": clean(base.get("target_auth_seq_id")),
            "label_seq_id": clean(base.get("target_label_seq_id")),
            "insertion_code": clean(base.get("target_insertion_code")),
            "component_id": component,
            "selected_altloc_original": selected_alt,
            "observed_heavy_atom_count": int(observed),
            "expected_heavy_atom_count": expected,
            "unresolved_heavy_atom_count": None if expected is None else max(0, expected - int(observed)),
            "completeness_status": "UNKNOWN_COMPONENT" if expected is None else ("COMPLETE" if observed >= expected else "INCOMPLETE"),
        })
    return pd.DataFrame(rows)


def append_validation_fields(base: dict, validation: dict | None) -> dict:
    fields = (
        "rsr", "rscc", "rsrz", "mean_occupancy", "mean_b_factor", "natoms_eds",
        "density_outlier", "geometry_outlier", "chirality_outlier", "clash_outlier",
        "bond_outlier_count", "angle_outlier_count", "chirality_outlier_count", "clash_outlier_count",
    )
    for field in fields:
        base[field] = None if validation is None else validation.get(field)
    return base


def metric_available(row: dict, names: tuple[str, ...]) -> bool:
    return all(finite(row.get(name)) is not None for name in names)


def quality_residue_row(context: dict, auth_index: dict, label_index: dict) -> dict:
    mapped, mapping_status, mapping_key = map_validation(
        auth_index, label_index,
        pdb_id=context.get("pdb_id"), model_id=context.get("model_id"),
        auth_asym_id=context.get("auth_asym_id"), label_asym_id=context.get("label_asym_id"),
        auth_seq_id=context.get("auth_seq_id"), label_seq_id=context.get("label_seq_id"),
        insertion_code=context.get("insertion_code"), component_id=context.get("component_id"),
        selected_alt=context.get("selected_altloc_original"),
    )
    row = dict(context)
    row.update({"validation_mapping_status": mapping_status, "validation_mapping_key": mapping_key})
    append_validation_fields(row, mapped)
    missing_metrics = not metric_available(row, ("rsr", "rscc", "mean_occupancy"))
    completeness = row.get("completeness_status")
    if mapped is None:
        quality = "VALIDATION_MAPPING_FAILED"
    elif missing_metrics or completeness == "UNKNOWN_COMPONENT":
        quality = "QUALITY_DATA_UNAVAILABLE"
    elif completeness == "INCOMPLETE":
        quality = "QUALITY_HARD_FAIL"
    elif (
        finite(row.get("rsrz")) is not None and finite(row.get("rsrz")) > 2.0
        or finite(row.get("rscc")) < 0.80
        or finite(row.get("rsr")) > 0.30
        or bool_value(row.get("density_outlier"))
        or bool_value(row.get("geometry_outlier"))
        or bool_value(row.get("chirality_outlier"))
        or bool_value(row.get("clash_outlier"))
    ):
        quality = "QUALITY_SOFT_FAIL"
    else:
        quality = "QUALITY_PASS"
    row["residue_quality_status"] = quality
    return row


def structural_gap_flags(unobserved: list[dict], binding: pd.DataFrame, pocket: pd.DataFrame) -> tuple[dict, list[dict]]:
    binding_positions = defaultdict(set)
    pocket_positions = defaultdict(set)
    for frame, target in ((binding, binding_positions), (pocket, pocket_positions)):
        for row in frame.to_dict("records"):
            label, auth, _, _ = residue_id_parts(row["protein_residue_id"])
            for value in (label, auth):
                try:
                    target[(row.get("chain_instance_id", ""), "position")].add(int(value))
                except ValueError:
                    pass
    by_source_chain = defaultdict(list)
    for row in unobserved:
        if row["unobserved_type"] == "unobserved_residue":
            by_source_chain[(row["pdb_id"], row["model_id"] or "1", row["auth_asym_id"], row["label_asym_id"])].append(row)
    flags = {}
    details = []
    chain_metadata = {}
    for frame in (binding, pocket):
        for row in frame.to_dict("records"):
            chain_metadata.setdefault(row.get("chain_instance_id", ""), row)
    for chain_id, meta in chain_metadata.items():
        candidates = by_source_chain.get((clean(meta.get("pdb_id")).lower(), clean(meta.get("model_id")) or "1", clean(meta.get("auth_asym_id")), clean(meta.get("label_asym_id"))), [])
        bind = binding_positions.get((chain_id, "position"), set())
        pocket_set = pocket_positions.get((chain_id, "position"), set())
        critical = False
        warning = False
        for gap in candidates:
            positions = []
            for field in ("label_seq_id", "auth_seq_id"):
                try:
                    positions.append(int(clean(gap.get(field))))
                except ValueError:
                    pass
            gap_critical = any((position - 1 in bind and position + 1 in bind) for position in positions)
            gap_warning = any(any(abs(position - observed) <= 1 for observed in pocket_set) for position in positions)
            if gap_critical or gap_warning:
                details.append({
                    "chain_instance_id": chain_id,
                    **gap,
                    "critical_interface_gap": gap_critical,
                    "adjacent_pocket_gap_warning": gap_warning,
                })
            critical = critical or gap_critical
            warning = warning or gap_warning
        flags[chain_id] = {"critical_interface_gap": critical, "adjacent_pocket_gap_warning": warning}
    return flags, details


def process_bucket(bucket: int) -> dict:
    started = time.time()
    bucket_out = OUT / f"bucket_id={bucket:03d}"
    marker = bucket_out / "_COMPLETE.json"
    if marker.exists():
        return json.loads(marker.read_text())
    bucket_out.mkdir(parents=True, exist_ok=True)

    pairs = read_bucket(P3 / "provisional_pairs", bucket)
    binding = read_bucket(P3 / "binding_residues", bucket)
    pocket = read_bucket(P3 / "pair_pocket_residues", bucket)
    chains = read_bucket(P3 / "contact_supported_chains", bucket)
    ligand_atoms = read_bucket(P2 / "prepared_ligand_assembly_atoms", bucket, [
        "filter_2_ligand_assembly_placement_id", "filter_2_source_ligand_instance_id", "pdb_id",
        "assembly_id", "model_id", "component_id", "entity_id", "label_asym_id", "auth_asym_id",
        "label_seq_id", "auth_seq_id", "insertion_code", "operator_path", "selected_altloc_original",
        "expected_heavy_atom_count", "observed_heavy_atom_count", "occupancy",
    ])
    topology = read_bucket(P2 / "ligand_topology_validation", bucket)
    receptor_atoms = read_bucket(P2 / "prepared_receptor_assembly_atoms", bucket, [
        "filter_1_chain_instance_id", "pdb_id", "model_id", "entity_id", "label_asym_id", "auth_asym_id",
        "label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id", "type_symbol", "label_atom_id",
        "selected_altloc_original",
    ])
    entry_validation = read_bucket(VAL / "entry_validation", bucket)
    residue_validation = read_bucket(VAL / "residue_validation", bucket)

    formal_placements = set(pairs["ligand_assembly_placement_id"])
    binding = binding[binding["ligand_assembly_placement_id"].isin(formal_placements)].copy()
    chains = chains[chains["ligand_assembly_placement_id"].isin(formal_placements)].copy()
    ligand_atoms = ligand_atoms[
        ligand_atoms["filter_2_ligand_assembly_placement_id"].isin(formal_placements)
    ].copy()
    formal_sources = set(ligand_atoms["filter_2_source_ligand_instance_id"])
    topology = topology[topology["source_ligand_instance_id"].isin(formal_sources)].copy()

    pdb_ids = sorted(set(pairs["pdb_id"].str.lower()))
    metadata_rows = []
    unobserved_rows = []
    for pdb_id in pdb_ids:
        metadata, unobserved = parse_structure_metadata(pdb_id)
        metadata_rows.append(metadata)
        unobserved_rows.extend(unobserved)
    metadata = {row["pdb_id"]: row for row in metadata_rows}
    entry_map = {clean(row["pdb_id"]).lower(): row for row in entry_validation.to_dict("records")}
    auth_index, label_index = validation_candidates(residue_validation)

    ligand = ligand_context(ligand_atoms, topology)
    ligand_map = {row["filter_2_ligand_assembly_placement_id"]: row for row in ligand.to_dict("records")}
    ligand_quality_rows = []
    for placement_id, context in ligand_map.items():
        mapped, mapping_status, mapping_key = map_validation(
            auth_index, label_index,
            pdb_id=context["pdb_id"], model_id=context["model_id"],
            auth_asym_id=context["auth_asym_id"], label_asym_id=context["label_asym_id"],
            auth_seq_id=context["auth_seq_id"], label_seq_id=context["label_seq_id"],
            insertion_code=context["insertion_code"], component_id=context["component_id"],
            selected_alt=context["selected_altloc_original"],
        )
        row = dict(context)
        row["ligand_assembly_placement_id"] = placement_id
        row["validation_mapping_status"] = mapping_status
        row["validation_mapping_key"] = mapping_key
        append_validation_fields(row, mapped)
        ligand_quality_rows.append(row)
    ligand_quality = pd.DataFrame(ligand_quality_rows)
    ligand_quality_map = {row["ligand_assembly_placement_id"]: row for row in ligand_quality_rows}

    needed = set()
    for frame in (binding, pocket):
        needed.update(zip(frame["chain_instance_id"], frame["protein_residue_id"]))
    residue_context = build_receptor_residue_context(receptor_atoms, needed)
    context_map = {(row["chain_instance_id"], row["protein_residue_id"]): row for row in residue_context.to_dict("records")}
    quality_context_map = {
        key: quality_residue_row(value, auth_index, label_index)
        for key, value in context_map.items()
    }

    def expand(frame: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for source in frame.to_dict("records"):
            quality = quality_context_map.get((source["chain_instance_id"], source["protein_residue_id"]))
            row = dict(source)
            if quality is None:
                row.update({"validation_mapping_status": "RECEPTOR_CONTEXT_NOT_FOUND", "residue_quality_status": "VALIDATION_MAPPING_FAILED"})
            else:
                row.update({key: value for key, value in quality.items() if key not in {"chain_instance_id", "protein_residue_id"}})
            rows.append(row)
        return pd.DataFrame(rows)

    binding_quality = expand(binding)
    pocket_quality = expand(pocket)
    gap_flags, gap_details = structural_gap_flags(unobserved_rows, binding_quality, pocket_quality)

    chain_quality_rows = []
    binding_by_placement_chain = binding_quality.groupby(["ligand_assembly_placement_id", "chain_instance_id"], sort=False)
    for source in chains.to_dict("records"):
        key = (source["ligand_assembly_placement_id"], source["chain_instance_id"])
        try:
            group = binding_by_placement_chain.get_group(key)
        except KeyError:
            group = pd.DataFrame()
        counts = Counter(group.get("residue_quality_status", pd.Series(dtype=str)).tolist())
        reliable = counts.get("QUALITY_PASS", 0)
        gap = gap_flags.get(source["chain_instance_id"], {})
        if gap.get("critical_interface_gap"):
            status = "CHAIN_QUALITY_HARD_FAIL"
        elif any(name in counts for name in ("VALIDATION_MAPPING_FAILED", "QUALITY_DATA_UNAVAILABLE")):
            status = "CHAIN_QUALITY_DATA_UNAVAILABLE"
        elif reliable >= 2:
            status = "QUALITY_SUPPORTED_CHAIN"
        else:
            status = "CHAIN_QUALITY_SUPPORT_LOST"
        chain_quality_rows.append({
            **source,
            "quality_pass_binding_residue_count": reliable,
            "quality_soft_fail_binding_residue_count": counts.get("QUALITY_SOFT_FAIL", 0),
            "quality_hard_fail_binding_residue_count": counts.get("QUALITY_HARD_FAIL", 0),
            "quality_unavailable_binding_residue_count": counts.get("VALIDATION_MAPPING_FAILED", 0) + counts.get("QUALITY_DATA_UNAVAILABLE", 0),
            "critical_interface_gap": bool(gap.get("critical_interface_gap")),
            "adjacent_pocket_gap_warning": bool(gap.get("adjacent_pocket_gap_warning")),
            "chain_quality_status": status,
        })
    chain_quality = pd.DataFrame(chain_quality_rows)

    pocket_by_pair = pocket_quality.groupby("pair_id", sort=False)
    chains_by_placement = chain_quality.groupby("ligand_assembly_placement_id", sort=False)
    pair_rows = []
    for pair in pairs.to_dict("records"):
        pair_id = pair["pair_id"]
        placement = pair["ligand_assembly_placement_id"]
        pdb_id = clean(pair["pdb_id"]).lower()
        reasons = []
        warnings = []
        meta = metadata.get(pdb_id, {})
        entry = entry_map.get(pdb_id)
        ligand_row = ligand_quality_map.get(placement)
        try:
            pocket_group = pocket_by_pair.get_group(pair_id)
        except KeyError:
            pocket_group = pd.DataFrame()
        try:
            chain_group = chains_by_placement.get_group(placement)
        except KeyError:
            chain_group = pd.DataFrame()

        method = meta.get("experimental_method_class", "unknown")
        terminal = ""
        entry_pass = ligand_pass = pocket_pass = chain_pass = False
        mapping_status = "VALIDATION_MAPPING_OK"
        if meta.get("mmcif_parse_status") != "PARSE_SUCCESS":
            terminal = "FILTER3_TECHNICAL_FAILURE"
            reasons.append("MMCIF_METADATA_PARSE_FAILED")
        elif method != "xray":
            terminal = "FILTER3_NON_XRAY_PROTOCOL_PENDING"
            reasons.append(f"METHOD_{method.upper()}_PENDING")
        elif entry is None or entry.get("parse_status") != "PARSE_SUCCESS":
            terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
            reasons.append("ENTRY_VALIDATION_UNAVAILABLE")
        else:
            entry_metrics = [finite(entry.get(name)) for name in ("resolution", "r_work", "r_free", "r_free_minus_r_work")]
            if any(value is None for value in entry_metrics):
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                reasons.append("ENTRY_QUALITY_METRIC_MISSING")
            else:
                resolution, r_work, r_free, r_gap = entry_metrics
                entry_pass = resolution <= 3.0 and r_work <= 0.40 and r_free <= 0.45 and abs(r_gap) <= 0.05
                if not entry_pass:
                    reasons.append("ENTRY_QUALITY_HARD_GATE_FAILED")

        if not terminal and ligand_row is None:
            terminal = "FILTER3_TECHNICAL_FAILURE"
            reasons.append("LIGAND_CONTEXT_MISSING")
        elif not terminal:
            if not str(ligand_row.get("validation_mapping_status", "")).startswith("MAPPED"):
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                mapping_status = "LIGAND_MAPPING_FAILED"
                reasons.append("LIGAND_VALIDATION_MAPPING_FAILED")
            elif not metric_available(ligand_row, ("rscc", "rsr", "mean_occupancy")):
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                reasons.append("LIGAND_QUALITY_METRIC_MISSING")
            else:
                missing = finite(ligand_row.get("missing_heavy_atom_count"))
                if missing is None:
                    terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                    reasons.append("LIGAND_COMPLETENESS_UNAVAILABLE")
                else:
                    ligand_pass = (
                        finite(ligand_row["rscc"]) >= 0.80
                        and finite(ligand_row["rsr"]) <= 0.30
                        and finite(ligand_row["mean_occupancy"]) >= 0.80
                        and missing == 0
                    )
                    if not ligand_pass:
                        reasons.append("LIGAND_QUALITY_HARD_GATE_FAILED")

        pocket_metrics = {"mean_rscc": None, "min_rscc": None, "mean_rsr": None, "max_rsr": None,
                          "mean_occupancy": None, "min_occupancy": None, "rsrz_outlier_count": None,
                          "unresolved_heavy_atom_count": None, "residue_count": len(pocket_group)}
        if not terminal:
            if pocket_group.empty:
                terminal = "FILTER3_TECHNICAL_FAILURE"
                reasons.append("POCKET_ROWS_MISSING")
            elif (pocket_group["validation_mapping_status"].astype(str).str.startswith("MAPPED").sum() != len(pocket_group)):
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                mapping_status = "POCKET_MAPPING_PARTIAL"
                reasons.append("POCKET_VALIDATION_MAPPING_PARTIAL")
            elif pocket_group[["rscc", "rsr", "mean_occupancy"]].isna().any().any():
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                reasons.append("POCKET_QUALITY_METRIC_MISSING")
            elif pocket_group["unresolved_heavy_atom_count"].isna().any():
                terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                reasons.append("POCKET_COMPLETENESS_UNAVAILABLE")
            else:
                pocket_metrics = {
                    "mean_rscc": float(pocket_group["rscc"].mean()),
                    "min_rscc": float(pocket_group["rscc"].min()),
                    "mean_rsr": float(pocket_group["rsr"].mean()),
                    "max_rsr": float(pocket_group["rsr"].max()),
                    "mean_occupancy": float(pocket_group["mean_occupancy"].mean()),
                    "min_occupancy": float(pocket_group["mean_occupancy"].min()),
                    "rsrz_outlier_count": int((pocket_group["rsrz"].fillna(-math.inf) > 2.0).sum()),
                    "unresolved_heavy_atom_count": int(pocket_group["unresolved_heavy_atom_count"].sum()),
                    "residue_count": len(pocket_group),
                }
                pocket_pass = (
                    pocket_metrics["mean_rscc"] >= 0.80
                    and pocket_metrics["mean_rsr"] <= 0.30
                    and pocket_metrics["mean_occupancy"] >= 0.80
                    and pocket_metrics["unresolved_heavy_atom_count"] == 0
                )
                if not pocket_pass:
                    reasons.append("POCKET_QUALITY_HARD_GATE_FAILED")

        supported_chain_count = 0
        soft_binding_count = 0
        hard_binding_count = 0
        unavailable_binding_count = 0
        important_gap_warning = False
        if not terminal:
            if chain_group.empty:
                terminal = "FILTER3_TECHNICAL_FAILURE"
                reasons.append("CHAIN_QUALITY_ROWS_MISSING")
            else:
                supported_chain_count = int((chain_group["chain_quality_status"] == "QUALITY_SUPPORTED_CHAIN").sum())
                soft_binding_count = int(chain_group["quality_soft_fail_binding_residue_count"].sum())
                hard_binding_count = int(chain_group["quality_hard_fail_binding_residue_count"].sum())
                unavailable_binding_count = int(chain_group["quality_unavailable_binding_residue_count"].sum())
                important_gap_warning = bool(chain_group["adjacent_pocket_gap_warning"].any())
                if unavailable_binding_count:
                    terminal = "FILTER3_VALIDATION_DATA_UNAVAILABLE"
                    mapping_status = "BINDING_MAPPING_PARTIAL"
                    reasons.append("BINDING_VALIDATION_MAPPING_PARTIAL")
                else:
                    chain_pass = supported_chain_count >= 1 and not bool(chain_group["critical_interface_gap"].any())
                    if not chain_pass:
                        reasons.append("DIRECT_BINDING_QUALITY_SUPPORT_LOST")

        geometry_warning_count = 0
        geometry_fatal = False
        if ligand_row is not None:
            geometry_warning_count = sum(bool_value(ligand_row.get(name)) for name in (
                "geometry_outlier", "chirality_outlier", "clash_outlier"
            ))
            geometry_fatal = (
                clean(ligand_row.get("mapping_status")) != "COMPLETE"
                or clean(ligand_row.get("topology_status")) != "TOPOLOGY_COMPLETE"
                or not bool_value(ligand_row.get("rdkit_parse_success"))
                or not bool_value(ligand_row.get("rdkit_sanitize_success"))
            )
            if geometry_warning_count:
                warnings.append("LIGAND_GEOMETRY_WARNING")

        if not terminal:
            if geometry_fatal:
                reasons.append("LIGAND_RAW_GEOMETRY_FATAL")
            hard_pass = entry_pass and ligand_pass and pocket_pass and chain_pass and not geometry_fatal
            if not hard_pass:
                terminal = "FILTER3_REJECT"
            else:
                high = (
                    finite(entry.get("resolution")) <= 2.5
                    and soft_binding_count == 0
                    and not important_gap_warning
                    and geometry_warning_count == 0
                )
                terminal = "FILTER3_HIGH_QUALITY" if high else "FILTER3_GOOD_QUALITY"
                if not high:
                    warnings.append("RETAINED_WITH_NONCRITICAL_WARNING")

        pair_rows.append({
            **pair,
            "experimental_method_raw": meta.get("experimental_method_raw", ""),
            "experimental_method_class": method,
            "validation_mapping_status": mapping_status,
            "entry_resolution": None if entry is None else entry.get("resolution"),
            "entry_r_work": None if entry is None else entry.get("r_work"),
            "entry_r_free": None if entry is None else entry.get("r_free"),
            "entry_r_free_minus_r_work": None if entry is None else entry.get("r_free_minus_r_work"),
            "ligand_rscc": None if ligand_row is None else ligand_row.get("rscc"),
            "ligand_rsr": None if ligand_row is None else ligand_row.get("rsr"),
            "ligand_mean_occupancy": None if ligand_row is None else ligand_row.get("mean_occupancy"),
            "ligand_unresolved_heavy_atom_count": None if ligand_row is None else ligand_row.get("missing_heavy_atom_count"),
            "pocket_residue_count": pocket_metrics["residue_count"],
            "pocket_mean_rscc": pocket_metrics["mean_rscc"],
            "pocket_min_rscc": pocket_metrics["min_rscc"],
            "pocket_mean_rsr": pocket_metrics["mean_rsr"],
            "pocket_max_rsr": pocket_metrics["max_rsr"],
            "pocket_mean_occupancy": pocket_metrics["mean_occupancy"],
            "pocket_min_occupancy": pocket_metrics["min_occupancy"],
            "pocket_rsrz_outlier_count": pocket_metrics["rsrz_outlier_count"],
            "pocket_unresolved_heavy_atom_count": pocket_metrics["unresolved_heavy_atom_count"],
            "quality_supported_chain_count": supported_chain_count,
            "binding_soft_fail_residue_count": soft_binding_count,
            "binding_hard_fail_residue_count": hard_binding_count,
            "binding_unavailable_residue_count": unavailable_binding_count,
            "critical_interface_gap": bool(chain_group["critical_interface_gap"].any()) if not chain_group.empty else False,
            "important_pocket_gap_warning": important_gap_warning,
            "ligand_geometry_warning_count": geometry_warning_count,
            "posebusters_status": "PENDING_SEPARATE_RAW_GEOMETRY_PASS",
            "entry_hard_gate_pass": entry_pass,
            "ligand_hard_gate_pass": ligand_pass,
            "pocket_hard_gate_pass": pocket_pass,
            "chain_support_gate_pass": chain_pass,
            "terminal_status_pre_posebusters": terminal,
            "reason_codes": ";".join(sorted(set(reasons))),
            "warning_codes": ";".join(sorted(set(warnings))),
            "rule_version": "filter3_quality_v1.0.0",
        })

    pair_results = pd.DataFrame(pair_rows)
    write_parquet(pair_results, bucket_out / "pair_quality_pre_posebusters.parquet")
    write_parquet(ligand_quality, bucket_out / "ligand_validation_mapping.parquet")
    write_parquet(binding_quality, bucket_out / "binding_residue_quality.parquet")
    write_parquet(pocket_quality, bucket_out / "pocket_residue_quality.parquet")
    write_parquet(chain_quality, bucket_out / "chain_quality_support.parquet")
    write_parquet(pd.DataFrame(metadata_rows), bucket_out / "entry_structure_metadata.parquet")
    if gap_details:
        write_parquet(pd.DataFrame(gap_details), bucket_out / "structural_gap_audit.parquet")

    result = {
        "status": "COMPLETED",
        "bucket_id": bucket,
        "pair_count": len(pair_results),
        "ligand_mapping_count": len(ligand_quality),
        "binding_residue_count": len(binding_quality),
        "pocket_residue_count": len(pocket_quality),
        "chain_quality_count": len(chain_quality),
        "terminal_status_pre_posebusters": dict(Counter(pair_results["terminal_status_pre_posebusters"])),
        "mapping_status": dict(Counter(pair_results["validation_mapping_status"])),
        "runtime_seconds": time.time() - started,
        "finished_at": utc(),
    }
    atomic_json(marker, result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bucket", type=int)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=255)
    args = parser.parse_args()
    buckets = [args.bucket] if args.bucket is not None else list(range(args.start, args.end + 1))
    started = time.time()
    totals = Counter()
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        future_map = {pool.submit(process_bucket, bucket): bucket for bucket in buckets}
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            totals.update({"pair_count": result["pair_count"], "binding_residue_count": result["binding_residue_count"], "pocket_residue_count": result["pocket_residue_count"]})
            for key, value in result["terminal_status_pre_posebusters"].items():
                totals[f"terminal::{key}"] += value
            progress = {
                "status": "RUNNING",
                "phase": "QUALITY_MAPPING_AND_PRECLASSIFICATION",
                "bucket_completed": len(results),
                "bucket_total": len(buckets),
                **dict(totals),
                "runtime_seconds": time.time() - started,
                "updated_at": utc(),
            }
            atomic_json(RUN / "status.json", progress)
            print(json.dumps(progress), flush=True)
    progress["status"] = "COMPLETED"
    progress["phase"] = "QUALITY_MAPPING_AND_PRECLASSIFICATION_COMPLETE"
    progress["finished_at"] = utc()
    atomic_json(RUN / "status.json", progress)


if __name__ == "__main__":
    main()
