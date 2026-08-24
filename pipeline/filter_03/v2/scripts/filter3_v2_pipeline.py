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
import yaml


STANDARD_HEAVY_ATOMS = {
    "ALA": "N CA C O CB", "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2", "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG", "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2", "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2", "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2", "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE", "MSE": "N CA C O CB CG SE CE",
    "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "N CA C O CB CG CD",
    "SER": "N CA C O CB OG", "THR": "N CA C O CB OG1 CG2",
    "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "N CA C O CB CG1 CG2",
    "SEC": "N CA C O CB SE",
}
STANDARD_HEAVY_ATOMS = {key: set(value.split()) for key, value in STANDARD_HEAVY_ATOMS.items()}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", ".", "?", "none", "false", "nan"} else text


def finite(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def bool_value(value) -> bool:
    return False if value is None or pd.isna(value) else bool(value)


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), tmp, compression="zstd")
    os.replace(tmp, path)


def read_bucket(root: Path, bucket: int, columns=None) -> pd.DataFrame:
    paths = [str(path) for path in sorted((root / f"bucket_id={bucket:03d}").glob("*.parquet"))]
    if not paths:
        return pd.DataFrame()
    return ds.dataset(paths, format="parquet").to_table(columns=columns).to_pandas(split_blocks=True, self_destruct=True)


def residue_parts(value: str) -> tuple[str, str, str, str]:
    parts = str(value).split("|") + [""] * 4
    return clean(parts[0]), clean(parts[1]), clean(parts[2]), clean(parts[3]).upper()


def numeric_position(label: str, auth: str):
    for value in (label, auth):
        try:
            return int(value)
        except (TypeError, ValueError):
            pass
    return None


def category_rows(block, prefix: str) -> list[dict]:
    category = block.get_mmcif_category(prefix)
    if not category:
        return []
    names = list(category)
    count = max((len(category[name]) for name in names), default=0)
    return [{name: clean(category[name][i]) if i < len(category[name]) else "" for name in names} for i in range(count)]


def parse_unobserved(path: Path, pdb_id: str) -> tuple[list[dict], str]:
    try:
        block = gemmi.cif.read_file(str(path)).sole_block()
        rows = []
        for row in category_rows(block, "_pdbx_unobs_or_zero_occ_residues."):
            rows.append({
                "pdb_id": pdb_id, "model_id": row.get("PDB_model_num", "") or "1",
                "auth_asym_id": row.get("auth_asym_id", ""), "label_asym_id": row.get("label_asym_id", ""),
                "auth_seq_id": row.get("auth_seq_id", ""), "label_seq_id": row.get("label_seq_id", ""),
                "component_id": row.get("auth_comp_id", row.get("label_comp_id", "")),
            })
        return rows, "PARSE_SUCCESS"
    except Exception as exc:
        return [], f"PARSE_FAILED:{type(exc).__name__}:{exc}"[:1000]


def observed_atom_map(atoms: pd.DataFrame, needed: set[tuple[str, str]]) -> dict:
    if atoms.empty or not needed:
        return {}
    needed_by_chain = defaultdict(set)
    for chain_id, residue_id in needed:
        needed_by_chain[chain_id].add(residue_parts(residue_id))
    result = defaultdict(set)
    for row in atoms.itertuples(index=False):
        chain = clean(row.filter_1_chain_instance_id)
        identity = (clean(row.label_seq_id), clean(row.auth_seq_id), clean(row.insertion_code), clean(row.label_comp_id).upper())
        if identity in needed_by_chain.get(chain, set()) and clean(row.type_symbol).upper() != "H":
            result[(chain, identity)].add(clean(row.label_atom_id).upper())
    return result


def pocket_completeness_rows(pocket: pd.DataFrame, binding: pd.DataFrame, atoms: pd.DataFrame, backbone: set[str]):
    binding_keys = set(zip(binding.get("ligand_assembly_placement_id", []), binding.get("chain_instance_id", []), binding.get("protein_residue_id", [])))
    needed = set(zip(pocket["chain_instance_id"], pocket["protein_residue_id"]))
    observed = observed_atom_map(atoms, needed)
    rows = []
    for row in pocket.to_dict("records"):
        identity = residue_parts(row["protein_residue_id"])
        component = identity[3]
        found = observed.get((row["chain_instance_id"], identity), set())
        expected = STANDARD_HEAVY_ATOMS.get(component)
        old_missing = finite(row.get("unresolved_heavy_atom_count"))
        backbone_missing = sorted(backbone - found)
        unknown = expected is None and old_missing not in (0, 0.0)
        side_missing = [] if expected is None else sorted((expected - backbone) - found)
        key = (row["ligand_assembly_placement_id"], row["chain_instance_id"], row["protein_residue_id"])
        rows.append({
            "pair_id": row["pair_id"], "ligand_assembly_placement_id": row["ligand_assembly_placement_id"],
            "chain_instance_id": row["chain_instance_id"], "protein_residue_id": row["protein_residue_id"],
            "component_id": component, "is_direct_binding_residue": key in binding_keys,
            "observed_heavy_atom_ids": ";".join(sorted(found)),
            "missing_backbone_heavy_atom_ids": ";".join(backbone_missing),
            "missing_backbone_heavy_atom_count": len(backbone_missing),
            "missing_sidechain_heavy_atom_ids": ";".join(side_missing),
            "missing_sidechain_heavy_atom_count": len(side_missing),
            "component_atom_template_status": "KNOWN_STANDARD" if expected is not None else "UNKNOWN_COMPONENT_TEMPLATE",
            "completeness_determinable": not unknown,
        })
    return pd.DataFrame(rows)


def build_gap_audit(pocket: pd.DataFrame, gap_candidates: pd.DataFrame, mmcif_root: Path) -> tuple[pd.DataFrame, set[str]]:
    if gap_candidates.empty:
        return pd.DataFrame(), set()
    candidate_pdbs = sorted(set(gap_candidates["pdb_id"].astype(str).str.lower()))
    unobserved = defaultdict(list)
    failed = set()
    for pdb_id in candidate_pdbs:
        rows, status = parse_unobserved(mmcif_root / pdb_id[1:3] / f"{pdb_id}.cif.gz", pdb_id)
        if status != "PARSE_SUCCESS":
            failed.add(pdb_id)
        for row in rows:
            key = (pdb_id, clean(row["model_id"]) or "1", clean(row["auth_asym_id"]), clean(row["label_asym_id"]))
            position = numeric_position(row["label_seq_id"], row["auth_seq_id"])
            if position is not None:
                unobserved[key].append((position, row["component_id"]))
    output = []
    for (pair_id, chain_id), group in pocket.groupby(["pair_id", "chain_instance_id"], sort=False):
        first = group.iloc[0]
        pdb_id = clean(first["pdb_id"]).lower()
        key = (pdb_id, clean(first["model_id"]) or "1", clean(first["auth_asym_id"]), clean(first["label_asym_id"]))
        missing = sorted(set(position for position, _ in unobserved.get(key, [])))
        if not missing:
            continue
        pocket_positions = set()
        for residue_id in group["protein_residue_id"]:
            label, auth, _, _ = residue_parts(residue_id)
            position = numeric_position(label, auth)
            if position is not None:
                pocket_positions.add(position)
        segments = []
        for position in missing:
            if not segments or position != segments[-1][-1] + 1:
                segments.append([position])
            else:
                segments[-1].append(position)
        for segment in segments:
            left, right = segment[0] - 1, segment[-1] + 1
            left_in, right_in = left in pocket_positions, right in pocket_positions
            if left_in or right_in:
                output.append({
                    "pair_id": pair_id, "chain_instance_id": chain_id, "pdb_id": pdb_id,
                    "segment_start": segment[0], "segment_end": segment[-1],
                    "left_flank_position": left, "right_flank_position": right,
                    "left_flank_in_6A_pocket": left_in, "right_flank_in_6A_pocket": right_in,
                    "critical_pocket_gap": left_in and right_in,
                    "pocket_gap_warning": left_in ^ right_in,
                })
    return pd.DataFrame(output), failed


def classify_terminal(state: dict) -> tuple[str, list[str], list[str]]:
    reasons = list(dict.fromkeys(state.get("reject", [])))
    warnings = list(dict.fromkeys(state.get("warnings", [])))
    if state.get("technical_failure"):
        return "FILTER3_TECHNICAL_FAILURE", [state["technical_failure"]], warnings
    if state.get("method") != "xray":
        return "FILTER3_NON_XRAY_PROTOCOL_PENDING", [f"METHOD_{str(state.get('method')).upper()}_PENDING"], warnings
    unavailable = list(dict.fromkeys(state.get("unavailable", [])))
    if unavailable:
        return "FILTER3_VALIDATION_DATA_UNAVAILABLE", unavailable, warnings
    if reasons:
        return "FILTER3_REJECT", reasons, warnings
    resolution = finite(state.get("resolution"))
    ligand_occupancy = finite(state.get("ligand_occupancy"))
    pocket_occupancy = finite(state.get("pocket_occupancy"))
    if ligand_occupancy is None:
        warnings.append("LIGAND_OCCUPANCY_UNAVAILABLE")
    elif ligand_occupancy < 0.80:
        warnings.append("LIGAND_OCCUPANCY_WARNING")
    if pocket_occupancy is None:
        warnings.append("POCKET_OCCUPANCY_UNAVAILABLE")
    elif pocket_occupancy < 0.80:
        warnings.append("POCKET_OCCUPANCY_WARNING")
    high = (
        resolution is not None and resolution <= 2.5
        and ligand_occupancy is not None and ligand_occupancy >= 0.80
        and pocket_occupancy is not None and pocket_occupancy >= 0.80
        and not warnings and not state.get("posebusters_warning")
    )
    return ("FILTER3_HIGH_QUALITY" if high else "FILTER3_GOOD_QUALITY"), [], sorted(set(warnings))


def process_bucket(bucket: int, config: dict, run_dir: Path, force: bool = False) -> dict:
    started = time.time()
    out = run_dir / "work/preclassification_batches" / f"bucket_id={bucket:03d}"
    marker = out / "_COMPLETE.json"
    if marker.exists() and not force:
        return json.loads(marker.read_text())
    out.mkdir(parents=True, exist_ok=True)
    p3 = Path(config["input"]["processing3_output"])
    p2 = Path(config["input"]["processing2_output"])
    evidence = Path(config["input"]["reusable_evidence_output"])
    mmcif = Path(config["input"]["mmcif_root"])

    pairs = read_bucket(p3 / "provisional_pairs", bucket)
    pair_ev = read_bucket(evidence / "filter3_pair_quality", bucket)
    ligand = read_bucket(evidence / "ligand_validation_mapping", bucket)
    pocket = read_bucket(evidence / "pocket_residue_quality", bucket)
    binding = read_bucket(evidence / "binding_residue_quality", bucket)
    gaps_old = read_bucket(evidence / "structural_gap_audit", bucket)
    pb = read_bucket(evidence / "posebusters_raw_geometry", bucket)
    atoms = read_bucket(p2 / "prepared_receptor_assembly_atoms", bucket, [
        "filter_1_chain_instance_id", "label_seq_id", "auth_seq_id", "insertion_code",
        "label_comp_id", "label_atom_id", "type_symbol",
    ])
    if pairs.empty:
        raise RuntimeError(f"missing P3 pairs for bucket {bucket}")
    if len(pair_ev) != len(pairs) or pair_ev["pair_id"].nunique() != len(pairs):
        raise RuntimeError(f"evidence accounting mismatch bucket={bucket} pairs={len(pairs)} evidence={len(pair_ev)}")

    backbone = set(config["pocket_completeness"]["backbone_atoms"])
    complete = pocket_completeness_rows(pocket, binding, atoms, backbone)
    gap_audit, gap_parse_failed = build_gap_audit(pocket, gaps_old, mmcif)
    gap_by_pair = defaultdict(lambda: {"critical": False, "warning": False})
    for row in gap_audit.to_dict("records"):
        gap_by_pair[row["pair_id"]]["critical"] |= bool(row["critical_pocket_gap"])
        gap_by_pair[row["pair_id"]]["warning"] |= bool(row["pocket_gap_warning"])

    pair_map = {row["pair_id"]: row for row in pair_ev.to_dict("records")}
    ligand_map = {row["ligand_assembly_placement_id"]: row for row in ligand.to_dict("records")}
    pb_map = {row["source_ligand_instance_id"]: row for row in pb.to_dict("records")}
    pocket_groups = {key: frame for key, frame in pocket.groupby("pair_id", sort=False)}
    complete_groups = {key: frame for key, frame in complete.groupby("pair_id", sort=False)}
    binding_groups = {key: frame for key, frame in binding.groupby("ligand_assembly_placement_id", sort=False)}

    pair_rows, binding_rows, pending_rows = [], [], []
    thresholds = config
    for pair in pairs.to_dict("records"):
        pair_id = pair["pair_id"]
        placement = pair["ligand_assembly_placement_id"]
        ev = pair_map[pair_id]
        lig = ligand_map.get(placement)
        pgroup = pocket_groups.get(pair_id, pd.DataFrame())
        cgroup = complete_groups.get(pair_id, pd.DataFrame())
        bgroup = binding_groups.get(placement, pd.DataFrame())
        state = {"method": clean(ev.get("experimental_method_class")) or "unknown", "technical_failure": False,
                 "unavailable": [], "reject": [], "warnings": [], "resolution": ev.get("entry_resolution"),
                 "ligand_occupancy": None if lig is None else lig.get("mean_occupancy"),
                 "pocket_occupancy": None, "posebusters_warning": False}
        if not clean(pair_id) or lig is None or pgroup.empty or cgroup.empty or bgroup.empty:
            state["technical_failure"] = "REQUIRED_FROZEN_CONTEXT_MISSING"
        if state["method"] == "xray" and not state["technical_failure"]:
            entry_metrics = [finite(ev.get(name)) for name in ("entry_resolution", "entry_r_work", "entry_r_free", "entry_r_free_minus_r_work")]
            if any(value is None for value in entry_metrics):
                state["unavailable"].append("ENTRY_QUALITY_METRIC_UNAVAILABLE")
            else:
                resolution, r_work, r_free, r_gap = entry_metrics
                if (resolution > thresholds["entry_thresholds"]["resolution_max_angstrom"]
                    or r_work > thresholds["entry_thresholds"]["r_work_max"]
                    or r_free > thresholds["entry_thresholds"]["r_free_max"]
                    or abs(r_gap) > thresholds["entry_thresholds"]["r_free_minus_r_work_abs_max"]):
                    state["reject"].append("ENTRY_QUALITY_HARD_FAIL")

            if not clean(lig.get("validation_mapping_status")).startswith("MAPPED"):
                state["unavailable"].append("LIGAND_VALIDATION_MAPPING_UNAVAILABLE")
            lrscc, lrsr = finite(lig.get("rscc")), finite(lig.get("rsr"))
            if lrscc is None or lrsr is None:
                state["unavailable"].append("LIGAND_DENSITY_METRIC_UNAVAILABLE")
            elif lrscc < thresholds["density_thresholds"]["ligand_rscc_min"] or lrsr > thresholds["density_thresholds"]["ligand_rsr_max"]:
                state["reject"].append("LIGAND_DENSITY_QUALITY_FAIL")

            mapped = pgroup["validation_mapping_status"].astype(str).str.startswith("MAPPED")
            if not bool(mapped.all()):
                state["unavailable"].append("POCKET_VALIDATION_MAPPING_UNAVAILABLE")
            prscc = pd.to_numeric(pgroup["rscc"], errors="coerce")
            prsr = pd.to_numeric(pgroup["rsr"], errors="coerce")
            if prscc.isna().any() or prsr.isna().any():
                state["unavailable"].append("POCKET_DENSITY_METRIC_UNAVAILABLE")
            else:
                if prscc.mean() < thresholds["density_thresholds"]["pocket_mean_rscc_min"] or prsr.mean() > thresholds["density_thresholds"]["pocket_mean_rsr_max"]:
                    state["reject"].append("POCKET_DENSITY_QUALITY_FAIL")
            occupancies = pd.to_numeric(pgroup["mean_occupancy"], errors="coerce")
            state["pocket_occupancy"] = None if occupancies.isna().any() else float(occupancies.mean())

            if not bool(cgroup["completeness_determinable"].all()):
                state["unavailable"].append("POCKET_COMPLETENESS_UNAVAILABLE")
            if int(cgroup["missing_backbone_heavy_atom_count"].sum()) > 0:
                state["reject"].append("POCKET_BACKBONE_INCOMPLETE")
            direct = cgroup[cgroup["is_direct_binding_residue"]]
            if int(direct["missing_sidechain_heavy_atom_count"].sum()) > 0:
                state["reject"].append("DIRECT_BINDING_SIDECHAIN_INCOMPLETE")
            nonbinding = cgroup[~cgroup["is_direct_binding_residue"]]
            nb_total = int(nonbinding["missing_sidechain_heavy_atom_count"].sum())
            nb_max = int(nonbinding["missing_sidechain_heavy_atom_count"].max()) if not nonbinding.empty else 0
            if nb_total >= thresholds["pocket_completeness"]["nonbinding_total_missing_sidechain_reject_min"] or nb_max >= thresholds["pocket_completeness"]["nonbinding_per_residue_missing_sidechain_reject_min"]:
                state["reject"].append("NONBINDING_POCKET_SIDECHAIN_INCOMPLETE")
            elif nb_total:
                state["warnings"].append("NONBINDING_POCKET_SIDECHAIN_WARNING")

            gap = gap_by_pair[pair_id]
            if clean(pair["pdb_id"]).lower() in gap_parse_failed:
                state["unavailable"].append("POCKET_GAP_AUDIT_UNAVAILABLE")
            elif gap["critical"]:
                state["reject"].append("CRITICAL_POCKET_GAP")
            elif gap["warning"]:
                state["warnings"].append("POCKET_GAP_WARNING")

            complete_by_residue = {(r["chain_instance_id"], r["protein_residue_id"]): r for r in cgroup.to_dict("records")}
            supported_by_chain = Counter()
            binding_unavailable = False
            binding_warning = False
            for residue in bgroup.to_dict("records"):
                comp = complete_by_residue.get((residue["chain_instance_id"], residue["protein_residue_id"]))
                mapped_ok = clean(residue.get("validation_mapping_status")).startswith("MAPPED")
                rscc, rsr = finite(residue.get("rscc")), finite(residue.get("rsr"))
                completeness_ok = bool(comp and comp["completeness_determinable"]
                                       and comp["missing_backbone_heavy_atom_count"] == 0
                                       and comp["missing_sidechain_heavy_atom_count"] == 0)
                supported = mapped_ok and rscc is not None and rsr is not None and completeness_ok and rscc >= thresholds["binding_support"]["residue_rscc_min"] and rsr <= thresholds["binding_support"]["residue_rsr_max"]
                if not mapped_ok or rscc is None or rsr is None or comp is None or not comp["completeness_determinable"]:
                    binding_unavailable = True
                    status = "BINDING_RESIDUE_QUALITY_UNAVAILABLE"
                elif supported:
                    supported_by_chain[residue["chain_instance_id"]] += 1
                    status = "QUALITY_SUPPORTED_BINDING_RESIDUE"
                else:
                    binding_warning = True
                    status = "BINDING_RESIDUE_QUALITY_WARNING"
                binding_rows.append({**residue, "pair_id": pair_id, "residue_structurally_complete": completeness_ok,
                                     "quality_supported": supported, "filter3_v2_binding_quality_status": status})
            if binding_unavailable:
                state["unavailable"].append("BINDING_RESIDUE_VALIDATION_UNAVAILABLE")
            supported_chains = sum(count >= thresholds["binding_support"]["supported_residues_per_chain_min"] for count in supported_by_chain.values())
            if supported_chains < thresholds["binding_support"]["supported_chains_per_pair_min"]:
                state["reject"].append("CHAIN_QUALITY_SUPPORT_LOST")
            if binding_warning:
                state["warnings"].append("BINDING_RESIDUE_QUALITY_WARNING")

            geometry_warning = any(bool_value(lig.get(name)) for name in ("geometry_outlier", "chirality_outlier", "clash_outlier"))
            if geometry_warning:
                state["warnings"].append("LIGAND_VALIDATION_GEOMETRY_WARNING")

            source_id = clean(lig.get("source_ligand_instance_id"))
            pbrow = pb_map.get(source_id)
            pre_status, pre_reasons, pre_warnings = classify_terminal(state)
            if pre_status in {"FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY"}:
                if pbrow is None:
                    final_status = "FILTER3_POSEBUSTERS_PENDING"
                    pending_rows.append({"bucket_id": bucket, "source_ligand_instance_id": source_id,
                                         "ligand_assembly_placement_id": placement, "pair_id": pair_id})
                    final_reasons, final_warnings = ["POSEBUSTERS_RESULT_PENDING"], pre_warnings
                elif clean(pbrow.get("posebusters_status")) != "COMPLETED":
                    final_status, final_reasons, final_warnings = "FILTER3_TECHNICAL_FAILURE", ["POSEBUSTERS_EXECUTION_FAILURE"], pre_warnings
                else:
                    fatal_chem = any(not bool_value(pbrow.get(name)) for name in ("sanitization", "all_atoms_connected", "no_radicals"))
                    fatal_clash = not bool_value(pbrow.get("internal_steric_clash"))
                    if fatal_chem: state["reject"].append("POSEBUSTERS_FATAL_CHEMISTRY")
                    if fatal_clash: state["reject"].append("POSEBUSTERS_FATAL_INTERNAL_CLASH")
                    pbwarning = any(not bool_value(pbrow.get(name)) for name in thresholds["posebusters"]["warning_checks"])
                    if pbwarning: state["warnings"].append("POSEBUSTERS_NONFATAL_WARNING")
                    state["posebusters_warning"] = pbwarning
                    final_status, final_reasons, final_warnings = classify_terminal(state)
            else:
                final_status, final_reasons, final_warnings = pre_status, pre_reasons, pre_warnings
        else:
            final_status, final_reasons, final_warnings = classify_terminal(state)

        pair_rows.append({
            **pair, "experimental_method_class": state["method"],
            "entry_resolution": ev.get("entry_resolution"), "entry_r_work": ev.get("entry_r_work"),
            "entry_r_free": ev.get("entry_r_free"), "entry_r_free_minus_r_work": ev.get("entry_r_free_minus_r_work"),
            "ligand_rscc": None if lig is None else lig.get("rscc"), "ligand_rsr": None if lig is None else lig.get("rsr"),
            "ligand_mean_occupancy": state.get("ligand_occupancy"),
            "pocket_mean_rscc": None if pgroup.empty else finite(pd.to_numeric(pgroup["rscc"], errors="coerce").mean()),
            "pocket_mean_rsr": None if pgroup.empty else finite(pd.to_numeric(pgroup["rsr"], errors="coerce").mean()),
            "pocket_mean_occupancy": state.get("pocket_occupancy"),
            "pocket_missing_backbone_heavy_atom_count": None if cgroup.empty else int(cgroup["missing_backbone_heavy_atom_count"].sum()),
            "direct_binding_missing_sidechain_heavy_atom_count": None if cgroup.empty else int(cgroup.loc[cgroup["is_direct_binding_residue"], "missing_sidechain_heavy_atom_count"].sum()),
            "nonbinding_pocket_missing_sidechain_heavy_atom_count": None if cgroup.empty else int(cgroup.loc[~cgroup["is_direct_binding_residue"], "missing_sidechain_heavy_atom_count"].sum()),
            "critical_pocket_gap": gap_by_pair[pair_id]["critical"], "pocket_gap_warning": gap_by_pair[pair_id]["warning"],
            "filter3_v2_terminal_status": final_status, "decision": "PASS" if final_status in {"FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY", "FILTER3_NON_XRAY_PROTOCOL_PENDING"} else ("REVIEW" if final_status in {"FILTER3_VALIDATION_DATA_UNAVAILABLE", "FILTER3_POSEBUSTERS_PENDING"} else ("FAIL" if final_status == "FILTER3_TECHNICAL_FAILURE" else "REJECT")),
            "destination": "ordinary_ground_truth" if final_status in {"FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY"} else ("method_protocol_pending" if final_status == "FILTER3_NON_XRAY_PROTOCOL_PENDING" else ("manual_review" if final_status in {"FILTER3_VALIDATION_DATA_UNAVAILABLE", "FILTER3_POSEBUSTERS_PENDING"} else ("technical_failure" if final_status == "FILTER3_TECHNICAL_FAILURE" else "excluded"))),
            "reason_codes": ";".join(sorted(set(final_reasons))), "warning_codes": ";".join(sorted(set(final_warnings))),
            "rsrz_used_for_membership": False, "ligand_completeness_inherited_from_processing2": True,
            "rule_version": "filter3_v2.0.0", "bucket_id": bucket,
        })

    pair_frame = pd.DataFrame(pair_rows)
    write_parquet(pair_frame, out / "pair_quality_pre_posebusters.parquet")
    write_parquet(complete, out / "pocket_completeness.parquet")
    write_parquet(pd.DataFrame(binding_rows), out / "binding_residue_quality_v2.parquet")
    if not gap_audit.empty: write_parquet(gap_audit, out / "structural_gap_audit_v2.parquet")
    if pending_rows: write_parquet(pd.DataFrame(pending_rows).drop_duplicates("source_ligand_instance_id"), out / "posebusters_pending_sources.parquet")
    result = {"status": "COMPLETED", "bucket_id": bucket, "pair_count": len(pair_frame),
              "terminal_status_counts": dict(Counter(pair_frame["filter3_v2_terminal_status"])),
              "pending_source_count": len({row["source_ligand_instance_id"] for row in pending_rows}),
              "runtime_seconds": time.time() - started, "finished_at": utc()}
    atomic_json(marker, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--bucket", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    run_dir = Path(args.run_dir)
    buckets = [args.bucket] if args.bucket is not None else list(range(config["runtime"]["bucket_start"], config["runtime"]["bucket_end"] + 1))
    workers = args.workers or config["runtime"]["workers"]
    started, results = time.time(), []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(process_bucket, bucket, config, run_dir, args.force): bucket for bucket in buckets}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            counts = Counter()
            for item in results: counts.update(item["terminal_status_counts"])
            progress = {"status": "RUNNING", "phase": "FILTER3_V2_PRECLASSIFICATION", "bucket_completed": len(results),
                        "bucket_total": len(buckets), "pair_count": sum(item["pair_count"] for item in results),
                        "pending_source_count": sum(item["pending_source_count"] for item in results),
                        "terminal_status_counts": dict(counts), "runtime_seconds": time.time() - started, "updated_at": utc()}
            atomic_json(run_dir / "status.json", progress)
            print(json.dumps(progress), flush=True)
    progress["status"] = "COMPLETED"
    progress["phase"] = "FILTER3_V2_PRECLASSIFICATION_COMPLETE"
    atomic_json(run_dir / "status.json", progress)


if __name__ == "__main__":
    main()
