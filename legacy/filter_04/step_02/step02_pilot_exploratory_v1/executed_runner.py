#!/usr/bin/env python3
"""Filter 4 Step 2: selected biological-assembly equivalence.

This stage is deliberately identity-only.  It consumes the frozen Step 1 lattice
instances and the exact, already-expanded assembly contexts consumed by frozen
Processing 2.  It never recomputes lattice proximity or physical contacts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import yaml
import gemmi


VERSION = "filter4_step2_v1.0.0"
SCHEMA_VERSION = "filter4_step2_schema_v1.0.0"
R_COLS = [f"R{i}{j}" for i in range(1, 4) for j in range(1, 4)]
T_COLS = ["tx", "ty", "tz"]
BA_R_COLS = [f"r{i}{j}" for i in range(1, 4) for j in range(1, 4)]
BA_T_COLS = ["t1", "t2", "t3"]
BLANKS = {"", ".", "?", "None", "False", "nan", "<NA>"}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in BLANKS else text


def sha256(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip", na_rep="")
    os.replace(tmp, path)


def parse_asym_ids(value: Any) -> tuple[str, ...]:
    text = clean(value).strip(";").replace("\n", "")
    return tuple(x.strip() for x in text.split(",") if x.strip())


def stable_rank(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def matrix_key(values: Iterable[float], digits: int = 9) -> tuple[float, ...]:
    return tuple(round(float(x), digits) for x in values)


def load_config(path: Path) -> dict:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    for key in ("step1_run", "processing2_run", "assembly_context", "mmcif_root"):
        if key not in cfg["input"]:
            raise RuntimeError(f"missing input.{key}")
    return cfg


def read_inputs(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    step1 = Path(cfg["input"]["step1_run"]) / "output"
    pairs = pd.read_csv(step1 / "01_pair_step1_inventory.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    neighbors = pd.read_csv(step1 / "02_crystal_neighbor_instances.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    sources = pd.read_csv(step1 / "03_source_object_inventory.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    for col in R_COLS + T_COLS + ["min_distance_to_ligand", "min_distance_to_pocket"]:
        if col in neighbors:
            neighbors[col] = pd.to_numeric(neighbors[col], errors="coerce")
    for col in ("touches_ligand_6A", "touches_pocket_6A"):
        neighbors[col] = neighbors[col].str.lower().eq("true")
    return pairs, neighbors, sources


def load_context(cfg: dict) -> pd.DataFrame:
    cols = ["pdb_id", "assembly_id", "assembly_gen_row_id", "oper_expression_raw",
            "asym_id_list_raw", "operator_path"] + BA_R_COLS + BA_T_COLS
    table = ds.dataset(cfg["input"]["assembly_context"], format="parquet").to_table(columns=cols)
    frame = table.to_pandas(split_blocks=True, self_destruct=True)
    frame["pdb_id"] = frame["pdb_id"].astype(str).str.lower()
    frame["assembly_id"] = frame["assembly_id"].astype(str)
    frame["asym_ids"] = frame["asym_id_list_raw"].map(parse_asym_ids)
    frame["operator_is_composed"] = frame["operator_path"].astype(str).str.contains("*", regex=False)
    frame["operator_is_identity"] = (
        np.max(np.abs(frame[BA_R_COLS].to_numpy(float).reshape(-1, 3, 3) - np.eye(3)), axis=(1, 2)) <= 1e-10
    ) & (np.max(np.abs(frame[BA_T_COLS].to_numpy(float)), axis=1) <= 1e-10)
    return frame


def choose_pilot(pairs: pd.DataFrame, neighbors: pd.DataFrame, sources: pd.DataFrame,
                 context: pd.DataFrame, count: int) -> pd.DataFrame:
    n = neighbors.copy()
    n["pdb_id"] = n["pdb_id"].str.lower()
    flags: dict[str, set[str]] = defaultdict(set)
    for pdb, group in n.groupby("pdb_id", sort=False):
        sym = pd.to_numeric(group["symmetry_operation_id"], errors="coerce").fillna(-1)
        hkl0 = (pd.to_numeric(group["cell_h"], errors="coerce").fillna(99).eq(0) &
                pd.to_numeric(group["cell_k"], errors="coerce").fillna(99).eq(0) &
                pd.to_numeric(group["cell_l"], errors="coerce").fillna(99).eq(0))
        if sym.nunique() > 1: flags[pdb].add("multiple_crystal_symops")
        if (sym.eq(0) & hkl0).any(): flags[pdb].add("same_asu_neighbor")
        if group["source_object_type"].eq("NONPOLYMER").any(): flags[pdb].add("nonpoly_neighbor")
        ligand = group["touches_ligand_6A"]; pocket = group["touches_pocket_6A"]
        if (ligand & ~pocket).any(): flags[pdb].add("ligand_only")
        if (~ligand & pocket).any(): flags[pdb].add("pocket_only")
        if (ligand & pocket).any(): flags[pdb].add("ligand_and_pocket")
    for (pdb, _aid), group in context.groupby(["pdb_id", "assembly_id"], sort=False):
        if group["operator_is_composed"].any(): flags[pdb].add("composed_ba_operator")
        if (~group["operator_is_identity"]).any(): flags[pdb].add("nonidentity_ba_operator")
        if group["operator_is_identity"].all(): flags[pdb].add("identity_only_ba")
        if len(group) == 1: flags[pdb].add("monomeric_context")
        if len(group) > 1: flags[pdb].add("multimeric_context")
    source_count = sources.groupby(["pdb_id", "label_asym_id"]).size()
    for (pdb, _), value in source_count.items():
        if value > 1: flags[str(pdb).lower()].add("multiple_source_objects_per_asym")

    available = sorted(set(n["pdb_id"]))
    strata = sorted({flag for values in flags.values() for flag in values})
    selected: list[str] = []
    reasons: dict[str, set[str]] = defaultdict(set)
    # Deterministic set cover with several representatives per critical stratum.
    target_per = max(3, min(12, count // max(1, len(strata))))
    for stratum in strata:
        candidates = sorted((p for p in available if stratum in flags[p]), key=lambda p: stable_rank(stratum + "|" + p))
        for pdb in candidates[:target_per]:
            if len(selected) >= count: break
            if pdb not in selected: selected.append(pdb)
            reasons[pdb].add(stratum)
    for pdb in sorted(available, key=lambda p: stable_rank("pilot_fill|" + p)):
        if len(selected) >= count: break
        if pdb not in selected:
            selected.append(pdb); reasons[pdb].add("deterministic_fill")
    rows = []
    for pdb in sorted(selected):
        rows.append({"pdb_id": pdb, "selection_reason": ";".join(sorted(reasons[pdb])),
                     "all_strata": ";".join(sorted(flags[pdb]))})
    return pd.DataFrame(rows)


def build_ba_inventory(sources: pd.DataFrame, pairs: pd.DataFrame, context: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    active = pairs[["pdb_id"]].drop_duplicates().copy()
    # assembly/model are encoded in candidate_pair_id but not separate Step 1 pair columns.
    pair_map = []
    for pid in pairs["candidate_pair_id"]:
        parts = str(pid).split("|")
        pair_map.append((str(pid), parts[2], parts[3]))
    pair_keys = pd.DataFrame(pair_map, columns=["candidate_pair_id", "assembly_id", "model_id"])
    pkeys = pairs[["candidate_pair_id", "pdb_id"]].merge(pair_keys, on="candidate_pair_id")
    needed = set(map(tuple, pkeys[["pdb_id", "assembly_id", "model_id"]].drop_duplicates().itertuples(index=False, name=None)))
    context_map: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in context.to_dict("records"):
        context_map[(clean(row["pdb_id"]).lower(), clean(row["assembly_id"]))].append(row)
    source_map: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in sources.to_dict("records"):
        source_map[(clean(row["pdb_id"]).lower(), clean(row["model_id"]) or "1")].append(row)

    raw_rows = []
    missing_context = []
    for pdb, aid, model in sorted(needed):
        contexts = context_map.get((pdb.lower(), aid), [])
        if not contexts:
            missing_context.append((pdb, aid, model)); continue
        for source in source_map.get((pdb.lower(), model), []):
            label = clean(source["label_asym_id"])
            for ctx in contexts:
                if label not in ctx["asym_ids"]:
                    continue
                rt = [ctx[c] for c in BA_R_COLS + BA_T_COLS]
                provenance = f'{ctx["assembly_gen_row_id"]}|{ctx["operator_path"]}'
                raw_rows.append({
                    "pdb_id": pdb, "assembly_id": aid, "model_id": model,
                    "source_object_id": source["source_object_id"], "source_object_type": source["source_object_type"],
                    "source_entity_id": source["entity_id"], "source_label_asym_id": source["label_asym_id"],
                    "source_auth_asym_id": source["auth_asym_id"], "source_comp_id": source["comp_id"],
                    "source_residue_id": source["residue_id"], "assembly_gen_row_id": ctx["assembly_gen_row_id"],
                    "assembly_operator_expression": ctx["oper_expression_raw"], "operator_path": ctx["operator_path"],
                    "operator_is_composed": bool(ctx["operator_is_composed"]), "operator_is_identity": bool(ctx["operator_is_identity"]),
                    "provenance": provenance,
                    **{f"R_ba_{i}{j}": float(ctx[f"r{i}{j}"]) for i in range(1,4) for j in range(1,4)},
                    **{f"t_ba_{i}": float(ctx[f"t{i}"]) for i in range(1,4)},
                    "_placement_key": matrix_key(rt),
                })
    raw = pd.DataFrame(raw_rows)
    if raw.empty:
        raise RuntimeError("BA inventory is empty")
    key_cols = ["pdb_id", "assembly_id", "model_id", "source_object_id", "_placement_key"]
    output = []
    for _, group in raw.groupby(key_cols, sort=False, dropna=False):
        first = group.iloc[0].to_dict()
        prov = sorted(set(group["provenance"]))
        first["assembly_gen_row_id"] = ";".join(sorted(set(group["assembly_gen_row_id"].astype(str))))
        first["assembly_operator_expression"] = ";".join(sorted(set(group["assembly_operator_expression"].astype(str))))
        first["operator_path"] = ";".join(sorted(set(group["operator_path"].astype(str))))
        first["operator_is_composed"] = bool(group["operator_is_composed"].any())
        first["operator_is_identity"] = bool(group["operator_is_identity"].all())
        first["provenance_multiplicity"] = len(prov)
        digest = hashlib.sha1("|".join(str(first[c]) for c in key_cols).encode()).hexdigest()[:16]
        first["ba_instance_id"] = f"BAI|{first['pdb_id']}|{first['assembly_id']}|{digest}"
        output.append(first)
    inventory = pd.DataFrame(output).drop(columns=["_placement_key", "provenance"])
    ordered = ["ba_instance_id", "pdb_id", "assembly_id", "model_id", "source_object_id", "source_object_type",
               "source_entity_id", "source_label_asym_id", "source_auth_asym_id", "source_comp_id", "source_residue_id",
               "assembly_gen_row_id", "assembly_operator_expression", "operator_path", "operator_is_composed",
               "operator_is_identity", "provenance_multiplicity"] + \
              [f"R_ba_{i}{j}" for i in range(1,4) for j in range(1,4)] + [f"t_ba_{i}" for i in range(1,4)]
    return inventory[ordered], {"missing_context_keys": missing_context, "raw_inventory_rows": len(raw),
                                "deduplicated_inventory_rows": len(inventory),
                                "collapsed_duplicate_provenances": len(raw) - len(inventory)}


def compare_neighbors(neighbors: pd.DataFrame, inventory: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rot_tol = float(cfg["equivalence"]["rotation_max_abs_tolerance"])
    trans_tol = float(cfg["equivalence"]["translation_max_abs_tolerance_angstrom"])
    review_factor = float(cfg["equivalence"]["near_threshold_factor"])
    ba_r_cols = [f"R_ba_{i}{j}" for i in range(1,4) for j in range(1,4)]
    ba_t_cols = [f"t_ba_{i}" for i in range(1,4)]
    index: dict[tuple[str, str, str, str], dict] = {}
    for key, group in inventory.groupby(["pdb_id", "assembly_id", "model_id", "source_object_id"], sort=False):
        index[tuple(map(str, key))] = {
            "ids": group["ba_instance_id"].to_numpy(), "r": group[ba_r_cols].to_numpy(float).reshape(-1,3,3),
            "t": group[ba_t_cols].to_numpy(float), "composed": group["operator_is_composed"].to_numpy(bool),
            "identity": group["operator_is_identity"].to_numpy(bool), "paths": group["operator_path"].to_numpy(),
        }
    results = []
    cache: dict[tuple, dict] = {}
    for row in neighbors.itertuples(index=False):
        key = (str(row.pdb_id).lower(), str(row.assembly_id), str(row.model_id), str(row.source_object_id))
        crystal_values = tuple(float(getattr(row, c)) for c in R_COLS + T_COLS)
        ckey = key + matrix_key(crystal_values, digits=10)
        result = cache.get(ckey)
        if result is None:
            group = index.get(key)
            if group is None:
                result = {"same_source_in_ba": False, "best_ba_instance_id": "", "best_ba_operator_path": "",
                          "best_ba_operator_is_composed": False, "best_ba_operator_is_identity": False,
                          "rotation_max_abs_diff": np.nan, "rotation_frobenius_diff": np.nan,
                          "translation_max_abs_diff": np.nan, "translation_norm_diff": np.nan,
                          "coordinate_rmsd": np.nan, "max_atom_deviation": np.nan,
                          "coordinate_verification_status": "NOT_APPLICABLE_SOURCE_ABSENT",
                          "equivalence_status": "EXTERNAL_CRYSTAL_INSTANCE",
                          "equivalence_reason": "SOURCE_OBJECT_NOT_IN_SELECTED_BA"}
            else:
                cr = np.array(crystal_values[:9]).reshape(3,3); ct = np.array(crystal_values[9:])
                dr = group["r"] - cr; dt = group["t"] - ct
                rmax = np.max(np.abs(dr), axis=(1,2)); rfrob = np.linalg.norm(dr, axis=(1,2))
                tmax = np.max(np.abs(dt), axis=1); tnorm = np.linalg.norm(dt, axis=1)
                score = np.sqrt(np.square(rfrob) + np.square(tnorm))
                best = int(np.argmin(score))
                within = rmax[best] <= rot_tol and tmax[best] <= trans_tol
                near = rmax[best] <= rot_tol * review_factor and tmax[best] <= trans_tol * review_factor
                if within:
                    status, reason, cv = "BA_EQUIVALENT", "FINAL_TRANSFORM_WITHIN_FROZEN_TOLERANCE", "NOT_REQUIRED_TRANSFORM_DECISIVE"
                elif near:
                    status, reason, cv = "BA_EQUIVALENCE_REVIEW", "NEAR_THRESHOLD_COORDINATE_VERIFICATION_REQUIRED", "PENDING"
                else:
                    status, reason, cv = "EXTERNAL_CRYSTAL_INSTANCE", "SAME_SOURCE_DIFFERENT_FINAL_PLACEMENT", "NOT_REQUIRED_TRANSFORM_DECISIVE"
                result = {"same_source_in_ba": True, "best_ba_instance_id": group["ids"][best],
                          "best_ba_operator_path": group["paths"][best],
                          "best_ba_operator_is_composed": bool(group["composed"][best]),
                          "best_ba_operator_is_identity": bool(group["identity"][best]),
                          "rotation_max_abs_diff": float(rmax[best]), "rotation_frobenius_diff": float(rfrob[best]),
                          "translation_max_abs_diff": float(tmax[best]), "translation_norm_diff": float(tnorm[best]),
                          "coordinate_rmsd": np.nan, "max_atom_deviation": np.nan,
                          "coordinate_verification_status": cv, "equivalence_status": status,
                          "equivalence_reason": reason}
            cache[ckey] = result
        results.append(result)
    return pd.concat([neighbors.reset_index(drop=True), pd.DataFrame(results)], axis=1)


def cif_text(value: Any) -> str:
    text = clean(value)
    return "" if text.lower() in {"false", "none", "null"} else text


def source_coordinates(block: gemmi.cif.Block, source_id: str) -> np.ndarray:
    parts = source_id.split("|")
    if parts[0] == "POLYMER" and len(parts) >= 5:
        _, model, entity, label, auth = parts[:5]
        nonpoly = False; comp = auth_seq = label_seq = ins = ""
    elif parts[0] == "NONPOLYMER" and len(parts) >= 9:
        _, model, entity, label, auth, comp, auth_seq, label_seq, ins = parts[:9]
        nonpoly = True
    else:
        return np.empty((0, 3), dtype=float)
    tags = ["_atom_site.pdbx_PDB_model_num", "_atom_site.label_entity_id", "_atom_site.label_asym_id",
            "_atom_site.auth_asym_id", "_atom_site.label_comp_id", "_atom_site.auth_comp_id",
            "_atom_site.auth_seq_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.type_symbol", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
    xyz = []
    try:
        for row in block.find(tags):
            values = [cif_text(row[i]) for i in range(len(tags))]
            rmodel, rentity, rlabel, rauth = values[:4]
            if (rmodel or "1") != (model or "1") or rentity != entity or rlabel != label or rauth != auth:
                continue
            if values[9].upper() in {"H", "D", "T"}: continue
            if nonpoly:
                rcomp = (values[4] or values[5]).upper()
                if rcomp != comp.upper() or values[6] != auth_seq or values[7] != label_seq or values[8] != ins:
                    continue
            xyz.append([float(values[10]), float(values[11]), float(values[12])])
    except Exception:
        return np.empty((0, 3), dtype=float)
    return np.asarray(xyz, dtype=float).reshape(-1, 3)


def coordinate_fallback(eq: pd.DataFrame, inventory: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, dict]:
    pending = eq["coordinate_verification_status"].eq("PENDING")
    meta = {"coordinate_fallback_requested_rows": int(pending.sum()), "coordinate_fallback_unique_cases": 0,
            "coordinate_fallback_resolved_equivalent": 0, "coordinate_fallback_resolved_external": 0,
            "coordinate_fallback_review": 0, "raw_mmcif_read": False}
    if not pending.any(): return eq, meta
    meta["raw_mmcif_read"] = True
    inv = inventory.set_index("ba_instance_id", drop=False)
    root = Path(cfg["input"]["mmcif_root"])
    rmsd_tol = float(cfg["equivalence"]["coordinate_rmsd_tolerance_angstrom"])
    max_tol = float(cfg["equivalence"]["max_atom_deviation_tolerance_angstrom"])
    cache: dict[tuple, tuple] = {}
    blocks: dict[str, gemmi.cif.Block | None] = {}
    for idx in eq.index[pending]:
        row = eq.loc[idx]; ba = inv.loc[row["best_ba_instance_id"]]
        case = (row["pdb_id"], row["source_object_id"], matrix_key(row[R_COLS + T_COLS], 10), row["best_ba_instance_id"])
        result = cache.get(case)
        if result is None:
            pdb = str(row["pdb_id"]).lower()
            if pdb not in blocks:
                path = root / f"{pdb}.cif.gz"
                if not path.exists(): path = root / pdb[1:3] / f"{pdb}.cif.gz"
                try: blocks[pdb] = gemmi.cif.read(str(path)).sole_block()
                except Exception: blocks[pdb] = None
            xyz = source_coordinates(blocks[pdb], str(row["source_object_id"])) if blocks[pdb] is not None else np.empty((0,3))
            if not len(xyz):
                result = (np.nan, np.nan, "BA_EQUIVALENCE_REVIEW", "COORDINATE_FALLBACK_SOURCE_ATOMS_UNAVAILABLE", "FAILED")
            else:
                cr = row[R_COLS].to_numpy(float).reshape(3,3); ct = row[T_COLS].to_numpy(float)
                br = ba[[f"R_ba_{i}{j}" for i in range(1,4) for j in range(1,4)]].to_numpy(float).reshape(3,3)
                bt = ba[[f"t_ba_{i}" for i in range(1,4)]].to_numpy(float)
                delta = (xyz @ cr.T + ct) - (xyz @ br.T + bt)
                dev = np.linalg.norm(delta, axis=1); rmsd = float(np.sqrt(np.mean(np.square(dev)))); maxdev = float(dev.max())
                if rmsd <= rmsd_tol and maxdev <= max_tol:
                    result = (rmsd, maxdev, "BA_EQUIVALENT", "COORDINATE_FALLBACK_WITHIN_FROZEN_TOLERANCE", "PASS_EQUIVALENT")
                else:
                    result = (rmsd, maxdev, "EXTERNAL_CRYSTAL_INSTANCE", "COORDINATE_FALLBACK_DIFFERENT_PLACEMENT", "PASS_EXTERNAL")
            cache[case] = result
        eq.loc[idx, ["coordinate_rmsd", "max_atom_deviation", "equivalence_status", "equivalence_reason",
                     "coordinate_verification_status"]] = result
    meta["coordinate_fallback_unique_cases"] = len(cache)
    meta["coordinate_fallback_resolved_equivalent"] = int(eq["coordinate_verification_status"].eq("PASS_EQUIVALENT").sum())
    meta["coordinate_fallback_resolved_external"] = int(eq["coordinate_verification_status"].eq("PASS_EXTERNAL").sum())
    meta["coordinate_fallback_review"] = int(eq["coordinate_verification_status"].eq("FAILED").sum())
    return eq, meta


def distribution(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["rotation_max_abs_diff", "rotation_frobenius_diff", "translation_max_abs_diff", "translation_norm_diff"]
    quantiles = [0, .001, .01, .05, .25, .5, .75, .95, .99, .999, 1]
    same = frame[frame["same_source_in_ba"]].copy()
    cohorts = [("ALL_SAME_SOURCE", same)] + [(status, group) for status, group in same.groupby("equivalence_status")]
    for cohort, group in cohorts:
        for metric in metrics:
            values = group[metric].dropna().to_numpy(float)
            if not len(values): continue
            qs = np.quantile(values, quantiles)
            for q, value in zip(quantiles, qs):
                rows.append({"cohort": cohort, "metric": metric, "statistic": f"q{q:g}", "value": value, "count": len(values)})
            exact = int(np.sum(values <= 1e-12))
            rows.append({"cohort": cohort, "metric": metric, "statistic": "count_le_1e-12", "value": exact, "count": len(values)})
            for low, high in [(-16,-12),(-12,-10),(-10,-8),(-8,-6),(-6,-4),(-4,-2),(-2,0),(0,2),(2,6)]:
                positive = values[values > 0]
                number = int(np.sum((np.log10(positive) >= low) & (np.log10(positive) < high))) if len(positive) else 0
                rows.append({"cohort": cohort, "metric": metric, "statistic": f"log10_bin_[{low},{high})", "value": number, "count": len(values)})
    return pd.DataFrame(rows)


def pair_inventory(pairs: pd.DataFrame, eq: pd.DataFrame) -> pd.DataFrame:
    aggs = []
    for pid, group in eq.groupby("candidate_pair_id", sort=False):
        status = group["equivalence_status"]
        external = status.eq("EXTERNAL_CRYSTAL_INSTANCE")
        review = status.eq("BA_EQUIVALENCE_REVIEW")
        lig = group["touches_ligand_6A"].astype(bool); pocket = group["touches_pocket_6A"].astype(bool)
        aggs.append({"candidate_pair_id": pid, "step1_neighbor_count_observed": len(group),
                     "ba_equivalent_instance_count": int(status.eq("BA_EQUIVALENT").sum()),
                     "external_crystal_instance_count": int(external.sum()), "review_instance_count": int(review.sum()),
                     "external_ligand_6A_count": int((external & lig).sum()),
                     "external_pocket_6A_count": int((external & pocket).sum()),
                     "external_both_6A_count": int((external & lig & pocket).sum()),
                     "has_external_crystal_neighbor": bool(external.any())})
    agg = pd.DataFrame(aggs)
    out = pairs.merge(agg, on="candidate_pair_id", how="left")
    count_cols = ["step1_neighbor_count_observed", "ba_equivalent_instance_count", "external_crystal_instance_count",
                  "review_instance_count", "external_ligand_6A_count", "external_pocket_6A_count", "external_both_6A_count"]
    out[count_cols] = out[count_cols].fillna(0).astype(np.int64)
    out["has_external_crystal_neighbor"] = out["has_external_crystal_neighbor"].fillna(False).astype(bool)
    status = []
    for row in out.itertuples():
        if row.step1_status == "NO_NEIGHBOR": value = "UPSTREAM_NO_NEIGHBOR"
        elif row.step1_status != "SUCCESS": value = "ERROR"
        elif row.review_instance_count: value = "BA_EQUIVALENCE_REVIEW"
        elif row.external_crystal_instance_count: value = "SUCCESS_EXTERNAL_NEIGHBOR"
        else: value = "SUCCESS_NO_EXTERNAL_NEIGHBOR"
        status.append(value)
    out["step2_status"] = status
    return out


def make_schema(frames: dict[str, pd.DataFrame]) -> dict:
    return {"schema_version": SCHEMA_VERSION, "datasets": {
        name: [{"column_name": c, "data_type": str(frame[c].dtype), "nullable": bool(frame[c].isna().any())}
               for c in frame.columns] for name, frame in frames.items()}}


def validate(mode: str, pairs_in: pd.DataFrame, neighbors_in: pd.DataFrame, inventory: pd.DataFrame,
             eq: pd.DataFrame, pair_out: pd.DataFrame, meta: dict, cfg: dict) -> dict:
    counts = Counter(eq["equivalence_status"])
    classified = sum(counts[x] for x in ("BA_EQUIVALENT", "EXTERNAL_CRYSTAL_INSTANCE", "BA_EQUIVALENCE_REVIEW"))
    duplicate_neighbors = int(eq.duplicated(["candidate_pair_id", "crystal_instance_key"]).sum())
    duplicate_ba = int(inventory["ba_instance_id"].duplicated().sum())
    expected_neighbors = int(cfg["validation"]["expected_full_instances"])
    expected_pairs = int(cfg["validation"]["expected_full_pairs"])
    controls = {
        "identity_ba_equivalent_present": bool(((eq["equivalence_status"] == "BA_EQUIVALENT") & eq["best_ba_operator_is_identity"]).any()),
        "nonidentity_ba_equivalent_present": bool(((eq["equivalence_status"] == "BA_EQUIVALENT") & ~eq["best_ba_operator_is_identity"]).any()),
        "composed_ba_case_present": bool(eq["best_ba_operator_is_composed"].any()),
        "same_source_different_placement_present": bool(eq["equivalence_reason"].eq("SAME_SOURCE_DIFFERENT_FINAL_PLACEMENT").any()),
        "source_not_in_ba_present": bool(eq["equivalence_reason"].eq("SOURCE_OBJECT_NOT_IN_SELECTED_BA").any()),
        "nonpoly_ba_equivalent_present": bool(((eq["equivalence_status"] == "BA_EQUIVALENT") & eq["source_object_type"].eq("NONPOLYMER")).any()),
    }
    p2_run = Path(cfg["input"]["processing2_run"])
    try:
        p2_frozen = json.loads((p2_run / "_FROZEN.json").read_text())
        p2_coordinate = json.loads((p2_run / "audit/assembly_coordinate_validation.json").read_text())
    except Exception:
        p2_frozen = {}; p2_coordinate = {}
    checks = {
        "input_instances_fully_accounted": len(eq) == len(neighbors_in) == classified,
        "input_pairs_fully_accounted": len(pair_out) == len(pairs_in),
        "silent_drop_zero": len(eq) - classified == 0,
        "duplicate_neighbor_rows_zero": duplicate_neighbors == 0,
        "duplicate_ba_instance_ids_zero": duplicate_ba == 0,
        "selected_assembly_context_complete": len(meta["missing_context_keys"]) == 0,
        "equivalence_reason_complete": bool(eq["equivalence_reason"].astype(str).str.len().gt(0).all()),
        "pair_neighbor_counts_match_step1": bool((pd.to_numeric(pair_out["n_unique_neighbor_instances"], errors="coerce").fillna(0).astype(int) == pair_out["step1_neighbor_count_observed"]).all()),
        "positive_and_negative_controls_present": all(controls.values()),
        "processing2_frozen_validation_pass": p2_frozen.get("status") == "FROZEN" and p2_frozen.get("validation_pass") is True,
        "processing2_coordinate_validation_pass": p2_coordinate.get("validation_pass") is True,
        "coordinate_fallback_pending_zero": not eq["coordinate_verification_status"].eq("PENDING").any(),
    }
    if mode == "full":
        checks["full_expected_instance_count"] = len(eq) == expected_neighbors
        checks["full_expected_pair_count"] = len(pair_out) == expected_pairs
    return {"run_mode": mode, "validated_at": utc(), "validation_pass": all(checks.values()), "checks": checks,
            "controls": controls, "counts": {"input_pairs": len(pairs_in), "input_instances": len(neighbors_in),
            "ba_inventory_instances": len(inventory), "ba_equivalent": counts["BA_EQUIVALENT"],
            "external": counts["EXTERNAL_CRYSTAL_INSTANCE"], "review": counts["BA_EQUIVALENCE_REVIEW"],
            "same_source": int(eq["same_source_in_ba"].sum()), "source_not_in_ba": int((~eq["same_source_in_ba"]).sum()),
            "duplicate_neighbor_rows": duplicate_neighbors, "duplicate_ba_instance_ids": duplicate_ba}}


def finalize(run: Path, frames: dict[str, pd.DataFrame], validation: dict, cfg: dict,
             elapsed: float, meta: dict, pilot_selection: pd.DataFrame | None) -> None:
    output = run / "output"
    paths = {
        "01_ba_instance_inventory.tsv.gz": frames["ba_inventory"],
        "02_step1_neighbor_ba_equivalence.tsv.gz": frames["equivalence"],
        "03_external_crystal_instances.tsv.gz": frames["external"],
        "04_pair_external_neighbor_inventory.tsv.gz": frames["pairs"],
        "05_transform_match_distribution.tsv.gz": frames["distribution"],
    }
    for name, frame in paths.items(): write_tsv_gz(frame, output / name)
    summary = pd.DataFrame([{**validation["counts"], "pdb_count": int(frames["pairs"]["pdb_id"].nunique()),
                             "runtime_seconds": elapsed, "validation_pass": validation["validation_pass"]}])
    summary.to_csv(output / "06_step2_summary.tsv", sep="\t", index=False)
    if pilot_selection is not None:
        pilot_selection.to_csv(run / "pilot_selection.tsv", sep="\t", index=False)
    atomic_json(run / "validation.json", validation)
    atomic_json(run / "input_provenance.json", {
        "step1_run": cfg["input"]["step1_run"], "step1_sha256sums": sha256(Path(cfg["input"]["step1_run"]) / "SHA256SUMS"),
        "processing2_run": cfg["input"]["processing2_run"],
        "processing2_frozen": str(Path(cfg["input"]["processing2_run"]) / "_FROZEN.json"),
        "assembly_context": cfg["input"]["assembly_context"],
        "assembly_context_role": "exact expanded final R/t consumed unchanged and coordinate-validated by frozen Processing 2",
        "raw_mmcif_read": bool(meta.get("raw_mmcif_read", False)), "mmcif_root_reserved_for_coordinate_fallback_only": cfg["input"]["mmcif_root"],
        "inventory_metadata": meta, "created_at": utc()})
    schema_frames = {"01_ba_instance_inventory": frames["ba_inventory"], "02_step1_neighbor_ba_equivalence": frames["equivalence"],
                     "03_external_crystal_instances": frames["external"], "04_pair_external_neighbor_inventory": frames["pairs"],
                     "05_transform_match_distribution": frames["distribution"], "06_step2_summary": summary}
    atomic_json(run / "output_schema.json", make_schema(schema_frames))
    tol = cfg["equivalence"]
    pc = Counter(frames["pairs"]["step2_status"])
    raw_note = "yes, coordinate fallback only" if meta.get("raw_mmcif_read") else "no; transform comparison was decisive"
    report = f"""# Filter 4 Step 2 — Selected Biological Assembly Equivalence\n\nRun: `{run.name}`  \nMode: `{validation['run_mode']}`  \nStatus: `{'PASS' if validation['validation_pass'] else 'FAIL'}`\n\n## Scientific boundary\n\nThis run classifies source-object/final-transform identity only. It does not compute 4 Å contacts, vdW contacts, packing severity, rejection, or benchmark exclusion. All Step 1 6 Å fields are inherited unchanged.\n\n## Inputs and identity\n\n- Step 1 frozen run: `{cfg['input']['step1_run']}`\n- Processing 2 frozen run: `{cfg['input']['processing2_run']}`\n- BA transform source: frozen `entry_assembly_context` consumed and coordinate-validated by Processing 2\n- Raw mmCIF reread: {raw_note}\n- Polymer source key: exact Step 1 `POLYMER|model|entity|label_asym|auth_asym`\n- Non-polymer source key: exact Step 1 residue-level `NONPOLYMER|model|entity|label_asym|auth_asym|comp|auth_seq|label_seq|ins`\n- Composed operators: already expanded by Processing 2 provenance; matching uses final Cartesian R/t, never operator IDs.\n\n## Frozen equivalence criterion\n\n- rotation max-absolute tolerance: `{tol['rotation_max_abs_tolerance']}`\n- translation max-absolute tolerance: `{tol['translation_max_abs_tolerance_angstrom']}` Å\n- near-threshold review factor: `{tol['near_threshold_factor']}`\n- coordinate fallback RMSD tolerance: `{tol['coordinate_rmsd_tolerance_angstrom']}` Å\n- coordinate fallback max-atom tolerance: `{tol['max_atom_deviation_tolerance_angstrom']}` Å\n\n## Counts\n\n- PDB: {frames['pairs']['pdb_id'].nunique():,}\n- pairs: {len(frames['pairs']):,}\n- Step 1 instances: {len(frames['equivalence']):,}\n- BA inventory instances: {len(frames['ba_inventory']):,}\n- BA_EQUIVALENT: {validation['counts']['ba_equivalent']:,}\n- EXTERNAL_CRYSTAL_INSTANCE: {validation['counts']['external']:,}\n- BA_EQUIVALENCE_REVIEW: {validation['counts']['review']:,}\n- same-source in BA: {validation['counts']['same_source']:,}\n- source not in BA: {validation['counts']['source_not_in_ba']:,}\n- silent drop: {len(frames['equivalence']) - sum(validation['counts'][x] for x in ('ba_equivalent','external','review')):,}\n- runtime: {elapsed:.1f} s\n\n## Pair statuses\n\n""" + "\n".join(f"- {k}: {v:,}" for k,v in sorted(pc.items())) + f"\n\n## Validation\n\n`{'PASS' if validation['validation_pass'] else 'FAIL'}`\n"
    (output / "07_step2_report.md").write_text(report, encoding="utf-8")

    files = [p for p in run.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "output_manifest.tsv"}]
    manifest_rows = []
    for p in sorted(files):
        rel = p.relative_to(run).as_posix()
        rows = ""
        if p.suffix == ".gz":
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh: rows = max(0, sum(1 for _ in fh) - 1)
        elif p.suffix in {".tsv", ".md", ".json", ".py", ".yaml"}: rows = ""
        manifest_rows.append({"relative_path": rel, "size_bytes": p.stat().st_size, "row_count": rows, "sha256": sha256(p),
                              "schema_version": SCHEMA_VERSION, "generated_by": VERSION})
    pd.DataFrame(manifest_rows).to_csv(run / "output_manifest.tsv", sep="\t", index=False)
    checksum_files = [p for p in run.rglob("*") if p.is_file() and p.name != "SHA256SUMS"]
    with (run / "SHA256SUMS").open("w", encoding="utf-8") as fh:
        for p in sorted(checksum_files): fh.write(f"{sha256(p)}  {p.relative_to(run).as_posix()}\n")
    if validation["validation_pass"]:
        atomic_json(run / "_FROZEN.json", {"status": "FROZEN", "run_id": run.name, "stage": VERSION,
                    "validation_pass": True, "frozen_at": utc(), "sha256sums_sha256": sha256(run / "SHA256SUMS")})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--mode", choices=["pilot", "full"], required=True)
    ap.add_argument("--pilot-count", type=int, default=300)
    args = ap.parse_args()
    start = time.time(); cfg = load_config(args.config); run = args.run_dir
    if run.exists() and any(run.iterdir()): raise SystemExit(f"run directory is not empty: {run}")
    run.mkdir(parents=True, exist_ok=True); (run / "output").mkdir()
    shutil.copy2(args.config, run / "config_snapshot.yaml"); shutil.copy2(Path(__file__), run / "executed_runner.py")
    pairs, neighbors, sources = read_inputs(cfg); context = load_context(cfg)
    pairs["pdb_id"] = pairs["pdb_id"].str.lower(); sources["pdb_id"] = sources["pdb_id"].str.lower()
    selection = None
    if args.mode == "pilot":
        selection = choose_pilot(pairs, neighbors, sources, context, args.pilot_count)
        keep = set(selection["pdb_id"])
        pairs = pairs[pairs["pdb_id"].isin(keep)].copy(); neighbors = neighbors[neighbors["pdb_id"].isin(keep)].copy()
        sources = sources[sources["pdb_id"].isin(keep)].copy(); context = context[context["pdb_id"].isin(keep)].copy()
    inventory, meta = build_ba_inventory(sources, pairs, context)
    eq = compare_neighbors(neighbors, inventory, cfg)
    eq, fallback_meta = coordinate_fallback(eq, inventory, cfg)
    meta.update(fallback_meta)
    pout = pair_inventory(pairs, eq)
    dist = distribution(eq)
    frames = {"ba_inventory": inventory, "equivalence": eq,
              "external": eq[eq["equivalence_status"].eq("EXTERNAL_CRYSTAL_INSTANCE")].copy(),
              "pairs": pout, "distribution": dist}
    validation = validate(args.mode, pairs, neighbors, inventory, eq, pout, meta, cfg)
    finalize(run, frames, validation, cfg, time.time() - start, meta, selection)
    print(json.dumps(validation, indent=2))
    if not validation["validation_pass"]: raise SystemExit(2)


if __name__ == "__main__":
    main()
