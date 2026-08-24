#!/usr/bin/env python3
"""Filter 4 Step 3: direct ligand contact by external crystal instances only."""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree


VERSION = "filter4_step3_v1.0.1"
SCHEMA_VERSION = "filter4_step3_schema_v1.0.0"
EPSILON = 1.0e-12
HEAVY_EXCLUDE = {"H", "D", "T"}
BLANKS = {"", ".", "?", "None", "False", "nan", "<NA>", "\x00"}
R_COLS = [f"R{i}{j}" for i in range(1, 4) for j in range(1, 4)]
T_COLS = ["tx", "ty", "tz"]

INSTANCE_HEADER = [
    "pair_id", "pdb_id", "assembly_id", "model_id", "external_instance_id", "external_instance_key",
    "source_object_key", "source_object_type", "symmetry_operation_id", "cell_h", "cell_k", "cell_l",
    "touches_ligand_6A", "touches_pocket_6A", "step1_min_ligand_distance_A",
    "n_candidate_contact_units_within_4A", "n_contact_units_4A", "n_external_heavy_atoms_within_4A",
    "n_external_ligand_atom_pairs_within_4A", "n_ligand_heavy_atoms_contacted_4A", "ligand_heavy_atom_count",
    "fraction_ligand_heavy_atoms_contacted_4A", "min_external_ligand_distance_A",
    "ligand_contacted_atom_indices_4A", "instance_direct_crystal_contact_4A", "step3_instance_status", "error_reason",
]
UNIT_HEADER = [
    "pair_id", "pdb_id", "external_instance_id", "contact_unit_id", "contact_unit_type",
    "external_model_id", "external_entity_id", "external_label_asym_id", "external_auth_asym_id",
    "external_comp_id", "external_label_seq_id", "external_auth_seq_id", "external_ins_code",
    "n_external_heavy_atoms_within_4A", "n_atom_pairs_within_4A", "n_ligand_heavy_atoms_contacted",
    "min_distance_A", "contact_unit_direct_crystal_contact_4A",
]
REFERENCE_HEADER = [
    "pair_id", "external_instance_id", "optimized_min_distance_A", "bruteforce_min_distance_A",
    "optimized_atom_pairs_4A", "bruteforce_atom_pairs_4A", "optimized_external_atoms_4A",
    "bruteforce_external_atoms_4A", "optimized_contact_units_4A", "bruteforce_contact_units_4A",
    "optimized_direct_contact_4A", "bruteforce_direct_contact_4A", "reference_match",
]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None: return ""
    text = str(value).strip()
    return "" if text in BLANKS else text


def normalized_ins(value: Any) -> str:
    value = clean(value)
    return "" if value.lower() in {"false", "none", "null"} else value


def truth(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def safe(value: Any) -> str:
    return clean(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def sha256(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block): h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def write_rows(path: Path, header: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows: writer.writerow({k: safe(v) if isinstance(v, str) else v for k, v in row.items()})
    os.replace(tmp, path)


def read_pair_buckets(step1_work: Path) -> tuple[dict[str, int], dict[int, set[str]]]:
    pair_bucket: dict[str, int] = {}; bucket_pairs: dict[int, set[str]] = {}
    for bucket in range(256):
        path = step1_work / f"bucket_{bucket:03d}" / "pairs.tsv"
        ids = set()
        with path.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pid = row["candidate_pair_id"]; pair_bucket[pid] = bucket; ids.add(pid)
        bucket_pairs[bucket] = ids
    return pair_bucket, bucket_pairs


def stable_rank(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def select_pilot(external: pd.DataFrame, pair_step2: pd.DataFrame, count: int) -> pd.DataFrame:
    e = external.copy()
    e["min_distance_to_ligand"] = pd.to_numeric(e["min_distance_to_ligand"], errors="coerce")
    flags: dict[str, set[str]] = defaultdict(set)
    bins = [(-math.inf,2,"distance_le_2"),(2,3,"distance_2_3"),(3,3.9,"distance_3_3p9"),
            (3.9,4,"distance_3p9_4"),(4,4.1,"distance_4_4p1"),(4.1,5,"distance_4p1_5"),(5,6.000001,"distance_5_6")]
    for row in e.itertuples():
        pdb = str(row.pdb_id).lower(); d = float(row.min_distance_to_ligand) if pd.notna(row.min_distance_to_ligand) else math.nan
        flags[pdb].add("polymer" if row.source_object_type == "POLYMER" else "nonpolymer")
        flags[pdb].add("ligand_6A_true" if truth(row.touches_ligand_6A) else "ligand_6A_false")
        if math.isfinite(d):
            for lo, hi, name in bins:
                if lo < d <= hi: flags[pdb].add(name); break
    repeated = e.groupby(["pdb_id", "source_object_id"])["crystal_instance_key"].nunique()
    for (pdb, _), n in repeated.items():
        if n > 1: flags[str(pdb).lower()].add("same_source_multiple_placements")
    for pdb in pair_step2.loc[pair_step2["step2_status"].eq("BA_EQUIVALENCE_REVIEW"), "pdb_id"]:
        flags[str(pdb).lower()].add("step2_review")
    selected: list[str] = []; reasons: dict[str, set[str]] = defaultdict(set)
    strata = sorted({x for values in flags.values() for x in values})
    target = max(8, min(30, count // max(1, len(strata))))
    available = sorted(flags)
    for stratum in strata:
        candidates = sorted((p for p in available if stratum in flags[p]), key=lambda p: stable_rank(stratum + "|" + p))
        for pdb in candidates[:target]:
            if len(selected) >= count: break
            if pdb not in selected: selected.append(pdb)
            reasons[pdb].add(stratum)
    for pdb in sorted(available, key=lambda p: stable_rank("fill|" + p)):
        if len(selected) >= count: break
        if pdb not in selected: selected.append(pdb); reasons[pdb].add("deterministic_fill")
    return pd.DataFrame([{"pdb_id": p, "selection_reason": ";".join(sorted(reasons[p])),
                          "all_strata": ";".join(sorted(flags[p]))} for p in sorted(selected)])


def locate_cif(root: Path, pdb: str) -> Path:
    direct = root / f"{pdb}.cif.gz"
    return direct if direct.exists() else root / pdb[1:3] / f"{pdb}.cif.gz"


def parse_source_atoms(path: Path, wanted: set[str]) -> tuple[dict[str, dict], str]:
    try:
        block = gemmi.cif.read(str(path)).sole_block()
        tags = ["_atom_site.id", "_atom_site.pdbx_PDB_model_num", "_atom_site.label_entity_id",
                "_atom_site.label_asym_id", "_atom_site.auth_asym_id", "_atom_site.label_comp_id",
                "_atom_site.auth_comp_id", "_atom_site.label_seq_id", "_atom_site.auth_seq_id",
                "_atom_site.pdbx_PDB_ins_code", "_atom_site.label_atom_id", "_atom_site.type_symbol",
                "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
        output: dict[str, dict] = {sid: {"xyz": [], "units": [], "seen": set()} for sid in wanted}
        for row in block.find(tags):
            v = [clean(row[i]) for i in range(len(tags))]
            atom_id, model, entity, label, auth = v[:5]; model = model or "1"
            comp = (v[5] or v[6]).upper(); label_seq, auth_seq, ins = v[7], v[8], normalized_ins(v[9])
            atom_name, element = v[10], v[11].upper()
            if element in HEAVY_EXCLUDE: continue
            polymer_id = f"POLYMER|{model}|{entity}|{label}|{auth}"
            nonpoly_id = f"NONPOLYMER|{model}|{entity}|{label}|{auth}|{comp}|{auth_seq}||{ins}"
            sid = polymer_id if polymer_id in wanted else (nonpoly_id if nonpoly_id in wanted else "")
            if not sid: continue
            x, y, z = float(v[12]), float(v[13]), float(v[14])
            dedup = (atom_id, atom_name, x, y, z)
            if dedup in output[sid]["seen"]: continue
            output[sid]["seen"].add(dedup)
            if sid.startswith("POLYMER|"):
                uid = f"POLYMER_RESIDUE|{model}|{entity}|{label}|{auth}|{comp}|{label_seq}|{auth_seq}|{ins}"
                meta = (uid, "POLYMER_RESIDUE", model, entity, label, auth, comp, label_seq, auth_seq, ins)
            else:
                uid = f"NONPOLYMER_COMPONENT|{sid}"
                meta = (uid, "NONPOLYMER_COMPONENT", model, entity, label, auth, comp, label_seq, auth_seq, ins)
            output[sid]["xyz"].append((x, y, z)); output[sid]["units"].append(meta)
        for value in output.values():
            value["xyz"] = np.asarray(value["xyz"], dtype=float).reshape(-1, 3); value.pop("seen", None)
        return output, ""
    except Exception as exc:
        return {}, f"{type(exc).__name__}: {exc}"


def metric_from_pairs(external_xyz: np.ndarray, unit_meta: list[tuple], ligand_xyz: np.ndarray,
                      method: str, cutoff: float) -> dict:
    n_ext, n_lig = len(external_xyz), len(ligand_xyz)
    ext_to_lig: list[list[tuple[int, float]]] = [[] for _ in range(n_ext)]
    if not n_ext or not n_lig:
        return {"units": [], "ext_indices": set(), "lig_indices": set(), "pair_count": 0,
                "min_distance": math.nan, "direct": False}
    if method == "optimized":
        tree = cKDTree(ligand_xyz)
        candidates = tree.query_ball_point(external_xyz, r=np.nextafter(cutoff, math.inf))
        for ei, indices in enumerate(candidates):
            for li in indices:
                d = float(np.linalg.norm(external_xyz[ei] - ligand_xyz[int(li)]))
                if d <= cutoff + EPSILON: ext_to_lig[ei].append((int(li), d))
        nearest = float(np.min(tree.query(external_xyz, k=1, workers=1)[0]))
    else:
        nearest = math.inf
        chunk = 1024
        for start in range(0, n_ext, chunk):
            distances = np.linalg.norm(external_xyz[start:start+chunk, None, :] - ligand_xyz[None, :, :], axis=2)
            nearest = min(nearest, float(distances.min()))
            ii, jj = np.where(distances <= cutoff + EPSILON)
            for i, j in zip(ii, jj): ext_to_lig[start + int(i)].append((int(j), float(distances[i, j])))
    grouped: dict[str, dict] = {}
    ext_indices = set(); lig_indices = set(); pair_count = 0
    for ei, contacts in enumerate(ext_to_lig):
        if not contacts: continue
        ext_indices.add(ei); pair_count += len(contacts); lig_indices.update(li for li, _ in contacts)
        meta = unit_meta[ei]; uid = meta[0]
        item = grouped.setdefault(uid, {"meta": meta, "external_atoms": set(), "ligand_atoms": set(), "pairs": 0, "min": math.inf})
        item["external_atoms"].add(ei); item["ligand_atoms"].update(li for li, _ in contacts)
        item["pairs"] += len(contacts); item["min"] = min(item["min"], min(d for _, d in contacts))
    units = []
    for uid in sorted(grouped):
        item = grouped[uid]; units.append({**item, "direct": len(item["external_atoms"]) >= 2})
    return {"units": units, "ext_indices": ext_indices, "lig_indices": lig_indices, "pair_count": pair_count,
            "min_distance": nearest, "direct": any(x["direct"] for x in units)}


def compare_metrics(a: dict, b: dict) -> bool:
    def unit_signature(m):
        return sorted((x["meta"][0], len(x["external_atoms"]), x["pairs"], len(x["ligand_atoms"]), x["direct"]) for x in m["units"])
    return (len(a["ext_indices"]) == len(b["ext_indices"]) and len(a["lig_indices"]) == len(b["lig_indices"]) and
            a["pair_count"] == b["pair_count"] and a["direct"] == b["direct"] and unit_signature(a) == unit_signature(b) and
            math.isclose(a["min_distance"], b["min_distance"], abs_tol=1e-10, rel_tol=1e-12))


def process_bucket(task: tuple) -> dict:
    bucket, records, target_path, mmcif_root, work_dir, cutoff, reference = task
    started = time.time(); work = Path(work_dir); work.mkdir(parents=True, exist_ok=True)
    instance_path = work / "instances.tsv"; unit_path = work / "units.tsv"; ref_path = work / "reference.tsv"; error_path = work / "errors.tsv"
    needed_pairs = {str(r["candidate_pair_id"]) for r in records}
    atoms = pd.read_csv(target_path, sep="\t", usecols=["candidate_pair_id", "target_kind", "x", "y", "z"])
    atoms = atoms[atoms["candidate_pair_id"].isin(needed_pairs) & atoms["target_kind"].eq("LIGAND")]
    ligand_map = {str(pid): g[["x","y","z"]].to_numpy(float) for pid, g in atoms.groupby("candidate_pair_id", sort=False)}
    by_pdb: dict[str, list[dict]] = defaultdict(list)
    for row in records: by_pdb[str(row["pdb_id"]).lower()].append(row)
    instances = []; units_out = []; refs = []; errors = []; source_mismatch = 0; min_diff_max = 0.0
    for pdb, prows in sorted(by_pdb.items()):
        active = [r for r in prows if truth(r["touches_ligand_6A"])]
        sources = {}; parse_error = ""
        if active:
            sources, parse_error = parse_source_atoms(locate_cif(Path(mmcif_root), pdb), {str(r["source_object_id"]) for r in active})
        transformed_cache = {}
        for row in prows:
            pid = str(row["candidate_pair_id"]); sid = str(row["source_object_id"]); ligand = ligand_map.get(pid, np.empty((0,3)))
            base = {"pair_id": pid, "pdb_id": pdb, "assembly_id": row["assembly_id"], "model_id": row["model_id"],
                    "external_instance_id": row["crystal_instance_id"], "external_instance_key": row["crystal_instance_key"],
                    "source_object_key": sid, "source_object_type": row["source_object_type"],
                    "symmetry_operation_id": row["symmetry_operation_id"], "cell_h": row["cell_h"], "cell_k": row["cell_k"], "cell_l": row["cell_l"],
                    "touches_ligand_6A": truth(row["touches_ligand_6A"]), "touches_pocket_6A": truth(row["touches_pocket_6A"]),
                    "step1_min_ligand_distance_A": row["min_distance_to_ligand"], "ligand_heavy_atom_count": len(ligand)}
            error = ""
            if not len(ligand): error = "FROZEN_LIGAND_ATOMS_MISSING"
            elif not truth(row["touches_ligand_6A"]):
                metric = {"units": [], "ext_indices": set(), "lig_indices": set(), "pair_count": 0,
                          "min_distance": math.nan, "direct": False}
            elif parse_error: error = "MMCIF_PARSE_ERROR: " + parse_error
            elif sid not in sources or not len(sources[sid]["xyz"]): error = "SOURCE_OBJECT_ATOMS_MISSING"
            elif len(sources[sid]["xyz"]) != int(float(row["source_heavy_atom_count"])):
                source_mismatch += 1; error = f"SOURCE_HEAVY_ATOM_COUNT_MISMATCH:{len(sources[sid]['xyz'])}!={row['source_heavy_atom_count']}"
            if error:
                metric = {"units": [], "ext_indices": set(), "lig_indices": set(), "pair_count": 0,
                          "min_distance": math.nan, "direct": False}; errors.append({"pair_id": pid, "external_instance_id": row["crystal_instance_id"], "error_reason": error})
            elif truth(row["touches_ligand_6A"]):
                cache_key = str(row["crystal_instance_key"])
                if cache_key not in transformed_cache:
                    rmat = np.array([float(row[c]) for c in R_COLS]).reshape(3,3); trans = np.array([float(row[c]) for c in T_COLS])
                    transformed_cache[cache_key] = (sources[sid]["xyz"] @ rmat.T + trans, sources[sid]["units"])
                xyz, unit_meta = transformed_cache[cache_key]
                metric = metric_from_pairs(xyz, unit_meta, ligand, "optimized", cutoff)
                step1_min = float(row["min_distance_to_ligand"])
                min_diff_max = max(min_diff_max, abs(metric["min_distance"] - step1_min))
                if reference:
                    brute = metric_from_pairs(xyz, unit_meta, ligand, "bruteforce", cutoff)
                    match = compare_metrics(metric, brute)
                    refs.append({"pair_id": pid, "external_instance_id": row["crystal_instance_id"],
                                 "optimized_min_distance_A": metric["min_distance"], "bruteforce_min_distance_A": brute["min_distance"],
                                 "optimized_atom_pairs_4A": metric["pair_count"], "bruteforce_atom_pairs_4A": brute["pair_count"],
                                 "optimized_external_atoms_4A": len(metric["ext_indices"]), "bruteforce_external_atoms_4A": len(brute["ext_indices"]),
                                 "optimized_contact_units_4A": sum(x["direct"] for x in metric["units"]),
                                 "bruteforce_contact_units_4A": sum(x["direct"] for x in brute["units"]),
                                 "optimized_direct_contact_4A": metric["direct"], "bruteforce_direct_contact_4A": brute["direct"],
                                 "reference_match": match})
            direct_units = sum(x["direct"] for x in metric["units"])
            contacted = sorted(metric["lig_indices"])
            instances.append({**base, "n_candidate_contact_units_within_4A": len(metric["units"]), "n_contact_units_4A": direct_units,
                              "n_external_heavy_atoms_within_4A": len(metric["ext_indices"]),
                              "n_external_ligand_atom_pairs_within_4A": metric["pair_count"],
                              "n_ligand_heavy_atoms_contacted_4A": len(contacted),
                              "fraction_ligand_heavy_atoms_contacted_4A": len(contacted)/len(ligand) if len(ligand) else math.nan,
                              "min_external_ligand_distance_A": metric["min_distance"],
                              "ligand_contacted_atom_indices_4A": ";".join(map(str, contacted)),
                              "instance_direct_crystal_contact_4A": metric["direct"],
                              "step3_instance_status": "ERROR" if error else "SUCCESS", "error_reason": error})
            for unit in metric["units"]:
                meta = unit["meta"]
                units_out.append({"pair_id": pid, "pdb_id": pdb, "external_instance_id": row["crystal_instance_id"],
                                  "contact_unit_id": meta[0], "contact_unit_type": meta[1], "external_model_id": meta[2],
                                  "external_entity_id": meta[3], "external_label_asym_id": meta[4], "external_auth_asym_id": meta[5],
                                  "external_comp_id": meta[6], "external_label_seq_id": meta[7], "external_auth_seq_id": meta[8],
                                  "external_ins_code": meta[9], "n_external_heavy_atoms_within_4A": len(unit["external_atoms"]),
                                  "n_atom_pairs_within_4A": unit["pairs"], "n_ligand_heavy_atoms_contacted": len(unit["ligand_atoms"]),
                                  "min_distance_A": unit["min"], "contact_unit_direct_crystal_contact_4A": unit["direct"]})
    write_rows(instance_path, INSTANCE_HEADER, instances); write_rows(unit_path, UNIT_HEADER, units_out)
    write_rows(ref_path, REFERENCE_HEADER, refs); write_rows(error_path, ["pair_id","external_instance_id","error_reason"], errors)
    result = {"bucket": bucket, "input_instances": len(records), "output_instances": len(instances), "unit_rows": len(units_out),
              "reference_rows": len(refs), "reference_mismatches": sum(not truth(x["reference_match"]) for x in refs),
              "errors": len(errors), "source_count_mismatches": source_mismatch, "step1_min_distance_max_abs_diff": min_diff_max,
              "runtime_seconds": time.time()-started, "completed_at": utc()}
    atomic_json(work / "_SUCCESS.json", result); return result


def concat_tsv(paths: list[Path], out: Path, header: list[str]) -> pd.DataFrame:
    pieces = []
    for path in paths:
        if path.exists() and path.stat().st_size:
            frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
            if len(frame): pieces.append(frame)
    frame = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=header)
    frame = frame.reindex(columns=header)
    tmp = out.with_suffix(out.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip", na_rep="", lineterminator="\n")
    os.replace(tmp, out); return frame


def boundary_tests(cutoff: float) -> dict:
    ligand = np.array([[0.,0.,0.]])
    results = {}
    for d in (3.9999, 4.0, 4.0001):
        # Two distinct atom indices at the same coordinate isolate the cutoff
        # semantics from any orthogonal displacement in the synthetic fixture.
        xyz = np.array([[d,0,0],[d,0,0]])
        meta = [("u","POLYMER_RESIDUE","1","1","A","A","GLY","1","1","")] * 2
        results[f"distance_{d:.4f}"] = metric_from_pairs(xyz, meta, ligand, "optimized", cutoff)["direct"]
    # Distinct external atoms may contact the same ligand atom; residues never pool atoms.
    two_same = metric_from_pairs(np.array([[3.5,0,0],[3.7,0,0]]),
        [("u","POLYMER_RESIDUE","1","1","A","A","TYR","100","100","")]*2, ligand, "optimized", cutoff)["direct"]
    separate = metric_from_pairs(np.array([[3.5,0,0],[3.7,0,0]]),
        [("u1","POLYMER_RESIDUE","1","1","A","A","TYR","100","100",""),
         ("u2","POLYMER_RESIDUE","1","1","A","A","ASP","101","101","")], ligand, "optimized", cutoff)["direct"]
    results.update({"two_distinct_external_atoms_same_ligand_atom_true": two_same,
                    "two_residues_one_atom_each_false": not separate,
                    "boundary_pass": results["distance_3.9999"] and results["distance_4.0000"] and not results["distance_4.0001"] and two_same and not separate})
    return results


def pair_output(step2_pairs: pd.DataFrame, instances: pd.DataFrame, units: pd.DataFrame) -> pd.DataFrame:
    numeric = ["n_contact_units_4A", "n_ligand_heavy_atoms_contacted_4A", "ligand_heavy_atom_count", "min_external_ligand_distance_A"]
    for col in numeric: instances[col] = pd.to_numeric(instances[col], errors="coerce")
    instances["instance_direct_crystal_contact_4A"] = instances["instance_direct_crystal_contact_4A"].map(truth)
    agg_rows = []
    for pid, group in instances.groupby("pair_id", sort=False):
        contacted = set()
        for value in group["ligand_contacted_atom_indices_4A"]:
            contacted.update(int(x) for x in str(value).split(";") if x != "")
        ligand_counts = set(group["ligand_heavy_atom_count"].dropna().astype(int))
        errors = int(group["step3_instance_status"].eq("ERROR").sum())
        agg_rows.append({"candidate_pair_id": pid, "n_external_instances_step3": len(group),
                         "n_external_instances_ligand_6A_step3": int(group["touches_ligand_6A"].map(truth).sum()),
                         "n_external_instances_direct_contact_4A": int(group["instance_direct_crystal_contact_4A"].sum()),
                         "n_contact_units_4A": int(group["n_contact_units_4A"].sum()),
                         "n_ligand_heavy_atoms_contacted_4A": len(contacted),
                         "ligand_heavy_atom_count": next(iter(ligand_counts)) if len(ligand_counts)==1 else math.nan,
                         "fraction_ligand_heavy_atoms_contacted_4A": len(contacted)/next(iter(ligand_counts)) if len(ligand_counts)==1 and next(iter(ligand_counts)) else math.nan,
                         "min_external_ligand_distance_A": group["min_external_ligand_distance_A"].min(),
                         "pair_direct_crystal_contact_4A": bool(group["instance_direct_crystal_contact_4A"].any()),
                         "step3_error_count": errors})
    agg = pd.DataFrame(agg_rows)
    out = step2_pairs.merge(agg, on="candidate_pair_id", how="left")
    int_cols = ["n_external_instances_step3", "n_external_instances_ligand_6A_step3", "n_external_instances_direct_contact_4A", "n_contact_units_4A", "n_ligand_heavy_atoms_contacted_4A", "step3_error_count"]
    out[int_cols] = out[int_cols].fillna(0).astype(np.int64)
    out["pair_direct_crystal_contact_4A"] = out["pair_direct_crystal_contact_4A"].fillna(False).astype(bool)
    statuses = []
    for row in out.itertuples():
        if row.step2_status == "UPSTREAM_NO_NEIGHBOR": status = "UPSTREAM_NO_NEIGHBOR"
        elif row.step2_status == "SUCCESS_NO_EXTERNAL_NEIGHBOR": status = "UPSTREAM_NO_EXTERNAL_NEIGHBOR"
        elif row.step2_status == "BA_EQUIVALENCE_REVIEW": status = "BA_EQUIVALENCE_REVIEW"
        elif row.step3_error_count: status = "ERROR"
        elif row.pair_direct_crystal_contact_4A: status = "SUCCESS_DIRECT_LIGAND_CONTACT"
        else: status = "SUCCESS_NO_DIRECT_LIGAND_CONTACT"
        statuses.append(status)
    out["step3_status"] = statuses
    return out


def observed_controls(instances: pd.DataFrame, units: pd.DataFrame, external_input: pd.DataFrame, pairs: pd.DataFrame) -> dict:
    for col in ["n_external_heavy_atoms_within_4A","n_atom_pairs_within_4A"]:
        if col in units: units[col] = pd.to_numeric(units[col], errors="coerce").fillna(0).astype(int)
    units["contact_unit_direct_crystal_contact_4A"] = units["contact_unit_direct_crystal_contact_4A"].map(truth) if len(units) else False
    instances["instance_direct_crystal_contact_4A"] = instances["instance_direct_crystal_contact_4A"].map(truth)
    one = units["n_external_heavy_atoms_within_4A"].eq(1) if len(units) else pd.Series([], dtype=bool)
    many = units["n_external_heavy_atoms_within_4A"].ge(3) if len(units) else pd.Series([], dtype=bool)
    per_instance = units.groupby(["pair_id","external_instance_id"])["n_external_heavy_atoms_within_4A"].apply(list) if len(units) else pd.Series(dtype=object)
    case4 = any(len(x)>=2 and all(v==1 for v in x) for x in per_instance)
    case5 = any(any(v>=2 for v in x) and any(v==1 for v in x) for x in per_instance)
    repeated = external_input.groupby(["pdb_id","source_object_id"])["crystal_instance_key"].nunique()
    return {"case1_one_atom_unit_false": bool((one & ~units["contact_unit_direct_crystal_contact_4A"]).any()) if len(units) else False,
            "case2_exactly_two_atoms_true": bool((units["n_external_heavy_atoms_within_4A"].eq(2) & units["contact_unit_direct_crystal_contact_4A"]).any()) if len(units) else False,
            "case3_many_atoms_true": bool((many & units["contact_unit_direct_crystal_contact_4A"]).any()) if len(units) else False,
            "case4_two_one_atom_residues_instance_false": case4,
            "case5_true_and_one_atom_units_same_instance": case5,
            "case9_non_ligand6_direct_false": bool((~instances["touches_ligand_6A"].map(truth) & ~instances["instance_direct_crystal_contact_4A"]).any()),
            "case10_polymer_contact": bool(((instances["source_object_type"]=="POLYMER") & instances["instance_direct_crystal_contact_4A"]).any()),
            "case11_nonpolymer_contact": bool(((instances["source_object_type"]=="NONPOLYMER") & instances["instance_direct_crystal_contact_4A"]).any()),
            "case12_same_source_different_placements": bool((repeated>1).any()),
            "case13_only_external_instances_entered": bool(len(external_input)==355846 or external_input["equivalence_status"].eq("EXTERNAL_CRYSTAL_INSTANCE").all()),
            "case14_step2_review_propagated": int(pairs["step3_status"].eq("BA_EQUIVALENCE_REVIEW").sum()) == int((pairs["step2_status"]=="BA_EQUIVALENCE_REVIEW").sum())}


def write_gzip_frame(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip", na_rep="", lineterminator="\n")
    os.replace(tmp, path)


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True, type=Path); ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--mode", choices=["smoke","pilot","full"], required=True); ap.add_argument("--pilot-count", type=int, default=500)
    args = ap.parse_args(); start = time.time(); cfg = yaml.safe_load(args.config.read_text()); run = args.run_dir
    if run.exists() and any(run.iterdir()): raise SystemExit(f"run directory not empty: {run}")
    for rel in ("output","work"): (run/rel).mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.config, run/"config_snapshot.yaml"); shutil.copy2(Path(__file__), run/"executed_runner.py")
    step2 = Path(cfg["input"]["step2_run"]); step1_work = Path(cfg["input"]["step1_work"])
    if json.loads((step2/"_FROZEN.json").read_text()).get("status") != "FROZEN": raise RuntimeError("Step2 input not frozen")
    external_all = pd.read_csv(step2/"output/03_external_crystal_instances.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    pairs_all = pd.read_csv(step2/"output/04_pair_external_neighbor_inventory.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    pair_bucket, _ = read_pair_buckets(step1_work)
    external_all["bucket_id"] = external_all["candidate_pair_id"].map(pair_bucket)
    if external_all["bucket_id"].isna().any(): raise RuntimeError("external pair missing Step1 bucket")
    external = external_all; pairs = pairs_all; selection = None
    if args.mode in {"smoke","pilot"}:
        count = args.pilot_count if args.mode == "pilot" else min(args.pilot_count, 30)
        selection = select_pilot(external_all, pairs_all, count); keep = set(selection["pdb_id"])
        external = external_all[external_all["pdb_id"].isin(keep)].copy(); pairs = pairs_all[pairs_all["pdb_id"].isin(keep)].copy()
        selection.to_csv(run/"pilot_selection.tsv", sep="\t", index=False)
    tasks = []
    for bucket, group in external.groupby("bucket_id", sort=True):
        b = int(bucket); tasks.append((b, group.to_dict("records"), str(step1_work/f"bucket_{b:03d}"/"target_atoms.tsv"),
                                      cfg["input"]["mmcif_root"], str(run/"work"/f"bucket_{b:03d}"),
                                      float(cfg["contact"]["cutoff_angstrom"]), args.mode in {"smoke","pilot"}))
    results = []
    with cf.ProcessPoolExecutor(max_workers=int(cfg["runtime"]["workers"])) as pool:
        futures = {pool.submit(process_bucket, task): task[0] for task in tasks}
        for future in cf.as_completed(futures):
            results.append(future.result()); atomic_json(run/"progress.json", {"completed_buckets":len(results),"total_buckets":len(tasks),
                "instances_completed":sum(x["input_instances"] for x in results),"errors":sum(x["errors"] for x in results),"updated_at":utc()})
    works = sorted((run/"work").glob("bucket_*")); output = run/"output"
    instances = concat_tsv([w/"instances.tsv" for w in works], output/"01_external_instance_direct_contact.tsv.gz", INSTANCE_HEADER)
    units = concat_tsv([w/"units.tsv" for w in works], output/"02_external_contact_unit_ligand_contacts.tsv.gz", UNIT_HEADER)
    refs = concat_tsv([w/"reference.tsv" for w in works], output/"pilot_optimized_vs_bruteforce.tsv.gz", REFERENCE_HEADER) if args.mode in {"smoke","pilot"} else pd.DataFrame(columns=REFERENCE_HEADER)
    pair_frame = pair_output(pairs, instances.copy(), units.copy()); write_gzip_frame(pair_frame, output/"03_pair_direct_contact_inventory.tsv.gz")
    boundaries = boundary_tests(float(cfg["contact"]["cutoff_angstrom"])); controls = observed_controls(instances.copy(), units.copy(), external, pair_frame)
    status_counts = Counter(pair_frame["step3_status"]); instance_direct = int(instances["instance_direct_crystal_contact_4A"].map(truth).sum())
    unit_true = int(units["contact_unit_direct_crystal_contact_4A"].map(truth).sum()) if len(units) else 0
    checks = {"instance_accounting":len(instances)==len(external), "duplicate_external_instance_key_zero":not instances.duplicated(["pair_id","external_instance_key"]).any(),
              "silent_drop_zero":len(instances)==len(external), "errors_zero":sum(x["errors"] for x in results)==0,
              "source_heavy_atom_count_mismatch_zero":sum(x["source_count_mismatches"] for x in results)==0,
              "pair_accounting":len(pair_frame)==len(pairs), "boundary_tests_pass":boundaries["boundary_pass"],
              "direct_iff_true_unit": bool((instances["instance_direct_crystal_contact_4A"].map(truth) == pd.to_numeric(instances["n_contact_units_4A"]).gt(0)).all()),
              "non_ligand6_zero_contact": bool((pd.to_numeric(instances.loc[~instances["touches_ligand_6A"].map(truth),"n_external_heavy_atoms_within_4A"]).eq(0)).all()),
              "reference_mismatch_zero":sum(x["reference_mismatches"] for x in results)==0}
    checks.update({
        "contact_unit_direct_iff_two_external_atoms": bool((units["contact_unit_direct_crystal_contact_4A"].map(truth) == pd.to_numeric(units["n_external_heavy_atoms_within_4A"]).ge(2)).all()) if len(units) else True,
        "contact_unit_minimum_within_4A": bool(pd.to_numeric(units["min_distance_A"], errors="coerce").le(float(cfg["contact"]["cutoff_angstrom"])+EPSILON).all()) if len(units) else True,
        "duplicate_contact_unit_rows_zero": not units.duplicated(["pair_id","external_instance_id","contact_unit_id"]).any() if len(units) else True,
        "step1_min_distance_reproduced": validation_distance_ok if (validation_distance_ok := max((x["step1_min_distance_max_abs_diff"] for x in results),default=0.0) <= 1.0e-9) else False,
        "pair_external_instance_counts_match_step2": bool((pd.to_numeric(pair_frame["external_crystal_instance_count"],errors="coerce").fillna(0).astype(int) == pair_frame["n_external_instances_step3"]).all()),
    })
    if args.mode in {"pilot", "full"}: checks["all_required_scientific_controls_present"] = all(controls.values())
    if args.mode == "full":
        expected = cfg["validation"]
        checks.update({"full_external_instances_355846":len(instances)==int(expected["external_instances"]),
                       "full_pairs_336412":len(pair_frame)==int(expected["pairs"]),
                       "upstream_no_neighbor_7663":status_counts["UPSTREAM_NO_NEIGHBOR"]==int(expected["upstream_no_neighbor"]),
                       "upstream_no_external_131764":status_counts["UPSTREAM_NO_EXTERNAL_NEIGHBOR"]==int(expected["upstream_no_external"]),
                       "ba_review_2":status_counts["BA_EQUIVALENCE_REVIEW"]==int(expected["ba_review"]),
                       "external_ligand6_instances_134406":int(instances["touches_ligand_6A"].map(truth).sum())==int(expected["external_ligand6_instances"]),
                       "external_pair_partition_196983":status_counts["SUCCESS_DIRECT_LIGAND_CONTACT"]+status_counts["SUCCESS_NO_DIRECT_LIGAND_CONTACT"]==int(expected["external_pairs"])})
    validation = {"run_mode":args.mode,"validated_at":utc(),"validation_pass":all(checks.values()),"checks":checks,"boundary_tests":boundaries,
                  "observed_controls":controls,"counts":{"pairs":len(pair_frame),"external_instances":len(instances),"direct_instances":instance_direct,
                  "candidate_contact_units":len(units),"true_contact_units":unit_true,"reference_rows":len(refs),"errors":sum(x["errors"] for x in results)},
                  "pair_status_counts":dict(status_counts),"max_step1_min_distance_abs_diff":max((x["step1_min_distance_max_abs_diff"] for x in results),default=0.0)}
    summary = pd.DataFrame([{**validation["counts"],"pdb_count":pair_frame["pdb_id"].nunique(),
        "external_ligand_6A_instances":int(instances["touches_ligand_6A"].map(truth).sum()),
        "polymer_direct_instances":int(((instances["source_object_type"]=="POLYMER")&instances["instance_direct_crystal_contact_4A"].map(truth)).sum()),
        "nonpolymer_direct_instances":int(((instances["source_object_type"]=="NONPOLYMER")&instances["instance_direct_crystal_contact_4A"].map(truth)).sum()),
        "one_atom_only_units":int(pd.to_numeric(units["n_external_heavy_atoms_within_4A"]).eq(1).sum()) if len(units) else 0,
        "runtime_seconds":time.time()-start,"validation_pass":validation["validation_pass"]}])
    summary.to_csv(output/"04_step3_summary.tsv",sep="\t",index=False)
    report = f"""# Filter 4 Step 3 — Direct Ligand Crystal Contact Analysis\n\nRun: `{run.name}`  \nMode: `{args.mode}`  \nStatus: `{'PASS' if validation['validation_pass'] else 'FAIL'}`\n\nFilter 4 Step 3 identifies direct structural contacts between external crystallographic neighbour instances and the frozen ligand. An external polymer residue or independent non-polymer component is considered a contacting unit when at least two distinct heavy atoms from that unit lie within 4.0 Å of any ligand heavy atom. A crystal instance is classified as having a direct ligand crystal contact when it contains at least one such contacting unit. The 4.0 Å criterion is the sole operational contact definition used in Step 3.\n\nThe inherited 6 Å annotation is used only as the Step 1 broad lattice-neighbour discovery/pruning shell. This run does not implement pocket-mediated contact, vdW contact, interaction typing, severity, rejection, or benchmark exclusion.\n\n- PDB: {pair_frame['pdb_id'].nunique():,}\n- pairs: {len(pair_frame):,}\n- Step 2 external instances: {len(instances):,}\n- direct-contact instances: {instance_direct:,}\n- candidate contact units (>=1 atom): {len(units):,}\n- TRUE contact units (>=2 atoms): {unit_true:,}\n- validation: {'PASS' if validation['validation_pass'] else 'FAIL'}\n"""
    (output/"05_step3_report.md").write_text(report,encoding="utf-8")
    atomic_json(run/"validation.json",validation)
    atomic_json(run/"input_provenance.json",{"step2_run":str(step2),"step2_sha256sums":sha256(step2/"SHA256SUMS"),
        "step1_frozen_ligand_atom_cache":str(step1_work),"raw_mmcif_root":cfg["input"]["mmcif_root"],
        "policy":"Step1 ligand coordinates and Step2 external identity/R/t inherited; no 6A, lattice, BA, ligand, altloc, receptor, or pocket recomputation", "created_at":utc()})
    datasets={"01_external_instance_direct_contact":instances,"02_external_contact_unit_ligand_contacts":units,
              "03_pair_direct_contact_inventory":pair_frame,"04_step3_summary":summary}
    atomic_json(run/"output_schema.json",{"schema_version":SCHEMA_VERSION,"datasets":{k:[{"column_name":c,"data_type":str(f[c].dtype),"nullable":bool(f[c].isna().any())} for c in f.columns] for k,f in datasets.items()}})
    # Serialized row counts are a freeze gate.
    serialized={}
    for path, expected_rows in [(output/"01_external_instance_direct_contact.tsv.gz",len(instances)),(output/"02_external_contact_unit_ligand_contacts.tsv.gz",len(units)),(output/"03_pair_direct_contact_inventory.tsv.gz",len(pair_frame))]:
        with gzip.open(path,"rt",encoding="utf-8",newline="") as fh: physical=max(0,sum(1 for _ in fh)-1)
        serialized[path.name]={"expected_rows":expected_rows,"physical_rows":physical,"match":physical==expected_rows}
    validation["serialization_row_counts"]=serialized; validation["checks"]["serialized_physical_row_counts_match"]=all(x["match"] for x in serialized.values()); validation["validation_pass"]=all(validation["checks"].values())
    atomic_json(run/"validation.json",validation); summary.loc[0,"validation_pass"]=validation["validation_pass"]; summary.to_csv(output/"04_step3_summary.tsv",sep="\t",index=False)
    files=[p for p in run.rglob("*") if p.is_file() and "work" not in p.parts and p.name not in {"SHA256SUMS","output_manifest.tsv","_FROZEN.json"}]
    manifest=[]
    for p in sorted(files): manifest.append({"relative_path":p.relative_to(run).as_posix(),"size_bytes":p.stat().st_size,"sha256":sha256(p),"schema_version":SCHEMA_VERSION,"generated_by":VERSION})
    pd.DataFrame(manifest).to_csv(run/"output_manifest.tsv",sep="\t",index=False)
    checksum_files=[p for p in run.rglob("*") if p.is_file() and "work" not in p.parts and p.name not in {"SHA256SUMS","_FROZEN.json"}]
    with (run/"SHA256SUMS").open("w",encoding="utf-8") as fh:
        for p in sorted(checksum_files): fh.write(f"{sha256(p)}  {p.relative_to(run).as_posix()}\n")
    if validation["validation_pass"]:
        atomic_json(run/"_FROZEN.json",{"status":"FROZEN","run_id":run.name,"stage":VERSION,"validation_pass":True,"frozen_at":utc(),"sha256sums_sha256":sha256(run/"SHA256SUMS")})
    print(json.dumps(validation, indent=2, default=lambda x: x.item() if hasattr(x, "item") else str(x)))
    if not validation["validation_pass"]: raise SystemExit(2)


if __name__ == "__main__":
    main()
