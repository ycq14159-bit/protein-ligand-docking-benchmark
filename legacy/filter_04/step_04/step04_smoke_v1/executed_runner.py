#!/usr/bin/env python3
"""Filter 4 Step 4: one-hop external-crystal -> frozen-binding-residue contact."""

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
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import yaml
from scipy.spatial import cKDTree


VERSION = "filter4_step4_v1.0.0"
SCHEMA_VERSION = "filter4_step4_schema_v1.0.0"
EPSILON = 1.0e-12
HEAVY_EXCLUDE = {"H", "D", "T"}
BLANKS = {"", ".", "?", "None", "False", "nan", "<NA>", "\x00"}
R_COLS = [f"R{i}{j}" for i in range(1, 4) for j in range(1, 4)]
T_COLS = ["tx", "ty", "tz"]

INSTANCE_HEADER = [
    "pair_id", "pdb_id", "assembly_id", "model_id", "external_instance_id", "external_instance_key",
    "source_object_key", "source_object_type", "symmetry_operation_id", "cell_h", "cell_k", "cell_l",
    "touches_ligand_6A", "touches_pocket_6A", "step1_min_pocket_distance_A",
    "n_frozen_binding_residues", "n_binding_residues_contacted_4A",
    "n_external_contact_units_contacting_binding_residues_4A",
    "n_external_heavy_atoms_in_binding_residue_contacts_4A",
    "min_external_binding_residue_distance_A", "instance_binding_residue_crystal_contact_4A",
    "step4_instance_status", "error_reason",
]
CONTACT_HEADER = [
    "pair_id", "pdb_id", "external_instance_id", "contact_unit_id", "contact_unit_type",
    "external_model_id", "external_entity_id", "external_label_asym_id", "external_auth_asym_id",
    "external_comp_id", "external_label_seq_id", "external_auth_seq_id", "external_ins_code",
    "binding_residue_id", "binding_chain_instance_id", "binding_protein_residue_id",
    "binding_model_id", "binding_entity_id", "binding_label_asym_id", "binding_auth_asym_id",
    "binding_comp_id", "binding_label_seq_id", "binding_auth_seq_id", "binding_ins_code",
    "n_external_heavy_atoms_within_4A", "n_atom_pairs_within_4A",
    "n_binding_residue_heavy_atoms_contacted", "min_distance_A",
    "external_to_binding_residue_contact_4A",
]
BRIDGED_HEADER = [
    "pair_id", "pdb_id", "binding_residue_id", "binding_chain_instance_id", "binding_protein_residue_id",
    "binding_model_id", "binding_entity_id", "binding_label_asym_id", "binding_auth_asym_id",
    "binding_comp_id", "binding_label_seq_id", "binding_auth_seq_id", "binding_ins_code",
    "binding_heavy_atom_count", "n_external_instances_contacting", "n_external_contact_units_contacting",
    "n_external_heavy_atoms_contributing", "min_external_distance_A", "crystal_bridged_binding_residue",
]
REFERENCE_HEADER = [
    "pair_id", "external_instance_id", "optimized_min_distance_A", "bruteforce_min_distance_A",
    "optimized_atom_pairs_4A", "bruteforce_atom_pairs_4A", "optimized_external_atoms_4A",
    "bruteforce_external_atoms_4A", "optimized_contact_combinations_4A",
    "bruteforce_contact_combinations_4A", "optimized_bridged_binding_residues",
    "bruteforce_bridged_binding_residues", "optimized_instance_contact_4A",
    "bruteforce_instance_contact_4A", "reference_match",
]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value: Any) -> str:
    if value is None:
        return ""
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
        while chunk := fh.read(block):
            h.update(chunk)
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
        for row in rows:
            writer.writerow({k: safe(v) if isinstance(v, str) else v for k, v in row.items()})
    os.replace(tmp, path)


def write_gzip_frame(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip", na_rep="", lineterminator="\n")
    os.replace(tmp, path)


def concat_tsv(paths: list[Path], out: Path, header: list[str]) -> pd.DataFrame:
    pieces = []
    with gzip.open(out, "wt", encoding="utf-8", newline="") as dst:
        dst.write("\t".join(header) + "\n")
        for path in paths:
            if not path.exists():
                continue
            frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
            if len(frame):
                frame = frame.reindex(columns=header)
                frame.to_csv(dst, sep="\t", index=False, header=False, lineterminator="\n")
                pieces.append(frame)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=header)


def bucket_path(root: Path, bucket: int) -> Path:
    plain, padded = root / f"bucket_id={bucket}", root / f"bucket_id={bucket:03d}"
    return plain if plain.exists() else padded


def read_bucket(root: Path, bucket: int, columns: list[str], filter_expression=None) -> pd.DataFrame:
    path = bucket_path(root, bucket)
    if not path.exists():
        return pd.DataFrame(columns=columns)
    return ds.dataset(path, format="parquet").to_table(columns=columns, filter=filter_expression).to_pandas(
        split_blocks=True, self_destruct=True
    )


def read_pair_buckets(step1_work: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    for bucket in range(256):
        path = step1_work / f"bucket_{bucket:03d}" / "pairs.tsv"
        with path.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                result[row["candidate_pair_id"]] = bucket
    return result


def stable_rank(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def select_pilot(external: pd.DataFrame, all_pairs: pd.DataFrame, count: int) -> pd.DataFrame:
    flags: dict[str, set[str]] = defaultdict(set)
    for row in external.itertuples():
        pdb = str(row.pdb_id).lower()
        flags[pdb].add("polymer" if row.source_object_type == "POLYMER" else "nonpolymer")
        flags[pdb].add("pocket_6A_true" if truth(row.touches_pocket_6A) else "pocket_6A_false")
        d = pd.to_numeric(pd.Series([row.min_distance_to_pocket]), errors="coerce").iloc[0]
        if pd.notna(d):
            flags[pdb].add("pocket_le_4" if float(d) <= 4.0 + EPSILON else "pocket_gt_4")
    repeated = external.groupby(["pdb_id", "source_object_id"])["crystal_instance_key"].nunique()
    for (pdb, _), n in repeated.items():
        if n > 1:
            flags[str(pdb).lower()].add("same_source_multiple_placements")
    for row in all_pairs.itertuples():
        if row.step3_status == "SUCCESS_DIRECT_LIGAND_CONTACT": flags[str(row.pdb_id).lower()].add("step3_direct_excluded")
        if row.step3_status == "BA_EQUIVALENCE_REVIEW": flags[str(row.pdb_id).lower()].add("ba_review")
    selected: list[str] = []
    reasons: dict[str, set[str]] = defaultdict(set)
    strata = sorted({x for values in flags.values() for x in values})
    quota = max(20, min(120, count // max(1, len(strata))))
    for stratum in strata:
        candidates = sorted((p for p in flags if stratum in flags[p]), key=lambda p: stable_rank(stratum + "|" + p))
        for pdb in candidates[:quota]:
            if len(selected) >= count: break
            if pdb not in selected: selected.append(pdb)
            reasons[pdb].add(stratum)
    for pdb in sorted(flags, key=lambda p: stable_rank("fill|" + p)):
        if len(selected) >= count: break
        if pdb not in selected: selected.append(pdb); reasons[pdb].add("deterministic_fill")
    return pd.DataFrame([{"pdb_id": p, "selection_reason": ";".join(sorted(reasons[p])),
                          "all_strata": ";".join(sorted(flags[p]))} for p in sorted(selected)])


def audit_binding_subset(pair_ids: set[str], pair_bucket: dict[str, int], p3_output: Path) -> dict:
    by_bucket: dict[int, set[str]] = defaultdict(set)
    for pid in pair_ids:
        by_bucket[pair_bucket[pid]].add(pid)
    binding_rows = pocket_rows = missing = duplicates = 0
    binding_pairs: set[str] = set(); examples = []
    for bucket, wanted in sorted(by_bucket.items()):
        placements = sorted(pid[3:] if pid.startswith("P3|") else pid for pid in wanted)
        bf = read_bucket(p3_output / "binding_residues", bucket,
                         ["ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id"],
                         ds.field("ligand_assembly_placement_id").isin(placements))
        bf["pair_id"] = "P3|" + bf["ligand_assembly_placement_id"].astype(str)
        bf = bf[bf["pair_id"].isin(wanted)]
        pf = read_bucket(p3_output / "pair_pocket_residues", bucket,
                         ["pair_id", "chain_instance_id", "protein_residue_id"], ds.field("pair_id").isin(sorted(wanted)))
        keys = list(map(tuple, bf[["pair_id", "chain_instance_id", "protein_residue_id"]].astype(str).to_numpy()))
        pocket_keys = set(map(tuple, pf[["pair_id", "chain_instance_id", "protein_residue_id"]].astype(str).to_numpy()))
        miss = [key for key in keys if key not in pocket_keys]
        binding_rows += len(bf); pocket_rows += len(pf); missing += len(miss)
        duplicates += int(bf.duplicated(["pair_id", "chain_instance_id", "protein_residue_id"]).sum())
        binding_pairs.update(bf["pair_id"])
        examples.extend(miss[:max(0, 10 - len(examples))])
    return {"step4_pair_count": len(pair_ids), "frozen_binding_residue_rows": binding_rows,
            "frozen_pocket_residue_rows": pocket_rows, "pairs_with_binding_residues": len(binding_pairs),
            "pairs_without_binding_residues": len(pair_ids - binding_pairs), "binding_not_in_pocket": missing,
            "duplicate_binding_keys": duplicates, "subset_pass": missing == 0, "examples": examples}


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
        output = {sid: {"xyz": [], "units": [], "seen": set()} for sid in wanted}
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
            x, y, z = float(v[12]), float(v[13]), float(v[14]); dedup = (atom_id, atom_name, x, y, z)
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


def residue_key(frame: pd.DataFrame) -> pd.Series:
    return (frame["label_seq_id"].map(clean) + "|" + frame["auth_seq_id"].map(clean) + "|" +
            frame["insertion_code"].map(normalized_ins) + "|" + frame["label_comp_id"].map(clean).str.upper())


def inverse_transform(xyz: np.ndarray, r: np.ndarray, t: np.ndarray) -> np.ndarray:
    return (xyz - t) @ np.linalg.inv(r).T


def load_binding_targets(bucket: int, pair_records: list[dict], cfg: dict) -> tuple[dict[str, list[dict]], list[dict], list[dict]]:
    p2 = Path(cfg["input"]["processing2_output"]); p3 = Path(cfg["input"]["processing3_output"])
    pairs = pd.DataFrame(pair_records); wanted = set(pairs["candidate_pair_id"].astype(str)); placements = sorted(pid[3:] for pid in wanted)
    binding = read_bucket(p3 / "binding_residues", bucket,
        ["ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id", "qualifying_atomic_contact_count"],
        ds.field("ligand_assembly_placement_id").isin(placements))
    binding["pair_id"] = "P3|" + binding["ligand_assembly_placement_id"].astype(str); binding = binding[binding["pair_id"].isin(wanted)]
    lig_cols = ["filter_2_ligand_assembly_placement_id", "r11", "r12", "r13", "r21", "r22", "r23", "r31", "r32", "r33", "t1", "t2", "t3"]
    ligand = read_bucket(p2 / "prepared_ligand_assembly_atoms", bucket, lig_cols,
                         ds.field("filter_2_ligand_assembly_placement_id").isin(placements))
    wanted_chains = sorted(set(binding["chain_instance_id"].dropna().astype(str)))
    rec_cols = ["filter_1_chain_instance_id", "model_id", "entity_id", "label_asym_id", "auth_asym_id",
                "label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id", "type_symbol", "source_atom_row_index",
                "Cartn_x", "Cartn_y", "Cartn_z"]
    receptor = read_bucket(p2 / "prepared_receptor_assembly_atoms", bucket, rec_cols,
                           ds.field("filter_1_chain_instance_id").isin(wanted_chains) if wanted_chains else None)
    receptor = receptor[~receptor["type_symbol"].astype(str).str.upper().isin(HEAVY_EXCLUDE)].copy()
    receptor["_residue_key"] = residue_key(receptor)
    rec_groups = {(a, b): v for (a, b), v in receptor.groupby(["filter_1_chain_instance_id", "_residue_key"], sort=False)}
    lig_groups = {str(k): v for k, v in ligand.groupby("filter_2_ligand_assembly_placement_id", sort=False)}
    targets: dict[str, list[dict]] = defaultdict(list); base_rows = []; errors = []
    for row in binding.to_dict("records"):
        pid = row["pair_id"]; placement = row["ligand_assembly_placement_id"]
        rg = rec_groups.get((row["chain_instance_id"], row["protein_residue_id"])); lg = lig_groups.get(str(placement))
        if rg is None or rg.empty or lg is None or lg.empty:
            errors.append({"pair_id": pid, "external_instance_id": "", "error_reason": "FROZEN_BINDING_RESIDUE_ATOMS_OR_TRANSFORM_MISSING"}); continue
        rg = rg.drop_duplicates(["filter_1_chain_instance_id", "source_atom_row_index"])
        first_l = lg.iloc[0]; rmat = np.array([[first_l[f"r{i}{j}"] for j in range(1,4)] for i in range(1,4)], float)
        trans = np.array([first_l["t1"], first_l["t2"], first_l["t3"]], float)
        xyz = inverse_transform(rg[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), rmat, trans)
        first = rg.iloc[0]; bid = f"{row['chain_instance_id']}|{row['protein_residue_id']}"
        meta = {"pair_id": pid, "pdb_id": pid[3:7].lower(), "binding_residue_id": bid,
                "binding_chain_instance_id": row["chain_instance_id"], "binding_protein_residue_id": row["protein_residue_id"],
                "binding_model_id": clean(first["model_id"]), "binding_entity_id": clean(first["entity_id"]),
                "binding_label_asym_id": clean(first["label_asym_id"]), "binding_auth_asym_id": clean(first["auth_asym_id"]),
                "binding_comp_id": clean(first["label_comp_id"]).upper(), "binding_label_seq_id": clean(first["label_seq_id"]),
                "binding_auth_seq_id": clean(first["auth_seq_id"]), "binding_ins_code": normalized_ins(first["insertion_code"]),
                "binding_heavy_atom_count": len(xyz)}
        targets[pid].append({"id": bid, "xyz": xyz, "meta": meta}); base_rows.append(meta)
    return targets, base_rows, errors


def metric_from_pairs(external_xyz: np.ndarray, unit_meta: list[tuple], targets: list[dict], method: str, cutoff: float) -> dict:
    if not len(external_xyz) or not targets:
        return {"contacts": [], "ext_indices": set(), "target_atoms": set(), "pair_count": 0,
                "min_distance": math.nan, "direct": False, "bridged": set(), "true_units": set(), "true_ext": set()}
    target_xyz = np.concatenate([x["xyz"] for x in targets], axis=0)
    atom_to_target = []
    for ti, target in enumerate(targets): atom_to_target.extend([(ti, ai) for ai in range(len(target["xyz"]))])
    ext_contacts: list[list[tuple[int, int, float]]] = [[] for _ in range(len(external_xyz))]
    if method == "optimized":
        tree = cKDTree(target_xyz); candidates = tree.query_ball_point(external_xyz, r=np.nextafter(cutoff, math.inf))
        nearest = float(np.min(tree.query(external_xyz, k=1, workers=1)[0]))
        for ei, indices in enumerate(candidates):
            for gi in indices:
                d = float(np.linalg.norm(external_xyz[ei] - target_xyz[int(gi)]))
                if d <= cutoff + EPSILON:
                    ti, ai = atom_to_target[int(gi)]; ext_contacts[ei].append((ti, ai, d))
    else:
        nearest = math.inf
        for start in range(0, len(external_xyz), 512):
            dist = np.linalg.norm(external_xyz[start:start+512, None, :] - target_xyz[None, :, :], axis=2)
            nearest = min(nearest, float(dist.min())); ii, jj = np.where(dist <= cutoff + EPSILON)
            for i, j in zip(ii, jj):
                ti, ai = atom_to_target[int(j)]; ext_contacts[start + int(i)].append((ti, ai, float(dist[i, j])))
    grouped: dict[tuple[str, int], dict] = {}; ext_indices = set(); target_atoms = set(); pair_count = 0
    for ei, contacts in enumerate(ext_contacts):
        if not contacts: continue
        ext_indices.add(ei); pair_count += len(contacts)
        for ti, ai, d in contacts:
            target_atoms.add((ti, ai)); meta = unit_meta[ei]; key = (meta[0], ti)
            item = grouped.setdefault(key, {"unit_meta": meta, "target": targets[ti], "external_atoms": set(),
                                            "binding_atoms": set(), "pairs": 0, "min": math.inf})
            item["external_atoms"].add(ei); item["binding_atoms"].add(ai); item["pairs"] += 1; item["min"] = min(item["min"], d)
    contacts = []
    for key in sorted(grouped):
        item = grouped[key]; item["direct"] = len(item["external_atoms"]) >= 2; contacts.append(item)
    true_contacts = [x for x in contacts if x["direct"]]
    return {"contacts": contacts, "ext_indices": ext_indices, "target_atoms": target_atoms, "pair_count": pair_count,
            "min_distance": nearest, "direct": bool(true_contacts),
            "bridged": {x["target"]["id"] for x in true_contacts}, "true_units": {x["unit_meta"][0] for x in true_contacts},
            "true_ext": set().union(*(x["external_atoms"] for x in true_contacts)) if true_contacts else set()}


def compare_metrics(a: dict, b: dict) -> bool:
    def signature(m):
        return sorted((x["unit_meta"][0], x["target"]["id"], len(x["external_atoms"]), x["pairs"],
                       len(x["binding_atoms"]), x["direct"]) for x in m["contacts"])
    return (len(a["ext_indices"]) == len(b["ext_indices"]) and len(a["target_atoms"]) == len(b["target_atoms"])
            and a["pair_count"] == b["pair_count"] and a["direct"] == b["direct"] and a["bridged"] == b["bridged"]
            and signature(a) == signature(b) and math.isclose(a["min_distance"], b["min_distance"], abs_tol=1e-10, rel_tol=1e-12))


def process_bucket(task: tuple) -> dict:
    bucket, records, pair_records, cfg, work_text, cutoff, reference, safe_prune = task
    started = time.time(); work = Path(work_text); work.mkdir(parents=True, exist_ok=True)
    targets, binding_base, target_errors = load_binding_targets(bucket, pair_records, cfg)
    by_pdb: dict[str, list[dict]] = defaultdict(list)
    for row in records: by_pdb[str(row["pdb_id"]).lower()].append(row)
    instances = []; contacts_out = []; refs = []; errors = list(target_errors); source_mismatch = 0
    binding_stats = {(x["pair_id"], x["binding_residue_id"]): {**x, "instances": set(), "units": set(), "atoms": set(), "min": math.inf} for x in binding_base}
    for pdb, prows in sorted(by_pdb.items()):
        active = [r for r in prows if (not safe_prune or truth(r["touches_pocket_6A"]))]
        sources = {}; parse_error = ""
        if active:
            sources, parse_error = parse_source_atoms(locate_cif(Path(cfg["input"]["mmcif_root"]), pdb), {str(r["source_object_id"]) for r in active})
        transformed_cache = {}
        for row in prows:
            pid = str(row["candidate_pair_id"]); sid = str(row["source_object_id"]); target = targets.get(pid, [])
            base = {"pair_id": pid, "pdb_id": pdb, "assembly_id": row["assembly_id"], "model_id": row["model_id"],
                    "external_instance_id": row["crystal_instance_id"], "external_instance_key": row["crystal_instance_key"],
                    "source_object_key": sid, "source_object_type": row["source_object_type"],
                    "symmetry_operation_id": row["symmetry_operation_id"], "cell_h": row["cell_h"], "cell_k": row["cell_k"], "cell_l": row["cell_l"],
                    "touches_ligand_6A": truth(row["touches_ligand_6A"]), "touches_pocket_6A": truth(row["touches_pocket_6A"]),
                    "step1_min_pocket_distance_A": row["min_distance_to_pocket"], "n_frozen_binding_residues": len(target)}
            error = ""; should_skip = safe_prune and not truth(row["touches_pocket_6A"])
            if not target: error = "FROZEN_BINDING_RESIDUE_ATOMS_MISSING"
            elif should_skip:
                metric = {"contacts": [], "ext_indices": set(), "target_atoms": set(), "pair_count": 0, "min_distance": math.nan,
                          "direct": False, "bridged": set(), "true_units": set(), "true_ext": set()}
            elif parse_error: error = "MMCIF_PARSE_ERROR: " + parse_error
            elif sid not in sources or not len(sources[sid]["xyz"]): error = "SOURCE_OBJECT_ATOMS_MISSING"
            elif len(sources[sid]["xyz"]) != int(float(row["source_heavy_atom_count"])):
                source_mismatch += 1; error = f"SOURCE_HEAVY_ATOM_COUNT_MISMATCH:{len(sources[sid]['xyz'])}!={row['source_heavy_atom_count']}"
            if error:
                metric = {"contacts": [], "ext_indices": set(), "target_atoms": set(), "pair_count": 0, "min_distance": math.nan,
                          "direct": False, "bridged": set(), "true_units": set(), "true_ext": set()}
                errors.append({"pair_id": pid, "external_instance_id": row["crystal_instance_id"], "error_reason": error})
            elif not should_skip:
                cache_key = str(row["crystal_instance_key"])
                if cache_key not in transformed_cache:
                    rmat = np.array([float(row[c]) for c in R_COLS]).reshape(3,3); trans = np.array([float(row[c]) for c in T_COLS])
                    transformed_cache[cache_key] = (sources[sid]["xyz"] @ rmat.T + trans, sources[sid]["units"])
                xyz, unit_meta = transformed_cache[cache_key]; metric = metric_from_pairs(xyz, unit_meta, target, "optimized", cutoff)
                if reference:
                    brute = metric_from_pairs(xyz, unit_meta, target, "bruteforce", cutoff); match = compare_metrics(metric, brute)
                    refs.append({"pair_id": pid, "external_instance_id": row["crystal_instance_id"],
                        "optimized_min_distance_A": metric["min_distance"], "bruteforce_min_distance_A": brute["min_distance"],
                        "optimized_atom_pairs_4A": metric["pair_count"], "bruteforce_atom_pairs_4A": brute["pair_count"],
                        "optimized_external_atoms_4A": len(metric["ext_indices"]), "bruteforce_external_atoms_4A": len(brute["ext_indices"]),
                        "optimized_contact_combinations_4A": sum(x["direct"] for x in metric["contacts"]),
                        "bruteforce_contact_combinations_4A": sum(x["direct"] for x in brute["contacts"]),
                        "optimized_bridged_binding_residues": len(metric["bridged"]), "bruteforce_bridged_binding_residues": len(brute["bridged"]),
                        "optimized_instance_contact_4A": metric["direct"], "bruteforce_instance_contact_4A": brute["direct"], "reference_match": match})
            instances.append({**base, "n_binding_residues_contacted_4A": len(metric["bridged"]),
                "n_external_contact_units_contacting_binding_residues_4A": len(metric["true_units"]),
                "n_external_heavy_atoms_in_binding_residue_contacts_4A": len(metric["true_ext"]),
                "min_external_binding_residue_distance_A": metric["min_distance"],
                "instance_binding_residue_crystal_contact_4A": metric["direct"], "step4_instance_status": "ERROR" if error else "SUCCESS", "error_reason": error})
            for item in metric["contacts"]:
                um = item["unit_meta"]; bm = item["target"]["meta"]
                contacts_out.append({"pair_id": pid, "pdb_id": pdb, "external_instance_id": row["crystal_instance_id"],
                    "contact_unit_id": um[0], "contact_unit_type": um[1], "external_model_id": um[2], "external_entity_id": um[3],
                    "external_label_asym_id": um[4], "external_auth_asym_id": um[5], "external_comp_id": um[6],
                    "external_label_seq_id": um[7], "external_auth_seq_id": um[8], "external_ins_code": um[9], **bm,
                    "n_external_heavy_atoms_within_4A": len(item["external_atoms"]), "n_atom_pairs_within_4A": item["pairs"],
                    "n_binding_residue_heavy_atoms_contacted": len(item["binding_atoms"]), "min_distance_A": item["min"],
                    "external_to_binding_residue_contact_4A": item["direct"]})
                if item["direct"]:
                    stat = binding_stats[(pid, bm["binding_residue_id"])]; iid = str(row["crystal_instance_id"])
                    stat["instances"].add(iid); stat["units"].add((iid, um[0])); stat["atoms"].update((iid, i) for i in item["external_atoms"])
                    stat["min"] = min(stat["min"], item["min"])
    bridged = []
    for key in sorted(binding_stats):
        x = binding_stats[key]; bridged.append({**{k:x[k] for k in BRIDGED_HEADER if k in x},
            "n_external_instances_contacting": len(x["instances"]), "n_external_contact_units_contacting": len(x["units"]),
            "n_external_heavy_atoms_contributing": len(x["atoms"]), "min_external_distance_A": x["min"] if x["instances"] else math.nan,
            "crystal_bridged_binding_residue": bool(x["instances"])})
    write_rows(work/"instances.tsv", INSTANCE_HEADER, instances); write_rows(work/"contacts.tsv", CONTACT_HEADER, contacts_out)
    write_rows(work/"bridged.tsv", BRIDGED_HEADER, bridged); write_rows(work/"reference.tsv", REFERENCE_HEADER, refs)
    write_rows(work/"errors.tsv", ["pair_id", "external_instance_id", "error_reason"], errors)
    return {"bucket": bucket, "input_instances": len(records), "output_instances": len(instances), "contact_rows": len(contacts_out),
            "binding_rows": len(bridged), "reference_rows": len(refs), "reference_mismatches": sum(not truth(x["reference_match"]) for x in refs),
            "errors": len(errors), "source_count_mismatches": source_mismatch, "runtime_seconds": time.time()-started}


def boundary_tests(cutoff: float) -> dict:
    target = [{"id":"b", "xyz":np.array([[0.,0.,0.]]), "meta":{}}]; results = {}
    for d in (3.9999, 4.0, 4.0001):
        xyz = np.array([[d,0,0],[d,0,0]]); meta = [("u","POLYMER_RESIDUE","1","1","A","A","GLY","1","1","")]*2
        results[f"distance_{d:.4f}"] = metric_from_pairs(xyz, meta, target, "optimized", cutoff)["direct"]
    same = metric_from_pairs(np.array([[3.5,0,0],[3.7,0,0]]), [("u","POLYMER_RESIDUE","1","1","A","A","TYR","100","100","")]*2, target, "optimized", cutoff)["direct"]
    separate = metric_from_pairs(np.array([[3.5,0,0],[3.7,0,0]]), [("u1","POLYMER_RESIDUE","1","1","A","A","TYR","100","100",""),("u2","POLYMER_RESIDUE","1","1","A","A","ASP","101","101","")], target, "optimized", cutoff)["direct"]
    results.update({"two_distinct_external_atoms_same_binding_atom_true": same, "two_external_units_one_atom_each_false": not separate,
                    "boundary_pass": results["distance_3.9999"] and results["distance_4.0000"] and not results["distance_4.0001"] and same and not separate})
    return results


def pair_output(step3_pairs: pd.DataFrame, instances: pd.DataFrame, bridged: pd.DataFrame, contacts: pd.DataFrame) -> pd.DataFrame:
    for col in ["n_external_contact_units_contacting_binding_residues_4A", "min_external_binding_residue_distance_A"]:
        instances[col] = pd.to_numeric(instances[col], errors="coerce")
    instances["instance_binding_residue_crystal_contact_4A"] = instances["instance_binding_residue_crystal_contact_4A"].map(truth)
    if len(bridged): bridged["crystal_bridged_binding_residue"] = bridged["crystal_bridged_binding_residue"].map(truth)
    bridged_groups = {pid: g for pid, g in bridged.groupby("pair_id", sort=False)} if len(bridged) else {}
    if len(contacts):
        true_frame = contacts[contacts["external_to_binding_residue_contact_4A"].map(truth)]
        contact_groups = {pid: g for pid, g in true_frame.groupby("pair_id", sort=False)}
    else:
        contact_groups = {}
    agg = []
    for pid, group in instances.groupby("pair_id", sort=False):
        bg = bridged_groups.get(pid, bridged.iloc[0:0]); pos = bg[bg["crystal_bridged_binding_residue"].map(truth)]
        true_contacts = contact_groups.get(pid, contacts.iloc[0:0])
        agg.append({"candidate_pair_id": pid, "binding_residue_count": len(bg), "n_crystal_bridged_binding_residues": len(pos),
            "fraction_binding_residues_crystal_bridged": len(pos)/len(bg) if len(bg) else math.nan,
            "n_external_instances_step4": len(group), "n_external_instances_contacting_binding_residues": int(group["instance_binding_residue_crystal_contact_4A"].sum()),
            "n_external_contact_units_contacting_binding_residues": len(true_contacts[["external_instance_id","contact_unit_id"]].drop_duplicates()) if len(true_contacts) else 0,
            "min_external_binding_residue_distance_A": group["min_external_binding_residue_distance_A"].min(),
            "pair_binding_residue_crystal_contact_4A": bool(len(pos)), "step4_error_count": int(group["step4_instance_status"].eq("ERROR").sum())})
    out = step3_pairs.merge(pd.DataFrame(agg), on="candidate_pair_id", how="left")
    ints = ["binding_residue_count","n_crystal_bridged_binding_residues","n_external_instances_step4","n_external_instances_contacting_binding_residues","n_external_contact_units_contacting_binding_residues","step4_error_count"]
    out[ints] = out[ints].fillna(0).astype(np.int64); out["pair_binding_residue_crystal_contact_4A"] = out["pair_binding_residue_crystal_contact_4A"].fillna(False).astype(bool)
    statuses=[]
    for row in out.itertuples():
        if row.step3_status == "UPSTREAM_NO_NEIGHBOR": status="UPSTREAM_NO_NEIGHBOR"
        elif row.step3_status == "UPSTREAM_NO_EXTERNAL_NEIGHBOR": status="UPSTREAM_NO_EXTERNAL_NEIGHBOR"
        elif row.step3_status == "SUCCESS_DIRECT_LIGAND_CONTACT": status="UPSTREAM_DIRECT_LIGAND_CONTACT_REJECT"
        elif row.step3_status == "BA_EQUIVALENCE_REVIEW": status="BA_EQUIVALENCE_REVIEW"
        elif row.step4_error_count: status="ERROR"
        elif row.pair_binding_residue_crystal_contact_4A: status="SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT"
        else: status="SUCCESS_NO_BINDING_RESIDUE_CRYSTAL_CONTACT"
        statuses.append(status)
    out["step4_status"] = statuses
    return out


def observed_controls(instances: pd.DataFrame, contacts: pd.DataFrame, bridged: pd.DataFrame, external: pd.DataFrame, pair_frame: pd.DataFrame) -> dict:
    for col in ["n_external_heavy_atoms_within_4A", "n_atom_pairs_within_4A"]:
        if col in contacts: contacts[col] = pd.to_numeric(contacts[col], errors="coerce").fillna(0).astype(int)
    if len(contacts): contacts["external_to_binding_residue_contact_4A"] = contacts["external_to_binding_residue_contact_4A"].map(truth)
    instances["instance_binding_residue_crystal_contact_4A"] = instances["instance_binding_residue_crystal_contact_4A"].map(truth)
    groups = contacts.groupby(["pair_id","external_instance_id","binding_residue_id"])["n_external_heavy_atoms_within_4A"].apply(list) if len(contacts) else pd.Series(dtype=object)
    per_binding_instances = contacts[contacts["external_to_binding_residue_contact_4A"]].groupby(["pair_id","binding_residue_id"])["external_instance_id"].nunique() if len(contacts) else pd.Series(dtype=int)
    per_instance_binding = contacts[contacts["external_to_binding_residue_contact_4A"]].groupby(["pair_id","external_instance_id"])["binding_residue_id"].nunique() if len(contacts) else pd.Series(dtype=int)
    repeated = external.groupby(["pdb_id","source_object_id"])["crystal_instance_key"].nunique()
    ordinary_negative = False
    if len(instances):
        d = pd.to_numeric(instances["step1_min_pocket_distance_A"], errors="coerce")
        ordinary_negative = bool((d.le(4.0+EPSILON) & ~instances["instance_binding_residue_crystal_contact_4A"]).any())
    return {
        "case1_one_atom_false": bool((contacts["n_external_heavy_atoms_within_4A"].eq(1) & ~contacts["external_to_binding_residue_contact_4A"]).any()) if len(contacts) else False,
        "case2_exactly_two_true": bool((contacts["n_external_heavy_atoms_within_4A"].eq(2) & contacts["external_to_binding_residue_contact_4A"]).any()) if len(contacts) else False,
        "case3_many_true": bool((contacts["n_external_heavy_atoms_within_4A"].ge(3) & contacts["external_to_binding_residue_contact_4A"]).any()) if len(contacts) else False,
        "case4_two_external_units_one_atom_each_false": any(len(x)>=2 and all(v==1 for v in x) for x in groups),
        "case5_true_and_one_atom_units": any(any(v>=2 for v in x) and any(v==1 for v in x) for x in groups),
        "case9_ordinary_pocket_negative": ordinary_negative,
        "case10_binding_residue_path_true": bool(instances["instance_binding_residue_crystal_contact_4A"].any()),
        "case11_instance_contacts_multiple_binding_residues": bool((per_instance_binding>=2).any()),
        "case12_multiple_instances_contact_same_binding_residue": bool((per_binding_instances>=2).any()),
        "case13_polymer_external_positive": bool(((instances["source_object_type"]=="POLYMER") & instances["instance_binding_residue_crystal_contact_4A"]).any()),
        "case14_nonpolymer_external_positive": bool(((instances["source_object_type"]=="NONPOLYMER") & instances["instance_binding_residue_crystal_contact_4A"]).any()),
        "case15_step3_direct_excluded": not set(external["candidate_pair_id"]) & set(pair_frame.loc[pair_frame["step4_status"].eq("UPSTREAM_DIRECT_LIGAND_CONTACT_REJECT"),"candidate_pair_id"]),
        "case16_only_external_instances_entered": external["equivalence_status"].eq("EXTERNAL_CRYSTAL_INSTANCE").all(),
        "case17_ba_review_propagated": int(pair_frame["step4_status"].eq("BA_EQUIVALENCE_REVIEW").sum()) == int(pair_frame["step3_status"].eq("BA_EQUIVALENCE_REVIEW").sum()),
        "same_source_different_placements_observed": bool((repeated>1).any()),
    }


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True,type=Path); ap.add_argument("--run-dir",required=True,type=Path)
    ap.add_argument("--mode",choices=["smoke","pilot","full"],required=True); ap.add_argument("--pilot-count",type=int,default=500)
    args=ap.parse_args(); started=time.time(); cfg=yaml.safe_load(args.config.read_text()); run=args.run_dir
    if run.exists() and any(run.iterdir()): raise SystemExit(f"run directory not empty: {run}")
    for rel in ("output","work"): (run/rel).mkdir(parents=True,exist_ok=True)
    shutil.copy2(args.config,run/"config_snapshot.yaml"); shutil.copy2(Path(__file__),run/"executed_runner.py")
    step3=Path(cfg["input"]["step3_run"]); step2=Path(cfg["input"]["step2_run"]); step1_work=Path(cfg["input"]["step1_work"])
    if json.loads((step3/"_FROZEN.json").read_text()).get("status")!="FROZEN": raise RuntimeError("Step3 input not frozen")
    p3_run=Path(cfg["input"]["processing3_output"]).parent
    p2_run=Path(cfg["input"]["processing2_output"]).parent
    if json.loads((p3_run/"_FROZEN.json").read_text()).get("status")!="FROZEN": raise RuntimeError("Processing3 binding source not frozen")
    if not (p2_run/"_FROZEN.json").exists(): raise RuntimeError("Processing2 coordinate source not frozen")
    step3_pairs_all=pd.read_csv(step3/"output/03_pair_direct_contact_inventory.tsv.gz",sep="\t",dtype=str,keep_default_na=False)
    step4_pair_ids=set(step3_pairs_all.loc[step3_pairs_all["step3_status"].eq("SUCCESS_NO_DIRECT_LIGAND_CONTACT"),"candidate_pair_id"])
    external_all=pd.read_csv(step2/"output/03_external_crystal_instances.tsv.gz",sep="\t",dtype=str,keep_default_na=False)
    external_all=external_all[external_all["candidate_pair_id"].isin(step4_pair_ids)].copy()
    pair_bucket=read_pair_buckets(step1_work); external_all["bucket_id"]=external_all["candidate_pair_id"].map(pair_bucket)
    if external_all["bucket_id"].isna().any(): raise RuntimeError("Step4 external pair missing Step1 bucket")
    external=external_all; pairs=step3_pairs_all; selection=None
    if args.mode in {"smoke","pilot"}:
        count=args.pilot_count if args.mode=="pilot" else min(args.pilot_count,30)
        selection=select_pilot(external_all,step3_pairs_all,count); keep=set(selection["pdb_id"])
        external=external_all[external_all["pdb_id"].isin(keep)].copy(); pairs=step3_pairs_all[step3_pairs_all["pdb_id"].isin(keep)].copy()
        selection.to_csv(run/"pilot_selection.tsv",sep="\t",index=False)
    scientific_pair_ids=set(external["candidate_pair_id"]); subset=audit_binding_subset(scientific_pair_ids,pair_bucket,Path(cfg["input"]["processing3_output"]))
    atomic_json(run/"binding_residue_subset_preflight.json",subset); safe_prune=bool(subset["subset_pass"] and cfg["policy"].get("use_pocket6_pruning_if_subset_pass",True))
    tasks=[]
    pair_records={pid:row for pid,row in ((r["candidate_pair_id"],r) for r in pairs.to_dict("records"))}
    for bucket,group in external.groupby("bucket_id",sort=True):
        b=int(bucket); pids=sorted(set(group["candidate_pair_id"])); tasks.append((b,group.to_dict("records"),[pair_records[x] for x in pids],cfg,str(run/"work"/f"bucket_{b:03d}"),float(cfg["contact"]["cutoff_angstrom"]),args.mode in {"smoke","pilot"},safe_prune))
    results=[]
    with cf.ProcessPoolExecutor(max_workers=int(cfg["runtime"]["workers"])) as pool:
        futures={pool.submit(process_bucket,t):t[0] for t in tasks}
        for future in cf.as_completed(futures):
            results.append(future.result()); atomic_json(run/"progress.json",{"completed_buckets":len(results),"total_buckets":len(tasks),
                "instances_completed":sum(x["input_instances"] for x in results),"errors":sum(x["errors"] for x in results),"updated_at":utc()})
    works=sorted((run/"work").glob("bucket_*")); output=run/"output"
    instances=concat_tsv([w/"instances.tsv" for w in works],output/"01_external_instance_binding_residue_contact.tsv.gz",INSTANCE_HEADER)
    contacts=concat_tsv([w/"contacts.tsv" for w in works],output/"02_external_contact_unit_binding_residue_contacts.tsv.gz",CONTACT_HEADER)
    bridged=concat_tsv([w/"bridged.tsv" for w in works],output/"03_crystal_bridged_binding_residues.tsv.gz",BRIDGED_HEADER)
    refs=concat_tsv([w/"reference.tsv" for w in works],output/"pilot_optimized_vs_bruteforce.tsv.gz",REFERENCE_HEADER) if args.mode in {"smoke","pilot"} else pd.DataFrame(columns=REFERENCE_HEADER)
    pair_frame=pair_output(pairs,instances.copy(),bridged.copy(),contacts.copy()); write_gzip_frame(pair_frame,output/"04_pair_binding_residue_contact_inventory.tsv.gz")
    boundaries=boundary_tests(float(cfg["contact"]["cutoff_angstrom"])); controls=observed_controls(instances.copy(),contacts.copy(),bridged.copy(),external,pair_frame)
    status_counts=Counter(pair_frame["step4_status"]); true_instances=int(instances["instance_binding_residue_crystal_contact_4A"].map(truth).sum())
    true_contacts=int(contacts["external_to_binding_residue_contact_4A"].map(truth).sum()) if len(contacts) else 0
    true_binding=int(bridged["crystal_bridged_binding_residue"].map(truth).sum()) if len(bridged) else 0
    checks={"step4_input_pair_accounting":len(scientific_pair_ids)==len(set(external["candidate_pair_id"])),"instance_accounting":len(instances)==len(external),
        "silent_drop_zero":len(instances)==len(external),"duplicate_external_instance_key_zero":not instances.duplicated(["pair_id","external_instance_key"]).any(),
        "binding_subset_pocket":subset["subset_pass"],"duplicate_frozen_binding_keys_zero":subset["duplicate_binding_keys"]==0,
        "frozen_binding_rows_accounting":len(bridged)==subset["frozen_binding_residue_rows"],"errors_zero":sum(x["errors"] for x in results)==0,
        "source_heavy_atom_count_mismatch_zero":sum(x["source_count_mismatches"] for x in results)==0,"boundary_tests_pass":boundaries["boundary_pass"],
        "reference_mismatch_zero":sum(x["reference_mismatches"] for x in results)==0,
        "contact_true_iff_two_external_atoms":bool((contacts["external_to_binding_residue_contact_4A"].map(truth)==pd.to_numeric(contacts["n_external_heavy_atoms_within_4A"]).ge(2)).all()) if len(contacts) else True,
        "contact_minimum_within_4A":bool(pd.to_numeric(contacts["min_distance_A"],errors="coerce").le(float(cfg["contact"]["cutoff_angstrom"])+EPSILON).all()) if len(contacts) else True,
        "duplicate_contact_combination_rows_zero":not contacts.duplicated(["pair_id","external_instance_id","contact_unit_id","binding_residue_id"]).any() if len(contacts) else True,
        "bridged_iff_positive_contact":bool((bridged["crystal_bridged_binding_residue"].map(truth)==pd.to_numeric(bridged["n_external_contact_units_contacting"]).gt(0)).all()) if len(bridged) else True,
        "instance_iff_bridged_residue":bool((instances["instance_binding_residue_crystal_contact_4A"].map(truth)==pd.to_numeric(instances["n_binding_residues_contacted_4A"]).gt(0)).all()),
        "non_pocket6_pruned_zero":bool((pd.to_numeric(instances.loc[~instances["touches_pocket_6A"].map(truth),"n_binding_residues_contacted_4A"]).eq(0)).all()) if safe_prune else True}
    if args.mode in {"pilot","full"}: checks["all_required_scientific_controls_present"]=all(controls.values())
    if args.mode=="full":
        exp=cfg["validation"]; checks.update({"full_pairs_336412":len(pair_frame)==int(exp["pairs"]),"step4_pairs_139403":len(scientific_pair_ids)==int(exp["step4_pairs"]),
            "direct_excluded_57580":status_counts["UPSTREAM_DIRECT_LIGAND_CONTACT_REJECT"]==int(exp["direct_pairs"]),
            "upstream_no_neighbor_7663":status_counts["UPSTREAM_NO_NEIGHBOR"]==int(exp["upstream_no_neighbor"]),
            "upstream_no_external_131764":status_counts["UPSTREAM_NO_EXTERNAL_NEIGHBOR"]==int(exp["upstream_no_external"]),
            "ba_review_2":status_counts["BA_EQUIVALENCE_REVIEW"]==int(exp["ba_review"]),
            "step4_pair_partition":status_counts["SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT"]+status_counts["SUCCESS_NO_BINDING_RESIDUE_CRYSTAL_CONTACT"]==int(exp["step4_pairs"])})
    validation={"run_mode":args.mode,"validated_at":utc(),"validation_pass":all(checks.values()),"checks":checks,"boundary_tests":boundaries,
        "binding_residue_subset_preflight":subset,"pocket6_pruning_enabled":safe_prune,"observed_controls":controls,
        "counts":{"pairs":len(pair_frame),"step4_pairs":len(scientific_pair_ids),"external_instances":len(instances),"contact_candidate_rows":len(contacts),
                  "true_contact_combinations":true_contacts,"frozen_binding_residues":len(bridged),"crystal_bridged_binding_residues":true_binding,
                  "contacting_instances":true_instances,"reference_rows":len(refs),"errors":sum(x["errors"] for x in results)},"pair_status_counts":dict(status_counts)}
    summary=pd.DataFrame([{**validation["counts"],"pdb_count":pair_frame["pdb_id"].nunique(),"step4_positive_pairs":status_counts["SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT"],
        "step4_negative_pairs":status_counts["SUCCESS_NO_BINDING_RESIDUE_CRYSTAL_CONTACT"],"polymer_contacting_instances":int(((instances["source_object_type"]=="POLYMER")&instances["instance_binding_residue_crystal_contact_4A"].map(truth)).sum()),
        "nonpolymer_contacting_instances":int(((instances["source_object_type"]=="NONPOLYMER")&instances["instance_binding_residue_crystal_contact_4A"].map(truth)).sum()),
        "one_atom_only_contact_combinations":int(pd.to_numeric(contacts["n_external_heavy_atoms_within_4A"]).eq(1).sum()) if len(contacts) else 0,
        "runtime_seconds":time.time()-started,"validation_pass":validation["validation_pass"]}])
    summary.to_csv(output/"05_step4_summary.tsv",sep="\t",index=False)
    report=f"""# Filter 4 Step 4 — Binding-Residue-Mediated Crystal Contact Analysis\n\nRun: `{run.name}`  \nMode: `{args.mode}`  \nStatus: `{'PASS' if validation['validation_pass'] else 'FAIL'}`\n\nFilter 4 Step 4 identifies the shortest residue-mediated structural pathway by which an external crystallographic neighbour may constrain the ligand-binding environment without directly contacting the ligand. Analysis is restricted to frozen ligand-binding residues. An external polymer residue or independent non-polymer component is considered to contact a frozen binding residue when at least two distinct heavy atoms from the same external contact unit lie within 4.0 Å of any heavy atom of that binding residue. A binding residue satisfying this condition is annotated as a crystal-bridged binding residue.\n\nNo propagation through non-binding pocket residues or multi-residue interaction networks is performed. Step 4 detects structural crystal-packing risk and does not by itself prove that the crystal neighbour caused a conformational change. No Step 4 exclusion policy is assigned.\n\n- complete Filter 4 pair inventory: {len(pair_frame):,}\n- Step 4 scientific pairs: {len(scientific_pair_ids):,}\n- Step 4 external instances: {len(instances):,}\n- positive pairs: {status_counts['SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT']:,}\n- crystal-bridged frozen binding residues: {true_binding:,}\n- validation: {'PASS' if validation['validation_pass'] else 'FAIL'}\n"""
    (output/"06_step4_report.md").write_text(report,encoding="utf-8"); atomic_json(run/"validation.json",validation)
    atomic_json(run/"input_provenance.json",{"step3_run":str(step3),"step3_sha256sums":sha256(step3/"SHA256SUMS"),"step2_run":str(step2),
        "processing3_frozen_binding_residues":str(Path(cfg["input"]["processing3_output"])/"binding_residues"),
        "processing3_frozen_pocket_residues":str(Path(cfg["input"]["processing3_output"])/"pair_pocket_residues"),
        "processing2_frozen_assembly_atoms":cfg["input"]["processing2_output"],"policy":"One-hop external contact unit to frozen binding residue only; no binding membership, pocket, lattice, BA, direct ligand, multi-hop, severity, or exclusion recomputation","created_at":utc()})
    datasets={"01_external_instance_binding_residue_contact":instances,"02_external_contact_unit_binding_residue_contacts":contacts,
              "03_crystal_bridged_binding_residues":bridged,"04_pair_binding_residue_contact_inventory":pair_frame,"05_step4_summary":summary}
    atomic_json(run/"output_schema.json",{"schema_version":SCHEMA_VERSION,"datasets":{k:[{"column_name":c,"data_type":str(f[c].dtype),"nullable":bool(f[c].isna().any())} for c in f.columns] for k,f in datasets.items()}})
    serialized={}
    for path,n in [(output/"01_external_instance_binding_residue_contact.tsv.gz",len(instances)),(output/"02_external_contact_unit_binding_residue_contacts.tsv.gz",len(contacts)),(output/"03_crystal_bridged_binding_residues.tsv.gz",len(bridged)),(output/"04_pair_binding_residue_contact_inventory.tsv.gz",len(pair_frame))]:
        with gzip.open(path,"rt",encoding="utf-8",newline="") as fh: physical=max(0,sum(1 for _ in fh)-1)
        serialized[path.name]={"expected_rows":n,"physical_rows":physical,"match":physical==n}
    validation["serialization_row_counts"]=serialized; validation["checks"]["serialized_physical_row_counts_match"]=all(x["match"] for x in serialized.values()); validation["validation_pass"]=all(validation["checks"].values())
    atomic_json(run/"validation.json",validation); summary.loc[0,"validation_pass"]=validation["validation_pass"]; summary.to_csv(output/"05_step4_summary.tsv",sep="\t",index=False)
    files=[p for p in run.rglob("*") if p.is_file() and "work" not in p.parts and p.name not in {"SHA256SUMS","output_manifest.tsv","_FROZEN.json"}]
    manifest=[{"relative_path":p.relative_to(run).as_posix(),"size_bytes":p.stat().st_size,"sha256":sha256(p),"schema_version":SCHEMA_VERSION,"generated_by":VERSION} for p in sorted(files)]
    pd.DataFrame(manifest).to_csv(run/"output_manifest.tsv",sep="\t",index=False)
    checksum_files=[p for p in run.rglob("*") if p.is_file() and "work" not in p.parts and p.name not in {"SHA256SUMS","_FROZEN.json"}]
    with (run/"SHA256SUMS").open("w",encoding="utf-8") as fh:
        for p in sorted(checksum_files): fh.write(f"{sha256(p)}  {p.relative_to(run).as_posix()}\n")
    if validation["validation_pass"]: atomic_json(run/"_FROZEN.json",{"status":"FROZEN","run_id":run.name,"stage":VERSION,"validation_pass":True,"frozen_at":utc(),"sha256sums_sha256":sha256(run/"SHA256SUMS")})
    print(json.dumps(validation,indent=2,default=lambda x:x.item() if hasattr(x,"item") else str(x)))
    if not validation["validation_pass"]: raise SystemExit(2)


if __name__ == "__main__":
    main()
