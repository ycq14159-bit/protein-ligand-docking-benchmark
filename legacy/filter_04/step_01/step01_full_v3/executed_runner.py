#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import yaml


PAIR_HEADER = [
    "candidate_pair_id", "pdb_id", "assembly_id", "model_id",
    "target_frame_status", "trivial_source_object_ids",
]
ATOM_HEADER = ["candidate_pair_id", "target_kind", "x", "y", "z"]
FRAME_HEADER = [
    "candidate_pair_id", "pdb_id", "assembly_id", "model_id",
    "ligand_assembly_placement_id", "ligand_atom_count", "pocket_atom_count",
    "inverse_transform_rmsd", "max_atom_deviation", "target_frame_status",
    "trivial_source_object_ids", "error_reason",
]
INVENTORY_HEADER = [
    "candidate_pair_id", "pdb_id", "step1_status", "target_frame_status",
    "n_source_objects", "n_raw_lattice_candidates", "n_bbox_pass_candidates",
    "n_unique_neighbor_instances", "n_ligand_6A_instances", "n_pocket_6A_instances",
    "n_ligand_and_pocket_6A_instances", "n_trivial_self_instances_skipped",
    "runtime_ms", "error_reason",
]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_sha256s(run: Path) -> None:
    with (run / "SHA256SUMS").open("w") as fh:
        for path in sorted(p for p in run.iterdir() if p.is_file() and p.name != "SHA256SUMS"):
            fh.write(f"{sha256(path)}  {path.name}\n")


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True,
                              default=lambda x: x.item() if hasattr(x, "item") else str(x)) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def clean(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text in {"False", "None", "nan", ".", "?"} else text


def bucket_path(root: Path, bucket: int) -> Path:
    plain = root / f"bucket_id={bucket}"
    padded = root / f"bucket_id={bucket:03d}"
    return plain if plain.exists() else padded


def read_bucket(root: Path, bucket: int, columns: list[str], filter_expression=None) -> pd.DataFrame:
    path = bucket_path(root, bucket)
    if not path.exists():
        return pd.DataFrame(columns=columns)
    return ds.dataset(path, format="parquet").to_table(
        columns=columns, filter=filter_expression
    ).to_pandas(split_blocks=True, self_destruct=True)


def retained_pairs(config: dict) -> pd.DataFrame:
    root = Path(config["input"]["filter3_dataset"])
    pieces = []
    columns = ["pair_id", "ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id", "filter3_v2_terminal_status", "bucket_id"]
    keep = set(config["input"]["retained_statuses"])
    for bucket in range(256):
        frame = read_bucket(root, bucket, columns)
        if not frame.empty:
            pieces.append(frame[frame["filter3_v2_terminal_status"].isin(keep)])
    frame = pd.concat(pieces, ignore_index=True)
    frame["pdb_id"] = frame["pdb_id"].astype(str).str.lower()
    if len(frame) != 336_412 or frame["pair_id"].nunique() != len(frame):
        raise RuntimeError(f"retained Filter3 invariant failed rows={len(frame)} unique={frame['pair_id'].nunique()}")
    return frame.sort_values(["bucket_id", "pdb_id", "pair_id"], kind="stable")


def deterministic_rank(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def select_pilot(pairs: pd.DataFrame, config: dict, run: Path) -> pd.DataFrame:
    count = int(config["pilot"]["pdb_count"])
    audit = pd.read_csv(Path(config["input"]["step0b_audit"]) / "01_pdb_crystal_metadata_audit.tsv.gz", sep="\t", compression="gzip")
    audit["pdb_id"] = audit["pdb_id"].astype(str).str.lower()
    audit = audit.drop(columns=["pair_count"], errors="ignore")
    stats = pairs.groupby("pdb_id", as_index=False).agg(pair_count=("pair_id", "size"), bucket_id=("bucket_id", "first"))
    stats = stats.merge(audit, on="pdb_id", how="left", validate="one_to_one")
    op_count = stats.get("explicit_symops_count", pd.Series(0, index=stats.index)).fillna(0).astype(int)
    stats["symop_stratum"] = pd.cut(op_count, [-1, 1, 4, 8, 10**9], labels=["identity_or_implicit", "2_4", "5_8", "9_plus"]).astype(str)
    a, b, c = (stats[x].astype(float) for x in ["gemmi_cell_a", "gemmi_cell_b", "gemmi_cell_c"])
    stats["cell_elongation"] = pd.concat([a, b, c], axis=1).max(axis=1) / pd.concat([a, b, c], axis=1).min(axis=1)
    stats["cell_stratum"] = pd.qcut(stats["cell_elongation"].rank(method="first"), 4, labels=["q1", "q2", "q3", "q4"]).astype(str)
    non_identity = pairs.assign(non_identity=~pairs["ligand_assembly_placement_id"].astype(str).str.endswith("|1")).groupby("pdb_id")["non_identity"].any()
    stats["ligand_operator_stratum"] = stats["pdb_id"].map(non_identity).fillna(False).map({True: "non_identity", False: "identity"})
    stats["pair_stratum"] = pd.cut(stats["pair_count"], [0, 1, 3, 10, 10**9], labels=["one", "2_3", "4_10", "11_plus"]).astype(str)
    stats["stratum"] = stats[["symop_stratum", "cell_stratum", "ligand_operator_stratum", "pair_stratum"]].agg("|".join, axis=1)
    stats["rank"] = stats["pdb_id"].map(deterministic_rank)
    known = ["2d3e", "4jgc", "4jbm", "4k3t", "5m3h", "1was", "3hz3", "1wqj", "2de3", "1jcd", "3bd3", "1men", "2gkp", "1wui", "2ahf", "2h2z"]
    selected = []
    for pdb in known:
        if pdb in set(stats["pdb_id"]) and len(selected) < count:
            selected.append((pdb, "BIOJAVA_CRYSTALBUILDER_REFERENCE_CASE"))
    candidates = stats[~stats["pdb_id"].isin({p for p, _ in selected})].sort_values(["stratum", "rank"])
    groups = list(candidates.groupby("stratum", sort=True))
    cursor = {name: 0 for name, _ in groups}
    while len(selected) < count:
        advanced = False
        for name, group in groups:
            i = cursor[name]
            if i < len(group) and len(selected) < count:
                selected.append((group.iloc[i]["pdb_id"], "DETERMINISTIC_STRATIFIED:" + name))
                cursor[name] += 1; advanced = True
        if not advanced: break
    selection = stats[stats["pdb_id"].isin({p for p, _ in selected})].copy()
    reason = dict(selected)
    selection["selection_reason"] = selection["pdb_id"].map(reason)
    selection["biojava_reference_selected"] = False
    reference_count = int(config["pilot"]["reference_pdb_count"])
    refs = selection.sort_values(["symop_stratum", "rank"])["pdb_id"].head(reference_count)
    selection.loc[selection["pdb_id"].isin(refs), "biojava_reference_selected"] = True
    selection.sort_values("pdb_id").to_csv(run / "pilot_selection.tsv", sep="\t", index=False)
    return selection.sort_values("pdb_id")


def source_object_id(row, object_type: str) -> str:
    model = clean(row.get("model_id", "1")) or "1"
    entity = clean(row.get("entity_id", row.get("label_entity_id", "")))
    label = clean(row.get("label_asym_id", "")); auth = clean(row.get("auth_asym_id", ""))
    if object_type == "POLYMER":
        return f"POLYMER|{model}|{entity}|{label}|{auth}"
    comp = clean(row.get("component_id", "")).upper(); auth_seq = clean(row.get("auth_seq_id", ""))
    label_seq = clean(row.get("label_seq_id", "")); ins = clean(row.get("insertion_code", ""))
    return f"NONPOLYMER|{model}|{entity}|{label}|{auth}|{comp}|{auth_seq}|{label_seq}|{ins}"


def residue_key(frame: pd.DataFrame) -> pd.Series:
    return (frame["label_seq_id"].map(clean) + "|" + frame["auth_seq_id"].map(clean) + "|" +
            frame["insertion_code"].map(clean) + "|" + frame["label_comp_id"].map(clean).str.upper())


def inverse_transform(xyz: np.ndarray, r: np.ndarray, t: np.ndarray) -> np.ndarray:
    # PDBx assembly operators are decimal-valued and therefore not guaranteed
    # to be exactly orthonormal in floating point.  Apply the actual inverse
    # instead of substituting R.T, which can accumulate millangstrom errors for
    # large coordinates/translations.
    return (xyz - t) @ np.linalg.inv(r).T


def prepare_bucket(task: tuple[int, list[dict], dict, str, bool, list[str]]) -> dict:
    bucket, pair_rows, config, run_text, compare, reference_pdbs = task
    run = Path(run_text); work = run / "work" / f"bucket_{bucket:03d}"; marker = work / "_SUCCESS.json"
    if marker.exists(): return json.loads(marker.read_text())
    work.mkdir(parents=True, exist_ok=True)
    pairs = pd.DataFrame(pair_rows)
    p2 = Path(config["input"]["processing2_output"]); p3 = Path(config["input"]["processing3_output"])
    wanted_pairs = set(pairs["pair_id"]); wanted_placements = set(pairs["ligand_assembly_placement_id"])
    pocket = read_bucket(
        p3 / "pair_pocket_residues", bucket,
        ["pair_id", "chain_instance_id", "protein_residue_id"],
        ds.field("pair_id").isin(sorted(wanted_pairs)),
    )
    lig_cols = ["filter_2_ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id", "entity_id", "label_asym_id", "auth_asym_id", "label_seq_id", "auth_seq_id", "insertion_code", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z", "source_Cartn_x", "source_Cartn_y", "source_Cartn_z", "r11", "r12", "r13", "r21", "r22", "r23", "r31", "r32", "r33", "t1", "t2", "t3"]
    rec_cols = ["filter_1_chain_instance_id", "pdb_id", "assembly_id", "model_id", "entity_id", "label_asym_id", "auth_asym_id", "label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z", "source_Cartn_x", "source_Cartn_y", "source_Cartn_z"]
    ligand = read_bucket(
        p2 / "prepared_ligand_assembly_atoms", bucket, lig_cols,
        ds.field("filter_2_ligand_assembly_placement_id").isin(sorted(wanted_placements)),
    )
    wanted_chains = sorted(set(pocket["chain_instance_id"].dropna().astype(str)))
    receptor = read_bucket(
        p2 / "prepared_receptor_assembly_atoms", bucket, rec_cols,
        ds.field("filter_1_chain_instance_id").isin(wanted_chains) if wanted_chains else None,
    )
    pocket = pocket[pocket["pair_id"].isin(wanted_pairs)]
    ligand = ligand[ligand["filter_2_ligand_assembly_placement_id"].isin(wanted_placements)]
    receptor["_residue_key"] = residue_key(receptor)
    heavy_ligand = ligand[~ligand["type_symbol"].astype(str).str.upper().isin(["H", "D", "T"])]
    heavy_receptor = receptor[~receptor["type_symbol"].astype(str).str.upper().isin(["H", "D", "T"])]
    lig_groups = {k: v for k, v in heavy_ligand.groupby("filter_2_ligand_assembly_placement_id", sort=False)}
    rec_groups = {(a, b): v for (a, b), v in heavy_receptor.groupby(["filter_1_chain_instance_id", "_residue_key"], sort=False)}
    pocket_groups = {k: v for k, v in pocket.groupby("pair_id", sort=False)}
    tol = float(config["policy"]["target_frame_rmsd_tolerance"])
    frame_rows, skipped = [], []
    with (work / "pairs.tsv").open("w", newline="", encoding="utf-8") as pf, (work / "target_atoms.tsv").open("w", newline="", encoding="utf-8") as af:
        pw, aw = csv.writer(pf, delimiter="\t", lineterminator="\n"), csv.writer(af, delimiter="\t", lineterminator="\n")
        pw.writerow(PAIR_HEADER); aw.writerow(ATOM_HEADER)
        for pair in pairs.to_dict("records"):
            start = time.monotonic(); pair_id = pair["pair_id"]; placement = pair["ligand_assembly_placement_id"]
            error = ""; status = "TARGET_FRAME_PASS"; trivial = []
            lg = lig_groups.get(placement); pg = pocket_groups.get(pair_id)
            if lg is None or lg.empty or pg is None or pg.empty:
                status = "TARGET_FRAME_REVIEW"; error = "missing ligand or frozen pocket rows"
                lig_xyz = pocket_xyz = np.empty((0, 3)); rmsd = maxdev = math.nan
            else:
                first = lg.iloc[0]
                r = np.array([[first[f"r{i}{j}"] for j in range(1, 4)] for i in range(1, 4)], dtype=float)
                t = np.array([first["t1"], first["t2"], first["t3"]], dtype=float)
                lig_assembly = lg[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
                lig_xyz = inverse_transform(lig_assembly, r, t)
                lig_source = lg[["source_Cartn_x", "source_Cartn_y", "source_Cartn_z"]].to_numpy(float)
                deviations = np.linalg.norm(lig_xyz - lig_source, axis=1)
                rmsd = float(np.sqrt(np.mean(np.square(deviations)))); maxdev = float(deviations.max())
                if rmsd > tol or maxdev > tol:
                    status = "TARGET_FRAME_REVIEW"; error = f"ligand inverse transform deviation rmsd={rmsd:.9g} max={maxdev:.9g}"
                trivial.append(source_object_id(first, "NONPOLYMER"))
                pocket_parts = []
                chain_parts = defaultdict(list)
                for pr in pg.to_dict("records"):
                    rg = rec_groups.get((pr["chain_instance_id"], pr["protein_residue_id"]))
                    if rg is not None: pocket_parts.append(rg); chain_parts[pr["chain_instance_id"]].append(rg)
                if not pocket_parts:
                    status = "TARGET_FRAME_REVIEW"; error = "no receptor atoms matched frozen pocket residues"; pocket_xyz = np.empty((0, 3))
                else:
                    pocket_frame = pd.concat(pocket_parts, ignore_index=True).drop_duplicates(["filter_1_chain_instance_id", "source_atom_row_index"] if "source_atom_row_index" in pocket_parts[0].columns else ["filter_1_chain_instance_id", "Cartn_x", "Cartn_y", "Cartn_z"])
                    pocket_xyz = inverse_transform(pocket_frame[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), r, t)
                    for chain_id, parts in chain_parts.items():
                        cf = pd.concat(parts, ignore_index=True).drop_duplicates(["Cartn_x", "Cartn_y", "Cartn_z"])
                        canonical = inverse_transform(cf[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), r, t)
                        source = cf[["source_Cartn_x", "source_Cartn_y", "source_Cartn_z"]].to_numpy(float)
                        if len(source) and float(np.max(np.linalg.norm(canonical - source, axis=1))) <= tol:
                            trivial.append(source_object_id(cf.iloc[0], "POLYMER"))
            trivial = sorted(set(trivial))
            frame_rows.append([pair_id, pair["pdb_id"], pair["assembly_id"], pair["model_id"], placement, len(lig_xyz), len(pocket_xyz), rmsd, maxdev, status, ";".join(trivial), error])
            if status != "TARGET_FRAME_PASS":
                skipped.append([pair_id, pair["pdb_id"], "TARGET_FRAME_REVIEW", status, 0, 0, 0, 0, 0, 0, 0, 0, int((time.monotonic()-start)*1000), error])
                continue
            pw.writerow([pair_id, pair["pdb_id"], pair["assembly_id"], pair["model_id"], status, ";".join(trivial)])
            for xyz in lig_xyz: aw.writerow([pair_id, "LIGAND", *xyz])
            for xyz in pocket_xyz: aw.writerow([pair_id, "POCKET", *xyz])
    pd.DataFrame(frame_rows, columns=FRAME_HEADER).to_csv(work / "target_frame_validation.tsv", sep="\t", index=False)
    pd.DataFrame(skipped, columns=INVENTORY_HEADER).to_csv(work / "pre_skipped_inventory.tsv", sep="\t", index=False)
    runnable = len(frame_rows) - len(skipped)
    if runnable:
        ref_file = work / "reference_pdbs.txt"; ref_file.write_text("\n".join(sorted(set(reference_pdbs) & set(pairs["pdb_id"]))) + "\n")
        command = [config["biojava"]["java_executable"], "-Xmx" + config["runtime"]["java_xmx"], "-jar", config["biojava"]["jar"],
                   "--pairs", str(work / "pairs.tsv"), "--atoms", str(work / "target_atoms.tsv"), "--mmcif-root", config["input"]["mmcif_root"],
                   "--out-dir", str(work / "java_output"), "--cell-mode", config["search"]["production_cell_mode"],
                   "--num-cells", str(config["search"]["exhaustive_num_cells"]), "--cutoff", str(config["search"]["cutoff_angstrom"]),
                   "--compare-auto-exhaustive", str(compare).lower(), "--reference-pdbs", str(ref_file)]
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (work / "java.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode: raise RuntimeError(f"Java bucket {bucket} failed exit={completed.returncode}: {completed.stdout[-4000:]}")
    result = {"bucket": bucket, "pair_count": len(pairs), "runnable_pair_count": runnable, "target_frame_review_count": len(skipped), "completed_at": utc()}
    atomic_json(marker, result)
    return result


def concatenate_tsv(inputs: list[Path], output: Path, header: list[str] | None = None) -> int:
    output.parent.mkdir(parents=True, exist_ok=True); count = 0
    with gzip.open(output, "wt", encoding="utf-8", newline="") as out:
        wrote_header = False
        for path in inputs:
            if not path.exists(): continue
            with path.open("r", encoding="utf-8", newline="") as inp:
                first = inp.readline()
                if not first: continue
                if not wrote_header:
                    out.write("\t".join(header) + "\n" if header else first); wrote_header = True
                for line in inp:
                    out.write(line); count += 1
        if not wrote_header and header: out.write("\t".join(header) + "\n")
    return count


def finalize(run: Path, config: dict, expected_pairs: int, pilot: bool, bucket_results: list[dict]) -> dict:
    works = sorted((run / "work").glob("bucket_*")); output = run / "output"; output.mkdir(parents=True, exist_ok=True)
    inv_inputs=[]
    for w in works:
        inv_inputs += [w / "pre_skipped_inventory.tsv", w / "java_output/pair_inventory.tsv"]
    pair_rows = concatenate_tsv(inv_inputs, output / "01_pair_step1_inventory.tsv.gz", INVENTORY_HEADER)
    hit_rows = concatenate_tsv([w / "java_output/hits.tsv" for w in works], output / "02_crystal_neighbor_instances.tsv.gz")
    source_rows = concatenate_tsv([w / "java_output/source_objects.tsv" for w in works], output / "03_source_object_inventory.tsv.gz")
    frame_rows = concatenate_tsv([w / "target_frame_validation.tsv" for w in works], output / "04_target_frame_validation.tsv.gz", FRAME_HEADER)
    lattice_rows = concatenate_tsv([w / "java_output/lattice_validation.tsv" for w in works], output / "05_lattice_reconstruction_validation.tsv.gz")
    inventory = pd.read_csv(output / "01_pair_step1_inventory.tsv.gz", sep="\t", compression="gzip")
    frame = pd.read_csv(output / "04_target_frame_validation.tsv.gz", sep="\t", compression="gzip")
    lattice = pd.read_csv(output / "05_lattice_reconstruction_validation.tsv.gz", sep="\t", compression="gzip") if lattice_rows else pd.DataFrame()
    checks = {
        "pair_accounting": len(inventory) == expected_pairs,
        "pair_id_unique": inventory["candidate_pair_id"].nunique() == len(inventory),
        "silent_drop_zero": len(inventory) == frame["candidate_pair_id"].nunique(),
        "target_frame_no_unexplained_failure": not inventory["step1_status"].eq("TARGET_FRAME_REVIEW").any(),
        "structure_parse_success_100pct": not inventory["step1_status"].eq("STRUCTURE_PARSE_ERROR").any(),
        "duplicate_instance_key_zero": True,
        "stored_transform_qc": lattice.empty or float(lattice["transform_max_deviation"].fillna(0).max()) <= float(config["policy"]["transform_tolerance"]),
        "auto_exhaustive_missed_zero": (not pilot) or (not lattice.empty and int(lattice["auto_exhaustive_missed_instances"].sum()) == 0),
        "biojava_reference_missed_zero": (not pilot) or (not lattice.empty and int(lattice["biojava_reference_missed_placements"].sum()) == 0),
        "lattice_validation_pass": lattice.empty or lattice["validation_status"].eq("PASS").all(),
    }
    # Exact duplicate audit is streamed to avoid materializing the full hits table.
    seen=set(); duplicate=0
    for chunk in pd.read_csv(output / "02_crystal_neighbor_instances.tsv.gz", sep="\t", compression="gzip", usecols=["candidate_pair_id","crystal_instance_key"], chunksize=250_000):
        for key in zip(chunk["candidate_pair_id"],chunk["crystal_instance_key"]):
            if key in seen: duplicate += 1
            else: seen.add(key)
    checks["duplicate_instance_key_zero"] = duplicate == 0
    status_counts = inventory["step1_status"].value_counts().to_dict()
    summary = pd.DataFrame([{ "run_id": run.name, "mode": "pilot" if pilot else "full", "pair_count": len(inventory),
        "pdb_count": int(inventory["pdb_id"].nunique()), "neighbor_instance_count": hit_rows,
        "ligand_shell_instance_count": int(inventory["n_ligand_6A_instances"].fillna(0).sum()),
        "pocket_shell_instance_count": int(inventory["n_pocket_6A_instances"].fillna(0).sum()),
        "target_frame_review_count": int(inventory["step1_status"].eq("TARGET_FRAME_REVIEW").sum()),
        "validation_pass": all(checks.values()) }])
    summary.to_csv(output / "06_step1_summary.tsv", sep="\t", index=False)
    validation = {"run_id": run.name, "validated_at": utc(), "validation_pass": all(checks.values()), "checks": checks,
                  "counts": {"pairs": pair_rows, "hits": hit_rows, "source_objects": source_rows, "target_frame_rows": frame_rows,
                             "lattice_rows": lattice_rows, "duplicate_pair_instance_keys": duplicate}, "step1_status_counts": status_counts}
    atomic_json(run / "validation.json", validation)
    report = f"""# Filter 4 Step 1 — Local Crystallographic Lattice Neighbor Search\n\nRun: `{run.name}`  \nMode: `{'pilot' if pilot else 'full'}`  \nBioJava: `7.2.5` (`5d047ab428b176a437131d31171c7f779caa239e`)  \nCell mode: `{config['search']['production_cell_mode']}`; pilot reference: exhaustive ±{config['search']['exhaustive_num_cells']}  \nCutoff: `{config['search']['cutoff_angstrom']} Å` discovery shell only\n\n- PDB: {inventory['pdb_id'].nunique():,}\n- Pairs: {len(inventory):,}\n- Unique neighbor instances: {hit_rows:,}\n- Ligand-shell instances: {int(inventory['n_ligand_6A_instances'].fillna(0).sum()):,}\n- Pocket-shell instances: {int(inventory['n_pocket_6A_instances'].fillna(0).sum()):,}\n- Target-frame review: {int(inventory['step1_status'].eq('TARGET_FRAME_REVIEW').sum()):,}\n- Structure parse errors: {int(inventory['step1_status'].eq('STRUCTURE_PARSE_ERROR').sum()):,}\n- Validation: {'PASS' if validation['validation_pass'] else 'FAIL'}\n\nThis run discovers crystallographic molecular instances only. It does not infer biological-assembly equivalence, packing contamination, physical interactions, severity, or pair eligibility.\n"""
    (output / "07_step1_report.md").write_text(report, encoding="utf-8")
    schemas={}
    for path in sorted(output.glob("0*")):
        if path.suffix==".gz":
            frame0=pd.read_csv(path,sep="\t",compression="gzip",nrows=100)
        elif path.suffix==".tsv": frame0=pd.read_csv(path,sep="\t",nrows=100)
        else: continue
        schemas[path.name]=[{"column_name":c,"inferred_dtype":str(frame0[c].dtype)} for c in frame0.columns]
    atomic_json(run / "output_schema.json", {"schema_version":"filter4_step1_schema_1.0.0","files":schemas})
    shutil.copy2(Path(__file__), run / "executed_runner.py")
    shutil.copy2(run / "config_snapshot.yaml", run / "executed_config.yaml")
    manifest=[]
    for path in sorted(list(output.iterdir())+[run/"validation.json",run/"output_schema.json",run/"executed_runner.py"]):
        if path.is_file(): manifest.append({"relative_path":str(path.relative_to(run)),"size_bytes":path.stat().st_size,"sha256":sha256(path)})
    pd.DataFrame(manifest).to_csv(run / "output_manifest.tsv",sep="\t",index=False)
    write_sha256s(run)
    return validation


def main() -> None:
    ap=argparse.ArgumentParser();ap.add_argument("--config",required=True);ap.add_argument("--run-dir",required=True);ap.add_argument("--mode",choices=["pilot","full"],required=True);ap.add_argument("--workers",type=int);args=ap.parse_args()
    config=yaml.safe_load(Path(args.config).read_text());run=Path(args.run_dir);run.mkdir(parents=True,exist_ok=True)
    shutil.copy2(args.config,run/"config_snapshot.yaml");atomic_json(run/"run_metadata.json",{"run_id":run.name,"mode":args.mode,"started_at":utc(),"status":"RUNNING"})
    pairs=retained_pairs(config);pilot=args.mode=="pilot"
    if pilot:
        selection=select_pilot(pairs,config,run);pairs=pairs[pairs["pdb_id"].isin(set(selection["pdb_id"]))].copy();reference_pdbs=selection.loc[selection["biojava_reference_selected"],"pdb_id"].tolist()
    else: reference_pdbs=[]
    tasks=[]
    for bucket,group in pairs.groupby("bucket_id",sort=True):tasks.append((int(bucket),group.to_dict("records"),config,str(run),pilot and bool(config["pilot"]["compare_auto_exhaustive"]),reference_pdbs))
    workers=args.workers or int(config["runtime"]["workers"]);results=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(prepare_bucket,t):t[0] for t in tasks}
        for future in as_completed(futures):
            result=future.result();results.append(result);atomic_json(run/"progress.json",{"completed_buckets":len(results),"total_buckets":len(tasks),"pairs_completed":sum(x["pair_count"] for x in results),"updated_at":utc()})
    validation=finalize(run,config,len(pairs),pilot,results)
    atomic_json(run/"run_metadata.json",{"run_id":run.name,"mode":args.mode,"completed_at":utc(),"status":"VALIDATED" if validation["validation_pass"] else "VALIDATION_FAILED"})
    # run_metadata is finalized after output generation, so refresh top-level
    # checksums only after its terminal state has been written.
    write_sha256s(run)
    if not validation["validation_pass"]:raise SystemExit(2)


if __name__=="__main__":main()
