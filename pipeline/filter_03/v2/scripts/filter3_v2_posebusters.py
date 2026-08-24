#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import sqlite3
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import yaml
from posebusters import PoseBusters
from rdkit import Chem
from rdkit.Geometry import Point3D


def utc(): return datetime.now(timezone.utc).isoformat()


def clean(value):
    if value is None: return ""
    text = str(value).strip()
    return "" if text.lower() in {"", ".", "?", "none", "false", "nan"} else text


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def write_parquet(frame, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), tmp, compression="zstd")
    os.replace(tmp, path)


def read_bucket(root, bucket, columns=None):
    paths = [str(path) for path in sorted((Path(root) / f"bucket_id={bucket:03d}").glob("*.parquet"))]
    if not paths: return pd.DataFrame()
    return ds.dataset(paths, format="parquet").to_table(columns=columns).to_pandas(split_blocks=True, self_destruct=True)


def buster():
    config = PoseBusters(config="mol", max_workers=0).config
    config["modules"] = [module for module in config["modules"] if module.get("name") != "Energy ratio"]
    return PoseBusters(config=config, max_workers=0, chunk_size=100)


def ccd_atom_properties(ccd, component_ids):
    result = {}
    connection = sqlite3.connect(f"file:{ccd}?mode=ro", uri=True)
    values = sorted(component_ids)
    for start in range(0, len(values), 500):
        chunk = values[start:start + 500]
        placeholders = ",".join("?" for _ in chunk)
        query = f"SELECT component_id, atom_id, charge, aromatic_flag FROM atoms WHERE component_id IN ({placeholders})"
        for component, atom_id, charge, aromatic in connection.execute(query, chunk):
            result[(component, atom_id)] = (int(charge or 0), str(aromatic).upper() == "Y")
    connection.close()
    return result


BOND_TYPES = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC,
              "DELO": Chem.BondType.AROMATIC}


def build_molecule(atoms, bonds, properties):
    try:
        atoms = atoms.drop_duplicates("label_atom_id").copy()
        rw, atom_index = Chem.RWMol(), {}
        component = clean(atoms.iloc[0]["component_id"]).upper()
        conformer = Chem.Conformer(len(atoms))
        for position, row in enumerate(atoms.to_dict("records")):
            atom_id, symbol = clean(row["label_atom_id"]), clean(row["type_symbol"])
            atom = Chem.Atom(symbol[:1].upper() + symbol[1:].lower())
            charge, aromatic = properties.get((component, atom_id), (0, False))
            atom.SetFormalCharge(charge); atom.SetIsAromatic(aromatic)
            index = rw.AddAtom(atom); atom_index[atom_id] = index
            conformer.SetAtomPosition(index, Point3D(float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])))
        for row in bonds.to_dict("records"):
            left, right = clean(row["atom_id_1"]), clean(row["atom_id_2"])
            if left not in atom_index or right not in atom_index: continue
            aromatic = clean(row.get("aromatic_flag")).upper() == "Y"
            bond_type = Chem.BondType.AROMATIC if aromatic else BOND_TYPES.get(clean(row.get("bond_order")).upper(), Chem.BondType.SINGLE)
            if rw.GetBondBetweenAtoms(atom_index[left], atom_index[right]) is None:
                rw.AddBond(atom_index[left], atom_index[right], bond_type)
                if bond_type == Chem.BondType.AROMATIC:
                    bond = rw.GetBondBetweenAtoms(atom_index[left], atom_index[right]); bond.SetIsAromatic(True)
                    rw.GetAtomWithIdx(atom_index[left]).SetIsAromatic(True); rw.GetAtomWithIdx(atom_index[right]).SetIsAromatic(True)
        molecule = rw.GetMol(); molecule.AddConformer(conformer)
        molecule.SetProp("_Name", clean(atoms.iloc[0]["filter_2_source_ligand_instance_id"]))
        return molecule, "BUILD_SUCCESS"
    except Exception as exc:
        return None, f"BUILD_FAILED:{type(exc).__name__}:{exc}"[:1000]


def normalize_report(report, source_id):
    if report.empty: return {"source_ligand_instance_id": source_id, "posebusters_status": "NO_REPORT"}
    normalized = {"source_ligand_instance_id": source_id, "posebusters_status": "COMPLETED"}
    for key, value in report.iloc[0].to_dict().items():
        if isinstance(key, tuple): key = "__".join(clean(part) for part in key if clean(part))
        key = clean(key).lower().replace(" ", "_")
        normalized[key] = value if isinstance(value, (bool, int, float, str)) or value is None else str(value)
    return normalized


def process_bucket(bucket, config, run_dir):
    started = time.time(); run_dir = Path(run_dir)
    out = run_dir / "work/posebusters_new_batches" / f"bucket_id={bucket:03d}"
    marker = out / "_COMPLETE.json"
    if marker.exists(): return json.loads(marker.read_text())
    pending_path = run_dir / "work/preclassification_batches" / f"bucket_id={bucket:03d}" / "posebusters_pending_sources.parquet"
    if not pending_path.exists():
        result = {"status": "COMPLETED", "bucket_id": bucket, "source_count": 0, "runtime_seconds": 0, "finished_at": utc()}
        atomic_json(marker, result); return result
    pending = pq.ParquetFile(pending_path).read().to_pandas()
    source_ids = set(pending["source_ligand_instance_id"])
    p2 = Path(config["input"]["processing2_output"])
    atoms = read_bucket(p2 / "prepared_ligand_assembly_atoms", bucket, [
        "filter_2_source_ligand_instance_id", "component_id", "label_atom_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"])
    bonds = read_bucket(p2 / "prepared_ligand_assembly_bonds", bucket)
    atoms = atoms[atoms["filter_2_source_ligand_instance_id"].isin(source_ids)].copy()
    bonds = bonds[bonds["filter_2_source_ligand_instance_id"].isin(source_ids)].copy()
    atom_groups = {key: frame for key, frame in atoms.groupby("filter_2_source_ligand_instance_id", sort=False)}
    bond_groups = {key: frame for key, frame in bonds.groupby("filter_2_source_ligand_instance_id", sort=False)}
    ccd = Path('/root/autodl-tmp/benchmark_1.0/processing_2_assembly_ready_structure_preparation/runs/20260810_full_01/input/ccd_active_snapshot.sqlite')
    properties = ccd_atom_properties(ccd, set(atoms["component_id"].astype(str).str.upper()))
    tool, rows = buster(), []
    for source_id in sorted(source_ids):
        if source_id not in atom_groups:
            rows.append({"source_ligand_instance_id": source_id, "posebusters_status": "SOURCE_ATOMS_MISSING"}); continue
        molecule, status = build_molecule(atom_groups[source_id], bond_groups.get(source_id, pd.DataFrame()), properties)
        if molecule is None:
            rows.append({"source_ligand_instance_id": source_id, "posebusters_status": status, "molecule_build_status": status}); continue
        try:
            row = normalize_report(tool.bust(molecule, full_report=True), source_id); row["molecule_build_status"] = status; rows.append(row)
        except Exception as exc:
            rows.append({"source_ligand_instance_id": source_id, "posebusters_status": f"EXECUTION_FAILED:{type(exc).__name__}:{exc}"[:1000], "molecule_build_status": status})
    frame = pd.DataFrame(rows)
    if not frame.empty: write_parquet(frame, out / "posebusters_results.parquet")
    result = {"status": "COMPLETED", "bucket_id": bucket, "source_count": len(source_ids),
              "posebusters_status_counts": dict(Counter(frame.get("posebusters_status", []))),
              "runtime_seconds": time.time() - started, "finished_at": utc()}
    atomic_json(marker, result); return result


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); parser.add_argument("--run-dir", required=True); parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(); config = yaml.safe_load(Path(args.config).read_text()); started = time.time(); results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_bucket, bucket, config, args.run_dir): bucket for bucket in range(256)}
        for future in concurrent.futures.as_completed(futures):
            result = future.result(); results.append(result)
            progress = {"status": "RUNNING", "phase": "FILTER3_V2_POSEBUSTERS_MISSING_ONLY", "bucket_completed": len(results), "bucket_total": 256,
                        "source_count": sum(item["source_count"] for item in results), "runtime_seconds": time.time() - started, "updated_at": utc()}
            atomic_json(Path(args.run_dir) / "posebusters_status.json", progress); print(json.dumps(progress), flush=True)
    progress["status"] = "COMPLETED"; progress["phase"] = "FILTER3_V2_POSEBUSTERS_COMPLETE"; atomic_json(Path(args.run_dir) / "posebusters_status.json", progress)


if __name__ == "__main__": main()
