#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
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
from posebusters import PoseBusters
from rdkit import Chem
from rdkit.Geometry import Point3D


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN = ROOT / "runs/20260812_full_01"
P2 = Path("/root/autodl-tmp/benchmark_1.0/processing_2_assembly_ready_structure_preparation/runs/20260810_full_01")
P2_OUT = P2 / "output"
CCD = P2 / "input/ccd_active_snapshot.sqlite"
QUALITY = RUN / "work/quality_batches"
OUT = RUN / "work/posebusters_batches"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean(value) -> str:
    if value is None:
        return ""
    value = str(value).strip()
    return "" if value.lower() in {"", ".", "?", "none", "false", "nan"} else value


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def bucket_paths(root: Path, bucket: int) -> list[str]:
    return [str(path) for path in sorted((root / f"bucket_id={bucket:03d}").glob("*.parquet"))]


def read_bucket(root: Path, bucket: int, columns=None) -> pd.DataFrame:
    paths = bucket_paths(root, bucket)
    if not paths:
        return pd.DataFrame()
    return ds.dataset(paths, format="parquet").to_table(columns=columns).to_pandas(split_blocks=True, self_destruct=True)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), temporary, compression="zstd")
    os.replace(temporary, path)


def buster() -> PoseBusters:
    config = PoseBusters(config="mol", max_workers=0).config
    config["modules"] = [
        module for module in config["modules"]
        if module.get("name") != "Energy ratio"
    ]
    return PoseBusters(config=config, max_workers=0, chunk_size=100)


def ccd_atom_properties(component_ids: set[str]) -> dict[tuple[str, str], tuple[int, bool]]:
    result = {}
    connection = sqlite3.connect(f"file:{CCD}?mode=ro", uri=True)
    values = sorted(component_ids)
    for start in range(0, len(values), 500):
        chunk = values[start:start + 500]
        placeholders = ",".join("?" for _ in chunk)
        query = f"SELECT component_id, atom_id, charge, aromatic_flag FROM atoms WHERE component_id IN ({placeholders})"
        for component, atom_id, charge, aromatic in connection.execute(query, chunk):
            result[(component, atom_id)] = (int(charge or 0), str(aromatic).upper() == "Y")
    connection.close()
    return result


BOND_TYPES = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "AROM": Chem.BondType.AROMATIC,
    "DELO": Chem.BondType.AROMATIC,
}


def build_molecule(atoms: pd.DataFrame, bonds: pd.DataFrame, properties: dict) -> tuple[Chem.Mol | None, str]:
    try:
        atoms = atoms.drop_duplicates("label_atom_id").copy()
        rw = Chem.RWMol()
        atom_index = {}
        component = clean(atoms.iloc[0]["component_id"]).upper()
        conformer = Chem.Conformer(len(atoms))
        for position, row in enumerate(atoms.to_dict("records")):
            atom_id = clean(row["label_atom_id"])
            symbol = clean(row["type_symbol"])
            symbol = symbol[:1].upper() + symbol[1:].lower()
            atom = Chem.Atom(symbol)
            charge, aromatic = properties.get((component, atom_id), (0, False))
            atom.SetFormalCharge(charge)
            atom.SetIsAromatic(aromatic)
            index = rw.AddAtom(atom)
            atom_index[atom_id] = index
            conformer.SetAtomPosition(index, Point3D(float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])))
        for row in bonds.to_dict("records"):
            left = clean(row["atom_id_1"])
            right = clean(row["atom_id_2"])
            if left not in atom_index or right not in atom_index:
                continue
            aromatic = clean(row.get("aromatic_flag")).upper() == "Y"
            bond_type = Chem.BondType.AROMATIC if aromatic else BOND_TYPES.get(clean(row.get("bond_order")).upper(), Chem.BondType.SINGLE)
            if rw.GetBondBetweenAtoms(atom_index[left], atom_index[right]) is None:
                rw.AddBond(atom_index[left], atom_index[right], bond_type)
                bond = rw.GetBondBetweenAtoms(atom_index[left], atom_index[right])
                if bond_type == Chem.BondType.AROMATIC:
                    bond.SetIsAromatic(True)
                    rw.GetAtomWithIdx(atom_index[left]).SetIsAromatic(True)
                    rw.GetAtomWithIdx(atom_index[right]).SetIsAromatic(True)
        molecule = rw.GetMol()
        molecule.AddConformer(conformer)
        molecule.SetProp("_Name", clean(atoms.iloc[0]["filter_2_source_ligand_instance_id"]))
        return molecule, "BUILD_SUCCESS"
    except Exception as exc:
        return None, f"BUILD_FAILED:{type(exc).__name__}:{exc}"[:1000]


def normalize_report(report: pd.DataFrame, source_id: str) -> dict:
    if report.empty:
        return {"source_ligand_instance_id": source_id, "posebusters_status": "NO_REPORT"}
    values = report.iloc[0].to_dict()
    normalized = {"source_ligand_instance_id": source_id, "posebusters_status": "COMPLETED"}
    for key, value in values.items():
        if isinstance(key, tuple):
            key = "__".join(clean(part) for part in key if clean(part))
        key = clean(key).lower().replace(" ", "_").replace("≤", "le")
        if isinstance(value, (bool, int, float, str)) or value is None:
            normalized[key] = value
        else:
            normalized[key] = str(value)
    return normalized


def process_bucket(bucket: int, limit: int | None = None) -> dict:
    started = time.time()
    bucket_out = OUT / f"bucket_id={bucket:03d}"
    marker = bucket_out / "_COMPLETE.json"
    if marker.exists() and limit is None:
        return json.loads(marker.read_text())
    pair_path = QUALITY / f"bucket_id={bucket:03d}" / "pair_quality_pre_posebusters.parquet"
    if not pair_path.exists():
        raise RuntimeError(f"quality bucket not complete: {bucket}")
    pairs = pq.read_table(pair_path).to_pandas()
    eligible = pairs[pairs["terminal_status_pre_posebusters"].isin(["FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY"])]
    placement_ids = set(eligible["ligand_assembly_placement_id"])
    atoms = read_bucket(P2_OUT / "prepared_ligand_assembly_atoms", bucket, [
        "filter_2_ligand_assembly_placement_id", "filter_2_source_ligand_instance_id", "component_id",
        "label_atom_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z",
    ])
    atoms = atoms[atoms["filter_2_ligand_assembly_placement_id"].isin(placement_ids)].copy()
    bonds = read_bucket(P2_OUT / "prepared_ligand_assembly_bonds", bucket)
    source_ids = list(dict.fromkeys(atoms["filter_2_source_ligand_instance_id"]))
    if limit is not None:
        source_ids = source_ids[:limit]
    source_set = set(source_ids)
    atoms = atoms[atoms["filter_2_source_ligand_instance_id"].isin(source_set)]
    bonds = bonds[bonds["filter_2_source_ligand_instance_id"].isin(source_set)]
    properties = ccd_atom_properties(set(atoms["component_id"].str.upper()))
    atom_groups = {key: group for key, group in atoms.groupby("filter_2_source_ligand_instance_id", sort=False)}
    bond_groups = {key: group for key, group in bonds.groupby("filter_2_source_ligand_instance_id", sort=False)}
    tool = buster()
    rows = []
    for source_id in source_ids:
        molecule, build_status = build_molecule(atom_groups[source_id], bond_groups.get(source_id, pd.DataFrame()), properties)
        if molecule is None:
            rows.append({"source_ligand_instance_id": source_id, "posebusters_status": build_status})
            continue
        try:
            report = tool.bust(molecule, full_report=True)
            row = normalize_report(report, source_id)
            row["molecule_build_status"] = build_status
            rows.append(row)
        except Exception as exc:
            rows.append({"source_ligand_instance_id": source_id, "posebusters_status": f"EXECUTION_FAILED:{type(exc).__name__}:{exc}"[:1000], "molecule_build_status": build_status})
    frame = pd.DataFrame(rows)
    output = bucket_out / ("posebusters_test.parquet" if limit is not None else "posebusters_results.parquet")
    write_parquet(frame, output)
    result = {
        "status": "COMPLETED",
        "bucket_id": bucket,
        "eligible_pair_count": len(eligible),
        "unique_source_ligand_count": len(source_ids),
        "posebusters_status": dict(Counter(frame["posebusters_status"])),
        "runtime_seconds": time.time() - started,
        "finished_at": utc(),
    }
    if limit is None:
        atomic_json(marker, result)
    else:
        atomic_json(bucket_out / "_TEST.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bucket", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=255)
    args = parser.parse_args()
    if args.limit is not None and args.bucket is None:
        raise SystemExit("--limit requires --bucket")
    buckets = [args.bucket] if args.bucket is not None else list(range(args.start, args.end + 1))
    started = time.time()
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_bucket, bucket, args.limit): bucket for bucket in buckets}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            results.append(result)
            progress = {
                "status": "RUNNING",
                "phase": "POSEBUSTERS_RAW_LIGAND_GEOMETRY",
                "bucket_completed": len(results),
                "bucket_total": len(buckets),
                "eligible_pair_count": sum(row["eligible_pair_count"] for row in results),
                "unique_source_ligand_count": sum(row["unique_source_ligand_count"] for row in results),
                "runtime_seconds": time.time() - started,
                "updated_at": utc(),
            }
            atomic_json(RUN / "posebusters_status.json", progress)
            print(json.dumps(progress), flush=True)
    progress["status"] = "COMPLETED"
    progress["phase"] = "POSEBUSTERS_RAW_LIGAND_GEOMETRY_COMPLETE"
    atomic_json(RUN / "posebusters_status.json", progress)


if __name__ == "__main__":
    main()
