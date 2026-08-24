#!/usr/bin/env python3
import hashlib
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

ROOT = Path("/root/autodl-tmp/benchmark_1.0/auxiliary_entry_work_packages")
DATASETS = [
    "entry_ligand_placements", "entry_receptor_chain_instances", "entry_assembly_context",
    "entry_ligand_source_atoms", "entry_receptor_source_atoms",
]


def _bucket(pdb_id: str) -> int:
    return int(hashlib.sha256(pdb_id.lower().encode()).hexdigest()[:8], 16) % 256


def _validate_pdb_id(pdb_id: str) -> str:
    value = pdb_id.strip().lower()
    if not re.fullmatch(r"[0-9a-z]{4,12}", value):
        raise ValueError(f"invalid PDB ID: {pdb_id!r}")
    return value


def load_entry_work_package(
    pdb_id: str,
    datasets: list[str] | None = None,
    columns: dict[str, list[str]] | None = None,
):
    pid = _validate_pdb_id(pdb_id)
    current = __import__("json").loads((ROOT / "CURRENT_BUILD.json").read_text())
    output = ROOT / current["relative_path"] / "output"
    manifest = pq.read_table(output / "entry_work_manifest.parquet", filters=[("pdb_id", "=", pid)])
    if manifest.num_rows != 1:
        raise KeyError(f"PDB ID not found exactly once in manifest: {pid}")
    row = manifest.to_pylist()[0]
    bid = _bucket(pid)
    if row["bucket_id"] != bid:
        raise RuntimeError(f"manifest bucket mismatch for {pid}: {row['bucket_id']} != {bid}")
    selected = DATASETS if datasets is None else datasets
    unknown = sorted(set(selected) - set(DATASETS))
    if unknown:
        raise ValueError(f"unknown datasets: {unknown}")
    result = {"manifest": manifest}
    for name in selected:
        bucket_path = output / name / f"bucket_id={bid:03d}"
        requested = None if columns is None else columns.get(name)
        if not bucket_path.exists():
            table = pa.table({"pdb_id": pa.array([], type=pa.string())})
        else:
            table = ds.dataset(bucket_path, format="parquet").to_table(columns=requested, filter=pc.field("pdb_id") == pid)
        expected = {
            "entry_ligand_placements": row["ligand_placement_count"],
            "entry_receptor_chain_instances": row["active_receptor_chain_instance_count"],
            "entry_ligand_source_atoms": row["ligand_source_atom_count"],
            "entry_receptor_source_atoms": row["receptor_source_atom_count"],
        }.get(name)
        if expected is not None and table.num_rows != expected:
            raise RuntimeError(f"{name} row mismatch for {pid}: {table.num_rows} != {expected}")
        result[name] = table
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_id")
    parser.add_argument("--datasets", nargs="*")
    args = parser.parse_args()
    package = load_entry_work_package(args.pdb_id, datasets=args.datasets)
    print({name: table.num_rows for name, table in package.items()})
