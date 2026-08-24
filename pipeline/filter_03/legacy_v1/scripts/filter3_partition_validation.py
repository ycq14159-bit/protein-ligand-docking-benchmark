#!/usr/bin/env python3
import json
import os
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN = ROOT / "runs/20260812_full_01"
P3 = Path(
    "/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/"
    "runs/20260811_full_01/output/provisional_pairs"
)
SOURCE = RUN / "work/normalized_validation/batches"
OUT = RUN / "work/normalized_validation/by_bucket"


def utc():
    return datetime.now(timezone.utc).isoformat()


def update_status(value):
    path = RUN / "status.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(tmp, path)


def main():
    pair_table = ds.dataset(P3, format="parquet", partitioning="hive").to_table(
        columns=["pdb_id", "bucket_id"]
    )
    pdb_to_bucket = {}
    conflicts = set()
    for pdb_id, bucket_id in zip(
        pair_table["pdb_id"].to_pylist(), pair_table["bucket_id"].to_pylist()
    ):
        previous = pdb_to_bucket.setdefault(pdb_id, int(bucket_id))
        if previous != int(bucket_id):
            conflicts.add(pdb_id)
    if conflicts:
        raise RuntimeError(f"PDB assigned to multiple buckets: {len(conflicts)}")

    OUT.mkdir(parents=True, exist_ok=True)
    entry_writers = {}
    residue_writers = {}
    totals = Counter()
    batches = sorted(SOURCE.glob("batch_*"))
    started = time.time()
    try:
        for index, batch in enumerate(batches, 1):
            for dataset_name, writers in (
                ("entry_validation", entry_writers),
                ("residue_validation", residue_writers),
            ):
                table = pq.read_table(batch / f"{dataset_name}.parquet")
                bucket_values = [pdb_to_bucket.get(pdb_id, -1) for pdb_id in table["pdb_id"].to_pylist()]
                if any(value < 0 for value in bucket_values):
                    totals[f"{dataset_name}_outside_input_rows"] += sum(value < 0 for value in bucket_values)
                table = table.append_column("bucket_id", pa.array(bucket_values, type=pa.int16()))
                for bucket_id in set(bucket_values):
                    if bucket_id < 0:
                        continue
                    subset = table.filter(pc.equal(table["bucket_id"], bucket_id))
                    writer = writers.get(bucket_id)
                    if writer is None:
                        directory = OUT / dataset_name / f"bucket_id={bucket_id:03d}"
                        directory.mkdir(parents=True, exist_ok=True)
                        writer = pq.ParquetWriter(
                            directory / "part-000000.parquet",
                            subset.schema,
                            compression="zstd",
                        )
                        writers[bucket_id] = writer
                    writer.write_table(subset)
                    totals[f"{dataset_name}_rows"] += subset.num_rows
            if index % 20 == 0 or index == len(batches):
                elapsed = time.time() - started
                update_status(
                    {
                        "status": "RUNNING",
                        "phase": "VALIDATION_REPARTITION",
                        "batch_done": index,
                        "batch_total": len(batches),
                        **dict(totals),
                        "runtime_seconds": elapsed,
                        "updated_at": utc(),
                    }
                )
    finally:
        for writer in list(entry_writers.values()) + list(residue_writers.values()):
            writer.close()
    status = {
        "status": "COMPLETED",
        "phase": "VALIDATION_REPARTITION_COMPLETE",
        "batch_done": len(batches),
        "batch_total": len(batches),
        "bucket_count": len(entry_writers),
        **dict(totals),
        "runtime_seconds": time.time() - started,
        "finished_at": utc(),
    }
    update_status(status)
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
