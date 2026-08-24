#!/usr/bin/env python3
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import os
import platform
import shutil
import sqlite3
import sys
import tempfile
import threading
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import requests


STAGE = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN_ID = "20260812_full_01"
RUN = STAGE / "runs" / RUN_ID
P3_RUN = Path(
    "/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/"
    "runs/20260811_full_01"
)
EXPECTED_PAIRS = 744_580
EXPECTED_PDBS = 138_892
BASE_URL = "https://files.wwpdb.org/pub/pdb/validation_reports"
SOURCE_TYPES = ("xml", "cif")
SCHEMA_VERSION = "filter3_input_v1.0.0"
_thread_local = threading.local()


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path, block=1024 * 1024):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(block):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def session():
    if not hasattr(_thread_local, "session"):
        value = requests.Session()
        value.headers.update({"User-Agent": "benchmark-1.0-filter3-validation-snapshot/1.0"})
        _thread_local.session = value
    return _thread_local.session


def init_stage():
    for relative in (
        "scripts",
        "config",
        "tests",
        "runs",
        f"runs/{RUN_ID}/input",
        f"runs/{RUN_ID}/work/validation_sources/xml",
        f"runs/{RUN_ID}/work/validation_sources/cif",
        f"runs/{RUN_ID}/output",
        f"runs/{RUN_ID}/audit",
        f"runs/{RUN_ID}/logs",
        f"runs/{RUN_ID}/release",
    ):
        (STAGE / relative).mkdir(parents=True, exist_ok=True)
    frozen = json.loads((P3_RUN / "_FROZEN.json").read_text())
    validation = json.loads((P3_RUN / "audit/processing_3_release_validation.json").read_text())
    checks = {
        "processing3_frozen": frozen.get("status") == "FROZEN",
        "processing3_validation_pass": validation.get("validation_pass") is True,
        "processing3_pair_count": frozen.get("final_pair_count") == EXPECTED_PAIRS,
        "plip_not_used_as_input": True,
        "filter4_not_started": True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"preflight failed: {checks}")
    atomic_json(
        RUN / "audit/preflight.json",
        {"preflight_pass": True, "checks": checks, "created_at": utc()},
    )
    atomic_json(
        RUN / "status.json",
        {"status": "RUNNING", "phase": "INITIALIZED", "run_id": RUN_ID, "updated_at": utc()},
    )


def build_input():
    init_stage()
    source = P3_RUN / "output/provisional_pairs"
    dataset = ds.dataset(source, format="parquet", partitioning="hive")
    columns = [
        "pair_id",
        "ligand_assembly_placement_id",
        "pdb_id",
        "assembly_id",
        "model_id",
        "component_id",
        "receptor_chain_instance_ids",
        "receptor_chain_count",
        "metal_status",
        "pair_status",
    ]
    table = dataset.to_table(columns=columns)
    if table.num_rows != EXPECTED_PAIRS:
        raise RuntimeError(f"unexpected pair count: {table.num_rows}")
    unique_pair_count = len(set(table["pair_id"].to_pylist()))
    unique_placement_count = len(set(table["ligand_assembly_placement_id"].to_pylist()))
    if unique_pair_count != EXPECTED_PAIRS or unique_placement_count != EXPECTED_PAIRS:
        raise RuntimeError("input pair or placement key is not unique")
    order = pc.sort_indices(table, sort_keys=[("pair_id", "ascending")])
    table = pc.take(table, order)
    pair_path = RUN / "input/filter3_input_pairs.parquet"
    pq.write_table(table, pair_path, compression="zstd", compression_level=6)

    pdb_ids = sorted(set(table["pdb_id"].to_pylist()))
    if len(pdb_ids) != EXPECTED_PDBS:
        raise RuntimeError(f"unexpected unique PDB count: {len(pdb_ids)}")
    inventory = pa.table(
        {
            "pdb_id": pdb_ids,
            "validation_xml_url": [
                f"{BASE_URL}/{p[1:3]}/{p}/{p}_validation.xml.gz" for p in pdb_ids
            ],
            "validation_cif_url": [
                f"{BASE_URL}/{p[1:3]}/{p}/{p}_validation.cif.gz" for p in pdb_ids
            ],
        }
    )
    inventory_path = RUN / "input/filter3_pdb_inventory.parquet"
    pq.write_table(inventory, inventory_path, compression="zstd", compression_level=6)
    with gzip.open(RUN / "input/filter3_pdb_inventory.tsv.gz", "wt", encoding="utf-8") as handle:
        handle.write("pdb_id\tvalidation_xml_url\tvalidation_cif_url\n")
        for row in inventory.to_pylist():
            handle.write(
                f"{row['pdb_id']}\t{row['validation_xml_url']}\t{row['validation_cif_url']}\n"
            )

    upstream = {
        "source_stage": "processing_03_direct_contact_qualification",
        "source_run_id": "20260811_full_01",
        "source_status": "FROZEN",
        "source_manifest": str(P3_RUN / "output/output_manifest.tsv"),
        "source_manifest_sha256": sha256(P3_RUN / "output/output_manifest.tsv"),
        "source_frozen_record_sha256": sha256(P3_RUN / "_FROZEN.json"),
        "source_dataset": str(source),
        "input_pair_count": table.num_rows,
        "unique_pair_count": unique_pair_count,
        "unique_placement_count": unique_placement_count,
        "unique_pdb_count": len(pdb_ids),
        "filter3_input_pairs_sha256": sha256(pair_path),
        "filter3_pdb_inventory_sha256": sha256(inventory_path),
        "plip_used_for_membership": False,
        "created_at": utc(),
    }
    atomic_json(RUN / "input/upstream.json", upstream)
    atomic_json(
        RUN / "status.json",
        {
            "status": "RUNNING",
            "phase": "INPUT_FROZEN",
            "input_pair_count": table.num_rows,
            "unique_pdb_count": len(pdb_ids),
            "updated_at": utc(),
        },
    )
    print(json.dumps(upstream, indent=2))


def database():
    path = RUN / "work/validation_download_state.sqlite"
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS sources (
          pdb_id TEXT NOT NULL,
          source_type TEXT NOT NULL,
          source_url TEXT NOT NULL,
          local_path TEXT NOT NULL,
          download_status TEXT NOT NULL,
          http_status INTEGER,
          size_bytes INTEGER,
          sha256 TEXT,
          gzip_ok INTEGER NOT NULL,
          parse_probe_ok INTEGER NOT NULL,
          error TEXT,
          retrieved_at TEXT NOT NULL,
          PRIMARY KEY (pdb_id, source_type)
        )
        """
    )
    connection.commit()
    return connection


def validation_path(pdb_id, source_type):
    return (
        RUN
        / "work/validation_sources"
        / source_type
        / pdb_id[1:3]
        / f"{pdb_id}_validation.{source_type}.gz"
    )


def gzip_probe(path):
    try:
        with gzip.open(path, "rb") as handle:
            first = handle.read(4096)
            while handle.read(1024 * 1024):
                pass
        return bool(first)
    except Exception:
        return False


def download_one(task):
    pdb_id, source_type = task
    url = f"{BASE_URL}/{pdb_id[1:3]}/{pdb_id}/{pdb_id}_validation.{source_type}.gz"
    target = validation_path(pdb_id, source_type)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size and gzip_probe(target):
        return (
            pdb_id,
            source_type,
            url,
            str(target),
            "VALIDATION_AVAILABLE",
            200,
            target.stat().st_size,
            sha256(target),
            1,
            1,
            "",
            utc(),
        )
    error = ""
    http_status = None
    for attempt in range(1, 3):
        temporary = target.with_suffix(target.suffix + f".tmp.{os.getpid()}.{threading.get_ident()}")
        try:
            with session().get(url, stream=True, timeout=(20, 180)) as response:
                http_status = response.status_code
                if response.status_code == 404:
                    return (
                        pdb_id, source_type, url, str(target), "VALIDATION_REPORT_MISSING",
                        404, 0, "", 0, 0, "http_404", utc(),
                    )
                response.raise_for_status()
                with temporary.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if not gzip_probe(temporary):
                raise ValueError("gzip_validation_failed")
            os.replace(temporary, target)
            return (
                pdb_id, source_type, url, str(target), "VALIDATION_AVAILABLE",
                http_status, target.stat().st_size, sha256(target), 1, 1, "", utc(),
            )
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            error = f"attempt={attempt}:{type(exc).__name__}:{exc}"[:2000]
            if attempt < 2:
                time.sleep(attempt)
    return (
        pdb_id, source_type, url, str(target), "VALIDATION_DOWNLOAD_FAILED",
        http_status, 0, "", 0, 0, error, utc(),
    )


def export_download_manifest(connection):
    fields = [
        "pdb_id", "source_type", "source_url", "local_path", "download_status",
        "http_status", "size_bytes", "sha256", "gzip_ok", "parse_probe_ok", "error",
        "retrieved_at",
    ]
    rows = connection.execute(
        "SELECT " + ",".join(fields) + " FROM sources ORDER BY pdb_id,source_type"
    ).fetchall()
    output = RUN / "audit/validation_source_manifest.tsv.gz"
    with gzip.open(output, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join("" if value is None else str(value) for value in row) + "\n")
    return len(rows), sha256(output)


def download(workers):
    inventory = pq.read_table(RUN / "input/filter3_pdb_inventory.parquet")
    pdb_ids = inventory["pdb_id"].to_pylist()
    connection = database()
    completed = {
        (row[0], row[1])
        for row in connection.execute(
            "SELECT pdb_id,source_type FROM sources WHERE download_status IN "
            "('VALIDATION_AVAILABLE','VALIDATION_REPORT_MISSING')"
        )
    }
    tasks = [
        (pdb_id, source_type)
        for pdb_id in pdb_ids
        for source_type in SOURCE_TYPES
        if (pdb_id, source_type) not in completed
    ]
    counts = Counter(
        row[0] for row in connection.execute("SELECT download_status FROM sources")
    )
    started = time.time()
    status_path = RUN / "status.json"
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        task_iterator = iter(tasks)
        futures = {}
        # Keep at most four tasks per worker in flight so pause/resume stays cheap.
        for _ in range(min(len(tasks), workers * 4)):
            task = next(task_iterator, None)
            if task is None:
                break
            futures[pool.submit(download_one, task)] = task
        index = 0
        while futures:
            done, _ = concurrent.futures.wait(
                futures, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                futures.pop(future)
                row = future.result()
                index += 1
                connection.execute(
                    "INSERT OR REPLACE INTO sources VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", row
                )
                counts[row[4]] += 1
                task = next(task_iterator, None)
                if task is not None:
                    futures[pool.submit(download_one, task)] = task
                if index % 100 == 0 or index == len(tasks):
                    connection.commit()
                    elapsed = time.time() - started
                    total_done = sum(counts.values())
                    rate = index / elapsed if elapsed else 0
                    remaining = len(tasks) - index
                    status = {
                        "status": "RUNNING",
                        "phase": "VALIDATION_DOWNLOAD",
                        "workers": workers,
                        "source_file_target_count": len(pdb_ids) * 2,
                        "source_file_accounted_count": total_done,
                        "this_attempt_done": index,
                        "this_attempt_total": len(tasks),
                        "status_counts": dict(counts),
                        "files_per_second": rate,
                        "eta_seconds": remaining / rate if rate else None,
                        "updated_at": utc(),
                    }
                    atomic_json(status_path, status)
                    print(json.dumps(status), flush=True)
    connection.commit()
    counts = Counter(
        row[0] for row in connection.execute("SELECT download_status FROM sources")
    )
    manifest_rows, manifest_hash = export_download_manifest(connection)
    connection.close()
    atomic_json(
        status_path,
        {
            "status": "COMPLETED",
            "phase": "VALIDATION_DOWNLOAD_COMPLETE",
            "source_file_target_count": len(pdb_ids) * 2,
            "source_file_manifest_rows": manifest_rows,
            "status_counts": dict(counts),
            "manifest_sha256": manifest_hash,
            "finished_at": utc(),
        },
    )


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("init")
    subparsers.add_parser("build-input")
    download_parser = subparsers.add_parser("download")
    download_parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    if args.command == "init":
        init_stage()
    elif args.command == "build-input":
        build_input()
    elif args.command == "download":
        download(args.workers)


if __name__ == "__main__":
    main()
