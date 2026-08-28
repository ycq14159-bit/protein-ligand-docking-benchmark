#!/usr/bin/env python3
"""Populate a resumable wwPDB validation XML cache for a fixed PDB-ID list."""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://files.rcsb.org/pub/pdb/validation_reports"
THREAD_LOCAL = threading.local()


def destination(root: Path, pdb_id: str) -> Path:
    return root / pdb_id[1:3] / pdb_id / f"{pdb_id}_validation.xml.gz"


def valid_gzip(path: Path) -> bool:
    try:
        if path.stat().st_size < 100:
            return False
        with gzip.open(path, "rb") as handle:
            head = handle.read(512)
        return b"<?xml" in head or b"ValidationReport" in head
    except Exception:
        return False


def session() -> requests.Session:
    if not hasattr(THREAD_LOCAL, "session"):
        value = requests.Session()
        value.headers["User-Agent"] = "benchmark-harmonized-validation-cache/1.0"
        THREAD_LOCAL.session = value
    return THREAD_LOCAL.session


def download_one(pdb_id: str, output_root: Path, retries: int) -> dict:
    target = destination(output_root, pdb_id)
    if valid_gzip(target):
        return {"pdb_id": pdb_id, "status": "ALREADY_CACHED", "bytes": target.stat().st_size}
    target.parent.mkdir(parents=True, exist_ok=True)
    url = f"{BASE_URL}/{pdb_id[1:3]}/{pdb_id}/{pdb_id}_validation.xml.gz"
    last = None
    for attempt in range(1, retries + 1):
        temporary = target.with_name(f"{target.name}.part.{os.getpid()}.{threading.get_ident()}")
        try:
            with session().get(url, timeout=(15, 90), stream=True) as response:
                if response.status_code == 404:
                    return {"pdb_id": pdb_id, "status": "NOT_AVAILABLE_404", "bytes": 0}
                response.raise_for_status()
                with temporary.open("wb") as handle:
                    for block in response.iter_content(1 << 20):
                        if block:
                            handle.write(block)
            if not valid_gzip(temporary):
                raise RuntimeError("downloaded payload is not a valid validation XML gzip")
            os.replace(temporary, target)
            return {"pdb_id": pdb_id, "status": "DOWNLOADED", "bytes": target.stat().st_size}
        except Exception as exc:
            last = f"{type(exc).__name__}: {exc}"
            temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(2 ** attempt, 20))
    return {"pdb_id": pdb_id, "status": "DOWNLOAD_FAILED", "bytes": 0, "detail": last}


def build_reuse_index(roots: list[Path]) -> dict[str, Path]:
    result = {}
    for root in roots:
        if not root.exists():
            continue
        for path in root.glob("*/*_validation.xml.gz"):
            pdb_id = path.name[:4].lower()
            result.setdefault(pdb_id, path)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-ids", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reuse-root", type=Path, action="append", default=[])
    parser.add_argument("--max-workers", type=int, default=16)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    pdb_ids = sorted({line.strip().lower() for line in args.pdb_ids.read_text().splitlines() if line.strip()})
    if any(len(value) != 4 for value in pdb_ids):
        raise RuntimeError("PDB-ID list contains a non-four-character value")
    reuse = build_reuse_index(args.reuse_root)
    rows = []
    pending = []
    for pdb_id in pdb_ids:
        target = destination(args.output_root, pdb_id)
        if valid_gzip(target):
            rows.append({"pdb_id": pdb_id, "status": "ALREADY_CACHED", "bytes": target.stat().st_size})
        elif pdb_id in reuse and valid_gzip(reuse[pdb_id]):
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(reuse[pdb_id], target)
                mode = "REUSED_HARDLINK"
            except OSError:
                shutil.copy2(reuse[pdb_id], target)
                mode = "REUSED_COPY"
            rows.append({"pdb_id": pdb_id, "status": mode, "bytes": target.stat().st_size})
        else:
            pending.append(pdb_id)

    print(f"population={len(pdb_ids)} reused_or_cached={len(rows)} download_pending={len(pending)}", flush=True)
    completed = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {
            pool.submit(download_one, pdb_id, args.output_root, args.retries): pdb_id
            for pdb_id in pending
        }
        for future in as_completed(futures):
            rows.append(future.result())
            completed += 1
            if completed % 1000 == 0 or completed == len(pending):
                print(f"download_completed={completed}/{len(pending)}", flush=True)

    frame = pd.DataFrame(rows).sort_values("pdb_id")
    if len(frame) != len(pdb_ids) or frame["pdb_id"].duplicated().any():
        raise RuntimeError("validation acquisition manifest does not close over input population")
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.manifest, sep="\t", index=False)
    print(frame["status"].value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
