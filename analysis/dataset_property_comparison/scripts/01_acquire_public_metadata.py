#!/usr/bin/env python3
"""Acquire only metadata needed for the harmonized six-dataset comparison.

PDBbind is intentionally excluded: it must be downloaded once through the
user's licensed browser session and ingested by 02_ingest_pdbbind.py.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from datetime import date, timedelta
from pathlib import Path

HIQBIND_REV = "2922a2a95cbb0dd4144a3dab572be90b77ef08c2"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, target: Path, retries: int = 6) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".partial")
    for attempt in range(retries):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "benchmark-metadata-acquisition/1.0"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                append = offset > 0 and response.status == 206
                if offset and not append:
                    offset = 0
                with partial.open("ab" if append else "wb") as out:
                    while True:
                        block = response.read(8 << 20)
                        if not block:
                            break
                        out.write(block)
            partial.replace(target)
            return
        except (OSError, urllib.error.URLError) as exc:
            if attempt + 1 == retries:
                raise RuntimeError(f"download failed: {url}") from exc
            time.sleep(min(60, 2**attempt))


def sources(root: Path):
    hiq = f"https://raw.githubusercontent.com/THGLab/HiQBind/{HIQBIND_REV}/figshare"
    yield "CROWN", "2026-06", "CROWN_metadata.parquet", (
        "https://zenodo.org/api/records/20825315/files/"
        "CROWN_metadata.parquet/content"
    ), root / "crown_202606" / "CROWN_metadata.parquet", "CC BY 4.0"
    for name in ("hiqbind_metadata.csv", "hiqbind_sm_metadata.csv", "README.md"):
        yield "HiQBind", f"Figshare-v3+git-{HIQBIND_REV}", name, (
            f"{hiq}/{name}"
        ), root / "hiqbind_v3" / name, "CC BY 4.0 / MIT code"
    yield "PLINDER", "2024-06/v2", "annotation_table.parquet", (
        "https://storage.googleapis.com/plinder/2024-06/v2/index/"
        "annotation_table.parquet"
    ), root / "plinder_2024-06_v2" / "annotation_table.parquet", "Apache-2.0"
    yield "BioLiP2/Q-BioLiP", "base-before-2026-01-02", "PL_annotation.csv", (
        "https://yanglab.qd.sdu.edu.cn/BioLiP2/DATA/application/PL/"
        "PL_annotation.csv"
    ), root / "biolip2_20260626" / "PL_annotation_base_before_20260102.csv", "official download"
    yield "BioLiP2/Q-BioLiP", "base-before-2026-01-02", "PL_nonredund_annotation.csv", (
        "https://yanglab.qd.sdu.edu.cn/BioLiP2/DATA/application/PL/"
        "PL_nonredund_annotation.csv"
    ), root / "biolip2_20260626" / "PL_nonredund_annotation_base_before_20260102.csv", "official download"
    yield "BioLiP2/Q-BioLiP", "dictionary-acquired-with-snapshot", "ligand.json", (
        "https://yanglab.qd.sdu.edu.cn/Q-BioLiP/ligand/ligand.json"
    ), root / "biolip2_20260626" / "ligand.json", "official download"
    current = date(2026, 1, 7)
    end = date(2026, 6, 24)
    while current <= end:
        stamp = current.strftime("%Y%m%d")
        name = f"Q-BioLiP-{stamp}.csv"
        yield "BioLiP2/Q-BioLiP", "weekly-through-2026-06-26", name, (
            f"https://yanglab.qd.sdu.edu.cn/BioLiP2/DATA/annotation/{name}"
        ), root / "biolip2_20260626" / "weekly" / name, "official download"
        current += timedelta(days=7)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--all", action="store_true", help="compatibility flag")
    args = parser.parse_args()
    manifest = []
    for dataset, version, filename, url, target, license_name in sources(args.root):
        status = "REUSED" if target.exists() and target.stat().st_size else "DOWNLOADED"
        if status == "DOWNLOADED":
            download(url, target)
        manifest.append({
            "dataset": dataset,
            "version": version,
            "filename": filename,
            "url": url,
            "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "file_size": target.stat().st_size,
            "sha256": sha256(target),
            "license": license_name,
            "status": status,
        })
    manifest_path = args.root / "manifests" / "download_manifest_public.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=manifest[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    print(json.dumps({"files": len(manifest), "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
