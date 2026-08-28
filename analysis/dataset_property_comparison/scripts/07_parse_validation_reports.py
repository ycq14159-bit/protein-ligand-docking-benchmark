#!/usr/bin/env python3
"""Normalize entry and candidate-ligand fields from cached wwPDB XML reports."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import math
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ENTRY_SCHEMA = pa.schema([
    ("pdb_id", pa.string()), ("resolution", pa.float64()),
    ("parse_status", pa.string()), ("parse_error", pa.string()),
])
RESIDUE_SCHEMA = pa.schema([
    ("pdb_id", pa.string()), ("model_id", pa.string()),
    ("auth_asym_id", pa.string()), ("label_asym_id", pa.string()),
    ("auth_seq_id", pa.string()), ("label_seq_id", pa.string()),
    ("insertion_code", pa.string()), ("component_id", pa.string()),
    ("alt_id", pa.string()), ("rsr", pa.float64()), ("rscc", pa.float64()),
])
CANDIDATE_COMPONENTS = set()


def init_worker(components):
    global CANDIDATE_COMPONENTS
    CANDIDATE_COMPONENTS = set(components)


def local_tag(element):
    return element.tag.split("}")[-1]


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value in {".", "?"} else value


def number(value):
    if value in (None, "", ".", "?", "NotAvailable", "Not available"):
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def parse_one(path_text: str, candidate_components: set[str]):
    path = Path(path_text)
    pdb_id = path.name[:4].lower()
    entry = {"pdb_id": pdb_id, "resolution": None, "parse_status": "PARSE_FAILED", "parse_error": ""}
    rows = []
    try:
        with gzip.open(path, "rb") as handle:
            root = ET.parse(handle).getroot()
        element = next((item for item in root.iter() if local_tag(item) == "Entry"), None)
        if element is None:
            raise ValueError("Entry element missing")
        entry.update({
            "resolution": number(element.get("PDB-resolution")),
            "parse_status": "PARSE_SUCCESS", "parse_error": "",
        })
        for item in root.iter():
            if local_tag(item) != "ModelledSubgroup":
                continue
            attrs = item.attrib
            component = clean(attrs.get("resname"))
            if component not in candidate_components:
                continue
            rows.append({
                "pdb_id": pdb_id, "model_id": clean(attrs.get("model")),
                "auth_asym_id": clean(attrs.get("chain")),
                "label_asym_id": clean(attrs.get("said")),
                "auth_seq_id": clean(attrs.get("resnum")),
                "label_seq_id": clean(attrs.get("seq")),
                "insertion_code": clean(attrs.get("icode")),
                "component_id": component, "alt_id": clean(attrs.get("altcode")),
                "rsr": number(attrs.get("rsr")), "rscc": number(attrs.get("rscc")),
            })
    except Exception as exc:
        entry["parse_error"] = f"{type(exc).__name__}: {exc}"[:1000]
    return entry, rows


def process_batch(batch_id, paths, output_root):
    batch = Path(output_root) / "batches" / f"batch_{batch_id:06d}"
    marker = batch / "_COMPLETE.json"
    if marker.exists():
        return json.loads(marker.read_text())
    batch.mkdir(parents=True, exist_ok=True)
    entries, residues = [], []
    for path in paths:
        entry, rows = parse_one(path, CANDIDATE_COMPONENTS)
        entries.append(entry)
        residues.extend(rows)
    pq.write_table(pa.Table.from_pylist(entries, schema=ENTRY_SCHEMA), batch / "entry_validation.parquet", compression="zstd")
    pq.write_table(pa.Table.from_pylist(residues, schema=RESIDUE_SCHEMA), batch / "ligand_candidate_validation.parquet", compression="zstd")
    summary = {
        "batch_id": batch_id, "pdb_count": len(entries),
        "parse_success": sum(row["parse_status"] == "PARSE_SUCCESS" for row in entries),
        "candidate_residue_rows": len(residues),
    }
    temporary = marker.with_suffix(".tmp")
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    os.replace(temporary, marker)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--prepared-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=200)
    args = parser.parse_args()
    components = set()
    for path in args.prepared_root.glob("*_properties*.parquet"):
        if path.name == "ours_properties.parquet":
            continue
        table = pd.read_parquet(path, columns=["ccd_id"])
        components.update(str(value).strip() for value in table["ccd_id"].dropna())
    paths = sorted(args.cache_root.glob("*/*/*_validation.xml.gz"))
    batches = [paths[i:i + args.batch_size] for i in range(0, len(paths), args.batch_size)]
    args.output_root.mkdir(parents=True, exist_ok=True)
    totals = {"pdb_count": 0, "parse_success": 0, "candidate_residue_rows": 0}
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers, initializer=init_worker, initargs=(sorted(components),)
    ) as pool:
        futures = [
            pool.submit(
                process_batch, index, [str(p) for p in batch], str(args.output_root)
            )
            for index, batch in enumerate(batches)
        ]
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            for key in totals:
                totals[key] += result[key]
            if completed % 50 == 0 or completed == len(futures):
                print(f"batches_completed={completed}/{len(futures)}", flush=True)
    summary = {
        **totals, "input_xml_files": len(paths), "candidate_component_count": len(components),
        "batch_count": len(batches),
    }
    (args.output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
