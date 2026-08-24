#!/usr/bin/env python3
import argparse
import concurrent.futures
import gzip
import json
import math
import os
import re
import time
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN = ROOT / "runs/20260812_full_01"
SOURCE = RUN / "work/validation_sources/xml"
OUT = RUN / "work/normalized_validation"
BATCH_SIZE = 100


ENTRY_SCHEMA = pa.schema(
    [
        ("pdb_id", pa.string()),
        ("validation_method_hint", pa.string()),
        ("resolution", pa.float64()),
        ("r_work", pa.float64()),
        ("r_free", pa.float64()),
        ("r_free_minus_r_work", pa.float64()),
        ("clashscore", pa.float64()),
        ("percent_rama_outliers", pa.float64()),
        ("percent_rota_outliers", pa.float64()),
        ("percent_rsrz_outliers", pa.float64()),
        ("bonds_rmsz", pa.float64()),
        ("angles_rmsz", pa.float64()),
        ("validation_schema", pa.string()),
        ("xml_creation_date", pa.string()),
        ("parse_status", pa.string()),
        ("parse_error", pa.string()),
    ]
)


RESIDUE_SCHEMA = pa.schema(
    [
        ("pdb_id", pa.string()),
        ("model_id", pa.string()),
        ("auth_asym_id", pa.string()),
        ("label_asym_id", pa.string()),
        ("auth_seq_id", pa.string()),
        ("label_seq_id", pa.string()),
        ("insertion_code", pa.string()),
        ("component_id", pa.string()),
        ("alt_id", pa.string()),
        ("entity_id", pa.string()),
        ("rsr", pa.float64()),
        ("rscc", pa.float64()),
        ("rsrz", pa.float64()),
        ("mean_occupancy", pa.float64()),
        ("mean_b_factor", pa.float64()),
        ("natoms_eds", pa.int64()),
        ("density_outlier", pa.bool_()),
        ("geometry_outlier", pa.bool_()),
        ("chirality_outlier", pa.bool_()),
        ("clash_outlier", pa.bool_()),
        ("bond_outlier_count", pa.int64()),
        ("angle_outlier_count", pa.int64()),
        ("chirality_outlier_count", pa.int64()),
        ("clash_outlier_count", pa.int64()),
    ]
)


def utc():
    return datetime.now(timezone.utc).isoformat()


def number(value):
    if value in (None, "", ".", "?", "NotAvailable", "Not available"):
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def integer(value):
    parsed = number(value)
    return int(parsed) if parsed is not None else None


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value in {".", "?"} else value


def local_tag(element):
    return element.tag.split("}")[-1]


def method_hint(entry):
    bins = clean(entry.get("percentilebins")).lower()
    attempted = clean(entry.get("attemptedValidationSteps")).lower()
    if "xray" in bins or "xtriage" in attempted or entry.get("PDB-R") is not None:
        return "xray"
    if "em" in bins or "em-" in attempted:
        return "cryo_em"
    if "nmr" in bins:
        return "nmr"
    return "unknown"


def child_outlier_counts(element):
    counts = Counter()
    for child in element.iter():
        if child is element:
            continue
        tag = local_tag(child).lower()
        if "bond" in tag and "outlier" in tag:
            counts["bond"] += 1
        if "angle" in tag and "outlier" in tag:
            counts["angle"] += 1
        if ("chiral" in tag or "stereo" in tag) and "outlier" in tag:
            counts["chirality"] += 1
        if "clash" in tag:
            counts["clash"] += 1
    return counts


def parse_one(path):
    pdb_id = path.name.split("_validation.")[0].lower()
    entry_row = {field.name: None for field in ENTRY_SCHEMA}
    entry_row.update({"pdb_id": pdb_id, "parse_status": "PARSE_FAILED", "parse_error": ""})
    residue_rows = []
    try:
        with gzip.open(path, "rb") as handle:
            root = ET.parse(handle).getroot()
        entry = next((e for e in root.iter() if local_tag(e) == "Entry"), None)
        if entry is None:
            raise ValueError("Entry element missing")
        r_work = number(entry.get("PDB-R"))
        r_free = number(entry.get("PDB-Rfree"))
        entry_row.update(
            {
                "validation_method_hint": method_hint(entry),
                "resolution": number(entry.get("PDB-resolution")),
                "r_work": r_work,
                "r_free": r_free,
                "r_free_minus_r_work": (r_free - r_work) if r_free is not None and r_work is not None else None,
                "clashscore": number(entry.get("clashscore")),
                "percent_rama_outliers": number(entry.get("percent-rama-outliers")),
                "percent_rota_outliers": number(entry.get("percent-rota-outliers")),
                "percent_rsrz_outliers": number(entry.get("percent-RSRZ-outliers")),
                "bonds_rmsz": number(entry.get("bonds_rmsz")),
                "angles_rmsz": number(entry.get("angles_rmsz")),
                "validation_schema": clean(root.attrib.get("{http://www.w3.org/2001/XMLSchema-instance}noNamespaceSchemaLocation")),
                "xml_creation_date": clean(entry.get("XMLcreationDate")),
                "parse_status": "PARSE_SUCCESS",
                "parse_error": "",
            }
        )
        for element in root.iter():
            if local_tag(element) != "ModelledSubgroup":
                continue
            attrs = element.attrib
            outliers = child_outlier_counts(element)
            rsrz = number(attrs.get("rsrz"))
            density_flag = rsrz is not None and rsrz > 2.0
            geometry_flag = bool(outliers["bond"] or outliers["angle"] or outliers["chirality"])
            residue_rows.append(
                {
                    "pdb_id": pdb_id,
                    "model_id": clean(attrs.get("model")),
                    "auth_asym_id": clean(attrs.get("chain")),
                    "label_asym_id": clean(attrs.get("said")),
                    "auth_seq_id": clean(attrs.get("resnum")),
                    "label_seq_id": clean(attrs.get("seq")),
                    "insertion_code": clean(attrs.get("icode")),
                    "component_id": clean(attrs.get("resname")),
                    "alt_id": clean(attrs.get("altcode")),
                    "entity_id": clean(attrs.get("ent")),
                    "rsr": number(attrs.get("rsr")),
                    "rscc": number(attrs.get("rscc")),
                    "rsrz": rsrz,
                    "mean_occupancy": number(attrs.get("avgoccu")),
                    "mean_b_factor": number(attrs.get("owab")),
                    "natoms_eds": integer(attrs.get("NatomsEDS")),
                    "density_outlier": density_flag,
                    "geometry_outlier": geometry_flag,
                    "chirality_outlier": bool(outliers["chirality"]),
                    "clash_outlier": bool(outliers["clash"]),
                    "bond_outlier_count": outliers["bond"],
                    "angle_outlier_count": outliers["angle"],
                    "chirality_outlier_count": outliers["chirality"],
                    "clash_outlier_count": outliers["clash"],
                }
            )
    except Exception as exc:
        entry_row["parse_error"] = f"{type(exc).__name__}:{exc}"[:2000]
    return entry_row, residue_rows


def process_batch(batch_id, paths):
    batch_dir = OUT / "batches" / f"batch_{batch_id:06d}"
    marker = batch_dir / "_COMPLETE.json"
    if marker.exists():
        return json.loads(marker.read_text())
    batch_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    residues = []
    for path_text in paths:
        entry, rows = parse_one(Path(path_text))
        entries.append(entry)
        residues.extend(rows)
    entry_table = pa.Table.from_pylist(entries, schema=ENTRY_SCHEMA)
    residue_table = pa.Table.from_pylist(residues, schema=RESIDUE_SCHEMA)
    pq.write_table(entry_table, batch_dir / "entry_validation.parquet", compression="zstd")
    pq.write_table(residue_table, batch_dir / "residue_validation.parquet", compression="zstd")
    result = {
        "status": "complete",
        "batch_id": batch_id,
        "pdb_count": len(entries),
        "parse_success": sum(row["parse_status"] == "PARSE_SUCCESS" for row in entries),
        "parse_failed": sum(row["parse_status"] != "PARSE_SUCCESS" for row in entries),
        "residue_rows": len(residues),
        "finished_at": utc(),
    }
    temporary = marker.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    os.replace(temporary, marker)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(SOURCE.rglob("*_validation.xml.gz"))
    batches = [paths[index : index + BATCH_SIZE] for index in range(0, len(paths), BATCH_SIZE)]
    started = time.time()
    totals = Counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {}
        iterator = iter(enumerate(batches))
        for _ in range(min(len(batches), args.workers * 2)):
            item = next(iterator, None)
            if item is None:
                break
            batch_id, batch = item
            futures[pool.submit(process_batch, batch_id, [str(p) for p in batch])] = batch_id
        done_count = 0
        while futures:
            done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                futures.pop(future)
                result = future.result()
                done_count += 1
                totals.update({
                    "pdb_count": result["pdb_count"],
                    "parse_success": result["parse_success"],
                    "parse_failed": result["parse_failed"],
                    "residue_rows": result["residue_rows"],
                })
                item = next(iterator, None)
                if item is not None:
                    batch_id, batch = item
                    futures[pool.submit(process_batch, batch_id, [str(p) for p in batch])] = batch_id
                if done_count % 10 == 0 or done_count == len(batches):
                    elapsed = time.time() - started
                    status = {
                        "status": "RUNNING",
                        "phase": "VALIDATION_PARSE",
                        "workers": args.workers,
                        "batch_done": done_count,
                        "batch_total": len(batches),
                        **dict(totals),
                        "runtime_seconds": elapsed,
                        "updated_at": utc(),
                    }
                    path = RUN / "status.json"
                    temporary = path.with_suffix(".json.tmp")
                    temporary.write_text(json.dumps(status, indent=2) + "\n")
                    os.replace(temporary, path)
                    print(json.dumps(status), flush=True)
    status.update({"status": "COMPLETED", "phase": "VALIDATION_PARSE_COMPLETE", "finished_at": utc()})
    path = RUN / "status.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(status, indent=2) + "\n")
    os.replace(temporary, path)


if __name__ == "__main__":
    main()
