#!/usr/bin/env python3
import argparse
import csv
import json
from collections import Counter
from pathlib import Path


MANIFEST_FIELDS = [
    "source",
    "remote_url",
    "remote_relative_path",
    "local_path",
    "pdb_id",
    "subdirectory",
    "expected_size",
    "downloaded_size",
    "http_status",
    "resumed",
    "gzip_ok",
    "parse_ok",
    "file_type",
    "status",
    "error",
]


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def index_tsv(path, key="pdb_id"):
    return {row[key].strip().lower(): row for row in read_tsv(path)}


def bool_text(value):
    text = str(value or "").strip().lower()
    if text in {"true", "1", "yes", "success", "ok"}:
        return "true"
    if text in {"false", "0", "no", "failed", "failure"}:
        return "false"
    return ""


def parse_status_text(value):
    text = str(value or "").strip().lower()
    if text == "success":
        return "true"
    if text in {"failed", "failure", "parse_failed"}:
        return "false"
    return ""


def combine_errors(*values):
    parts = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in parts:
            parts.append(text)
    return " | ".join(parts)


def remote_relative(url, marker):
    if marker in url:
        return url.split(marker, 1)[1]
    return ""


def write_manifest(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize(path):
    counts = Counter()
    gzip_counts = Counter()
    parse_counts = Counter()
    ids = set()
    rows = 0
    for row in read_tsv(path):
        rows += 1
        ids.add(row["pdb_id"])
        counts[row["status"] or "(empty)"] += 1
        gzip_counts[row["gzip_ok"] or "(empty)"] += 1
        parse_counts[row["parse_ok"] or "(empty)"] += 1
    return {
        "path": str(path),
        "rows": rows,
        "unique_pdb_ids": len(ids),
        "duplicate_pdb_id_rows": rows - len(ids),
        "status_counts": dict(counts),
        "gzip_ok_counts": dict(gzip_counts),
        "parse_ok_counts": dict(parse_counts),
    }


def build_standard_source(runtime_path, validation_path, parse_statuses, source, file_type, marker):
    validation = index_tsv(validation_path)
    output = []
    for row in read_tsv(runtime_path):
        pdb_id = row["pdb_id"].strip().lower()
        val = validation.get(pdb_id, {})
        url = row.get("url", "") or val.get("url", "")
        parse_status = parse_statuses.get(pdb_id, "")
        parse_ok = parse_status_text(parse_status)
        validation_status = str(val.get("validation_status", "")).strip().lower()
        status = row.get("status", "")
        extra_error = ""
        if validation_status and validation_status != "ok":
            status = "failed_validation"
        elif parse_ok == "false":
            status = "download_succeeded_parse_failed"
            extra_error = f"source_parse_status={parse_status}"
        output.append(
            {
                "source": source,
                "remote_url": url,
                "remote_relative_path": remote_relative(url, marker),
                "local_path": val.get("local_path", "") or row.get("local_path", ""),
                "pdb_id": pdb_id,
                "subdirectory": pdb_id[1:3],
                "expected_size": row.get("expected_size", "") or val.get("file_size", ""),
                "downloaded_size": val.get("file_size", "") or row.get("downloaded_size", ""),
                "http_status": row.get("http_status", ""),
                "resumed": bool_text(row.get("resumed", "")),
                "gzip_ok": bool_text(val.get("gzip_ok", "") or row.get("gzip_ok", "")),
                "parse_ok": parse_ok,
                "file_type": file_type,
                "status": status,
                "error": combine_errors(row.get("error"), val.get("error"), extra_error),
            }
        )
    return output


def build_bundle(inventory_path, qc_path, processed_path):
    qc = index_tsv(qc_path)
    processed = index_tsv(processed_path)
    output = []
    for inv in read_tsv(inventory_path):
        pdb_id = inv["pdb_id"].strip().lower()
        subdir = pdb_id[1:3]
        rel = f"{subdir}/{pdb_id}/{pdb_id}-pdb-bundle.tar.gz"
        url = f"https://files.wwpdb.org/pub/pdb/compatible/pdb_bundle/{rel}"
        qc_row = qc.get(pdb_id, {})
        proc = processed.get(pdb_id, {})
        gzip_ok = "true" if qc_row.get("qc_status", "").strip().lower() == "success" else "false"
        parse_ok = bool_text(proc.get("parse_ok", ""))
        status = proc.get("status", "") or qc_row.get("qc_status", "") or inv.get("inventory_status", "")
        output.append(
            {
                "source": "pdb_bundle",
                "remote_url": url,
                "remote_relative_path": rel,
                "local_path": inv.get("bundle_path", "") or proc.get("source_path", ""),
                "pdb_id": pdb_id,
                "subdirectory": subdir,
                "expected_size": inv.get("file_size", "") or proc.get("file_size", ""),
                "downloaded_size": inv.get("file_size", "") or proc.get("file_size", ""),
                "http_status": "",
                "resumed": "",
                "gzip_ok": gzip_ok,
                "parse_ok": parse_ok,
                "file_type": "pdb_bundle_tar_gz",
                "status": status,
                "error": combine_errors(inv.get("error"), qc_row.get("error"), proc.get("error")),
            }
        )
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive_root", required=True)
    parser.add_argument("--project_root", required=True)
    parser.add_argument("--out_dir", required=True)
    args = parser.parse_args()

    archive = Path(args.archive_root)
    project = Path(args.project_root)
    out_dir = Path(args.out_dir)

    parse_by_source = {"mmcif": {}, "pdb": {}}
    unified = project / "data_stage1_v2_unified_full/unified_structure_manifest.tsv"
    for row in read_tsv(unified):
        pdb_id = row["pdb_id"].strip().lower()
        parse_by_source["mmcif"][pdb_id] = row.get("mmcif_parse_status", "")
        parse_by_source["pdb"][pdb_id] = row.get("legacy_pdb_parse_status", "")

    mmcif_path = out_dir / "mmcif_final_download_manifest.tsv"
    pdb_path = out_dir / "pdb_final_download_manifest.tsv"
    bundle_path = out_dir / "pdb_bundle_final_download_manifest.tsv"

    write_manifest(
        mmcif_path,
        build_standard_source(
            archive / "manifests/mmcif_download_runtime_restart_20260711_112900.tsv",
            archive / "manifests/final_mmcif_validation.tsv",
            parse_by_source["mmcif"],
            "mmcif",
            "mmcif_gz",
            "/divided/mmCIF/",
        ),
    )
    write_manifest(
        pdb_path,
        build_standard_source(
            archive / "manifests/pdb_download_runtime_restart_20260711_112900.tsv",
            archive / "manifests/final_pdb_validation.tsv",
            parse_by_source["pdb"],
            "pdb",
            "legacy_pdb_gz",
            "/divided/pdb/",
        ),
    )
    bundle_manifests = archive / "processed/pdb_bundle_v1/manifests"
    write_manifest(
        bundle_path,
        build_bundle(
            bundle_manifests / "pdb_bundle_inventory.tsv",
            bundle_manifests / "pdb_bundle_tar_qc.tsv",
            bundle_manifests / "pdb_bundle_processed_manifest.tsv",
        ),
    )

    summary = {
        "manifest_fields": MANIFEST_FIELDS,
        "size_semantics": {
            "expected_size": "Runtime expected_size when present; otherwise final validated local file size.",
            "downloaded_size": "Final local file size from validation/inventory.",
        },
        "manifests": [summarize(path) for path in (mmcif_path, pdb_path, bundle_path)],
    }
    summary_path = out_dir / "normalized_download_manifest_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
