#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

MANIFESTS = Path("/root/autodl-tmp/pdb_archive_v2/manifests")


def main() -> int:
    summary_path = MANIFESTS / "final_archive_validation_summary.json"
    missing_path = MANIFESTS / "final_missing_or_corrupt.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = []
    if missing_path.exists() and missing_path.stat().st_size > 2:
        with missing_path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    repair = summary.get("bundle_6q9e_repair", {})
    if repair and not repair.get("final_qc_ok"):
        rows = [r for r in rows if not (r.get("source") == "pdb_bundle" and r.get("pdb_id") == "6q9e")]
        rows.append(
            {
                "pdb_id": "6q9e",
                "source": "pdb_bundle",
                "url": repair.get("official_url", ""),
                "local_path": repair.get("local_path", ""),
                "exists": True,
                "file_size": repair.get("local_size_after") or repair.get("local_size_before") or "",
                "zero_byte": False,
                "gzip_ok": False,
                "html_or_error_page": False,
                "filename_valid": True,
                "subdir_matches_pdb_id": True,
                "duplicate_url": False,
                "duplicate_local_path": False,
                "validation_status": "bundle_tar_qc_failed_after_redownload",
                "error": repair.get("final_qc_error", ""),
            }
        )
    fields = [
        "pdb_id",
        "source",
        "url",
        "local_path",
        "exists",
        "file_size",
        "zero_byte",
        "gzip_ok",
        "html_or_error_page",
        "filename_valid",
        "subdir_matches_pdb_id",
        "duplicate_url",
        "duplicate_local_path",
        "validation_status",
        "error",
    ]
    with missing_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    summary["missing_or_corrupt_count"] = len(rows)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"missing_or_corrupt_count": len(rows), "missing_path": str(missing_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
