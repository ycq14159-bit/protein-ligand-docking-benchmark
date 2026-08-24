#!/usr/bin/env python3
"""Validate the Git/data separation contract without accessing scientific data."""

from __future__ import annotations

from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[2]
BLOCKED = {".parquet", ".cif", ".mmcif", ".pdb", ".sdf", ".mol2", ".sqlite", ".gz", ".jar", ".pem", ".key"}
REQUIRED = [
    "README.md",
    ".gitignore",
    "manifests/data_locations.yaml",
    "manifests/frozen_runs.yaml",
    "manifests/source_import.tsv",
    "docs/DATA_MANAGEMENT.md",
    "docs/VERSIONING.md",
    "docs/REPRODUCIBILITY.md",
]


def main() -> None:
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    blocked = [str(path.relative_to(ROOT)) for path in ROOT.rglob("*") if path.is_file() and ".git" not in path.parts and path.suffix.lower() in BLOCKED]
    result = {
        "validation_pass": not missing and not blocked,
        "missing_required_metadata": missing,
        "blocked_data_or_credential_files": blocked,
    }
    print(json.dumps(result, indent=2))
    if not result["validation_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
