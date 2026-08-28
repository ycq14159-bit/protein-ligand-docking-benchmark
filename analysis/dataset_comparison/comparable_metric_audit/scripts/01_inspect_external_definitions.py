from __future__ import annotations

import json
import platform
import sys

import pandas as pd
import pyarrow
import rdkit

from audit_common import EXT, QC, sha256, write_tsv


def main() -> None:
    files = [
        ("CROWN current metadata", EXT / "crown" / "CROWN_metadata_2026-06.parquet"),
        ("CATH-Plus v4.4 domain list", EXT / "cath" / "cath-domain-list-v4_4_0.txt"),
        ("PLINDER 2024-06/v2 annotation", EXT / "plinder" / "annotation_table_2024-06_v2.parquet"),
        ("PLINDER 2024-04/v1 diagnostic annotation", EXT / "plinder" / "annotation_table_2024-04_v1.parquet"),
        ("PLINDER v0.2.0 classifier", EXT / "plinder" / "ligand_utils.py"),
        ("PLINDER v0.2.0 artifact list", EXT / "plinder" / "artifacts_badlist.csv"),
    ]
    rows = []
    for label, path in files:
        rows.append({"source": label, "path": str(path), "exists": path.exists(),
                     "size_bytes": path.stat().st_size if path.exists() else "",
                     "sha256": sha256(path) if path.exists() else ""})
    write_tsv(QC / "external_file_inventory.tsv", rows)
    env = {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "pandas": pd.__version__, "pyarrow": pyarrow.__version__, "rdkit": rdkit.__version__,
    }
    (QC / "environment.json").write_text(json.dumps(env, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

