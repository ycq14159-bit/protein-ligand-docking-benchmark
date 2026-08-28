from __future__ import annotations

import csv
import hashlib
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1]
REFS = AUDIT / "references"
EXT = AUDIT / "external"
OUT = AUDIT / "output"
QC = AUDIT / "qc"
for _p in (REFS, EXT, OUT, QC, AUDIT / "logs"):
    _p.mkdir(parents=True, exist_ok=True)

CROWN_REPORTED_CATH = 2041
PLINDER_REPORTED_TOTAL = 649915
PLINDER_REPORTED_ION = 22728
PLINDER_REPORTED_ARTIFACT = 18626
OURS_N = 91860


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def calibration_status(path: Path) -> str:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[-1]["status"]

