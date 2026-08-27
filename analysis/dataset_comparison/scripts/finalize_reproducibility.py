#!/usr/bin/env python3
"""Capture/compare deterministic output hashes and write the final checksum set."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EXCLUDE = {"hashes.sha256", "reproducibility_qc.tsv"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inventory(root: Path) -> dict[str, str]:
    result = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        if path.name in EXCLUDE or rel.startswith("logs/") or path.suffix.lower() == ".pdf":
            continue
        result[rel] = digest(path)
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--capture", type=Path)
    ap.add_argument("--compare", type=Path)
    args = ap.parse_args()
    current = inventory(args.root)
    if args.capture:
        args.capture.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return
    if not args.compare:
        ap.error("one of --capture or --compare is required")
    baseline = json.loads(args.compare.read_text(encoding="utf-8"))
    names = sorted(set(baseline) | set(current))
    rows = ["path\tfirst_sha256\tsecond_sha256\tmatch"]
    for name in names:
        rows.append(f"{name}\t{baseline.get(name, '')}\t{current.get(name, '')}\t{baseline.get(name) == current.get(name)}")
    all_match = baseline == current
    rows.append(f"__ALL__\t\t\t{all_match}")
    (args.root / "qc" / "reproducibility_qc.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    if not all_match:
        raise SystemExit("Deterministic rerun hash comparison failed")
    final = inventory(args.root)
    (args.root / "qc" / "hashes.sha256").write_text(
        "".join(f"{value}  {name}\n" for name, value in sorted(final.items())), encoding="utf-8")


if __name__ == "__main__":
    main()
