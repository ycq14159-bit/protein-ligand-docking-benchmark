#!/usr/bin/env python3
"""Validate the three official, browser-acquired PDBbind v2020 archives."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

EXPECTED = {
    "PDBbind_v2020_plain_text_index.tar.gz": None,
    "PDBbind_v2020_other_PL.tar.gz": 14127,
    "PDBbind_v2020_refined.tar.gz": 5316,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pdb_ids(tf: tarfile.TarFile) -> set[str]:
    found = set()
    for member in tf.getmembers():
        for token in Path(member.name).parts:
            value = token.lower()
            if len(value) == 4 and value[0].isdigit() and value.isalnum():
                found.add(value)
    return found


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = {"dataset": "PDBbind", "version": "v2020", "archives": {}}
    sets = {}
    for filename, expected in EXPECTED.items():
        path = args.root / filename
        if not path.exists():
            raise FileNotFoundError(path)
        with tarfile.open(path, "r:gz") as tf:
            ids = pdb_ids(tf)
            members = len(tf.getmembers())
        result["archives"][filename] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "tar_members": members,
            "pdb_ids": len(ids),
        }
        if expected is not None and len(ids) != expected:
            raise RuntimeError(f"{filename}: expected {expected} PDB IDs, found {len(ids)}")
        sets[filename] = ids
    other = sets["PDBbind_v2020_other_PL.tar.gz"]
    refined = sets["PDBbind_v2020_refined.tar.gz"]
    result["population"] = {
        "other": len(other),
        "refined": len(refined),
        "intersection": len(other & refined),
        "union": len(other | refined),
        "expected_union": 19443,
    }
    result["status"] = "PASS" if not (other & refined) and len(other | refined) == 19443 else "FAIL"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["population"], indent=2))


if __name__ == "__main__":
    main()
