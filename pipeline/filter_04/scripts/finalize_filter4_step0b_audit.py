#!/usr/bin/env python3
"""Post-run metadata finalizer; does not alter PDB audit rows or classifications."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_04_crystal_packing_influence/audit/crystal_metadata_audit_v1")
EXCLUDED = {"SHA256SUMS", "output_manifest.tsv", "runtime.log"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


audit = pd.read_csv(ROOT / "01_pdb_crystal_metadata_audit.tsv.gz", sep="\t")
columns = [
    "cell_6_fields_complete", "any_spacegroup_identifier_present",
    "sg_hm_alt_present", "sg_hall_present", "sg_it_number_present",
    "symmetry_hm_present", "symmetry_full_hm_present",
    "symmetry_hall_present", "symmetry_it_number_present",
    "explicit_symops_present", "fract_matrix_9_complete",
    "fract_vector_3_complete", "fract_transform_complete",
    "gemmi_parse_success", "gemmi_cell_is_crystal",
    "gemmi_spacegroup_resolvable", "gemmi_ready",
]
coverage = {
    "unique_pdb_count": len(audit),
    "pair_count_sum": int(audit["pair_count"].sum()),
    "counts": {column: int(audit[column].fillna(False).astype(bool).sum()) for column in columns},
    "generated_at": datetime.now(timezone.utc).isoformat(),
}
write_json(ROOT / "tag_coverage_summary.json", coverage)
write_json(ROOT / "postrun_finalize.json", {
    "completed_at": datetime.now(timezone.utc).isoformat(),
    "formal_record_tables_modified": False,
    "classification_modified": False,
    "runtime_log_excluded_from_formal_hashes": True,
    "reason": "runtime.log is an externally redirected stream and changes after in-process hashing",
})

files = sorted(path for path in ROOT.iterdir() if path.is_file() and path.name not in EXCLUDED)
manifest = pd.DataFrame([
    {"relative_path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)}
    for path in files
])
manifest.to_csv(ROOT / "output_manifest.tsv", sep="\t", index=False)

hash_files = sorted(files + [ROOT / "output_manifest.tsv"])
with (ROOT / "SHA256SUMS").open("w", encoding="ascii") as handle:
    for path in hash_files:
        handle.write(f"{sha256(path)}  {path.name}\n")

print(json.dumps(coverage, indent=2, sort_keys=True))
