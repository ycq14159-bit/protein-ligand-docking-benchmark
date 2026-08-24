#!/usr/bin/env python3
import json
from pathlib import Path

run = Path("/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v2/runs/20260804_full_01")
batch_root = run / "work/batches"
complete_files = sorted(batch_root.glob("batch_*/complete.json"))
all_ids = []
batch_ids = []
missing_table_files = []
tables = [
    "source_instances", "source_exclusions", "inorganic_review", "ccd_review",
    "provisional_source_ligands", "ligand_assembly_logical_placements",
    "no_retained_assembly_mapping", "water_exclusion_summary",
    "context_exclusion_summary", "entries",
]
for complete in complete_files:
    payload = json.loads(complete.read_text())
    batch_ids.append(payload["batch_id"])
    all_ids.extend(payload["pdb_ids"])
    for table in tables:
        path = complete.parent / f"{table}.tsv.gz"
        if not path.is_file() or path.stat().st_size == 0:
            missing_table_files.append(str(path))

expected_batches = (248037 + 199) // 200
result = {
    "status_json": json.loads((run / "status.json").read_text()),
    "expected_batches": expected_batches,
    "complete_batches": len(complete_files),
    "batch_ids_unique": len(set(batch_ids)),
    "batch_id_min": min(batch_ids) if batch_ids else None,
    "batch_id_max": max(batch_ids) if batch_ids else None,
    "completed_pdb_rows": len(all_ids),
    "completed_pdb_unique": len(set(all_ids)),
    "duplicate_completed_pdb": len(all_ids) - len(set(all_ids)),
    "missing_formal_batch_files": len(missing_table_files),
    "temporary_or_partial_files": len(list(batch_root.rglob("*.tmp"))) + len(list(batch_root.rglob("*.partial"))),
}
result["validation_pass"] = all([
    result["status_json"].get("status") == "COMPLETED",
    result["complete_batches"] == expected_batches,
    result["batch_ids_unique"] == expected_batches,
    result["batch_id_min"] == 0,
    result["batch_id_max"] == expected_batches - 1,
    result["completed_pdb_rows"] == 248037,
    result["completed_pdb_unique"] == 248037,
    result["missing_formal_batch_files"] == 0,
    result["temporary_or_partial_files"] == 0,
])
(run / "audit/batch_completion_validation.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["validation_pass"] else 1)
