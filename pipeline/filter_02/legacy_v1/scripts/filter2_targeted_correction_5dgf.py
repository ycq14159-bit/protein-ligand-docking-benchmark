from __future__ import annotations

import csv
import gzip
import json
import os
from collections import Counter
from pathlib import Path

import filter2_pipeline as p


PDB_ID = "5dgf"


def replace_rows(path: Path, fields: list[str], replacement: list[dict]) -> dict:
    tmp = Path(str(path) + ".correction.tmp")
    old = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as source, gzip.open(tmp, "wt", encoding="utf-8", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        inserted = False
        for row in reader:
            belongs_to_target = row.get("pdb_id") == PDB_ID or row.get("source_component_instance_id", "").startswith(PDB_ID + "|")
            if belongs_to_target:
                old.append(row)
                if not inserted:
                    writer.writerows(replacement)
                    inserted = True
            else:
                writer.writerow(row)
        if not inserted:
            raise RuntimeError(f"{PDB_ID} not found in {path}")
    os.replace(tmp, path)
    return {"path": str(path), "old_rows": len(old), "new_rows": len(replacement)}


def main():
    p.load_globals()
    paths = {row["pdb_id"]: row["mmcif_path"] for row in p.iter_tsv(p.OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz")}
    parsed = p.parse_entry((PDB_ID, paths[PDB_ID]))
    old_route_counts = Counter()
    for row in p.iter_tsv(p.OUT / "full/filter_2_sources.tsv.gz"):
        if row["pdb_id"] == PDB_ID:
            old_route_counts[row["filter_2_route"]] += 1
    new_route_counts = Counter(row["filter_2_route"] for row in parsed["sources"])
    if any(row["filter_2_route"] == "ordinary_small_molecule_candidate" and row["polymer_context"] != "independent_nonpolymer" for row in parsed["sources"]):
        raise RuntimeError("Targeted correction still leaves polymer-context ordinary rows")

    results = []
    for key, fields in p.TABLE_FIELDS.items():
        results.append(replace_rows(p.OUT / f"full/filter_2_{key}.tsv.gz", fields, parsed[key]))

    batch_dir = None
    for candidate in sorted((p.OUT / "checkpoints/batches").glob("batch_*")):
        if any(row["pdb_id"] == PDB_ID for row in p.iter_tsv(candidate / "entries.tsv.gz")):
            batch_dir = candidate
            break
    if batch_dir is None:
        raise RuntimeError("Checkpoint batch for 5dgf not found")
    for key, fields in p.TABLE_FIELDS.items():
        results.append(replace_rows(batch_dir / f"{key}.tsv.gz", fields, parsed[key]))

    audit = {
        "pdb_id": PDB_ID,
        "reason": "other_polymer context was not explicitly routed before CCD small-organic fallback",
        "rule_fix": "other_polymer -> polymer_or_modified_residue; unresolved contexts -> unresolved_review",
        "old_route_counts": dict(old_route_counts),
        "new_route_counts": dict(new_route_counts),
        "affected_source_instance_ids": [row["source_component_instance_id"] for row in parsed["sources"] if row["polymer_context"] == "other_polymer"],
        "rewrites": results,
        "checkpoint_batch": str(batch_dir),
    }
    target = p.OUT / "validation/filter_2_targeted_correction_5dgf.json"
    target.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
