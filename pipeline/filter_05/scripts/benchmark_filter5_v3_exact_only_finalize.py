#!/usr/bin/env python3
"""Freeze Benchmark Filter 5 v3 after removing all near-similarity deduplication."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


EXPECTED = 91860


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), temporary, compression="zstd")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    stage = root / "filter_05_benchmark_redundancy_reduction"
    old_run = stage / "runs/20260827_exact_near_v2_01"
    run = stage / "runs" / args.run_id
    if run.exists():
        raise RuntimeError(f"refusing to overwrite {run}")
    (run / "final").mkdir(parents=True)
    (run / "validation").mkdir(parents=True)

    old_frozen = json.loads((old_run / "_FROZEN.json").read_text())
    if old_frozen.get("status") != "FROZEN":
        raise RuntimeError("source Filter 5 v2 run is not frozen")

    source_inventory = old_run / "step02_exact_redundancy/output/01_step2_pair_inventory.parquet"
    inventory = pq.read_table(source_inventory).to_pandas()
    if len(inventory) != EXPECTED or inventory["pair_id"].nunique() != EXPECTED:
        raise RuntimeError("source exact-equivalence inventory does not close to 91,860 unique pairs")

    inventory["filter5_v3_final_status"] = ""
    redundant = inventory["step2_exact_role"].eq("STEP2_EXACT_REDUNDANT")
    representative = inventory["step2_exact_role"].eq("STEP2_EXACT_REPRESENTATIVE")
    review = (~redundant & ~representative) & (
        inventory["step1_pair_status"].eq("STEP1_REVIEW")
        | inventory["step2_pair_status"].eq("STEP2_PAIR_MAPPING_REVIEW")
    )
    unique = ~redundant & ~representative & ~review
    inventory.loc[redundant, "filter5_v3_final_status"] = "BENCHMARK_F5_EXACT_REDUNDANT"
    inventory.loc[representative, "filter5_v3_final_status"] = "BENCHMARK_F5_EXACT_REPRESENTATIVE"
    inventory.loc[review, "filter5_v3_final_status"] = "BENCHMARK_F5_REVIEW_RETAIN"
    inventory.loc[unique, "filter5_v3_final_status"] = "BENCHMARK_F5_EXACT_UNIQUE_RETAIN"
    inventory["filter5_v3_membership"] = inventory["filter5_v3_final_status"].ne(
        "BENCHMARK_F5_EXACT_REDUNDANT"
    )
    inventory["near_similarity_deduplication_applied"] = False
    inventory["filter5_policy_version"] = "benchmark_filter5_v3_exact_only"

    retained = inventory[inventory["filter5_v3_membership"]].copy()
    removed = inventory[~inventory["filter5_v3_membership"]].copy()
    counts = inventory["filter5_v3_final_status"].value_counts().sort_index().to_dict()
    quality = retained["filter3_quality_class"].value_counts().sort_index().to_dict()

    # Recompute the essential exact-group representative invariant.
    grouped = inventory[inventory["exact_redundancy_group_id"].ne("")]
    group_check = grouped.groupby("exact_redundancy_group_id").agg(
        members=("pair_id", "size"),
        representatives=("filter5_v3_final_status", lambda s: int((s == "BENCHMARK_F5_EXACT_REPRESENTATIVE").sum())),
    )
    checks = {
        "formal_input_is_91860": bool(len(inventory) == EXPECTED),
        "input_pair_id_unique": bool(inventory["pair_id"].is_unique),
        "status_complete": bool(inventory["filter5_v3_final_status"].ne("").all()),
        "partition_closure": bool(len(retained) + len(removed) == EXPECTED),
        "retained_pair_id_unique": bool(retained["pair_id"].is_unique),
        "one_representative_per_exact_group": bool(group_check["representatives"].eq(1).all()),
        "near_similarity_deduplication_absent": bool((~inventory["near_similarity_deduplication_applied"]).all()),
        "exact_removed_matches_source": bool(len(removed) == 26698),
        "formal_output_is_65162": bool(len(retained) == 65162),
    }
    if not all(checks.values()):
        raise RuntimeError(f"validation failed: {checks}")

    full_path = run / "final/01_benchmark_filter5_v3_full_inventory.parquet"
    retained_path = run / "final/02_benchmark_filter5_v3_retained_cases.parquet"
    removed_path = run / "final/03_benchmark_filter5_v3_exact_redundant_pairs.parquet"
    write_parquet(inventory.sort_values("pair_id"), full_path)
    write_parquet(retained.sort_values("pair_id"), retained_path)
    write_parquet(removed.sort_values("pair_id"), removed_path)

    now = datetime.now(timezone.utc).isoformat()
    summary = {
        "stage": "Benchmark Filter 5 v3: Strict Exact-Equivalent Redocking-Case Deduplication",
        "policy_version": "benchmark_filter5_v3_exact_only",
        "status": "PASS",
        "generated_at": now,
        "formal_input": EXPECTED,
        "input_quality": {
            "FILTER3_HIGH_QUALITY": 7674,
            "FILTER3_GOOD_QUALITY": 84186,
        },
        "exact_groups": 12417,
        "exact_redundant_removed": len(removed),
        "formal_output_retained": len(retained),
        "formal_output_quality": quality,
        "final_status_counts": counts,
        "representative_rank": [
            "resolution ASC", "R_free ASC", "abs_Rfree_Rwork ASC", "R_work ASC", "pair_id ASC"
        ],
        "near_similarity_step": "REMOVED_BY_SCIENTIFIC_POLICY",
        "removed_logic": [
            "receptor similarity clustering",
            "ligand Tanimoto similarity clustering",
            "6A pocket-environment similarity clustering",
        ],
    }
    validation = {
        "status": "PASS",
        "validated_at": now,
        "checks": checks,
        "formal_input": EXPECTED,
        "formal_output": len(retained),
        "duplicate_pair_id": int(inventory["pair_id"].duplicated().sum()),
        "missing_pair_id": int(inventory["pair_id"].astype(str).eq("").sum()),
        "final_status_counts": counts,
        "formal_output_quality": quality,
    }
    provenance = {
        "created_at": now,
        "run_id": args.run_id,
        "source_frozen_run": str(old_run),
        "source_frozen_marker_sha256": sha256(old_run / "_FROZEN.json"),
        "source_exact_inventory": str(source_inventory),
        "source_exact_inventory_sha256": sha256(source_inventory),
        "reused_scientific_evidence": "Filter 5 v2 Steps 1-2 only; scientific logic unchanged",
        "excluded_from_new_policy": "All former Step 3 near-similarity evidence and decisions",
        "old_v2_run_preserved": True,
    }
    write_json(run / "final/summary.json", summary)
    write_json(run / "validation/validation.json", validation)
    write_json(run / "provenance.json", provenance)

    manifest_rows = []
    for path in sorted(p for p in run.rglob("*") if p.is_file() and p.name not in {"_FROZEN.json", "manifest.tsv.gz"}):
        manifest_rows.append((str(path.relative_to(run)), path.stat().st_size, sha256(path)))
    manifest_path = run / "manifest.tsv.gz"
    with gzip.open(manifest_path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("relative_path\tsize_bytes\tsha256\n")
        for rel, size, digest in manifest_rows:
            handle.write(f"{rel}\t{size}\t{digest}\n")

    frozen = {
        "status": "FROZEN",
        "stage": "benchmark_filter5_v3_exact_only",
        "frozen_at": now,
        "formal_input_pairs": EXPECTED,
        "formal_output_pairs": len(retained),
        "formal_status_counts": counts,
        "formal_output_quality": quality,
        "manifest_sha256": sha256(manifest_path),
        "summary_sha256": sha256(run / "final/summary.json"),
        "validation_sha256": sha256(run / "validation/validation.json"),
    }
    write_json(run / "_FROZEN.json", frozen)
    write_json(stage / "CURRENT_RUN.json", {
        "run_id": args.run_id,
        "run_path": str(run),
        "status": "FROZEN",
        "frozen_marker_sha256": sha256(run / "_FROZEN.json"),
        "policy_version": "benchmark_filter5_v3_exact_only",
    })
    print(json.dumps(frozen, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
