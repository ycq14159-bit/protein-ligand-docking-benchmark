#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import processing4_pipeline as p4  # noqa: E402


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def failure_category(reason: object) -> str:
    s = str(reason)
    if "no stereo-preserving descriptor mapping" in s:
        return "NO_STEREO_PRESERVING_MAPPING"
    if "ligand graph identity mismatch" in s:
        return "POST_GENERATION_GRAPH_STEREO_MISMATCH"
    if "ETKDGv3 returned" in s:
        return "ETKDG_EMBEDDING_FAILED"
    if "UFF parameters unavailable" in s:
        return "UFF_PARAMETERS_UNAVAILABLE"
    return s.split(":", 1)[0]


def verify_bucket(run_dir: str, bid: int, case_ids: list[str]) -> dict:
    run = Path(run_dir)
    errors = []
    checked_files = 0
    checked_bytes = 0
    for cid in case_ids:
        root = run / f"output/cases/bucket_{bid:03d}" / cid
        try:
            marker = json.loads((root / "_SUCCESS.json").read_text(encoding="utf-8"))
            if marker.get("status") != "P4_DOCKING_READY":
                raise RuntimeError("success marker status mismatch")
            expected = marker.get("sha256", {})
            if set(expected) != set(p4.READY_FILES):
                raise RuntimeError("success marker file set mismatch")
            for name in p4.READY_FILES:
                path = root / name
                actual = p4.sha256_file(path)
                checked_files += 1
                checked_bytes += path.stat().st_size
                if actual != expected[name]:
                    raise RuntimeError(f"sha256 mismatch: {name}")
            metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
            if metadata.get("processing4_status") != "P4_DOCKING_READY":
                raise RuntimeError("metadata status mismatch")
        except Exception as exc:
            errors.append({"case_id": cid, "error": f"{type(exc).__name__}: {exc}"})
    return {"bucket_id": bid, "cases": len(case_ids), "checked_files": checked_files,
            "checked_bytes": checked_bytes, "errors": errors}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve()
    out = run / "review_rescue_v1"
    statuses = pd.concat(
        [pq.read_table(p).to_pandas() for p in sorted((run / "work/buckets").glob("bucket_*.parquet"))],
        ignore_index=True,
    )
    inventory = pq.read_table(run / "input/case_inventory.parquet").to_pandas()
    audit = pq.read_table(run / "review_audit_v2/case_audit.parquet").to_pandas()
    transitions = pq.read_table(out / "status_transitions.parquet").to_pandas()
    merged = statuses.merge(
        inventory[["case_id", "component_id", "filter3_quality_class", "bucket_id"]],
        on="case_id", how="left", validate="one_to_one",
    )
    merged = merged.merge(
        audit[["case_id", "audit_v2_reason", "raw_graph_id"]],
        on="case_id", how="left", validate="one_to_one",
    )
    tcols = transitions[["case_id", "rescue_rule_applied", "rescue_rule_version",
                         "staging_status", "staging_reason", "previous_status", "previous_reason"]].copy()
    merged = merged.merge(tcols, on="case_id", how="left", validate="one_to_one")
    merged["final_reason_v2"] = ""
    merged.loc[merged["status"].eq("P4_DOCKING_READY"), "final_reason_v2"] = "BASELINE_DOCKING_READY"
    rescued = merged["rescue_rule_applied"].fillna(False).astype(bool)
    merged.loc[rescued, "final_reason_v2"] = "RESCUED|" + merged.loc[rescued, "audit_v2_reason"].astype(str)
    review = merged["status"].eq("P4_PREPARATION_REVIEW")
    merged.loc[review, "final_reason_v2"] = merged.loc[review, "audit_v2_reason"].fillna("OTHER_PREPARATION_REVIEW")
    attempted_failed = review & merged["staging_status"].eq("RESCUE_STAGING_FAILED")
    merged.loc[attempted_failed, "final_reason_v2"] = (
        merged.loc[attempted_failed, "audit_v2_reason"].astype(str) +
        "|RESCUE_REJECTED_" + merged.loc[attempted_failed, "staging_reason"].map(failure_category)
    )
    start_failed = merged["status"].eq("P4_LIGAND_START_GENERATION_FAILED")
    merged.loc[start_failed, "final_reason_v2"] = merged.loc[start_failed, "reason"].map(failure_category)
    pq.write_table(pa.Table.from_pandas(merged, preserve_index=False),
                   out / "final_case_status_with_audit.parquet", compression="zstd")

    status_by_quality = (
        merged.groupby(["filter3_quality_class", "status"]).size().rename("pairs").reset_index()
    )
    status_by_quality.to_csv(out / "final_status_by_filter3_quality.tsv", sep="\t", index=False)
    reason_summary = (
        merged.groupby(["status", "final_reason_v2"], dropna=False)
        .agg(pairs=("case_id", "size"), unique_ccd=("component_id", "nunique"),
             unique_audited_graphs=("raw_graph_id", "nunique"))
        .reset_index().sort_values(["status", "pairs"], ascending=[True, False])
    )
    reason_summary.to_csv(out / "final_reason_v2_summary.tsv", sep="\t", index=False)
    ligand_distribution = (
        merged.groupby(["status", "component_id"], dropna=False).size().rename("pairs").reset_index()
        .sort_values(["status", "pairs"], ascending=[True, False])
    )
    ligand_distribution.to_csv(out / "final_ligand_component_distribution.tsv.gz", sep="\t", index=False, compression="gzip")

    ready = merged[merged["status"].eq("P4_DOCKING_READY")]
    by_bucket = {int(bid): group["case_id"].astype(str).tolist() for bid, group in ready.groupby("bucket_id")}
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(verify_bucket, str(run), bid, ids): bid for bid, ids in by_bucket.items()}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps({k: v for k, v in result.items() if k != "errors"}), flush=True)
    errors = [item for result in results for item in result["errors"]]
    hash_report = {
        "status": "PASS" if not errors and sum(x["cases"] for x in results) == 149521 else "FAIL",
        "ready_cases_expected": 149521,
        "ready_cases_checked": sum(x["cases"] for x in results),
        "delivery_files_checked": sum(x["checked_files"] for x in results),
        "delivery_bytes_checked": sum(x["checked_bytes"] for x in results),
        "hash_error_count": len(errors), "hash_error_examples": errors[:20],
        "validated_at": utc(),
    }
    p4.atomic_json(out / "final_ready_hash_validation.json", hash_report)
    payload = {
        "created_at": utc(),
        "input_total": int(len(merged)),
        "status_counts": {str(k): int(v) for k, v in merged["status"].value_counts().items()},
        "status_percent": {str(k): round(100.0 * int(v) / len(merged), 6)
                           for k, v in merged["status"].value_counts().items()},
        "status_by_filter3_quality": status_by_quality.to_dict("records"),
        "unique_ccd_by_status": {str(k): int(v) for k, v in merged.groupby("status")["component_id"].nunique().items()},
        "rescued_cases": int(rescued.sum()),
        "rescued_by_audit_reason": {str(k): int(v) for k, v in merged[rescued]["audit_v2_reason"].value_counts().items()},
        "rescue_attempt_rejected": int(attempted_failed.sum()),
        "hash_validation": hash_report,
    }
    p4.atomic_json(out / "final_status_summary.json", payload)
    print(json.dumps(payload, indent=2), flush=True)
    if hash_report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
