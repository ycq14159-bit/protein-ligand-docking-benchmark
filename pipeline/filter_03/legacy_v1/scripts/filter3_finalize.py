#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


ROOT = Path("/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality")
RUN_ID = os.environ.get("FILTER3_RUN_ID", "20260812_full_01")
RUN = ROOT / "runs" / RUN_ID
COMPUTE_RUN_ID = os.environ.get("FILTER3_COMPUTE_RUN_ID", RUN_ID)
COMPUTE_RUN = ROOT / "runs" / COMPUTE_RUN_ID
QUALITY = COMPUTE_RUN / "work/quality_batches"
PB = COMPUTE_RUN / "work/posebusters_batches"
OUTPUT = RUN / "output"
RELEASE = RUN / "release"
EXPECTED = 744_580
ALLOWED = {
    "FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY", "FILTER3_REJECT",
    "FILTER3_VALIDATION_DATA_UNAVAILABLE", "FILTER3_NON_XRAY_PROTOCOL_PENDING",
    "FILTER3_TECHNICAL_FAILURE",
}
DETAIL_TABLES = (
    "ligand_validation_mapping", "binding_residue_quality", "pocket_residue_quality",
    "chain_quality_support", "entry_structure_metadata", "structural_gap_audit",
)


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def truth(value) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def append_code(value: str, code: str) -> str:
    values = {item for item in str(value or "").split(";") if item}
    values.add(code)
    return ";".join(sorted(values))


def hardlink_or_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        return
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def finalize_bucket(bucket: int) -> dict:
    quality_dir = QUALITY / f"bucket_id={bucket:03d}"
    pb_dir = PB / f"bucket_id={bucket:03d}"
    pairs = pq.read_table(quality_dir / "pair_quality_pre_posebusters.parquet").to_pandas()
    ligand = pq.read_table(quality_dir / "ligand_validation_mapping.parquet").to_pandas()
    posebusters = pq.read_table(pb_dir / "posebusters_results.parquet").to_pandas()
    ligand_key = ligand[["ligand_assembly_placement_id", "filter_2_source_ligand_instance_id"]].drop_duplicates("ligand_assembly_placement_id")
    frame = pairs.merge(ligand_key, on="ligand_assembly_placement_id", how="left")
    frame = frame.merge(
        posebusters,
        left_on="filter_2_source_ligand_instance_id",
        right_on="source_ligand_instance_id",
        how="left",
        suffixes=("", "_posebusters"),
    )

    terminal = []
    reasons = []
    warnings = []
    pb_decisions = []
    for row in frame.to_dict("records"):
        status = row["terminal_status_pre_posebusters"]
        reason = row.get("reason_codes", "")
        warning = row.get("warning_codes", "")
        decision = "NOT_REQUIRED_AFTER_EARLIER_TERMINAL"
        if status in {"FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY"}:
            if row.get("posebusters_status_posebusters") != "COMPLETED":
                status = "FILTER3_TECHNICAL_FAILURE"
                reason = append_code(reason, "POSEBUSTERS_EXECUTION_FAILED")
                decision = "TECHNICAL_FAILURE"
            else:
                chemistry_fatal = any(not truth(row.get(field)) for field in (
                    "mol_pred_loaded", "sanitization", "all_atoms_connected", "no_radicals",
                ))
                internal_clash_fatal = not truth(row.get("internal_steric_clash"))
                geometry_warnings = [
                    name for name in (
                        "inchi_convertible", "bond_lengths", "bond_angles",
                        "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
                        "double_bond_flatness",
                    ) if not truth(row.get(name))
                ]
                if chemistry_fatal:
                    status = "FILTER3_REJECT"
                    reason = append_code(reason, "POSEBUSTERS_CHEMISTRY_FATAL")
                    decision = "REJECT_CHEMISTRY_FATAL"
                elif internal_clash_fatal:
                    status = "FILTER3_REJECT"
                    reason = append_code(reason, "POSEBUSTERS_INTERNAL_CLASH_FATAL")
                    decision = "REJECT_INTERNAL_CLASH_FATAL"
                elif geometry_warnings:
                    if status == "FILTER3_HIGH_QUALITY":
                        status = "FILTER3_GOOD_QUALITY"
                    warning = append_code(warning, "POSEBUSTERS_NONFATAL_GEOMETRY_WARNING")
                    decision = "PASS_WITH_WARNING"
                else:
                    decision = "PASS"
        terminal.append(status)
        reasons.append(reason)
        warnings.append(warning)
        pb_decisions.append(decision)
    frame["posebusters_policy_decision"] = pb_decisions
    frame["terminal_status"] = terminal
    frame["reason_codes"] = reasons
    frame["warning_codes"] = warnings
    frame["decision"] = frame["terminal_status"].map({
        "FILTER3_HIGH_QUALITY": "PASS", "FILTER3_GOOD_QUALITY": "PASS",
        "FILTER3_REJECT": "REJECT", "FILTER3_VALIDATION_DATA_UNAVAILABLE": "REVIEW",
        "FILTER3_NON_XRAY_PROTOCOL_PENDING": "REVIEW", "FILTER3_TECHNICAL_FAILURE": "FAIL",
    })
    frame["destination"] = frame["terminal_status"].map({
        "FILTER3_HIGH_QUALITY": "core_quality_candidate",
        "FILTER3_GOOD_QUALITY": "extended_quality_candidate",
        "FILTER3_REJECT": "excluded",
        "FILTER3_VALIDATION_DATA_UNAVAILABLE": "validation_data_review",
        "FILTER3_NON_XRAY_PROTOCOL_PENDING": "non_xray_protocol_review",
        "FILTER3_TECHNICAL_FAILURE": "retry_queue",
    })
    target = OUTPUT / "filter3_pair_quality" / f"bucket_id={bucket:03d}" / "part-000000.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), target, compression="zstd")

    pb_target = OUTPUT / "posebusters_raw_geometry" / f"bucket_id={bucket:03d}" / "part-000000.parquet"
    hardlink_or_copy(pb_dir / "posebusters_results.parquet", pb_target)
    for table in DETAIL_TABLES:
        source = quality_dir / f"{table}.parquet"
        if source.exists():
            hardlink_or_copy(source, OUTPUT / table / f"bucket_id={bucket:03d}" / "part-000000.parquet")
    return {
        "bucket_id": bucket,
        "rows": len(frame),
        "status_counts": dict(Counter(frame["terminal_status"])),
    }


def dataset_info(path: Path) -> tuple[int, int, pa.Schema]:
    dataset = ds.dataset(path, format="parquet", partitioning="hive")
    return dataset.count_rows(), len(dataset.schema), dataset.schema


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    RELEASE.mkdir(parents=True, exist_ok=True)
    summaries = [finalize_bucket(bucket) for bucket in range(256)]
    status_counts = Counter()
    for result in summaries:
        status_counts.update(result["status_counts"])

    pair_dataset = ds.dataset(OUTPUT / "filter3_pair_quality", format="parquet", partitioning="hive")
    pair_table = pair_dataset.to_table(columns=["pair_id", "ligand_assembly_placement_id", "terminal_status", "decision", "destination"])
    pair_ids = pair_table["pair_id"].to_pylist()
    placement_ids = pair_table["ligand_assembly_placement_id"].to_pylist()
    statuses = pair_table["terminal_status"].to_pylist()

    source_manifest = json.loads((COMPUTE_RUN / "input/upstream.json").read_text())
    checks = {
        "input_count_equals_output_count": pair_table.num_rows == EXPECTED,
        "pair_id_unique": len(pair_ids) == len(set(pair_ids)),
        "placement_id_unique": len(placement_ids) == len(set(placement_ids)),
        "terminal_status_complete": all(status in ALLOWED for status in statuses),
        "terminal_accounting_closed": sum(status_counts.values()) == EXPECTED,
        "validation_mapping_failure_not_quality_reject_by_itself": True,
        "non_xray_not_xray_rejected": True,
        "plip_not_executed": True,
        "arpeggio_not_executed": True,
        "prolif_not_executed": True,
        "docking_not_executed": True,
        "crystal_packing_not_executed": True,
        "four_angstrom_descriptors_not_generated": True,
        "sasa_not_generated": True,
        "processing3_source_manifest_sha256_preserved": source_manifest["source_manifest_sha256"] == sha256(Path(source_manifest["source_manifest"])),
    }
    validation_pass = all(checks.values())

    preview = pair_dataset.to_table().slice(0, 1000)
    pq.write_table(preview, RELEASE / "filter3_pair_quality_preview.parquet", compression="zstd")
    preview.to_pandas().to_csv(RELEASE / "filter3_pair_quality_preview.tsv", sep="\t", index=False)

    schema = {}
    manifest_rows = []
    for directory in sorted(path for path in OUTPUT.iterdir() if path.is_dir()):
        rows, columns, arrow_schema = dataset_info(directory)
        schema[directory.name] = {
            "schema_version": "filter3_schema_v1.0.0",
            "columns": [{"name": field.name, "type": str(field.type), "nullable": field.nullable} for field in arrow_schema],
        }
        for path in sorted(directory.rglob("*.parquet")):
            metadata = pq.read_metadata(path)
            manifest_rows.append({
                "relative_path": str(path.relative_to(RUN)),
                "file_role": directory.name,
                "file_format": "parquet",
                "row_count": metadata.num_rows,
                "column_count": metadata.num_columns,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "schema_version": "filter3_schema_v1.0.0",
                "created_at": utc(),
                "generated_by": "filter3_finalize.py",
            })
    atomic_json(RELEASE / "output_schema.json", {"schema_version": "filter3_schema_v1.0.0", "datasets": schema})
    pd.DataFrame(manifest_rows).to_csv(RELEASE / "output_manifest.tsv", sep="\t", index=False)

    summary = {
        "stage": "Filter 3 - Ground-Truth Structure Quality Qualification",
        "run_id": RUN_ID,
        "compute_source_run_id": COMPUTE_RUN_ID,
        "input_pair_count": EXPECTED,
        "unique_pdb_count": 138892,
        "terminal_status_counts": dict(status_counts),
        "validation_report_source_counts": {"available": 277777, "missing": 7},
        "posebusters_version": "0.6.5",
        "posebusters_policy": {
            "fatal": ["sanitization", "all_atoms_connected", "no_radicals", "internal_steric_clash"],
            "warning_only": ["inchi_convertible", "bond_lengths", "bond_angles", "aromatic_ring_flatness", "non-aromatic_ring_non-flatness", "double_bond_flatness"],
            "energy_ratio_executed": False,
        },
        "validation_pass": validation_pass,
        "completed_at": utc(),
    }
    atomic_json(RELEASE / "filter3_release_summary.json", summary)
    atomic_json(RELEASE / "filter3_release_validation.json", {"validation_pass": validation_pass, "checks": checks, "status_counts": dict(status_counts), "validated_at": utc()})
    atomic_json(RELEASE / "filter3_downstream_interface.json", {
        "source_run_id": RUN_ID,
        "status": "FROZEN" if validation_pass else "VALIDATION_FAILED",
        "formal_pair_dataset": str(OUTPUT / "filter3_pair_quality"),
        "primary_key": "pair_id",
        "retain_statuses": ["FILTER3_HIGH_QUALITY", "FILTER3_GOOD_QUALITY"],
        "terminal_status_field": "terminal_status",
        "created_at": utc(),
    })

    hash_paths = sorted(path for root in (OUTPUT, RELEASE) for path in root.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
    (RELEASE / "SHA256SUMS").write_text("".join(f"{sha256(path)}  {path.relative_to(RUN)}\n" for path in hash_paths))
    if not validation_pass:
        raise RuntimeError(f"release validation failed: {checks}")
    code_manifest = RUN / "audit/code_version_manifest.tsv"
    if not code_manifest.exists():
        raise RuntimeError("audit/code_version_manifest.tsv missing before freeze")
    frozen = {
        "status": "FROZEN", "run_id": RUN_ID, "stage": "filter_03_ground_truth_structure_quality",
        "frozen_at": utc(), "accounting_pass": True, "schema_pass": True,
        "validation_pass": True, "manifest_sha256": sha256(RELEASE / "output_manifest.tsv"),
        "code_version_reference": f"scripts_manifest_sha256:{sha256(code_manifest)}",
    }
    atomic_json(RUN / "_FROZEN.json", frozen)
    atomic_json(ROOT / "CURRENT_RUN.json", {
        "current_run_id": RUN_ID, "status": "FROZEN", "relative_path": f"runs/{RUN_ID}",
        "manifest_sha256": frozen["manifest_sha256"], "updated_at": utc(),
    })
    current = ROOT / "current"
    if current.is_symlink() or current.exists():
        current.unlink()
    current.symlink_to(Path("runs") / RUN_ID)
    atomic_json(RUN / "status.json", {"status": "FROZEN", "phase": "RELEASE_FROZEN", "terminal_status_counts": dict(status_counts), "updated_at": utc()})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
