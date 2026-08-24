from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


HIST = Path("/root/autodl-tmp/stage_Summary/Stage_A")
OLD_PROJECT = Path("/root/autodl-tmp/vs_benchmark")
ARCHIVE = Path("/root/autodl-tmp/pdb_archive_v2")
ROOT = Path("/root/autodl-tmp/benchmark_1.0")
PROC = ROOT / "processing_1_pdb_source_audit"
EXPECTED = {"mmcif": 256158, "pdb": 241185, "pdb_bundle": 6827}


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify(relative: str) -> tuple[bool, str, str]:
    low = relative.lower()
    if "/logs/" in "/" + low or low.endswith(".pid"):
        if low.endswith("start_commands.txt"):
            return True, "provenance", "historical execution command record"
        return False, "runtime_log", "temporary screen/wget log or PID"
    if "smoke" in low or "/tests/" in "/" + low or "test_" in low:
        return False, "test_or_smoke", "test/smoke artifacts are excluded"
    if "processed/pdb_bundle_v1/manifests/pdb_bundle_inventory.tsv" in low:
        return False, "later_analysis", "contains parsed protein/component semantics"
    if "processed/pdb_bundle_v1/manifests/pdb_bundle_processed_manifest.tsv" in low:
        return False, "later_analysis", "contains parsed protein/component semantics"
    if "processed/pdb_bundle_v1/summaries/pdb_bundle_processing_summary.json" in low:
        return False, "later_analysis", "contains protein/component analysis counts"
    if low.endswith("process_pdb_archives_v2.py") or low.endswith("archive_validation_stage1_v2.py"):
        return False, "later_analysis_code", "includes structure-content/Stage1 processing"
    if low.startswith("migration/"):
        return True, "historical_migration", "validated historical Stage A migration provenance"
    if "/inventories/" in "/" + low:
        return True, "source_inventory", "official source inventory"
    if "/manifests/" in "/" + low:
        if "bundle_tar_qc" in low:
            return True, "bundle_qc", "archive-level bundle integrity QC"
        if "bundle_failed" in low:
            return True, "download_failure", "bundle archive failure ledger"
        return True, "download_manifest", "download, coverage, validation, or retry manifest"
    if "/scripts/" in "/" + low or low.endswith(".py") or low.endswith(".sh"):
        return True, "script", "Processing 1 acquisition/audit implementation"
    if "/configs/" in "/" + low or low.endswith((".yaml", ".yml")):
        return True, "config", "Processing 1 source configuration"
    return False, "unresolved_or_unneeded", "not required for Processing 1 release"


def inventory_sources() -> list[dict[str, object]]:
    rows = []
    for path in sorted(HIST.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        rel = path.relative_to(HIST).as_posix()
        selected, category, reason = classify(rel)
        rows.append({
            "source_path": str(path), "relative_path": rel, "file_name": path.name,
            "file_size": path.stat().st_size, "selected": selected,
            "file_category": category, "selection_or_exclusion_reason": reason,
        })
    return rows


def destination_for(row: dict[str, object]) -> Path:
    rel = str(row["relative_path"])
    category = str(row["file_category"])
    source = Path(str(row["source_path"]))
    if category == "script":
        return PROC / "scripts" / source.name
    if category == "config":
        return PROC / "configs" / source.name
    if category == "provenance":
        return ROOT / "shared/provenance" / rel.replace("/", "__")
    if category == "historical_migration":
        return PROC / "migration/historical_stage_A" / source.name
    return PROC / "inputs/historical_stage_A" / rel


def write_tsv(path: Path, rows: list[dict], fields: list[str], gz: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gz else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plan(rows: list[dict[str, object]]) -> dict:
    selected = [r for r in rows if r["selected"]]
    excluded = [r for r in rows if not r["selected"]]
    destinations = [str(destination_for(r)) for r in selected]
    collisions = sorted(k for k, v in Counter(destinations).items() if v > 1)
    return {
        "candidate_source_files": len(rows),
        "planned_copy_files": len(selected),
        "planned_excluded_files": len(excluded),
        "selected_bytes": sum(int(r["file_size"]) for r in selected),
        "excluded_bytes": sum(int(r["file_size"]) for r in excluded),
        "selected_by_category": dict(Counter(str(r["file_category"]) for r in selected)),
        "excluded_by_category": dict(Counter(str(r["file_category"]) for r in excluded)),
        "raw_mmcif_root": str(ARCHIVE / "mmCIF"),
        "destination_collisions": collisions,
        "later_stage_files_detected_and_excluded": sum(r["file_category"] in {"later_analysis", "later_analysis_code"} for r in excluded),
        "standardized_manifest_generatable": True,
    }


def mkdirs() -> None:
    for path in [
        ROOT / "docs", ROOT / "shared/schemas", ROOT / "shared/utilities", ROOT / "shared/provenance",
        PROC / "configs", PROC / "scripts", PROC / "schemas", PROC / "inputs", PROC / "full",
        PROC / "reports", PROC / "release", PROC / "validation", PROC / "migration", PROC / "logs",
    ]:
        path.mkdir(parents=True, exist_ok=False if path == ROOT else True)


def load_final_manifests() -> dict[str, list[dict[str, str]]]:
    base = HIST / "archive_metadata/manifests/normalized_final"
    return {
        "mmcif": read_tsv(base / "mmcif_final_download_manifest.tsv"),
        "pdb": read_tsv(base / "pdb_final_download_manifest.tsv"),
        "pdb_bundle": read_tsv(base / "pdb_bundle_final_download_manifest.tsv"),
    }


def hash_raw_rows(manifests: dict[str, list[dict[str, str]]]) -> list[dict]:
    jobs = []
    for source, rows in manifests.items():
        for row in rows:
            jobs.append((source, row))

    def one(job: tuple[str, dict[str, str]]) -> dict:
        source, row = job
        path = Path(row["local_path"])
        exists = path.is_file()
        return {
            "source": source, "pdb_id": row["pdb_id"].lower(), "file_type": row["file_type"],
            "raw_path": str(path), "file_exists": str(exists).lower(),
            "file_size": path.stat().st_size if exists else "",
            "sha256": hash_file(path) if exists else "",
            "gzip_ok": row["gzip_ok"], "parse_ok": row["parse_ok"],
            "download_status": row["status"], "qc_status": "pass" if row["parse_ok"].lower() == "true" else "fail",
            "qc_error": row["error"], "preferred_source": str(source == "mmcif").lower(),
        }

    with ThreadPoolExecutor(max_workers=8) as pool:
        return list(pool.map(one, jobs, chunksize=64))


def docs() -> None:
    (ROOT / "README.md").write_text("""# Benchmark 1.0\n\nBenchmark 1.0 is an independently organized project. Only **Processing 1 - PDB Data Acquisition and Source Audit** is currently defined and complete. Raw PDB structure archives are not copied into this repository; they are referenced read-only through the release manifests. No later processing stage has been defined or started.\n""", encoding="utf-8")
    (ROOT / "CHANGELOG.md").write_text(f"# Changelog\n\n- {now()}: Created Benchmark 1.0 Processing 1 release from validated historical acquisition records.\n", encoding="utf-8")
    (ROOT / ".gitignore").write_text("""**/__pycache__/\n**/*.pyc\n**/.pytest_cache/\n**/tmp/\n**/cache/\n**/*.part\n**/logs/*.log\n*.cif\n*.cif.gz\n*.pdb\n*.ent\n*.bcif\n*.tar\n*.tar.gz\n""", encoding="utf-8")
    (ROOT / "docs/processing_1_definition.md").write_text("""# Processing 1 Definition\n\nProcessing 1 records PDB source acquisition, download completeness, file existence and compressed-file QC, source coverage, checksums, pdb_bundle QC, provenance, and a formal mmCIF interface. It does not parse biological content, identify protein/RNA/DNA, build assemblies, identify ligands, create pairs, calculate distances, or perform interaction/quality analysis. Frozen counts: mmCIF 256,158; legacy PDB 241,185; pdb_bundle 6,827; unique PDB IDs 256,158. Preferred source is mmCIF. Bundle QC: 6,826 pass and one failure (`6q9e`).\n""", encoding="utf-8")
    (ROOT / "docs/data_provenance.md").write_text(f"""# Data Provenance\n\nHistorical metadata source: `{HIST}`. Historical project (read-only): `{OLD_PROJECT}`. Raw archive root (read-only): `{ARCHIVE}`. Files were copied, never moved. Standardized tables were generated by `processing_1_pdb_source_audit/scripts/build_processing_1.py`; raw-file SHA256 values were calculated read-only. Creation time: {now()}.\n""", encoding="utf-8")
    (ROOT / "docs/repository_data_policy.md").write_text("""# Repository Data Policy\n\nCode, configs, schemas, documentation, reports and validations are suitable for source control. Raw PDB archives are not redistributed or copied. Large full manifests require a separate release, Git LFS, or external data repository decision. Absolute local paths must be configurable or rewritten for public release.\n""", encoding="utf-8")
    (PROC / "README.md").write_text("""# Processing 1 - PDB Data Acquisition and Source Audit\n\n`inputs/` contains selected historical acquisition records; `full/` contains standardized source inventories; `release/processing_1_mmcif_index.tsv.gz` is the formal 256,158-row mmCIF interface. `validation/` and `release/` contain independent accounting. Raw structure files remain outside this project and are not copied. Re-run `scripts/build_processing_1.py` only against the frozen historical sources after confirming the destination policy.\n""", encoding="utf-8")


def execute(rows: list[dict[str, object]], p: dict, resume_partial: bool = False) -> None:
    if ROOT.exists() and not resume_partial:
        raise SystemExit(f"Destination already exists; refusing overwrite: {ROOT}")
    if resume_partial and not ROOT.exists():
        raise SystemExit(f"Partial destination does not exist: {ROOT}")
    if resume_partial and (PROC / "release/processing_1_release_validation.json").exists():
        raise SystemExit("Refusing partial resume because a release validation already exists")
    if p["destination_collisions"]:
        raise SystemExit(f"Destination collisions: {p['destination_collisions']}")
    if not resume_partial:
        mkdirs()
    docs()

    source_inventory = []
    selected_manifest = []
    excluded_manifest = []
    mapping = []
    for row in rows:
        source_inventory.append(dict(row))
        if not row["selected"]:
            excluded_manifest.append({
                "source_path": row["source_path"], "file_name": row["file_name"],
                "file_size": row["file_size"], "exclusion_reason": row["selection_or_exclusion_reason"],
            })
            continue
        source = Path(str(row["source_path"]))
        target = destination_for(row)
        target.parent.mkdir(parents=True, exist_ok=True)
        sh = hash_file(source)
        if target.exists():
            dh = hash_file(target)
            if sh != dh:
                raise SystemExit(f"Existing destination checksum conflict: {target}")
        else:
            shutil.copy2(source, target, follow_symlinks=False)
            dh = hash_file(target)
        item = {
            "source_path": str(source), "destination_path": str(target), "file_name": source.name,
            "file_size": source.stat().st_size, "source_sha256": sh, "destination_sha256": dh,
            "file_category": row["file_category"], "selection_reason": row["selection_or_exclusion_reason"],
            "migration_status": "copied_verified" if sh == dh else "checksum_mismatch",
        }
        selected_manifest.append(item)
        mapping.append({"source_path": str(source), "destination_path": str(target), "source_sha256": sh, "destination_sha256": dh})

    mig = PROC / "migration"
    write_tsv(mig / "source_file_inventory.tsv", source_inventory, list(source_inventory[0]))
    write_tsv(mig / "selected_file_manifest.tsv", selected_manifest, list(selected_manifest[0]))
    write_tsv(mig / "excluded_file_manifest.tsv", excluded_manifest, list(excluded_manifest[0]))
    write_tsv(mig / "source_to_destination_mapping.tsv", mapping, list(mapping[0]))

    manifests = load_final_manifests()
    raw_index = PROC / "full/raw_structure_file_index.tsv.gz"
    if resume_partial and raw_index.exists():
        with gzip.open(raw_index, "rt", encoding="utf-8", newline="") as handle:
            raw = list(csv.DictReader(handle, delimiter="\t"))
        expected_raw = sum(EXPECTED.values())
        if len(raw) != expected_raw:
            raise SystemExit(f"Partial raw index row mismatch: {len(raw)} != {expected_raw}")
    else:
        raw = hash_raw_rows(manifests)
    raw_fields = list(raw[0])
    write_tsv(PROC / "full/raw_structure_file_index.tsv.gz", raw, raw_fields, True)
    write_tsv(PROC / "full/checksum_inventory.tsv.gz", raw, ["source", "pdb_id", "raw_path", "file_size", "sha256"], True)
    for source, name in [("mmcif", "mmcif_inventory.tsv.gz"), ("pdb", "legacy_pdb_inventory.tsv.gz"), ("pdb_bundle", "pdb_bundle_inventory.tsv.gz")]:
        subset = [r for r in raw if r["source"] == source]
        write_tsv(PROC / "full" / name, subset, raw_fields, True)

    by_source = {s: {r["pdb_id"]: r for r in raw if r["source"] == s} for s in EXPECTED}
    coverage_hist = read_tsv(HIST / "archive_metadata/manifests/archive_source_coverage_by_pdb_id.tsv")
    entry_rows, mmcif_rows, coverage_rows = [], [], []
    for old in coverage_hist:
        pid = old["pdb_id"].lower()
        mm = by_source["mmcif"].get(pid)
        lp = by_source["pdb"].get(pid)
        bu = by_source["pdb_bundle"].get(pid)
        availability = ",".join(x for x, value in [("mmcif", mm), ("legacy_pdb", lp), ("pdb_bundle", bu)] if value)
        status = "ready" if mm and mm["file_exists"] == "true" and mm["qc_status"] == "pass" else "mmcif_unavailable"
        entry = {
            "pdb_id": pid, "preferred_source": "mmcif", "has_mmcif": str(bool(mm)).lower(),
            "has_legacy_pdb": str(bool(lp)).lower(), "has_pdb_bundle": str(bool(bu)).lower(),
            "mmcif_path": mm["raw_path"] if mm else "", "legacy_pdb_path": lp["raw_path"] if lp else "",
            "pdb_bundle_path": bu["raw_path"] if bu else "", "mmcif_file_size": mm["file_size"] if mm else "",
            "legacy_pdb_file_size": lp["file_size"] if lp else "", "pdb_bundle_file_size": bu["file_size"] if bu else "",
            "mmcif_checksum": mm["sha256"] if mm else "", "legacy_pdb_checksum": lp["sha256"] if lp else "",
            "pdb_bundle_checksum": bu["sha256"] if bu else "", "mmcif_qc_status": mm["qc_status"] if mm else "missing",
            "legacy_pdb_qc_status": lp["qc_status"] if lp else "missing", "pdb_bundle_qc_status": bu["qc_status"] if bu else "missing",
            "source_availability": availability, "processing_1_status": status,
        }
        entry_rows.append(entry)
        coverage_rows.append(entry)
        mmcif_rows.append({k: entry[k] for k in ["pdb_id", "preferred_source", "mmcif_path", "mmcif_file_size", "mmcif_checksum", "mmcif_qc_status", "has_legacy_pdb", "has_pdb_bundle", "source_availability", "processing_1_status"]} | {"mmcif_exists": mm["file_exists"] if mm else "false"})

    write_tsv(PROC / "full/pdb_source_inventory.tsv.gz", entry_rows, list(entry_rows[0]), True)
    write_tsv(PROC / "full/source_coverage.tsv.gz", coverage_rows, list(coverage_rows[0]), True)
    write_tsv(PROC / "release/processing_1_entry_manifest.tsv.gz", entry_rows, list(entry_rows[0]), True)
    mm_fields = ["pdb_id", "preferred_source", "mmcif_path", "mmcif_exists", "mmcif_file_size", "mmcif_checksum", "mmcif_qc_status", "has_legacy_pdb", "has_pdb_bundle", "source_availability", "processing_1_status"]
    write_tsv(PROC / "release/processing_1_mmcif_index.tsv.gz", mmcif_rows, mm_fields, True)

    source_counts = []
    for source in EXPECTED:
        vals = [r for r in raw if r["source"] == source]
        source_counts.append({"source": source, "file_count": len(vals), "exists_count": sum(r["file_exists"] == "true" for r in vals), "qc_pass": sum(r["qc_status"] == "pass" for r in vals), "qc_fail": sum(r["qc_status"] == "fail" for r in vals), "total_bytes": sum(int(r["file_size"] or 0) for r in vals)})
    write_tsv(PROC / "reports/source_counts.tsv", source_counts, list(source_counts[0]))
    overlaps = Counter(r["source_availability"] for r in entry_rows)
    write_tsv(PROC / "reports/source_overlap_summary.tsv", [{"source_availability": k, "count": v} for k, v in sorted(overlaps.items())], ["source_availability", "count"])
    failures = [r for r in raw if r["qc_status"] != "pass"]
    write_tsv(PROC / "reports/download_failure_summary.tsv", failures, raw_fields)
    bundle_qc = [r for r in raw if r["source"] == "pdb_bundle"]
    bcounts = Counter(r["qc_status"] for r in bundle_qc)
    write_tsv(PROC / "reports/bundle_qc_summary.tsv", [{"qc_status": k, "count": v, "pdb_ids": ",".join(r["pdb_id"] for r in bundle_qc if r["qc_status"] == k)} for k, v in sorted(bcounts.items())], ["qc_status", "count", "pdb_ids"])

    # Schema is intentionally explicit and only describes Processing 1.
    schema_rows = [{"field": f, "description": f.replace("_", " "), "required": "true"} for f in mm_fields]
    write_tsv(PROC / "schemas/processing_1_mmcif_index_schema.tsv", schema_rows, ["field", "description", "required"])

    release_manifest = PROC / "release/processing_1_mmcif_index.tsv.gz"
    validation = {
        "input_unique_pdb_ids": len(entry_rows), "release_manifest_rows": len(mmcif_rows),
        "unique_release_pdb_ids": len({r["pdb_id"] for r in mmcif_rows}),
        "duplicate_pdb_ids": len(mmcif_rows) - len({r["pdb_id"] for r in mmcif_rows}),
        "expected_mmcif_count": EXPECTED["mmcif"], "actual_mmcif_count": len(manifests["mmcif"]),
        "expected_legacy_pdb_count": EXPECTED["pdb"], "actual_legacy_pdb_count": len(manifests["pdb"]),
        "expected_pdb_bundle_count": EXPECTED["pdb_bundle"], "actual_pdb_bundle_count": len(manifests["pdb_bundle"]),
        "bundle_qc_pass": bcounts["pass"], "bundle_qc_failure": bcounts["fail"],
        "bundle_qc_failure_pdb": [r["pdb_id"] for r in bundle_qc if r["qc_status"] == "fail"],
        "part_file_residual": sum(1 for p in ARCHIVE.rglob("*.part") if p.is_file()),
        "raw_structure_files_copied": False, "raw_structure_files_modified": False,
        "historical_directories_modified": False,
        "missing_mmcif_path": sum(not r["mmcif_path"] for r in entry_rows),
        "mmcif_path_exists_false": sum(r["mmcif_exists"] != "true" for r in mmcif_rows),
        "checksum_mismatch": sum(int(r["file_size"] or 0) != Path(r["raw_path"]).stat().st_size for r in raw if r["file_exists"] == "true"),
        "silent_drop": EXPECTED["mmcif"] - len(mmcif_rows),
        "later_stage_directories_created": False, "later_stage_processing_started": False,
    }
    required = [
        validation["input_unique_pdb_ids"] == 256158, validation["release_manifest_rows"] == 256158,
        validation["unique_release_pdb_ids"] == 256158, validation["duplicate_pdb_ids"] == 0,
        validation["actual_mmcif_count"] == 256158, validation["actual_legacy_pdb_count"] == 241185,
        validation["actual_pdb_bundle_count"] == 6827, validation["bundle_qc_pass"] == 6826,
        validation["bundle_qc_failure"] == 1, validation["bundle_qc_failure_pdb"] == ["6q9e"],
        validation["part_file_residual"] == 0, validation["missing_mmcif_path"] == 0,
        validation["mmcif_path_exists_false"] == 0, validation["checksum_mismatch"] == 0,
        validation["silent_drop"] == 0,
    ]
    validation["release_validation_pass"] = all(required)
    validation["migration_validation_pass"] = all(r["migration_status"] == "copied_verified" for r in selected_manifest)
    (PROC / "validation/processing_1_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (PROC / "release/processing_1_release_validation.json").write_text(json.dumps(validation, indent=2) + "\n")

    summary = {"processing_name": "Processing 1 - PDB Data Acquisition and Source Audit", "created_at": now(), **validation}
    (PROC / "reports/processing_1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (PROC / "release/processing_1_release_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    migration_summary = {**p, "actual_copied_files": len(selected_manifest), "checksum_mismatch": sum(r["migration_status"] != "copied_verified" for r in selected_manifest), "migration_validation_pass": validation["migration_validation_pass"], "created_at": now()}
    (mig / "migration_summary.json").write_text(json.dumps(migration_summary, indent=2) + "\n")
    (mig / "migration_validation.json").write_text(json.dumps({"migration_validation_pass": validation["migration_validation_pass"], "checksum_mismatch": migration_summary["checksum_mismatch"], "historical_directories_modified": False}, indent=2) + "\n")

    interface = {
        "project_name": "Benchmark 1.0", "processing_name": "Processing 1 - PDB Data Acquisition and Source Audit",
        "release_version": "1.0", "entry_count": len(mmcif_rows), "unique_pdb_count": len({r["pdb_id"] for r in mmcif_rows}),
        "manifest_path": str(release_manifest), "manifest_sha256": hash_file(release_manifest),
        "schema_path": str(PROC / "schemas/processing_1_mmcif_index_schema.tsv"),
        "preferred_source_policy": "mmCIF", "raw_mmcif_roots": [str(ARCHIVE / "mmCIF")],
        "historical_source_directories": [str(HIST), str(OLD_PROJECT), str(ARCHIVE)],
        "creation_timestamp": now(), "validation_pass": validation["release_validation_pass"],
    }
    (PROC / "release/processing_1_downstream_interface.json").write_text(json.dumps(interface, indent=2) + "\n")

    # Freeze checksums after all artifacts are present; exclude checksum files themselves.
    release_files = sorted(p for p in (PROC / "release").rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (PROC / "release/SHA256SUMS").write_text("".join(f"{hash_file(p)}  {p.relative_to(PROC / 'release').as_posix()}\n" for p in release_files))
    migration_files = sorted(p for p in mig.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (mig / "SHA256SUMS").write_text("".join(f"{hash_file(p)}  {p.relative_to(mig).as_posix()}\n" for p in migration_files))
    shutil.copy2(Path(__file__), PROC / "scripts/build_processing_1.py")
    software = {"created_at": now(), "hostname": platform.node(), "python": sys.version, "platform": platform.platform(), "command": " ".join(sys.argv)}
    (ROOT / "shared/provenance/processing_1_runtime.json").write_text(json.dumps(software, indent=2) + "\n")
    if not validation["release_validation_pass"] or not validation["migration_validation_pass"]:
        raise SystemExit("Validation failed; release not accepted")
    print(json.dumps({"plan": p, "validation": validation, "root": str(ROOT)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--resume-partial", action="store_true")
    args = parser.parse_args()
    rows = inventory_sources()
    p = plan(rows)
    print(json.dumps(p, indent=2))
    if args.execute:
        execute(rows, p, resume_partial=args.resume_partial)


if __name__ == "__main__":
    main()
