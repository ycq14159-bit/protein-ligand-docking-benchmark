#!/usr/bin/env python3
"""Filter 4 Step 0B: read-only crystallographic lattice metadata audit.

Audit-class precedence:
1. MMCIF_PARSE_ERROR: raw mmCIF cannot be parsed by Gemmi.
2. DIRECT_READY: all six raw cell tags are valid, at least one raw space-group
   identifier is present, and Gemmi yields a crystal cell plus resolvable group.
3. GEMMI_READY_TAG_INCOMPLETE: Gemmi is ready although direct raw tags are incomplete.
4. INCONSISTENT_METADATA: direct raw cell and space-group metadata look complete,
   but Gemmi cannot produce a fully usable crystallographic result.
5. CELL_PROBLEM / SPACEGROUP_PROBLEM / CELL_AND_SPACEGROUP_PROBLEM according to
   the two independently tested Gemmi readiness dimensions.

No metadata recovery, symmetry expansion, neighbour search, or pair filtering is
performed by this program.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import tarfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import pandas as pd
import pyarrow.dataset as ds
import yaml


ALLOWED_CLASSES = {
    "DIRECT_READY",
    "GEMMI_READY_TAG_INCOMPLETE",
    "CELL_PROBLEM",
    "SPACEGROUP_PROBLEM",
    "CELL_AND_SPACEGROUP_PROBLEM",
    "MMCIF_PARSE_ERROR",
    "INCONSISTENT_METADATA",
}

CELL_TAGS = {
    "cell_length_a": "_cell.length_a",
    "cell_length_b": "_cell.length_b",
    "cell_length_c": "_cell.length_c",
    "cell_angle_alpha": "_cell.angle_alpha",
    "cell_angle_beta": "_cell.angle_beta",
    "cell_angle_gamma": "_cell.angle_gamma",
}

SPACEGROUP_TAGS = {
    "sg_hm_alt": "_space_group.name_H-M_alt",
    "sg_hall": "_space_group.name_Hall",
    "sg_it_number": "_space_group.IT_number",
    "symmetry_hm": "_symmetry.space_group_name_H-M",
    "symmetry_full_hm": "_symmetry.pdbx_full_space_group_name_H-M",
    "symmetry_hall": "_symmetry.space_group_name_Hall",
    "symmetry_it_number": "_symmetry.Int_Tables_number",
}

FRACT_MATRIX_TAGS = [f"_atom_sites.fract_transf_matrix[{i}][{j}]" for i in range(1, 4) for j in range(1, 4)]
FRACT_VECTOR_TAGS = [f"_atom_sites.fract_transf_vector[{i}]" for i in range(1, 4)]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def clean(value) -> str:
    return "" if value is None else str(value).strip()


def scalar_status(value, kind: str) -> tuple[str, float | None, str]:
    text = clean(value)
    if not text:
        return "MISSING", None, text
    if text == "?":
        return "UNKNOWN_QUESTION_MARK", None, text
    if text == ".":
        return "NOT_APPLICABLE_DOT", None, text
    try:
        number = float(text)
    except ValueError:
        return "NON_NUMERIC", None, text
    if not math.isfinite(number):
        return "INVALID_VALUE", None, text
    valid = number > 0 if kind == "length" else 0 < number < 180
    return ("PRESENT_VALID" if valid else "INVALID_VALUE"), number, text


def numeric_complete(block, tags: list[str]) -> tuple[bool, int]:
    valid = 0
    for tag in tags:
        text = clean(block.find_value(tag))
        if text in {"", "?", "."}:
            continue
        try:
            if math.isfinite(float(text)):
                valid += 1
        except ValueError:
            pass
    return valid == len(tags), valid


def values(block, tag: str) -> list[str]:
    try:
        return [clean(value) for value in block.find_values(tag) if clean(value) not in {"", "?", "."}]
    except Exception:
        return []


def spacegroup_resolvable(name: str) -> bool:
    if not clean(name):
        return False
    try:
        return gemmi.find_spacegroup_by_name(clean(name)) is not None
    except Exception:
        return False


def cryst1_from_lines(lines) -> dict:
    result = {
        "cryst1_present": False,
        "cryst1_cell_valid": False,
        "cryst1_spacegroup": "",
        "cryst1_spacegroup_nonempty": False,
        "cryst1_error": "",
    }
    for raw in lines:
        line = raw.decode("ascii", errors="replace") if isinstance(raw, bytes) else raw
        if not line.startswith("CRYST1"):
            continue
        result["cryst1_present"] = True
        try:
            a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
            alpha, beta, gamma = float(line[33:40]), float(line[40:47]), float(line[47:54])
            result["cryst1_cell_valid"] = (
                all(math.isfinite(x) for x in (a, b, c, alpha, beta, gamma))
                and a > 0 and b > 0 and c > 0
                and all(0 < x < 180 for x in (alpha, beta, gamma))
            )
            result.update({
                "cryst1_a": a, "cryst1_b": b, "cryst1_c": c,
                "cryst1_alpha": alpha, "cryst1_beta": beta, "cryst1_gamma": gamma,
            })
        except Exception as exc:
            result["cryst1_error"] = f"{type(exc).__name__}: {exc}"
        result["cryst1_spacegroup"] = clean(line[55:66])
        result["cryst1_spacegroup_nonempty"] = bool(result["cryst1_spacegroup"])
        break
    return result


def local_cryst1(pdb_id: str, legacy_root: str, bundle_root: str) -> dict:
    legacy = Path(legacy_root) / pdb_id[1:3] / f"pdb{pdb_id}.ent.gz"
    bundle = Path(bundle_root) / pdb_id / f"{pdb_id}-pdb-bundle.tar.gz"
    base = {
        "local_legacy_pdb_available": legacy.exists(),
        "local_bundle_available": bundle.exists(),
        "local_pdb_available": legacy.exists() or bundle.exists(),
        "cryst1_source": "",
        "cryst1_present": False,
        "cryst1_cell_valid": False,
        "cryst1_spacegroup": "",
        "cryst1_spacegroup_nonempty": False,
        "cryst1_error": "",
    }
    if legacy.exists():
        try:
            with gzip.open(legacy, "rt", encoding="ascii", errors="replace") as handle:
                parsed = cryst1_from_lines(handle)
            base.update(parsed)
            if parsed["cryst1_present"]:
                base["cryst1_source"] = "legacy_pdb"
                return base
        except Exception as exc:
            base["cryst1_error"] = f"legacy {type(exc).__name__}: {exc}"
    if bundle.exists():
        try:
            with tarfile.open(bundle, "r:gz") as archive:
                for member in archive.getmembers():
                    if not member.isfile() or not member.name.lower().endswith((".pdb", ".ent")):
                        continue
                    handle = archive.extractfile(member)
                    if handle is None:
                        continue
                    parsed = cryst1_from_lines(handle)
                    if parsed["cryst1_present"]:
                        base.update(parsed)
                        base["cryst1_source"] = f"pdb_bundle:{member.name}"
                        return base
        except Exception as exc:
            base["cryst1_error"] = f"bundle {type(exc).__name__}: {exc}"
    return base


def audit_one(task: tuple[str, int, str, str, str]) -> dict:
    pdb_id, pair_count, mmcif_root, legacy_root, bundle_root = task
    path = Path(mmcif_root) / pdb_id[1:3] / f"{pdb_id}.cif.gz"
    row = {
        "pdb_id": pdb_id,
        "pair_count": pair_count,
        "mmcif_path": str(path),
        "mmcif_exists": path.exists(),
        "mmcif_parse_success": False,
        "mmcif_parse_error": "",
    }
    for field in CELL_TAGS:
        row[field] = None
        row[f"{field}_raw"] = ""
        row[f"{field}_status"] = "MISSING"
    for field in SPACEGROUP_TAGS:
        row[f"{field}_value"] = ""
        row[f"{field}_present"] = False
    row.update({
        "cell_6_fields_complete": False,
        "any_spacegroup_identifier_present": False,
        "space_group_symop_present": False,
        "space_group_symop_count": 0,
        "symmetry_equiv_present": False,
        "symmetry_equiv_count": 0,
        "explicit_symops_present": False,
        "explicit_symops_count": 0,
        "fract_matrix_9_complete": False,
        "fract_matrix_valid_count": 0,
        "fract_vector_3_complete": False,
        "fract_vector_valid_count": 0,
        "fract_transform_complete": False,
        "gemmi_parse_success": False,
        "gemmi_parse_error": "",
        "gemmi_cell_a": None, "gemmi_cell_b": None, "gemmi_cell_c": None,
        "gemmi_cell_alpha": None, "gemmi_cell_beta": None, "gemmi_cell_gamma": None,
        "gemmi_cell_is_crystal": False,
        "gemmi_spacegroup_hm": "",
        "gemmi_spacegroup_nonempty": False,
        "gemmi_spacegroup_resolvable": False,
        "gemmi_ready": False,
    })
    block = None
    if not path.exists():
        row["mmcif_parse_error"] = "MISSING_FROZEN_MMCIF"
    else:
        try:
            block = gemmi.cif.read_file(str(path)).sole_block()
            row["mmcif_parse_success"] = True
            for field, tag in CELL_TAGS.items():
                kind = "length" if "length" in field else "angle"
                status, number, raw = scalar_status(block.find_value(tag), kind)
                row[field], row[f"{field}_raw"], row[f"{field}_status"] = number, raw, status
            row["cell_6_fields_complete"] = all(row[f"{field}_status"] == "PRESENT_VALID" for field in CELL_TAGS)
            for field, tag in SPACEGROUP_TAGS.items():
                value = clean(block.find_value(tag))
                present = value not in {"", "?", "."}
                row[f"{field}_value"], row[f"{field}_present"] = value, present
            row["any_spacegroup_identifier_present"] = any(row[f"{field}_present"] for field in SPACEGROUP_TAGS)
            modern, old = values(block, "_space_group_symop.operation_xyz"), values(block, "_symmetry_equiv.pos_as_xyz")
            row.update({
                "space_group_symop_present": bool(modern), "space_group_symop_count": len(modern),
                "symmetry_equiv_present": bool(old), "symmetry_equiv_count": len(old),
                "explicit_symops_present": bool(modern or old), "explicit_symops_count": len(modern) + len(old),
            })
            matrix_complete, matrix_count = numeric_complete(block, FRACT_MATRIX_TAGS)
            vector_complete, vector_count = numeric_complete(block, FRACT_VECTOR_TAGS)
            row.update({
                "fract_matrix_9_complete": matrix_complete, "fract_matrix_valid_count": matrix_count,
                "fract_vector_3_complete": vector_complete, "fract_vector_valid_count": vector_count,
                "fract_transform_complete": matrix_complete and vector_complete,
            })
        except Exception as exc:
            row["mmcif_parse_error"] = f"{type(exc).__name__}: {exc}"
        try:
            structure = gemmi.read_structure(str(path))
            row["gemmi_parse_success"] = True
            cell = structure.cell
            row.update({
                "gemmi_cell_a": cell.a, "gemmi_cell_b": cell.b, "gemmi_cell_c": cell.c,
                "gemmi_cell_alpha": cell.alpha, "gemmi_cell_beta": cell.beta, "gemmi_cell_gamma": cell.gamma,
                "gemmi_cell_is_crystal": bool(cell.is_crystal()),
                "gemmi_spacegroup_hm": clean(structure.spacegroup_hm),
            })
            row["gemmi_spacegroup_nonempty"] = bool(row["gemmi_spacegroup_hm"])
            row["gemmi_spacegroup_resolvable"] = spacegroup_resolvable(row["gemmi_spacegroup_hm"])
            row["gemmi_ready"] = row["gemmi_cell_is_crystal"] and row["gemmi_spacegroup_resolvable"]
        except Exception as exc:
            row["gemmi_parse_error"] = f"{type(exc).__name__}: {exc}"

    if not row["mmcif_parse_success"] or not row["gemmi_parse_success"]:
        audit_class = "MMCIF_PARSE_ERROR"
    elif row["gemmi_ready"]:
        audit_class = "DIRECT_READY" if row["cell_6_fields_complete"] and row["any_spacegroup_identifier_present"] else "GEMMI_READY_TAG_INCOMPLETE"
    elif row["cell_6_fields_complete"] and row["any_spacegroup_identifier_present"]:
        audit_class = "INCONSISTENT_METADATA"
    elif not row["gemmi_cell_is_crystal"] and row["gemmi_spacegroup_resolvable"]:
        audit_class = "CELL_PROBLEM"
    elif row["gemmi_cell_is_crystal"] and not row["gemmi_spacegroup_resolvable"]:
        audit_class = "SPACEGROUP_PROBLEM"
    else:
        audit_class = "CELL_AND_SPACEGROUP_PROBLEM"
    row["audit_class"] = audit_class

    if not row["gemmi_ready"]:
        row.update(local_cryst1(pdb_id, legacy_root, bundle_root))
    else:
        row.update({
            "local_legacy_pdb_available": False, "local_bundle_available": False,
            "local_pdb_available": False, "cryst1_source": "", "cryst1_present": False,
            "cryst1_cell_valid": False, "cryst1_spacegroup": "",
            "cryst1_spacegroup_nonempty": False, "cryst1_error": "",
        })
    row["cryst1_recovery_candidate"] = bool(row["cryst1_cell_valid"] and row["cryst1_spacegroup_nonempty"])
    return row


def report_table(audit: pd.DataFrame, pair_total: int) -> list[dict]:
    total_pdb = len(audit)
    rows = []
    for name in sorted(ALLOWED_CLASSES):
        subset = audit[audit["audit_class"] == name]
        n_pairs = int(subset["pair_count"].sum())
        rows.append({
            "audit_class": name,
            "n_pdb": len(subset),
            "pct_pdb": 100.0 * len(subset) / total_pdb if total_pdb else 0.0,
            "n_pairs": n_pairs,
            "pct_pairs": 100.0 * n_pairs / pair_total if pair_total else 0.0,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    input_cfg, runtime_cfg = config["input"], config["runtime"]
    output = Path(config["output"]["audit_dir"])
    if args.limit:
        output = Path("/tmp/filter4_step0b_preflight")
    output.mkdir(parents=True, exist_ok=True)
    started = time.time()

    dataset_path = Path(input_cfg["filter3_dataset"])
    # bucket_id is stored both in the files and in Hive-style directory names,
    # with different integer widths. Read the authoritative in-file schema only.
    dataset = ds.dataset(dataset_path, format="parquet")
    required = ["pair_id", "pdb_id", "filter3_v2_terminal_status"]
    missing_columns = sorted(set(required) - set(dataset.schema.names))
    if missing_columns:
        raise RuntimeError(f"Missing Filter 3 columns: {missing_columns}")
    frame = dataset.to_table(columns=required).to_pandas()
    retained = frame[frame["filter3_v2_terminal_status"].isin(input_cfg["retained_statuses"])].copy()
    expected = int(input_cfg["expected_pair_count"])
    if not args.limit and len(retained) != expected:
        raise RuntimeError(f"Expected {expected} retained pairs, observed {len(retained)}")
    if retained["pair_id"].isna().any() or retained["pdb_id"].isna().any():
        raise RuntimeError("Retained Filter 3 row has missing pair_id or pdb_id")
    if retained["pair_id"].duplicated().any():
        raise RuntimeError("Duplicate pair_id in retained Filter 3 input")
    retained["pdb_id"] = retained["pdb_id"].astype(str).str.lower()
    counts = retained.groupby("pdb_id", as_index=False).agg(pair_count=("pair_id", "size")).sort_values("pdb_id")
    if args.limit:
        counts = counts.head(args.limit).copy()
        retained = retained[retained["pdb_id"].isin(counts["pdb_id"])].copy()
    pair_total = len(retained)

    manifest = Path(input_cfg["filter3_release_manifest"])
    input_snapshot = {
        "created_at": utc(),
        "filter3_dataset": str(dataset_path),
        "filter3_release_manifest": str(manifest),
        "filter3_release_manifest_sha256": sha256(manifest),
        "selection_field": "filter3_v2_terminal_status",
        "selection_values": input_cfg["retained_statuses"],
        "pair_count": pair_total,
        "unique_pdb_count": len(counts),
        "config_path": str(config_path),
        "config_sha256": sha256(config_path),
    }
    atomic_json(output / "input_snapshot.json", input_snapshot)

    tasks = [
        (row.pdb_id, int(row.pair_count), input_cfg["mmcif_root"], input_cfg["legacy_pdb_root"], input_cfg["pdb_bundle_root"])
        for row in counts.itertuples(index=False)
    ]
    results = []
    workers = min(int(runtime_cfg["workers"]), os.cpu_count() or 1)
    progress_every = int(runtime_cfg["progress_every_pdb"])
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for index, result in enumerate(executor.map(audit_one, tasks, chunksize=8), start=1):
            results.append(result)
            if index % progress_every == 0 or index == len(tasks):
                atomic_json(output / "progress.json", {
                    "status": "RUNNING", "processed_pdb": index, "total_pdb": len(tasks),
                    "processed_pairs": sum(item["pair_count"] for item in results),
                    "class_counts": dict(Counter(item["audit_class"] for item in results)),
                    "runtime_seconds": time.time() - started, "updated_at": utc(),
                })

    audit = pd.DataFrame(results).sort_values("pdb_id").reset_index(drop=True)
    impact = audit[["pdb_id", "pair_count", "audit_class", "gemmi_ready", "cryst1_recovery_candidate"]].copy()
    problems = audit[~audit["gemmi_ready"]].copy()
    audit.to_csv(output / "01_pdb_crystal_metadata_audit.tsv.gz", sep="\t", index=False, compression="gzip")
    impact.to_csv(output / "02_pair_impact_by_crystal_metadata.tsv.gz", sep="\t", index=False, compression="gzip")
    problems.to_csv(output / "03_problem_pdbs.tsv", sep="\t", index=False)

    stats = report_table(audit, pair_total)
    not_ready = audit[~audit["gemmi_ready"]]
    checks = {
        "filter3_retained_pair_count": pair_total == (expected if not args.limit else pair_total),
        "pair_count_per_pdb_sum": int(audit["pair_count"].sum()) == pair_total,
        "pdb_audit_row_count": len(audit) == len(counts),
        "pdb_id_unique": not audit["pdb_id"].duplicated().any(),
        "audit_class_complete": not audit["audit_class"].isna().any(),
        "audit_class_allowed": set(audit["audit_class"]) <= ALLOWED_CLASSES,
        "each_pdb_exactly_one_class": len(audit) == audit["pdb_id"].nunique(),
        "no_pair_level_filtering": True,
        "no_recovery_executed": True,
        "no_symmetry_mates_generated": True,
        "no_network_download": True,
        "filter3_release_manifest_unchanged": sha256(manifest) == input_snapshot["filter3_release_manifest_sha256"],
    }
    validation = {
        "validation_pass": all(checks.values()),
        "checks": checks,
        "pair_count": pair_total,
        "unique_pdb_count": len(audit),
        "audit_class_statistics": stats,
        "gemmi_ready_pdb_count": int(audit["gemmi_ready"].sum()),
        "gemmi_not_ready_pdb_count": len(not_ready),
        "not_ready_with_explicit_symops": int(not_ready["explicit_symops_present"].sum()),
        "not_ready_with_fract_matrix": int(not_ready["fract_matrix_9_complete"].sum()),
        "not_ready_with_fract_transform": int(not_ready["fract_transform_complete"].sum()),
        "not_ready_with_local_cryst1_candidate": int(not_ready["cryst1_recovery_candidate"].sum()),
        "completed_at": utc(),
        "runtime_seconds": time.time() - started,
    }
    atomic_json(output / "validation.json", validation)

    lines = [
        "# Filter 4 Step 0B - Crystallographic Lattice Metadata Audit",
        "",
        f"- Filter 3 HIGH + GOOD pairs: **{pair_total:,}**",
        f"- Unique PDB entries: **{len(audit):,}**",
        f"- Gemmi READY PDB: **{int(audit['gemmi_ready'].sum()):,}** ({100*float(audit['gemmi_ready'].mean()):.6f}%)",
        f"- Gemmi NOT READY PDB: **{len(not_ready):,}**",
        "",
        "## Audit classes",
        "",
        "| Audit class | n PDB | % PDB | n pairs | % pairs |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in stats:
        lines.append(f"| {item['audit_class']} | {item['n_pdb']:,} | {item['pct_pdb']:.6f} | {item['n_pairs']:,} | {item['pct_pairs']:.6f} |")
    lines += [
        "",
        "## Gemmi NOT READY recovery signals (audit only)",
        "",
        f"- Explicit symmetry operations present: {int(not_ready['explicit_symops_present'].sum()):,}",
        f"- Complete 9-element fractional matrix: {int(not_ready['fract_matrix_9_complete'].sum()):,}",
        f"- Complete matrix + vector transform: {int(not_ready['fract_transform_complete'].sum()):,}",
        f"- Local valid CRYST1 recovery candidate: {int(not_ready['cryst1_recovery_candidate'].sum()):,}",
        "",
        "No recovery, symmetry-mate construction, neighbour search, or pair exclusion was performed.",
    ]
    (output / "04_crystal_metadata_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    schema = {"schema_version": "1.0.0", "primary_key": "pdb_id", "columns": [{"name": name, "dtype": str(dtype)} for name, dtype in audit.dtypes.items()]}
    atomic_json(output / "output_schema.json", schema)
    atomic_json(output / "progress.json", {
        "status": "COMPLETED", "processed_pdb": len(audit), "total_pdb": len(audit),
        "processed_pairs": pair_total, "class_counts": dict(Counter(audit["audit_class"])),
        "runtime_seconds": time.time() - started, "updated_at": utc(),
    })
    manifest_rows = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name not in {"output_manifest.tsv", "SHA256SUMS", "runtime.log"}:
            manifest_rows.append({"relative_path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(output / "output_manifest.tsv", sep="\t", index=False)
    with (output / "SHA256SUMS").open("w", encoding="ascii") as handle:
        for path in sorted(output.iterdir()):
            if path.is_file() and path.name not in {"SHA256SUMS", "runtime.log"}:
                handle.write(f"{sha256(path)}  {path.name}\n")
    if not validation["validation_pass"]:
        raise SystemExit(2)
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
