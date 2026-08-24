#!/usr/bin/env python3
"""Strict, staged rescue for representation-equivalent Processing 4 reviews.

The stage command writes only under review_rescue_v1/staged_cases.  It does not
alter baseline case directories or bucket status files.  The apply command is
intentionally separate and requires a passing staging validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import processing4_pipeline as p4  # noqa: E402
import processing4_review_audit_v2 as audit  # noqa: E402

RDLogger.DisableLog("rdApp.error")
RESCUE_VERSION = "processing4_review_rescue_v1.0.0"
ALLOWED_REASONS = {
    "AROMATIC_KEKULE_REPRESENTATION_EQUIVALENT",
    "ATOM_MAPPING_ONLY_MISMATCH",
    "CCD_DESCRIPTOR_ONLY_WARNING",
}
STAGE_READY = "RESCUE_STAGING_READY"
STAGE_FAILED = "RESCUE_STAGING_FAILED"
_MAPPING_NAME_CACHE: dict[str, tuple[list[str], int]] = {}
_START_CACHE: dict[str, tuple[Chem.Mol, int]] = {}
_ORIGINAL_INDEPENDENT_START = p4.independent_start


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def compatible_edge(a: dict[str, Any], b: dict[str, Any]) -> bool:
    if a["order"] == b["order"]:
        return True
    orders = {a["order"], b["order"]}
    return "AROM" in orders and orders <= {"AROM", "SING", "DOUB"}


def stereo_mapping_ok(descriptor: Chem.Mol, raw: Chem.Mol,
                      mapping: dict[int, int], ccd: dict[str, Any]) -> tuple[bool, str]:
    """Require mapped descriptor stereo to agree with CCD atom/bond declarations."""
    obj = Chem.Mol(descriptor)
    Chem.AssignStereochemistry(obj, cleanIt=True, force=True)
    atom_decl = {str(r[0]): str(r[4]).upper() for r in ccd["atoms"]}
    bond_decl = {frozenset((str(r[0]), str(r[1]))): str(r[4]).upper() for r in ccd["bonds"]}
    represented_atom_names = set()
    for atom in obj.GetAtoms():
        raw_atom = raw.GetAtomWithIdx(mapping[atom.GetIdx()])
        name = raw_atom.GetProp("_CCD_ATOM_ID")
        declared = atom_decl.get(name, "N")
        cip = atom.GetProp("_CIPCode").upper() if atom.HasProp("_CIPCode") else ""
        if declared in {"R", "S"}:
            represented_atom_names.add(name)
            if cip != declared:
                return False, f"CCD atom stereo {name}:{declared} != descriptor:{cip or 'UNSPECIFIED'}"
    expected_atoms = {name for name, value in atom_decl.items() if value in {"R", "S"} and
                      any(a.GetProp("_CCD_ATOM_ID") == name for a in raw.GetAtoms())}
    if represented_atom_names != expected_atoms:
        return False, "not all observed CCD atom stereo declarations were represented"

    for bond in obj.GetBonds():
        ra = raw.GetAtomWithIdx(mapping[bond.GetBeginAtomIdx()]).GetProp("_CCD_ATOM_ID")
        rb = raw.GetAtomWithIdx(mapping[bond.GetEndAtomIdx()]).GetProp("_CCD_ATOM_ID")
        declared = bond_decl.get(frozenset((ra, rb)), "N")
        stereo = bond.GetStereo().name.upper()
        normalized = "E" if stereo.endswith("E") else "Z" if stereo.endswith("Z") else "N"
        if declared in {"E", "Z"} and normalized != declared:
            return False, f"CCD bond stereo {ra}-{rb}:{declared} != descriptor:{normalized}"
    return True, ""


def strict_descriptor_mapping(descriptor: Chem.Mol, raw: Chem.Mol,
                              ccd: dict[str, Any], max_maps: int = 512) -> tuple[dict[int, int], int]:
    params = Chem.AdjustQueryParameters()
    params.makeBondsGeneric = True
    params.adjustDegree = True
    query = Chem.AdjustQueryProperties(Chem.Mol(descriptor), params)
    matches = raw.GetSubstructMatches(
        query, useChirality=False, uniquify=False, maxMatches=max_maps
    )
    valid = []
    examined = 0
    rejection_examples = []
    for match in matches:
        examined += 1
        if len(match) != descriptor.GetNumAtoms() or raw.GetNumAtoms() != descriptor.GetNumAtoms():
            continue
        mapping = {idx: int(raw_idx) for idx, raw_idx in enumerate(match)}
        charge_ok = all(
            descriptor.GetAtomWithIdx(i).GetFormalCharge() == raw.GetAtomWithIdx(j).GetFormalCharge()
            for i, j in mapping.items()
        )
        bond_ok = all(
            compatible_edge(
                {"order": audit.bond_label(bond)},
                {"order": audit.bond_label(raw.GetBondBetweenAtoms(
                    mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]))},
            )
            for bond in descriptor.GetBonds()
        )
        if not charge_ok or not bond_ok:
            continue
        ok, reason = stereo_mapping_ok(descriptor, raw, mapping, ccd)
        if ok:
            names = tuple(raw.GetAtomWithIdx(mapping[i]).GetProp("_CCD_ATOM_ID")
                          for i in range(descriptor.GetNumAtoms()))
            valid.append((names, dict(mapping)))
        elif len(rejection_examples) < 3:
            rejection_examples.append(reason)
        if examined >= max_maps:
            break
    if not valid:
        detail = "; ".join(rejection_examples) or "no charge/bond-compatible graph isomorphism"
        raise RuntimeError(f"no stereo-preserving descriptor mapping after {examined} maps: {detail}")
    valid.sort(key=lambda item: item[0])
    return valid[0][1], examined


def rescued_frozen_ligand(atom_rows: pd.DataFrame, bond_rows: pd.DataFrame,
                          ccd: dict[str, Any]) -> tuple[Chem.Mol, dict[str, Any]]:
    raw, raw_info = audit.ccd_raw_molecule(atom_rows, bond_rows, ccd)
    sanitized, full_ok, partial_ok, sanitize_error = audit.sanitize_probe(raw)
    if not partial_ok or sanitized is None:
        raise RuntimeError(f"authoritative raw graph invalid: {sanitize_error}")
    descriptor = audit.descriptor_molecule(ccd)
    coord = {
        str(r.label_atom_id): (str(r.type_symbol).upper(), float(r.Cartn_x),
                               float(r.Cartn_y), float(r.Cartn_z))
        for r in atom_rows.itertuples() if str(r.type_symbol).upper() != "H"
    }
    if descriptor is None:
        if audit.declared_stereo_count(ccd) != 0 or not full_ok:
            raise RuntimeError("descriptor missing and authoritative raw graph is not unambiguous")
        mol = Chem.Mol(sanitized)
        mapping_examined = 0
        method = "authoritative_raw_graph_descriptor_warning"
    else:
        descriptor = Chem.RemoveHs(descriptor)
        cache_key = raw_info["raw_graph_id"] + ":" + hashlib.sha256(
            str(ccd.get("descriptor", "")).encode()
        ).hexdigest()[:24]
        if cache_key in _MAPPING_NAME_CACHE:
            mapped_names, mapping_examined = _MAPPING_NAME_CACHE[cache_key]
            raw_by_name = {a.GetProp("_CCD_ATOM_ID"): a.GetIdx() for a in raw.GetAtoms()}
            mapping = {idx: raw_by_name[name] for idx, name in enumerate(mapped_names)}
        else:
            mapping, mapping_examined = strict_descriptor_mapping(descriptor, raw, ccd)
            mapped_names = [raw.GetAtomWithIdx(mapping[idx]).GetProp("_CCD_ATOM_ID")
                            for idx in range(descriptor.GetNumAtoms())]
            _MAPPING_NAME_CACHE[cache_key] = (mapped_names, mapping_examined)
        mol = Chem.Mol(descriptor)
        for didx, ridx in mapping.items():
            name = raw.GetAtomWithIdx(ridx).GetProp("_CCD_ATOM_ID")
            mol.GetAtomWithIdx(didx).SetProp("_CCD_ATOM_ID", name)
        method = "descriptor_representation_mapped_to_authoritative_raw_graph"
    conf = Chem.Conformer(mol.GetNumAtoms())
    atom_order = []
    for idx, atom in enumerate(mol.GetAtoms()):
        if not atom.HasProp("_CCD_ATOM_ID"):
            raise RuntimeError(f"mapped atom {idx} lacks CCD atom id")
        name = atom.GetProp("_CCD_ATOM_ID")
        if name not in coord:
            raise RuntimeError(f"mapped CCD atom {name} lacks frozen coordinates")
        element, x, y, z = coord[name]
        if atom.GetSymbol().upper() != element:
            raise RuntimeError(f"mapped element mismatch at {name}")
        conf.SetAtomPosition(idx, (x, y, z))
        atom_order.append(name)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if mol.GetNumHeavyAtoms() != len(coord):
        raise RuntimeError("rescue heavy atom mapping does not close")
    if Chem.GetFormalCharge(mol) != raw_info["raw_formal_charge"]:
        raise RuntimeError("rescue formal charge changed relative to authoritative raw graph")
    return mol, {
        "ccd_atom_order": atom_order,
        "ccd_inchikey": ccd.get("inchikey", ""),
        "rescue_graph_construction_method": method,
        "rescue_mapping_count_examined": mapping_examined,
        "authoritative_raw_graph_id": raw_info["raw_graph_id"],
    }


def cached_independent_start(graph_mol: Chem.Mol, config: dict[str, Any]) -> tuple[Chem.Mol, int]:
    key = "|".join([
        audit.canonical(graph_mol, True), str(config["etkdg_random_seed"]),
        str(config["etkdg_enforce_chirality"]), str(config["uff_max_iterations"]),
    ])
    if key not in _START_CACHE:
        mol, code = _ORIGINAL_INDEPENDENT_START(graph_mol, config)
        _START_CACHE[key] = (Chem.Mol(mol), int(code))
    mol, code = _START_CACHE[key]
    return Chem.Mol(mol), int(code)


def annotate_success(case_dir: Path, audit_row: pd.Series) -> None:
    metadata_path = case_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["review_rescue"] = {
        "previous_status": str(audit_row["previous_status"]),
        "previous_reason": str(audit_row["previous_reason"]),
        "audit_v2_reason": str(audit_row["audit_v2_reason"]),
        "rescue_rule_applied": True,
        "rescue_rule_version": RESCUE_VERSION,
        "new_status": "P4_DOCKING_READY",
        "new_reason": "STRICT_REPRESENTATION_EQUIVALENCE_RESCUE",
    }
    metadata["created_at"] = utc()
    p4.atomic_json(metadata_path, metadata)
    sums = {name: p4.sha256_file(case_dir / name) for name in p4.READY_FILES}
    p4.atomic_json(case_dir / "_SUCCESS.json", {
        "status": "P4_DOCKING_READY", "sha256": sums, "created_at": utc(),
        "rescue_rule_version": RESCUE_VERSION,
    })


def stage_bucket(run_dir: str, bucket_id: int) -> list[dict[str, Any]]:
    run = Path(run_dir)
    audit_frame = pq.read_table(
        run / "review_audit_v2/case_audit.parquet", filters=[("bucket_id", "=", bucket_id)]
    ).to_pandas()
    audit_frame = audit_frame[audit_frame["audit_v2_reason"].isin(ALLOWED_REASONS)]
    if audit_frame.empty:
        return []
    inventory = pq.read_table(
        run / "input/case_inventory.parquet", filters=[("bucket_id", "=", bucket_id)]
    ).to_pandas()
    rows = audit_frame.merge(inventory, on=["case_id", "pair_id", "bucket_id", "component_id",
                                             "filter3_quality_class"], how="left", validate="one_to_one")
    p2out, p3out = p4.P2_RUN / "output", p4.P3_RUN / "output"
    lig_atoms = p4.read_partition(p2out / "prepared_ligand_assembly_atoms", bucket_id)
    lig_bonds = p4.read_partition(p2out / "prepared_ligand_assembly_bonds", bucket_id)
    rec_atoms = p4.read_partition(p2out / "prepared_receptor_assembly_atoms", bucket_id)
    pocket = p4.read_partition(p3out / "pair_pocket_residues", bucket_id)
    binding = p4.read_partition(p3out / "binding_residues", bucket_id)
    ccd = p4.CCDStore(p4.P2_RUN / "input/ccd_active_snapshot.sqlite")
    config = p4.read_config(run / "input/config_snapshot.json")
    stage_root = run / f"review_rescue_v1/staged_cases/bucket_{bucket_id:03d}"
    stage_root.mkdir(parents=True, exist_ok=True)
    p4.build_frozen_ligand = rescued_frozen_ligand
    p4.independent_start = cached_independent_start
    results = []
    for _, row in rows.iterrows():
        cid = str(row["case_id"])
        target = stage_root / cid
        tmp = stage_root / ("." + cid + f".tmp-{os.getpid()}")
        if tmp.exists():
            shutil.rmtree(tmp)
        if target.exists():
            shutil.rmtree(target)
        tmp.mkdir()
        placement = str(row["ligand_assembly_placement_id"])
        chains = [x for x in str(row["receptor_chain_instance_ids"]).split(",") if x]
        try:
            p4.create_case(
                tmp, row,
                rec_atoms[rec_atoms["filter_1_chain_instance_id"].astype(str).isin(chains)],
                lig_atoms[lig_atoms["filter_2_ligand_assembly_placement_id"].astype(str).eq(placement)],
                lig_bonds[lig_bonds["filter_2_ligand_assembly_placement_id"].astype(str).eq(placement)],
                pocket[pocket["pair_id"].astype(str).eq(str(row["pair_id"]))] if not pocket.empty else pocket,
                binding[binding["ligand_assembly_placement_id"].astype(str).eq(placement)] if not binding.empty else binding,
                ccd, config,
            )
            annotate_success(tmp, row)
            status, error = STAGE_READY, ""
        except Exception as exc:
            status, error = STAGE_FAILED, f"{type(exc).__name__}: {exc}"[:4000]
            p4.atomic_json(tmp / "_RESCUE_FAILED.json", {
                "status": status, "error": error,
                "traceback": traceback.format_exc(limit=10), "created_at": utc(),
            })
        os.replace(tmp, target)
        results.append({
            "case_id": cid, "pair_id": str(row["pair_id"]), "bucket_id": bucket_id,
            "component_id": str(row["component_id"]),
            "filter3_quality_class": str(row["filter3_quality_class"]),
            "previous_status": str(row["previous_status"]),
            "previous_reason": str(row["previous_reason"]),
            "audit_v2_reason": str(row["audit_v2_reason"]),
            "rescue_rule_version": RESCUE_VERSION,
            "staging_status": status, "staging_reason": error,
        })
    return results


def validate_staging(run: Path, frame: pd.DataFrame) -> dict[str, Any]:
    errors = []
    for row in frame[frame["staging_status"].eq(STAGE_READY)].itertuples():
        case_dir = run / f"review_rescue_v1/staged_cases/bucket_{int(row.bucket_id):03d}" / row.case_id
        success_path = case_dir / "_SUCCESS.json"
        try:
            success = json.loads(success_path.read_text())
            for name in p4.READY_FILES:
                path = case_dir / name
                if not path.is_file() or p4.sha256_file(path) != success["sha256"].get(name):
                    raise RuntimeError(f"missing/hash mismatch {name}")
        except Exception as exc:
            errors.append({"case_id": row.case_id, "error": str(exc)})
    expected = int(frame.shape[0])
    report = {
        "status": "PASS" if expected == 1450 and frame["case_id"].is_unique and not errors else "FAIL",
        "rescue_rule_version": RESCUE_VERSION,
        "expected_candidates": 1450,
        "staged_candidates": expected,
        "unique_candidates": int(frame["case_id"].nunique()),
        "staging_status_counts": {str(k): int(v) for k, v in frame["staging_status"].value_counts().items()},
        "ready_by_audit_reason": {
            str(k): int(v) for k, v in frame[frame["staging_status"].eq(STAGE_READY)]["audit_v2_reason"].value_counts().items()
        },
        "failed_by_audit_reason": {
            str(k): int(v) for k, v in frame[frame["staging_status"].eq(STAGE_FAILED)]["audit_v2_reason"].value_counts().items()
        },
        "ready_file_error_count": len(errors),
        "ready_file_error_examples": errors[:20],
        "validated_at": utc(),
    }
    p4.atomic_json(run / "review_rescue_v1/staging_validation.json", report)
    return report


def stage_command(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    if (run / "_FROZEN.json").exists():
        raise RuntimeError("refusing to stage against a frozen run")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(stage_bucket, str(run), bid): bid for bid in range(256)}
        for future in as_completed(futures):
            bid = futures[future]
            rows = future.result()
            results.extend(rows)
            print(json.dumps({"bucket_id": bid, "candidates": len(rows),
                              "ready": sum(x["staging_status"] == STAGE_READY for x in rows)}), flush=True)
    frame = pd.DataFrame(results).sort_values("case_id")
    out = run / "review_rescue_v1"
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), out / "staging_results.parquet", compression="zstd")
    frame.to_csv(out / "staging_results.tsv.gz", sep="\t", index=False, compression="gzip")
    report = validate_staging(run, frame)
    print(json.dumps(report, indent=2), flush=True)
    if report["status"] != "PASS":
        raise SystemExit(2)


def verify_staged_success(case_dir: Path) -> dict[str, Any]:
    success = json.loads((case_dir / "_SUCCESS.json").read_text(encoding="utf-8"))
    if success.get("status") != "P4_DOCKING_READY" or success.get("rescue_rule_version") != RESCUE_VERSION:
        raise RuntimeError("staged success marker status/version mismatch")
    for name in p4.READY_FILES:
        path = case_dir / name
        if not path.is_file() or p4.sha256_file(path) != success.get("sha256", {}).get(name):
            raise RuntimeError(f"staged hash mismatch: {name}")
    return success


def ensure_pre_apply_snapshot(run: Path) -> Path:
    snapshot = run / "review_rescue_v1/pre_apply_status_buckets"
    snapshot.mkdir(parents=True, exist_ok=True)
    source_files = sorted((run / "work/buckets").glob("bucket_*.parquet"))
    if len(source_files) != 256:
        raise RuntimeError(f"expected 256 baseline status buckets, found {len(source_files)}")
    manifest = []
    for source in source_files:
        target = snapshot / source.name
        if not target.exists():
            shutil.copy2(source, target)
        manifest.append({"file": source.name, "sha256": p4.sha256_file(target),
                         "size_bytes": target.stat().st_size})
    p4.atomic_json(snapshot / "snapshot_manifest.json", {
        "created_at": utc(), "files": manifest,
        "baseline_validation_sha256": p4.sha256_file(run / "validation/validation.json"),
    })
    return snapshot


def apply_command(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    if (run / "_FROZEN.json").exists():
        raise RuntimeError("refusing to apply rescue to a frozen run")
    validation = json.loads((run / "review_rescue_v1/staging_validation.json").read_text())
    if validation.get("status") != "PASS" or validation.get("staging_status_counts", {}).get(STAGE_READY) != 1181:
        raise RuntimeError("staging validation is not the expected PASS with 1181 ready cases")
    ensure_pre_apply_snapshot(run)
    state_path = run / "review_rescue_v1/apply_state.json"
    if state_path.exists():
        state = json.loads(state_path.read_text())
        apply_id = state["apply_id"]
    else:
        apply_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        state = {"apply_id": apply_id, "started_at": utc(), "status": "IN_PROGRESS",
                 "rescue_rule_version": RESCUE_VERSION}
        p4.atomic_json(state_path, state)

    staged = pq.read_table(run / "review_rescue_v1/staging_results.parquet").to_pandas()
    apply_bucket_root = run / "review_rescue_v1/apply_buckets"
    apply_bucket_root.mkdir(parents=True, exist_ok=True)
    for bid in range(256):
        subset = staged[staged["bucket_id"].eq(bid)].copy()
        if subset.empty:
            continue
        status_path = run / f"work/buckets/bucket_{bid:03d}.parquet"
        statuses = pq.read_table(status_path).to_pandas()
        transitions = []
        for row in subset.itertuples():
            mask = statuses["case_id"].astype(str).eq(str(row.case_id))
            if int(mask.sum()) != 1:
                raise RuntimeError(f"status row closure failed for {row.case_id}")
            current = statuses.loc[mask].iloc[0]
            transition = {
                "case_id": str(row.case_id), "pair_id": str(row.pair_id), "bucket_id": bid,
                "component_id": str(row.component_id),
                "filter3_quality_class": str(row.filter3_quality_class),
                "previous_status": str(row.previous_status),
                "previous_reason": str(row.previous_reason),
                "audit_v2_reason": str(row.audit_v2_reason),
                "rescue_rule_version": RESCUE_VERSION,
                "staging_status": str(row.staging_status),
                "staging_reason": str(row.staging_reason),
            }
            if row.staging_status != STAGE_READY:
                transition.update({
                    "rescue_rule_applied": False,
                    "new_status": str(current["status"]),
                    "new_reason": str(current["reason"]),
                    "application_status": "NOT_APPLIED_STAGING_FAILED",
                })
                transitions.append(transition)
                continue
            target = run / f"output/cases/bucket_{bid:03d}" / str(row.case_id)
            stage_dir = run / f"review_rescue_v1/staged_cases/bucket_{bid:03d}" / str(row.case_id)
            already_applied = False
            if str(current["status"]) == "P4_DOCKING_READY" and target.exists():
                success = verify_staged_success(target)
                already_applied = success.get("rescue_rule_version") == RESCUE_VERSION
            if not already_applied:
                if str(current["status"]) != "P4_PREPARATION_REVIEW":
                    raise RuntimeError(f"unexpected pre-apply status for {row.case_id}: {current['status']}")
                if str(current["reason"]) != str(row.previous_reason):
                    raise RuntimeError(f"previous reason changed for {row.case_id}")
                verify_staged_success(stage_dir)
                backup_root = run / f"work/superseded_nonready/review_rescue_v1/bucket_{bid:03d}"
                backup_root.mkdir(parents=True, exist_ok=True)
                backup = backup_root / (str(row.case_id) + f".superseded-{apply_id}")
                if target.exists() and not backup.exists():
                    os.replace(target, backup)
                if stage_dir.exists():
                    os.replace(stage_dir, target)
                verify_staged_success(target)
            metadata = json.loads((target / "metadata.json").read_text())
            statuses.loc[mask, "status"] = "P4_DOCKING_READY"
            statuses.loc[mask, "reason"] = "STRICT_REPRESENTATION_EQUIVALENCE_RESCUE|" + str(row.audit_v2_reason)
            statuses.loc[mask, "receptor_atoms"] = metadata.get("receptor_heavy_atom_count")
            statuses.loc[mask, "ligand_heavy_atoms"] = metadata.get("ligand", {}).get("heavy_atom_count")
            transition.update({
                "rescue_rule_applied": True,
                "new_status": "P4_DOCKING_READY",
                "new_reason": "STRICT_REPRESENTATION_EQUIVALENCE_RESCUE|" + str(row.audit_v2_reason),
                "application_status": "ALREADY_APPLIED" if already_applied else "APPLIED",
            })
            transitions.append(transition)
        p4.atomic_parquet(status_path, statuses)
        p4.atomic_parquet(apply_bucket_root / f"bucket_{bid:03d}.parquet", pd.DataFrame(transitions))
        print(json.dumps({"bucket_id": bid, "candidates": len(subset),
                          "applied": sum(x["rescue_rule_applied"] for x in transitions)}), flush=True)

    parts = sorted(apply_bucket_root.glob("bucket_*.parquet"))
    transition_frame = pd.concat([pq.read_table(path).to_pandas() for path in parts], ignore_index=True)
    transition_frame = transition_frame.sort_values("case_id")
    if len(transition_frame) != 1450 or not transition_frame["case_id"].is_unique:
        raise RuntimeError("transition provenance universe does not close")
    p4.atomic_parquet(run / "review_rescue_v1/status_transitions.parquet", transition_frame)
    transition_frame.to_csv(run / "review_rescue_v1/status_transitions.tsv.gz", sep="\t", index=False, compression="gzip")
    applied = int(transition_frame["rescue_rule_applied"].sum())
    if applied != 1181:
        raise RuntimeError(f"expected 1181 applied rescues, got {applied}")
    state.update({"status": "APPLIED", "finished_at": utc(), "applied_rescues": applied,
                  "not_applied": int(len(transition_frame) - applied)})
    p4.atomic_json(state_path, state)
    p4.atomic_json(run / "review_rescue_v1/apply_summary.json", state)
    print(json.dumps(state, indent=2), flush=True)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("stage")
    q.add_argument("--run-dir", required=True)
    q.add_argument("--workers", type=int, default=8)
    q.set_defaults(func=stage_command)
    q = sub.add_parser("apply")
    q.add_argument("--run-dir", required=True)
    q.set_defaults(func=apply_command)
    return p


if __name__ == "__main__":
    ns = parser().parse_args()
    ns.func(ns)
