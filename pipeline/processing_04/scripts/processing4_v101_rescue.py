#!/usr/bin/env python3
"""Processing 4 v1.0.1: non-ready-only rescue over the frozen v1.0.0 universe."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import gemmi
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import processing4_pipeline as p4  # noqa: E402
import processing4_review_audit_v2 as audit  # noqa: E402

RDLogger.DisableLog("rdApp.error")
VERSION = "processing4_v1.0.1"
BASELINE_READY = "P4_DOCKING_READY"
REVIEW = "P4_PREPARATION_REVIEW"
START_FAILED = "P4_LIGAND_START_GENERATION_FAILED"
_ORIGINAL_BUILD = p4.build_frozen_ligand


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def graph_key_local(mol: Chem.Mol, stereo_source: Chem.Mol | None = None) -> str:
    """P4-local hard gate: composition/connectivity/charge, not stereochemistry."""
    obj = Chem.RemoveHs(Chem.Mol(mol))
    obj.RemoveAllConformers()
    return Chem.MolToSmiles(obj, canonical=True, isomericSmiles=False)


def compatible_bond(a: str, b: str) -> bool:
    if a == b:
        return True
    values = {a, b}
    return "AROM" in values and values <= {"AROM", "SING", "DOUB"}


def relaxed_descriptor_mapping(descriptor: Chem.Mol, raw: Chem.Mol) -> dict[int, int] | None:
    left, right = audit.mol_graph(descriptor), audit.mol_graph(raw)
    if left.number_of_nodes() != right.number_of_nodes() or left.number_of_edges() != right.number_of_edges():
        return None
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        left, right, node_match=lambda a, b: a["element"] == b["element"]
    )
    valid: list[tuple[tuple[str, ...], dict[int, int]]] = []
    for examined, mapping in enumerate(matcher.isomorphisms_iter(), start=1):
        if examined > 4096:
            break
        if any(left.nodes[i]["charge"] != right.nodes[j]["charge"] for i, j in mapping.items()):
            continue
        ok = True
        for a, b, attr in left.edges(data=True):
            other = right.edges[mapping[a], mapping[b]]
            if not compatible_bond(attr["order"], other["order"]):
                ok = False
                break
        if not ok:
            continue
        names = tuple(raw.GetAtomWithIdx(mapping[i]).GetProp("_CCD_ATOM_ID") for i in range(descriptor.GetNumAtoms()))
        valid.append((names, dict(mapping)))
    if not valid:
        return None
    valid.sort(key=lambda item: item[0])
    return valid[0][1]


def attach_native_coordinates(mol: Chem.Mol, mapping: dict[int, int] | None,
                              raw: Chem.Mol, atom_rows: pd.DataFrame) -> tuple[Chem.Mol, list[str]]:
    coords = {
        str(r.label_atom_id): (str(r.type_symbol).upper(), float(r.Cartn_x), float(r.Cartn_y), float(r.Cartn_z))
        for r in atom_rows.itertuples() if str(r.type_symbol).upper() != "H"
    }
    obj = Chem.Mol(mol)
    names: list[str] = []
    for idx, atom in enumerate(obj.GetAtoms()):
        raw_idx = mapping[idx] if mapping is not None else idx
        raw_atom = raw.GetAtomWithIdx(raw_idx)
        if not raw_atom.HasProp("_CCD_ATOM_ID"):
            raise RuntimeError("P4 local graph atom lacks CCD atom identifier")
        name = raw_atom.GetProp("_CCD_ATOM_ID")
        if name not in coords or atom.GetSymbol().upper() != coords[name][0]:
            raise RuntimeError(f"P4 local atom/coordinate mapping failed at {name}")
        atom.SetProp("_CCD_ATOM_ID", name)
        names.append(name)
    if len(names) != len(coords) or len(set(names)) != len(names):
        raise RuntimeError("P4 local heavy-atom mapping does not close")
    conf = Chem.Conformer(obj.GetNumAtoms())
    for idx, name in enumerate(names):
        _element, x, y, z = coords[name]
        conf.SetAtomPosition(idx, (x, y, z))
    obj.RemoveAllConformers()
    obj.AddConformer(conf, assignId=True)
    return obj, names


def local_frozen_ligand(atom_rows: pd.DataFrame, bond_rows: pd.DataFrame,
                        ccd: dict[str, Any]) -> tuple[Chem.Mol, dict[str, Any]]:
    """Keep the normal graph when possible; otherwise use a P4-local usable graph without stereo veto."""
    try:
        mol, info = _ORIGINAL_BUILD(atom_rows, bond_rows, ccd)
        info["p4_local_graph_method"] = "v1.0.0_frozen_graph_without_cross_source_stereo_veto"
        return mol, info
    except Exception as original_error:
        raw, raw_info = audit.ccd_raw_molecule(atom_rows, bond_rows, ccd)
        if raw_info["missing_ccd_atom_count"] or raw_info["unsupported_bond_count"]:
            raise RuntimeError("P4_LOCAL_GRAPH_UNUSABLE: incomplete atom or bond identity") from original_error
        sanitized, full_ok, partial_ok, sanitize_error = audit.sanitize_probe(raw)
        descriptor = audit.descriptor_molecule(ccd)
        mapping = None
        method = ""
        if descriptor is not None:
            descriptor = Chem.RemoveHs(descriptor)
            mapping = relaxed_descriptor_mapping(descriptor, raw)
        if descriptor is not None and mapping is not None:
            base = Chem.Mol(descriptor)
            method = "descriptor_graph_relaxed_mapping_no_stereo_veto"
        elif full_ok and sanitized is not None:
            base = Chem.Mol(sanitized)
            mapping = None
            method = "sanitized_current_frozen_raw_graph"
        else:
            raise RuntimeError(f"P4_LOCAL_GRAPH_UNUSABLE: {sanitize_error or 'no complete graph mapping'}") from original_error
        if base.GetNumAtoms() > 1 and len(Chem.GetMolFrags(base)) != 1:
            raise RuntimeError("P4_LOCAL_GRAPH_UNUSABLE: disconnected molecular graph")
        mol, names = attach_native_coordinates(base, mapping, raw, atom_rows)
        if Chem.GetFormalCharge(mol) != int(raw_info["raw_formal_charge"]):
            raise RuntimeError("P4_LOCAL_GRAPH_UNUSABLE: formal charge changed")
        return mol, {
            "ccd_atom_order": names,
            "ccd_inchikey": ccd.get("inchikey", ""),
            "p4_local_graph_method": method,
            "p4_local_raw_graph_id": raw_info["raw_graph_id"],
            "p4_local_original_error": f"{type(original_error).__name__}: {original_error}"[:2000],
            "p4_local_raw_full_sanitize_ok": bool(full_ok),
            "p4_local_raw_partial_sanitize_ok": bool(partial_ok),
        }


def geometry_qc(mol: Chem.Mol) -> dict[str, float]:
    if mol.GetNumConformers() != 1:
        raise RuntimeError("generated start conformer count is not one")
    conf = mol.GetConformer()
    xyz = np.asarray([tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=float)
    if xyz.size == 0 or not np.isfinite(xyz).all():
        raise RuntimeError("generated start coordinates are empty or non-finite")
    minimum = math.inf
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            minimum = min(minimum, float(np.linalg.norm(xyz[i] - xyz[j])))
    if minimum < 0.35:
        raise RuntimeError(f"generated start has overlapping heavy atoms: {minimum:.6f} A")
    maximum_bond = 0.0
    minimum_bond = math.inf
    for bond in mol.GetBonds():
        length = float(np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()]))
        minimum_bond = min(minimum_bond, length)
        maximum_bond = max(maximum_bond, length)
    if minimum_bond < 0.45 or maximum_bond > 3.50:
        raise RuntimeError(f"generated start has abnormal bond geometry: {minimum_bond:.6f}-{maximum_bond:.6f} A")
    return {"minimum_heavy_atom_distance_angstrom": minimum if math.isfinite(minimum) else 0.0,
            "minimum_bond_length_angstrom": minimum_bond if math.isfinite(minimum_bond) else 0.0,
            "maximum_bond_length_angstrom": maximum_bond}


def independent_start_v101(graph_mol: Chem.Mol, config: dict[str, Any]) -> tuple[Chem.Mol, int]:
    base = Chem.Mol(graph_mol)
    base.RemoveAllConformers()
    if base.GetNumConformers() != 0:
        raise RuntimeError("coordinate removal failed")
    with_h = Chem.AddHs(base)
    primary = AllChem.ETKDGv3()
    primary.randomSeed = int(config["etkdg_random_seed"])
    primary.enforceChirality = bool(config["etkdg_enforce_chirality"])
    primary.numThreads = 1
    code = AllChem.EmbedMolecule(with_h, primary)
    fallback_used = False
    seed = int(config["etkdg_random_seed"])
    method = "ETKDGv3_UFF"
    if code != 0:
        with_h = Chem.AddHs(base)
        fallback = AllChem.ETKDGv3()
        fallback.randomSeed = int(config["etkdg_fallback_random_seed"])
        fallback.enforceChirality = bool(config["etkdg_enforce_chirality"])
        fallback.numThreads = 1
        fallback.useRandomCoords = True
        code = AllChem.EmbedMolecule(with_h, fallback)
        fallback_used = True
        seed = int(config["etkdg_fallback_random_seed"])
        method = "ETKDGv3_random_coords_UFF"
    if code != 0:
        raise RuntimeError(f"all deterministic ETKDGv3 protocols failed; final code {code}")
    uff_available = bool(AllChem.UFFHasAllMoleculeParams(with_h))
    if uff_available:
        uff_code = int(AllChem.UFFOptimizeMolecule(with_h, maxIters=int(config["uff_max_iterations"])))
    else:
        uff_code = -999
        method = "ETKDGv3_no_forcefield_fallback" if not fallback_used else "ETKDGv3_random_coords_no_forcefield_fallback"
    start = Chem.RemoveHs(with_h)
    if start.GetNumAtoms() != base.GetNumAtoms():
        raise RuntimeError("RemoveHs changed heavy atom count")
    qc = geometry_qc(start)
    start.SetProp("P4_START_CONFORMER_METHOD", method)
    start.SetProp("P4_FALLBACK_USED", json.dumps(fallback_used))
    start.SetProp("P4_UFF_AVAILABLE", json.dumps(uff_available))
    start.SetProp("P4_GENERATION_SEED", str(seed))
    start.SetProp("P4_GEOMETRY_QC", json.dumps(qc, sort_keys=True))
    return start, uff_code


def finish_ready_case(case_dir: Path, row: pd.Series) -> dict[str, Any]:
    ref = p4.load_sdf_one(case_dir / "ligand_reference.sdf")
    start = p4.load_sdf_one(case_dir / "ligand_start.sdf")
    isomeric = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(ref)), canonical=True, isomericSmiles=True)
    (case_dir / "ligand.smi").write_text(f"{isomeric}\t{row['case_id']}\n", encoding="utf-8")
    smi = Chem.MolFromSmiles(isomeric)
    if smi is None:
        raise RuntimeError("P4-local canonical isomeric SMILES parse-back failed")
    keys = {graph_key_local(ref), graph_key_local(start), graph_key_local(smi)}
    charges = {Chem.GetFormalCharge(ref), Chem.GetFormalCharge(start), Chem.GetFormalCharge(smi)}
    compositions = {tuple(sorted(a.GetSymbol() for a in mol.GetAtoms())) for mol in (ref, start, smi)}
    if len(keys) != 1 or len(charges) != 1 or len(compositions) != 1:
        raise RuntimeError("P4-local graph composition/connectivity/formal-charge consistency failed")
    qc = geometry_qc(start)
    method = start.GetProp("P4_START_CONFORMER_METHOD") if start.HasProp("P4_START_CONFORMER_METHOD") else "ETKDGv3_UFF"
    fallback_used = json.loads(start.GetProp("P4_FALLBACK_USED")) if start.HasProp("P4_FALLBACK_USED") else False
    uff_available = json.loads(start.GetProp("P4_UFF_AVAILABLE")) if start.HasProp("P4_UFF_AVAILABLE") else True
    generation_seed = int(start.GetProp("P4_GENERATION_SEED")) if start.HasProp("P4_GENERATION_SEED") else 24301
    metadata_path = case_dir / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["processing4_status"] = BASELINE_READY
    metadata["previous_processing4_status"] = str(row["previous_status"])
    metadata["original_review_reason"] = str(row["previous_reason"])
    metadata["rescue_version"] = VERSION
    metadata["rescue_applied"] = True
    metadata["stereo_cross_source_validation"] = "ignored_as_non_gating_in_v1.0.1"
    metadata["historical_warning"] = str(row.get("historical_warning", ""))
    metadata["ligand_start_generation_method"] = method
    metadata["fallback_used"] = bool(fallback_used)
    metadata["UFF_available"] = bool(uff_available)
    metadata["ETKDG_seed"] = generation_seed
    metadata["start_geometry_qc"] = qc
    metadata["native_pose_leakage_control"] = "coordinate_free_graph_only; no native-coordinate fallback"
    metadata["ligand"]["canonical_isomeric_smiles"] = isomeric
    metadata["created_at"] = utc()
    p4.atomic_json(metadata_path, metadata)
    sums = {name: sha256(case_dir / name) for name in p4.READY_FILES}
    p4.atomic_json(case_dir / "_SUCCESS.json", {
        "status": BASELINE_READY, "sha256": sums, "created_at": utc(), "rescue_version": VERSION
    })
    return {"method": method, "fallback_used": bool(fallback_used), "uff_available": bool(uff_available)}


def configure_worker() -> None:
    p4.build_frozen_ligand = local_frozen_ligand
    p4.graph_key = graph_key_local
    p4.independent_start = independent_start_v101


def process_bucket(run_dir: str, base_dir: str, bid: int) -> dict[str, Any]:
    configure_worker()
    run, base = Path(run_dir), Path(base_dir)
    inv = pq.read_table(run / "input/nonready_inventory.parquet", filters=[("bucket_id", "=", bid)]).to_pandas()
    if inv.empty:
        return {"bucket_id": bid, "cases": 0, "ready": 0, "review": 0, "start_failed": 0}
    config = p4.read_config(run / "input/config_snapshot.json")
    p2out, p3out = p4.P2_RUN / "output", p4.P3_RUN / "output"
    lig_atoms = p4.read_partition(p2out / "prepared_ligand_assembly_atoms", bid)
    lig_bonds = p4.read_partition(p2out / "prepared_ligand_assembly_bonds", bid)
    rec_atoms = p4.read_partition(p2out / "prepared_receptor_assembly_atoms", bid)
    pocket = p4.read_partition(p3out / "pair_pocket_residues", bid)
    binding = p4.read_partition(p3out / "binding_residues", bid)
    ccd = p4.CCDStore(p4.P2_RUN / "input/ccd_active_snapshot.sqlite")
    target_root = run / f"output/rescued_cases/bucket_{bid:03d}"
    target_root.mkdir(parents=True, exist_ok=True)
    results = []
    for _, row in inv.iterrows():
        cid = str(row["case_id"])
        target = target_root / cid
        marker = target / "_SUCCESS.json"
        review_marker = target / "_REVIEW.json"
        if marker.exists() or review_marker.exists():
            status_obj = json.loads((marker if marker.exists() else review_marker).read_text())
            results.append({"case_id": cid, "pair_id": str(row["pair_id"]), "bucket_id": bid,
                            "status": status_obj["status"], "reason": status_obj.get("error", "RESUMED"),
                            "previous_status": str(row["previous_status"]), "previous_reason": str(row["previous_reason"])})
            continue
        tmp = target_root / ("." + cid + f".tmp-{os.getpid()}")
        if tmp.exists():
            shutil.rmtree(tmp)
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
            extra = finish_ready_case(tmp, row)
            status, reason = BASELINE_READY, "P4_LOCAL_DOCKING_READINESS_RESCUE"
        except Exception as exc:
            text = f"{type(exc).__name__}: {exc}"[:4000]
            if isinstance(exc, p4.LigandStartFailure) or "ETKDG" in text:
                status = START_FAILED
            else:
                status = REVIEW
            reason = text
            extra = {}
            p4.atomic_json(tmp / "metadata.json", {
                "case_id": cid, "pair_id": str(row["pair_id"]), "processing4_status": status,
                "previous_processing4_status": str(row["previous_status"]),
                "original_review_reason": str(row["previous_reason"]), "rescue_version": VERSION,
                "rescue_applied": True, "error": text,
                "stereo_cross_source_validation": "ignored_as_non_gating_in_v1.0.1",
                "native_pose_leakage_control": "coordinate_free_graph_only; no native-coordinate fallback",
                "created_at": utc(),
            })
            p4.atomic_json(tmp / "_REVIEW.json", {
                "status": status, "error": text, "traceback": traceback.format_exc(limit=10),
                "rescue_version": VERSION, "created_at": utc(),
            })
        os.replace(tmp, target)
        results.append({"case_id": cid, "pair_id": str(row["pair_id"]), "bucket_id": bid,
                        "status": status, "reason": reason,
                        "previous_status": str(row["previous_status"]), "previous_reason": str(row["previous_reason"]), **extra})
    frame = pd.DataFrame(results).sort_values("case_id")
    p4.atomic_parquet(run / f"work/buckets/bucket_{bid:03d}.parquet", frame)
    counts = frame.status.value_counts()
    return {"bucket_id": bid, "cases": len(frame), "ready": int(counts.get(BASELINE_READY, 0)),
            "review": int(counts.get(REVIEW, 0)), "start_failed": int(counts.get(START_FAILED, 0))}


def prepare(args: argparse.Namespace) -> None:
    run, base = Path(args.run_dir).resolve(), Path(args.base_run).resolve()
    if run.exists() and any(run.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty run: {run}")
    if json.loads((base / "_FROZEN.json").read_text())["status"] != "FROZEN":
        raise RuntimeError("baseline is not frozen")
    for path in [run / "input", run / "output/rescued_cases", run / "work/buckets", run / "logs", run / "validation"]:
        path.mkdir(parents=True, exist_ok=True)
    inv = pd.read_parquet(base / "input/case_inventory.parquet")
    statuses = pd.read_parquet(base / "output/processing4_case_inventory.parquet")
    merged = inv.merge(statuses[["case_id", "status", "reason"]], on="case_id", validate="one_to_one")
    merged = merged.rename(columns={"status": "previous_status", "reason": "previous_reason"})
    nonready = merged[merged.previous_status.ne(BASELINE_READY)].copy()
    if len(merged) != 158226 or len(nonready) != 8705:
        raise RuntimeError(f"unexpected baseline universe: total={len(merged)} nonready={len(nonready)}")
    audit_path = base / "review_audit_v2/case_audit.parquet"
    if audit_path.exists():
        hist = pd.read_parquet(audit_path, columns=["case_id", "audit_v2_reason"]).drop_duplicates("case_id")
        nonready = nonready.merge(hist, on="case_id", how="left", validate="one_to_one")
        nonready["historical_warning"] = nonready.audit_v2_reason.fillna("")
    else:
        nonready["historical_warning"] = ""
    if args.limit:
        nonready = nonready.sort_values(["previous_status", "previous_reason", "case_id"]).groupby("previous_status", group_keys=False).head(args.limit)
    p4.atomic_parquet(run / "input/nonready_inventory.parquet", nonready)
    p4.atomic_parquet(run / "input/full_case_inventory.parquet", merged)
    config = json.loads((base / "input/config_snapshot.json").read_text())
    config.update({"schema_version": "processing4_schema_v1.0.1", "policy_version": "processing4_policy_v1.0.1",
                   "etkdg_fallback_random_seed": 24302,
                   "cross_source_stereo_validation": "warning_only",
                   "baseline_ready_reference_mode": "frozen_external_reference"})
    p4.atomic_json(run / "input/config_snapshot.json", config)
    provenance = {
        "stage": VERSION, "created_at": utc(), "baseline_run": str(base),
        "baseline_frozen_marker_sha256": sha256(base / "_FROZEN.json"),
        "baseline_output_manifest_sha256": sha256(base / "output_manifest.parquet"),
        "total_universe": 158226, "baseline_ready_immutable": 149521,
        "nonready_rescue_universe": int(len(nonready)),
        "policy": "P4-local docking usability; cross-source stereochemistry is non-gating",
        "native_pose_leakage_prohibition": "unchanged",
    }
    p4.atomic_json(run / "input/provenance.json", provenance)
    (run / "CHANGELOG.md").write_text(
        "# Processing 4 v1.0.1\n\n"
        "- Cross-source stereochemistry disagreement no longer blocks P4 readiness.\n"
        "- P4 validates only local docking representation usability.\n"
        "- Deterministic ETKDG random-coordinate fallback added (seed 24302).\n"
        "- UFF-unavailable conformers may be retained after geometry QC.\n"
        "- Native-pose leakage prohibition unchanged.\n"
        "- Existing 149,521 READY cases remain in the frozen v1.0.0 baseline and are not regenerated.\n",
        encoding="utf-8",
    )
    p4.atomic_json(run / "status.json", {"status": "PREPARED", "stage": VERSION,
                                         "nonready_cases": len(nonready), "created_at": utc()})
    print(json.dumps({"run": str(run), "nonready_cases": len(nonready)}, indent=2))


def run_rescue(args: argparse.Namespace) -> None:
    run, base = Path(args.run_dir).resolve(), Path(args.base_run).resolve()
    inv = pd.read_parquet(run / "input/nonready_inventory.parquet", columns=["bucket_id"])
    buckets = sorted(set(int(x) for x in inv.bucket_id))
    summaries = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_bucket, str(run), str(base), bid): bid for bid in buckets}
        for future in as_completed(futures):
            result = future.result()
            summaries.append(result)
            print(json.dumps(result), flush=True)
    p4.atomic_json(run / "run_summary.json", {"stage": VERSION, "finished_at": utc(),
                                               "buckets": sorted(summaries, key=lambda x: x["bucket_id"])})


def validate_rescued_case(item: tuple[Path, pd.Series]) -> str | None:
    case_dir, row = item
    try:
        success = json.loads((case_dir / "_SUCCESS.json").read_text())
        for name in p4.READY_FILES:
            path = case_dir / name
            if not path.is_file() or sha256(path) != success["sha256"].get(name):
                raise RuntimeError(f"missing/hash mismatch: {name}")
        gemmi.read_structure(str(case_dir / "receptor.cif"))
        gemmi.read_structure(str(case_dir / "receptor.pdb"))
        ref, start = p4.load_sdf_one(case_dir / "ligand_reference.sdf"), p4.load_sdf_one(case_dir / "ligand_start.sdf")
        smi = Chem.MolFromSmiles((case_dir / "ligand.smi").read_text().split("\t", 1)[0])
        if smi is None or len({graph_key_local(ref), graph_key_local(start), graph_key_local(smi)}) != 1:
            raise RuntimeError("local graph validation failed")
        geometry_qc(start)
        site = json.loads((case_dir / "site.json").read_text())
        center = list(site["site_center"].values())
        sizes = list(site["search_box"].values())
        if not all(math.isfinite(float(x)) for x in center + sizes) or not all(float(x) > 0 for x in sizes):
            raise RuntimeError("site validation failed")
        metadata = json.loads((case_dir / "metadata.json").read_text())
        if metadata.get("native_pose_leakage_control") != "coordinate_free_graph_only; no native-coordinate fallback":
            raise RuntimeError("native-pose leakage control metadata mismatch")
        return None
    except Exception as exc:
        return f"{row['case_id']}: {type(exc).__name__}: {exc}"


def validate_baseline_case(item: tuple[Path, pd.Series]) -> str | None:
    base, row = item
    case_dir = base / f"output/cases/bucket_{int(row['bucket_id']):03d}" / str(row["case_id"])
    try:
        success = json.loads((case_dir / "_SUCCESS.json").read_text())
        for name in p4.READY_FILES:
            path = case_dir / name
            if not path.is_file() or sha256(path) != success["sha256"].get(name):
                raise RuntimeError(f"baseline missing/hash mismatch: {name}")
        return None
    except Exception as exc:
        return f"{row['case_id']}: {type(exc).__name__}: {exc}"


def validate(args: argparse.Namespace) -> None:
    run, base = Path(args.run_dir).resolve(), Path(args.base_run).resolve()
    full = pd.read_parquet(run / "input/full_case_inventory.parquet")
    base_ready = full[full.previous_status.eq(BASELINE_READY)][["case_id", "pair_id", "bucket_id"]].copy()
    parts = sorted((run / "work/buckets").glob("bucket_*.parquet"))
    rescued = pd.concat([pd.read_parquet(path) for path in parts], ignore_index=True) if parts else pd.DataFrame()
    expected_nonready = pd.read_parquet(run / "input/nonready_inventory.parquet")
    errors = []
    if len(rescued) != len(expected_nonready) or not rescued.case_id.is_unique or set(rescued.case_id) != set(expected_nonready.case_id):
        errors.append("nonready status universe does not close")
    rescued_ready = rescued[rescued.status.eq(BASELINE_READY)]
    items = [(run / f"output/rescued_cases/bucket_{int(r.bucket_id):03d}" / str(r.case_id), pd.Series(r._asdict()))
             for r in rescued_ready.itertuples(index=False)]
    rescued_file_errors = []
    with ThreadPoolExecutor(max_workers=args.hash_workers) as pool:
        for value in pool.map(validate_rescued_case, items):
            if value: rescued_file_errors.append(value)
    baseline_errors = []
    if not args.skip_baseline_full_hash:
        base_items = [(base, row) for _, row in base_ready.iterrows()]
        with ThreadPoolExecutor(max_workers=args.hash_workers) as pool:
            for value in pool.map(validate_baseline_case, base_items):
                if value: baseline_errors.append(value)
    final_rows = pd.concat([
        pd.DataFrame({"case_id": base_ready.case_id, "pair_id": base_ready.pair_id,
                      "status": BASELINE_READY, "reason": "UNCHANGED_FROZEN_V1.0.0_BASELINE",
                      "case_source": "processing4_v1.0.0_frozen_baseline"}),
        rescued[["case_id", "pair_id", "status", "reason"]].assign(case_source="processing4_v1.0.1_rescue"),
    ], ignore_index=True).sort_values("case_id")
    counts = {str(k): int(v) for k, v in final_rows.status.value_counts().items()}
    closure_ok = len(final_rows) == 158226 and final_rows.case_id.is_unique and set(final_rows.case_id) == set(full.case_id)
    status_ok = sum(counts.values()) == 158226
    report = {
        "status": "PASS" if closure_ok and status_ok and not errors and not rescued_file_errors and not baseline_errors else "FAIL",
        "stage": VERSION, "expected_cases": 158226, "actual_records": len(final_rows),
        "missing": len(set(full.case_id) - set(final_rows.case_id)),
        "extra": len(set(final_rows.case_id) - set(full.case_id)),
        "duplicate": int(final_rows.case_id.duplicated().sum()), "status_counts": counts,
        "baseline_ready_checked": 0 if args.skip_baseline_full_hash else len(base_ready),
        "baseline_ready_file_hash_errors": len(baseline_errors),
        "rescued_ready_checked": len(rescued_ready), "rescued_ready_file_errors": len(rescued_file_errors),
        "error_examples": (errors + baseline_errors[:10] + rescued_file_errors[:10])[:30],
        "validated_at": utc(),
    }
    p4.atomic_parquet(run / "output/processing4_case_inventory.parquet", final_rows)
    p4.atomic_parquet(run / "output/rescue_status_transitions.parquet", rescued.sort_values("case_id"))
    p4.atomic_json(run / "validation/validation.json", report)
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(2)


def freeze(args: argparse.Namespace) -> None:
    run, base = Path(args.run_dir).resolve(), Path(args.base_run).resolve()
    validation = json.loads((run / "validation/validation.json").read_text())
    if validation.get("status") != "PASS" or validation.get("baseline_ready_checked") != 149521:
        raise RuntimeError("full validation including baseline hashes is not PASS")
    manifest = []
    for path in sorted((run / "output").rglob("*")):
        if path.is_file():
            manifest.append({"relative_path": str(path.relative_to(run)), "size_bytes": path.stat().st_size,
                             "sha256": sha256(path)})
    p4.atomic_parquet(run / "output_manifest.parquet", pd.DataFrame(manifest))
    marker = {"status": "FROZEN", "stage": VERSION, "frozen_at": utc(),
              "baseline_run": str(base), "baseline_ready_reference_count": 149521,
              "validation_sha256": sha256(run / "validation/validation.json"),
              "output_manifest_sha256": sha256(run / "output_manifest.parquet")}
    p4.atomic_json(run / "_FROZEN.json", marker)
    print(json.dumps(marker, indent=2))


def parser() -> argparse.ArgumentParser:
    default_root = Path(os.environ.get("PROCESSING4_BENCHMARK_ROOT", Path(__file__).resolve().parents[2])).resolve()
    default_base = default_root / "processing_04_docking_ready_case_construction/runs/p4_full_v1"
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--run-dir", required=True); q.add_argument("--base-run", default=str(default_base))
    q.add_argument("--limit", type=int, default=0); q.set_defaults(func=prepare)
    q = sub.add_parser("run")
    q.add_argument("--run-dir", required=True); q.add_argument("--base-run", default=str(default_base))
    q.add_argument("--workers", type=int, default=8); q.set_defaults(func=run_rescue)
    q = sub.add_parser("validate")
    q.add_argument("--run-dir", required=True); q.add_argument("--base-run", default=str(default_base))
    q.add_argument("--hash-workers", type=int, default=16); q.add_argument("--skip-baseline-full-hash", action="store_true")
    q.set_defaults(func=validate)
    q = sub.add_parser("freeze")
    q.add_argument("--run-dir", required=True); q.add_argument("--base-run", default=str(default_base)); q.set_defaults(func=freeze)
    return p


if __name__ == "__main__":
    ns = parser().parse_args()
    ns.func(ns)
