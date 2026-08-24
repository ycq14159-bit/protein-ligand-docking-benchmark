#!/usr/bin/env python3
"""Read-only reason audit for Processing 4 PREPARATION_REVIEW cases.

This script never edits case directories or bucket status files.  It classifies
the existing review universe and emits evidence/provenance under a separate
review_audit_v2 directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import rdkit
from rdkit import Chem, RDLogger

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import processing4_pipeline as p4  # noqa: E402

RDLogger.DisableLog("rdApp.error")
AUDIT_VERSION = "processing4_review_reason_audit_v2.0.0"
REVIEW_STATUS = "P4_PREPARATION_REVIEW"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def original_category(reason: str) -> str:
    if "ligand graph identity mismatch" in reason:
        return "GRAPH_IDENTITY_OR_STEREO_MISMATCH"
    if "KekulizeException" in reason or "Can't kekulize" in reason:
        return "KEKULIZATION_FAILED"
    if "CCD descriptor graph cannot be mapped" in reason:
        return "CCD_DESCRIPTOR_MAPPING_FAILED"
    if "canonical descriptor SMILES missing or unparseable" in reason:
        return "CCD_DESCRIPTOR_MISSING_OR_BAD"
    if "frozen ligand atoms/bonds missing" in reason:
        return "UPSTREAM_FROZEN_GRAPH_MISSING"
    return "OTHER_PREPARATION_REVIEW"


def bond_label(bond: Chem.Bond) -> str:
    if bond.GetIsAromatic() or bond.GetBondType() == Chem.BondType.AROMATIC:
        return "AROM"
    value = float(bond.GetBondTypeAsDouble())
    if value == 1.0:
        return "SING"
    if value == 2.0:
        return "DOUB"
    if value == 3.0:
        return "TRIP"
    return f"OTHER:{value:g}"


def mol_graph(mol: Chem.Mol) -> nx.Graph:
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(
            atom.GetIdx(),
            element=atom.GetSymbol().upper(),
            charge=int(atom.GetFormalCharge()),
            aromatic=bool(atom.GetIsAromatic()),
        )
    for bond in mol.GetBonds():
        graph.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
            order=bond_label(bond), aromatic=bool(bond.GetIsAromatic()),
        )
    return graph


def graph_comparison(left: Chem.Mol, right: Chem.Mol, max_maps: int = 2048) -> dict[str, Any]:
    """Compare element/connectivity, charge, and bond order across graph isomorphisms."""
    lg, rg = mol_graph(left), mol_graph(right)
    if lg.number_of_nodes() != rg.number_of_nodes() or lg.number_of_edges() != rg.number_of_edges():
        return {"topology_equal": False, "mapping_count_examined": 0}
    if Counter(nx.get_node_attributes(lg, "element").values()) != Counter(nx.get_node_attributes(rg, "element").values()):
        return {"topology_equal": False, "mapping_count_examined": 0}
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        lg, rg, node_match=lambda a, b: a["element"] == b["element"]
    )
    best = None
    examined = 0
    for mapping in matcher.isomorphisms_iter():
        examined += 1
        charge_mismatch = 0
        exact_bond_mismatch = 0
        aromatic_compatible_mismatch = 0
        nonaromatic_bond_mismatch = 0
        for li, ri in mapping.items():
            if lg.nodes[li]["charge"] != rg.nodes[ri]["charge"]:
                charge_mismatch += 1
        for la, lb, attr in lg.edges(data=True):
            ra, rb = mapping[la], mapping[lb]
            other = rg.edges[ra, rb]
            if attr["order"] == other["order"]:
                continue
            exact_bond_mismatch += 1
            orders = {attr["order"], other["order"]}
            if "AROM" in orders and orders <= {"AROM", "SING", "DOUB"}:
                aromatic_compatible_mismatch += 1
            else:
                nonaromatic_bond_mismatch += 1
        score = (charge_mismatch, nonaromatic_bond_mismatch,
                 exact_bond_mismatch - aromatic_compatible_mismatch, exact_bond_mismatch)
        candidate = {
            "topology_equal": True,
            "charge_mismatch_atoms": charge_mismatch,
            "bond_order_mismatch_edges": exact_bond_mismatch,
            "aromatic_compatible_mismatch_edges": aromatic_compatible_mismatch,
            "nonaromatic_bond_mismatch_edges": nonaromatic_bond_mismatch,
            "mapping_count_examined": examined,
        }
        if best is None or score < best[0]:
            best = (score, candidate)
        if score == (0, 0, 0, 0) or examined >= max_maps:
            break
    if best is None:
        return {"topology_equal": False, "mapping_count_examined": examined}
    best[1]["mapping_count_examined"] = examined
    return best[1]


def ccd_raw_molecule(atom_rows: pd.DataFrame, bond_rows: pd.DataFrame,
                     ccd: dict[str, Any]) -> tuple[Chem.Mol, dict[str, Any]]:
    observed = {
        str(r.label_atom_id): str(r.type_symbol).upper()
        for r in atom_rows.itertuples() if str(r.type_symbol).upper() != "H"
    }
    ccd_atom_by_name = {str(r[0]): r for r in ccd["atoms"] if str(r[1]).upper() != "H"}
    missing_ccd_atoms = sorted(set(observed) - set(ccd_atom_by_name))
    rw = Chem.RWMol()
    index: dict[str, int] = {}
    for name in sorted(observed):
        if name not in ccd_atom_by_name:
            continue
        _name, element, charge, aromatic, _stereo, _ordinal = ccd_atom_by_name[name]
        atom = Chem.Atom(str(element).title())
        atom.SetFormalCharge(int(charge))
        atom.SetBoolProp("_P4_RAW_AROMATIC", str(aromatic).upper() == "Y")
        atom.SetProp("_CCD_ATOM_ID", name)
        if str(aromatic).upper() == "Y":
            atom.SetIsAromatic(True)
        index[name] = rw.AddAtom(atom)
    unsupported_bonds = []
    seen = set()
    for row in bond_rows.sort_values("bond_index").itertuples():
        a, b = str(row.atom_id_1), str(row.atom_id_2)
        if a not in index or b not in index:
            continue
        edge = tuple(sorted((a, b)))
        if edge in seen:
            continue
        seen.add(edge)
        order = str(row.bond_order).upper()
        aromatic = str(row.aromatic_flag).upper() == "Y"
        types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                 "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
        bt = Chem.BondType.AROMATIC if aromatic else types.get(order)
        if bt is None:
            unsupported_bonds.append(f"{a}-{b}:{order}")
            continue
        rw.AddBond(index[a], index[b], bt)
        if aromatic:
            rw.GetAtomWithIdx(index[a]).SetIsAromatic(True)
            rw.GetAtomWithIdx(index[b]).SetIsAromatic(True)
    mol = rw.GetMol()
    signature = {
        "atoms": sorted((n, observed[n], int(ccd_atom_by_name[n][2]), str(ccd_atom_by_name[n][3]))
                        for n in observed if n in ccd_atom_by_name),
        "bonds": sorted((str(r.atom_id_1), str(r.atom_id_2), str(r.bond_order), str(r.aromatic_flag))
                        for r in bond_rows.itertuples()
                        if str(r.atom_id_1) in observed and str(r.atom_id_2) in observed),
    }
    graph_id = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:24]
    info = {
        "raw_graph_id": graph_id,
        "observed_heavy_atoms": len(observed),
        "raw_bond_count": mol.GetNumBonds(),
        "raw_aromatic_atom_count": sum(a.GetIsAromatic() for a in mol.GetAtoms()),
        "raw_aromatic_bond_count": sum(b.GetIsAromatic() for b in mol.GetBonds()),
        "raw_formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "missing_ccd_atom_count": len(missing_ccd_atoms),
        "missing_ccd_atom_ids": ",".join(missing_ccd_atoms[:30]),
        "unsupported_bond_count": len(unsupported_bonds),
    }
    return mol, info


def sanitize_probe(raw: Chem.Mol) -> tuple[Chem.Mol | None, bool, bool, str]:
    full = Chem.Mol(raw)
    try:
        Chem.SanitizeMol(full)
        return full, True, True, ""
    except Exception as exc:
        full_error = f"{type(exc).__name__}: {exc}"
    partial = Chem.Mol(raw)
    try:
        ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        Chem.SanitizeMol(partial, sanitizeOps=ops)
        return partial, False, True, full_error
    except Exception as exc:
        return None, False, False, full_error + " | no-kekule: " + f"{type(exc).__name__}: {exc}"


def canonical(mol: Chem.Mol, isomeric: bool) -> str:
    obj = Chem.RemoveHs(Chem.Mol(mol))
    obj.RemoveAllConformers()
    return Chem.MolToSmiles(obj, canonical=True, isomericSmiles=isomeric)


def stereo_profile(mol: Chem.Mol) -> dict[str, Any]:
    obj = Chem.Mol(mol)
    Chem.AssignStereochemistry(obj, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(obj, includeUnassigned=True, includeCIP=True)
    bonds = [b.GetStereo().name for b in obj.GetBonds()
             if b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)]
    return {
        "assigned_atom_stereo": sum(label != "?" for _, label in centers),
        "unassigned_atom_stereo": sum(label == "?" for _, label in centers),
        "assigned_bond_stereo": len(bonds),
        "cip_labels": ",".join(sorted(label for _, label in centers)),
        "bond_stereo_labels": ",".join(sorted(bonds)),
    }


def load_case_molecules(case_dir: Path) -> tuple[Chem.Mol | None, Chem.Mol | None, Chem.Mol | None]:
    try:
        smi = (case_dir / "ligand.smi").read_text(encoding="utf-8").split("\t", 1)[0].strip()
        m_smi = Chem.MolFromSmiles(smi)
    except Exception:
        m_smi = None
    def sdf(name: str) -> Chem.Mol | None:
        path = case_dir / name
        if not path.exists():
            return None
        try:
            return next((m for m in Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True) if m is not None), None)
        except Exception:
            return None
    return m_smi, sdf("ligand_reference.sdf"), sdf("ligand_start.sdf")


def classify_existing_graph_mismatch(case_dir: Path) -> tuple[str, dict[str, Any]]:
    smi, ref, start = load_case_molecules(case_dir)
    evidence: dict[str, Any] = {}
    if smi is None or ref is None or start is None:
        return "OTHER_GRAPH_VALIDATION_REVIEW", {"file_parse_complete": False}
    mols = {"smi": smi, "reference": ref, "start": start}
    evidence["file_parse_complete"] = True
    for name, mol in mols.items():
        evidence[f"{name}_nonisomeric_smiles"] = canonical(mol, False)
        evidence[f"{name}_isomeric_smiles"] = canonical(mol, True)
        evidence[f"{name}_formal_charge"] = int(Chem.GetFormalCharge(mol))
        for key, value in stereo_profile(mol).items():
            evidence[f"{name}_{key}"] = value
    noniso = {evidence[f"{name}_nonisomeric_smiles"] for name in mols}
    iso = {evidence[f"{name}_isomeric_smiles"] for name in mols}
    charges = {evidence[f"{name}_formal_charge"] for name in mols}
    if len(charges) != 1:
        return "FORMAL_CHARGE_MISMATCH", evidence
    if len(noniso) != 1:
        comparisons = [graph_comparison(smi, ref), graph_comparison(smi, start)]
        evidence["pair_graph_comparisons"] = json.dumps(comparisons, sort_keys=True)
        if not all(x.get("topology_equal") for x in comparisons):
            return "GRAPH_TOPOLOGY_TRUE_MISMATCH", evidence
        if any(x.get("charge_mismatch_atoms", 0) for x in comparisons):
            return "FORMAL_CHARGE_MISMATCH", evidence
        if all(x.get("nonaromatic_bond_mismatch_edges", 0) == 0 and
               x.get("bond_order_mismatch_edges", 0) == x.get("aromatic_compatible_mismatch_edges", 0)
               for x in comparisons):
            return "AROMATIC_KEKULE_REPRESENTATION_EQUIVALENT", evidence
        return "BOND_ORDER_TRUE_MISMATCH", evidence
    if len(iso) == 1:
        return "STEREO_REPRESENTATION_EQUIVALENT", evidence
    profiles = [stereo_profile(m) for m in mols.values()]
    assigned_atoms = {x["assigned_atom_stereo"] for x in profiles}
    assigned_bonds = {x["assigned_bond_stereo"] for x in profiles}
    unassigned = {x["unassigned_atom_stereo"] for x in profiles}
    if len(assigned_atoms) > 1 or len(assigned_bonds) > 1 or max(unassigned) > min(unassigned):
        return "STEREO_ONE_SIDE_UNSPECIFIED", evidence
    return "STEREO_TRUE_CONFLICT", evidence


def descriptor_molecule(ccd: dict[str, Any]) -> Chem.Mol | None:
    if not ccd.get("descriptor"):
        return None
    try:
        mol = Chem.MolFromSmiles(ccd["descriptor"])
        return Chem.RemoveHs(mol) if mol is not None else None
    except Exception:
        return None


def declared_stereo_count(ccd: dict[str, Any]) -> int:
    values = [str(r[4]).upper() for r in ccd["atoms"]] + [str(r[4]).upper() for r in ccd["bonds"]]
    return sum(v not in {"", "N", "NONE", "?"} for v in values)


def classify_raw_review(category: str, raw: Chem.Mol, sanitized: Chem.Mol | None,
                        full_ok: bool, partial_ok: bool, descriptor: Chem.Mol | None,
                        ccd: dict[str, Any], config: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    evidence: dict[str, Any] = {
        "raw_full_sanitize_ok": full_ok,
        "raw_no_kekule_sanitize_ok": partial_ok,
        "ccd_descriptor_parse_ok": descriptor is not None,
        "ccd_declared_stereo_count": declared_stereo_count(ccd),
    }
    if category == "UPSTREAM_FROZEN_GRAPH_MISSING":
        return "UPSTREAM_FROZEN_GRAPH_INCOMPLETE", evidence
    if not partial_ok or sanitized is None:
        return "TRUE_INVALID_FROZEN_GRAPH", evidence
    if descriptor is not None:
        comparison = graph_comparison(raw, descriptor)
        evidence.update({f"descriptor_{k}": v for k, v in comparison.items()})
        if not comparison.get("topology_equal"):
            return "GRAPH_TOPOLOGY_TRUE_MISMATCH", evidence
        if comparison.get("charge_mismatch_atoms", 0):
            return "FORMAL_CHARGE_MISMATCH", evidence
        bond_bad = comparison.get("nonaromatic_bond_mismatch_edges", 0)
        bond_all = comparison.get("bond_order_mismatch_edges", 0)
        bond_arom = comparison.get("aromatic_compatible_mismatch_edges", 0)
        if bond_bad:
            return "BOND_ORDER_TRUE_MISMATCH", evidence
        if not full_ok and (bond_all == bond_arom or bond_all == 0):
            return "AROMATIC_KEKULE_REPRESENTATION_EQUIVALENT", evidence
        if category == "CCD_DESCRIPTOR_MAPPING_FAILED" and bond_all > 0 and bond_all == bond_arom:
            return "AROMATIC_KEKULE_REPRESENTATION_EQUIVALENT", evidence
        if category == "CCD_DESCRIPTOR_MAPPING_FAILED" and bond_all == 0:
            return "ATOM_MAPPING_ONLY_MISMATCH", evidence
        if category == "KEKULIZATION_FAILED":
            return "TRUE_INVALID_FROZEN_GRAPH", evidence
        if category == "CCD_DESCRIPTOR_MAPPING_FAILED":
            return "CCD_DESCRIPTOR_TRUE_IDENTITY_CONFLICT", evidence

    if category == "CCD_DESCRIPTOR_MISSING_OR_BAD" and full_ok:
        raw_copy = Chem.Mol(sanitized)
        try:
            smi = canonical(raw_copy, True)
            parsed = Chem.MolFromSmiles(smi)
            roundtrip_ok = parsed is not None and canonical(parsed, True) == smi
        except Exception:
            roundtrip_ok = False
        evidence["authoritative_graph_smiles_roundtrip_ok"] = roundtrip_ok
        start_ok = False
        start_reason = ""
        if roundtrip_ok and evidence["ccd_declared_stereo_count"] == 0:
            try:
                start, _code = p4.independent_start(raw_copy, config)
                start_ok = canonical(start, True) == smi
            except Exception as exc:
                start_reason = f"{type(exc).__name__}: {exc}"
        evidence["audit_etkdg_uff_graph_roundtrip_ok"] = start_ok
        evidence["audit_etkdg_uff_error"] = start_reason
        if roundtrip_ok and start_ok and evidence["ccd_declared_stereo_count"] == 0:
            return "CCD_DESCRIPTOR_ONLY_WARNING", evidence
        if evidence["ccd_declared_stereo_count"] > 0:
            return "STEREO_ONE_SIDE_UNSPECIFIED", evidence
        return "OTHER_GRAPH_VALIDATION_REVIEW", evidence
    if category == "CCD_DESCRIPTOR_MISSING_OR_BAD":
        return "TRUE_INVALID_FROZEN_GRAPH", evidence
    return "OTHER_GRAPH_VALIDATION_REVIEW", evidence


def audit_bucket(run_dir: str, bucket_id: int) -> list[dict[str, Any]]:
    run = Path(run_dir)
    status_path = run / f"work/buckets/bucket_{bucket_id:03d}.parquet"
    statuses = pq.read_table(status_path).to_pandas()
    statuses = statuses[statuses["status"].eq(REVIEW_STATUS)].copy()
    if statuses.empty:
        return []
    inventory = pq.read_table(
        run / "input/case_inventory.parquet", filters=[("bucket_id", "=", bucket_id)]
    ).to_pandas()
    rows = statuses.merge(inventory, on=["case_id", "pair_id"], how="left", validate="one_to_one")
    atoms = p4.read_partition(p4.P2_RUN / "output/prepared_ligand_assembly_atoms", bucket_id)
    bonds = p4.read_partition(p4.P2_RUN / "output/prepared_ligand_assembly_bonds", bucket_id)
    store = p4.CCDStore(p4.P2_RUN / "input/ccd_active_snapshot.sqlite")
    config = p4.read_config(run / "input/config_snapshot.json")
    output = []
    for row in rows.itertuples():
        reason = str(row.reason)
        category = original_category(reason)
        placement = str(row.ligand_assembly_placement_id)
        component = str(row.component_id)
        lig_atoms = atoms[atoms["filter_2_ligand_assembly_placement_id"].astype(str).eq(placement)] if not atoms.empty else atoms
        lig_bonds = bonds[bonds["filter_2_ligand_assembly_placement_id"].astype(str).eq(placement)] if not bonds.empty else bonds
        base: dict[str, Any] = {
            "case_id": str(row.case_id), "pair_id": str(row.pair_id),
            "bucket_id": bucket_id, "component_id": component,
            "filter3_quality_class": str(row.filter3_quality_class),
            "previous_status": str(row.status), "previous_reason": reason,
            "original_reason_category": category,
            "audit_version": AUDIT_VERSION,
            "rescue_rule_applied": False, "new_status": str(row.status),
        }
        if category == "UPSTREAM_FROZEN_GRAPH_MISSING":
            base.update({
                "audit_v2_reason": "UPSTREAM_FROZEN_GRAPH_INCOMPLETE",
                "raw_graph_id": "", "observed_heavy_atoms": int(len(lig_atoms)),
                "raw_bond_count": int(len(lig_bonds)),
            })
            output.append(base)
            continue
        try:
            ccd = store.component(component)
            raw, raw_info = ccd_raw_molecule(lig_atoms, lig_bonds, ccd)
            sanitized, full_ok, partial_ok, sanitize_error = sanitize_probe(raw)
            base.update(raw_info)
            base["raw_sanitize_error"] = sanitize_error
            base["ccd_descriptor_sha256"] = hashlib.sha256(str(ccd.get("descriptor", "")).encode()).hexdigest()[:24]
            if category == "GRAPH_IDENTITY_OR_STEREO_MISMATCH":
                case_dir = run / f"output/cases/bucket_{bucket_id:03d}" / str(row.case_id)
                audit_reason, evidence = classify_existing_graph_mismatch(case_dir)
            else:
                audit_reason, evidence = classify_raw_review(
                    category, raw, sanitized, full_ok, partial_ok,
                    descriptor_molecule(ccd), ccd, config,
                )
            base["audit_v2_reason"] = audit_reason
            base.update(evidence)
        except Exception as exc:
            base["audit_v2_reason"] = "OTHER_GRAPH_VALIDATION_REVIEW"
            base["audit_exception"] = f"{type(exc).__name__}: {exc}"
            base.setdefault("raw_graph_id", "")
        output.append(base)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve()
    out = run / "review_audit_v2"
    out.mkdir(parents=True, exist_ok=True)
    buckets = list(range(256))
    all_rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(audit_bucket, str(run), bid): bid for bid in buckets}
        for future in as_completed(futures):
            bid = futures[future]
            rows = future.result()
            all_rows.extend(rows)
            print(json.dumps({"bucket_id": bid, "review_rows": len(rows)}), flush=True)
    frame = pd.DataFrame(all_rows).sort_values(["audit_v2_reason", "component_id", "case_id"])
    if len(frame) != 9223 or not frame["case_id"].is_unique:
        raise RuntimeError(f"audit universe closure failed: rows={len(frame)}, unique={frame['case_id'].nunique()}")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), out / "case_audit.parquet", compression="zstd")
    frame.to_csv(out / "case_audit.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = (
        frame.groupby(["original_reason_category", "audit_v2_reason"], dropna=False)
        .agg(pairs=("case_id", "size"), unique_ligand_graphs=("raw_graph_id", "nunique"),
             unique_ccd=("component_id", "nunique"))
        .reset_index().sort_values(["original_reason_category", "pairs"], ascending=[True, False])
    )
    summary.to_csv(out / "reason_summary.tsv", sep="\t", index=False)
    concentration = (
        frame.groupby(["original_reason_category", "audit_v2_reason", "component_id"], dropna=False)
        .agg(pairs=("case_id", "size"), unique_ligand_graphs=("raw_graph_id", "nunique"))
        .reset_index().sort_values("pairs", ascending=False)
    )
    concentration.to_csv(out / "component_concentration.tsv", sep="\t", index=False)
    payload = {
        "audit_version": AUDIT_VERSION,
        "created_at": utc(),
        "read_only_audit": True,
        "baseline_ready_untouched": 148340,
        "review_cases_audited": int(len(frame)),
        "audit_reason_counts": {str(k): int(v) for k, v in frame["audit_v2_reason"].value_counts().items()},
        "original_reason_counts": {str(k): int(v) for k, v in frame["original_reason_category"].value_counts().items()},
        "rescue_status_changes_applied": 0,
        "rdkit_version": rdkit.__version__,
    }
    (out / "audit_summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    provenance = {
        "audit_version": AUDIT_VERSION, "created_at": utc(),
        "source_run": str(run),
        "source_validation_sha256": p4.sha256_file(run / "validation/validation.json"),
        "source_case_inventory_sha256": p4.sha256_file(run / "input/case_inventory.parquet"),
        "constraints": [
            "No existing P4_DOCKING_READY case was read for regeneration",
            "No bucket status or case directory was modified",
            "No frozen ligand identity was changed",
            "No native-coordinate fallback was used",
            "P4_LIGAND_START_GENERATION_FAILED cases were excluded",
        ],
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
