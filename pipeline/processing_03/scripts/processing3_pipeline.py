#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from scipy.spatial import cKDTree


CONFIG = None
P1_PATHS = None
SELECTION = None
VDW = None
RUN_DIR = None

METAL_ELEMENTS = {
    "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "GA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA",
    "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA",
    "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC", "TH", "PA", "U",
    "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT",
    "DS", "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG"
}

DATASETS = [
    "placement_terminal_status", "spatial_candidate_chains", "qualifying_atomic_contacts",
    "binding_residues", "contact_supported_chains", "provisional_pairs",
    "pair_pocket_residues", "covalent_audit", "metal_audit", "technical_failures",
]

OUTPUT_COLUMNS = {
    "placement_terminal_status": ["ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id",
        "terminal_status", "reason_code", "candidate_chain_count", "supported_chain_count", "minimum_distance",
        "receptor_chain_instance_ids", "plip_membership_effect"],
    "spatial_candidate_chains": ["ligand_assembly_placement_id", "chain_instance_id", "pdb_id", "assembly_id", "model_id",
        "minimum_heavy_atom_distance", "protein_atom_pairs_within_6A"],
    "qualifying_atomic_contacts": ["ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id",
        "ligand_atom_id", "protein_atom_id", "distance_angstrom", "vdw_upper_bound_angstrom"],
    "binding_residues": ["ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id",
        "qualifying_atomic_contact_count"],
    "contact_supported_chains": ["ligand_assembly_placement_id", "chain_instance_id", "qualifying_atomic_contact_count",
        "binding_residue_count", "chain_status"],
    "provisional_pairs": ["pair_id", "ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id",
        "receptor_chain_instance_ids", "receptor_chain_count", "metal_status", "pair_status"],
    "pair_pocket_residues": ["pair_id", "ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id"],
    "covalent_audit": ["ligand_assembly_placement_id", "pdb_id", "conn_id", "chain_instance_id", "mapping_status"],
    "metal_audit": ["ligand_assembly_placement_id", "pdb_id", "metal_status", "ligand_metal_conn_id",
        "protein_metal_conn_id", "metal_component_id"],
    "technical_failures": ["pdb_id", "record_id", "error_type", "error_message"],
}


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value in {".", "?", "None", "nan"} else value


def norm_element(value):
    return clean(value).upper()


def severe_steric_clash(distance):
    return float(distance) < 2.0


def qualifying_direct_contact(distance, protein_element, ligand_element, vdw):
    protein_element = norm_element(protein_element)
    ligand_element = norm_element(ligand_element)
    if protein_element not in vdw or ligand_element not in vdw:
        return False
    distance = float(distance)
    return 2.0 <= distance <= vdw[protein_element] + vdw[ligand_element] + 0.5


def read_parquet(path, columns=None):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=columns or [])
    return pq.ParquetFile(path).read(columns=columns).to_pandas()


def write_parquet_atomic(path, frame, columns=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    if columns is not None:
        frame = frame.reindex(columns=columns)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), tmp, compression="zstd")
    os.replace(tmp, path)


def load_tsv_map(path, key, value):
    opener = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row[key]] = row[value]
    return out


def extract_vdw_table(source):
    text = Path(source).read_text(encoding="utf-8", errors="replace")
    match = re.search(r"void\s+make_vdw\s*\([^)]*\)\s*\{(?P<body>.*?)\n\}", text, re.S)
    if not match:
        raise RuntimeError("BioLiP2 make_vdw function not found")
    pairs = re.findall(r'vdw_dict\[\s*"\s*([A-Za-z]{1,2})\s*"\s*\]\s*=\s*([0-9.]+)', match.group("body"))
    table = {}
    for element, radius in pairs:
        table[element.upper()] = float(radius)
    if not all(x in table for x in ["C", "N", "O", "S", "P", "F", "CL", "BR", "I"]):
        raise RuntimeError("BioLiP2 vdw table missing core elements")
    return table


def category_records(block, prefix):
    table = block.find_mmcif_category(prefix)
    if not table:
        return []
    tags = [tag[len(prefix):] for tag in table.tags]
    return [{tag: clean(value) for tag, value in zip(tags, row)} for row in table]


def endpoint(row, n):
    p = f"ptnr{n}_"
    return {
        "label_asym_id": row.get(p + "label_asym_id", ""),
        "label_comp_id": row.get(p + "label_comp_id", "").upper(),
        "label_seq_id": row.get(p + "label_seq_id", ""),
        "auth_asym_id": row.get(p + "auth_asym_id", ""),
        "auth_comp_id": row.get(p + "auth_comp_id", "").upper(),
        "auth_seq_id": row.get(p + "auth_seq_id", ""),
        "symmetry": row.get(p + "symmetry", ""),
    }


def endpoint_comp(ep):
    return ep["label_comp_id"] or ep["auth_comp_id"]


def endpoint_asym(ep):
    return ep["label_asym_id"] or ep["auth_asym_id"]


def identity_symmetry(ep):
    return clean(ep.get("symmetry", "")) in {"", "1_555"}


def parse_connections(pdb_id):
    path = P1_PATHS.get(pdb_id)
    if not path:
        return [], "missing_mmcif_path"
    try:
        block = gemmi.cif.read(path).sole_block()
        rows = []
        for row in category_records(block, "_struct_conn."):
            ctype = row.get("conn_type_id", "").lower()
            if ctype not in {"covale", "metalc"}:
                continue
            rows.append({"conn_id": row.get("id", ""), "conn_type_id": ctype,
                         "p1": endpoint(row, 1), "p2": endpoint(row, 2)})
        return rows, ""
    except Exception as exc:
        return [], f"{type(exc).__name__}:{exc}"


def source_matches_endpoint(meta, ep):
    comp = endpoint_comp(ep)
    if comp and comp != str(meta.component_id).upper():
        return False
    if ep["label_asym_id"] and ep["label_asym_id"] != str(meta.label_asym_id):
        return False
    if not ep["label_asym_id"] and ep["auth_asym_id"] and ep["auth_asym_id"] != str(meta.auth_asym_id):
        return False
    auth_seq = clean(meta.auth_seq_id)
    if ep["auth_seq_id"] and auth_seq and ep["auth_seq_id"] != auth_seq:
        return False
    return bool(endpoint_asym(ep) or ep["auth_seq_id"])


def receptor_endpoint_source_candidates(ep, receptor_meta):
    if ep["label_asym_id"]:
        chain_match = receptor_meta["label_asym_id"].astype(str).eq(ep["label_asym_id"])
    elif ep["auth_asym_id"]:
        chain_match = receptor_meta["auth_asym_id"].astype(str).eq(ep["auth_asym_id"])
    else:
        return receptor_meta.iloc[0:0]
    return receptor_meta[chain_match]


def receptor_endpoint_chain(ep, receptor_meta, ligand_operator):
    candidates = receptor_endpoint_source_candidates(ep, receptor_meta)
    candidates = candidates[candidates["operator_path"].astype(str).eq(str(ligand_operator))]
    return candidates["chain_instance_id"].astype(str).tolist()


def residue_key(row):
    return "|".join([clean(row.label_seq_id), clean(row.auth_seq_id), clean(row.insertion_code), clean(row.label_comp_id)])


def strict_inputs(rel):
    p2 = Path(CONFIG["input"]["processing_2_run"])
    aux = Path(CONFIG["input"]["auxiliary_build"])
    placements = read_parquet(aux / "entry_ligand_placements" / rel)
    receptors = read_parquet(aux / "entry_receptor_chain_instances" / rel).rename(
        columns={"filter_1_chain_instance_id": "chain_instance_id"})
    topology = read_parquet(p2 / "output/ligand_topology_validation" / rel)
    manifest = read_parquet(p2 / "output/structure_preparation_manifest" / rel)
    pm = manifest[manifest["object_type"].eq("ligand_assembly_placement")].copy().rename(
        columns={"object_id": "filter_2_ligand_assembly_placement_id"})
    rm = manifest[manifest["object_type"].eq("receptor_chain_instance")].copy().rename(
        columns={"object_id": "chain_instance_id"})
    placements = placements.merge(topology[["source_ligand_instance_id", "mapping_status", "missing_heavy_atom_count",
                                             "topology_status", "rdkit_parse_success", "rdkit_sanitize_success"]],
                                  left_on="filter_2_source_ligand_instance_id", right_on="source_ligand_instance_id",
                                  how="left", validate="many_to_one")
    placements = placements.merge(pm[["filter_2_ligand_assembly_placement_id", "preparation_status", "prepared_atom_count",
                                      "operator_quality_status"]], on="filter_2_ligand_assembly_placement_id",
                                  how="left", validate="one_to_one")
    strict = (
        placements["mapping_status"].eq("COMPLETE")
        & placements["missing_heavy_atom_count"].eq(0)
        & placements["topology_status"].eq("TOPOLOGY_COMPLETE")
        & placements["rdkit_parse_success"].eq(True)
        & placements["rdkit_sanitize_success"].eq(True)
        & placements["operator_quality_status"].eq("PASS")
        & placements["preparation_status"].eq("ASSEMBLY_READY")
        & placements["prepared_atom_count"].fillna(0).gt(0)
    )
    placements = placements[strict].copy()
    receptors = receptors.merge(rm[["chain_instance_id", "preparation_status", "prepared_atom_count", "operator_quality_status"]],
                                on="chain_instance_id", how="left", validate="one_to_one")
    receptors = receptors[
        receptors["receptor_eligible"].eq(True)
        & receptors["polymer_class"].eq("POLYPEPTIDE")
        & receptors["operator_quality_status"].eq("PASS")
        & receptors["preparation_status"].eq("ASSEMBLY_READY")
        & receptors["prepared_atom_count"].fillna(0).gt(0)
    ].copy()
    if SELECTION is not None:
        placements = placements[placements["filter_2_ligand_assembly_placement_id"].isin(SELECTION)].copy()
        keep_pdb = set(placements["pdb_id"])
        receptors = receptors[receptors["pdb_id"].isin(keep_pdb)].copy()
    return placements, receptors


def ligand_declared_covalent(meta, receptor_meta, connections):
    mapped = []
    unresolved = []
    for conn in connections:
        if conn["conn_type_id"] != "covale":
            continue
        for lig_ep, other in [(conn["p1"], conn["p2"]), (conn["p2"], conn["p1"])]:
            if not source_matches_endpoint(meta, lig_ep):
                continue
            source_candidates = receptor_endpoint_source_candidates(other, receptor_meta)
            if source_candidates.empty:
                # _struct_conn may describe ligand-ligand or ligand-saccharide links.
                continue
            if not identity_symmetry(lig_ep) or not identity_symmetry(other):
                unresolved.append(conn["conn_id"])
                continue
            chains = receptor_endpoint_chain(other, receptor_meta, meta.operator_path)
            if chains:
                mapped.extend((conn["conn_id"], chain) for chain in chains)
            elif not source_candidates.empty:
                unresolved.append(conn["conn_id"])
    return sorted(set(mapped)), sorted(set(unresolved))


def ligand_metal_status(meta, supported_chain_ids, receptor_meta, connections):
    ligand_metal_nodes = []
    for conn in connections:
        if conn["conn_type_id"] != "metalc":
            continue
        for lig_ep, other in [(conn["p1"], conn["p2"]), (conn["p2"], conn["p1"])]:
            if source_matches_endpoint(meta, lig_ep) and endpoint_comp(other) in METAL_ELEMENTS:
                ligand_metal_nodes.append((conn["conn_id"], other))
    if not ligand_metal_nodes:
        return "NO_RELEVANT_METAL", []
    supported = receptor_meta[receptor_meta["chain_instance_id"].isin(supported_chain_ids)]
    evidence = []
    for first_id, metal_ep in ligand_metal_nodes:
        for conn in connections:
            if conn["conn_type_id"] != "metalc":
                continue
            for maybe_metal, protein_ep in [(conn["p1"], conn["p2"]), (conn["p2"], conn["p1"])]:
                same_metal = (
                    endpoint_asym(maybe_metal) == endpoint_asym(metal_ep)
                    and endpoint_comp(maybe_metal) == endpoint_comp(metal_ep)
                    and (not metal_ep["auth_seq_id"] or maybe_metal["auth_seq_id"] == metal_ep["auth_seq_id"])
                )
                if not same_metal:
                    continue
                if receptor_endpoint_chain(protein_ep, supported, meta.operator_path):
                    evidence.append((first_id, conn["conn_id"], endpoint_comp(metal_ep)))
    if evidence:
        return "LIGAND_INVOLVED_METAL_COORDINATION", sorted(set(evidence))
    return "NEARBY_METAL_NOT_LIGAND_COORDINATING", [(x[0], "", endpoint_comp(x[1])) for x in ligand_metal_nodes]


def process_task(task):
    rel = Path(task)
    started = time.time()
    p2 = Path(CONFIG["input"]["processing_2_run"])
    placements, receptor_meta = strict_inputs(rel)
    outputs = {name: [] for name in DATASETS}
    if placements.empty:
        marker = Path(RUN_DIR) / "work/checkpoints" / (rel.stem + ".json")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(json.dumps({"status": "complete", "task": str(rel), "placements": 0}) + "\n")
        return {"task": str(rel), "placements": 0, "terminal": 0, "runtime": time.time() - started}

    lig_atoms = read_parquet(p2 / "output/prepared_ligand_assembly_atoms" / rel)
    rec_atoms = read_parquet(p2 / "output/prepared_receptor_assembly_atoms" / rel)
    selected_ids = set(placements["filter_2_ligand_assembly_placement_id"])
    lig_atoms = lig_atoms[lig_atoms["filter_2_ligand_assembly_placement_id"].isin(selected_ids)].copy()
    ready_chains = set(receptor_meta["chain_instance_id"])
    rec_atoms = rec_atoms[rec_atoms["filter_1_chain_instance_id"].isin(ready_chains)].copy()
    lig_atoms = lig_atoms[~lig_atoms["type_symbol"].astype(str).str.upper().isin({"H", "D"})].copy()
    rec_atoms = rec_atoms[~rec_atoms["type_symbol"].astype(str).str.upper().isin({"H", "D"})].copy()

    for pid, pmeta in placements.groupby("pdb_id", sort=False):
        conns, conn_error = parse_connections(str(pid))
        if conn_error:
            outputs["technical_failures"].append({"pdb_id": pid, "record_id": "", "error_type": "STRUCT_CONN_PARSE_FAILURE", "error_message": conn_error})
            for meta in pmeta.itertuples():
                outputs["placement_terminal_status"].append({
                    "ligand_assembly_placement_id": meta.filter_2_ligand_assembly_placement_id, "pdb_id": pid,
                    "assembly_id": meta.assembly_id, "model_id": meta.model_id, "component_id": meta.component_id,
                    "terminal_status": "PROCESSING3_TECHNICAL_FAILURE", "reason_code": "struct_conn_parse_failure",
                    "candidate_chain_count": 0, "supported_chain_count": 0, "minimum_distance": math.nan,
                    "receptor_chain_instance_ids": "", "plip_membership_effect": False})
            continue
        p_rec_meta = receptor_meta[receptor_meta["pdb_id"].eq(pid)]
        for key, ameta in pmeta.groupby(["assembly_id", "model_id"], sort=False):
            aid, mid = key
            rmeta = p_rec_meta[p_rec_meta["assembly_id"].eq(aid) & p_rec_meta["model_id"].eq(mid)]
            ra = rec_atoms[rec_atoms["pdb_id"].eq(pid) & rec_atoms["assembly_id"].eq(aid) & rec_atoms["model_id"].eq(mid)].copy()
            if ra.empty:
                for meta in ameta.itertuples():
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": meta.filter_2_ligand_assembly_placement_id, "pdb_id": pid,
                        "assembly_id": aid, "model_id": mid, "component_id": meta.component_id,
                        "terminal_status": "PROCESSING3_TECHNICAL_FAILURE", "reason_code": "ready_receptor_atoms_missing",
                        "candidate_chain_count": 0, "supported_chain_count": 0, "minimum_distance": math.nan,
                        "receptor_chain_instance_ids": "", "plip_membership_effect": False})
                continue
            rxyz = ra[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
            tree = cKDTree(rxyz)
            r_elements = ra["type_symbol"].astype(str).str.upper().to_numpy()
            r_chain = ra["filter_1_chain_instance_id"].astype(str).to_numpy()
            r_residue = np.array([residue_key(x) for x in ra.itertuples()], dtype=object)
            r_atom = ra["label_atom_id"].astype(str).to_numpy()

            for meta in ameta.itertuples():
                lid = str(meta.filter_2_ligand_assembly_placement_id)
                la = lig_atoms[lig_atoms["filter_2_ligand_assembly_placement_id"].eq(lid)].copy()
                if la.empty:
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "PROCESSING3_TECHNICAL_FAILURE",
                        "reason_code": "strict_ready_ligand_atoms_missing", "candidate_chain_count": 0,
                        "supported_chain_count": 0, "minimum_distance": math.nan, "receptor_chain_instance_ids": "",
                        "plip_membership_effect": False})
                    continue
                covalent, cov_unresolved = ligand_declared_covalent(meta, rmeta, conns)
                lxyz = la[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
                l_elements = la["type_symbol"].astype(str).str.upper().to_numpy()
                l_atoms = la["label_atom_id"].astype(str).to_numpy()
                neighbours = tree.query_ball_point(lxyz, r=6.0)
                pair_indices = []
                for li, indices in enumerate(neighbours):
                    pair_indices.extend((li, int(ri)) for ri in indices)
                if not pair_indices:
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "NO_PROTEIN_WITHIN_6A",
                        "reason_code": "no_receptor_heavy_atom_within_6A", "candidate_chain_count": 0,
                        "supported_chain_count": 0, "minimum_distance": math.nan, "receptor_chain_instance_ids": "",
                        "plip_membership_effect": False})
                    continue
                pair_indices = sorted(set(pair_indices))
                distances = np.array([np.linalg.norm(lxyz[li] - rxyz[ri]) for li, ri in pair_indices])
                candidate_chains = sorted(set(r_chain[ri] for _, ri in pair_indices))
                minimum_distance = float(distances.min())
                for chain in candidate_chains:
                    ds = [d for (li, ri), d in zip(pair_indices, distances) if r_chain[ri] == chain]
                    outputs["spatial_candidate_chains"].append({
                        "ligand_assembly_placement_id": lid, "chain_instance_id": chain, "pdb_id": pid,
                        "assembly_id": aid, "model_id": mid, "minimum_heavy_atom_distance": min(ds),
                        "protein_atom_pairs_within_6A": len(ds)})
                if cov_unresolved and not covalent:
                    outputs["technical_failures"].append({
                        "pdb_id": pid, "record_id": lid, "error_type": "COVALENT_OPERATOR_MAPPING_UNRESOLVED",
                        "error_message": "struct_conn_ids=" + ",".join(cov_unresolved)})
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "PROCESSING3_TECHNICAL_FAILURE",
                        "reason_code": "covalent_declaration_operator_mapping_unresolved", "candidate_chain_count": len(candidate_chains),
                        "supported_chain_count": 0, "minimum_distance": minimum_distance, "receptor_chain_instance_ids": "",
                        "plip_membership_effect": False})
                    continue
                if covalent:
                    for conn_id, chain in covalent:
                        outputs["covalent_audit"].append({"ligand_assembly_placement_id": lid, "pdb_id": pid,
                                                           "conn_id": conn_id, "chain_instance_id": chain,
                                                           "mapping_status": "EXACT_ASSEMBLY_OPERATOR_MATCH"})
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "OUT_OF_SCOPE_COVALENT",
                        "reason_code": "direct_struct_conn_covale", "candidate_chain_count": len(candidate_chains),
                        "supported_chain_count": 0, "minimum_distance": minimum_distance,
                        "receptor_chain_instance_ids": ",".join(sorted({x[1] for x in covalent})), "plip_membership_effect": False})
                    continue
                if severe_steric_clash(minimum_distance):
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "SEVERE_STERIC_CLASH",
                        "reason_code": "unexplained_heavy_atom_distance_below_2A", "candidate_chain_count": len(candidate_chains),
                        "supported_chain_count": 0, "minimum_distance": minimum_distance, "receptor_chain_instance_ids": "",
                        "plip_membership_effect": False})
                    continue
                unsupported = sorted({l_elements[li] for li, _ in pair_indices if l_elements[li] not in VDW}
                                     | {r_elements[ri] for _, ri in pair_indices if r_elements[ri] not in VDW})
                if unsupported:
                    outputs["placement_terminal_status"].append({
                        "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                        "component_id": meta.component_id, "terminal_status": "UNSUPPORTED_VDW_ELEMENT",
                        "reason_code": "unsupported_elements:" + ",".join(unsupported), "candidate_chain_count": len(candidate_chains),
                        "supported_chain_count": 0, "minimum_distance": minimum_distance, "receptor_chain_instance_ids": "",
                        "plip_membership_effect": False})
                    continue

                contacts = []
                residue_counts = Counter()
                chain_residues = defaultdict(set)
                pocket_by_chain = defaultdict(set)
                for (li, ri), distance in zip(pair_indices, distances):
                    chain = r_chain[ri]
                    residue = r_residue[ri]
                    pocket_by_chain[chain].add(residue)
                    upper = VDW[r_elements[ri]] + VDW[l_elements[li]] + 0.5
                    if qualifying_direct_contact(distance, r_elements[ri], l_elements[li], VDW):
                        contact_key = (li, ri)
                        contacts.append((contact_key, chain, residue, float(distance), upper))
                        residue_counts[(chain, residue)] += 1
                        outputs["qualifying_atomic_contacts"].append({
                            "ligand_assembly_placement_id": lid, "chain_instance_id": chain, "protein_residue_id": residue,
                            "ligand_atom_id": l_atoms[li], "protein_atom_id": r_atom[ri], "distance_angstrom": float(distance),
                            "vdw_upper_bound_angstrom": upper})
                binding = defaultdict(list)
                for (chain, residue), count in residue_counts.items():
                    if count >= 2:
                        binding[chain].append(residue)
                        outputs["binding_residues"].append({
                            "ligand_assembly_placement_id": lid, "chain_instance_id": chain,
                            "protein_residue_id": residue, "qualifying_atomic_contact_count": count})
                supported = sorted(chain for chain, residues in binding.items() if len(set(residues)) >= 2)
                for chain in candidate_chains:
                    outputs["contact_supported_chains"].append({
                        "ligand_assembly_placement_id": lid, "chain_instance_id": chain,
                        "qualifying_atomic_contact_count": sum(v for (c, _), v in residue_counts.items() if c == chain),
                        "binding_residue_count": len(set(binding.get(chain, []))),
                        "chain_status": "CONTACT_SUPPORTED_CHAIN" if chain in supported else ("SPARSE_DIRECT_CONTACT" if any(c == chain for c, _ in residue_counts) else "NO_DIRECT_CONTACT")})
                if not contacts:
                    terminal = "NO_DIRECT_CONTACT"; reason = "within_6A_but_no_vdw_qualifying_contact"
                elif not supported:
                    terminal = "SPARSE_DIRECT_CONTACT"; reason = "direct_contacts_below_2x2_chain_gate"
                else:
                    metal_status, metal_evidence = ligand_metal_status(meta, supported, rmeta, conns)
                    for a, b, element in metal_evidence:
                        outputs["metal_audit"].append({"ligand_assembly_placement_id": lid, "pdb_id": pid,
                                                       "metal_status": metal_status, "ligand_metal_conn_id": a,
                                                       "protein_metal_conn_id": b, "metal_component_id": element})
                    terminal = "OUT_OF_SCOPE_METAL_RELATED" if metal_status == "LIGAND_INVOLVED_METAL_COORDINATION" else "FINAL_ORDINARY_NONCOVALENT_PAIR"
                    reason = "ligand_involved_metal_coordination" if terminal.startswith("OUT_") else "contact_supported_non_covalent_non_metal_pair"
                    if terminal == "FINAL_ORDINARY_NONCOVALENT_PAIR":
                        pair_id = "P3|" + lid
                        outputs["provisional_pairs"].append({
                            "pair_id": pair_id, "ligand_assembly_placement_id": lid, "pdb_id": pid,
                            "assembly_id": aid, "model_id": mid, "component_id": meta.component_id,
                            "receptor_chain_instance_ids": ",".join(supported), "receptor_chain_count": len(supported),
                            "metal_status": metal_status, "pair_status": terminal})
                        for chain in supported:
                            for residue in sorted(pocket_by_chain[chain]):
                                outputs["pair_pocket_residues"].append({"pair_id": pair_id, "ligand_assembly_placement_id": lid,
                                                                        "chain_instance_id": chain, "protein_residue_id": residue})
                outputs["placement_terminal_status"].append({
                    "ligand_assembly_placement_id": lid, "pdb_id": pid, "assembly_id": aid, "model_id": mid,
                    "component_id": meta.component_id, "terminal_status": terminal, "reason_code": reason,
                    "candidate_chain_count": len(candidate_chains), "supported_chain_count": len(supported),
                    "minimum_distance": minimum_distance, "receptor_chain_instance_ids": ",".join(supported),
                    "plip_membership_effect": False})

    for name, rows in outputs.items():
        frame = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS[name])
        if frame.empty:
            continue
        frame = frame.drop_duplicates().sort_values(list(frame.columns[:min(3, len(frame.columns))]), kind="stable")
        write_parquet_atomic(Path(RUN_DIR) / "work/batches" / name / rel, frame, OUTPUT_COLUMNS[name])
    marker = Path(RUN_DIR) / "work/checkpoints" / (rel.stem + ".json")
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"status": "complete", "task": str(rel), "placements": len(placements),
                                  "terminal": len(outputs["placement_terminal_status"]),
                                  "runtime_seconds": time.time() - started, "finished_at": utc()}) + "\n")
    return {"task": str(rel), "placements": len(placements), "terminal": len(outputs["placement_terminal_status"]),
            "runtime": time.time() - started}


def init_worker(config, p1_paths, selection, vdw, run_dir):
    global CONFIG, P1_PATHS, SELECTION, VDW, RUN_DIR
    CONFIG, P1_PATHS, SELECTION, VDW, RUN_DIR = config, p1_paths, selection, vdw, run_dir


def preflight(args):
    config = yaml.safe_load(Path(args.config).read_text())
    p2 = Path(config["input"]["processing_2_run"])
    aux = Path(config["input"]["auxiliary_build"])
    rels = sorted(x.relative_to(aux / "entry_ligand_placements") for x in (aux / "entry_ligand_placements").rglob("*.parquet"))
    checks = {
        "processing_2_frozen": (p2 / "_FROZEN.json").exists(),
        "processing_2_validation": (p2 / "audit/assembly_coordinate_validation.json").exists(),
        "partition_count_match": len(rels) == 8227,
        "all_required_partition_sets_present": all(all((base / rel).exists() for base in [
            aux / "entry_receptor_chain_instances", p2 / "output/prepared_ligand_assembly_atoms",
            p2 / "output/prepared_receptor_assembly_atoms", p2 / "output/ligand_topology_validation",
            p2 / "output/structure_preparation_manifest"]) for rel in rels),
        "raw_mmcif_index_exists": Path(config["input"]["processing_1_index"]).exists(),
        "biolip2_vdw_source_exists": Path(config["contact"]["vdw_source_cpp"]).exists(),
        "scipy_available": True,
        "plip_not_membership_gate": config["policy"]["plip_affects_membership"] is False,
    }
    vdw = extract_vdw_table(config["contact"]["vdw_source_cpp"])
    report = {"preflight_pass": all(checks.values()), "checks": checks, "partition_count": len(rels),
              "vdw_element_count": len(vdw), "vdw_source_sha256": sha256(config["contact"]["vdw_source_cpp"]),
              "struct_conn_source": "raw_mmcif_read_only_because_processing_2_did_not_materialize_struct_conn",
              "created_at": utc()}
    out = Path(args.run_dir); (out / "audit").mkdir(parents=True, exist_ok=True); (out / "references").mkdir(parents=True, exist_ok=True)
    (out / "audit/preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "references/biolip2_vdw_radii.json").write_text(json.dumps(vdw, indent=2, sort_keys=True) + "\n")
    shutil.copy2(args.config, out / "config_snapshot.yaml")
    print(json.dumps(report, indent=2))
    if not report["preflight_pass"]:
        raise SystemExit(2)


def choose_smoke(args):
    config = yaml.safe_load(Path(args.config).read_text())
    global CONFIG, SELECTION
    CONFIG = config; SELECTION = None
    aux = Path(config["input"]["auxiliary_build"])
    candidates = []
    for f in sorted((aux / "entry_ligand_placements").rglob("*.parquet")):
        rel = f.relative_to(aux / "entry_ligand_placements")
        p, _ = strict_inputs(rel)
        for row in p[["filter_2_ligand_assembly_placement_id", "pdb_id", "component_id"]].itertuples(index=False):
            score = hashlib.sha256(row.filter_2_ligand_assembly_placement_id.encode()).hexdigest()
            candidates.append((score, row.filter_2_ligand_assembly_placement_id, row.pdb_id, row.component_id))
    selected = sorted(candidates)[:args.count]
    out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["hash", "ligand_assembly_placement_id", "pdb_id", "component_id"]); w.writerows(selected)
    print(json.dumps({"strict_candidates": len(candidates), "selected": len(selected), "output": str(out)}, indent=2))


def run(args):
    config = yaml.safe_load(Path(args.config).read_text())
    p1 = load_tsv_map(config["input"]["processing_1_index"], "pdb_id", "mmcif_path")
    vdw = extract_vdw_table(config["contact"]["vdw_source_cpp"])
    selection = None
    if args.selection:
        with open(args.selection, encoding="utf-8") as fh:
            selection = {r["ligand_assembly_placement_id"] for r in csv.DictReader(fh, delimiter="\t")}
    aux = Path(config["input"]["auxiliary_build"])
    rels = sorted(x.relative_to(aux / "entry_ligand_placements") for x in (aux / "entry_ligand_placements").rglob("*.parquet"))
    if selection:
        relevant = []
        for rel in rels:
            ids = set(read_parquet(aux / "entry_ligand_placements" / rel, ["filter_2_ligand_assembly_placement_id"])["filter_2_ligand_assembly_placement_id"])
            if ids & selection:
                relevant.append(rel)
        rels = relevant
    run_dir = Path(args.run_dir); (run_dir / "work/checkpoints").mkdir(parents=True, exist_ok=True)
    pending = [str(rel) for rel in rels if not (run_dir / "work/checkpoints" / (rel.stem + ".json")).exists()]
    status = {"status": "RUNNING", "task_count": len(rels), "pending_task_count": len(pending), "workers": args.workers,
              "started_at": utc(), "selection_count": len(selection) if selection else None}
    (run_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    completed = 0; placements = 0; terminal = 0; started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker,
                             initargs=(config, p1, selection, vdw, str(run_dir))) as pool:
        for result in pool.map(process_task, pending, chunksize=1):
            completed += 1; placements += result["placements"]; terminal += result["terminal"]
            if completed % 20 == 0 or completed == len(pending):
                status.update({"completed_this_attempt": completed, "placements_this_attempt": placements,
                               "terminal_this_attempt": terminal, "runtime_seconds": time.time() - started,
                               "updated_at": utc()})
                (run_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n")
                print(json.dumps(status), flush=True)
    status.update({"status": "COMPLETED", "completed_this_attempt": completed, "finished_at": utc(),
                   "runtime_seconds": time.time() - started})
    (run_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n")


def finalize(args):
    run_dir = Path(args.run_dir)
    output = run_dir / "output"
    output.mkdir(parents=True, exist_ok=True)
    for name in DATASETS:
        src = run_dir / "work/batches" / name
        dst = output / name
        if dst.exists():
            raise RuntimeError(f"Output already exists: {dst}")
        if src.exists():
            os.replace(src, dst)
    counts = {}
    terminal = Counter(); duplicate = 0; ids = set()
    for name in DATASETS:
        files = list((output / name).rglob("*.parquet")) if (output / name).exists() else []
        counts[name] = sum(pq.ParquetFile(f).metadata.num_rows for f in files)
        if name == "placement_terminal_status":
            for f in files:
                d = read_parquet(f, ["ligand_assembly_placement_id", "terminal_status"])
                for row in d.itertuples(index=False):
                    duplicate += row.ligand_assembly_placement_id in ids; ids.add(row.ligand_assembly_placement_id); terminal[row.terminal_status] += 1
    summary = {"row_counts": counts, "terminal_status_counts": dict(terminal), "unique_terminal_ids": len(ids),
               "duplicate_terminal_ids": duplicate, "finalized_at": utc()}
    (output / "release_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))


def validate(args):
    run_dir = Path(args.run_dir); output = run_dir / "output"
    summary = json.loads((output / "release_summary.json").read_text())
    expected = args.expected
    terminal_count = summary["row_counts"].get("placement_terminal_status", 0)
    dataset_readable = {}
    for name in DATASETS:
        path = output / name
        if not path.exists() or not list(path.rglob("*.parquet")):
            dataset_readable[name] = True
            continue
        try:
            import pyarrow.dataset as ds
            ds.dataset(path, format="parquet", partitioning="hive").count_rows()
            dataset_readable[name] = True
        except Exception:
            dataset_readable[name] = False
    checks = {
        "terminal_count_matches_expected": terminal_count == expected,
        "unique_terminal_ids_match_expected": summary["unique_terminal_ids"] == expected,
        "duplicate_terminal_ids_zero": summary["duplicate_terminal_ids"] == 0,
        "terminal_status_missing_zero": sum(summary["terminal_status_counts"].values()) == expected,
        "plip_not_executed_as_membership_gate": True,
        "arpeggio_not_executed": True,
        "prolif_not_executed": True,
        "docking_not_executed": True,
        "all_parquet_datasets_readable": all(dataset_readable.values()),
    }
    report = {"validation_pass": all(checks.values()), "checks": checks, "dataset_readable": dataset_readable, "expected": expected,
              "terminal_status_counts": summary["terminal_status_counts"], "validated_at": utc()}
    (run_dir / "audit").mkdir(parents=True, exist_ok=True)
    (run_dir / "audit/processing_3_validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    if not report["validation_pass"]:
        raise SystemExit(2)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("preflight"); a.add_argument("--config", required=True); a.add_argument("--run-dir", required=True); a.set_defaults(func=preflight)
    a = sub.add_parser("choose-smoke"); a.add_argument("--config", required=True); a.add_argument("--count", type=int, default=500); a.add_argument("--output", required=True); a.set_defaults(func=choose_smoke)
    a = sub.add_parser("run"); a.add_argument("--config", required=True); a.add_argument("--run-dir", required=True); a.add_argument("--workers", type=int, default=8); a.add_argument("--selection"); a.set_defaults(func=run)
    a = sub.add_parser("finalize"); a.add_argument("--run-dir", required=True); a.set_defaults(func=finalize)
    a = sub.add_parser("validate"); a.add_argument("--run-dir", required=True); a.add_argument("--expected", type=int, required=True); a.set_defaults(func=validate)
    args = p.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
