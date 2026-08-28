#!/usr/bin/env python3
"""Map PDBbind ligand SDF coordinates to deposited PDB ligand instances."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import io
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from scipy.optimize import linear_sum_assignment

RDLogger.DisableLog("rdApp.*")


def pdb_element(line: str) -> str:
    value = line[76:78].strip().upper() if len(line) >= 78 else ""
    if value:
        return value
    name = line[12:16].strip().upper().lstrip("0123456789")
    return name[:2] if name[:2] in {"CL", "BR"} else name[:1]


def deposited_candidates(path: Path):
    residues = defaultdict(lambda: defaultdict(list))
    with gzip.open(path, "rt", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("HETATM") or len(line) < 54:
                continue
            element = pdb_element(line)
            if element in {"H", "D"}:
                continue
            component = line[17:20].strip()
            if component in {"HOH", "DOD"}:
                continue
            physical = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), component)
            alt = line[16:17].strip()
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            residues[physical][alt].append((element, xyz))
    result = []
    for physical, alts in residues.items():
        common = alts.get("", [])
        explicit = sorted(key for key in alts if key)
        if not explicit:
            result.append((*physical, "", common))
        else:
            for alt in explicit:
                result.append((*physical, alt, common + alts[alt]))
    return result


def sdf_coordinates(block: str):
    molecule = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=False)
    if molecule is None or molecule.GetNumConformers() == 0:
        return None
    conformer = molecule.GetConformer()
    elements, coords = [], []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        point = conformer.GetAtomPosition(atom.GetIdx())
        elements.append(atom.GetSymbol().upper())
        coords.append((point.x, point.y, point.z))
    return elements, np.asarray(coords, dtype=float)


def compare(source_elements, source_xyz, candidate_atoms):
    candidate_elements = [item[0] for item in candidate_atoms]
    if len(source_elements) != len(candidate_elements) or Counter(source_elements) != Counter(candidate_elements):
        return None
    candidate_xyz = np.asarray([item[1] for item in candidate_atoms], dtype=float)
    distances = np.linalg.norm(source_xyz[:, None, :] - candidate_xyz[None, :, :], axis=2)
    for i, source_element in enumerate(source_elements):
        for j, candidate_element in enumerate(candidate_elements):
            if source_element != candidate_element:
                distances[i, j] = 1e6
    row, col = linear_sum_assignment(distances)
    selected = distances[row, col]
    if np.any(selected >= 1e5):
        return None
    return float(np.sqrt(np.mean(selected ** 2))), float(np.max(selected))


def map_one(item, pdb_root: Path, rmsd_cutoff: float, max_cutoff: float):
    pdb_id, block = item
    structure_path = pdb_root / f"pdb{pdb_id}.ent.gz"
    if not structure_path.exists():
        structure_path = pdb_root / pdb_id[1:3] / f"pdb{pdb_id}.ent.gz"
    base = {"pdb_id": pdb_id, "auth_asym_id": None, "auth_seq_id": None,
            "insertion_code": None, "component_id": None, "alt_id": None,
            "coordinate_rmsd_A": None, "max_atom_distance_A": None, "match_count": 0}
    if not structure_path.exists():
        return {**base, "mapping_status": "DEPOSITED_PDB_FILE_MISSING"}
    source = sdf_coordinates(block)
    if source is None:
        return {**base, "mapping_status": "PDBBIND_SDF_COORDINATES_UNREADABLE"}
    elements, xyz = source
    matches = []
    try:
        for chain, resnum, icode, component, alt, atoms in deposited_candidates(structure_path):
            score = compare(elements, xyz, atoms)
            if score is None:
                continue
            rmsd, maximum = score
            if rmsd <= rmsd_cutoff and maximum <= max_cutoff:
                matches.append((rmsd, maximum, chain, resnum, icode, component, alt))
    except Exception:
        return {**base, "mapping_status": "DEPOSITED_PDB_PARSE_FAILED"}
    if not matches:
        return {**base, "mapping_status": "NO_COORDINATE_MATCH"}
    matches.sort()
    if len(matches) != 1:
        return {**base, "mapping_status": "AMBIGUOUS_COORDINATE_MATCH", "match_count": len(matches)}
    rmsd, maximum, chain, resnum, icode, component, alt = matches[0]
    return {
        **base, "auth_asym_id": chain, "auth_seq_id": resnum,
        "insertion_code": icode, "component_id": component, "alt_id": alt,
        "coordinate_rmsd_A": rmsd, "max_atom_distance_A": maximum,
        "match_count": 1, "mapping_status": "EXACT_COORDINATE_MATCH",
    }


def sdf_items(archives):
    for archive in archives:
        with tarfile.open(archive, "r:gz") as tf:
            for member in tf:
                if member.isfile() and member.name.endswith("_ligand.sdf"):
                    pdb_id = Path(member.name).parent.name.lower()
                    block = tf.extractfile(member).read().decode("utf-8", errors="replace")
                    yield pdb_id, block


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--pdb-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--rmsd-cutoff", type=float, default=0.10)
    parser.add_argument("--max-cutoff", type=float, default=0.20)
    parser.add_argument("--batch-size", type=int, default=500)
    args = parser.parse_args()
    root = args.external_root / "pdbbind_v2020"
    archives = [root / "PDBbind_v2020_refined.tar.gz", root / "PDBbind_v2020_other_PL.tar.gz"]
    rows, batch = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for item in sdf_items(archives):
            batch.append(item)
            if len(batch) < args.batch_size:
                continue
            rows.extend(pool.map(lambda value: map_one(value, args.pdb_root, args.rmsd_cutoff, args.max_cutoff), batch))
            print(f"mapped={len(rows)}", flush=True)
            batch = []
        if batch:
            rows.extend(pool.map(lambda value: map_one(value, args.pdb_root, args.rmsd_cutoff, args.max_cutoff), batch))
    frame = pd.DataFrame(rows)
    if len(frame) != 19443 or frame["pdb_id"].nunique() != 19443:
        raise RuntimeError(f"PDBbind mapping population mismatch: {len(frame)} rows")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(args.output, index=False)
    print(frame["mapping_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
