#!/usr/bin/env python3
"""Processing 4 v2: frozen docking-ready case construction.

This stage is deliberately a serializer/preparation stage. All case, receptor,
ligand, and pocket identities come from frozen upstream tables.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import sqlite3
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import gemmi
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = Path(
    os.environ.get("PROCESSING4_BENCHMARK_ROOT")
    or os.environ.get("BENCHMARK_DATA_ROOT")
    or Path.cwd()
).expanduser().resolve()
P2_COORD_RUN = BENCHMARK / "processing_2_assembly_ready_structure_preparation/runs/20260810_full_01"
P2_FORMAL_RUN = BENCHMARK / "processing_2_assembly_ready_structure_preparation/runs/20260826_validation_provenance_required_01"
P3_DETAIL_RUN = BENCHMARK / "processing_03_direct_contact_qualification/runs/20260811_full_01"
P3_FORMAL_RUN = BENCHMARK / "processing_03_direct_contact_qualification/runs/20260826_p2_validation_provenance_01"
F4_RUN = BENCHMARK / "filter_04_crystal_packing_influence/database_runs/20260826_filter3_118255_strict_posebusters_01"
F5_RUN = BENCHMARK / "filter_05_benchmark_redundancy_reduction/runs/20260828_exact_only_v3_01"
F5_RETAINED = F5_RUN / "final/02_benchmark_filter5_v3_retained_cases.parquet"
OLD_P4_RUN = BENCHMARK / (
    "processing_04_docking_ready_case_construction/runs/"
    "p4_benchmark_f5v2_30187_gt3_v3_0_1_01"
)
CCD_IDEAL_CACHE = BENCHMARK / "processing_04_docking_ready_case_construction/references/ccd_ideal_heavy_cache_v20260711_v2.sqlite"
DEFAULT_CONFIG = REPOSITORY_ROOT / "pipeline/processing_04/config/processing4_benchmark_f5v3_exact_only.json"
EXPECTED_CASES = 65162
CHAIN_CODES = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
READY_FILES = (
    "receptor.cif", "receptor.pdb", "ligand_reference.sdf", "ligand.smi",
    "ligand_start.sdf", "site.json", "metadata.json",
)


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), tmp, compression="zstd")
    os.replace(tmp, path)


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def stable_bucket(pdb_id: str, count: int = 256) -> int:
    return int(hashlib.sha256(pdb_id.lower().encode()).hexdigest()[:8], 16) % count


def case_id(pair_id: str, pdb_id: str) -> str:
    return f"P4_{pdb_id.lower()}_{hashlib.sha256(pair_id.encode()).hexdigest()[:20]}"


def read_config(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "policy_version", "schema_version", "site_definition_version",
        "etkdg_random_seed", "uff_max_iterations", "bucket_count",
    }
    missing = required - payload.keys()
    if missing:
        raise ValueError(f"config missing keys: {sorted(missing)}")
    return payload


def dataset_table(path: Path, columns: list[str] | None = None) -> pa.Table:
    if not path.exists():
        raise FileNotFoundError(path)
    return ds.dataset(str(path), format="parquet").to_table(columns=columns)


def prepare_run(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    if (run / "_FROZEN.json").exists():
        raise RuntimeError("refusing to modify frozen run")
    config_path = Path(args.config).resolve()
    config = read_config(config_path)
    run.mkdir(parents=True, exist_ok=True)
    (run / "input").mkdir(exist_ok=True)
    (run / "output/cases").mkdir(parents=True, exist_ok=True)
    (run / "work/buckets").mkdir(parents=True, exist_ok=True)
    (run / "logs").mkdir(exist_ok=True)

    retained_full = pq.read_table(F5_RETAINED).to_pandas()
    retained = retained_full[["pair_id", "filter3_quality_class", "normalized_ccd_isomeric_smiles",
                              "filter5_v3_final_status"]].copy()
    retained = retained.rename(columns={"normalized_ccd_isomeric_smiles": "ligand_exact_smiles",
                                        "filter5_v3_final_status": "filter5_final_status"}).fillna("")
    if len(retained) != EXPECTED_CASES or not retained["pair_id"].is_unique:
        raise RuntimeError(f"unexpected Filter 5 retained inventory: {len(retained)}")
    if args.mode == "smoke":
        retained["_score"] = retained["pair_id"].map(
            lambda x: hashlib.sha256(("processing4-smoke-v1|" + x).encode()).hexdigest()
        )
        retained = retained.sort_values("_score").head(args.limit).drop(columns="_score")
    elif args.mode != "full":
        raise ValueError(args.mode)

    pairs = pq.read_table(P3_FORMAL_RUN / "output/ordinary_noncovalent_pairs_with_plip.parquet").to_pandas()
    pairs = pairs[pairs["pair_id"].isin(set(retained["pair_id"]))].copy()
    merged = retained.merge(pairs, on="pair_id", how="left", validate="one_to_one")
    missing = merged["ligand_assembly_placement_id"].isna()
    if missing.any():
        raise RuntimeError(f"{int(missing.sum())} retained pairs missing from frozen P3")
    merged["bucket_id"] = merged["pdb_id"].map(lambda x: stable_bucket(str(x), config["bucket_count"]))
    merged["case_id"] = [case_id(p, q) for p, q in zip(merged["pair_id"], merged["pdb_id"])]
    if not merged["case_id"].is_unique:
        raise RuntimeError("case_id collision")
    merged = merged.sort_values(["bucket_id", "case_id"]).reset_index(drop=True)
    atomic_parquet(run / "input/case_inventory.parquet", merged)
    shutil.copy2(config_path, run / "input/config_snapshot.json")
    provenance = {
        "stage": "processing4_benchmark_v3.1.0_filter5_v3_exact_only", "created_at": utc(), "mode": args.mode,
        "expected_case_count": int(len(merged)), "filter5_retained_count": EXPECTED_CASES,
        "inputs": {
            "filter5_retained": str(F5_RETAINED),
            "filter5_retained_sha256": sha256_file(F5_RETAINED),
            "filter5_frozen_marker": str(F5_RUN / "_FROZEN.json"),
            "filter4_frozen_marker": str(F4_RUN / "_FROZEN.json"),
            "processing3_formal_frozen_marker": str(P3_FORMAL_RUN / "_FROZEN.json"),
            "processing2_formal_frozen_marker": str(P2_FORMAL_RUN / "_FROZEN.json"),
            "processing3_detail_coordinate_source": str(P3_DETAIL_RUN / "_FROZEN.json"),
            "processing2_coordinate_source": str(P2_COORD_RUN / "_FROZEN.json"),
            "ccd_ideal_cache": str(CCD_IDEAL_CACHE),
            "ccd_ideal_cache_sha256": sha256_file(CCD_IDEAL_CACHE),
        },
        "versions": {"python": sys.version, "rdkit": rdkit.__version__,
                     "pyarrow": pa.__version__, "gemmi": gemmi.__version__},
        "prohibited_reruns": ["PLIP", "Arpeggio", "pocket_detection", "crystal_packing", "Filter5_deduplication"],
        "native_coordinate_fallback": "FORBIDDEN",
        "membership_rule": "ligand heavy_atom_count > 3; preparation inconsistencies are warnings, not scientific filters",
    }
    reused = 0
    regenerated_for_heavy_ideal = 0
    if OLD_P4_RUN.exists():
        old_status = pq.read_table(OLD_P4_RUN / "output/processing4_case_inventory.parquet").to_pandas()
        old_status = old_status[(old_status.status == "P4_DOCKING_READY") &
                                (pd.to_numeric(old_status.ligand_heavy_atoms, errors="coerce") > 3)]
        old_by_pair = old_status.set_index("pair_id")
        ideal_connection = sqlite3.connect(f"file:{CCD_IDEAL_CACHE}?mode=ro", uri=True)
        heavy_ideal_components = {str(x[0]) for x in ideal_connection.execute(
            "select distinct component_id from ideal_atoms"
        )}
        ideal_connection.close()
        for row in merged.itertuples(index=False):
            if row.pair_id not in old_by_pair.index:
                continue
            old = old_by_pair.loc[row.pair_id]
            if str(old.ligand_start_generation_method) != "CCD_IDEAL" and str(row.component_id) in heavy_ideal_components:
                regenerated_for_heavy_ideal += 1
                continue
            src = OLD_P4_RUN / f"output/cases/bucket_{int(row.bucket_id):03d}" / str(old.case_id)
            dst = run / f"output/cases/bucket_{int(row.bucket_id):03d}" / str(row.case_id)
            if not (src / "_SUCCESS.json").is_file() or dst.exists():
                continue
            dst.mkdir(parents=True)
            for source in src.iterdir():
                if source.is_file():
                    os.link(source, dst / source.name)
            reused += 1
    provenance["reused_prior_ready_cases"] = reused
    provenance["regenerated_for_new_heavy_atom_ideal_policy"] = regenerated_for_heavy_ideal
    provenance["reuse_source"] = str(OLD_P4_RUN)
    atomic_json(run / "input/provenance.json", provenance)
    atomic_json(run / "status.json", {"status": "PREPARED", "mode": args.mode,
                                      "case_count": len(merged), "created_at": utc()})
    print(json.dumps({"run_dir": str(run), "cases": len(merged),
                      "buckets": int(merged["bucket_id"].nunique()), "reused_ready": reused,
                      "regenerated_for_heavy_ideal": regenerated_for_heavy_ideal}, indent=2))


class CCDStore:
    def __init__(self, path: Path, ideal_cache: Path):
        self.conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        self.ideal_conn = sqlite3.connect(f"file:{ideal_cache}?mode=ro", uri=True)
        self.cache: dict[str, dict[str, Any]] = {}

    def component(self, component_id: str) -> dict[str, Any]:
        if component_id in self.cache:
            return self.cache[component_id]
        comp = self.conn.execute(
            "select formal_charge,descriptor_smiles_canonical,descriptor_smiles,descriptor_inchikey "
            "from components where component_id=?", (component_id,)
        ).fetchone()
        if comp is None:
            raise ValueError(f"CCD component missing: {component_id}")
        atoms = self.conn.execute(
            "select atom_id,element,charge,aromatic_flag,stereo_config,ordinal from atoms "
            "where component_id=? order by ordinal", (component_id,)
        ).fetchall()
        bonds = self.conn.execute(
            "select atom_id_1,atom_id_2,bond_order,aromatic_flag,stereo_config,ordinal from bonds "
            "where component_id=? order by ordinal", (component_id,)
        ).fetchall()
        out = {"formal_charge": int(comp[0]), "descriptor": clean_descriptor(comp[1] or comp[2]),
               "inchikey": comp[3] or "", "atoms": atoms, "bonds": bonds}
        self.cache[component_id] = out
        return out

    def ideal_coordinates(self, component_id: str) -> dict[str, tuple[str, float, float, float]]:
        rows = self.ideal_conn.execute(
            "select atom_id,element,x,y,z from ideal_atoms where component_id=? order by ordinal",
            (component_id,),
        ).fetchall()
        return {
            str(atom_id): (str(element), float(x), float(y), float(z))
            for atom_id, element, x, y, z in rows
        }


def clean_descriptor(value: str) -> str:
    value = (value or "").strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        try:
            return json.loads(value)
        except Exception:
            return value[1:-1]
    return value


def connectivity_isomorphism(candidate: Chem.Mol, mapped: Chem.Mol) -> tuple[int, ...] | None:
    """Map candidate atoms to observed CCD atoms using element+connectivity only."""
    params = Chem.AdjustQueryParameters.NoAdjustments()
    params.makeBondsGeneric = True
    query = Chem.AdjustQueryProperties(Chem.Mol(candidate), params)
    match = mapped.GetSubstructMatch(query, useChirality=False)
    if not match or len(match) != candidate.GetNumAtoms():
        return None
    return tuple(int(x) for x in match)


def build_frozen_ligand(atom_rows: pd.DataFrame, bond_rows: pd.DataFrame,
                        ccd: dict[str, Any], authoritative_smiles: str = "") -> tuple[Chem.Mol, dict[str, Any]]:
    coord_by_name = {
        str(r.label_atom_id): (str(r.type_symbol), float(r.Cartn_x), float(r.Cartn_y), float(r.Cartn_z))
        for r in atom_rows.itertuples()
    }
    observed = set(coord_by_name)
    observed_heavy = {name for name, values in coord_by_name.items() if values[0].upper() != "H"}
    ccd_atoms = [r for r in ccd["atoms"] if r[0] in observed_heavy and str(r[1]).upper() != "H"]
    rw = Chem.RWMol()
    ccd_index: dict[str, int] = {}
    aromatic_names: set[str] = set()
    for name, element, charge, aromatic, _stereo, _ordinal in ccd_atoms:
        atom = Chem.Atom(str(element).title())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("_CCD_ATOM_ID", str(name))
        if str(aromatic).upper() == "Y":
            atom.SetIsAromatic(True)
            aromatic_names.add(str(name))
        ccd_index[str(name)] = rw.AddAtom(atom)
    bond_types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                  "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for r in bond_rows.sort_values("bond_index").itertuples():
        a, b = str(r.atom_id_1), str(r.atom_id_2)
        if a not in ccd_index or b not in ccd_index:
            continue
        aromatic = str(r.aromatic_flag).upper() == "Y"
        bt = Chem.BondType.AROMATIC if aromatic else bond_types.get(str(r.bond_order).upper())
        if bt is None:
            raise ValueError(f"unsupported bond order {r.bond_order}")
        rw.AddBond(ccd_index[a], ccd_index[b], bt)
        if aromatic:
            rw.GetAtomWithIdx(ccd_index[a]).SetIsAromatic(True)
            rw.GetAtomWithIdx(ccd_index[b]).SetIsAromatic(True)
    mapped = rw.GetMol()
    mapped.UpdatePropertyCache(strict=False)
    representation_warnings = []
    try:
        Chem.SanitizeMol(mapped)
    except Exception as exc:
        representation_warnings.append(f"OBSERVED_GRAPH_SANITIZE_WARNING[{type(exc).__name__}: {exc}]")

    candidates = []
    if authoritative_smiles:
        candidates.append(("FILTER5_EXACT_SMILES", authoritative_smiles))
    if ccd["descriptor"] and ccd["descriptor"] != authoritative_smiles:
        candidates.append(("CCD_DESCRIPTOR", ccd["descriptor"]))
    mol = None
    match = None
    graph_source = ""
    for source, smiles in candidates:
        candidate = Chem.MolFromSmiles(smiles)
        if candidate is None:
            representation_warnings.append(f"{source}_PARSE_WARNING")
            continue
        candidate = Chem.RemoveHs(candidate)
        if candidate.GetNumAtoms() != mapped.GetNumAtoms():
            representation_warnings.append(f"{source}_ATOM_COUNT_WARNING")
            continue
        candidate_match = connectivity_isomorphism(candidate, mapped)
        if candidate_match is None:
            representation_warnings.append(f"{source}_CONNECTIVITY_MAPPING_WARNING")
            continue
        mol, match, graph_source = Chem.Mol(candidate), candidate_match, source
        break
    if mol is None or match is None:
        # Last-resort serialization uses the upstream frozen observed graph itself.
        # This is not a membership gate; any representation limitation is recorded.
        try:
            Chem.SanitizeMol(mapped, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            mol, match, graph_source = Chem.Mol(mapped), tuple(range(mapped.GetNumAtoms())), "UPSTREAM_FROZEN_GRAPH"
            representation_warnings.append("DESCRIPTOR_MAPPING_UNAVAILABLE_USED_UPSTREAM_FROZEN_GRAPH")
        except Exception as exc:
            raise ValueError(f"upstream frozen graph cannot be serialized: {type(exc).__name__}: {exc}")
    conf = Chem.Conformer(mol.GetNumAtoms())
    atom_map: list[str] = []
    for qidx, midx in enumerate(match):
        name = mapped.GetAtomWithIdx(midx).GetProp("_CCD_ATOM_ID")
        element, x, y, z = coord_by_name[name]
        if mol.GetAtomWithIdx(qidx).GetSymbol().upper() != element.upper():
            raise ValueError(f"element mismatch at {name}")
        mol.GetAtomWithIdx(qidx).SetProp("_CCD_ATOM_ID", name)
        conf.SetAtomPosition(qidx, (x, y, z))
        atom_map.append(name)
    mol.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if Chem.GetFormalCharge(mol) != int(ccd["formal_charge"]):
        raise ValueError("formal charge mismatch against CCD")
    if len(atom_map) != len(observed_heavy):
        raise ValueError("heavy atom mapping does not close")
    return mol, {"ccd_atom_order": atom_map, "ccd_inchikey": ccd["inchikey"],
                 "graph_source": graph_source, "representation_warnings": representation_warnings}


def graph_key(mol: Chem.Mol, stereo_source: Chem.Mol | None = None) -> str:
    copy = Chem.RemoveHs(Chem.Mol(mol))
    copy.RemoveAllConformers()
    if stereo_source is not None:
        source = Chem.RemoveHs(Chem.Mol(stereo_source))
        source.RemoveAllConformers()
        match = copy.GetSubstructMatch(source, useChirality=False)
        if not match or len(match) != source.GetNumAtoms():
            raise RuntimeError("candidate graph cannot be mapped to frozen stereo source")
        for source_idx, candidate_idx in enumerate(match):
            source_atom = source.GetAtomWithIdx(source_idx)
            if source_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                copy.GetAtomWithIdx(candidate_idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for source_bond in source.GetBonds():
            a = match[source_bond.GetBeginAtomIdx()]
            b = match[source_bond.GetEndAtomIdx()]
            candidate_bond = copy.GetBondBetweenAtoms(a, b)
            if source_bond.GetStereo() in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
                candidate_bond.SetStereo(Chem.BondStereo.STEREONONE)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=True)


def write_sdf(path: Path, mol: Chem.Mol, case: str, role: str) -> None:
    obj = Chem.Mol(mol)
    obj.SetProp("_Name", case)
    obj.SetProp("PROCESSING4_ROLE", role)
    writer = Chem.SDWriter(str(path))
    writer.write(obj)
    writer.close()
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"empty SDF: {path}")


def independent_start(graph_mol: Chem.Mol, config: dict[str, Any]) -> tuple[Chem.Mol, int | None, str]:
    base = Chem.Mol(graph_mol)
    base.RemoveAllConformers()
    if base.GetNumConformers() != 0:
        raise RuntimeError("coordinate removal failed")
    with_h = Chem.AddHs(base)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(config["etkdg_random_seed"])
    params.enforceChirality = bool(config["etkdg_enforce_chirality"])
    params.numThreads = 1
    code = AllChem.EmbedMolecule(with_h, params)
    method = "ETKDGv3_UFF_FALLBACK"
    if code != 0:
        params.useRandomCoords = True
        params.randomSeed = int(config["etkdg_random_seed"])
        code = AllChem.EmbedMolecule(with_h, params)
        method = "ETKDGv3_RANDOM_COORDS_FALLBACK"
    if code != 0:
        raise RuntimeError(f"ETKDGv3 deterministic fallbacks returned {code}")
    uff_code = None
    if AllChem.UFFHasAllMoleculeParams(with_h):
        uff_code = AllChem.UFFOptimizeMolecule(with_h, maxIters=int(config["uff_max_iterations"]))
    else:
        method += "_NO_UFF"
    start = Chem.RemoveHs(with_h)
    if start.GetNumAtoms() != base.GetNumAtoms():
        raise RuntimeError("RemoveHs changed heavy atom count")
    return start, None if uff_code is None else int(uff_code), method


def ccd_ideal_start(
    graph_mol: Chem.Mol,
    component_id: str,
    ccd_store: CCDStore,
) -> Chem.Mol:
    """Put frozen-graph atoms at wwPDB CCD ideal coordinates by CCD atom_id."""
    ideal = ccd_store.ideal_coordinates(component_id)
    if not ideal:
        raise RuntimeError("CCD ideal coordinates unavailable")
    start = Chem.Mol(graph_mol)
    start.RemoveAllConformers()
    conf = Chem.Conformer(start.GetNumAtoms())
    xyz: list[tuple[float, float, float]] = []
    for index, atom in enumerate(start.GetAtoms()):
        if not atom.HasProp("_CCD_ATOM_ID"):
            raise RuntimeError(f"frozen atom {index} lacks CCD atom_id")
        atom_id = atom.GetProp("_CCD_ATOM_ID")
        if atom_id not in ideal:
            raise RuntimeError(f"CCD ideal atom mapping missing: {atom_id}")
        element, x, y, z = ideal[atom_id]
        if atom.GetSymbol().upper() != element.upper():
            raise RuntimeError(f"CCD ideal element mismatch: {atom_id}")
        if not all(math.isfinite(v) for v in (x, y, z)):
            raise RuntimeError(f"CCD ideal non-finite coordinate: {atom_id}")
        conf.SetAtomPosition(index, (x, y, z))
        xyz.append((x, y, z))
    arr = np.asarray(xyz, dtype=float)
    if len(arr) != start.GetNumHeavyAtoms() or len(arr) < 3:
        raise RuntimeError("CCD ideal heavy-atom mapping does not close")
    delta = arr[:, None, :] - arr[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances += np.eye(len(arr)) * 999.0
    if float(distances.min()) <= 0.10:
        raise RuntimeError("CCD ideal geometry contains overlapping heavy atoms")
    for bond in start.GetBonds():
        length = float(np.linalg.norm(arr[bond.GetBeginAtomIdx()] - arr[bond.GetEndAtomIdx()]))
        if not math.isfinite(length) or length <= 0.20 or length >= 6.0:
            raise RuntimeError(f"CCD ideal bonded geometry invalid: {length}")
    start.AddConformer(conf, assignId=True)
    return start


def parse_int(value: Any, fallback: int) -> int:
    text = str(value).strip()
    try:
        return int(float(text))
    except Exception:
        return fallback


def normalized_icode(value: Any) -> str:
    text = str(value or "").strip()
    if text in {"", ".", "?", "False", "None", "nan"}:
        return " "
    return text[0]


def build_receptor_structure(rec_rows: pd.DataFrame, chain_ids: list[str]) -> tuple[gemmi.Structure, dict[str, str]]:
    if len(chain_ids) > len(CHAIN_CODES):
        raise ValueError("receptor has more than 62 chain instances; PDB compatibility unavailable")
    structure = gemmi.Structure()
    structure.name = "processing4_receptor"
    model = gemmi.Model("1")
    chain_map: dict[str, str] = {}
    serial = 1
    for chain_pos, chain_instance in enumerate(chain_ids):
        rows = rec_rows[rec_rows["filter_1_chain_instance_id"].astype(str) == chain_instance].copy()
        if rows.empty:
            raise ValueError(f"missing receptor chain instance {chain_instance}")
        out_chain = CHAIN_CODES[chain_pos]
        chain_map[chain_instance] = out_chain
        chain = gemmi.Chain(out_chain)
        residue_cols = ["label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id"]
        for res_pos, (_key, group) in enumerate(rows.groupby(residue_cols, sort=False, dropna=False), start=1):
            first = group.iloc[0]
            residue = gemmi.Residue()
            residue.name = str(first["label_comp_id"] or first["auth_comp_id"])
            seqnum = parse_int(first["auth_seq_id"], parse_int(first["label_seq_id"], res_pos))
            residue.seqid = gemmi.SeqId(seqnum, normalized_icode(first["insertion_code"]))
            label_seq = parse_int(first["label_seq_id"], res_pos)
            if label_seq > 0:
                residue.label_seq = label_seq
            for row in group.itertuples():
                atom = gemmi.Atom()
                atom.name = str(row.label_atom_id)
                atom.element = gemmi.Element(str(row.type_symbol).title())
                atom.pos = gemmi.Position(float(row.Cartn_x), float(row.Cartn_y), float(row.Cartn_z))
                atom.occ = float(row.occupancy) if pd.notna(row.occupancy) else 1.0
                atom.b_iso = float(row.B_iso_or_equiv) if pd.notna(row.B_iso_or_equiv) else 0.0
                atom.serial = serial
                serial += 1
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    return structure, chain_map


def structure_heavy_coordinates(path: Path) -> np.ndarray:
    st = gemmi.read_structure(str(path))
    coords = []
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element.name.upper() != "H":
                        coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
    return np.asarray(coords, dtype=float)


def receptor_source_coordinates(rec_rows: pd.DataFrame, chain_ids: list[str]) -> np.ndarray:
    coords = []
    for cid in chain_ids:
        rows = rec_rows[rec_rows["filter_1_chain_instance_id"].astype(str) == cid]
        for r in rows.itertuples():
            if str(r.type_symbol).upper() != "H":
                coords.append((float(r.Cartn_x), float(r.Cartn_y), float(r.Cartn_z)))
    return np.asarray(coords, dtype=float)


def native_site(atom_rows: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    heavy = atom_rows[atom_rows["type_symbol"].astype(str).str.upper() != "H"]
    xyz = heavy[["Cartn_x", "Cartn_y", "Cartn_z"]].astype(float).to_numpy()
    if len(xyz) == 0 or not np.isfinite(xyz).all():
        raise ValueError("native ligand heavy coordinates missing/non-finite")
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    center = (lo + hi) / 2.0
    extent = hi - lo
    margin = float(config["site_margin_angstrom_per_side"])
    minimum = float(config["site_minimum_dimension_angstrom"])
    size = np.maximum(extent + 2.0 * margin, minimum)
    if not np.all(size > 0) or np.any(lo < center - size / 2.0 - 1e-8) or np.any(hi > center + size / 2.0 + 1e-8):
        raise RuntimeError("site box does not cover native ligand")
    return {
        "site_center": {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])},
        "search_box": {"size_x": float(size[0]), "size_y": float(size[1]), "size_z": float(size[2])},
        "definition": "native_ligand_based",
        "center_policy": config["site_center"],
        "margin_angstrom_per_side": margin,
        "minimum_dimension_angstrom": minimum,
        "site_definition_version": config["site_definition_version"],
    }


def load_sdf_one(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True)
    mols = [m for m in supplier if m is not None]
    if len(mols) != 1:
        raise RuntimeError(f"expected one parseable molecule in {path.name}")
    return mols[0]


def create_case(tmp: Path, case_row: pd.Series, rec_rows: pd.DataFrame,
                lig_atoms: pd.DataFrame, lig_bonds: pd.DataFrame,
                pocket: pd.DataFrame, binding: pd.DataFrame,
                ccd_store: CCDStore, config: dict[str, Any]) -> dict[str, Any]:
    pair_id = str(case_row["pair_id"])
    cid = str(case_row["case_id"])
    chain_ids = [x for x in str(case_row["receptor_chain_instance_ids"]).split(",") if x]
    if len(chain_ids) != int(case_row["receptor_chain_count"]):
        raise ValueError("receptor chain count does not close")
    placement = str(case_row["ligand_assembly_placement_id"])
    component = str(case_row["component_id"])
    if lig_atoms.empty:
        raise ValueError("frozen ligand atoms missing")
    if set(lig_atoms["filter_2_ligand_assembly_placement_id"].astype(str)) != {placement}:
        raise ValueError("ligand placement reconstruction mismatch")
    ligand_heavy_atom_count = int(
        lig_atoms.loc[lig_atoms["type_symbol"].astype(str).str.upper() != "H", "label_atom_id"].nunique()
    )
    if ligand_heavy_atom_count < int(config["minimum_ligand_heavy_atom_count"]):
        raise TooSmallForPoseDocking(ligand_heavy_atom_count)
    if lig_bonds.empty:
        raise ValueError("frozen ligand bonds missing for multi-atom ligand")
    present_chains = set(rec_rows["filter_1_chain_instance_id"].astype(str))
    if present_chains != set(chain_ids):
        raise ValueError("receptor chain instance reconstruction mismatch")

    structure, chain_map = build_receptor_structure(rec_rows, chain_ids)
    structure.make_mmcif_document().write_file(str(tmp / "receptor.cif"))
    structure.write_pdb(str(tmp / "receptor.pdb"))
    source_xyz = receptor_source_coordinates(rec_rows, chain_ids)
    cif_xyz = structure_heavy_coordinates(tmp / "receptor.cif")
    pdb_xyz = structure_heavy_coordinates(tmp / "receptor.pdb")
    if source_xyz.shape != cif_xyz.shape or source_xyz.shape != pdb_xyz.shape:
        raise RuntimeError("receptor parsed heavy atom count mismatch")
    cif_delta = float(np.max(np.abs(source_xyz - cif_xyz))) if len(source_xyz) else 0.0
    pdb_delta = float(np.max(np.abs(source_xyz - pdb_xyz))) if len(source_xyz) else 0.0
    if cif_delta > float(config["coordinate_tolerance_cif_angstrom"]):
        raise RuntimeError(f"CIF receptor coordinate delta {cif_delta}")
    if pdb_delta > float(config["coordinate_tolerance_pdb_angstrom"]):
        raise RuntimeError(f"PDB receptor coordinate delta {pdb_delta}")

    mol_ref, ligand_info = build_frozen_ligand(
        lig_atoms, lig_bonds, ccd_store.component(component), str(case_row.get("ligand_exact_smiles", ""))
    )
    if mol_ref.GetNumHeavyAtoms() != ligand_heavy_atom_count:
        raise RuntimeError("frozen graph heavy-atom count disagrees with prepared coordinates")
    frozen_key = graph_key(mol_ref)
    frozen_stereo_source = Chem.MolFromSmiles(frozen_key)
    if frozen_stereo_source is None:
        raise RuntimeError("frozen canonical isomeric SMILES cannot be parsed")
    write_sdf(tmp / "ligand_reference.sdf", mol_ref, cid, "experimental_native_pose")
    (tmp / "ligand.smi").write_text(f"{frozen_key}\t{cid}\n", encoding="utf-8")
    graph_only = Chem.Mol(mol_ref)
    graph_only.RemoveAllConformers()
    ideal_error = ""
    uff_code: int | None = None
    try:
        mol_start = ccd_ideal_start(graph_only, component, ccd_store)
        start_method = "CCD_IDEAL"
        start_role = "wwpdb_ccd_ideal_coordinates_on_frozen_graph"
    except Exception as exc:
        ideal_error = f"{type(exc).__name__}: {exc}"
        try:
            mol_start, uff_code, start_method = independent_start(graph_only, config)
            start_role = "independent_etkdgv3_uff_fallback_start"
        except Exception as fallback_exc:
            raise LigandStartFailure(
                f"CCD_IDEAL_FAILED[{ideal_error}]; ETKDG_UFF_FAILED[{type(fallback_exc).__name__}: {fallback_exc}]"
            ) from fallback_exc
    write_sdf(tmp / "ligand_start.sdf", mol_start, cid, start_role)

    ref_parsed = load_sdf_one(tmp / "ligand_reference.sdf")
    start_parsed = load_sdf_one(tmp / "ligand_start.sdf")
    smi_text = (tmp / "ligand.smi").read_text(encoding="utf-8").split("\t", 1)[0]
    smi_parsed = Chem.MolFromSmiles(smi_text)
    if smi_parsed is None:
        raise RuntimeError("ligand.smi parse failed")
    keys = {
        graph_key(ref_parsed, frozen_stereo_source),
        graph_key(start_parsed, frozen_stereo_source),
        graph_key(smi_parsed, frozen_stereo_source),
        frozen_key,
    }
    graph_validation_warning = ""
    if len(keys) != 1:
        graph_validation_warning = f"SERIALIZATION_GRAPH_REPRESENTATION_WARNING[{sorted(keys)}]"
    input_native = np.asarray([
        (float(r.Cartn_x), float(r.Cartn_y), float(r.Cartn_z))
        for r in lig_atoms.set_index("label_atom_id").loc[ligand_info["ccd_atom_order"]].itertuples()
    ])
    ref_conf = ref_parsed.GetConformer()
    parsed_native = np.asarray([tuple(ref_conf.GetAtomPosition(i)) for i in range(ref_parsed.GetNumAtoms())])
    native_delta = float(np.max(np.abs(input_native - parsed_native)))
    if native_delta > 0.00011:
        raise RuntimeError(f"ligand reference coordinate delta {native_delta}")

    site = native_site(lig_atoms, config)
    site["case_id"] = cid
    site["pair_id"] = pair_id
    site["pocket_residue_ids"] = sorted(set(
        str(r.chain_instance_id) + "|" + str(r.protein_residue_id) for r in pocket.itertuples()
    ))
    binding_ids = sorted(set(
        str(r.chain_instance_id) + "|" + str(r.protein_residue_id) for r in binding.itertuples()
    ))
    site["binding_residue_ids"] = binding_ids
    atomic_json(tmp / "site.json", site)
    metadata = {
        "case_id": cid, "pair_id": pair_id, "pdb_id": str(case_row["pdb_id"]),
        "assembly_id": str(case_row["assembly_id"]), "model_id": str(case_row["model_id"]),
        "receptor_chain_instances": chain_ids, "receptor_chain_file_mapping": chain_map,
        "receptor_heavy_atom_count": int(len(source_xyz)),
        "receptor_coordinate_max_delta_cif_angstrom": cif_delta,
        "receptor_coordinate_max_delta_pdb_angstrom": pdb_delta,
        "ligand": {"placement_id": placement, "component_id": component,
                   "formal_charge": int(Chem.GetFormalCharge(mol_ref)),
                   "heavy_atom_count": mol_ref.GetNumHeavyAtoms(), "canonical_isomeric_smiles": frozen_key,
                   "experimental_coordinate_max_delta_angstrom": native_delta, **ligand_info},
        "ligand_start_generation_method": start_method,
        "ligand_start_generation_primary_error": ideal_error,
        "ligand_start_coordinate_source": (
            "frozen_wwPDB_components.cif_ideal_coordinates"
            if start_method == "CCD_IDEAL" else "coordinate_free_frozen_graph_ETKDGv3_UFF"
        ),
        "binding_residue_ids": binding_ids,
        "pocket_residue_ids": site["pocket_residue_ids"],
        "receptor_format_version": config["receptor_format_version"],
        "ligand_format_version": config["ligand_format_version"],
        "rdkit_version": rdkit.__version__, "ETKDG_version": config["etkdg_version"],
        "ETKDG_seed": config["etkdg_random_seed"], "ETKDG_enforce_chirality": config["etkdg_enforce_chirality"],
        "UFF_parameters": {"max_iterations": config["uff_max_iterations"], "return_code": uff_code},
        "site_definition_version": config["site_definition_version"],
        "upstream_filter5_status": str(case_row["filter5_final_status"]),
        "upstream_filter3_quality_class": str(case_row["filter3_quality_class"]),
        "upstream_filter4_status": "PASS",
        "processing4_status": "P4_DOCKING_READY",
        "native_pose_leakage_control": "CCD ideal coordinates or coordinate-free graph only; native-coordinate fallback forbidden",
        "representation_warnings": ligand_info.get("representation_warnings", []) +
                                   ([graph_validation_warning] if graph_validation_warning else []),
        "created_at": utc(),
    }
    atomic_json(tmp / "metadata.json", metadata)
    for filename in READY_FILES:
        if not (tmp / filename).is_file() or (tmp / filename).stat().st_size == 0:
            raise RuntimeError(f"missing/empty required file {filename}")
    sums = {f: sha256_file(tmp / f) for f in READY_FILES}
    atomic_json(tmp / "_SUCCESS.json", {"status": "P4_DOCKING_READY", "sha256": sums, "created_at": utc()})
    return {"case_id": cid, "pair_id": pair_id, "status": "P4_DOCKING_READY", "reason": "",
            "receptor_atoms": int(len(source_xyz)), "ligand_heavy_atoms": int(mol_ref.GetNumHeavyAtoms()),
            "ligand_start_generation_method": start_method}


class LigandStartFailure(RuntimeError):
    pass


class TooSmallForPoseDocking(RuntimeError):
    def __init__(self, heavy_atom_count: int):
        self.heavy_atom_count = int(heavy_atom_count)
        super().__init__(f"ligand heavy_atom_count={heavy_atom_count} is below pose-docking minimum")


def bucket_paths(run: Path, bid: int) -> tuple[Path, Path]:
    return run / f"output/cases/bucket_{bid:03d}", run / f"work/buckets/bucket_{bid:03d}.parquet"


def read_partition(root: Path, bid: int) -> pd.DataFrame:
    path = root / f"bucket_id={bid:03d}"
    if not path.exists():
        return pd.DataFrame()
    return dataset_table(path).to_pandas()


def run_bucket(run_dir: str, bid: int, force: bool = False,
               retry_nonready: bool = False) -> dict[str, Any]:
    started = time.time()
    run = Path(run_dir)
    config = read_config(run / "input/config_snapshot.json")
    inventory = pq.read_table(run / "input/case_inventory.parquet",
                              filters=[("bucket_id", "=", bid)]).to_pandas()
    if inventory.empty:
        return {"bucket_id": bid, "cases": 0, "ready": 0, "review": 0, "seconds": 0.0}
    case_root, status_path = bucket_paths(run, bid)
    case_root.mkdir(parents=True, exist_ok=True)

    p2out = P2_COORD_RUN / "output"
    p3out = P3_DETAIL_RUN / "output"
    lig_atoms = read_partition(p2out / "prepared_ligand_assembly_atoms", bid)
    lig_bonds = read_partition(p2out / "prepared_ligand_assembly_bonds", bid)
    rec_atoms = read_partition(p2out / "prepared_receptor_assembly_atoms", bid)
    pocket = read_partition(p3out / "pair_pocket_residues", bid)
    binding = read_partition(p3out / "binding_residues", bid)
    ccd = CCDStore(P2_COORD_RUN / "input/ccd_active_snapshot.sqlite", CCD_IDEAL_CACHE)
    results: list[dict[str, Any]] = []
    for _, row in inventory.iterrows():
        cid = str(row["case_id"])
        target = case_root / cid
        if (target / "_SUCCESS.json").exists() and not force:
            prior_metadata = json.loads((target / "metadata.json").read_text(encoding="utf-8"))
            results.append({"case_id": cid, "pair_id": row["pair_id"], "status": "P4_DOCKING_READY",
                            "reason": "REUSED_VALIDATED_PRIOR_SUCCESS",
                            "receptor_atoms": prior_metadata.get("receptor_heavy_atom_count"),
                            "ligand_heavy_atoms": prior_metadata.get("ligand", {}).get("heavy_atom_count"),
                            "ligand_start_generation_method": prior_metadata.get("ligand_start_generation_method", "")})
            continue
        if target.exists() and not force and not retry_nonready:
            results.append({"case_id": cid, "pair_id": row["pair_id"], "status": "P4_PREPARATION_REVIEW",
                            "reason": "EXISTING_NON_SUCCESS_CASE_REQUIRES_FORCE", "receptor_atoms": None,
                            "ligand_heavy_atoms": None, "ligand_start_generation_method": ""})
            continue
        tmp = case_root / ("." + cid + f".tmp-{os.getpid()}")
        if tmp.exists():
            shutil.rmtree(tmp)
        tmp.mkdir()
        placement = str(row["ligand_assembly_placement_id"])
        chains = [x for x in str(row["receptor_chain_instance_ids"]).split(",") if x]
        try:
            result = create_case(
                tmp, row,
                rec_atoms[rec_atoms["filter_1_chain_instance_id"].astype(str).isin(chains)],
                lig_atoms[lig_atoms["filter_2_ligand_assembly_placement_id"].astype(str) == placement],
                lig_bonds[lig_bonds["filter_2_ligand_assembly_placement_id"].astype(str) == placement],
                pocket[pocket["pair_id"].astype(str) == str(row["pair_id"])] if not pocket.empty else pocket,
                binding[binding["ligand_assembly_placement_id"].astype(str) == placement] if not binding.empty else binding,
                ccd, config,
            )
        except Exception as exc:
            if isinstance(exc, TooSmallForPoseDocking):
                status = "P4_OUT_OF_SCOPE_TOO_SMALL_FOR_POSE_DOCKING"
            elif isinstance(exc, LigandStartFailure):
                status = "P4_LIGAND_START_GENERATION_FAILED"
            else:
                status = "P4_PREPARATION_REVIEW"
            error = f"{type(exc).__name__}: {exc}"[:4000]
            metadata = {"case_id": cid, "pair_id": str(row["pair_id"]),
                        "pdb_id": str(row["pdb_id"]), "component_id": str(row["component_id"]),
                        "processing4_status": status, "error": error, "created_at": utc()}
            if isinstance(exc, TooSmallForPoseDocking):
                metadata["ligand_heavy_atom_count"] = exc.heavy_atom_count
                metadata["scope_interpretation"] = "scientifically valid database ligand; outside heavy-atom RMSD pose-docking scope"
            atomic_json(tmp / "metadata.json", metadata)
            marker_name = "_OUT_OF_SCOPE.json" if isinstance(exc, TooSmallForPoseDocking) else "_REVIEW.json"
            atomic_json(tmp / marker_name, {"status": status, "error": error,
                                           "traceback": traceback.format_exc(limit=8), "created_at": utc()})
            result = {"case_id": cid, "pair_id": row["pair_id"], "status": status, "reason": error,
                      "receptor_atoms": None,
                      "ligand_heavy_atoms": exc.heavy_atom_count if isinstance(exc, TooSmallForPoseDocking) else None,
                      "ligand_start_generation_method": ""}
        if target.exists() and (force or retry_nonready):
            backup_root = run / f"work/superseded_nonready/bucket_{bid:03d}"
            backup_root.mkdir(parents=True, exist_ok=True)
            backup = backup_root / (target.name + ".superseded-" + datetime.now().strftime("%Y%m%d%H%M%S"))
            os.replace(target, backup)
        os.replace(tmp, target)
        results.append(result)
    atomic_parquet(status_path, pd.DataFrame(results))
    ready = sum(r["status"] == "P4_DOCKING_READY" for r in results)
    return {"bucket_id": bid, "cases": len(results), "ready": ready,
            "review": len(results) - ready, "seconds": time.time() - started}


def run_command(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    if (run / "_FROZEN.json").exists():
        raise RuntimeError("refusing to execute frozen run")
    inv = pq.read_table(run / "input/case_inventory.parquet", columns=["bucket_id"]).to_pandas()
    buckets = sorted(set(int(x) for x in inv["bucket_id"]))
    if args.bucket is not None:
        buckets = [int(args.bucket)]
    summaries = []
    if args.workers == 1 or len(buckets) == 1:
        for bid in buckets:
            summary = run_bucket(str(run), bid, args.force, args.retry_nonready)
            summaries.append(summary)
            print(json.dumps(summary), flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_bucket, str(run), bid, args.force, args.retry_nonready): bid for bid in buckets}
            for future in as_completed(futures):
                summary = future.result()
                summaries.append(summary)
                print(json.dumps(summary), flush=True)
    atomic_json(run / "run_summary.json", {"status": "EXECUTED", "finished_at": utc(),
                                           "buckets": sorted(summaries, key=lambda x: x["bucket_id"])})


def validate_command(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    config = read_config(run / "input/config_snapshot.json")
    inventory = pq.read_table(run / "input/case_inventory.parquet").to_pandas()
    status_files = sorted((run / "work/buckets").glob("bucket_*.parquet"))
    statuses = pd.concat([pq.read_table(p).to_pandas() for p in status_files], ignore_index=True) if status_files else pd.DataFrame()
    duplicated = int(statuses["case_id"].duplicated().sum()) if not statuses.empty else 0
    missing = sorted(set(inventory["case_id"]) - set(statuses.get("case_id", [])))
    extras = sorted(set(statuses.get("case_id", [])) - set(inventory["case_id"]))
    counts = statuses["status"].value_counts().to_dict() if not statuses.empty else {}
    ready_file_errors = []
    nonready_contract_errors = []
    case_to_bucket = dict(zip(inventory["case_id"].astype(str), inventory["bucket_id"].astype(int)))
    old_ready_status = pq.read_table(
        OLD_P4_RUN / "output/processing4_case_inventory.parquet",
        columns=["case_id", "status"],
    ).to_pandas()
    old_ready_case_ids = set(old_ready_status.loc[
        old_ready_status["status"].eq("P4_DOCKING_READY"), "case_id"
    ].astype(str))
    # A restart also labels successes created earlier in this same run as reused.
    # Only cases that are members of the frozen old P4 READY set may inherit its manifest lineage.
    reused_case_ids = set(statuses.loc[
        statuses["reason"].eq("REUSED_VALIDATED_PRIOR_SUCCESS")
        & statuses["case_id"].astype(str).isin(old_ready_case_ids), "case_id"
    ].astype(str))
    old_marker = json.loads((OLD_P4_RUN / "_FROZEN.json").read_text())
    old_validation = json.loads((OLD_P4_RUN / "validation/validation.json").read_text())
    reuse_lineage_valid = (
        old_marker.get("status") == "FROZEN"
        and old_validation.get("status") == "PASS"
        and old_marker.get("validation_sha256") == sha256_file(OLD_P4_RUN / "validation/validation.json")
        and old_marker.get("output_manifest_sha256") == sha256_file(OLD_P4_RUN / "output_manifest.parquet")
    )
    for row in statuses[statuses["status"] == "P4_DOCKING_READY"].itertuples() if not statuses.empty else []:
        if str(row.case_id) in reused_case_ids:
            continue
        bid = case_to_bucket[str(row.case_id)]
        root = run / f"output/cases/bucket_{bid:03d}" / row.case_id
        absent = [f for f in READY_FILES if not (root / f).is_file()]
        try:
            success = json.loads((root / "_SUCCESS.json").read_text())
            metadata = json.loads((root / "metadata.json").read_text())
            site = json.loads((root / "site.json").read_text())
            hash_bad = [f for f in READY_FILES if f not in absent and success.get("sha256", {}).get(f) != sha256_file(root / f)]
            sizes = site["search_box"]
            box_bad = any(float(sizes[key]) < float(config["site_minimum_dimension_angstrom"])
                          for key in ("size_x", "size_y", "size_z"))
            method_bad = metadata.get("ligand_start_generation_method") not in {
                "CCD_IDEAL", "ETKDGv3_UFF_FALLBACK",
                "ETKDGv3_UFF_FALLBACK_NO_UFF", "ETKDGv3_RANDOM_COORDS_FALLBACK",
                "ETKDGv3_RANDOM_COORDS_FALLBACK_NO_UFF",
            }
            status_bad = metadata.get("processing4_status") != "P4_DOCKING_READY"
        except Exception as exc:
            hash_bad, box_bad, method_bad, status_bad = [], True, True, True
            absent.append(f"validation_exception:{type(exc).__name__}:{exc}")
        if absent or hash_bad or box_bad or method_bad or status_bad:
            ready_file_errors.append({"case_id": row.case_id, "missing": absent,
                                      "hash_bad": hash_bad, "box_bad": box_bad,
                                      "method_bad": method_bad, "status_bad": status_bad})
    for row in statuses[statuses["status"] != "P4_DOCKING_READY"].itertuples() if not statuses.empty else []:
        bid = case_to_bucket[str(row.case_id)]
        root = run / f"output/cases/bucket_{bid:03d}" / row.case_id
        try:
            metadata = json.loads((root / "metadata.json").read_text())
            if metadata.get("processing4_status") != row.status:
                raise RuntimeError("metadata status mismatch")
            if row.status == "P4_OUT_OF_SCOPE_TOO_SMALL_FOR_POSE_DOCKING":
                if not (root / "_OUT_OF_SCOPE.json").is_file():
                    raise RuntimeError("out-of-scope marker missing")
                if int(metadata.get("ligand_heavy_atom_count", 999)) >= int(config["minimum_ligand_heavy_atom_count"]):
                    raise RuntimeError("out-of-scope heavy atom count is not below threshold")
            elif not (root / "_REVIEW.json").is_file():
                raise RuntimeError("review/failure marker missing")
        except Exception as exc:
            nonready_contract_errors.append({"case_id": row.case_id, "error": f"{type(exc).__name__}: {exc}"})
    allowed_statuses = {
        "P4_DOCKING_READY", "P4_OUT_OF_SCOPE_TOO_SMALL_FOR_POSE_DOCKING",
        "P4_LIGAND_START_GENERATION_FAILED", "P4_PREPARATION_REVIEW",
    }
    unexpected_statuses = sorted(set(counts) - allowed_statuses)
    scientific_non_scope_statuses = sorted(set(counts) - {
        "P4_DOCKING_READY", "P4_OUT_OF_SCOPE_TOO_SMALL_FOR_POSE_DOCKING"
    })
    report = {
        "status": "PASS" if not missing and not extras and duplicated == 0 and not ready_file_errors
        and not nonready_contract_errors and not unexpected_statuses and not scientific_non_scope_statuses
        and reuse_lineage_valid else "FAIL",
        "expected_cases": len(inventory), "status_rows": len(statuses), "status_counts": counts,
        "missing_case_count": len(missing), "extra_case_count": len(extras), "duplicate_case_count": duplicated,
        "ready_file_error_count": len(ready_file_errors),
        "nonready_contract_error_count": len(nonready_contract_errors),
        "unexpected_statuses": unexpected_statuses,
        "scientific_non_scope_statuses": scientific_non_scope_statuses,
        "validation_mode": "FULL_FOR_REGENERATED_AND_NONREADY_PLUS_FROZEN_MANIFEST_LINEAGE_FOR_HARDLINK_REUSE",
        "reused_ready_case_count": len(reused_case_ids),
        "regenerated_ready_case_count": int((statuses["status"] == "P4_DOCKING_READY").sum()) - len(reused_case_ids),
        "reuse_lineage_valid": reuse_lineage_valid,
        "reuse_source_run": str(OLD_P4_RUN),
        "missing_examples": missing[:20],
        "ready_file_error_examples": ready_file_errors[:20], "validated_at": utc(),
        "nonready_contract_error_examples": nonready_contract_errors[:20],
        "site_definition_version": config["site_definition_version"],
        "minimum_ligand_heavy_atom_count": config["minimum_ligand_heavy_atom_count"],
    }
    atomic_json(run / "validation/validation.json", report)
    if len(statuses):
        atomic_parquet(run / "output/processing4_case_inventory.parquet", statuses.sort_values("case_id"))
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(2)


def freeze_command(args: argparse.Namespace) -> None:
    run = Path(args.run_dir).resolve()
    validation = json.loads((run / "validation/validation.json").read_text())
    if validation["status"] != "PASS":
        raise RuntimeError("validation is not PASS")
    config = read_config(run / "input/config_snapshot.json")
    if str(config["policy_version"]).endswith("-pilot") and not args.allow_pilot:
        raise RuntimeError("pilot policy cannot be formally frozen without --allow-pilot")
    statuses = pq.read_table(run / "output/processing4_case_inventory.parquet").to_pandas()
    upstream = pq.read_table(
        run / "input/case_inventory.parquet",
        columns=["pair_id", "filter3_quality_class"],
    ).to_pandas()
    joined = statuses.merge(upstream, on="pair_id", how="left", validate="one_to_one")
    release_summary = {
        "stage": "Processing 4: Docking-Ready Case Construction",
        "policy_version": config["policy_version"],
        "formal_input_total": int(len(joined)),
        "formal_input_high": int((joined["filter3_quality_class"] == "FILTER3_HIGH_QUALITY").sum()),
        "formal_input_good": int((joined["filter3_quality_class"] == "FILTER3_GOOD_QUALITY").sum()),
        "status_counts": {str(k): int(v) for k, v in joined["status"].value_counts().items()},
        "quality_by_status": {
            str(status): {str(k): int(v) for k, v in frame["filter3_quality_class"].value_counts().items()}
            for status, frame in joined.groupby("status", sort=True)
        },
        "ligand_start_generation_method_counts": {
            str(k): int(v) for k, v in joined.loc[
                joined["status"] == "P4_DOCKING_READY", "ligand_start_generation_method"
            ].value_counts().items()
        },
        "site_definition_version": config["site_definition_version"],
        "minimum_ligand_heavy_atom_count": int(config["minimum_ligand_heavy_atom_count"]),
        "native_coordinate_fallback": "FORBIDDEN",
        "created_at": utc(),
    }
    atomic_json(run / "output/release_summary.json", release_summary)
    manifest = []
    old_ready_status = pq.read_table(
        OLD_P4_RUN / "output/processing4_case_inventory.parquet",
        columns=["case_id", "status"],
    ).to_pandas()
    old_ready_case_ids = set(old_ready_status.loc[
        old_ready_status["status"].eq("P4_DOCKING_READY"), "case_id"
    ].astype(str))
    reused_case_ids = set(statuses.loc[
        statuses["reason"].eq("REUSED_VALIDATED_PRIOR_SUCCESS")
        & statuses["case_id"].astype(str).isin(old_ready_case_ids), "case_id"
    ].astype(str))
    old_manifest = pq.read_table(OLD_P4_RUN / "output_manifest.parquet").to_pandas()
    for row in old_manifest.itertuples(index=False):
        parts = Path(str(row.relative_path)).parts
        if len(parts) >= 5 and parts[0] == "output" and parts[1] == "cases" and parts[3] in reused_case_ids:
            manifest.append({"relative_path": str(row.relative_path), "size_bytes": int(row.size_bytes),
                             "sha256": str(row.sha256)})
    for path in sorted((run / "output").rglob("*")):
        if path.is_file():
            relative = path.relative_to(run)
            parts = relative.parts
            if len(parts) >= 5 and parts[0] == "output" and parts[1] == "cases" and parts[3] in reused_case_ids:
                continue
            manifest.append({"relative_path": str(relative), "size_bytes": path.stat().st_size,
                             "sha256": sha256_file(path)})
    manifest_frame = pd.DataFrame(manifest).drop_duplicates("relative_path", keep="last").sort_values("relative_path")
    atomic_parquet(run / "output_manifest.parquet", manifest_frame)
    marker = {"status": "FROZEN", "stage": "processing4_benchmark_v3.1.0_filter5_v3_exact_only", "frozen_at": utc(),
              "formal_input_pairs": int(len(joined)),
              "formal_ready_pairs": int((joined["status"] == "P4_DOCKING_READY").sum()),
              "formal_status_counts": release_summary["status_counts"],
              "release_summary_sha256": sha256_file(run / "output/release_summary.json"),
              "validation_sha256": sha256_file(run / "validation/validation.json"),
              "output_manifest_sha256": sha256_file(run / "output_manifest.parquet")}
    atomic_json(run / "_FROZEN.json", marker)
    atomic_json(run.parent.parent / "CURRENT_RUN.json", {
        "current_run_id": run.name,
        "relative_path": str(run.relative_to(run.parent.parent)),
        "status": "FROZEN",
        "updated_at": marker["frozen_at"],
        "output_manifest_sha256": marker["output_manifest_sha256"],
    })
    print(json.dumps(marker, indent=2))


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("prepare-run")
    q.add_argument("--run-dir", required=True)
    q.add_argument("--config", default=str(DEFAULT_CONFIG))
    q.add_argument("--mode", choices=["smoke", "full"], required=True)
    q.add_argument("--limit", type=int, default=100)
    q.set_defaults(func=prepare_run)
    q = sub.add_parser("run")
    q.add_argument("--run-dir", required=True)
    q.add_argument("--workers", type=int, default=1,
                   help="Concurrent bucket readers; keep at 1-2 on the shared data volume")
    q.add_argument("--bucket", type=int)
    q.add_argument("--force", action="store_true")
    q.add_argument("--retry-nonready", action="store_true",
                   help="Keep successful cases and atomically rebuild existing review/failure cases")
    q.set_defaults(func=run_command)
    q = sub.add_parser("validate")
    q.add_argument("--run-dir", required=True)
    q.set_defaults(func=validate_command)
    q = sub.add_parser("freeze")
    q.add_argument("--run-dir", required=True)
    q.add_argument("--allow-pilot", action="store_true")
    q.set_defaults(func=freeze_command)
    return p


if __name__ == "__main__":
    ns = parser().parse_args()
    ns.func(ns)
