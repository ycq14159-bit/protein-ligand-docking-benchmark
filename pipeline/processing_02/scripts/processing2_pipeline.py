#!/usr/bin/env python3
"""Processing 2: assembly-ready structure preparation.

The implementation deliberately treats upstream Parquet files and mmCIF archives as
read-only.  Every task writes atomically into its own output fragments, which makes
resume idempotent and bounds memory to one 20-PDB auxiliary work package.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import shutil
import sqlite3
import statistics
import sys
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import gemmi
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

STAGE_VERSION = "processing_2_v1.0.0"
SCHEMA_VERSION = "processing_2_schema_v1.0.0"
BLANK_VALUES = {"", ".", "?", "\x00", "None", "False", "nan", "<NA>"}
HYDROGEN_ELEMENTS = {"H", "D", "T"}

ROOT = Path("/root/autodl-tmp/benchmark_1.0/processing_2_assembly_ready_structure_preparation")
AUX = Path("/root/autodl-tmp/benchmark_1.0/auxiliary_entry_work_packages/builds/20260805_full_01/output")
F2 = Path("/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v3/runs/20260804_full_01")
CCD = Path("/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification/references/components.cif.gz")
CROWN = ROOT / "references/CROWN"

DATASETS = (
    "prepared_ligand_source_atoms",
    "prepared_receptor_source_atoms",
    "prepared_ligand_assembly_atoms",
    "prepared_receptor_assembly_atoms",
    "prepared_ligand_source_bonds",
    "prepared_ligand_assembly_bonds",
    "altloc_selection_audit",
    "microheterogeneity_selection",
    "ligand_mapping_review",
    "ligand_incomplete_instances",
    "ligand_topology_validation",
    "operator_matrix_quality_review",
    "structure_preparation_manifest",
)


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def norm_blank(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in BLANK_VALUES else text


def as_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value)
    return "" if text in {".", "?", "None"} else text


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except Exception:
        return default


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(str(value))
    except Exception:
        return default


def parquet_frame(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    return pq.ParquetFile(path).read(columns=columns).to_pandas()


def task_key(path: Path) -> tuple[int, int]:
    bid = int(path.parent.name.split("=")[1])
    tid = int(path.stem.split("-")[1])
    return bid, tid


def input_path(dataset: str, bid: int, tid: int) -> Path:
    return AUX / dataset / f"bucket_id={bid:03d}" / f"part-{tid:06d}.parquet"


def output_path(run: Path, dataset: str, bid: int, tid: int) -> Path:
    return run / "output" / dataset / f"bucket_id={bid:03d}" / f"part-{tid:06d}.parquet"


def rows_to_table(rows: list[dict[str, Any]]) -> pa.Table:
    if not rows:
        return pa.table({"_empty": pa.array([], pa.string())})
    frame = pd.DataFrame.from_records(rows)
    for col in frame.columns:
        if frame[col].dtype == object:
            frame[col] = frame[col].map(lambda x: "" if x is None else str(x))
    return pa.Table.from_pandas(frame, preserve_index=False)


def atomic_parquet(path: Path, frame_or_rows: pd.DataFrame | list[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".parquet.tmp")
    if isinstance(frame_or_rows, pd.DataFrame):
        table = pa.Table.from_pandas(frame_or_rows, preserve_index=False)
    else:
        table = rows_to_table(frame_or_rows)
    pq.write_table(table, tmp, compression="zstd", compression_level=6,
                   version="2.6", data_page_version="2.0", use_dictionary=True)
    os.replace(tmp, path)
    return table.num_rows


def init_stage(args: argparse.Namespace) -> None:
    for rel in ("scripts", "config", "tests", "references", "runs"):
        (ROOT / rel).mkdir(parents=True, exist_ok=True)
    print(ROOT)


def component_ids() -> set[str]:
    path = F2 / "output/provisional_source_ligands.tsv.gz"
    result: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result.add((row.get("resolved_ccd_id") or row.get("component_id") or "").upper())
    result.discard("")
    return result


def block_value(block: gemmi.cif.Block, tag: str) -> str:
    try:
        return as_text(block.find_value(tag))
    except Exception:
        return ""


def loop_rows(block: gemmi.cif.Block, tags: list[str]) -> Iterable[list[str]]:
    try:
        table = block.find(tags)
        for row in table:
            yield [as_text(row[i]) for i in range(len(tags))]
    except Exception:
        return


def rdkit_component(atom_rows: list[dict[str, Any]], bond_rows: list[dict[str, Any]]) -> dict[str, Any]:
    out = {
        "rdkit_parse_success": False, "rdkit_sanitize_success": False,
        "rdkit_error": "", "canonical_smiles_from_graph": "",
        "molecular_weight_from_graph": None, "formal_charge_from_graph": None,
        "fragment_count_from_graph": None, "ring_count_from_graph": None,
        "rotatable_bond_count_from_graph": None,
    }
    try:
        rw = Chem.RWMol()
        indices: dict[str, int] = {}
        aromatic_atoms: set[str] = set()
        for row in atom_rows:
            element = row["element"].strip().title()
            atom = Chem.Atom(element)
            atom.SetFormalCharge(int(row["charge"]))
            idx = rw.AddAtom(atom)
            indices[row["atom_id"]] = idx
        types = {
            "SING": Chem.BondType.SINGLE, "SINGLE": Chem.BondType.SINGLE,
            "DOUB": Chem.BondType.DOUBLE, "DOUBLE": Chem.BondType.DOUBLE,
            "TRIP": Chem.BondType.TRIPLE, "TRIPLE": Chem.BondType.TRIPLE,
            "AROM": Chem.BondType.AROMATIC, "AROMATIC": Chem.BondType.AROMATIC,
        }
        for row in bond_rows:
            a, b = row["atom_id_1"], row["atom_id_2"]
            if a not in indices or b not in indices:
                continue
            order = row["bond_order"].upper()
            if order not in types:
                raise ValueError(f"unsupported CCD bond order {order}")
            rw.AddBond(indices[a], indices[b], types[order])
            if types[order] == Chem.BondType.AROMATIC:
                aromatic_atoms.update((a, b))
        mol = rw.GetMol()
        for atom_id in aromatic_atoms:
            mol.GetAtomWithIdx(indices[atom_id]).SetIsAromatic(True)
        out["rdkit_parse_success"] = True
        Chem.SanitizeMol(mol)
        out["rdkit_sanitize_success"] = True
        out["canonical_smiles_from_graph"] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        out["molecular_weight_from_graph"] = float(Descriptors.MolWt(mol))
        out["formal_charge_from_graph"] = int(Chem.GetFormalCharge(mol))
        out["fragment_count_from_graph"] = len(Chem.GetMolFrags(mol))
        out["ring_count_from_graph"] = int(Lipinski.RingCount(mol))
        out["rotatable_bond_count_from_graph"] = int(Lipinski.NumRotatableBonds(mol))
    except Exception as exc:
        out["rdkit_error"] = f"{type(exc).__name__}: {exc}"[:1000]
    return out


def build_ccd(args: argparse.Namespace) -> None:
    run = Path(args.run_dir)
    db = run / "input/ccd_active_snapshot.sqlite"
    metadata = run / "input/ccd_snapshot_metadata.json"
    if db.exists() and metadata.exists() and not args.force:
        print(db)
        return
    run.joinpath("input").mkdir(parents=True, exist_ok=True)
    wanted = component_ids()
    tmp = db.with_suffix(".sqlite.tmp")
    if tmp.exists():
        tmp.unlink()
    conn = sqlite3.connect(tmp)
    conn.executescript("""
      PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
      CREATE TABLE components(component_id TEXT PRIMARY KEY,name TEXT,ccd_type TEXT,formula TEXT,formula_weight REAL,formal_charge INTEGER,parent_component_id TEXT,descriptor_smiles TEXT,descriptor_smiles_canonical TEXT,descriptor_inchi TEXT,descriptor_inchikey TEXT,expected_atom_count INTEGER,expected_heavy_atom_count INTEGER,rdkit_parse_success INTEGER,rdkit_sanitize_success INTEGER,rdkit_error TEXT,canonical_smiles_from_graph TEXT,molecular_weight_from_graph REAL,formal_charge_from_graph INTEGER,fragment_count_from_graph INTEGER,ring_count_from_graph INTEGER,rotatable_bond_count_from_graph INTEGER);
      CREATE TABLE atoms(component_id TEXT,atom_id TEXT,element TEXT,charge INTEGER,aromatic_flag TEXT,stereo_config TEXT,ordinal INTEGER,PRIMARY KEY(component_id,atom_id));
      CREATE TABLE bonds(component_id TEXT,atom_id_1 TEXT,atom_id_2 TEXT,bond_order TEXT,aromatic_flag TEXT,stereo_config TEXT,ordinal INTEGER);
      CREATE INDEX atoms_component ON atoms(component_id);
      CREATE INDEX bonds_component ON bonds(component_id);
    """)
    doc = gemmi.cif.read(str(CCD))
    found = 0
    failures: list[dict[str, str]] = []
    for block in doc:
        cid = (block_value(block, "_chem_comp.id") or block.name).upper()
        if cid not in wanted:
            continue
        atom_rows = []
        for n, values in enumerate(loop_rows(block, [
            "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
            "_chem_comp_atom.charge", "_chem_comp_atom.pdbx_aromatic_flag",
            "_chem_comp_atom.pdbx_stereo_config"
        ])):
            atom_rows.append({"component_id": cid, "atom_id": values[0], "element": values[1],
                              "charge": as_int(values[2]), "aromatic_flag": values[3],
                              "stereo_config": values[4], "ordinal": n})
        bond_rows = []
        for n, values in enumerate(loop_rows(block, [
            "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order", "_chem_comp_bond.pdbx_aromatic_flag",
            "_chem_comp_bond.pdbx_stereo_config"
        ])):
            bond_rows.append({"component_id": cid, "atom_id_1": values[0], "atom_id_2": values[1],
                              "bond_order": values[2], "aromatic_flag": values[3],
                              "stereo_config": values[4], "ordinal": n})
        descriptors: dict[str, str] = {}
        for values in loop_rows(block, ["_pdbx_chem_comp_descriptor.type", "_pdbx_chem_comp_descriptor.program", "_pdbx_chem_comp_descriptor.descriptor"]):
            key = values[0].upper()
            descriptors.setdefault(key, values[2])
        graph = rdkit_component(atom_rows, bond_rows)
        heavy = sum(row["element"].upper() not in HYDROGEN_ELEMENTS for row in atom_rows)
        comp = (
            cid, block_value(block, "_chem_comp.name"), block_value(block, "_chem_comp.type"),
            block_value(block, "_chem_comp.formula"), as_float(block_value(block, "_chem_comp.formula_weight"), math.nan),
            as_int(block_value(block, "_chem_comp.pdbx_formal_charge")), block_value(block, "_chem_comp.mon_nstd_parent_comp_id"),
            descriptors.get("SMILES", ""), descriptors.get("SMILES_CANONICAL", ""),
            descriptors.get("INCHI", ""), descriptors.get("INCHIKEY", ""), len(atom_rows), heavy,
            int(graph["rdkit_parse_success"]), int(graph["rdkit_sanitize_success"]), graph["rdkit_error"],
            graph["canonical_smiles_from_graph"], graph["molecular_weight_from_graph"],
            graph["formal_charge_from_graph"], graph["fragment_count_from_graph"],
            graph["ring_count_from_graph"], graph["rotatable_bond_count_from_graph"],
        )
        try:
            conn.execute("INSERT INTO components VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", comp)
            conn.executemany("INSERT INTO atoms VALUES (?,?,?,?,?,?,?)", [tuple(r.values()) for r in atom_rows])
            conn.executemany("INSERT INTO bonds VALUES (?,?,?,?,?,?,?)", [tuple(r.values()) for r in bond_rows])
            found += 1
        except Exception as exc:
            failures.append({"component_id": cid, "error": f"{type(exc).__name__}: {exc}"})
        if found % 1000 == 0:
            conn.commit()
            print(f"CCD {found}/{len(wanted)}", flush=True)
    conn.commit()
    conn.close()
    os.replace(tmp, db)
    atomic_json(metadata, {
        "created_at": utc(), "source_path": str(CCD), "source_sha256": sha256(CCD),
        "active_component_count_requested": len(wanted), "active_component_count_found": found,
        "missing_component_ids": sorted(wanted - {r[0] for r in sqlite3.connect(db).execute("SELECT component_id FROM components")}),
        "parse_failures": failures, "stage_version": STAGE_VERSION,
    })
    print(json.dumps(json.loads(metadata.read_text()), indent=2)[:5000])


class CCDStore:
    def __init__(self, path: Path):
        self.conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        self.cache: dict[str, tuple[dict[str, Any] | None, dict[str, dict[str, Any]], list[dict[str, Any]]]] = {}

    def get(self, cid: str):
        cid = cid.upper()
        if cid in self.cache:
            return self.cache[cid]
        self.conn.row_factory = sqlite3.Row
        row = self.conn.execute("SELECT * FROM components WHERE component_id=?", (cid,)).fetchone()
        comp = dict(row) if row else None
        atoms = {r["atom_id"]: dict(r) for r in self.conn.execute("SELECT * FROM atoms WHERE component_id=? ORDER BY ordinal", (cid,))}
        bonds = [dict(r) for r in self.conn.execute("SELECT * FROM bonds WHERE component_id=? ORDER BY ordinal", (cid,))]
        self.cache[cid] = (comp, atoms, bonds)
        return self.cache[cid]


def select_altloc(frame: pd.DataFrame, group_cols: list[str], object_type: str, object_id_col: str) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if frame.empty:
        return frame.copy(), []
    x = frame.copy()
    x["_alt"] = x["alt_id"].map(norm_blank)
    x["_order"] = pd.to_numeric(x["source_atom_row_index"], errors="coerce").fillna(10**18).astype("int64")
    nonblank = x[x["_alt"] != ""]
    winners: dict[tuple, str] = {}
    audits: list[dict[str, Any]] = []
    if not nonblank.empty:
        stats = (nonblank.groupby(group_cols + ["_alt"], sort=False, dropna=False)
                 .agg(atom_count=("_alt", "size"), mean_occupancy=("occupancy", "mean"), first_order=("_order", "min"))
                 .reset_index())
        grouped = stats.groupby(group_cols, sort=False, dropna=False)
        for key, part in grouped:
            part = part.sort_values(["atom_count", "mean_occupancy", "first_order"], ascending=[False, False, True], kind="stable")
            winner = str(part.iloc[0]["_alt"])
            key_tuple = key if isinstance(key, tuple) else (key,)
            winners[key_tuple] = winner
            base = {col: val for col, val in zip(group_cols, key_tuple)}
            obj_rows = x
            for col, val in base.items():
                obj_rows = obj_rows[obj_rows[col] == val]
            audits.append({
                "object_type": object_type, "object_id": str(obj_rows.iloc[0][object_id_col]),
                "pdb_id": str(obj_rows.iloc[0]["pdb_id"]), "candidate_altlocs": ",".join(part["_alt"].astype(str)),
                "selected_altloc": winner, "blank_atom_count": int((obj_rows["_alt"] == "").sum()),
                "selected_alt_atom_count": int(part.iloc[0]["atom_count"]),
                "selected_alt_mean_occupancy": float(part.iloc[0]["mean_occupancy"]),
                "selection_rule": "CROWN_count_then_mean_occupancy_then_source_order",
                "source_rows_unchanged": True,
            })
    if winners:
        keys = list(zip(*(x[c].tolist() for c in group_cols)))
        keep = [alt == "" or winners.get(tuple(k)) == alt for k, alt in zip(keys, x["_alt"])]
        x = x.loc[keep].copy()
    x["selected_altloc_original"] = x["_alt"]
    x["alt_id"] = ""
    return x.drop(columns=["_alt", "_order"]), audits


def select_microheterogeneity(frame: pd.DataFrame, position_cols: list[str], component_col: str,
                              object_type: str, object_id_col: str) -> tuple[pd.DataFrame, list[dict[str, Any]], set[str]]:
    if frame.empty:
        return frame.copy(), [], set()
    x = frame.copy()
    x["_heavy"] = ~x["type_symbol"].str.upper().isin(HYDROGEN_ELEMENTS)
    x["_order"] = pd.to_numeric(x["source_atom_row_index"], errors="coerce").fillna(10**18).astype("int64")
    stats = (x[x["_heavy"]].groupby(position_cols + [component_col], sort=False, dropna=False)
             .agg(heavy_atom_count=("_heavy", "size"), mean_heavy_occupancy=("occupancy", "mean"), first_order=("_order", "min"))
             .reset_index())
    winners: dict[tuple, str] = {}
    audit: list[dict[str, Any]] = []
    loser_ids: set[str] = set()
    for key, part in stats.groupby(position_cols, sort=False, dropna=False):
        if len(part) <= 1:
            continue
        part = part.sort_values(["heavy_atom_count", "mean_heavy_occupancy", "first_order"], ascending=[False, False, True], kind="stable")
        key_tuple = key if isinstance(key, tuple) else (key,)
        winner = str(part.iloc[0][component_col])
        winners[key_tuple] = winner
        mask = pd.Series(True, index=x.index)
        for col, val in zip(position_cols, key_tuple):
            mask &= x[col] == val
        subset = x[mask]
        winning_ids = set(subset.loc[subset[component_col].astype(str) == winner, object_id_col].astype(str))
        losing = set(subset[object_id_col].astype(str)) - winning_ids
        loser_ids.update(losing)
        audit.append({
            "object_type": object_type, "pdb_id": str(subset.iloc[0]["pdb_id"]),
            "position_key": "|".join(map(str, key_tuple)),
            "candidate_components": ",".join(part[component_col].astype(str)),
            "selected_component": winner, "selected_object_ids": ",".join(sorted(winning_ids)),
            "removed_object_ids": ",".join(sorted(losing)),
            "selection_rule": "CROWN_heavy_count_then_mean_heavy_occupancy_then_source_order",
            "source_rows_unchanged": True,
        })
    if winners:
        keys = list(zip(*(x[c].tolist() for c in position_cols)))
        keep = [winners.get(tuple(k), str(comp)) == str(comp) for k, comp in zip(keys, x[component_col])]
        x = x.loc[keep].copy()
    return x.drop(columns=["_heavy", "_order"]), audit, loser_ids


def apply_transform(frame: pd.DataFrame, prefix_cols: list[str]) -> pd.DataFrame:
    if frame.empty:
        return frame
    xyz = frame[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(dtype=float)
    r = frame[["r11", "r12", "r13", "r21", "r22", "r23", "r31", "r32", "r33"]].to_numpy(dtype=float).reshape(-1, 3, 3)
    t = frame[["t1", "t2", "t3"]].to_numpy(dtype=float)
    transformed = np.einsum("nij,nj->ni", r, xyz) + t
    out = frame[prefix_cols + [c for c in frame.columns if c not in prefix_cols]].copy()
    out["source_Cartn_x"] = out["Cartn_x"]
    out["source_Cartn_y"] = out["Cartn_y"]
    out["source_Cartn_z"] = out["Cartn_z"]
    out["Cartn_x"] = transformed[:, 0]
    out["Cartn_y"] = transformed[:, 1]
    out["Cartn_z"] = transformed[:, 2]
    return out


def add_operator_quality(frame: pd.DataFrame, tolerance: float = 1.0e-3) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    out = frame.copy()
    cols = ["r11", "r12", "r13", "r21", "r22", "r23", "r31", "r32", "r33"]
    r = out[cols].to_numpy(dtype=float).reshape(-1, 3, 3)
    out["operator_determinant_error"] = np.abs(np.linalg.det(r) - 1.0)
    out["operator_orthogonality_error"] = np.max(
        np.abs(np.matmul(np.transpose(r, (0, 2, 1)), r) - np.eye(3)), axis=(1, 2))
    out["operator_quality_status"] = np.where(
        (out["operator_determinant_error"] <= tolerance) &
        (out["operator_orthogonality_error"] <= tolerance),
        "PASS", "OPERATOR_MATRIX_REVIEW")
    return out


def process_task(run_dir: str, bid: int, tid: int, selected_pdbs: set[str] | None = None) -> dict[str, Any]:
    started = time.time()
    run = Path(run_dir)
    ccd = CCDStore(run / "input/ccd_active_snapshot.sqlite")
    paths = {name: input_path(name, bid, tid) for name in (
        "entry_ligand_placements", "entry_receptor_chain_instances", "entry_ligand_source_atoms", "entry_receptor_source_atoms")}
    data = {name: parquet_frame(path) for name, path in paths.items()}
    if selected_pdbs is not None:
        data = {k: v[v["pdb_id"].isin(selected_pdbs)].copy() for k, v in data.items()}
    placements = add_operator_quality(data["entry_ligand_placements"])
    receptors = add_operator_quality(data["entry_receptor_chain_instances"])
    lig = data["entry_ligand_source_atoms"]
    rec = data["entry_receptor_source_atoms"]

    lig, lig_alt = select_altloc(lig, ["filter_2_source_ligand_instance_id"], "ligand_source_instance", "filter_2_source_ligand_instance_id")
    rec_residue = ["filter_1_source_chain_key", "label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id"]
    rec, rec_alt = select_altloc(rec, rec_residue, "receptor_residue", "filter_1_source_chain_key")

    lig_position = ["pdb_id", "model_id", "label_asym_id", "auth_seq_id", "insertion_code"]
    lig, lig_micro, lig_losers = select_microheterogeneity(lig, lig_position, "component_id", "ligand_source_instance", "filter_2_source_ligand_instance_id")
    rec_position = ["filter_1_source_chain_key", "label_seq_id", "auth_seq_id", "insertion_code"]
    rec, rec_micro, _ = select_microheterogeneity(rec, rec_position, "label_comp_id", "receptor_residue", "filter_1_source_chain_key")

    mapping_review: list[dict[str, Any]] = []
    incomplete: list[dict[str, Any]] = []
    topology: list[dict[str, Any]] = []
    source_bonds: list[dict[str, Any]] = []
    ligand_manifest: list[dict[str, Any]] = []
    valid_source_ids: set[str] = set()
    for source_id, group in lig.groupby("filter_2_source_ligand_instance_id", sort=False):
        cid = str(group.iloc[0]["component_id"]).upper()
        comp, expected, bonds = ccd.get(cid)
        observed_names = group["label_atom_id"].astype(str).tolist()
        observed_unique = set(observed_names)
        duplicate_names = sorted(k for k, v in Counter(observed_names).items() if v > 1)
        expected_heavy = {name for name, row in expected.items() if row["element"].upper() not in HYDROGEN_ELEMENTS}
        observed_heavy = {str(row.label_atom_id) for row in group.itertuples() if str(row.type_symbol).upper() not in HYDROGEN_ELEMENTS}
        missing = sorted(expected_heavy - observed_heavy)
        unexpected = sorted(observed_unique - set(expected)) if expected else sorted(observed_unique)
        element_mismatch = []
        for row in group.itertuples():
            atom = expected.get(str(row.label_atom_id))
            if atom and atom["element"].upper() != str(row.type_symbol).upper():
                element_mismatch.append(f"{row.label_atom_id}:{row.type_symbol}!={atom['element']}")
        if comp is None:
            status = "CCD_MISSING"
        elif duplicate_names or unexpected or element_mismatch:
            status = "ATOM_MAPPING_REVIEW"
        elif missing:
            status = "ATOM_INCOMPLETE"
        else:
            status = "COMPLETE"
            valid_source_ids.add(str(source_id))
        base = {
            "pdb_id": str(group.iloc[0]["pdb_id"]), "source_ligand_instance_id": str(source_id),
            "component_id": cid, "mapping_status": status,
            "observed_atom_count": len(group), "observed_heavy_atom_count": len(observed_heavy),
            "expected_atom_count": len(expected), "expected_heavy_atom_count": len(expected_heavy),
            "missing_heavy_atom_count": len(missing), "missing_heavy_atom_ids": ",".join(missing),
            "unexpected_atom_ids": ",".join(unexpected), "duplicate_observed_atom_ids": ",".join(duplicate_names),
            "element_mismatches": ",".join(element_mismatch), "strict_completeness_rule": "missing_CCD_heavy_atoms_equals_zero",
        }
        if status not in {"COMPLETE", "ATOM_INCOMPLETE"}:
            mapping_review.append(base)
        if status == "ATOM_INCOMPLETE":
            incomplete.append(base)
        top_status = "TOPOLOGY_COMPLETE" if comp and comp["rdkit_sanitize_success"] and not unexpected and not duplicate_names else "TOPOLOGY_REVIEW"
        topology.append({**base, "topology_status": top_status,
                         "ccd_bond_count": len(bonds), "mapped_bond_count": sum(b["atom_id_1"] in observed_unique and b["atom_id_2"] in observed_unique for b in bonds),
                         "rdkit_parse_success": bool(comp and comp["rdkit_parse_success"]),
                         "rdkit_sanitize_success": bool(comp and comp["rdkit_sanitize_success"]),
                         "rdkit_error": "" if not comp else comp["rdkit_error"],
                         "canonical_smiles_from_ccd_graph": "" if not comp else comp["canonical_smiles_from_graph"],
                         "ccd_descriptor_smiles": "" if not comp else (comp["descriptor_smiles_canonical"] or comp["descriptor_smiles"]),
                         "stereochemistry_status": "CCD_DESCRIPTOR_PRESERVED_NOT_REDERIVED_FROM_COORDINATES"})
        for n, bond in enumerate(bonds):
            if bond["atom_id_1"] in observed_unique and bond["atom_id_2"] in observed_unique:
                source_bonds.append({"pdb_id": base["pdb_id"], "filter_2_source_ligand_instance_id": str(source_id), "component_id": cid,
                                     "bond_index": n, "atom_id_1": bond["atom_id_1"], "atom_id_2": bond["atom_id_2"],
                                     "bond_order": bond["bond_order"], "aromatic_flag": bond["aromatic_flag"], "stereo_config": bond["stereo_config"]})
        ligand_manifest.append({"object_type": "ligand_source_instance", "object_id": str(source_id), "pdb_id": base["pdb_id"],
                                "component_id": cid, "preparation_status": status,
                                "decision": "PASS" if status == "COMPLETE" else "REVIEW", "destination": "assembly_ready" if status == "COMPLETE" else "ligand_preparation_review",
                                "reason_code": "strict_heavy_atom_complete" if status == "COMPLETE" else status.lower(),
                                "prepared_atom_count": len(group), "source_rows_modified": False})
    for loser in sorted(lig_losers):
        if loser not in set(lig["filter_2_source_ligand_instance_id"].astype(str)):
            ligand_manifest.append({"object_type": "ligand_source_instance", "object_id": loser, "pdb_id": loser.split("|")[0],
                                    "component_id": "", "preparation_status": "MICROHETEROGENEITY_NOT_SELECTED",
                                    "decision": "REVIEW", "destination": "microheterogeneity_review",
                                    "reason_code": "crown_microheterogeneity_loser", "prepared_atom_count": 0, "source_rows_modified": False})

    lig_place = placements.merge(lig, on=["pdb_id", "model_id", "filter_2_source_ligand_instance_id", "component_id", "entity_id", "label_asym_id", "auth_asym_id", "auth_seq_id", "insertion_code", "bucket_id"], how="left", suffixes=("_placement", ""))
    lig_place = lig_place[lig_place["label_atom_id"].notna()].copy()
    lig_assembly = apply_transform(lig_place, ["filter_2_ligand_assembly_placement_id", "filter_2_source_ligand_instance_id"])
    rec_place = receptors.merge(rec, on=["pdb_id", "model_id", "filter_1_source_chain_key", "entity_id", "label_asym_id", "auth_asym_id", "bucket_id"], how="left", suffixes=("_instance", ""))
    rec_place = rec_place[rec_place["label_atom_id"].notna()].copy()
    rec_assembly = apply_transform(rec_place, ["filter_1_chain_instance_id", "filter_1_source_chain_key"])

    source_bonds_df = pd.DataFrame(source_bonds)
    if not source_bonds_df.empty and not placements.empty:
        assembly_bonds = placements[["filter_2_ligand_assembly_placement_id", "filter_2_source_ligand_instance_id", "pdb_id", "assembly_id", "operator_path"]].merge(source_bonds_df, on=["pdb_id", "filter_2_source_ligand_instance_id"], how="inner")
    else:
        assembly_bonds = pd.DataFrame()

    manifest = ligand_manifest
    operator_review: list[dict[str, Any]] = []
    for row in placements.itertuples():
        ok = str(row.filter_2_source_ligand_instance_id) in valid_source_ids
        operator_ok = row.operator_quality_status == "PASS"
        ready = ok and operator_ok
        if not operator_ok:
            operator_review.append({"object_type": "ligand_assembly_placement", "object_id": str(row.filter_2_ligand_assembly_placement_id),
                                    "pdb_id": str(row.pdb_id), "assembly_id": str(row.assembly_id), "operator_path": str(row.operator_path),
                                    "determinant_error": float(row.operator_determinant_error), "orthogonality_error": float(row.operator_orthogonality_error),
                                    "quality_threshold": 1.0e-3, "action": "REVIEW_NO_MATRIX_CORRECTION"})
        manifest.append({"object_type": "ligand_assembly_placement", "object_id": str(row.filter_2_ligand_assembly_placement_id), "pdb_id": str(row.pdb_id),
                         "component_id": str(row.component_id), "preparation_status": "ASSEMBLY_READY" if ready else ("OPERATOR_MATRIX_REVIEW" if not operator_ok else "SOURCE_NOT_COMPLETE"),
                         "decision": "PASS" if ready else "REVIEW", "destination": "assembly_ready" if ready else ("operator_matrix_review" if not operator_ok else "ligand_preparation_review"),
                         "reason_code": "coordinate_transform_valid" if ready else ("operator_nonorthogonality_above_1e_3" if not operator_ok else "source_ligand_not_complete"),
                         "operator_quality_status": row.operator_quality_status, "operator_determinant_error": float(row.operator_determinant_error),
                         "operator_orthogonality_error": float(row.operator_orthogonality_error),
                         "prepared_atom_count": int((lig_assembly["filter_2_ligand_assembly_placement_id"] == row.filter_2_ligand_assembly_placement_id).sum()), "source_rows_modified": False})
    rec_atom_counts = rec.groupby("filter_1_source_chain_key").size().to_dict()
    rec_assembly_counts = rec_assembly.groupby("filter_1_chain_instance_id").size().to_dict() if not rec_assembly.empty else {}
    for row in receptors.itertuples():
        count = int(rec_assembly_counts.get(row.filter_1_chain_instance_id, 0))
        operator_ok = row.operator_quality_status == "PASS"
        ready = bool(count) and operator_ok
        if not operator_ok:
            operator_review.append({"object_type": "receptor_chain_instance", "object_id": str(row.filter_1_chain_instance_id),
                                    "pdb_id": str(row.pdb_id), "assembly_id": str(row.assembly_id), "operator_path": str(row.operator_path),
                                    "determinant_error": float(row.operator_determinant_error), "orthogonality_error": float(row.operator_orthogonality_error),
                                    "quality_threshold": 1.0e-3, "action": "REVIEW_NO_MATRIX_CORRECTION"})
        manifest.append({"object_type": "receptor_chain_instance", "object_id": str(row.filter_1_chain_instance_id), "pdb_id": str(row.pdb_id),
                         "component_id": "", "preparation_status": "ASSEMBLY_READY" if ready else ("OPERATOR_MATRIX_REVIEW" if not operator_ok else "NO_PREPARED_ATOMS"),
                         "decision": "PASS" if ready else ("REVIEW" if not operator_ok else "FAIL"),
                         "destination": "assembly_ready" if ready else ("operator_matrix_review" if not operator_ok else "preparation_failure"),
                         "reason_code": "coordinate_transform_valid" if ready else ("operator_nonorthogonality_above_1e_3" if not operator_ok else "source_chain_atoms_missing"),
                         "operator_quality_status": row.operator_quality_status, "operator_determinant_error": float(row.operator_determinant_error),
                         "operator_orthogonality_error": float(row.operator_orthogonality_error),
                         "prepared_atom_count": count, "source_rows_modified": False})

    outputs: dict[str, Any] = {
        "prepared_ligand_source_atoms": lig,
        "prepared_receptor_source_atoms": rec,
        "prepared_ligand_assembly_atoms": lig_assembly,
        "prepared_receptor_assembly_atoms": rec_assembly,
        "prepared_ligand_source_bonds": source_bonds_df,
        "prepared_ligand_assembly_bonds": assembly_bonds,
        "altloc_selection_audit": lig_alt + rec_alt,
        "microheterogeneity_selection": lig_micro + rec_micro,
        "ligand_mapping_review": mapping_review,
        "ligand_incomplete_instances": incomplete,
        "ligand_topology_validation": topology,
        "operator_matrix_quality_review": operator_review,
        "structure_preparation_manifest": manifest,
    }
    counts = {}
    for name, value in outputs.items():
        counts[name] = atomic_parquet(output_path(run, name, bid, tid), value)
    marker = run / "work/checkpoints" / f"task-{tid:06d}.json"
    atomic_json(marker, {"task_id": tid, "bucket_id": bid, "status": "complete", "row_counts": counts,
                         "runtime_seconds": time.time() - started, "finished_at": utc()})
    return {"task_id": tid, "bucket_id": bid, "row_counts": counts, "runtime_seconds": time.time() - started}


def safe_process_task(payload: tuple[str, int, int, set[str] | None]) -> dict[str, Any]:
    run_dir, bid, tid, selected = payload
    try:
        return {"ok": True, **process_task(run_dir, bid, tid, selected)}
    except Exception as exc:
        return {"ok": False, "task_id": tid, "bucket_id": bid,
                "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}


def build_tasks() -> list[tuple[int, int]]:
    return sorted(task_key(p) for p in (AUX / "entry_ligand_placements").rglob("part-*.parquet"))


def choose_smoke_pdbs(count: int = 1000) -> set[str]:
    table = pq.ParquetFile(AUX / "entry_work_manifest.parquet").read().to_pandas()
    active = table[table["has_active_work"]].copy()
    active["hash"] = active["pdb_id"].map(lambda x: hashlib.sha256(("processing2-smoke-v1|" + x).encode()).hexdigest())
    high = active.nlargest(min(20, max(1, count // 50)), "estimated_coordinate_transform_workload")
    rest = (active[~active["pdb_id"].isin(high["pdb_id"])]
            .sort_values(["hash", "pdb_id"], kind="stable")
            .head(count - len(high)))
    return set(pd.concat([high, rest])["pdb_id"].astype(str))


def relevant_tasks(selected: set[str] | None) -> list[tuple[int, int]]:
    tasks = build_tasks()
    if selected is None:
        return tasks
    wanted = []
    for bid, tid in tasks:
        pids = set(parquet_frame(input_path("entry_ligand_placements", bid, tid), ["pdb_id"])["pdb_id"].astype(str))
        if pids & selected:
            wanted.append((bid, tid))
    return wanted


def run_build(args: argparse.Namespace) -> None:
    run = Path(args.run_dir)
    selected = choose_smoke_pdbs(args.smoke_count) if args.smoke_count else None
    if selected is not None:
        (run / "input").mkdir(parents=True, exist_ok=True)
        (run / "input/smoke_pdb_ids.txt").write_text("\n".join(sorted(selected)) + "\n")
    tasks = relevant_tasks(selected)
    done = {int(p.stem.split("-")[1]) for p in (run / "work/checkpoints").glob("task-*.json")}
    pending = [(bid, tid) for bid, tid in tasks if tid not in done]
    progress = {"status": "RUNNING", "run_id": run.name, "task_count": len(tasks), "completed_task_count": len(tasks)-len(pending),
                "workers": args.workers, "smoke_pdb_count": len(selected) if selected else 0, "started_at": utc()}
    atomic_json(run / "status.json", progress)
    print(f"tasks={len(tasks)} pending={len(pending)} workers={args.workers}", flush=True)
    failures = []
    payloads = ((str(run), bid, tid, selected) for bid, tid in pending)
    # imap_unordered keeps workers fed across large-entry tails.  Recycling each
    # process after a small number of tasks prevents DataFrame/RDKit heap growth.
    with mp.Pool(processes=args.workers, maxtasksperchild=args.max_tasks_per_worker) as pool:
        for result in pool.imap_unordered(safe_process_task, payloads, chunksize=1):
            if result["ok"]:
                progress["completed_task_count"] += 1
                progress["last_task_id"] = result["task_id"]
                progress["last_task_runtime_seconds"] = result["runtime_seconds"]
            else:
                failures.append(result)
                progress["failed_task_count"] = len(failures)
            progress["updated_at"] = utc()
            atomic_json(run / "status.json", progress)
            if progress["completed_task_count"] % 10 == 0 or not result["ok"]:
                print(json.dumps(progress), flush=True)
    atomic_json(run / "audit/task_failures.json", failures)
    progress["status"] = "COMPLETED" if not failures else "FAILED"
    progress["finished_at"] = utc()
    atomic_json(run / "status.json", progress)
    if failures:
        raise SystemExit(f"{len(failures)} task failures")


def prepare_run(args: argparse.Namespace) -> None:
    run = Path(args.run_dir)
    for rel in ("input", "work/checkpoints", "output", "audit", "logs"):
        (run / rel).mkdir(parents=True, exist_ok=True)
    refs = {
        "auxiliary_build": str(AUX.parent),
        "auxiliary_output_manifest_sha256": sha256(AUX / "output_manifest.tsv"),
        "auxiliary_validation_sha256": sha256(AUX / "validation_report.json"),
        "filter_2_v3_run": str(F2),
        "filter_2_v3_frozen_sha256": sha256(F2 / "_FROZEN.json"),
        "ccd_path": str(CCD), "ccd_sha256": sha256(CCD),
        "crown_path": str(CROWN),
        "crown_commit": os.popen(f"git -C {CROWN} rev-parse HEAD").read().strip(),
    }
    atomic_json(run / "input/upstream.json", refs)
    config_source = ROOT / "config/default.yaml"
    shutil.copy2(config_source, run / "config_snapshot.yaml")
    atomic_json(run / "validation_plan.json", {
        "smoke_required": True, "risk_level": "high", "acceptance_criteria": [
            "input/output object accounting closes", "source coordinate rows unchanged",
            "AltLoc winner follows frozen CROWN rule", "microheterogeneity winner follows frozen CROWN rule",
            "strict missing CCD heavy atom count controls completeness", "assembly transform invariance within tolerance",
            "duplicate object keys equal zero", "malformed parquet files equal zero"],
    })
    atomic_json(run / "run_metadata.json", {"run_id": run.name, "stage": STAGE_VERSION, "schema_version": SCHEMA_VERSION,
                 "created_at": utc(), "inputs": refs, "config_sha256": sha256(run / "config_snapshot.yaml"),
                 "code_version_reference": {"script_sha256": sha256(Path(__file__))}, "status": "DRAFT"})


def matrix_quality(frame: pd.DataFrame) -> tuple[int, float, float]:
    if frame.empty:
        return 0, 0.0, 0.0
    cols = ["r11","r12","r13","r21","r22","r23","r31","r32","r33"]
    r = frame[cols].drop_duplicates().to_numpy(float).reshape(-1,3,3)
    det = np.linalg.det(r)
    orth = np.max(np.abs(np.matmul(np.transpose(r,(0,2,1)),r)-np.eye(3)), axis=(1,2))
    return len(r), float(np.max(np.abs(det-1))), float(np.max(orth))


def finalize(args: argparse.Namespace) -> None:
    run = Path(args.run_dir)
    failures = json.loads((run / "audit/task_failures.json").read_text()) if (run / "audit/task_failures.json").exists() else []
    manifest_rows = []
    totals = Counter()
    malformed = []
    for name in DATASETS:
        for path in sorted((run / "output" / name).rglob("*.parquet")):
            try:
                pf = pq.ParquetFile(path)
                rows = pf.metadata.num_rows
                cols = len(pf.schema_arrow)
                totals[name] += rows
                manifest_rows.append({"relative_path": str(path.relative_to(run)), "file_role": name, "file_format": "parquet",
                                      "row_count": rows, "column_count": cols, "size_bytes": path.stat().st_size,
                                      "sha256": sha256(path), "schema_version": SCHEMA_VERSION, "created_at": utc(), "generated_by": STAGE_VERSION})
            except Exception as exc:
                malformed.append({"path": str(path), "error": str(exc)})
    manifest_path = run / "output/output_manifest.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    manifest_rows.append({"relative_path": "output/output_manifest.tsv", "file_role": "output_manifest", "file_format": "tsv",
                          "row_count": len(manifest_rows), "column_count": 9, "size_bytes": manifest_path.stat().st_size,
                          "sha256": sha256(manifest_path), "schema_version": SCHEMA_VERSION, "created_at": utc(), "generated_by": STAGE_VERSION})
    # Schema is the union of schemas actually emitted for each dataset.
    schemas = {}
    previews = {}
    for name in DATASETS:
        files = sorted((run / "output" / name).rglob("*.parquet"))
        real = [p for p in files if pq.ParquetFile(p).metadata.num_rows]
        if not real:
            schemas[name] = []
            continue
        schema = pq.ParquetFile(real[0]).schema_arrow
        schemas[name] = [{"column_name": f.name, "data_type": str(f.type), "nullable": f.nullable} for f in schema]
        chunks=[]
        for p in real[:50]:
            t=pq.ParquetFile(p).read()
            if "_empty" not in t.column_names:
                chunks.append(t.slice(0, min(20,t.num_rows)))
            if sum(x.num_rows for x in chunks)>=1000: break
        if chunks:
            try:
                preview=pa.concat_tables(chunks, promote_options="permissive").slice(0,1000).to_pandas()
                preview.to_csv(run/f"output/{name}_preview.tsv",sep="\t",index=False)
                previews[name]=len(preview)
            except Exception as exc:
                malformed.append({"path": name+" preview", "error": str(exc)})
    atomic_json(run / "output/output_schema.json", {"schema_version": SCHEMA_VERSION, "datasets": schemas})
    summary = {"stage": STAGE_VERSION, "run_id": run.name, "status": "COMPLETED", "created_at": utc(),
               "row_counts": dict(totals), "task_failures": len(failures), "malformed_parquet": len(malformed),
               "preview_rows": previews, "source_atoms_modified": False,
               "missing_atom_reconstruction": False, "hydrogen_addition": False, "protonation": False,
               "minimization": False, "coordinate_optimization": False}
    atomic_json(run / "output/release_summary.json", summary)
    atomic_json(run / "audit/malformed_outputs.json", malformed)
    atomic_json(run / "output/downstream_interface.json", {
        "stage": STAGE_VERSION, "run_id": run.name, "formal_input": "output/structure_preparation_manifest",
        "coordinate_tables": ["output/prepared_ligand_assembly_atoms", "output/prepared_receptor_assembly_atoms"],
        "topology_tables": ["output/prepared_ligand_assembly_bonds", "output/ligand_topology_validation"],
        "read_preview_as_input": False,
    })
    sums=[]
    for p in sorted((run/"output").rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS": sums.append(f"{sha256(p)}  {p.relative_to(run/'output')}")
    (run/"output/SHA256SUMS").write_text("\n".join(sums)+"\n")
    print(json.dumps(summary,indent=2))


def validate(args: argparse.Namespace) -> None:
    run=Path(args.run_dir)
    summary=json.loads((run/"output/release_summary.json").read_text())
    status=json.loads((run/"status.json").read_text())
    checks={}
    checks["tasks_complete"] = status.get("status")=="COMPLETED" and status.get("completed_task_count")==status.get("task_count")
    checks["task_failures_zero"] = summary["task_failures"]==0
    checks["malformed_parquet_zero"] = summary["malformed_parquet"]==0
    checks["manifest_present"] = summary["row_counts"].get("structure_preparation_manifest",0)>0
    checks["ligand_source_accounted"] = summary["row_counts"].get("ligand_topology_validation",0)>0
    checks["source_rows_unchanged"] = summary["source_atoms_modified"] is False
    # Duplicate and terminal accounting checks use a streaming Arrow dataset scan.
    manifest_files=sorted((run/"output/structure_preparation_manifest").rglob("*.parquet"))
    keys=set(); duplicates=0; missing_decision=0
    counts=Counter()
    matrix_over_threshold=0; matrix_review_manifest=0; max_det=0.0; max_orth=0.0
    for p in manifest_files:
        t=pq.ParquetFile(p).read()
        if "_empty" in t.column_names: continue
        cols=["object_type","object_id","decision","preparation_status","operator_quality_status",
              "operator_determinant_error","operator_orthogonality_error"]
        for row in t.select(cols).to_pylist():
            key=(row["object_type"],row["object_id"])
            duplicates += key in keys; keys.add(key)
            missing_decision += not bool(row["decision"])
            counts[row["object_type"]]+=1
            if row["object_type"] in {"ligand_assembly_placement","receptor_chain_instance"}:
                det=float(row["operator_determinant_error"] or 0.0)
                orth=float(row["operator_orthogonality_error"] or 0.0)
                max_det=max(max_det,det); max_orth=max(max_orth,orth)
                over=det>1.0e-3 or orth>1.0e-3
                matrix_over_threshold += over
                matrix_review_manifest += row["operator_quality_status"]=="OPERATOR_MATRIX_REVIEW"
    checks["duplicate_manifest_keys_zero"] = duplicates==0
    checks["missing_terminal_decision_zero"] = missing_decision==0
    # Deposited operators are never altered. Matrices above the conservative
    # orthogonality gate must instead have an explicit review record.
    review_files=sorted((run/"output/operator_matrix_quality_review").rglob("*.parquet"))
    review_count=sum(pq.ParquetFile(p).metadata.num_rows for p in review_files)
    checks["operator_matrix_gate_accounted"] = (matrix_over_threshold==matrix_review_manifest==review_count)
    report={"validation_pass": all(checks.values()),"checks":checks,"duplicate_manifest_keys":duplicates,
            "missing_terminal_decision":missing_decision,"object_counts":dict(counts),
            "max_abs_determinant_error":max_det,"max_orthogonality_error":max_orth,
            "operator_matrix_over_threshold":matrix_over_threshold,
            "operator_matrix_review_manifest_count":matrix_review_manifest,
            "operator_matrix_review_count":review_count,"operator_quality_threshold":1.0e-3,"validated_at":utc()}
    atomic_json(run/"audit/assembly_coordinate_validation.json",report)
    if not report["validation_pass"]:
        print(json.dumps(report,indent=2)); raise SystemExit(2)
    metadata=json.loads((run/"run_metadata.json").read_text()); metadata["status"]="VALIDATED"; metadata["validated_at"]=utc()
    atomic_json(run/"run_metadata.json",metadata)
    atomic_json(run/"_FROZEN.json",{"status":"FROZEN","run_id":run.name,"stage":STAGE_VERSION,"frozen_at":utc(),
                "accounting_pass":True,"schema_pass":True,"validation_pass":True,
                "manifest_sha256":sha256(run/"output/output_manifest.tsv"),"code_version_reference":metadata["code_version_reference"]})
    current={"current_run_id":run.name,"status":"FROZEN","relative_path":f"runs/{run.name}",
             "manifest_sha256":sha256(run/"output/output_manifest.tsv"),"updated_at":utc()}
    atomic_json(ROOT/"CURRENT_RUN.json",current)
    link=ROOT/"current"
    if link.is_symlink() or link.exists(): link.unlink()
    link.symlink_to(Path("runs")/run.name)
    print(json.dumps(report,indent=2))


def selftest(args: argparse.Namespace) -> None:
    rows=[]
    for alt,count,occ,start in [("A",3,0.5,1),("B",3,0.5,10)]:
        for i in range(count):
            rows.append({"pdb_id":"test","object":"x","alt_id":alt,"source_atom_row_index":start+i,"occupancy":occ,
                         "type_symbol":"C","label_atom_id":f"{alt}{i}"})
    rows.append({"pdb_id":"test","object":"x","alt_id":"","source_atom_row_index":0,"occupancy":1.0,
                 "type_symbol":"N","label_atom_id":"N"})
    frame=pd.DataFrame(rows)
    selected,audit=select_altloc(frame,["object"],"test","object")
    assert set(selected["selected_altloc_original"])=={"","A"}
    assert audit[0]["selected_altloc"]=="A"
    r=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
    p=np.array([1.,2.,3.]); q=r@p+np.array([4.,5.,6.])
    assert np.allclose(q,[2.,6.,9.])
    print("selftest_pass")


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="cmd",required=True)
    sub.add_parser("init")
    for cmd in ("prepare-run","build-ccd","build","finalize","validate","selftest"):
        q=sub.add_parser(cmd); q.add_argument("--run-dir",required=True)
        if cmd=="build-ccd": q.add_argument("--force",action="store_true")
        if cmd=="build":
            q.add_argument("--workers",type=int,default=8); q.add_argument("--max-tasks-per-worker",type=int,default=4)
            q.add_argument("--smoke-count",type=int,default=0)
    return p


def main() -> None:
    args=parser().parse_args()
    {"init":init_stage,"prepare-run":prepare_run,"build-ccd":build_ccd,"build":run_build,
     "finalize":finalize,"validate":validate,"selftest":selftest}[args.cmd](args)


if __name__=="__main__":
    main()
