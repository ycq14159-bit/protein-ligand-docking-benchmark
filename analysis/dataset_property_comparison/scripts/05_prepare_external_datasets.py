#!/usr/bin/env python3
"""Prepare CROWN, HiQBind, PLINDER and Ours using one RDKit definition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

OUT_COLUMNS = [
    "dataset", "dataset_version", "source_entry_id", "pdb_id", "ligand_id",
    "ccd_id", "ligand_chain", "ligand_resnum", "ligand_altloc", "source_smiles",
    "canonical_isomeric_smiles", "heavy_atoms", "rotatable_bonds", "hba", "qed",
    "experimental_method", "resolution", "rsr", "rscc", "descriptor_status",
    "resolution_status", "rsr_status", "rscc_status", "source_file", "source_row_id",
]


def descriptors(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return (None, None, None, None, None, "MISSING_SMILES")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (None, None, None, None, None, "PARSE_FAILED")
        return (
            Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            int(mol.GetNumHeavyAtoms()),
            int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
            int(rdMolDescriptors.CalcNumHBA(mol)),
            float(QED.qed(mol)),
            "OK",
        )
    except Exception:
        return (None, None, None, None, None, "PARSE_FAILED")


def add_descriptors(frame, smiles_column):
    smiles = frame[smiles_column].tolist()
    keys = [value if isinstance(value, str) else None for value in smiles]
    cache = {value: descriptors(value) for value in set(keys)}
    values = [cache[value] for value in keys]
    desc = pd.DataFrame(values, columns=[
        "canonical_isomeric_smiles", "heavy_atoms", "rotatable_bonds", "hba", "qed",
        "descriptor_status",
    ])
    for col in desc:
        frame[col] = desc[col].values
    return frame


def status(series):
    return np.where(series.notna(), "AVAILABLE", "NO_DATA")


def finish(frame, path):
    for col in OUT_COLUMNS:
        if col not in frame:
            frame[col] = None
    path.parent.mkdir(parents=True, exist_ok=True)
    frame[OUT_COLUMNS].to_parquet(path, index=False)


def crown(external, prepared):
    source = external / "crown_202606" / "CROWN_metadata.parquet"
    d = pd.read_parquet(source)
    d = d.rename(columns={
        "basename": "source_entry_id", "pdb_id": "pdb_id", "lig_name": "ccd_id",
        "SMILES": "source_smiles", "resolution": "resolution", "ligand_rsr": "rsr",
        "ligand_rscc": "rscc",
    })
    d["dataset"] = "CROWN"
    d["dataset_version"] = "2026-06"
    d["ligand_id"] = d["source_entry_id"]
    d["experimental_method"] = "X-RAY DIFFRACTION"
    d["source_file"] = str(source)
    d["source_row_id"] = np.arange(len(d), dtype=np.int64)
    d = add_descriptors(d, "source_smiles")
    d["resolution_status"] = status(d["resolution"])
    d["rsr_status"] = status(d["rsr"])
    d["rscc_status"] = status(d["rscc"])
    finish(d, prepared / "crown_properties.parquet")


def hiqbind(external, prepared):
    source = external / "hiqbind_v3" / "hiqbind_sm_metadata.csv"
    d = pd.read_csv(source)
    d["source_row_id"] = np.arange(len(d), dtype=np.int64)
    d["source_entry_id"] = (
        d["PDBID"].astype(str).str.lower() + ":" + d["Ligand Chain"].astype(str) + ":" +
        d["Ligand Residue Number"].astype(str) + ":row" + d["source_row_id"].astype(str)
    )
    d = d.rename(columns={
        "PDBID": "pdb_id", "Ligand Name": "ccd_id", "Ligand Chain": "ligand_chain",
        "Ligand Residue Number": "ligand_resnum", "Ligand SMILES": "source_smiles",
        "Resolution": "resolution",
    })
    d["dataset"] = "HiQBind"
    d["dataset_version"] = "Figshare-v3-small-molecule"
    d["ligand_id"] = d["source_entry_id"]
    d["experimental_method"] = None
    d["rsr"] = np.nan
    d["rscc"] = np.nan
    d["source_file"] = str(source)
    d = add_descriptors(d, "source_smiles")
    d["resolution_status"] = status(d["resolution"])
    d["rsr_status"] = "VALIDATION_MAPPING_REQUIRED"
    d["rscc_status"] = "VALIDATION_MAPPING_REQUIRED"
    finish(d, prepared / "hiqbind_properties.parquet")


def pick_instance_metric(instance, asym_id, keys, values):
    if keys is None or values is None:
        return np.nan
    target = f"{int(instance)}.{asym_id}"
    try:
        matches = [i for i, key in enumerate(keys) if str(key) == target]
        if len(matches) != 1:
            return np.nan
        value = values[matches[0]]
        return float(value) if value is not None else np.nan
    except Exception:
        return np.nan


def plinder(external, prepared):
    source = external / "plinder_2024-06_v2" / "annotation_table.parquet"
    cols = [
        "entry_pdb_id", "entry_determination_method", "entry_resolution", "system_id",
        "ligand_id", "ligand_ccd_code", "ligand_asym_id", "ligand_auth_id",
        "ligand_instance", "ligand_is_proper", "ligand_rdkit_canonical_smiles",
        "ligand_smiles", "system_ligand_chains_asym_id",
        "system_ligand_chains_validation_average_rsr",
        "system_ligand_chains_validation_average_rscc",
    ]
    d = pd.read_parquet(source, columns=cols)
    d = d[d["ligand_is_proper"]].copy().reset_index(drop=True)
    d["source_row_id"] = np.arange(len(d), dtype=np.int64)
    shared = zip(
        d["ligand_instance"], d["ligand_asym_id"],
        d["system_ligand_chains_asym_id"],
        d["system_ligand_chains_validation_average_rsr"],
    )
    d["rsr"] = [pick_instance_metric(*row) for row in shared]
    shared = zip(
        d["ligand_instance"], d["ligand_asym_id"],
        d["system_ligand_chains_asym_id"],
        d["system_ligand_chains_validation_average_rscc"],
    )
    d["rscc"] = [pick_instance_metric(*row) for row in shared]
    d["source_smiles"] = d["ligand_rdkit_canonical_smiles"].where(
        d["ligand_rdkit_canonical_smiles"].notna(), d["ligand_smiles"]
    )
    d = d.rename(columns={
        "entry_pdb_id": "pdb_id", "entry_determination_method": "experimental_method",
        "entry_resolution": "resolution", "ligand_ccd_code": "ccd_id",
        "ligand_auth_id": "ligand_chain",
    })
    d["dataset"] = "PLINDER"
    d["dataset_version"] = "2024-06/v2"
    d["source_entry_id"] = d["ligand_id"]
    d["source_file"] = str(source)
    d = add_descriptors(d, "source_smiles")
    d["resolution_status"] = status(d["resolution"])
    d["rsr_status"] = status(d["rsr"])
    d["rscc_status"] = status(d["rscc"])
    finish(d, prepared / "plinder_properties.parquet")


def ours(repo, prepared):
    source = repo / "analysis" / "dataset_property_comparison" / "prepared" / "ours_properties.parquet"
    d = pd.read_parquet(source)
    d["dataset_version"] = "20260826_Filter4_PASS"
    d["source_entry_id"] = d["pair_id"]
    d["ligand_id"] = d["ligand_assembly_placement_id"]
    d["source_smiles"] = d["normalized_ccd_isomeric_smiles"]
    d["canonical_isomeric_smiles"] = d["normalized_ccd_isomeric_smiles"]
    d["descriptor_status"] = d["descriptor_parse_status"]
    d["resolution_status"] = status(d["resolution"])
    d["rsr_status"] = status(d["rsr"])
    d["rscc_status"] = status(d["rscc"])
    d["source_file"] = str(source)
    d["source_row_id"] = np.arange(len(d), dtype=np.int64)
    finish(d, prepared / "ours_properties_harmonized.parquet")


def biolip2(external, prepared):
    root = external / "biolip2_20260626"
    ligand_rows = json.loads((root / "ligand.json").read_text(encoding="utf-8"))
    ligand_table = pd.DataFrame(ligand_rows).rename(
        columns={"ligid": "ccd_id", "smiles": "source_smiles"}
    )
    ligand_table["source_smiles"] = ligand_table["source_smiles"].map(
        lambda value: value if isinstance(value, str) else None
    )
    if ligand_table["ccd_id"].duplicated().any():
        raise RuntimeError("BioLiP2 ligand dictionary contains duplicate ligand IDs")
    ligand_table["ccd_id"] = ligand_table["ccd_id"].astype(str).str.strip()
    ligand_ids = set(ligand_table["ccd_id"])
    non_small_labels = {
        "DNA", "RNA", "DNA+RNA", "DNAH", "RNAH", "NUC", "NUCH", "K-MER",
        "KMER", "PEPTIDE", "PROTEIN", "LI", "NA", "K", "RB", "CS", "MG",
        "CA", "SR", "BA", "AL", "GA", "IN", "TL", "SN", "PB", "BI", "SC",
        "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "Y", "ZR",
        "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "HF", "TA", "W",
        "RE", "OS", "IR", "PT", "AU", "HG", "LA", "CE", "PR", "ND", "SM",
        "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "F", "CL",
        "BR", "IOD",
    }

    base_file = root / "PL_annotation_base_before_20260102.csv"
    base = pd.read_csv(base_file, dtype=str)
    base_ccd = base["Ligand_ID"].astype(str).str.strip()
    base = base.loc[
        base_ccd.isin(ligand_ids) & ~base_ccd.str.upper().isin(non_small_labels)
    ].copy()
    base["source_entry_id"] = "BASE:" + base["Ligand_file"]
    base["pdb_id"] = base["Assembly_ID"].str[:4].str.lower()
    base["ligand_id"] = base["Ligand_file"]
    base["ccd_id"] = base["Ligand_ID"].str.strip()
    prefix = base["Assembly_ID"] + "_" + base["ccd_id"] + "_"
    base["ligand_chain"] = [
        value[len(pre):] if value.startswith(pre) else None
        for value, pre in zip(base["Ligand_file"], prefix)
    ]
    base["resolution"] = np.nan
    base["source_file"] = str(base_file)
    base["snapshot_segment"] = "OFFICIAL_BASE_PROTEIN_SMALL_MOLECULE"

    weekly_frames = []
    for weekly_file in sorted((root / "weekly").glob("Q-BioLiP-*.csv")):
        weekly = pd.read_csv(weekly_file, dtype=str)
        weekly["weekly_source_file"] = str(weekly_file)
        weekly_frames.append(weekly)
    weekly = pd.concat(weekly_frames, ignore_index=True)
    weekly_ccd = weekly["Ligand ID"].astype(str).str.strip()
    keep = (
        weekly["Relevant"].astype(str).str.lower().eq("yes")
        & weekly_ccd.isin(ligand_ids)
        & ~weekly_ccd.str.upper().isin(non_small_labels)
    )
    weekly = weekly.loc[keep].copy()
    weekly["source_entry_id"] = "WEEKLY:" + weekly["Ligand Detail"]
    weekly["pdb_id"] = weekly["PDB ID"].str.lower()
    weekly["ligand_id"] = weekly["Ligand Detail"]
    weekly["ccd_id"] = weekly["Ligand ID"].str.strip()
    prefix = weekly["Assembly ID"] + "_" + weekly["ccd_id"] + "_"
    weekly["ligand_chain"] = [
        value[len(pre):] if value.startswith(pre) else None
        for value, pre in zip(weekly["Ligand Detail"], prefix)
    ]
    weekly["resolution"] = pd.to_numeric(weekly["Resolution (Å)"], errors="coerce")
    weekly["source_file"] = weekly["weekly_source_file"]
    weekly["snapshot_segment"] = "WEEKLY_RELEVANT_SMALL_MOLECULE"

    cols = [
        "source_entry_id", "pdb_id", "ligand_id", "ccd_id", "ligand_chain",
        "resolution", "source_file", "snapshot_segment",
    ]
    d = pd.concat([base[cols], weekly[cols]], ignore_index=True)
    if d["source_entry_id"].duplicated().any():
        raise RuntimeError("BioLiP2 snapshot source entry IDs are not unique")
    d = d.merge(ligand_table[["ccd_id", "source_smiles"]], on="ccd_id", how="left", validate="many_to_one")
    d["dataset"] = "BioLiP2/Q-BioLiP"
    d["dataset_version"] = "snapshot-through-2026-06-26"
    d["experimental_method"] = None
    d["rsr"] = np.nan
    d["rscc"] = np.nan
    d["source_row_id"] = np.arange(len(d), dtype=np.int64)
    d = add_descriptors(d, "source_smiles")
    d["resolution_status"] = status(d["resolution"])
    d["rsr_status"] = "VALIDATION_MAPPING_REQUIRED"
    d["rscc_status"] = "VALIDATION_MAPPING_REQUIRED"
    finish(d, prepared / "biolip2_properties.parquet")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--datasets", nargs="+", default=["crown", "hiqbind", "plinder", "ours"])
    args = parser.parse_args()
    for name in args.datasets:
        globals()[name](args.external_root, args.prepared) if name != "ours" else ours(args.repo_root, args.prepared)


if __name__ == "__main__":
    main()
