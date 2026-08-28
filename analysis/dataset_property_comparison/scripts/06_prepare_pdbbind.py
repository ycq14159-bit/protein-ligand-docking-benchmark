#!/usr/bin/env python3
"""Build the PDBbind v2020 comparison table from the official archives.

The archives remain outside Git.  One released ligand SDF is read per PDBbind
complex; no structure archive is redownloaded from another source.
"""

from __future__ import annotations

import argparse
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def parse_index(index_archive: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(index_archive, "r:gz") as tf:
        member = tf.getmember("index/INDEX_general_PL_data.2020")
        text = tf.extractfile(member).read().decode("utf-8", errors="replace")
    pattern = re.compile(
        r"^(?P<pdb>[0-9][A-Za-z0-9]{3})\s+(?P<resolution>\S+)\s+"
        r"(?P<year>\d{4})\s+(?P<pk>\S+)\s+(?P<affinity>\S+)\s+"
        r"//\s+(?P<reference>\S+)\s+\((?P<ligand>.*)\)\s*$"
    )
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        item = match.groupdict()
        try:
            item["resolution"] = float(item["resolution"])
        except ValueError:
            item["resolution"] = np.nan
        rows.append(item)
    result = pd.DataFrame(rows)
    if len(result) != 19443 or result["pdb"].nunique() != 19443:
        raise RuntimeError(
            f"Unexpected PDBbind index population: rows={len(result)}, "
            f"unique_pdb={result['pdb'].nunique()}"
        )
    return result


def mol_properties(block: str, source_format: str):
    try:
        if source_format == "SDF":
            mol = Chem.MolFromMolBlock(
                block, sanitize=True, removeHs=False, strictParsing=False
            )
        else:
            mol = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
        if mol is None:
            return (None, None, None, None, None, None, f"{source_format}_PARSE_FAILED")
        heavy = Chem.RemoveHs(mol)
        smiles = Chem.MolToSmiles(heavy, canonical=True, isomericSmiles=True)
        return (
            smiles,
            smiles,
            int(heavy.GetNumHeavyAtoms()),
            int(rdMolDescriptors.CalcNumRotatableBonds(heavy)),
            int(rdMolDescriptors.CalcNumHBA(heavy)),
            float(QED.qed(heavy)),
            "OK",
        )
    except Exception:
        return (None, None, None, None, None, None, f"{source_format}_PARSE_FAILED")


def read_ligands(archives: list[Path]) -> pd.DataFrame:
    rows = []
    for archive in archives:
        released_set = "refined" if "refined" in archive.name else "general-minus-refined"
        ligand_records = {}
        with tarfile.open(archive, "r:gz") as tf:
            for member in tf:
                if not member.isfile() or not member.name.endswith(("_ligand.sdf", "_ligand.mol2")):
                    continue
                pdb_id = Path(member.name).parent.name.lower()
                handle = tf.extractfile(member)
                block = handle.read().decode("utf-8", errors="replace")
                if member.name.endswith("_ligand.mol2"):
                    ligand_records.setdefault(pdb_id, {})["MOL2"] = (
                        member.name, mol_properties(block, "MOL2")
                    )
                    continue
                ligand_records.setdefault(pdb_id, {})["SDF"] = (
                    member.name, mol_properties(block, "SDF")
                )
        for pdb_id, record in ligand_records.items():
            source_member, props = record.get("SDF", (None, (None,) * 6 + ("SDF_MISSING",)))
            source_format = "SDF"
            if props[-1] != "OK" and "MOL2" in record:
                    source_member, props = record["MOL2"]
                    source_format = "MOL2_FALLBACK"
            if props[-1] != "OK" and source_format == "MOL2_FALLBACK":
                props = (*props[:-1], "BOTH_OFFICIAL_FORMATS_PARSE_FAILED")
            rows.append((pdb_id, released_set, source_member, source_format, *props))
    result = pd.DataFrame(rows, columns=[
        "pdb_id", "pdbbind_subset", "archive_member", "descriptor_source_format", "source_smiles",
        "canonical_isomeric_smiles", "heavy_atoms", "rotatable_bonds", "hba",
        "qed", "descriptor_status",
    ])
    if len(result) != 19443 or result["pdb_id"].nunique() != 19443:
        raise RuntimeError(
            f"Unexpected archive population: rows={len(result)}, "
            f"unique_pdb={result['pdb_id'].nunique()}"
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.external_root / "pdbbind_v2020"
    index_path = root / "PDBbind_v2020_plain_text_index.tar.gz"
    archives = [
        root / "PDBbind_v2020_refined.tar.gz",
        root / "PDBbind_v2020_other_PL.tar.gz",
    ]
    index = parse_index(index_path).rename(columns={"pdb": "pdb_id", "ligand": "ccd_id"})
    ligands = read_ligands(archives)
    merged = index.merge(ligands, on="pdb_id", how="outer", validate="one_to_one", indicator=True)
    if (merged["_merge"] != "both").any():
        raise RuntimeError(merged["_merge"].value_counts().to_dict())
    merged = merged.drop(columns="_merge")
    merged["dataset"] = "PDBbind"
    merged["dataset_version"] = "v2020-official"
    merged["source_entry_id"] = merged["pdb_id"]
    merged["ligand_id"] = merged["pdb_id"]
    merged["experimental_method"] = "X-RAY DIFFRACTION"
    merged["rsr"] = np.nan
    merged["rscc"] = np.nan
    merged["resolution_status"] = np.where(merged["resolution"].notna(), "AVAILABLE", "NO_DATA")
    merged["rsr_status"] = "VALIDATION_MAPPING_REQUIRED"
    merged["rscc_status"] = "VALIDATION_MAPPING_REQUIRED"
    merged["source_file"] = ";".join(str(p) for p in [index_path, *archives])
    merged["source_row_id"] = np.arange(len(merged), dtype=np.int64)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
