#!/usr/bin/env python3
"""Build and audit the Q-BioLiP snapshot through 2026-06-26.

The script preserves all weekly rows and records a provisional small-molecule
classification. It does not silently discard ambiguous ligand classes.
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import pandas as pd

POLYMER_LABELS = {"DNA", "RNA", "k-mer", "PEPTIDE", "PROTEIN"}
METAL_OR_SINGLE_ION = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "AL", "GA",
    "IN", "TL", "SN", "PB", "BI", "SC", "TI", "V", "CR", "MN", "FE",
    "CO", "NI", "CU", "ZN", "Y", "ZR", "NB", "MO", "TC", "RU", "RH",
    "PD", "AG", "CD", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU",
    "HG", "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB", "DY", "HO",
    "ER", "TM", "YB", "LU", "F", "CL", "BR", "IOD",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()

    ligand_rows = json.loads((args.root / "ligand.json").read_text(encoding="utf-8"))
    ligand_ids = {str(row["ligid"]).strip() for row in ligand_rows}
    weekly = []
    for filename in sorted(glob.glob(str(args.root / "weekly" / "Q-BioLiP-*.csv"))):
        frame = pd.read_csv(filename, dtype=str).fillna("")
        frame["source_file"] = Path(filename).name
        weekly.append(frame)
    frame = pd.concat(weekly, ignore_index=True)
    frame["ligand_dictionary_match"] = frame["Ligand ID"].isin(ligand_ids)
    frame["explicit_polymer_label"] = frame["Ligand ID"].str.upper().isin(POLYMER_LABELS)
    frame["single_atom_metal_or_ion_label"] = frame["Ligand ID"].str.upper().isin(METAL_OR_SINGLE_ION)
    frame["provisional_small_molecule"] = (
        frame["ligand_dictionary_match"]
        & ~frame["explicit_polymer_label"]
        & ~frame["single_atom_metal_or_ion_label"]
    )
    frame["snapshot_through"] = "2026-06-26"
    frame.to_parquet(args.output, index=False)

    base = pd.read_csv(args.root / "PL_annotation_base_before_20260102.csv", dtype=str)
    audit = {
        "snapshot_through": "2026-06-26",
        "base_official_PL_rows": int(len(base)),
        "weekly_files": len(weekly),
        "weekly_all_rows": int(len(frame)),
        "weekly_unique_qbiolip_id": int(frame["Q-BioLiP ID"].nunique()),
        "weekly_ligand_dictionary_match": int(frame["ligand_dictionary_match"].sum()),
        "weekly_explicit_polymer": int(frame["explicit_polymer_label"].sum()),
        "weekly_single_atom_metal_or_ion": int(frame["single_atom_metal_or_ion_label"].sum()),
        "weekly_provisional_small_molecule": int(frame["provisional_small_molecule"].sum()),
        "status": "REVIEW_REQUIRED",
        "note": "Population is preserved; provisional class is not yet the formal adapter membership rule.",
    }
    args.audit.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
