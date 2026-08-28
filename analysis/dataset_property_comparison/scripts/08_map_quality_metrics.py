#!/usr/bin/env python3
"""Map normalized wwPDB ligand validation records onto harmonized datasets."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def text(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def seq(value):
    value = text(value)
    if value.endswith(".0"):
        value = value[:-2]
    try:
        return str(int(value))
    except ValueError:
        return value


def load_validation(root: Path):
    entries = pd.concat(
        [pd.read_parquet(path) for path in sorted(root.glob("batches/*/entry_validation.parquet"))],
        ignore_index=True,
    )
    residues = pd.concat(
        [pd.read_parquet(path) for path in sorted(root.glob("batches/*/ligand_candidate_validation.parquet"))],
        ignore_index=True,
    )
    entries = entries.sort_values("pdb_id").drop_duplicates("pdb_id", keep="last")
    residues["pdb_key"] = residues["pdb_id"].map(text).str.lower()
    residues["ccd_key"] = residues["component_id"].map(text).str.upper()
    residues["chain_key"] = residues["auth_asym_id"].map(text)
    residues["resnum_key"] = residues["auth_seq_id"].map(seq)
    residues["icode_key"] = residues["insertion_code"].map(text)
    residues["alt_key"] = residues["alt_id"].map(text)
    return entries, residues


def load_validation_targeted(root: Path, pdb_ids: set[str], components: set[str]):
    entries = pd.concat(
        [pd.read_parquet(path) for path in sorted(root.glob("batches/*/entry_validation.parquet"))],
        ignore_index=True,
    ).sort_values("pdb_id").drop_duplicates("pdb_id", keep="last")
    selected = []
    columns = [
        "pdb_id", "model_id", "auth_asym_id", "label_asym_id", "auth_seq_id",
        "label_seq_id", "insertion_code", "component_id", "alt_id", "rsr", "rscc",
    ]
    for path in sorted(root.glob("batches/*/ligand_candidate_validation.parquet")):
        frame = pd.read_parquet(path, columns=columns)
        frame = frame[
            frame["pdb_id"].isin(pdb_ids) & frame["component_id"].isin(components)
        ]
        if len(frame):
            selected.append(frame)
    residues = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame(columns=columns)
    residues["pdb_key"] = residues["pdb_id"].map(text).str.lower()
    residues["ccd_key"] = residues["component_id"].map(text).str.upper()
    residues["chain_key"] = residues["auth_asym_id"].map(text)
    residues["resnum_key"] = residues["auth_seq_id"].map(seq)
    residues["icode_key"] = residues["insertion_code"].map(text)
    residues["alt_key"] = residues["alt_id"].map(text)
    return entries, residues


def apply_entry_resolution(frame, entries):
    resolution = entries.set_index("pdb_id")["resolution"]
    missing = frame["resolution"].isna()
    frame.loc[missing, "resolution"] = frame.loc[missing, "pdb_key"].map(resolution)
    frame["resolution_status"] = np.where(frame["resolution"].notna(), "AVAILABLE", "NO_DATA")


def map_candidates(frame, residues, left_keys, right_keys, eligible=None):
    frame = frame.copy().reset_index(drop=True)
    frame["_row_id"] = np.arange(len(frame), dtype=np.int64)
    if eligible is None:
        eligible = pd.Series(True, index=frame.index)
    left = frame.loc[eligible, ["_row_id", *left_keys]]
    right_columns = list(dict.fromkeys([*right_keys, "rsr", "rscc", "model_id", "alt_key"]))
    merged = left.merge(
        residues[right_columns], left_on=left_keys, right_on=right_keys,
        how="left", sort=False,
    )
    matched = merged[merged["model_id"].notna()].copy()
    counts = matched.groupby("_row_id").size()
    exact_ids = set(counts[counts == 1].index)
    multiple_ids = set(counts[counts > 1].index)
    status = pd.Series("VALIDATION_RESIDUE_NOT_FOUND", index=frame.index, dtype=object)
    status.loc[~eligible] = frame.loc[~eligible, "validation_mapping_status"]
    status.loc[list(exact_ids)] = "EXACT_INSTANCE_MATCH"
    status.loc[list(multiple_ids)] = "AMBIGUOUS_VALIDATION_INSTANCE"
    exact = matched[matched["_row_id"].isin(exact_ids)].set_index("_row_id")
    frame["rsr"] = np.nan
    frame["rscc"] = np.nan
    if len(exact):
        frame.loc[exact.index, "rsr"] = exact["rsr"]
        frame.loc[exact.index, "rscc"] = exact["rscc"]
    metrics_missing = status.eq("EXACT_INSTANCE_MATCH") & (frame["rsr"].isna() | frame["rscc"].isna())
    status.loc[metrics_missing] = "EXACT_INSTANCE_METRICS_MISSING"
    frame["validation_mapping_status"] = status
    frame["rsr_status"] = np.where(frame["rsr"].notna(), "AVAILABLE", status)
    frame["rscc_status"] = np.where(frame["rscc"].notna(), "AVAILABLE", status)
    return frame.drop(columns="_row_id")


def common_keys(frame):
    frame["pdb_key"] = frame["pdb_id"].map(text).str.lower()
    frame["ccd_key"] = frame["ccd_id"].map(text).str.upper()
    frame["chain_key"] = frame["ligand_chain"].map(text) if "ligand_chain" in frame else ""
    if "ligand_resnum" in frame:
        frame["resnum_key"] = frame["ligand_resnum"].map(seq)
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-root", type=Path, required=True)
    parser.add_argument("--validation-root", type=Path, required=True)
    args = parser.parse_args()
    hiq_output = args.prepared_root / "hiqbind_properties_quality.parquet"
    bio_output = args.prepared_root / "biolip2_properties_quality.parquet"
    if hiq_output.exists() and bio_output.exists():
        coordinate_map = pd.read_parquet(args.prepared_root / "pdbbind_ligand_instance_mapping.parquet")
        exact = coordinate_map[coordinate_map["mapping_status"].eq("EXACT_COORDINATE_MATCH")]
        entries, residues = load_validation_targeted(
            args.validation_root,
            set(exact["pdb_id"].astype(str).str.lower()),
            set(exact["component_id"].dropna().astype(str).str.upper()),
        )
    else:
        entries, residues = load_validation(args.validation_root)
    report_ok = set(entries.loc[entries.parse_status.eq("PARSE_SUCCESS"), "pdb_id"])

    # HiQBind: PDB + CCD + author chain + author residue number.
    output = hiq_output
    if not output.exists():
        path = args.prepared_root / "hiqbind_properties.parquet"
        hiq = common_keys(pd.read_parquet(path))
        hiq["validation_mapping_status"] = np.where(
            hiq["pdb_key"].isin(report_ok), "PENDING", "VALIDATION_REPORT_UNAVAILABLE"
        )
        apply_entry_resolution(hiq, entries)
        hiq = map_candidates(
            hiq, residues,
            ["pdb_key", "ccd_key", "chain_key", "resnum_key"],
            ["pdb_key", "ccd_key", "chain_key", "resnum_key"],
            eligible=hiq["validation_mapping_status"].eq("PENDING"),
        )
        hiq.to_parquet(output, index=False)

    # BioLiP2: released identity has PDB + CCD + ligand chain but no residue number.
    output = bio_output
    if not output.exists():
        path = args.prepared_root / "biolip2_properties.parquet"
        bio = common_keys(pd.read_parquet(path))
        bio["validation_mapping_status"] = np.where(
            bio["pdb_key"].isin(report_ok), "PENDING", "VALIDATION_REPORT_UNAVAILABLE"
        )
        apply_entry_resolution(bio, entries)
        bio = map_candidates(
            bio, residues,
            ["pdb_key", "ccd_key", "chain_key"],
            ["pdb_key", "ccd_key", "chain_key"],
            eligible=bio["validation_mapping_status"].eq("PENDING"),
        )
        bio.to_parquet(output, index=False)

    # PDBbind: first require unique heavy-atom coordinate mapping to a deposited residue.
    pdbbind = common_keys(pd.read_parquet(args.prepared_root / "pdbbind_properties.parquet"))
    coord = pd.read_parquet(args.prepared_root / "pdbbind_ligand_instance_mapping.parquet")
    coord = coord.rename(columns={
        "auth_asym_id": "mapped_chain", "auth_seq_id": "mapped_resnum",
        "insertion_code": "mapped_icode", "component_id": "mapped_ccd", "alt_id": "mapped_alt",
        "mapping_status": "coordinate_mapping_status",
    })
    pdbbind = pdbbind.merge(coord, on="pdb_id", how="left", validate="one_to_one")
    pdbbind["pdb_key"] = pdbbind["pdb_id"].map(text).str.lower()
    pdbbind["mapped_ccd_key"] = pdbbind["mapped_ccd"].map(text).str.upper()
    pdbbind["mapped_chain_key"] = pdbbind["mapped_chain"].map(text)
    pdbbind["mapped_resnum_key"] = pdbbind["mapped_resnum"].map(seq)
    pdbbind["mapped_icode_key"] = pdbbind["mapped_icode"].map(text)
    pdbbind["mapped_alt_key"] = pdbbind["mapped_alt"].map(text)
    pdbbind["validation_mapping_status"] = pdbbind["coordinate_mapping_status"]
    exact_coordinate = pdbbind["coordinate_mapping_status"].eq("EXACT_COORDINATE_MATCH")
    pdbbind.loc[exact_coordinate & ~pdbbind["pdb_key"].isin(report_ok), "validation_mapping_status"] = "VALIDATION_REPORT_UNAVAILABLE"
    eligible = exact_coordinate & pdbbind["pdb_key"].isin(report_ok)
    apply_entry_resolution(pdbbind, entries)
    pdbbind = map_candidates(
        pdbbind, residues,
        ["pdb_key", "mapped_ccd_key", "mapped_chain_key", "mapped_resnum_key", "mapped_icode_key", "mapped_alt_key"],
        ["pdb_key", "ccd_key", "chain_key", "resnum_key", "icode_key", "alt_key"],
        eligible=eligible,
    )
    pdbbind.to_parquet(args.prepared_root / "pdbbind_properties_quality.parquet", index=False)

    # Already released or frozen exact metrics are retained with explicit provenance.
    for name, provenance in [
        ("crown", "RELEASE_WWPDB_DERIVED_METRICS"),
        ("plinder", "RELEASE_EXACT_INSTANCE_WWPDB_METRICS"),
        ("ours_properties_harmonized", "FROZEN_EXACT_INSTANCE_WWPDB_METRICS"),
    ]:
        source = args.prepared_root / (name + "_properties.parquet" if name in {"crown", "plinder"} else name + ".parquet")
        frame = common_keys(pd.read_parquet(source))
        frame["validation_mapping_status"] = provenance
        frame["rsr_status"] = np.where(frame["rsr"].notna(), "AVAILABLE", "NO_DATA")
        frame["rscc_status"] = np.where(frame["rscc"].notna(), "AVAILABLE", "NO_DATA")
        frame.to_parquet(args.prepared_root / (name + "_quality.parquet"), index=False)


if __name__ == "__main__":
    main()
