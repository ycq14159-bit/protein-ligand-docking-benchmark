#!/usr/bin/env python3
"""Create formal population, descriptor, mapping, missingness and summary QC."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

FILES = {
    "PDBbind": "pdbbind_properties_quality.parquet",
    "HiQBind": "hiqbind_properties_quality.parquet",
    "BioLiP2": "biolip2_properties_quality.parquet",
    "PLINDER": "plinder_quality.parquet",
    "CROWN": "crown_quality.parquet",
    "Ours": "ours_properties_harmonized_quality.parquet",
}
LOCKED_RAW = {
    "PDBbind": 19443, "HiQBind": 31572, "BioLiP2": 710466,
    "PLINDER": 1357906, "CROWN": 141261, "Ours": 91860,
}
METRICS = ["heavy_atoms", "rotatable_bonds", "hba", "qed", "resolution", "rsr", "rscc"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-root", type=Path, required=True)
    parser.add_argument("--qc-root", type=Path, required=True)
    args = parser.parse_args()
    args.qc_root.mkdir(parents=True, exist_ok=True)
    populations, descriptors, missingness, summaries, mappings = [], [], [], [], []
    for dataset, filename in FILES.items():
        frame = pd.read_parquet(args.prepared_root / filename)
        key = "source_entry_id"
        row = {
            "dataset": dataset, "locked_raw_N": LOCKED_RAW[dataset],
            "formal_mode_b_N": len(frame), "unique_source_entry_id_N": frame[key].nunique(dropna=True),
            "duplicate_source_entry_id_N": int(frame.duplicated(key).sum()),
            "missing_source_entry_id_N": int(frame[key].isna().sum()),
            "unique_pdb_N": frame["pdb_id"].nunique(dropna=True),
            "descriptor_N": int(pd.to_numeric(frame["heavy_atoms"], errors="coerce").notna().sum()),
            "resolution_N": int(pd.to_numeric(frame["resolution"], errors="coerce").notna().sum()),
            "rsr_N": int(pd.to_numeric(frame["rsr"], errors="coerce").notna().sum()),
            "rscc_N": int(pd.to_numeric(frame["rscc"], errors="coerce").notna().sum()),
        }
        populations.append(row)
        for status, count in frame["descriptor_status"].value_counts(dropna=False).items():
            descriptors.append({"dataset": dataset, "descriptor_status": str(status), "N": int(count)})
        for status, count in frame["validation_mapping_status"].value_counts(dropna=False).items():
            mappings.append({"dataset": dataset, "validation_mapping_status": str(status), "N": int(count)})
        for metric in METRICS:
            series = pd.to_numeric(frame[metric], errors="coerce")
            available = series.dropna()
            missingness.append({
                "dataset": dataset, "metric": metric, "total_N": len(frame),
                "available_N": len(available), "missing_N": int(series.isna().sum()),
                "missing_percent": float(series.isna().mean() * 100),
            })
            if len(available):
                summaries.append({
                    "dataset": dataset, "metric": metric, "N": len(available),
                    "mean": float(available.mean()), "sd": float(available.std()),
                    "min": float(available.min()), "q1": float(available.quantile(.25)),
                    "median": float(available.median()), "q3": float(available.quantile(.75)),
                    "max": float(available.max()),
                })
    population = pd.DataFrame(populations)
    if (population["formal_mode_b_N"] != population["unique_source_entry_id_N"]).any():
        raise RuntimeError("At least one formal adapter has a non-unique source entry key")
    population.to_csv(args.qc_root / "population_qc.tsv", sep="\t", index=False)
    pd.DataFrame(descriptors).to_csv(args.qc_root / "descriptor_qc.tsv", sep="\t", index=False)
    pd.DataFrame(mappings).to_csv(args.qc_root / "validation_mapping_qc.tsv", sep="\t", index=False)
    pd.DataFrame(missingness).to_csv(args.qc_root / "missingness_qc.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.qc_root / "summary_statistics.tsv", sep="\t", index=False)
    (args.qc_root / "qc_summary.json").write_text(json.dumps({
        "status": "PASS", "datasets": len(FILES),
        "formal_total_rows_across_datasets": int(population.formal_mode_b_N.sum()),
    }, indent=2) + "\n")
    print(population.to_string(index=False))


if __name__ == "__main__":
    main()
