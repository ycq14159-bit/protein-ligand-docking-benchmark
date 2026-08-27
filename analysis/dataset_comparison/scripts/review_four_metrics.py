#!/usr/bin/env python3
"""Read-only focused review of four potentially non-comparable metrics."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--analysis-root", required=True, type=Path)
    args = ap.parse_args()
    root, ar = args.data_root, args.analysis_root
    qc = ar / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    entry_path = ar / "data" / "ours_entry_properties.parquet"
    cols = ["pair_id", "pdb_id", "assembly_id", "resolved_ccd_id",
            "source_ligand_instance_id", "observed_heavy_atom_count",
            "previous_terminal_route", "previous_reason_code",
            "historical_artifact_provenance_proxy", "one_heavy_atom_proxy",
            "unresolved_pocket_atom_count",
            "pocket_missing_backbone_heavy_atom_count",
            "direct_binding_missing_sidechain_heavy_atom_count",
            "nonbinding_pocket_missing_sidechain_heavy_atom_count"]
    entry = pq.read_table(entry_path, columns=cols).to_pandas()

    f2p = root / "filter_2_ligand_qualification_v4/runs/20260825_dual_source_strict_01/output/01_source_membership.tsv.gz"
    f2cols = ["source_ligand_instance_id", "resolved_ccd_id", "ccd_name", "ccd_type",
              "formula", "formal_charge", "observed_heavy_atom_count", "is_suspicious",
              "biolip_list_match", "overlap_class", "q_match_level", "q_relevance_class",
              "dual_source_status", "dual_source_pass", "previous_terminal_route",
              "previous_reason_code", "terminal_route", "reason_code"]
    f2 = pd.read_csv(f2p, sep="\t", usecols=f2cols, low_memory=False)
    f2 = f2[f2.source_ligand_instance_id.isin(set(entry.source_ligand_instance_id))]
    formal = entry[["pair_id", "pdb_id", "assembly_id", "source_ligand_instance_id"]].merge(
        f2, on="source_ligand_instance_id", how="left", validate="many_to_one")

    artifact_pattern = r"artifact|additive|solvent|simple_inorganic"
    art = formal[(formal.previous_terminal_route.fillna("").str.contains(artifact_pattern, case=False, regex=True)) |
                 (formal.previous_reason_code.fillna("").str.contains(artifact_pattern, case=False, regex=True))].copy()
    art.sort_values("pair_id").to_csv(qc / "artifact_proxy_3526_pair_audit.tsv.gz", sep="\t", index=False,
                                      compression={"method": "gzip", "mtime": 0}, lineterminator="\n")
    art_summary = (art.groupby(["previous_terminal_route", "previous_reason_code", "is_suspicious",
                                "overlap_class", "q_match_level", "q_relevance_class",
                                "dual_source_status", "dual_source_pass"], dropna=False)
                   .agg(entries=("pair_id", "size"), unique_sources=("source_ligand_instance_id", "nunique"),
                        unique_ccd=("resolved_ccd_id", "nunique"))
                   .reset_index().sort_values("entries", ascending=False))
    art_summary.to_csv(qc / "artifact_proxy_provenance_summary.tsv", sep="\t", index=False, lineterminator="\n")
    (art.groupby(["resolved_ccd_id", "ccd_name"], dropna=False)
        .agg(entries=("pair_id", "size"), unique_pdb=("pdb_id", "nunique"),
             unique_sources=("source_ligand_instance_id", "nunique"))
        .reset_index().sort_values(["entries", "resolved_ccd_id"], ascending=[False, True])
        .to_csv(qc / "artifact_proxy_top_ccd.tsv", sep="\t", index=False, lineterminator="\n"))

    ions = formal[pd.to_numeric(formal.observed_heavy_atom_count, errors="coerce").eq(1)].copy()
    ions.sort_values("pair_id").to_csv(qc / "one_heavy_atom_28_pair_audit.tsv", sep="\t", index=False, lineterminator="\n")
    (ions.groupby(["resolved_ccd_id", "ccd_name", "formula", "formal_charge", "ccd_type"], dropna=False)
         .agg(entries=("pair_id", "size"), unique_pdb=("pdb_id", "nunique"),
              unique_sources=("source_ligand_instance_id", "nunique"))
         .reset_index().sort_values(["entries", "resolved_ccd_id"], ascending=[False, True])
         .to_csv(qc / "one_heavy_atom_identity_summary.tsv", sep="\t", index=False, lineterminator="\n"))

    missing = entry[entry.unresolved_pocket_atom_count.gt(0)].copy()
    f3_path = root / "filter_03_ground_truth_structure_quality_database/runs/20260826_processing3_176900_strict_posebusters_02/output/filter3_retained_pairs.tsv.gz"
    f3 = pd.read_csv(f3_path, sep="\t", usecols=["pair_id", "benchmark_filter3_terminal_status", "reason_codes", "warning_codes"])
    missing = missing.merge(f3, on="pair_id", how="left", validate="one_to_one")
    missing.sort_values("pair_id").to_csv(qc / "unresolved_pocket_317_pair_audit.tsv.gz", sep="\t", index=False,
                                          compression={"method": "gzip", "mtime": 0}, lineterminator="\n")
    miss_cols = ["pocket_missing_backbone_heavy_atom_count",
                 "direct_binding_missing_sidechain_heavy_atom_count",
                 "nonbinding_pocket_missing_sidechain_heavy_atom_count"]
    summary = []
    for col in miss_cols:
        summary.append({"missing_atom_category": col, "affected_entries": int(missing[col].gt(0).sum()),
                        "missing_atoms": int(missing[col].sum())})
    pd.DataFrame(summary).to_csv(qc / "unresolved_pocket_category_summary.tsv", sep="\t", index=False, lineterminator="\n")
    (missing.groupby(["benchmark_filter3_terminal_status", "warning_codes"], dropna=False)
            .size().rename("entries").reset_index().sort_values("entries", ascending=False)
            .to_csv(qc / "unresolved_pocket_status_summary.tsv", sep="\t", index=False, lineterminator="\n"))

    cath_path = root / "protein_provenance_annotation/references/sifts_20260826/flatfiles/pdb_chain_cath_uniprot.tsv.gz"
    with pd.io.common.get_handle(cath_path, "r", compression="gzip") as handle:
        header = handle.handle.readline().strip()
        columns = handle.handle.readline().strip()
    cath_report = [
        "metric\tvalue\tinterpretation",
        f"reported_unique_CATH_ID\t58346\tDistinct SIFTS CATH_ID strings among participating receptor chains",
        f"SIFTS_file_comment\t{header}\tFrozen reference version line",
        f"SIFTS_columns\t{columns}\tOnly a CATH domain identifier is supplied; no C/A/T/H classification code",
        "CROWN_comparable_classification_count\tNA\tCannot derive classification IDs from the frozen SIFTS bridge alone",
    ]
    (qc / "cath_metric_review.tsv").write_text("\n".join(cath_report) + "\n", encoding="utf-8")

    checks = [
        ("artifact_proxy_entries", len(art), 3526),
        ("artifact_proxy_all_dual_source_pass", int(art.dual_source_pass.fillna(False).sum()), len(art)),
        ("one_heavy_atom_entries", len(ions), 28),
        ("unresolved_pocket_entries", len(missing), 317),
        ("unresolved_direct_binding_entries", int(missing.direct_binding_missing_sidechain_heavy_atom_count.gt(0).sum()), 0),
        ("unresolved_backbone_entries", int(missing.pocket_missing_backbone_heavy_atom_count.gt(0).sum()), 0),
    ]
    out = pd.DataFrame(checks, columns=["check", "observed", "expected"])
    out["pass"] = out.observed == out.expected
    out.to_csv(qc / "four_metric_review_validation.tsv", sep="\t", index=False, lineterminator="\n")
    if not out["pass"].all():
        raise RuntimeError("Focused metric review validation failed")


if __name__ == "__main__":
    main()
