#!/usr/bin/env python3
"""Generate the paper-primary harmonized six-dataset CROWN-style figure."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
from scipy.stats import gaussian_kde

FILES = {
    "PDBbind": "pdbbind_properties_quality.parquet",
    "HiQBind": "hiqbind_properties_quality.parquet",
    "BioLiP2": "biolip2_properties_quality.parquet",
    "PLINDER": "plinder_quality.parquet",
    "CROWN": "crown_quality.parquet",
    "Ours": "ours_properties_harmonized_quality.parquet",
}
COLORS = {
    "PDBbind": "#4C78A8", "HiQBind": "#F58518", "BioLiP2": "#54A24B",
    "PLINDER": "#E45756", "CROWN": "#B279A2", "Ours": "#222222",
}
PANELS = [
    ("heavy_atoms", "(a) Heavy atoms", (0, 80)),
    ("rotatable_bonds", "(b) Rotatable bonds", (0, 25)),
    ("hba", "(c) H-bond acceptors", (0, 20)),
    ("qed", "(d) QED", (0, 1)),
    ("resolution", "(e) Resolution (Å)", (0.5, 4.0)),
    ("rsr", "(f) Ligand RSR", (0, 0.5)),
    ("rscc", "(g) Ligand RSCC", (0.5, 1.0)),
]


def density(values, limits, points=400):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    grid = np.linspace(*limits, points)
    if len(values) < 2 or np.all(values == values[0]):
        return grid, np.zeros_like(grid)
    estimator = gaussian_kde(values, bw_method="scott")
    estimator.set_bandwidth(estimator.factor * 0.85)
    return grid, estimator(grid)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    frames = {dataset: pd.read_parquet(args.prepared_root / filename) for dataset, filename in FILES.items()}
    counts = {dataset: len(frame) for dataset, frame in frames.items()}

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
        "axes.labelsize": 8.5, "legend.fontsize": 8.2, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(13.2, 6.9), constrained_layout=True)
    outer = fig.add_gridspec(2, 4, height_ratios=[1, 1.13])
    axes = []
    for index in range(4):
        axes.append((fig.add_subplot(outer[0, index]), None))
    for index in range(3):
        nested = outer[1, index].subgridspec(2, 1, height_ratios=[4.2, 1.1], hspace=0.07)
        axes.append((fig.add_subplot(nested[0]), fig.add_subplot(nested[1])))
    legend_axis = fig.add_subplot(outer[1, 3])
    display_rows = []
    legend_handles = []
    legend_labels = []
    for panel_index, ((metric, title, limits), (axis, missing_axis)) in enumerate(zip(PANELS, axes)):
        for dataset, frame in frames.items():
            values = pd.to_numeric(frame[metric], errors="coerce").dropna().to_numpy()
            grid, curve = density(values, limits)
            width = 2.2 if dataset == "Ours" else 1.45
            line, = axis.plot(grid, curve, color=COLORS[dataset], lw=width, alpha=0.98)
            if panel_index == 0:
                legend_handles.append(line)
                legend_labels.append(f"{dataset} (n = {counts[dataset]:,})")
            within = int(((values >= limits[0]) & (values <= limits[1])).sum())
            display_rows.append({
                "dataset": dataset, "metric": metric, "available_N": len(values),
                "within_display_range_N": within, "outside_display_range_N": len(values) - within,
                "x_min": limits[0], "x_max": limits[1],
            })
        axis.set_title(title, loc="left", fontweight="semibold")
        axis.set_xlim(limits)
        axis.set_ylabel("Density")
        axis.grid(axis="y", color="#D9D9D9", lw=0.55, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        if missing_axis is not None:
            missing = [frames[name][metric].isna().mean() * 100 for name in FILES]
            x = np.arange(len(FILES))
            missing_axis.bar(x, missing, color=[COLORS[name] for name in FILES], width=0.72)
            missing_axis.set_ylim(0, 100)
            missing_axis.set_ylabel("No data\n(%)", fontsize=7.3)
            missing_axis.set_xticks(x, ["PDB", "HiQ", "Bio", "PLI", "CROWN", "Ours"], fontsize=6.8)
            missing_axis.set_yticks([0, 50, 100])
            missing_axis.tick_params(axis="y", labelsize=6.5)
            missing_axis.spines[["top", "right"]].set_visible(False)
            missing_axis.grid(axis="y", color="#E5E5E5", lw=0.45)
        else:
            axis.set_xlabel(title.split(") ", 1)[-1])

    legend_axis.axis("off")
    legend_axis.legend(legend_handles, legend_labels, loc="center left", frameon=False, handlelength=3)
    legend_axis.text(
        0, 0.19, "Gaussian KDE\nbw_adjust = 0.85\nEntry-weighted",
        transform=legend_axis.transAxes, color="#555555", linespacing=1.5,
    )
    stem = args.output_root / "figure_harmonized_comparison"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(display_rows).to_csv(args.output_root.parent / "qc" / "display_range_qc.tsv", sep="\t", index=False)
    try:
        commit = subprocess.check_output(["git", "-C", str(args.repo_root), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        commit = "UNAVAILABLE"
    provenance = {
        "figure": "figure_harmonized_comparison", "mode": "B_harmonized_fixed_version",
        "datasets": counts, "rdkit_version": rdkit.__version__,
        "descriptor_definitions": ["GetNumHeavyAtoms", "CalcNumRotatableBonds", "CalcNumHBA", "QED.qed"],
        "kde": {"implementation": "scipy.stats.gaussian_kde", "base_bandwidth": "Scott", "bw_adjust": 0.85},
        "axis_limits": {metric: limits for metric, _, limits in PANELS},
        "validation_sources": {
            "PDBbind": "wwPDB XML after unique deposited-coordinate mapping",
            "HiQBind": "wwPDB XML exact PDB+CCD+chain+residue mapping",
            "BioLiP2": "wwPDB XML unique PDB+CCD+chain mapping",
            "PLINDER": "released exact-instance wwPDB validation arrays",
            "CROWN": "released wwPDB-derived ligand validation metrics",
            "Ours": "frozen exact-instance wwPDB validation provenance",
        },
        "git_commit_before_figure_commit": commit,
    }
    (args.output_root / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
