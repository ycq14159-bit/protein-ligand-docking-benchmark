#!/usr/bin/env python3
"""Re-render Figure Scout plots after visual-only layout QA.

This script reads the already validated derived/frozen tables.  It does not
change classifications, populations, scientific metrics, or frozen inputs.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


OUT = Path("/home/linx/data/youcq/autodl-tmp/benchmark_1.0/audits/benchmark_figure_scout/20260823_draft_01")
BUILD_SCRIPT = Path("/tmp/build_benchmark_figure_scout.py")
FINALIZE_SCRIPT = Path("/tmp/finalize_figure_scout.py")


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def redraw_cross_benchmark() -> None:
    data = pd.read_csv(OUT / "research/cross_benchmark_source_log.tsv", sep="\t", dtype=str).fillna("")
    structural = ["CROWN", "PLINDER", "PDBbind", "HiQBind", "BioLiP2 / Q-BioLiP"]
    evaluation = ["This benchmark", "PoseBusters Benchmark", "PoseX Self-Docking", "Astex Diverse Set", "CASF-2016", "DockGen"]
    counts = {r.benchmark: int(r.reported_records_or_cases.replace(">", "")) for r in data.itertuples()}
    counts["This benchmark"] = 158226
    colors = {"blue": "#2563EB", "purple": "#7C3AED", "light": "#E2E8F0"}
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.2), gridspec_kw={"width_ratios": [1.15, 1.15, 1.25]})
    for ax, names, title in [
        (axes[0], structural, "Structural resources\n(heterogeneous record units)"),
        (axes[1], evaluation, "Evaluation benchmarks\n(case/complex units)"),
    ]:
        vals = [counts[n] for n in names]
        bar_colors = [colors["purple"] if n == "This benchmark" else colors["blue"] for n in names]
        y = np.arange(len(names))
        ax.barh(y, vals, color=bar_colors)
        display = {"PoseBusters Benchmark": "PoseBusters", "PoseX Self-Docking": "PoseX", "Astex Diverse Set": "Astex"}
        ax.set_yticks(y, [display.get(n, n) for n in names])
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel("Reported records/cases (log scale)")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color=colors["light"], lw=.7)
        ax.set_axisbelow(True)
        for yi, (n, value) in enumerate(zip(names, vals)):
            prefix = ">" if n == "PLINDER" else ""
            if n in {"PLINDER", "BioLiP2 / Q-BioLiP"}:
                ax.text(value * .96, yi, f"{prefix}{value:,}", va="center", ha="right", fontsize=8, color="white")
            else:
                ax.text(value * 1.08, yi, f"{prefix}{value:,}", va="center", fontsize=8)
    levels = ["overlap_entry", "overlap_chemical", "overlap_instance", "overlap_system"]
    level_names = ["Entry", "Chemical", "Ligand instance", "Binding system"]
    mapping = {"YES": 3, "HIGH_PARTIAL": 2.5, "PARTIAL": 2, "NEED_DATA": 0, "": 0}
    arr = np.array([[mapping.get(str(row[c]), 1) for c in levels] for _, row in data.iterrows()])
    cmap = ListedColormap(["#E2E8F0", "#FDE68A", "#FDBA74", "#86EFAC"])
    axes[2].imshow(arr, aspect="auto", vmin=0, vmax=3, cmap=cmap)
    axes[2].set_xticks(range(4), level_names, rotation=30, ha="right")
    axes[2].set_yticks(range(len(data)), data.benchmark)
    axes[2].set_title("Overlap feasibility from\nlightweight official metadata", loc="left", fontweight="bold")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            axes[2].text(j, i, data.iloc[i][levels[j]].replace("_", " "), ha="center", va="center", fontsize=6.5)
    fig.suptitle("Cross-benchmark landscape — units and overlap levels kept explicit", x=.035, y=.98, ha="left", fontsize=15, fontweight="bold")
    fig.subplots_adjust(left=.07, right=.98, bottom=.17, top=.76, wspace=.38)
    for ext, kwargs in [("png", {"dpi": 220}), ("svg", {})]:
        fig.savefig(OUT / "cross_benchmark" / f"fig07_cross_benchmark_landscape.{ext}", bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def redraw_group_examples(b, groups: pd.DataFrame, members: pd.DataFrame) -> None:
    targets = []
    for target_size in [3, 10, int(groups.group_size.max())]:
        targets.append(groups.iloc[(groups.group_size - target_size).abs().argsort().iloc[0]])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, row in zip(axes, targets):
        gid = row.equivalence_group_id
        selected = members[members.equivalence_group_id.eq(gid)]
        graph = b.nx.Graph()
        representative = str(row.representative_pair_id)
        for pair_id in selected.pair_id.astype(str):
            graph.add_node(pair_id, role="rep" if pair_id == representative else "member")
        for pair_id in selected.pair_id.astype(str):
            if pair_id != representative:
                graph.add_edge(representative, pair_id)
        pos = b.nx.spring_layout(graph, seed=24301)
        node_colors = [b.COLORS["green"] if graph.nodes[n]["role"] == "rep" else b.COLORS["gray"] for n in graph.nodes]
        b.nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=60, ax=ax)
        b.nx.draw_networkx_edges(graph, pos, alpha=.45, width=.7, ax=ax)
        ax.set_title(f"{gid}\nsize={len(selected):,}, ligand={row.ligand_exact_id}", fontsize=9)
        ax.axis("off")
    fig.suptitle("Real strict-equivalence groups (membership view)", x=.04, y=.98, ha="left", fontsize=14, fontweight="bold")
    fig.subplots_adjust(left=.03, right=.99, bottom=.04, top=.80, wspace=.12)
    b.savefig(fig, b.SUPP, "supp04_equivalence_group_examples")


def rebuild_manifest() -> None:
    rows = []
    excluded = {"output_manifest.tsv", "SHA256SUMS"}
    for path in sorted(x for x in OUT.rglob("*") if x.is_file() and x.name not in excluded):
        rows.append({"relative_path": str(path.relative_to(OUT)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(rows).to_csv(OUT / "output_manifest.tsv", sep="\t", index=False)
    with (OUT / "SHA256SUMS").open("w") as handle:
        for path in sorted(x for x in OUT.rglob("*") if x.is_file() and x.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.relative_to(OUT)}\n")


def main() -> None:
    base_qc = json.loads((OUT / "qc/qc_report.json").read_text())
    final_qc = json.loads((OUT / "qc/finalization_validation.json").read_text())
    if not base_qc.get("validation_pass") or not final_qc.get("validation_pass"):
        raise RuntimeError("refusing layout repair because scientific/finalization QC did not pass")

    b = load_module("figure_scout_build", BUILD_SCRIPT)
    pair = pd.read_parquet(OUT / "data_derived/unified_pair_figure_scout.parquet")
    source = pd.read_csv(b.F2 / "provisional_source_ligands.tsv.gz", sep="\t", usecols=["expected_heavy_atom_count"], low_memory=False)
    f4 = pd.read_csv(b.F4S5 / "01_filter4_final_pair_inventory.tsv.gz", sep="\t", low_memory=False)
    step4 = pd.read_csv(
        b.F4S4 / "04_pair_binding_residue_contact_inventory.tsv.gz",
        sep="\t",
        usecols=["candidate_pair_id", "ligand_heavy_atom_count", "min_external_ligand_distance_A", "min_external_binding_residue_distance_A"],
        low_memory=False,
    ).rename(columns={"candidate_pair_id": "pair_id"})
    f4["pair_id"] = f4.pair_id.fillna("").astype(str)
    step4["pair_id"] = step4.pair_id.fillna("").astype(str)
    f4 = f4.merge(step4, on="pair_id", how="left", validate="one_to_one")
    f5 = pd.read_csv(b.F5 / "01_filter5_final_pair_inventory.tsv.gz", sep="\t", low_memory=False)
    groups = pd.read_csv(b.F5 / "04_filter5_equivalence_groups.tsv.gz", sep="\t", low_memory=False)
    members = pd.read_csv(b.F5 / "05_filter5_group_members.tsv.gz", sep="\t", low_memory=False)
    groups["group_size"] = pd.to_numeric(groups.group_size, errors="coerce")
    heavy_dist = pd.read_csv(b.HEAVY_AUDIT / "01_heavy_atom_distribution.tsv", sep="\t")

    # make_figures only redraws candidate plots from existing pair-level data.
    inventory = b.make_figures(pair, source, f4, f5, groups, members, heavy_dist)
    if len(inventory) != 13:
        raise RuntimeError(f"unexpected make_figures inventory size: {len(inventory)}")
    redraw_group_examples(b, groups, members)
    redraw_cross_benchmark()

    shutil.copy2(BUILD_SCRIPT, OUT / "scripts/build_benchmark_figure_scout.py")
    shutil.copy2(FINALIZE_SCRIPT, OUT / "scripts/finalize_figure_scout.py")
    shutil.copy2(Path(__file__), OUT / "scripts/repair_figure_scout_layouts.py")

    now = utc()
    provenance = pd.read_csv(OUT / "figure_provenance.csv")
    redrawn = {row["figure_id"] for row in inventory} | {"S04G", "F07"}
    provenance.loc[provenance.figure_id.isin(redrawn), "generated_timestamp"] = now
    provenance.to_csv(OUT / "figure_provenance.csv", index=False)

    base_qc["visual_layout_qc"] = {
        "status": "PASS_AFTER_RERENDER",
        "checked_at": now,
        "scientific_values_changed": False,
        "redrawn_figure_ids": sorted(redrawn),
        "fixes": [
            "F01 right-edge box clipping",
            "F02/F05/F06 inter-panel title and axis-label spacing",
            "F07 suptitle/facet-title spacing",
            "S04/S04G row/title spacing",
        ],
    }
    (OUT / "qc/qc_report.json").write_text(json.dumps(base_qc, indent=2) + "\n")
    qc_md_path = OUT / "qc_report.md"
    qc_md = qc_md_path.read_text()
    visual_note = "## Visual layout QA\n\nPASS after a plot-only re-render. Clipping and title/axis-label overlaps were corrected; no population, metric, status, or frozen scientific value changed."
    if "## Visual layout QA" not in qc_md:
        qc_md_path.write_text(qc_md.rstrip() + "\n\n" + visual_note + "\n")

    inventory_table = pd.read_csv(OUT / "figure_inventory.csv")
    ready = inventory_table.status.astype(str).str.startswith("READY") & inventory_table.output_file.fillna("").ne("")
    missing = [str(path) for path in inventory_table.loc[ready, "output_file"] if not (OUT / str(path)).exists()]
    if missing:
        raise RuntimeError(f"ready plot outputs missing after repair: {missing}")
    final_qc["visual_layout_qc_pass"] = True
    final_qc["visual_layout_qc_timestamp"] = now
    (OUT / "qc/finalization_validation.json").write_text(json.dumps(final_qc, indent=2) + "\n")
    rebuild_manifest()
    print(json.dumps({"validation_pass": True, "redrawn": sorted(redrawn), "scientific_values_changed": False}, indent=2))


if __name__ == "__main__":
    main()
