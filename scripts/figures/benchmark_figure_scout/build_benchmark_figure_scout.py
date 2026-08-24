#!/usr/bin/env python3
"""Read-only full benchmark Figure Scout over frozen benchmark_1.0 results."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import shutil
import textwrap

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch
import networkx as nx

from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


ROOT = Path("/home/linx/data/youcq/autodl-tmp/benchmark_1.0")
RUN_ID = "20260823_draft_01"
OUT = ROOT / "audits/benchmark_figure_scout" / RUN_ID
DERIVED = OUT / "data_derived"
MAIN = OUT / "main_candidates"
SUPP = OUT / "supplementary_candidates"
CROSS = OUT / "cross_benchmark"
RESEARCH = OUT / "research"
SCRIPTS = OUT / "scripts"
LOGS = OUT / "logs"
QC = OUT / "qc"
FUTURE = OUT / "future_docking"

F1 = ROOT / "filter_1_protein_receptor_qualification/release"
F2 = ROOT / "filter_2_ligand_qualification_v3/runs/20260804_full_01/output"
P2 = ROOT / "processing_2_assembly_ready_structure_preparation/runs/20260810_full_01/output"
P3 = ROOT / "processing_03_direct_contact_qualification/runs/20260811_full_01/output"
F3 = ROOT / "filter_03_ground_truth_structure_quality_v2/runs/20260814_full_01/output"
F4ROOT = ROOT / "filter_04_crystal_packing_influence"
F4S1 = F4ROOT / "step_01_lattice_neighbor_search/runs/step01_full_v3/output"
F4S2 = F4ROOT / "step_02_biological_assembly_equivalence/runs/step02_full_v3/output"
F4S3 = F4ROOT / "step_03_direct_ligand_crystal_contact/runs/step03_full_v1/output"
F4S4 = F4ROOT / "step_04_binding_residue_mediated_crystal_contact/runs/step04_full_v1/output"
F4S5 = F4ROOT / "step_05_final_crystal_packing_decision/runs/step05_full_v1/output"
F5 = ROOT / "filter_05_equivalent_redocking_case/step_05_strict_equivalent_grouping_and_representative_selection/runs/step05_full_v1/output"
HEAVY_AUDIT = ROOT / "audits/ligand_minimum_pose_complexity_census/20260823_full_02"

COLORS = {
    "blue": "#2563EB", "orange": "#F59E0B", "green": "#16A34A", "red": "#DC2626",
    "purple": "#7C3AED", "gray": "#64748B", "light": "#E2E8F0", "dark": "#0F172A",
}
POP_ORDER = ["F3 HIGH+GOOD", "F4 PASS", "F5 FINAL"]
QUALITY_ORDER = ["HIGH", "GOOD"]
HEAVY_BINS = ["<3", "3–5", "6–10", "11–20", "21–40", ">40", "missing"]


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def norm(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str)


def as_num(frame: pd.DataFrame, columns: list[str]) -> None:
    for col in columns:
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def heavy_bin(values: pd.Series) -> pd.Categorical:
    x = pd.to_numeric(values, errors="coerce")
    labels = np.select(
        [x.lt(3), x.between(3, 5), x.between(6, 10), x.between(11, 20), x.between(21, 40), x.gt(40)],
        HEAVY_BINS[:-1], default="missing"
    )
    return pd.Categorical(labels, categories=HEAVY_BINS, ordered=True)


def ecdf(values: pd.Series | np.ndarray):
    x = np.sort(pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy())
    if not len(x):
        return x, x
    return x, np.arange(1, len(x) + 1) / len(x)


def savefig(fig: plt.Figure, folder: Path, stem: str) -> tuple[str, str]:
    png = folder / f"{stem}.png"
    svg = folder / f"{stem}.svg"
    fig.savefig(png, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return str(png.relative_to(OUT)), str(svg.relative_to(OUT))


def style_ax(ax, grid=True):
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color=COLORS["light"], linewidth=.7, alpha=.8)
        ax.set_axisbelow(True)


def annotate_bars(ax, fmt="{:,}", rotation=0):
    for patch in ax.patches:
        value = patch.get_height()
        if np.isfinite(value):
            ax.annotate(fmt.format(value), (patch.get_x() + patch.get_width()/2, value),
                        xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                        fontsize=8, rotation=rotation)


def qsummary(frame: pd.DataFrame, population: str, metrics: list[str]) -> list[dict]:
    rows = []
    for metric in metrics:
        vals = pd.to_numeric(frame[metric], errors="coerce") if metric in frame else pd.Series(dtype=float)
        a = vals.dropna()
        rows.append({
            "population": population, "metric": metric, "n_total": len(frame), "n_available": len(a),
            "n_missing": len(frame)-len(a), "available_fraction": len(a)/len(frame) if len(frame) else np.nan,
            "median": a.median() if len(a) else np.nan, "q1": a.quantile(.25) if len(a) else np.nan,
            "q3": a.quantile(.75) if len(a) else np.nan, "min": a.min() if len(a) else np.nan,
            "max": a.max() if len(a) else np.nan,
        })
    return rows


def aggregate_pocket(path: Path) -> pd.DataFrame:
    parts = []
    for file in sorted(path.rglob("*.parquet")):
        d = pq.read_table(file, columns=["pair_id", "chain_instance_id", "protein_residue_id"]).to_pandas()
        g = d.groupby("pair_id", sort=False).agg(
            pocket_residue_count=("protein_residue_id", "size"),
            pocket_chain_count=("chain_instance_id", "nunique"),
        ).reset_index()
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    if out.pair_id.duplicated().any():
        out = out.groupby("pair_id", as_index=False).agg(
            pocket_residue_count=("pocket_residue_count", "sum"),
            pocket_chain_count=("pocket_chain_count", "max"),
        )
    return out


def aggregate_contacts(path: Path) -> pd.DataFrame:
    parts = []
    for file in sorted(path.rglob("*.parquet")):
        d = pq.read_table(file, columns=[
            "ligand_assembly_placement_id", "protein_residue_id", "ligand_atom_id", "distance_angstrom"
        ]).to_pandas()
        g = d.groupby("ligand_assembly_placement_id", sort=False).agg(
            qualifying_contact_count=("distance_angstrom", "size"),
            native_min_contact_distance=("distance_angstrom", "min"),
            contacted_ligand_atom_count=("ligand_atom_id", "nunique"),
            contacted_protein_residue_count=("protein_residue_id", "nunique"),
        ).reset_index()
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    if out.ligand_assembly_placement_id.duplicated().any():
        out = out.groupby("ligand_assembly_placement_id", as_index=False).agg(
            qualifying_contact_count=("qualifying_contact_count", "sum"),
            native_min_contact_distance=("native_min_contact_distance", "min"),
            contacted_ligand_atom_count=("contacted_ligand_atom_count", "max"),
            contacted_protein_residue_count=("contacted_protein_residue_count", "max"),
        )
    return out


def aggregate_binding_quality(path: Path) -> pd.DataFrame:
    parts = []
    for file in sorted(path.rglob("*.parquet")):
        d = pq.read_table(file, columns=["pair_id", "rsrz", "rscc", "rsr", "qualifying_atomic_contact_count"]).to_pandas()
        g = d.groupby("pair_id", sort=False).agg(
            direct_binding_residue_rows=("pair_id", "size"),
            direct_binding_rsrz_available=("rsrz", "count"),
            direct_binding_rsrz_median=("rsrz", "median"),
            direct_binding_rsrz_max=("rsrz", "max"),
            direct_binding_rscc_min=("rscc", "min"),
            direct_binding_rsr_max=("rsr", "max"),
            direct_binding_contact_count=("qualifying_atomic_contact_count", "sum"),
        ).reset_index()
        parts.append(g)
    out = pd.concat(parts, ignore_index=True)
    if out.pair_id.duplicated().any():
        raise RuntimeError("binding-residue pair rows unexpectedly span parquet partitions")
    return out


def build_rdkit_descriptors(topology: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    t = topology[["component_id", "canonical_smiles_from_ccd_graph"]].copy()
    t["component_id"] = norm(t.component_id)
    t["canonical_smiles_from_ccd_graph"] = norm(t.canonical_smiles_from_ccd_graph)
    conflicts = t.groupby("component_id").canonical_smiles_from_ccd_graph.nunique().sort_values(ascending=False)
    conflict_components = conflicts[conflicts.gt(1)]
    unique = t.sort_values(["component_id", "canonical_smiles_from_ccd_graph"]).drop_duplicates("component_id")
    rows = []
    for row in unique.itertuples(index=False):
        smi = row.canonical_smiles_from_ccd_graph
        result = {
            "component_id": row.component_id, "canonical_smiles": smi, "descriptor_status": "MISSING_SMILES",
            "rdkit_heavy_atom_count": np.nan, "rotatable_bond_count": np.nan, "ring_count": np.nan,
            "aromatic_ring_count": np.nan, "hbd_count": np.nan, "hba_count": np.nan,
            "fraction_csp3": np.nan, "murcko_scaffold_smiles": "",
        }
        if smi:
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    result.update({
                        "descriptor_status": "PASS",
                        "rdkit_heavy_atom_count": mol.GetNumHeavyAtoms(),
                        "rotatable_bond_count": Lipinski.NumRotatableBonds(mol),
                        "ring_count": rdMolDescriptors.CalcNumRings(mol),
                        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(mol),
                        "hbd_count": rdMolDescriptors.CalcNumHBD(mol),
                        "hba_count": rdMolDescriptors.CalcNumHBA(mol),
                        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
                        "murcko_scaffold_smiles": MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
                    })
                else:
                    result["descriptor_status"] = "RDKIT_PARSE_FAILED"
            except Exception as exc:
                result["descriptor_status"] = f"RDKIT_ERROR:{type(exc).__name__}"
        rows.append(result)
    return pd.DataFrame(rows), {
        "unique_components": len(unique), "components_with_multiple_frozen_graph_smiles": int(len(conflict_components)),
        "conflict_component_examples": conflict_components.head(20).index.tolist(),
    }


def make_figures(pair: pd.DataFrame, source: pd.DataFrame, f4: pd.DataFrame,
                 f5: pd.DataFrame, groups: pd.DataFrame, members: pd.DataFrame,
                 heavy_dist: pd.DataFrame) -> list[dict]:
    inventory = []
    hg = pair[pair.f3_hg].copy()
    f4pass = pair[pair.f4_pass].copy()
    final = pair[pair.f5_retained].copy()

    # Figure 1: separate lanes and pair attrition.
    fig = plt.figure(figsize=(13, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], wspace=.26, hspace=.34)
    ax = fig.add_subplot(gs[0, :])
    ax.axis("off")
    lanes = [
        ("Receptor lane", [("F1 entries", 248037), ("Assemblies", 360611), ("Source polypeptides", 1073451), ("Assembly receptor chains", 2145537), ("P2 coord-ready", 834222)]),
        ("Ligand lane", [("F2 source", 852968), ("Mapped source", 851966), ("Placements", 1151324), ("Active assembly keys", 234975), ("P2 complete source", 746509)]),
        ("Pair lane", [("P3 candidate pairs", 744580), ("F3 HIGH+GOOD", 336412), ("F4 PASS", 241545), ("F5 retained", 158226)]),
    ]
    lane_colors = [COLORS["blue"], COLORS["orange"], COLORS["purple"]]
    for li, (label, nodes) in enumerate(lanes):
        y = 2.3 - li
        ax.text(.01, y, label, va="center", ha="left", fontweight="bold", color=lane_colors[li], transform=ax.transData)
        # Keep the rounded boxes fully inside the axes; the previous .93 endpoint
        # clipped the rightmost box after export.
        xs = np.linspace(.20, .90, len(nodes))
        for i, ((name, value), x) in enumerate(zip(nodes, xs)):
            box = FancyBboxPatch((x-.075, y-.25), .15, .5, boxstyle="round,pad=0.02,rounding_size=0.02",
                                 edgecolor=lane_colors[li], facecolor="white", linewidth=1.5)
            ax.add_patch(box)
            ax.text(x, y+.06, name, ha="center", va="center", fontsize=8.5)
            ax.text(x, y-.09, f"{value:,}", ha="center", va="center", fontsize=10, fontweight="bold")
            if i:
                ax.annotate("", xy=(x-.078, y), xytext=(xs[i-1]+.078, y),
                            arrowprops=dict(arrowstyle="->", color=lane_colors[li], lw=1.5))
    ax.set_xlim(0, 1); ax.set_ylim(-.1, 2.8)
    ax.set_title("Benchmark construction landscape — statistical units remain separate", loc="left", fontsize=15, fontweight="bold")
    ax = fig.add_subplot(gs[1, 0])
    stages = ["F3 HIGH+GOOD", "F4 PASS", "F5 retained"]
    counts = [336412, 241545, 158226]
    ax.bar(stages, counts, color=[COLORS["blue"], COLORS["green"], COLORS["purple"]])
    annotate_bars(ax); ax.set_ylabel("Pair count"); ax.set_title("Pair-level attrition", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=12); style_ax(ax)
    ax = fig.add_subplot(gs[1, 1])
    state_order = ["F5_RETAIN_UNIQUE", "F5_RETAIN_REPRESENTATIVE", "F5_REVIEW_RETAIN", "F5_REDUNDANT_EQUIVALENT_CASE"]
    vc = f5.filter5_final_status.value_counts().reindex(state_order).fillna(0)
    labels = ["Retain unique", "Representative", "Review retain", "Redundant"]
    ax.bar(labels, vc.values, color=[COLORS["blue"], COLORS["green"], COLORS["orange"], COLORS["gray"]])
    annotate_bars(ax); ax.set_ylabel("Pair count"); ax.set_title("Filter 5 final-state composition", loc="left", fontweight="bold")
    ax.tick_params(axis="x", rotation=15); style_ax(ax)
    p = savefig(fig, MAIN, "fig01_construction_landscape")
    inventory.append(inv("F01", "Benchmark construction landscape", "How was the benchmark constructed without mixing units?",
                         "Frozen F1/F2/P2/P3/F3/F4/F5 ledgers", "counts by receptor/ligand/pair lane", "READY", p[0], "KEEP"))

    # Figure 2: selection pressure and heavy atom conditional rates.
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    populations = [("F2 source", source.expected_heavy_atom_count, COLORS["gray"]),
                   ("F3 HIGH+GOOD", hg.expected_heavy_atom_count, COLORS["blue"]),
                   ("F4 PASS", f4pass.expected_heavy_atom_count, COLORS["green"]),
                   ("F5 FINAL", final.expected_heavy_atom_count, COLORS["purple"])]
    for label, vals, color in populations:
        x, y = ecdf(vals); axes[0,0].plot(x, y, label=f"{label} (N={len(vals):,})", color=color, lw=1.7)
    axes[0,0].set_xlim(0, 80); axes[0,0].set_xlabel("CCD expected heavy-atom count"); axes[0,0].set_ylabel("ECDF")
    axes[0,0].set_title("Ligand-complexity shift", loc="left", fontweight="bold"); axes[0,0].legend(fontsize=8); style_ax(axes[0,0])
    heavy_rate = hg.groupby("heavy_atom_bin", observed=False).agg(
        denominator=("pair_id", "size"), rejected=("f4_reject", "sum")
    ).reset_index(); heavy_rate["rate"] = heavy_rate.rejected / heavy_rate.denominator
    axes[0,1].bar(heavy_rate.heavy_atom_bin.astype(str), heavy_rate.rate*100, color=COLORS["red"])
    axes[0,1].set_ylabel("F4 reject rate (%)"); axes[0,1].set_xlabel("Heavy-atom bin")
    axes[0,1].set_title("P(F4 REJECT | ligand size)", loc="left", fontweight="bold"); axes[0,1].tick_params(axis="x", rotation=20); style_ax(axes[0,1])
    f5_rate = f4pass.groupby("heavy_atom_bin", observed=False).agg(
        denominator=("pair_id", "size"), redundant=("f5_redundant", "sum")
    ).reset_index(); f5_rate["rate"] = f5_rate.redundant / f5_rate.denominator
    axes[1,0].bar(f5_rate.heavy_atom_bin.astype(str), f5_rate.rate*100, color=COLORS["gray"])
    axes[1,0].set_ylabel("F5 redundant rate (%)"); axes[1,0].set_xlabel("Heavy-atom bin")
    axes[1,0].set_title("P(F5 REDUNDANT | ligand size)", loc="left", fontweight="bold"); axes[1,0].tick_params(axis="x", rotation=20); style_ax(axes[1,0])
    qrows=[]
    for q in ["HIGH", "GOOD"]:
        sub=hg[hg.f3_quality.eq(q)]; p4rate=sub.f4_pass.mean(); p5rate=sub.f5_retained.mean()
        qrows.append((q,p4rate,p5rate))
    x=np.arange(2); w=.34
    axes[1,1].bar(x-w/2,[r[1]*100 for r in qrows],w,label="Survive F4",color=COLORS["green"])
    axes[1,1].bar(x+w/2,[r[2]*100 for r in qrows],w,label="Reach F5 retained",color=COLORS["purple"])
    axes[1,1].set_xticks(x,[r[0] for r in qrows]); axes[1,1].set_ylabel("Share of F3 class (%)"); axes[1,1].set_ylim(0,100)
    axes[1,1].set_title("HIGH vs GOOD downstream survival", loc="left", fontweight="bold"); axes[1,1].legend(); style_ax(axes[1,1])
    fig.suptitle("Selection pressure and distribution shift", x=.06, ha="left", fontsize=15, fontweight="bold")
    fig.subplots_adjust(left=.09, right=.98, bottom=.08, top=.91, hspace=.42, wspace=.20)
    p=savefig(fig, MAIN, "fig02_selection_pressure")
    inventory.append(inv("F02", "Selection pressure and distribution shift", "Which ligand sizes and F3 classes are preferentially removed?",
                         "F2 source; F3 HIGH+GOOD; F4 PASS; F5 retained", "heavy atoms, F4 rejection, F5 redundancy, F3 class", "READY", p[0], "KEEP"))
    heavy_rate.to_csv(DERIVED/"heavy_atom_bin_f4_reject_rate.tsv",sep="\t",index=False)
    f5_rate.to_csv(DERIVED/"heavy_atom_bin_f5_redundant_rate.tsv",sep="\t",index=False)

    # Figure 2b: shared descriptor ECDFs.
    metrics=[("formula_weight","Molecular weight (Da)",(0,1000)),("rotatable_bond_count","Rotatable bonds",(0,25)),
             ("ring_count","Ring count",(0,10)),("pocket_residue_count","6 Å pocket residues",(0,100)),
             ("receptor_total_declared_length","Total receptor sequence length",(0,2000)),("direct_binding_residue_rows","Direct-binding residues",(0,30))]
    fig,axes=plt.subplots(2,3,figsize=(14,8.2))
    for ax,(metric,label,xlim) in zip(axes.flat,metrics):
        for pop,df,color in [("F3 H+G",hg,COLORS["blue"]),("F4 PASS",f4pass,COLORS["green"]),("F5 FINAL",final,COLORS["purple"])]:
            x,y=ecdf(df[metric]); ax.plot(x,y,label=f"{pop} n={len(x):,}",color=color,lw=1.5)
        ax.set_xlim(*xlim); ax.set_xlabel(label); ax.set_ylabel("ECDF"); ax.set_title(label,loc="left",fontweight="bold"); style_ax(ax)
    axes[0,0].legend(fontsize=7)
    fig.suptitle("Shared-property shifts across pair-level populations",x=.06,ha="left",fontsize=15,fontweight="bold")
    p=savefig(fig,MAIN,"fig02b_descriptor_shifts")
    inventory.append(inv("F02B","Shared-property distribution shifts","Does filtering shift chemistry, pocket size or receptor size?",
                         "F3 HIGH+GOOD, F4 PASS, F5 retained","MW, rotatable bonds, rings, pocket/receptor/binding size","READY",p[0],"MERGE"))

    # Figure 3: crystal packing landscape.
    fig,axes=plt.subplots(2,3,figsize=(14.5,8.5))
    summary=f4.groupby(["filter4_decision","filter4_reason"],dropna=False).size().reset_index(name="count")
    reason_order=["NO_CRYSTALLOGRAPHIC_NEIGHBOR","NO_EXTERNAL_CRYSTAL_NEIGHBOR","EXTERNAL_NEIGHBOR_NO_RELEVANT_CONTACT",
                  "DIRECT_LIGAND_CRYSTAL_CONTACT","BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT","BA_EQUIVALENCE_UNRESOLVED"]
    colors=[COLORS["blue"],COLORS["green"],"#0D9488",COLORS["red"],COLORS["orange"],COLORS["purple"]]
    s=summary.set_index("filter4_reason").reindex(reason_order).fillna(0)
    axes[0,0].barh([x.replace("_"," ").title() for x in reason_order][::-1],s["count"].values[::-1],color=colors[::-1])
    axes[0,0].set_xlabel("Pairs"); axes[0,0].set_title("Final decision routes",loc="left",fontweight="bold"); style_ax(axes[0,0],grid=False)
    decision=f4.filter4_decision.value_counts().reindex(["PASS","REJECT","REVIEW"]).fillna(0)
    axes[0,1].bar(decision.index,decision.values,color=[COLORS["green"],COLORS["red"],COLORS["purple"]]); annotate_bars(axes[0,1])
    axes[0,1].set_ylabel("Pairs"); axes[0,1].set_title("PASS / REJECT / REVIEW",loc="left",fontweight="bold"); style_ax(axes[0,1])
    neigh=[799269,355846,2]
    axes[0,2].bar(["BA-equivalent","External","Review"],neigh,color=[COLORS["blue"],COLORS["orange"],COLORS["purple"]]); annotate_bars(axes[0,2])
    axes[0,2].set_ylabel("Neighbour instances"); axes[0,2].set_title("1,155,117 neighbour instances",loc="left",fontweight="bold"); axes[0,2].tick_params(axis="x",rotation=12); style_ax(axes[0,2])
    clipped=f4.n_external_instances.clip(upper=20)
    axes[1,0].hist(clipped,bins=np.arange(-.5,21.5,1),color=COLORS["gray"]); axes[1,0].set_yscale("log")
    axes[1,0].set_xlabel("External neighbours per pair (20 = ≥20)"); axes[1,0].set_ylabel("Pairs, log scale"); axes[1,0].set_title("External-neighbour burden",loc="left",fontweight="bold"); style_ax(axes[1,0])
    direct=f4[f4.filter4_reason.eq("DIRECT_LIGAND_CRYSTAL_CONTACT")].fraction_ligand_heavy_atoms_contacted_4A.dropna()
    bridge=f4[f4.filter4_reason.eq("BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT")].fraction_binding_residues_crystal_bridged.dropna()
    for vals,label,color in [(direct,"Direct ligand contact",COLORS["red"]),(bridge,"Pocket-mediated",COLORS["orange"])]:
        x,y=ecdf(vals); axes[1,1].plot(x,y,label=f"{label} (n={len(x):,})",color=color,lw=1.8)
    axes[1,1].set_xlabel("Fraction affected"); axes[1,1].set_ylabel("ECDF"); axes[1,1].set_title("Crystal-contact severity",loc="left",fontweight="bold"); axes[1,1].legend(fontsize=8); style_ax(axes[1,1])
    qtab=pair[pair.f3_hg].groupby(["f3_quality","filter4_decision"]).size().unstack(fill_value=0).reindex(QUALITY_ORDER)
    qrate=qtab.div(qtab.sum(axis=1),axis=0)*100
    bottom=np.zeros(len(qrate))
    for decision_name,color in [("PASS",COLORS["green"]),("REJECT",COLORS["red"]),("REVIEW",COLORS["purple"])]:
        vals=qrate.get(decision_name,pd.Series(0,index=qrate.index)); axes[1,2].bar(qrate.index,vals,bottom=bottom,label=decision_name,color=color); bottom+=vals.to_numpy()
    axes[1,2].set_ylim(0,100); axes[1,2].set_ylabel("Within-class percentage"); axes[1,2].set_title("Crystal decision by F3 quality",loc="left",fontweight="bold"); axes[1,2].legend(fontsize=8); style_ax(axes[1,2])
    fig.suptitle("Crystal-packing landscape",x=.055,ha="left",fontsize=15,fontweight="bold")
    p=savefig(fig,MAIN,"fig03_crystal_packing_landscape")
    inventory.append(inv("F03","Crystal-packing landscape","How often and how severely does the crystal lattice affect the binding site?",
                         "336,412 F4 input pairs; 1,155,117 neighbour instances","F4 route, neighbour class/count, affected fractions, F3 class","READY",p[0],"KEEP"))

    # Figure 3b: distances and severity.
    fig,axes=plt.subplots(1,3,figsize=(14,4.2))
    for vals,label,color in [(f4.min_external_ligand_distance_A,"Ligand",COLORS["red"]),(f4.min_external_binding_residue_distance_A,"Binding residues",COLORS["orange"])]:
        x,y=ecdf(vals); axes[0].plot(x,y,label=f"{label} (n={len(x):,})",color=color)
    axes[0].set_xlim(0,10); axes[0].set_xlabel("Minimum external distance (Å)"); axes[0].set_ylabel("ECDF"); axes[0].legend(); axes[0].set_title("Minimum distances",loc="left",fontweight="bold"); style_ax(axes[0])
    axes[1].hexbin(f4.ligand_heavy_atom_count,f4.fraction_ligand_heavy_atoms_contacted_4A,gridsize=45,mincnt=1,bins="log",cmap="viridis")
    axes[1].set_xlim(0,100); axes[1].set_xlabel("Ligand heavy atoms"); axes[1].set_ylabel("Fraction ligand atoms contacted"); axes[1].set_title("Direct-contact severity vs ligand size",loc="left",fontweight="bold"); style_ax(axes[1],grid=False)
    axes[2].hexbin(f4.binding_residue_count,f4.n_crystal_bridged_binding_residues,gridsize=40,mincnt=1,bins="log",cmap="magma")
    axes[2].set_xlim(0,40); axes[2].set_ylim(0,30); axes[2].set_xlabel("Binding residues"); axes[2].set_ylabel("Crystal-bridged residues"); axes[2].set_title("Pocket-mediated burden",loc="left",fontweight="bold"); style_ax(axes[2],grid=False)
    fig.suptitle("Crystal-contact severity diagnostics",x=.06,ha="left",fontsize=14,fontweight="bold")
    p=savefig(fig,SUPP,"supp03_crystal_severity")
    inventory.append(inv("S03","Crystal-contact severity diagnostics","Which rejected cases have the most severe packing?",
                         "F4 pairs with available contact metrics","minimum distance, contact fractions, bridged residues","READY",p[0],"SUPPLEMENTARY"))

    # Figure 4: equivalence and redundancy.
    fig,axes=plt.subplots(2,3,figsize=(14.5,8.6))
    axes[0,0].bar(["F4 PASS","F5 retained"],[241545,158226],color=[COLORS["green"],COLORS["purple"]]); annotate_bars(axes[0,0]); axes[0,0].set_ylabel("Pairs"); axes[0,0].set_title("Strict-equivalence attrition",loc="left",fontweight="bold"); style_ax(axes[0,0])
    sizes=pd.to_numeric(groups.group_size,errors="coerce").dropna().astype(int)
    axes[0,1].hist(np.minimum(sizes,50),bins=np.arange(1.5,51.5),color=COLORS["blue"]); axes[0,1].set_yscale("log")
    axes[0,1].set_xlabel("Equivalent-group size (50 = ≥50)"); axes[0,1].set_ylabel("Groups, log scale"); axes[0,1].set_title("Group-size distribution",loc="left",fontweight="bold"); style_ax(axes[0,1])
    ux=np.sort(sizes.unique()); ccdf=np.array([(sizes>=v).mean() for v in ux]); axes[0,2].loglog(ux,ccdf,color=COLORS["purple"],lw=1.7)
    axes[0,2].set_xlabel("Group size"); axes[0,2].set_ylabel("P(group size ≥ x)"); axes[0,2].set_title("Group-size CCDF",loc="left",fontweight="bold"); style_ax(axes[0,2])
    redund=(sizes-1).sort_values(ascending=False).to_numpy(); cum=np.cumsum(redund)/redund.sum(); share=np.arange(1,len(redund)+1)/len(redund)
    axes[1,0].plot(share*100,cum*100,color=COLORS["red"]); axes[1,0].plot([0,100],[0,100],color=COLORS["gray"],lw=.8)
    axes[1,0].set_xlabel("Largest groups included (%)"); axes[1,0].set_ylabel("Redundant cases explained (%)"); axes[1,0].set_title("Redundancy concentration",loc="left",fontweight="bold"); style_ax(axes[1,0])
    f5pass=pair[pair.f4_pass].copy(); rtab=f5pass.groupby("heavy_atom_bin",observed=False).agg(n=("pair_id","size"),redundant=("f5_redundant","sum")).reset_index(); rtab["rate"]=rtab.redundant/rtab.n
    axes[1,1].bar(rtab.heavy_atom_bin.astype(str),rtab.rate*100,color=COLORS["gray"]); axes[1,1].set_ylabel("Redundant (%)"); axes[1,1].set_xlabel("Heavy-atom bin"); axes[1,1].tick_params(axis="x",rotation=20); axes[1,1].set_title("Redundancy vs ligand complexity",loc="left",fontweight="bold"); style_ax(axes[1,1])
    recurrence_bins=pd.qcut(f5pass.pdb_recurrence.rank(method="first"),q=6,labels=["Q1","Q2","Q3","Q4","Q5","Q6"])
    pr=f5pass.groupby(recurrence_bins,observed=False).agg(n=("pair_id","size"),redundant=("f5_redundant","sum"),median_pdb_recurrence=("pdb_recurrence","median")).reset_index(); pr["rate"]=pr.redundant/pr.n
    axes[1,2].plot(np.arange(len(pr)),pr.rate*100,marker="o",color=COLORS["orange"]); axes[1,2].set_xticks(np.arange(len(pr)),[f"{q}\nmed {m:.0f}" for q,m in zip(pr.iloc[:,0],pr.median_pdb_recurrence)])
    axes[1,2].set_ylabel("Redundant (%)"); axes[1,2].set_xlabel("PDB-recurrence quantile"); axes[1,2].set_title("Redundancy vs PDB recurrence",loc="left",fontweight="bold"); style_ax(axes[1,2])
    fig.suptitle("Filter 5 strict equivalence and redundancy",x=.055,ha="left",fontsize=15,fontweight="bold")
    p=savefig(fig,MAIN,"fig04_filter5_redundancy")
    inventory.append(inv("F04","Filter 5 equivalence and redundancy","How is redundancy distributed and concentrated?",
                         "241,545 F4 PASS pairs; 32,188 multi-member groups","group size, CCDF, concentration, ligand/PDB recurrence","READY",p[0],"KEEP"))

    # Representative-vs-redundant quality bias.
    rep=pair[pair.f5_status.eq("F5_RETAIN_REPRESENTATIVE")]
    red=pair[pair.f5_redundant]
    fig,axes=plt.subplots(2,3,figsize=(14,9))
    bias_metrics=[("expected_heavy_atom_count","Heavy atoms",(0,80)),("entry_resolution","Resolution (Å)",(0.5,5)),("ligand_rscc","Ligand RSCC",(.4,1.02)),
                  ("ligand_rsr","Ligand RSR",(0,.6)),("pocket_residue_count","Pocket residues",(0,100)),("receptor_total_declared_length","Receptor length",(0,2000))]
    for ax,(metric,label,xlim) in zip(axes.flat,bias_metrics):
        for df,l,c in [(rep,"Representative",COLORS["green"]),(red,"Redundant",COLORS["gray"])]:
            x,y=ecdf(df[metric]); ax.plot(x,y,label=f"{l} n={len(x):,}",color=c)
        ax.set_xlim(*xlim); ax.set_xlabel(label); ax.set_ylabel("ECDF"); ax.set_title(label,loc="left",fontweight="bold"); style_ax(ax)
    axes[0,0].legend(fontsize=8); fig.suptitle("Representative-selection bias scout",x=.06,ha="left",fontsize=15,fontweight="bold")
    fig.subplots_adjust(left=.08, right=.98, bottom=.08, top=.90, hspace=.42, wspace=.20)
    p=savefig(fig,SUPP,"supp04_representative_bias")
    inventory.append(inv("S04","Representative-selection bias scout","Does representative selection shift size or quality?",
                         "F5 representatives vs redundant members","heavy atoms, resolution, RSCC, RSR, pocket and receptor size","READY",p[0],"SUPPLEMENTARY"))

    # Figure 5: final benchmark characterization.
    char_metrics=[("expected_heavy_atom_count","Ligand heavy atoms",(0,80)),("formula_weight","Molecular weight (Da)",(0,1000)),("rotatable_bond_count","Rotatable bonds",(0,25)),
                  ("ring_count","Ring count",(0,10)),("pocket_residue_count","6 Å pocket residues",(0,100)),("qualifying_contact_count","Qualifying atomic contacts",(0,100)),
                  ("entry_resolution","Resolution (Å)",(.5,5)),("ligand_rscc","Ligand RSCC",(.4,1.02)),("ligand_rsr","Ligand RSR",(0,.6))]
    fig,axes=plt.subplots(3,3,figsize=(14,12.4))
    for ax,(metric,label,xlim) in zip(axes.flat,char_metrics):
        x,y=ecdf(final[metric]); ax.plot(x,y,color=COLORS["purple"],lw=1.8)
        vals=pd.to_numeric(final[metric],errors="coerce").dropna(); med=vals.median() if len(vals) else np.nan
        if np.isfinite(med): ax.axvline(med,color=COLORS["orange"],lw=1,ls="--",label=f"median {med:.2g}")
        ax.set_xlim(*xlim); ax.set_xlabel(label); ax.set_ylabel("ECDF"); ax.set_title(f"{label} (n={len(vals):,})",loc="left",fontweight="bold"); ax.legend(fontsize=7); style_ax(ax)
    fig.suptitle("Final retained benchmark characterization (N=158,226)",x=.055,ha="left",fontsize=15,fontweight="bold")
    fig.subplots_adjust(left=.08, right=.98, bottom=.07, top=.92, hspace=.55, wspace=.20)
    p=savefig(fig,MAIN,"fig05_final_characterization")
    inventory.append(inv("F05","Final benchmark characterization","What protein, ligand, complex and experimental space does the final benchmark cover?",
                         "158,226 F5-retained cases","ligand descriptors, pocket/contact size, resolution, RSCC, RSR","READY",p[0],"KEEP"))

    # Identity recurrence supplementary.
    fig,axes=plt.subplots(1,3,figsize=(14,4.2))
    compfreq=final.component_id.value_counts(); pdbfreq=final.pdb_id.value_counts(); scaff=final.murcko_scaffold_smiles.replace("",np.nan).value_counts()
    for ax,series,title in [(axes[0],compfreq,"Ligand CCD frequency"),(axes[1],pdbfreq,"PDB frequency"),(axes[2],scaff,"Murcko scaffold frequency")]:
        vals=np.sort(series.to_numpy()); y=np.arange(1,len(vals)+1)/len(vals); ax.loglog(vals,y,color=COLORS["purple"]); ax.set_xlabel("Cases per identity"); ax.set_ylabel("CDF of identities"); ax.set_title(title,loc="left",fontweight="bold"); style_ax(ax)
    fig.suptitle("Identity recurrence in the final benchmark",x=.06,ha="left",fontsize=14,fontweight="bold")
    p=savefig(fig,SUPP,"supp05_identity_recurrence")
    inventory.append(inv("S05","Identity recurrence","How concentrated are cases by PDB, ligand identity and scaffold?",
                         "158,226 F5-retained cases","CCD/PDB/Murcko recurrence","READY",p[0],"SUPPLEMENTARY"))

    # Figure 5c: chemical-identity concentration (potential benchmark-composition issue).
    fig,axes=plt.subplots(2,2,figsize=(13.5,8.3))
    top=compfreq.head(15).sort_values()
    axes[0,0].barh(top.index,top.values,color=COLORS["purple"]); axes[0,0].set_xlabel("Final retained cases"); axes[0,0].set_title("Most frequent CCD identities",loc="left",fontweight="bold"); style_ax(axes[0,0],grid=False)
    ordered=compfreq.sort_values(ascending=False).to_numpy(); cumulative=np.cumsum(ordered)/ordered.sum(); share=np.arange(1,len(ordered)+1)/len(ordered)
    axes[0,1].plot(share*100,cumulative*100,color=COLORS["red"]); axes[0,1].plot([0,100],[0,100],color=COLORS["gray"],lw=.8)
    axes[0,1].set_xlabel("Most frequent ligand identities included (%)"); axes[0,1].set_ylabel("Final cases explained (%)"); axes[0,1].set_title("Ligand-identity concentration curve",loc="left",fontweight="bold"); style_ax(axes[0,1])
    share_rows=[]
    for label,df in [("F3 HIGH+GOOD",hg),("F4 PASS",f4pass),("F5 FINAL",final)]:
        vc=df.component_id.value_counts(); share_rows.append((label,vc.head(3).sum()/len(df)*100,vc.head(6).sum()/len(df)*100,vc.head(30).sum()/len(df)*100))
    x=np.arange(3); w=.24
    for i,(name,color) in enumerate([("Top 3",COLORS["red"]),("Top 6",COLORS["orange"]),("Top 30",COLORS["blue"])]):
        axes[1,0].bar(x+(i-1)*w,[r[i+1] for r in share_rows],w,label=name,color=color)
    axes[1,0].set_xticks(x,[r[0] for r in share_rows],rotation=12); axes[1,0].set_ylabel("Population explained (%)"); axes[1,0].set_title("Concentration survives filtering",loc="left",fontweight="bold"); axes[1,0].legend(); style_ax(axes[1,0])
    rep_component=pair.set_index("pair_id").component_id.to_dict(); g=groups[["representative_pair_id","group_size"]].copy(); g["component_id"]=g.representative_pair_id.map(rep_component)
    gtop=g.groupby("component_id").agg(groups=("group_size","size"),redundant_cases=("group_size",lambda v:int((v-1).sum())),largest_group=("group_size","max")).sort_values("redundant_cases",ascending=False).head(12).sort_values("redundant_cases")
    axes[1,1].barh(gtop.index,gtop.redundant_cases,color=COLORS["gray"]); axes[1,1].set_xlabel("Redundant cases in multi-member groups"); axes[1,1].set_title("Identities driving Filter 5 redundancy",loc="left",fontweight="bold"); style_ax(axes[1,1],grid=False)
    fig.suptitle("Ligand-identity concentration — priority composition audit",x=.055,ha="left",fontsize=15,fontweight="bold")
    p=savefig(fig,MAIN,"fig05c_ligand_identity_concentration")
    pd.DataFrame(share_rows,columns=["population","top3_share_pct","top6_share_pct","top30_share_pct"]).to_csv(DERIVED/"ligand_identity_concentration.tsv",sep="\t",index=False)
    gtop.reset_index().to_csv(DERIVED/"ligand_identity_redundancy_drivers.tsv",sep="\t",index=False)
    inventory.append(inv("F05C","Ligand-identity concentration","Is the final corpus dominated by a small number of CCD identities?",
                         "F3 HIGH+GOOD, F4 PASS and 158,226 F5-retained cases","top CCD share, concentration curve, group redundancy drivers","READY",p[0],"KEEP"))

    # Figure 6: quality map.
    quality_metrics=["entry_resolution","entry_r_work","entry_r_free","entry_r_free_minus_r_work","ligand_rscc","ligand_rsr","ligand_mean_occupancy","pocket_mean_rscc","pocket_mean_rsr","pocket_mean_occupancy","direct_binding_rsrz_median"]
    strata=[("F3 HIGH",pair.f3_status.eq("FILTER3_HIGH_QUALITY")),("F3 GOOD",pair.f3_status.eq("FILTER3_GOOD_QUALITY")),("F4 REJECT",pair.f4_reject),("F5 FINAL",pair.f5_retained)]
    avail=np.zeros((len(strata),len(quality_metrics)))
    den=np.zeros_like(avail,dtype=int)
    for i,(_,mask) in enumerate(strata):
        sub=pair[mask]; den[i,:]=len(sub)
        for j,m in enumerate(quality_metrics): avail[i,j]=sub[m].notna().mean()*100 if len(sub) else np.nan
    fig,axes=plt.subplots(2,2,figsize=(14,10.5))
    im=axes[0,0].imshow(avail,aspect="auto",vmin=0,vmax=100,cmap="viridis")
    axes[0,0].set_xticks(range(len(quality_metrics)),[m.replace("entry_","").replace("ligand_","").replace("pocket_","pocket ") for m in quality_metrics],rotation=52,ha="right",fontsize=8)
    axes[0,0].set_yticks(range(len(strata)),[s[0] for s in strata]); axes[0,0].set_title("Metric availability (%)",loc="left",fontweight="bold")
    for i in range(avail.shape[0]):
        for j in range(avail.shape[1]): axes[0,0].text(j,i,f"{avail[i,j]:.0f}",ha="center",va="center",fontsize=7,color="white" if avail[i,j]<60 else "black")
    fig.colorbar(im,ax=axes[0,0],fraction=.035,pad=.02)
    for metric,label,xlim,ax in [("ligand_rscc","Ligand RSCC",(.4,1.02),axes[0,1]),("entry_resolution","Resolution (Å)",(.5,5),axes[1,0])]:
        for q,c in [("HIGH",COLORS["blue"]),("GOOD",COLORS["orange"])]:
            x,y=ecdf(hg.loc[hg.f3_quality.eq(q),metric]); ax.plot(x,y,label=q,color=c)
        ax.set_xlim(*xlim); ax.set_xlabel(label); ax.set_ylabel("ECDF"); ax.legend(); ax.set_title(f"HIGH vs GOOD: {label}",loc="left",fontweight="bold"); style_ax(ax)
    hb=axes[1,1].hexbin(hg.entry_resolution,hg.ligand_rscc,gridsize=55,mincnt=1,bins="log",cmap="viridis")
    axes[1,1].set_xlim(.5,5); axes[1,1].set_ylim(.4,1.02); axes[1,1].set_xlabel("Resolution (Å)"); axes[1,1].set_ylabel("Ligand RSCC"); axes[1,1].set_title("RSCC × resolution density",loc="left",fontweight="bold"); fig.colorbar(hb,ax=axes[1,1],label="log count")
    fig.suptitle("Ground-truth quality map",x=.055,ha="left",fontsize=15,fontweight="bold")
    fig.subplots_adjust(left=.08, right=.95, bottom=.07, top=.91, hspace=.72, wspace=.25)
    p=savefig(fig,MAIN,"fig06_ground_truth_quality")
    inventory.append(inv("F06","Ground-truth quality map","How complete and distinct are the frozen quality strata?",
                         "F3/F4/F5 frozen pair strata","availability, resolution, R factors, RSCC/RSR/occupancy, binding-residue RSRZ","READY",p[0],"KEEP"))

    # Figure 8: descriptive removal heatmaps.
    fig,axes=plt.subplots(1,3,figsize=(15.5,4.8))
    pocket_bins=pd.cut(hg.pocket_residue_count,[-np.inf,10,20,30,40,60,np.inf],labels=["≤10","11–20","21–30","31–40","41–60",">60"])
    def rate_matrix(df,row,col,outcome):
        num=df.assign(_row=row,_col=col).groupby(["_row","_col"],observed=False)[outcome].agg(["sum","count"])
        return (num["sum"]/num["count"]*100).unstack(),num["count"].unstack()
    mat,nmat=rate_matrix(hg,hg.heavy_atom_bin,pocket_bins,"f4_reject")
    draw_heat(axes[0],mat,nmat,"P(F4 REJECT)","Heavy-atom bin","Pocket residues")
    f4pbins=pd.cut(f4pass.pocket_residue_count,[-np.inf,10,20,30,40,60,np.inf],labels=["≤10","11–20","21–30","31–40","41–60",">60"])
    mat2,nmat2=rate_matrix(f4pass,f4pass.heavy_atom_bin,f4pbins,"f5_redundant")
    draw_heat(axes[1],mat2,nmat2,"P(F5 REDUNDANT)","Heavy-atom bin","Pocket residues")
    qbins=pd.cut(hg.ligand_rscc,[0,.8,.85,.9,.95,1.01],labels=["<.80",".80–.85",".85–.90",".90–.95","≥.95"])
    rbins=pd.cut(hg.entry_resolution,[0,1.5,2,2.5,3,10],labels=["≤1.5","1.5–2.0","2.0–2.5","2.5–3.0",">3.0"])
    mat3,nmat3=rate_matrix(hg,qbins,rbins,"f5_retained")
    draw_heat(axes[2],mat3,nmat3,"P(final retained)","Ligand RSCC","Resolution (Å)")
    fig.suptitle("What predicts downstream removal? Descriptive rates only",x=.04,ha="left",fontsize=15,fontweight="bold")
    p=savefig(fig,MAIN,"fig08_downstream_removal_predictors")
    inventory.append(inv("F08","Downstream removal predictors","Which measured strata are associated with F4 rejection and F5 redundancy?",
                         "F3 HIGH+GOOD and F4 PASS denominators, separately","ligand size × pocket size; RSCC × resolution","READY",p[0],"KEEP"))

    # Heavy-atom audit supplementary, reused from final full_02 data.
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    hd=heavy_dist.copy(); axes[0].bar(hd.heavy_atom_bin.astype(str),hd.source_ligand_count,color=COLORS["blue"]); axes[0].set_yscale("log"); axes[0].tick_params(axis="x",rotation=30); axes[0].set_ylabel("F2 source ligands, log scale"); axes[0].set_title("Frozen CCD heavy-atom census",loc="left",fontweight="bold"); style_ax(axes[0])
    flow=[8543,4768,2590,2470,1019]; labs=["F2 placements","F3 input","F4 input","F4 PASS/F5 input","Current P4 input"]
    axes[1].plot(range(len(flow)),flow,marker="o",color=COLORS["orange"]); axes[1].set_xticks(range(len(flow)),labs,rotation=25,ha="right"); axes[1].set_ylabel("Cases with <3 heavy atoms"); axes[1].set_title("Known minimum-complexity lineage",loc="left",fontweight="bold"); style_ax(axes[1])
    fig.suptitle("Ligand minimum pose complexity — recorded issue, not a rule change",x=.05,ha="left",fontsize=14,fontweight="bold")
    p=savefig(fig,SUPP,"supp01_heavy_atom_census")
    inventory.append(inv("S01","Ligand minimum-complexity census","How many very-small ligands survive downstream?",
                         "Frozen full_02 heavy-atom audit","CCD expected heavy atoms and official lineage counts","READY",p[0],"SUPPLEMENTARY"))

    return inventory


def inv(fid,title,question,population,metrics,status,output,recommendation):
    return {"figure_id":fid,"title":title,"question":question,"population":population,"metrics":metrics,
            "status":status,"source":"frozen benchmark_1.0 + derived read-only joins","output_file":output,"recommendation":recommendation}


def draw_heat(ax, mat, nmat, title, ylabel, xlabel):
    arr=mat.to_numpy(dtype=float)
    im=ax.imshow(arr,aspect="auto",vmin=0,vmax=np.nanpercentile(arr,95) if np.isfinite(arr).any() else 100,cmap="magma")
    ax.set_xticks(range(len(mat.columns)),[str(x) for x in mat.columns],rotation=35,ha="right")
    ax.set_yticks(range(len(mat.index)),[str(x) for x in mat.index]); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title+" (%)",loc="left",fontweight="bold")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if np.isfinite(arr[i,j]):
                ax.text(j,i,f"{arr[i,j]:.1f}\nn={int(nmat.iloc[i,j]):,}",ha="center",va="center",fontsize=6.5,color="white")


def main():
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing Figure Scout: {OUT}")
    for path in [DERIVED,MAIN,SUPP,CROSS,RESEARCH,SCRIPTS,LOGS,QC,FUTURE]: path.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlepad":7,"figure.dpi":120})
    log=[]
    def note(msg): log.append(f"{utc()}\t{msg}"); print(msg,flush=True)

    note("load Filter 2 source and placement ledgers")
    source_cols=["pdb_id","source_ligand_instance_id","resolved_ccd_id","component_id","formula_weight","formal_charge","expected_heavy_atom_count","fragment_count","element_set","carbon_atom_count","contains_metal"]
    source=pd.read_csv(F2/"provisional_source_ligands.tsv.gz",sep="\t",usecols=source_cols,low_memory=False)
    source["source_ligand_instance_id"]=norm(source.source_ligand_instance_id); source["resolved_ccd_id"]=norm(source.resolved_ccd_id)
    as_num(source,["formula_weight","formal_charge","expected_heavy_atom_count","fragment_count","carbon_atom_count"])
    source["heavy_atom_bin"]=heavy_bin(source.expected_heavy_atom_count)
    placements=pd.read_csv(F2/"ligand_assembly_logical_placements.tsv.gz",sep="\t",usecols=["assembly_ligand_placement_id","source_ligand_instance_id"],low_memory=False).rename(columns={"assembly_ligand_placement_id":"ligand_assembly_placement_id"})
    placements["ligand_assembly_placement_id"]=norm(placements.ligand_assembly_placement_id); placements["source_ligand_instance_id"]=norm(placements.source_ligand_instance_id)

    note("derive RDKit descriptors from frozen Processing 2 CCD graph SMILES")
    topo=ds.dataset(str(P2/"ligand_topology_validation"),format="parquet").to_table(columns=["component_id","canonical_smiles_from_ccd_graph"]).to_pandas()
    descriptors,descriptor_qc=build_rdkit_descriptors(topo)
    descriptors.to_parquet(DERIVED/"ccd_rdkit_descriptors.parquet",index=False,compression="zstd")
    del topo

    note("load Processing 3 and Filter 3 pair ledgers")
    p3=ds.dataset(str(P3/"provisional_pairs"),format="parquet").to_table(columns=["pair_id","ligand_assembly_placement_id","receptor_chain_instance_ids","receptor_chain_count"]).to_pandas()
    f3_cols=["pair_id","ligand_assembly_placement_id","pdb_id","assembly_id","model_id","component_id","receptor_chain_instance_ids","receptor_chain_count","experimental_method_class","entry_resolution","entry_r_work","entry_r_free","entry_r_free_minus_r_work","ligand_rscc","ligand_rsr","ligand_mean_occupancy","pocket_mean_rscc","pocket_mean_rsr","pocket_mean_occupancy","pocket_missing_backbone_heavy_atom_count","direct_binding_missing_sidechain_heavy_atom_count","nonbinding_pocket_missing_sidechain_heavy_atom_count","critical_pocket_gap","pocket_gap_warning","filter3_v2_terminal_status","decision","reason_codes","warning_codes"]
    pair=ds.dataset(str(F3/"filter3_pair_quality_v2"),format="parquet").to_table(columns=f3_cols).to_pandas()
    for c in ["pair_id","ligand_assembly_placement_id","pdb_id","component_id"]: pair[c]=norm(pair[c])

    note("aggregate 6 Å pocket residues, atomic contacts and binding-residue RSRZ")
    pockets=aggregate_pocket(P3/"pair_pocket_residues"); contacts=aggregate_contacts(P3/"qualifying_atomic_contacts"); bq=aggregate_binding_quality(F3/"binding_residue_quality_v2")
    pockets.to_parquet(DERIVED/"pair_pocket_aggregates.parquet",index=False,compression="zstd")
    contacts.to_parquet(DERIVED/"placement_contact_aggregates.parquet",index=False,compression="zstd")
    bq.to_parquet(DERIVED/"pair_binding_quality_aggregates.parquet",index=False,compression="zstd")

    note("join official placement/source lineage")
    pmeta=placements.merge(source[["source_ligand_instance_id","resolved_ccd_id","formula_weight","formal_charge","expected_heavy_atom_count","fragment_count","element_set","carbon_atom_count","contains_metal"]],on="source_ligand_instance_id",how="left",validate="many_to_one",indicator="_source_join")
    pair=pair.merge(pmeta,on="ligand_assembly_placement_id",how="left",validate="one_to_one",indicator="_placement_join")
    pair=pair.merge(descriptors,left_on="resolved_ccd_id",right_on="component_id",how="left",validate="many_to_one",suffixes=("","_descriptor"),indicator="_descriptor_join")
    pair=pair.merge(pockets,on="pair_id",how="left",validate="one_to_one",indicator="_pocket_join")
    pair=pair.merge(contacts,on="ligand_assembly_placement_id",how="left",validate="one_to_one",indicator="_contact_join")
    pair=pair.merge(bq,on="pair_id",how="left",validate="one_to_one",indicator="_binding_quality_join")

    note("derive receptor sequence-length aggregates from frozen Filter 1 chain ledger")
    chains=pd.read_csv(F1/"filter_1_receptor_chain_instances.tsv.gz",sep="\t",usecols=["chain_instance_id","declared_sequence_length","observed_residue_count"])
    chains["chain_instance_id"]=norm(chains.chain_instance_id); cmap=chains.set_index("chain_instance_id")[["declared_sequence_length","observed_residue_count"]]
    exploded=pair[["pair_id","receptor_chain_instance_ids"]].copy(); exploded["chain_instance_id"]=exploded.receptor_chain_instance_ids.fillna("").str.split(","); exploded=exploded.explode("chain_instance_id")
    exploded=exploded.join(cmap,on="chain_instance_id")
    cagg=exploded.groupby("pair_id",as_index=False).agg(receptor_total_declared_length=("declared_sequence_length","sum"),receptor_max_chain_length=("declared_sequence_length","max"),receptor_total_observed_residues=("observed_residue_count","sum"),receptor_chain_metadata_available=("declared_sequence_length","count"))
    pair=pair.merge(cagg,on="pair_id",how="left",validate="one_to_one")

    note("load and join detailed Filter 4 and final Filter 5 state")
    f4=pd.read_csv(F4S5/"01_filter4_final_pair_inventory.tsv.gz",sep="\t",low_memory=False)
    f4["pair_id"]=norm(f4.pair_id); as_num(f4,["n_external_instances","n_external_ligand_6A","n_external_pocket_6A","n_direct_contact_instances","n_direct_contact_units","n_ligand_heavy_atoms_contacted_4A","fraction_ligand_heavy_atoms_contacted_4A","binding_residue_count","n_crystal_bridged_binding_residues","fraction_binding_residues_crystal_bridged","n_binding_residue_contacting_external_instances"])
    step4=pd.read_csv(F4S4/"04_pair_binding_residue_contact_inventory.tsv.gz",sep="\t",usecols=["candidate_pair_id","ligand_heavy_atom_count","min_external_ligand_distance_A","min_external_binding_residue_distance_A"],low_memory=False).rename(columns={"candidate_pair_id":"pair_id"})
    step4["pair_id"]=norm(step4.pair_id); as_num(step4,["ligand_heavy_atom_count","min_external_ligand_distance_A","min_external_binding_residue_distance_A"])
    f4=f4.merge(step4,on="pair_id",how="left",validate="one_to_one",indicator="_step4_detail_join")
    f5=pd.read_csv(F5/"01_filter5_final_pair_inventory.tsv.gz",sep="\t",low_memory=False); f5["pair_id"]=norm(f5.pair_id); as_num(f5,["group_size"])
    groups=pd.read_csv(F5/"04_filter5_equivalence_groups.tsv.gz",sep="\t",low_memory=False); members=pd.read_csv(F5/"05_filter5_group_members.tsv.gz",sep="\t",low_memory=False)
    as_num(groups,["group_size"])
    f4keep=f4[["pair_id","filter4_decision","filter4_reason","reject_stage","n_external_instances","n_external_ligand_6A","n_external_pocket_6A","n_direct_contact_instances","n_direct_contact_units","n_ligand_heavy_atoms_contacted_4A","fraction_ligand_heavy_atoms_contacted_4A","binding_residue_count","n_crystal_bridged_binding_residues","fraction_binding_residues_crystal_bridged","n_binding_residue_contacting_external_instances","ligand_heavy_atom_count","min_external_ligand_distance_A","min_external_binding_residue_distance_A"]]
    pair=pair.merge(f4keep,on="pair_id",how="left",validate="one_to_one",indicator="_f4_join")
    pair=pair.merge(f5[["pair_id","filter3_quality_class","step1_block_id","step2_exact_site_group_id","equivalence_group_id","group_size","representative_pair_id","filter5_final_status","filter5_final_reason"]],on="pair_id",how="left",validate="one_to_one",indicator="_f5_join")
    pair["f3_status"]=pair.filter3_v2_terminal_status; pair["f3_quality"]=np.select([pair.f3_status.eq("FILTER3_HIGH_QUALITY"),pair.f3_status.eq("FILTER3_GOOD_QUALITY")],["HIGH","GOOD"],default="OTHER")
    pair["f3_hg"]=pair.f3_quality.isin(QUALITY_ORDER); pair["f4_pass"]=pair.filter4_decision.eq("PASS"); pair["f4_reject"]=pair.filter4_decision.eq("REJECT")
    pair["f5_status"]=pair.filter5_final_status; pair["f5_redundant"]=pair.f5_status.eq("F5_REDUNDANT_EQUIVALENT_CASE")
    pair["f5_retained"]=pair.f5_status.isin(["F5_RETAIN_UNIQUE","F5_RETAIN_REPRESENTATIVE","F5_REVIEW_RETAIN"])
    pair["heavy_atom_bin"]=heavy_bin(pair.expected_heavy_atom_count)
    pair["pdb_recurrence"]=pair.pdb_id.map(pair[pair.f4_pass].pdb_id.value_counts()).fillna(0)
    pair["ligand_identity_recurrence"]=pair.component_id.map(pair[pair.f4_pass].component_id.value_counts()).fillna(0)
    pair.to_parquet(DERIVED/"unified_pair_figure_scout.parquet",index=False,compression="zstd")
    pair[pair.f5_retained].to_parquet(DERIVED/"final_retained_characterization.parquet",index=False,compression="zstd")

    note("write quality availability and distribution summaries")
    metrics=["expected_heavy_atom_count","formula_weight","rotatable_bond_count","ring_count","aromatic_ring_count","hbd_count","hba_count","fraction_csp3","pocket_residue_count","receptor_chain_count","receptor_total_declared_length","direct_binding_residue_rows","qualifying_contact_count","entry_resolution","entry_r_work","entry_r_free","entry_r_free_minus_r_work","ligand_rscc","ligand_rsr","ligand_mean_occupancy","pocket_mean_rscc","pocket_mean_rsr","pocket_mean_occupancy","direct_binding_rsrz_median"]
    summaries=[]
    for name,sub in [("F3 HIGH+GOOD",pair[pair.f3_hg]),("F4 PASS",pair[pair.f4_pass]),("F5 FINAL",pair[pair.f5_retained]),("F4 REJECT",pair[pair.f4_reject])]: summaries+=qsummary(sub,name,metrics)
    pd.DataFrame(summaries).to_csv(DERIVED/"metric_availability_and_summary.tsv",sep="\t",index=False)

    note("generate candidate figures")
    heavy_dist=pd.read_csv(HEAVY_AUDIT/"01_heavy_atom_distribution.tsv",sep="\t")
    inventory=make_figures(pair,source,f4,f5,groups,members,heavy_dist)

    note("select real F4 and F5 examples")
    examples=[]
    for label,mask,sortcol,asc in [
        ("clean",f4.filter4_reason.eq("NO_CRYSTALLOGRAPHIC_NEIGHBOR"),"n_external_instances",True),
        ("direct_ligand_contact",f4.filter4_reason.eq("DIRECT_LIGAND_CRYSTAL_CONTACT"),"fraction_ligand_heavy_atoms_contacted_4A",False),
        ("pocket_mediated",f4.filter4_reason.eq("BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT"),"fraction_binding_residues_crystal_bridged",False),
    ]:
        x=f4[mask].sort_values(sortcol,ascending=asc).head(5).copy(); x.insert(0,"example_class",label); examples.append(x)
    pd.concat(examples,ignore_index=True).to_csv(DERIVED/"real_f4_structure_example_candidates.tsv",sep="\t",index=False)
    inventory.append(inv("F03M","Crystal-packing molecular examples","What do clean, direct and pocket-mediated cases look like structurally?","15 real frozen F4 cases selected for rendering","native structures and crystal mates","NEED_RENDERER","","NEED_DATA"))

    # Real equivalence-group examples as membership graphs.
    targets=[]
    for target_size in [3,10,int(groups.group_size.max())]:
        row=groups.iloc[(groups.group_size-target_size).abs().argsort().iloc[0]]; targets.append(row)
    fig,axes=plt.subplots(1,3,figsize=(15,4.8))
    group_rows=[]
    for ax,row in zip(axes,targets):
        gid=row.equivalence_group_id; m=members[members.equivalence_group_id.eq(gid)]
        g=nx.Graph(); rep=str(row.representative_pair_id)
        for pid in m.pair_id.astype(str): g.add_node(pid,role="rep" if pid==rep else "member")
        for pid in m.pair_id.astype(str):
            if pid!=rep: g.add_edge(rep,pid)
        pos=nx.spring_layout(g,seed=24301); cols=[COLORS["green"] if g.nodes[n]["role"]=="rep" else COLORS["gray"] for n in g.nodes]
        nx.draw_networkx_nodes(g,pos,node_color=cols,node_size=60,ax=ax); nx.draw_networkx_edges(g,pos,alpha=.45,width=.7,ax=ax)
        ax.set_title(f"{gid}\nsize={len(m):,}, ligand={row.ligand_exact_id}",fontsize=9); ax.axis("off")
        for x in m.pair_id: group_rows.append({"equivalence_group_id":gid,"pair_id":x,"representative_pair_id":rep,"group_size":len(m)})
    fig.suptitle("Real strict-equivalence groups (membership view)",x=.04,ha="left",fontsize=14,fontweight="bold")
    fig.subplots_adjust(left=.03, right=.99, bottom=.04, top=.80, wspace=.12)
    p=savefig(fig,SUPP,"supp04_equivalence_group_examples")
    pd.DataFrame(group_rows).to_csv(DERIVED/"equivalence_group_example_members.tsv",sep="\t",index=False)
    inventory.append(inv("S04G","Real equivalence-group examples","What do small, medium and largest frozen groups look like?","Three real F5 groups","group membership and representative","READY_MEMBERSHIP_VIEW",p[0],"SUPPLEMENTARY"))

    # Upstream preparation scout.
    manifest=ds.dataset(str(P2/"structure_preparation_manifest"),format="parquet").to_table(columns=["object_type","preparation_status","decision","reason_code"]).to_pandas()
    prep=manifest.groupby(["object_type","preparation_status","decision"],dropna=False).size().reset_index(name="count")
    prep.to_csv(DERIVED/"processing2_preparation_status.tsv",sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(10,4.8)); plot=prep[prep.object_type.isin(["ligand_source_instance","receptor_chain_instance"])]
    labels=(plot.object_type.str.replace("_"," ")+"\n"+plot.preparation_status.str.replace("_"," ")).tolist()
    ax.barh(labels[::-1],plot["count"].to_numpy()[::-1],color=[COLORS["green"] if x=="PASS" else COLORS["orange"] for x in plot.decision][::-1]); ax.set_xlabel("Objects"); ax.set_title("Processing 2 preparation outcomes",loc="left",fontweight="bold"); style_ax(ax,grid=False)
    p=savefig(fig,SUPP,"supp02_processing2_preparation")
    inventory.append(inv("S02","Processing 2 preparation outcomes","Where do incomplete/mapping-review objects arise before pairs?","P2 source ligand and receptor ledgers","complete, missing-heavy-atom, mapping review, coordinate-ready","READY",p[0],"SUPPLEMENTARY"))

    # Placeholder inventory rows for cross benchmark and future docking.
    inventory.append(inv("F07","Cross-benchmark landscape","How does scale, scope and overlap compare to public benchmarks?","Official public benchmark metadata","scale, methodology, entry/chemical/instance/system overlap","RESEARCH_IN_PROGRESS","","NEED_DATA"))
    inventory.append(inv("F09","Future docking evaluation","How do methods perform under frozen strata?","No docking predictions currently available","RMSD, PB-valid, failures, runtime","BLOCKED_FUTURE_DATA","","FUTURE"))
    pd.DataFrame(inventory).to_csv(OUT/"figure_inventory.csv",index=False)

    note("write future docking schema and analysis skeleton")
    docking_cols=["pair_id","case_id","method","method_version","seed","rank","rmsd_A","pb_valid","runtime_seconds","run_status","failure_code","crystal_packing_stratum","quality_stratum","pocket_completeness_stratum"]
    pd.DataFrame(columns=docking_cols).to_csv(FUTURE/"future_docking_predictions_schema.csv",index=False)
    (FUTURE/"README.md").write_text("# Future docking figure\n\nStatus: `BLOCKED_FUTURE_DATA`. No performance values were simulated. Populate the schema only with real predictions.\n")
    (FUTURE/"analyze_future_docking.py").write_text("#!/usr/bin/env python3\n# Skeleton only: validate real prediction rows, then compute Top-1/Best-of-N RMSD and PB-valid metrics with explicit denominators.\nraise SystemExit('BLOCKED_FUTURE_DATA: provide real docking predictions')\n")

    note("run strict QC")
    f1e=pd.read_csv(F1/"filter_1_receptor_qualified_entries.tsv.gz",sep="\t",usecols=["pdb_id"])
    f1a=pd.read_csv(F1/"filter_1_receptor_qualified_assemblies.tsv.gz",sep="\t",usecols=["pdb_id","assembly_id"])
    f1c=pd.read_csv(F1/"filter_1_receptor_chain_instances.tsv.gz",sep="\t",usecols=["chain_instance_id"])
    f1_source=pd.read_csv(ROOT/"filter_1_protein_receptor_qualification/full/filter_1_source_chains.tsv.gz",sep="\t",usecols=["is_polypeptide"])
    auxiliary=json.loads((ROOT/"auxiliary_entry_work_packages/builds/20260805_full_01/output/build_summary.json").read_text())
    prep_counts=prep.groupby(["object_type","preparation_status"])["count"].sum().to_dict()
    checks={
        "f1_entries_248037":len(f1e)==248037,"f1_assemblies_360611":len(f1a)==360611,
        "f1_source_polypeptide_chains_1073451":int(truth(f1_source.is_polypeptide).sum())==1073451,
        "f1_receptor_assembly_chain_instances_2145537":len(f1c)==2145537,
        "f2_source_852968":len(source)==852968,"f2_placements_1151324":len(placements)==1151324,
        "active_assembly_keys_234975":int(auxiliary["active_assembly_key_count"])==234975,
        "active_receptor_instances_834548":int(auxiliary["active_receptor_chain_instances"])==834548,
        "p2_complete_source_746509":int(prep_counts.get(("ligand_source_instance","COMPLETE"),-1))==746509,
        "p2_incomplete_source_100354":int(prep_counts.get(("ligand_source_instance","ATOM_INCOMPLETE"),-1))==100354,
        "p2_mapping_review_5103":int(prep_counts.get(("ligand_source_instance","ATOM_MAPPING_REVIEW"),-1))==5103,
        "p2_coordinate_ready_receptors_834222":int(prep_counts.get(("receptor_chain_instance","ASSEMBLY_READY"),-1))==834222,
        "p3_pairs_744580":len(pair)==744580,"f3_high_120297":int(pair.f3_status.eq("FILTER3_HIGH_QUALITY").sum())==120297,
        "f3_good_216115":int(pair.f3_status.eq("FILTER3_GOOD_QUALITY").sum())==216115,"f3_hg_336412":int(pair.f3_hg.sum())==336412,
        "f4_total_336412":len(f4)==336412,"f4_pass_241545":int(f4.filter4_decision.eq("PASS").sum())==241545,
        "f4_reject_94865":int(f4.filter4_decision.eq("REJECT").sum())==94865,"f4_review_2":int(f4.filter4_decision.eq("REVIEW").sum())==2,
        "f5_total_241545":len(f5)==241545,"f5_redundant_83319":int(f5.filter5_final_status.eq("F5_REDUNDANT_EQUIVALENT_CASE").sum())==83319,
        "f5_retained_158226":int(f5.filter5_final_status.ne("F5_REDUNDANT_EQUIVALENT_CASE").sum())==158226,
        "f5_groups_32188":len(groups)==32188,"f5_group_members_close":int(groups.group_size.sum())==len(members),
        "source_key_unique":source.source_ligand_instance_id.is_unique,"placement_key_unique":placements.ligand_assembly_placement_id.is_unique,
        "pair_key_unique":pair.pair_id.is_unique,"f4_key_unique":f4.pair_id.is_unique,"f5_key_unique":f5.pair_id.is_unique,
        "placement_source_join_missing_zero":int(pmeta._source_join.ne("both").sum())==0,
        "f3_placement_join_missing_zero":int(pair._placement_join.ne("both").sum())==0,
        "f4_join_expected_only_hg":int(pair.loc[pair.f3_hg,"_f4_join"].ne("both").sum())==0 and int(pair.loc[~pair.f3_hg,"_f4_join"].eq("both").sum())==0,
        "f5_join_expected_only_f4pass":int(pair.loc[pair.f4_pass,"_f5_join"].ne("both").sum())==0 and int(pair.loc[~pair.f4_pass,"_f5_join"].eq("both").sum())==0,
        "heavy_atom_audit_pass":json.loads((HEAVY_AUDIT/"validation_report.json").read_text())["validation_pass"] is True,
    }
    qc={"run_id":RUN_ID,"generated_at":utc(),"validation_pass":bool(all(checks.values())),"checks":{k:bool(v) for k,v in checks.items()},
        "descriptor_qc":descriptor_qc,"join_diagnostics":{
            "descriptor_missing_pairs":int(pair._descriptor_join.ne("both").sum()),"pocket_aggregate_missing_pairs":int(pair._pocket_join.ne("both").sum()),
            "contact_aggregate_missing_pairs":int(pair._contact_join.ne("both").sum()),"binding_quality_missing_pairs":int(pair._binding_quality_join.ne("both").sum()),
        },"findings":[
            {"id":"UNIT_LABEL_F1_1073451","severity":"HIGH","status":"CORRECTED_IN_FIGURE_SCOUT","detail":"1,073,451 is deposited/source polypeptide chains, not assembly chain instances; receptor assembly-chain instances are 2,145,537."},
            {"id":"P2_SPARSE_PARQUET_SCHEMA","severity":"HIGH","status":"AVOIDED","detail":"Some sparse P2 datasets include zero-row _empty schema fragments; default dataset schema discovery can hide non-empty rows. Figure Scout uses the unified structure_preparation_manifest."},
            {"id":"F4_STEP4_VALIDATION_FALSE","severity":"MEDIUM","status":"EXPLICIT_REVIEW_ROUTE","detail":"Step 4 summary records validation_pass=False because two BA-equivalence unresolved cases flow explicitly to final REVIEW; no silent repair was applied."},
            {"id":"F5_LOG_HASH_MISMATCH","severity":"MEDIUM","status":"SCIENTIFIC_OUTPUTS_INTACT","detail":"F5 SHA256SUMS mismatch is limited to logs/run.log written after freeze checksum creation; all scientific outputs match."},
            {"id":"F5_FREEZE_ANCHOR_WEAK","severity":"MEDIUM","status":"DOCUMENTED","detail":"F5 _FROZEN.json does not hash-anchor its manifest/SHA256SUMS, weaker than earlier stages."},
            {"id":"POSEBUSTERS_APPROVED_NOT_FROZEN","severity":"HIGH","status":"SEPARATED","detail":"156,621 is an approved counterfactual target awaiting versioned rerun; authoritative frozen F5 remains 158,226 and is used for all current figures."},
            {"id":"LIGAND_IDENTITY_CONCENTRATION","severity":"HIGH","status":"FIGURE_ADDED","detail":"SO4, GOL and EDO account for 55,372/158,226 final cases (35.0%); the six most frequent CCD identities (adding PO4, ACT and PEG) account for 42.44%. Heavy-atom >=3 would not resolve this composition issue."},
            {"id":"F4_CONDITIONAL_SEVERITY_DENOMINATORS","severity":"HIGH","status":"EXPLICIT","detail":"F4 direct-contact and binding-residue severity metrics are conditional on Step 3/4 eligibility and missing values were not filled with zero."},
            {"id":"LIGAND_RSRZ_UNAVAILABLE","severity":"MEDIUM","status":"DROP_PANEL","detail":"Ligand RSRZ is available for only 22/336,412 F3 HIGH+GOOD pairs; binding/pocket-residue RSRZ has high residue-level coverage and is kept with its own denominator."},
            {"id":"POSEBUSTERS_SPLIT_PROVENANCE","severity":"HIGH","status":"DOCUMENTED","detail":"Complete frozen PoseBusters evidence is split across 224,994 old-source and 58,943 v2-new-source rows (283,937 unique source ligands). Reading only v2 posebusters_new_results would create false missingness."},
            {"id":"F5_IDENTITY_MISSING_29","severity":"MEDIUM","status":"DENOMINATOR_REQUIRED","detail":"29/241,545 F5 input pairs have unresolved ligand_exact_id under explicit F5_STEP1_REVIEW; chemical-identity analyses use denominator 241,516 when this field is required."},
            {"id":"F5_EDGE_VS_GROUP_GRAPH","severity":"MEDIUM","status":"SEPARATED","detail":"The 2,529,689 strict-edge table includes 206,891 edges not used within final groups; group density must not be computed from the entire edge table without final-group membership filtering."},
            {"id":"LEGACY_ABSOLUTE_PATHS","severity":"LOW","status":"REMAPPED","detail":"Some provenance records retain /root/autodl-tmp paths; Figure Scout uses the verified private-server root."}
        ]}
    (QC/"qc_report.json").write_text(json.dumps(qc,indent=2,default=str)+"\n")

    source_paths=[F1/"filter_1_release_summary.json",F2/"provisional_source_ligands.tsv.gz",F2/"ligand_assembly_logical_placements.tsv.gz",P2/"release_summary.json",P3/"release_summary.json",F3.parent/"_FROZEN.json",F4S5.parent/"_FROZEN.json",F5.parent/"_FROZEN.json",HEAVY_AUDIT/"validation_report.json"]
    provenance={"run_id":RUN_ID,"mode":"READ_ONLY_EXPLORATORY_FIGURE_SCOUT","created_at":utc(),"formal_result_changes":False,
                "inputs":[{"path":str(p),"size_bytes":p.stat().st_size,"sha256":sha256(p)} for p in source_paths if p.exists()]}
    (OUT/"input_provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")

    # Baseline narrative reports; research supplement is finalized separately.
    (OUT/"README.md").write_text(f"""# Benchmark Figure Scout\n\nRun: `{RUN_ID}`\n\nA read-only exploratory visualization audit over frozen `benchmark_1.0`. No status, representative, cutoff, or frozen population was changed.\n\nGenerated candidate figures: {len(inventory)} inventory entries. QC: **{'PASS' if qc['validation_pass'] else 'FAIL'}**.\n\nThe existing `<3 heavy atoms` issue is visualized as a recorded supplementary scout and is not treated as a rule change.\n""")
    (OUT/"figure_review.md").write_text("""# Figure review — preliminary\n\n## Strongest main-text candidates\n\n- F01 keeps receptor, ligand and pair units separate while showing pair attrition.\n- F03 makes crystal packing a distinctive scientific contribution rather than a generic cleaning step.\n- F04 shows that redundancy is highly group-structured and should not be described as simple random downsampling.\n- F06 gives quality provenance and denominator-aware missingness.\n- F08 is descriptive, not causal, and should be merged with F02 if the paper has a tight figure budget.\n\n## Supplementary\n\nThe Processing 2 preparation ledger, ligand minimum-complexity census, detailed crystal severity, representative-bias scout and identity recurrence belong in supplementary material unless one becomes a central unexpected result.\n\n## Needs additional work\n\nCross-benchmark overlap needs official manifests and mapping at explicitly named granularity. Crystal molecular examples have been selected from frozen cases but require a reproducible crystal-mate renderer. No docking performance figure is currently possible.\n\nThe construction panel must label 1,073,451 as **source polypeptide chains**. It must not be called assembly chain instances; the latter count is 2,145,537.\n""")
    (OUT/"data_gap_report.md").write_text("""# Data gap report — preliminary\n\n## Available now\n\nFrozen ligand descriptors, receptor chain lengths, pocket/contact sizes, experimental quality metrics, crystal-packing severity and Filter 5 grouping support Figures F01–F06 and F08.\n\n## Derivable from existing frozen data\n\nRDKit descriptors and Murcko scaffolds were derived only from frozen Processing 2 CCD-graph SMILES. Pocket/contact counts were aggregated from frozen Processing 3 rows.\n\n## Missing or external\n\nSpecies, protein/domain family and PDB release year are not present in the unified frozen pair ledger and need entry/SIFTS metadata. Cross-benchmark fine-grained overlap needs official manifests. Crystal molecular panels need a renderer.\n\nProcessing 2 sparse review datasets contain zero-row `_empty` schema fragments that can mislead automatic PyArrow schema discovery. Use the frozen manifest or an explicitly unified schema; never infer a zero review count from the first fragment.\n\nThe approved PoseBusters-strict 156,621 target is not yet a frozen rerun. Current figures therefore use the authoritative 158,226 frozen F5 population and keep the future target separate.\n\n## Must wait\n\nDocking RMSD, PB-validity, failure rate and runtime require real docking predictions and remain `BLOCKED_FUTURE_DATA`.\n""")
    (OUT/"qc_report.md").write_text(f"""# QC report\n\nOverall: **{'PASS' if qc['validation_pass'] else 'FAIL'}**.\n\nAll requested headline counts, primary keys and expected cross-stage joins were checked. See `qc/qc_report.json` for machine-readable checks and explicit missingness.\n\nImportant non-silent findings are listed in `qc/qc_report.json`: the F1 source-vs-assembly chain unit correction, Processing 2 sparse-Parquet schema trap, two explicit F4 review cases, an F5 non-scientific log checksum mismatch, weaker F5 freeze anchoring, legacy absolute paths, and the separation between frozen 158,226 and the approved-but-not-rerun 156,621 target.\n""")
    (LOGS/"build.log").write_text("\n".join(log)+"\n")
    shutil.copy2(Path(__file__),SCRIPTS/Path(__file__).name)

    # Figure provenance one row per candidate, including denominators and sources.
    prov=[]
    for row in inventory:
        prov.append({**row,"n_total":"see population","n_metric_available":"see data_derived/metric_availability_and_summary.tsv","n_missing":"explicit in metric summary",
                     "denominator":row["population"],"source_files":"input_provenance.json","source_columns":row["metrics"],"filters_applied":"frozen statuses only","script":f"scripts/{Path(__file__).name}","generated_timestamp":utc()})
    pd.DataFrame(prov).to_csv(OUT/"figure_provenance.csv",index=False)

    # Final output manifest (cross-benchmark research may append in finalization).
    rows=[]
    for p in sorted(x for x in OUT.rglob("*") if x.is_file()):
        rows.append({"relative_path":str(p.relative_to(OUT)),"size_bytes":p.stat().st_size,"sha256":sha256(p)})
    pd.DataFrame(rows).to_csv(OUT/"output_manifest.tsv",sep="\t",index=False)
    note(f"complete: validation={qc['validation_pass']} outputs={len(rows)}")
    print(json.dumps({"output":str(OUT),"validation_pass":qc["validation_pass"],"figure_inventory":len(inventory),"descriptor_qc":descriptor_qc},indent=2))
    if not qc["validation_pass"]: raise SystemExit(2)


if __name__=="__main__":
    main()
