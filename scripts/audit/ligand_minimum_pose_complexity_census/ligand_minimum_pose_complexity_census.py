#!/usr/bin/env python3
"""Read-only Filter 2 v3 heavy-atom census and downstream lineage audit."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq


ROOT = Path("/home/linx/data/youcq/autodl-tmp/benchmark_1.0")
RUN_ID = "20260823_full_02"
OUT = ROOT / "audits/ligand_minimum_pose_complexity_census" / RUN_ID

F2 = ROOT / "filter_2_ligand_qualification_v3/runs/20260804_full_01"
F2_LEGACY = ROOT / "filter_2_ligand_qualification"
F3 = ROOT / "filter_03_ground_truth_structure_quality_v2/runs/20260814_full_01"
F4 = ROOT / "filter_04_crystal_packing_influence/step_05_final_crystal_packing_decision/runs/step05_full_v1"
F5 = ROOT / "filter_05_equivalent_redocking_case/step_05_strict_equivalent_grouping_and_representative_selection/runs/step05_full_v1"
P4 = ROOT / "processing_04_docking_ready_case_construction/runs/p4_full_v1_0_1"
RECON = ROOT / "reconciliation/filter3_filter4_filter5_p4_20260822"

SOURCE_PATH = F2 / "output/provisional_source_ligands.tsv.gz"
NOMAP_PATH = F2 / "output/ligand_assembly_no_mapping.tsv.gz"
PLACEMENT_PATH = F2 / "output/ligand_assembly_logical_placements.tsv.gz"
CCD_PATH = F2_LEGACY / "references/ccd_component_cache.tsv.gz"
F3_PATH = F3 / "output/filter3_pair_quality_v2"
F4_PATH = F4 / "output/01_filter4_final_pair_inventory.tsv.gz"
F5_FINAL_PATH = F5 / "output/01_filter5_final_pair_inventory.tsv.gz"
F5_GROUP_PATH = F5 / "output/04_filter5_equivalence_groups.tsv.gz"
F5_MEMBER_PATH = F5 / "output/05_filter5_group_members.tsv.gz"
P4_INPUT_PATH = P4 / "input/full_case_inventory.parquet"
P4_STATUS_PATH = P4 / "output/processing4_case_inventory.parquet"


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False)


def write_json(payload: dict, name: str) -> None:
    (OUT / name).write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def normalize_string(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str)


def heavy_bin(value) -> str:
    if pd.isna(value):
        return "missing_or_invalid"
    try:
        v = int(value)
    except (TypeError, ValueError):
        return "missing_or_invalid"
    if v < 0:
        return "missing_or_invalid"
    if v <= 5:
        return str(v)
    if v <= 10:
        return "6-10"
    if v <= 20:
        return "11-20"
    if v <= 40:
        return "21-40"
    return ">40"


BIN_ORDER = ["0", "1", "2", "3", "4", "5", "6-10", "11-20", "21-40", ">40", "missing_or_invalid"]


def charge_category(value) -> str:
    if pd.isna(value):
        return "missing_or_unknown"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "missing_or_unknown"
    if v < 0:
        return "negative"
    if v > 0:
        return "positive"
    return "zero"


def bool_text(value) -> str:
    if pd.isna(value):
        return "missing"
    return str(value).strip().lower()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"refusing to overwrite existing audit directory: {OUT}")
    OUT.mkdir(parents=True)

    source_cols = [
        "pdb_id", "selected_model_id", "component_id", "source_ligand_instance_id",
        "resolved_ccd_id", "ccd_identity_status", "ccd_name", "ccd_type", "formula",
        "formula_weight", "formal_charge", "standard_total_atom_count",
        "expected_heavy_atom_count", "element_set", "carbon_atom_count", "fragment_count",
        "contains_metal", "descriptor_availability", "terminal_route", "decision",
    ]
    source = pd.read_csv(SOURCE_PATH, sep="\t", usecols=source_cols, low_memory=False)
    for col in ["pdb_id", "selected_model_id", "component_id", "source_ligand_instance_id", "resolved_ccd_id"]:
        source[col] = normalize_string(source[col])
    for col in ["formula_weight", "formal_charge", "standard_total_atom_count", "expected_heavy_atom_count", "carbon_atom_count", "fragment_count"]:
        source[col] = pd.to_numeric(source[col], errors="coerce")
    source["heavy_atom_bin"] = source.expected_heavy_atom_count.map(heavy_bin)
    source["heavy_class"] = np.select(
        [source.expected_heavy_atom_count.lt(3), source.expected_heavy_atom_count.ge(3)],
        ["<3", ">=3"],
        default="missing",
    )

    placements = pd.read_csv(PLACEMENT_PATH, sep="\t", usecols=[
        "pdb_id", "assembly_id", "selected_model_id", "source_ligand_instance_id",
        "component_id", "assembly_ligand_placement_id", "mapping_status",
    ], low_memory=False)
    placements = placements.rename(columns={"assembly_ligand_placement_id": "ligand_assembly_placement_id"})
    for col in ["pdb_id", "assembly_id", "selected_model_id", "source_ligand_instance_id", "component_id", "ligand_assembly_placement_id"]:
        placements[col] = normalize_string(placements[col])
    placements["active_assembly_key"] = placements.pdb_id + "|" + placements.assembly_id + "|" + placements.selected_model_id

    nomap = pd.read_csv(NOMAP_PATH, sep="\t", usecols=["source_ligand_instance_id", "mapping_status"])
    nomap["source_ligand_instance_id"] = normalize_string(nomap.source_ligand_instance_id)

    source_ids = set(source.source_ligand_instance_id)
    mapped_ids = set(placements.source_ligand_instance_id)
    nomap_ids = set(nomap.source_ligand_instance_id)
    source["mapped"] = source.source_ligand_instance_id.isin(mapped_ids)
    source["no_mapping"] = source.source_ligand_instance_id.isin(nomap_ids)

    pmeta = placements[["source_ligand_instance_id", "ligand_assembly_placement_id", "active_assembly_key"]].merge(
        source[["source_ligand_instance_id", "expected_heavy_atom_count", "heavy_atom_bin", "heavy_class", "resolved_ccd_id"]],
        on="source_ligand_instance_id", how="left", validate="many_to_one", indicator=True,
    )

    # Distribution ledger.
    distribution_rows = []
    for bin_name in BIN_ORDER:
        s = source[source.heavy_atom_bin.eq(bin_name)]
        p = pmeta[pmeta.heavy_atom_bin.eq(bin_name)]
        distribution_rows.append({
            "heavy_atom_bin": bin_name,
            "source_ligand_count": len(s),
            "mapped_source_ligand_count": int(s.mapped.sum()),
            "no_mapping_source_ligand_count": int(s.no_mapping.sum()),
            "ligand_placement_count": len(p),
            "unique_pdb_count": s.pdb_id.nunique(),
            "unique_active_assembly_key_count": p.active_assembly_key.nunique(),
            "unique_ccd_component_count": s.resolved_ccd_id.replace("", np.nan).nunique(),
        })
    dist = pd.DataFrame(distribution_rows)
    write_tsv(dist, "01_heavy_atom_distribution.tsv")

    lt3_source = source[source.expected_heavy_atom_count.lt(3)].copy()
    ge3_source = source[source.expected_heavy_atom_count.ge(3)].copy()
    lt3_ids = set(lt3_source.source_ligand_instance_id)
    lt3_pmeta = pmeta[pmeta.source_ligand_instance_id.isin(lt3_ids)].copy()
    lt3_placement_ids = set(lt3_pmeta.ligand_assembly_placement_id)
    dist_summary = {
        "source_total": len(source),
        "mapped_source_total": len(mapped_ids),
        "no_mapping_source_total": len(nomap_ids),
        "placement_total": len(placements),
        "heavy_atom_bins": {r["heavy_atom_bin"]: int(r["source_ligand_count"]) for r in distribution_rows},
        "lt3": {
            "source_ligands": len(lt3_source),
            "mapped_source_ligands": int(lt3_source.mapped.sum()),
            "no_mapping_source_ligands": int(lt3_source.no_mapping.sum()),
            "placements": len(lt3_pmeta),
            "unique_pdb": lt3_source.pdb_id.nunique(),
            "unique_active_assembly_key": lt3_pmeta.active_assembly_key.nunique(),
            "unique_ccd": lt3_source.resolved_ccd_id.replace("", np.nan).nunique(),
        },
        "ge3": {"source_ligands": len(ge3_source)},
    }
    write_json(dist_summary, "02_heavy_atom_distribution_summary.json")

    # PDB and active-assembly co-occurrence: a small ligand never excludes its whole container.
    cooccurrence_rows = []
    for scope, frame, key in [
        ("PDB", source, "pdb_id"),
        ("ACTIVE_ASSEMBLY", pmeta, "active_assembly_key"),
    ]:
        flags = frame.groupby(key).heavy_class.agg(
            has_lt3=lambda values: bool((values == "<3").any()),
            has_ge3=lambda values: bool((values == ">=3").any()),
            has_missing=lambda values: bool((values == "missing").any()),
        ).reset_index()
        flags["composition"] = np.select(
            [flags.has_lt3 & flags.has_ge3, flags.has_lt3, flags.has_ge3],
            ["LT3_AND_GE3", "LT3_ONLY", "GE3_ONLY"],
            default="MISSING_ONLY",
        )
        for composition, group in flags.groupby("composition"):
            cooccurrence_rows.append({
                "scope": scope,
                "composition": composition,
                "container_count": len(group),
                "has_missing_heavy_atom_count": int(group.has_missing.sum()),
            })
    cooccurrence = pd.DataFrame(cooccurrence_rows).sort_values(["scope", "composition"])
    write_tsv(cooccurrence, "02b_pdb_active_assembly_cooccurrence.tsv")

    # Component-level census using Filter 2 v3 frozen metadata, with source/placement counts.
    component_fields = [
        "resolved_ccd_id", "ccd_name", "ccd_type", "formula", "formula_weight", "formal_charge",
        "standard_total_atom_count", "expected_heavy_atom_count", "element_set", "carbon_atom_count",
        "fragment_count", "contains_metal", "descriptor_availability",
    ]
    comp_base = source.sort_values("source_ligand_instance_id").drop_duplicates("resolved_ccd_id")[component_fields].copy()
    comp_counts = source.groupby("resolved_ccd_id", dropna=False).agg(
        source_ligand_count=("source_ligand_instance_id", "size"),
        mapped_source_ligand_count=("mapped", "sum"),
        unique_pdb_count=("pdb_id", "nunique"),
    ).reset_index()
    placement_counts = pmeta.groupby("resolved_ccd_id", dropna=False).agg(
        placement_count=("ligand_assembly_placement_id", "size"),
        unique_active_assembly_key_count=("active_assembly_key", "nunique"),
    ).reset_index()
    comps = comp_base.merge(comp_counts, on="resolved_ccd_id", validate="one_to_one").merge(
        placement_counts, on="resolved_ccd_id", how="left", validate="one_to_one"
    )
    comps["placement_count"] = comps.placement_count.fillna(0).astype(int)
    comps["unique_active_assembly_key_count"] = comps.unique_active_assembly_key_count.fillna(0).astype(int)
    comps = comps.rename(columns={
        "resolved_ccd_id": "component_id", "ccd_name": "component_name", "ccd_type": "component_type",
        "fragment_count": "connected_component_count",
    })
    component_order = [
        "component_id", "component_name", "component_type", "formula", "formula_weight", "formal_charge",
        "standard_total_atom_count", "expected_heavy_atom_count", "element_set", "carbon_atom_count",
        "contains_metal", "connected_component_count", "descriptor_availability", "source_ligand_count",
        "mapped_source_ligand_count", "placement_count", "unique_pdb_count", "unique_active_assembly_key_count",
    ]
    comps = comps[component_order]
    lt3_comps = comps[comps.expected_heavy_atom_count.lt(3)].sort_values(["source_ligand_count", "component_id"], ascending=[False, True])
    write_tsv(lt3_comps, "03_lt3_component_full.tsv")
    for count, filename in [(1, "04_heavy_atom_1_components.tsv"), (2, "05_heavy_atom_2_components.tsv")]:
        frame = comps[comps.expected_heavy_atom_count.eq(count)].sort_values(["source_ligand_count", "component_id"], ascending=[False, True])
        write_tsv(frame, filename)
        write_tsv(frame.head(50), f"top_50_heavy_atom_{count}_components.tsv")
    for count, filename in [(3, "06_heavy_atom_3_top100.tsv"), (4, "07_heavy_atom_4_top100.tsv"), (5, "08_heavy_atom_5_top100.tsv")]:
        frame = comps[comps.expected_heavy_atom_count.eq(count)].sort_values(["source_ligand_count", "component_id"], ascending=[False, True]).head(100)
        write_tsv(frame, filename)

    # Element, type, charge, MW and connectivity descriptions.
    element_rows = []
    for h in [1, 2, 3]:
        subset = source[source.expected_heavy_atom_count.eq(h)]
        for element, group in subset.groupby(subset.element_set.fillna("missing_or_unknown")):
            element_rows.append({
                "expected_heavy_atom_count": h, "element_set": element,
                "source_ligand_count": len(group), "placement_count": len(pmeta[pmeta.source_ligand_instance_id.isin(set(group.source_ligand_instance_id))]),
                "unique_ccd_count": group.resolved_ccd_id.nunique(), "unique_pdb_count": group.pdb_id.nunique(),
                "carbon_zero_source_count": int(group.carbon_atom_count.fillna(-1).eq(0).sum()),
                "carbon_ge1_source_count": int(group.carbon_atom_count.ge(1).sum()),
            })
    element_df = pd.DataFrame(element_rows).sort_values(["expected_heavy_atom_count", "source_ligand_count"], ascending=[True, False])
    write_tsv(element_df[element_df.expected_heavy_atom_count.isin([1, 2])], "09_lt3_element_set_distribution.tsv")
    write_tsv(element_df, "09b_heavy_atom_1_2_3_element_set_distribution.tsv")

    type_rows = []
    for component_type, group in lt3_source.groupby(lt3_source.ccd_type.fillna("missing_or_unknown")):
        gids = set(group.source_ligand_instance_id)
        type_rows.append({
            "component_type": component_type, "source_ligand_count": len(group),
            "placement_count": len(pmeta[pmeta.source_ligand_instance_id.isin(gids)]),
            "unique_ccd_count": group.resolved_ccd_id.nunique(), "unique_pdb_count": group.pdb_id.nunique(),
        })
    write_tsv(pd.DataFrame(type_rows).sort_values("source_ligand_count", ascending=False), "10_lt3_component_type_distribution.tsv")

    charge_rows = []
    lt3_source["charge_category"] = lt3_source.formal_charge.map(charge_category)
    lt3_source["exact_charge"] = lt3_source.formal_charge.map(lambda v: "missing_or_unknown" if pd.isna(v) else str(int(v)) if float(v).is_integer() else str(v))
    for (category, exact), group in lt3_source.groupby(["charge_category", "exact_charge"], dropna=False):
        gids = set(group.source_ligand_instance_id)
        charge_rows.append({
            "charge_category": category, "exact_formal_charge": exact, "source_ligand_count": len(group),
            "placement_count": len(pmeta[pmeta.source_ligand_instance_id.isin(gids)]),
            "unique_ccd_count": group.resolved_ccd_id.nunique(), "unique_pdb_count": group.pdb_id.nunique(),
        })
    write_tsv(pd.DataFrame(charge_rows).sort_values(["charge_category", "source_ligand_count"], ascending=[True, False]), "11_lt3_formal_charge_distribution.tsv")

    mw_rows = []
    for label, subset in [("heavy_atom_1", source[source.expected_heavy_atom_count.eq(1)]), ("heavy_atom_2", source[source.expected_heavy_atom_count.eq(2)]), ("heavy_atom_ge3", ge3_source)]:
        values = subset.formula_weight.dropna().astype(float)
        mw_rows.append({
            "group": label, "count": len(values), "missing": len(subset) - len(values),
            "min": values.min() if len(values) else np.nan, "p25": values.quantile(.25) if len(values) else np.nan,
            "median": values.median() if len(values) else np.nan, "p75": values.quantile(.75) if len(values) else np.nan,
            "p90": values.quantile(.90) if len(values) else np.nan, "p95": values.quantile(.95) if len(values) else np.nan,
            "max": values.max() if len(values) else np.nan,
        })
    write_tsv(pd.DataFrame(mw_rows), "12_heavy_atom_mw_summary.tsv")
    write_tsv(lt3_comps[["component_id", "formula", "formula_weight", "expected_heavy_atom_count", "source_ligand_count", "placement_count"]], "12b_lt3_formula_weight_detail.tsv")

    connected_rows = []
    for label, subset in [("heavy_atom_lt3", lt3_source), ("heavy_atom_ge3", ge3_source)]:
        cats = pd.Series(np.where(subset.fragment_count.eq(1), "connected_component_count_1", np.where(subset.fragment_count.gt(1), "connected_component_count_gt1", "missing_or_invalid")), index=subset.index)
        for category, group in subset.groupby(cats):
            gids = set(group.source_ligand_instance_id)
            connected_rows.append({
                "heavy_atom_group": label, "connected_component_category": category,
                "source_ligand_count": len(group), "placement_count": len(pmeta[pmeta.source_ligand_instance_id.isin(gids)]),
                "unique_ccd_count": group.resolved_ccd_id.nunique(), "unique_pdb_count": group.pdb_id.nunique(),
            })
    write_tsv(pd.DataFrame(connected_rows), "13_connected_component_census.tsv")
    write_tsv(comps[comps.expected_heavy_atom_count.ge(3) & comps.connected_component_count.gt(1)].sort_values("source_ligand_count", ascending=False), "13b_ge3_multicomponent_ccd_full.tsv")

    unsupported_rows = []
    for label, subset in [("heavy_atom_lt3", lt3_source), ("heavy_atom_ge3", ge3_source)]:
        for desc, group in subset.groupby(subset.descriptor_availability.fillna("missing_or_unknown")):
            unsupported_rows.append({
                "heavy_atom_group": label, "field": "descriptor_availability", "value": desc,
                "source_ligand_count": len(group), "unique_ccd_count": group.resolved_ccd_id.nunique(),
            })
    write_tsv(pd.DataFrame(unsupported_rows), "13c_unsupported_chemistry_preaudit.tsv")

    # Pair lineage: official placement key only.
    placement_heavy = pmeta[["ligand_assembly_placement_id", "source_ligand_instance_id", "expected_heavy_atom_count", "heavy_class", "resolved_ccd_id"]].copy()
    placement_heavy = placement_heavy.rename(columns={"resolved_ccd_id": "ccd_component_id"})
    f3 = ds.dataset(str(F3_PATH), format="parquet").to_table(columns=[
        "pair_id", "ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id",
        "filter3_v2_terminal_status", "reason_codes", "warning_codes",
    ]).to_pandas()
    for col in ["pair_id", "ligand_assembly_placement_id", "pdb_id", "assembly_id", "model_id", "component_id"]:
        f3[col] = normalize_string(f3[col])
    f3j = f3.merge(placement_heavy, on="ligand_assembly_placement_id", how="left", validate="many_to_one", indicator="_placement_join")
    f3lt = f3j[f3j.expected_heavy_atom_count.lt(3)].copy()
    entered_lt3_placements = set(f3lt.ligand_assembly_placement_id)
    filter3_flow_rows = [{
        "scope": "mapped_lt3_placement_ledger", "status": "ENTERED_FILTER3",
        "pair_count": len(f3lt), "unique_placement_count": len(entered_lt3_placements), "unique_pdb_count": f3lt.pdb_id.nunique(),
    }, {
        "scope": "mapped_lt3_placement_ledger", "status": "DID_NOT_ENTER_FILTER3",
        "pair_count": 0, "unique_placement_count": len(lt3_placement_ids - entered_lt3_placements),
        "unique_pdb_count": lt3_pmeta[lt3_pmeta.ligand_assembly_placement_id.isin(lt3_placement_ids - entered_lt3_placements)].source_ligand_instance_id.map(dict(zip(source.source_ligand_instance_id, source.pdb_id))).nunique(),
    }]
    for status, group in f3lt.groupby("filter3_v2_terminal_status"):
        filter3_flow_rows.append({
            "scope": "filter3_pairs", "status": status, "pair_count": len(group),
            "unique_placement_count": group.ligand_assembly_placement_id.nunique(), "unique_pdb_count": group.pdb_id.nunique(),
        })
    write_tsv(pd.DataFrame(filter3_flow_rows), "14_lt3_filter3_flow.tsv")
    write_tsv(pd.DataFrame(filter3_flow_rows), "heavy_atom_lt3_filter3_flow.tsv")
    write_tsv(f3lt[["pair_id", "ligand_assembly_placement_id", "source_ligand_instance_id", "pdb_id", "assembly_id", "model_id", "component_id", "expected_heavy_atom_count", "filter3_v2_terminal_status", "reason_codes", "warning_codes"]], "14b_lt3_filter3_pairs.tsv")

    pair_heavy = f3j[["pair_id", "ligand_assembly_placement_id", "source_ligand_instance_id", "expected_heavy_atom_count", "heavy_class", "ccd_component_id"]].copy()
    f4 = pd.read_csv(F4_PATH, sep="\t", low_memory=False)
    f4["pair_id"] = normalize_string(f4.pair_id)
    f4j = f4.merge(pair_heavy, on="pair_id", how="left", validate="one_to_one", indicator="_f3_join")
    f4lt = f4j[f4j.expected_heavy_atom_count.lt(3)].copy()
    f4_flow = f4lt.groupby(["filter4_decision", "filter4_reason", "reject_stage"], dropna=False).agg(
        pair_count=("pair_id", "size"), unique_placement_count=("ligand_assembly_placement_id", "nunique")
    ).reset_index().sort_values("pair_count", ascending=False)
    write_tsv(f4_flow, "15_lt3_filter4_flow.tsv")
    write_tsv(f4lt[["pair_id", "ligand_assembly_placement_id", "source_ligand_instance_id", "expected_heavy_atom_count", "filter4_decision", "filter4_reason", "reject_stage"]], "15b_lt3_filter4_pairs.tsv")

    f5 = pd.read_csv(F5_FINAL_PATH, sep="\t", low_memory=False)
    f5["pair_id"] = normalize_string(f5.pair_id)
    f5j = f5.merge(pair_heavy, on="pair_id", how="left", validate="one_to_one", indicator="_f3_join")
    f5lt = f5j[f5j.expected_heavy_atom_count.lt(3)].copy()
    f5_flow = f5lt.groupby("filter5_final_status", dropna=False).agg(
        pair_count=("pair_id", "size"), unique_placement_count=("ligand_assembly_placement_id", "nunique"),
        unique_ccd_count=("ccd_component_id", "nunique")
    ).reset_index().sort_values("pair_count", ascending=False)
    write_tsv(f5_flow, "16_lt3_filter5_flow.tsv")

    members = pd.read_csv(F5_MEMBER_PATH, sep="\t", dtype=str).fillna("")
    members["pair_id"] = normalize_string(members.pair_id)
    member_heavy = members.merge(pair_heavy[["pair_id", "expected_heavy_atom_count", "ccd_component_id"]], on="pair_id", how="left", validate="many_to_one", indicator="_heavy_join")
    groups = pd.read_csv(F5_GROUP_PATH, sep="\t", dtype=str).fillna("")
    group_meta = groups.set_index("equivalence_group_id").to_dict("index")
    f5_index = f5.set_index("pair_id")
    rep_rows = []
    for gid, group in member_heavy.groupby("equivalence_group_id", sort=False):
        rep_pair = str(group.representative_pair_id.iloc[0])
        rep = group[group.pair_id.eq(rep_pair)]
        if rep.empty or not pd.to_numeric(rep.expected_heavy_atom_count, errors="coerce").lt(3).all():
            continue
        h = pd.to_numeric(group.expected_heavy_atom_count, errors="coerce")
        ge3 = int(h.ge(3).sum())
        lt3 = int(h.lt(3).sum())
        frow = f5_index.loc[rep_pair]
        meta = group_meta.get(gid, {})
        rep_rows.append({
            "equivalence_group_id": gid, "step1_block_id": frow.get("step1_block_id", ""),
            "step2_exact_site_group_id": meta.get("step2_exact_site_group_id", frow.get("step2_exact_site_group_id", "")),
            "current_representative_pair_id": rep_pair, "representative_component_id": rep.ccd_component_id.iloc[0],
            "representative_heavy_atom_count": int(rep.expected_heavy_atom_count.iloc[0]),
            "group_size": len(group), "ge3_surviving_member_count": ge3, "lt3_member_count": lt3,
            "impact_class": "A_RESELECT_REPRESENTATIVE" if ge3 > 0 else "B_GROUP_FULLY_REMOVED",
        })
    rep_impact = pd.DataFrame(rep_rows)
    if not rep_impact.empty:
        rep_impact = rep_impact.sort_values(["impact_class", "equivalence_group_id"])
    write_tsv(rep_impact, "17_filter5_lt3_representative_impact.tsv")
    write_tsv(rep_impact, "filter5_lt3_representative_impact.tsv")

    redundant = member_heavy[member_heavy.member_role.eq("REDUNDANT") & pd.to_numeric(member_heavy.expected_heavy_atom_count, errors="coerce").lt(3)].copy()
    rep_heavy_map = dict(zip(member_heavy.pair_id, pd.to_numeric(member_heavy.expected_heavy_atom_count, errors="coerce")))
    redundant["representative_heavy_atom_count"] = redundant.representative_pair_id.map(rep_heavy_map)
    redundant["representative_ge3"] = redundant.representative_heavy_atom_count.ge(3)
    redundant_impact = redundant[["equivalence_group_id", "pair_id", "ccd_component_id", "expected_heavy_atom_count", "representative_pair_id", "representative_heavy_atom_count", "representative_ge3"]].sort_values(["equivalence_group_id", "pair_id"])
    write_tsv(redundant_impact, "18_filter5_lt3_redundant_impact.tsv")

    # P4 current frozen flow and approved 156,621 target-membership reconciliation.
    p4input = pq.read_table(P4_INPUT_PATH, columns=["pair_id", "case_id", "component_id", "ligand_assembly_placement_id", "filter5_final_status"]).to_pandas()
    p4status = pq.read_table(P4_STATUS_PATH).to_pandas()
    p4input["pair_id"] = normalize_string(p4input.pair_id)
    p4status["pair_id"] = normalize_string(p4status.pair_id)
    p4j = p4input.merge(p4status, on=["pair_id", "case_id"], how="left", validate="one_to_one", indicator="_status_join").merge(
        pair_heavy, on="pair_id", how="left", validate="one_to_one", indicator="_f3_join"
    )
    p4lt = p4j[p4j.expected_heavy_atom_count.lt(3)].copy()
    p4lt["etkdg_attempted"] = ~p4lt.reason.fillna("").str.contains("frozen ligand atoms/bonds missing", regex=False)
    p4lt["etkdg_success"] = p4lt.status.eq("P4_DOCKING_READY")
    p4lt["etkdg_failure"] = p4lt.status.eq("P4_LIGAND_START_GENERATION_FAILED")
    p4_flow_rows = []
    for status, group in p4lt.groupby("status"):
        p4_flow_rows.append({
            "scope": "CURRENT_FROZEN_P4_158226", "status": status, "case_count": len(group),
            "etkdg_attempted": int(group.etkdg_attempted.sum()), "etkdg_success": int(group.etkdg_success.sum()),
            "etkdg_failure": int(group.etkdg_failure.sum()), "other_preparation_state": int((~group.etkdg_success & ~group.etkdg_failure).sum()),
        })
    candidate_path = RECON / "counterfactual_filter5_inventory.tsv"
    strict_target_ids = set()
    if candidate_path.exists():
        candidate = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
        strict_target_ids = set(candidate.loc[candidate.candidate_filter5_status.ne("F5_REDUNDANT_EQUIVALENT_CASE"), "pair_id"])
        strict_existing = p4lt[p4lt.pair_id.isin(strict_target_ids)]
        p4_flow_rows.append({
            "scope": "APPROVED_PB_STRICT_TARGET_MEMBERSHIP_156621_EXISTING_LT3_ONLY", "status": "EXISTING_IN_P4",
            "case_count": len(strict_existing), "etkdg_attempted": int(strict_existing.etkdg_attempted.sum()),
            "etkdg_success": int(strict_existing.etkdg_success.sum()), "etkdg_failure": int(strict_existing.etkdg_failure.sum()),
            "other_preparation_state": int((~strict_existing.etkdg_success & ~strict_existing.etkdg_failure).sum()),
        })
    write_tsv(pd.DataFrame(p4_flow_rows), "19_lt3_processing4_flow.tsv")

    failures = p4j[p4j.status.eq("P4_LIGAND_START_GENERATION_FAILED")].copy()
    source_component_lookup = comps.set_index("component_id")
    failure_rows = []
    for row in failures.itertuples(index=False):
        ccd = str(row.ccd_component_id)
        meta = source_component_lookup.loc[ccd] if ccd in source_component_lookup.index else pd.Series(dtype=object)
        failure_rows.append({
            "pair_id": row.pair_id, "pdb_id": str(row.pair_id).split("|")[1], "component_id": ccd,
            "expected_heavy_atom_count": row.expected_heavy_atom_count, "formula": meta.get("formula", ""),
            "formula_weight": meta.get("formula_weight", np.nan), "formal_charge": meta.get("formal_charge", np.nan),
            "element_set": meta.get("element_set", ""), "carbon_atom_count": meta.get("carbon_atom_count", np.nan),
            "connected_component_count": meta.get("connected_component_count", np.nan),
            "filter5_status": row.filter5_final_status, "p4_status": row.status, "p4_reason": row.reason,
        })
    failure_audit = pd.DataFrame(failure_rows).sort_values(["component_id", "pair_id"])
    write_tsv(failure_audit, "20_p4_etkdg_failures_chemistry_audit.tsv")
    write_tsv(failure_audit, "p4_etkdg_failures_chemistry_audit.tsv")

    retained_statuses = {"F5_RETAIN_UNIQUE", "F5_RETAIN_REPRESENTATIVE", "F5_REVIEW_RETAIN"}
    current_retained_lt3 = f5lt[f5lt.filter5_final_status.isin(retained_statuses)]
    impact_metrics = [
        ("filter2", "would_exclude_source_ligands", len(lt3_source), "source_ligand_instance_id"),
        ("filter2", "would_exclude_mapped_source_ligands", int(lt3_source.mapped.sum()), "source_ligand_instance_id"),
        ("filter2", "would_remove_placements", len(lt3_pmeta), "ligand_assembly_placement_id"),
        ("filter2", "would_affect_pdbs", lt3_source.pdb_id.nunique(), "pdb_id"),
        ("filter2", "would_affect_active_assembly_keys", lt3_pmeta.active_assembly_key.nunique(), "pdb_id+assembly_id+model_id"),
        ("filter3", "would_remove_filter3_inputs", len(f3lt), "pair_id"),
        ("filter4", "would_remove_filter4_inputs", len(f4lt), "pair_id"),
        ("filter4", "would_remove_filter4_pass", int(f4lt.filter4_decision.eq("PASS").sum()), "pair_id"),
        ("filter5", "would_remove_filter5_inputs", len(f5lt), "pair_id"),
        ("filter5", "would_remove_current_retained_cases_before_reselection", len(current_retained_lt3), "pair_id; not a projected final count"),
        ("filter5", "current_representatives_affected", len(rep_impact), "equivalence_group_id"),
        ("filter5", "groups_requiring_representative_reselection", int(rep_impact.impact_class.eq("A_RESELECT_REPRESENTATIVE").sum()) if not rep_impact.empty else 0, "equivalence_group_id"),
        ("filter5", "groups_entirely_removed", int(rep_impact.impact_class.eq("B_GROUP_FULLY_REMOVED").sum()) if not rep_impact.empty else 0, "equivalence_group_id"),
        ("processing4", "would_remove_current_frozen_p4_inputs", len(p4lt), "pair_id"),
        ("processing4", "current_etkdg_failures_heavy_atom_lt3", int(failure_audit.expected_heavy_atom_count.lt(3).sum()), "pair_id"),
        ("processing4", "current_etkdg_failures_heavy_atom_ge3", int(failure_audit.expected_heavy_atom_count.ge(3).sum()), "pair_id"),
    ]
    impact = pd.DataFrame(impact_metrics, columns=["stage", "metric", "value", "unit_or_note"])
    write_tsv(impact, "21_hypothetical_ge3_impact_summary.tsv")
    impact_json = {stage: dict(zip(group.metric, group.value.astype(int))) for stage, group in impact.groupby("stage")}
    impact_json["warning"] = "Filter 5 final count must not be calculated by simple subtraction; representative reselection was not executed."
    write_json(impact_json, "22_hypothetical_ge3_impact_summary.json")

    # Strict validation and provenance.
    checks = {
        "filter2_source_total_852968": len(source) == 852_968,
        "filter2_source_id_unique": source.source_ligand_instance_id.is_unique,
        "heavy_atom_bins_sum_852968": int(dist.source_ligand_count.sum()) == 852_968,
        "mapped_source_total_851966": len(mapped_ids) == 851_966,
        "no_mapping_source_total_1002": len(nomap_ids) == 1_002,
        "mapped_plus_nomap_closes_source": mapped_ids.isdisjoint(nomap_ids) and mapped_ids | nomap_ids == source_ids,
        "placement_total_1151324": len(placements) == 1_151_324,
        "placement_id_unique": placements.ligand_assembly_placement_id.is_unique,
        "placement_mapping_status_all_mapped": bool(placements.mapping_status.eq("mapped").all()),
        "placement_source_join_missing_zero": int(pmeta._merge.ne("both").sum()) == 0,
        "filter3_pair_id_unique": f3.pair_id.is_unique,
        "filter3_placement_join_missing_zero": int(f3j._placement_join.ne("both").sum()) == 0,
        "filter4_pair_id_unique": f4.pair_id.is_unique,
        "filter4_f3_join_missing_zero": int(f4j._f3_join.ne("both").sum()) == 0,
        "filter5_pair_id_unique": f5.pair_id.is_unique,
        "filter5_f3_join_missing_zero": int(f5j._f3_join.ne("both").sum()) == 0,
        "filter5_group_member_pair_unique": members.pair_id.is_unique,
        "filter5_group_member_heavy_join_missing_zero": int(member_heavy._heavy_join.ne("both").sum()) == 0,
        "p4_input_pair_id_unique": p4input.pair_id.is_unique,
        "p4_status_pair_id_unique": p4status.pair_id.is_unique,
        "p4_status_join_missing_zero": int(p4j._status_join.ne("both").sum()) == 0,
        "p4_f3_join_missing_zero": int(p4j._f3_join.ne("both").sum()) == 0,
        "p4_input_count_158226": len(p4input) == 158_226,
        "p4_etkdg_failure_count_37": len(failures) == 37,
        "unexpected_duplicate_zero": all([source.source_ligand_instance_id.is_unique, placements.ligand_assembly_placement_id.is_unique, f3.pair_id.is_unique, f4.pair_id.is_unique, f5.pair_id.is_unique, p4input.pair_id.is_unique]),
        "ambiguous_join_zero": all([
            source.source_ligand_instance_id.is_unique,
            placements.ligand_assembly_placement_id.is_unique,
            f3.pair_id.is_unique,
            f4.pair_id.is_unique,
            f5.pair_id.is_unique,
            members.pair_id.is_unique,
            p4input.pair_id.is_unique,
            p4status.pair_id.is_unique,
        ]),
        "pdb_cooccurrence_accounting_closes": int(cooccurrence.loc[cooccurrence.scope.eq("PDB"), "container_count"].sum()) == source.pdb_id.nunique(),
        "active_assembly_cooccurrence_accounting_closes": int(cooccurrence.loc[cooccurrence.scope.eq("ACTIVE_ASSEMBLY"), "container_count"].sum()) == pmeta.active_assembly_key.nunique(),
        "silent_drop_zero": all([
            len(source) == len(source_ids),
            len(placements) == len(pmeta),
            len(f3) == len(f3j),
            len(f4) == len(f4j),
            len(f5) == len(f5j),
            len(members) == len(member_heavy),
            len(p4input) == len(p4j),
        ]),
    }
    validation = {
        "run_id": RUN_ID, "validated_at": utc(), "validation_pass": all(checks.values()),
        "checks": checks,
        "counts": {
            "source": len(source), "mapped_source": len(mapped_ids), "no_mapping_source": len(nomap_ids),
            "placements": len(placements), "filter3_pairs": len(f3), "filter4_pairs": len(f4),
            "filter5_pairs": len(f5), "p4_cases": len(p4input), "p4_etkdg_failures": len(failures),
        },
        "join_diagnostics": {
            "placement_source_missing": int(pmeta._merge.ne("both").sum()),
            "filter3_placement_missing": int(f3j._placement_join.ne("both").sum()),
            "filter4_f3_missing": int(f4j._f3_join.ne("both").sum()),
            "filter5_f3_missing": int(f5j._f3_join.ne("both").sum()),
            "group_member_heavy_missing": int(member_heavy._heavy_join.ne("both").sum()),
            "p4_status_missing": int(p4j._status_join.ne("both").sum()),
            "p4_f3_missing": int(p4j._f3_join.ne("both").sum()),
        },
    }
    write_json(validation, "validation_report.json")

    provenance_paths = [
        SOURCE_PATH, NOMAP_PATH, PLACEMENT_PATH, CCD_PATH, F3 / "_FROZEN.json", F4 / "_FROZEN.json",
        F5 / "_FROZEN.json", P4 / "_FROZEN.json", F5_FINAL_PATH, F5_GROUP_PATH, F5_MEMBER_PATH,
        P4_INPUT_PATH, P4_STATUS_PATH,
    ]
    provenance = {
        "audit_run_id": RUN_ID, "created_at": utc(), "mode": "READ_ONLY_CENSUS_AND_DOWNSTREAM_IMPACT_AUDIT",
        "scientific_rule_changes": False, "formal_membership_changes": False,
        "heavy_atom_field": "Filter 2 v3 provisional_source_ligands.expected_heavy_atom_count",
        "lineage_keys": ["source_ligand_instance_id", "ligand_assembly_placement_id", "pair_id", "equivalence_group_id"],
        "inputs": [{"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)} for path in provenance_paths],
    }
    write_json(provenance, "input_provenance.json")

    readme = f"""# Filter 2 v3 ligand heavy-atom complexity census

Run: `{RUN_ID}`

This directory contains a read-only census and downstream impact audit. It does not change Filter 2/3/4/5, Processing 4, representative selection, or any frozen membership.

The only complexity field used is frozen Filter 2 v3 `expected_heavy_atom_count` derived from the CCD standard definition. Current observed atom count, AltLoc-selected atom count and RDKit heavy-atom count were not used as census substitutes.

Current downstream statistics refer to the existing frozen lineage. The separately approved 156,621-case PoseBusters-strict target is shown only as an additional membership-reconciliation scope where available; it is not represented as an already frozen P4 run.

Validation: **{'PASS' if validation['validation_pass'] else 'FAIL'}**.
"""
    (OUT / "README.md").write_text(readme)

    # Output manifest after all requested outputs exist.
    manifest = []
    for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in {"output_manifest.tsv", "SHA256SUMS"}):
        rows = None
        if path.suffix == ".tsv":
            with path.open("rb") as handle:
                rows = max(sum(1 for _ in handle) - 1, 0)
        manifest.append({"relative_path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path), "row_count": rows})
    write_tsv(pd.DataFrame(manifest), "output_manifest.tsv")
    with (OUT / "SHA256SUMS").open("w") as handle:
        for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.name}\n")

    print(json.dumps({
        "output": str(OUT), "validation_pass": validation["validation_pass"],
        "distribution": dist_summary, "impact": impact_json,
        "lt3_representatives": len(rep_impact),
        "groups_reselection": int(rep_impact.impact_class.eq("A_RESELECT_REPRESENTATIVE").sum()) if not rep_impact.empty else 0,
        "groups_removed": int(rep_impact.impact_class.eq("B_GROUP_FULLY_REMOVED").sum()) if not rep_impact.empty else 0,
    }, indent=2, default=str))

    if not validation["validation_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
