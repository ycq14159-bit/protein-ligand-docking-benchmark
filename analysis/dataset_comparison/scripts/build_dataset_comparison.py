#!/usr/bin/env python3
"""Build the auditable CROWN-style dataset comparison statistics.

The script is read-only with respect to the frozen benchmark root.  It writes
only below --output-root and derives every value in the Ours column from the
current, explicitly frozen final database-construction release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


FINAL_F4_RUN = "20260826_filter3_118255_strict_posebusters_01"
PROVENANCE_RUN = "20260826_filter4_pass_01"
F2_RUN = "20260825_dual_source_strict_01"
P2_RUN = "20260826_validation_provenance_required_01"
F3_STRICT_RUN = "20260826_processing3_176900_strict_posebusters_02"
F5_IDENTITY_RUN = "20260826_filter4_91860_rmsd1A_lexquality_v2"
F3_V2_RUN = "20260814_full_01"

STANDARD_PROTEIN_COMPONENTS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
    "VAL", "SEC",
}

PROPERTIES = [
    ("Dataset scope", "Total entries"),
    ("Dataset scope", "Unique PDB-CCD pairs"),
    ("Dataset scope", "Unique PDB IDs"),
    ("Dataset scope", "Unique UniProt IDs"),
    ("Dataset scope", "Unique CATH IDs"),
    ("Dataset scope", "Unique species"),
    ("Dataset scope", "Affinity annotations"),
    ("Ligand diversity", "Unique CCD IDs"),
    ("Ligand diversity", "Unique Murcko scaffolds"),
    ("Ligand diversity", "Ion ligands"),
    ("Ligand diversity", "Covalent ligands"),
    ("Ligand diversity", "Artifact ligands"),
    ("Structure issues", "Missing bonds"),
    ("Structure issues", "Steric overlaps"),
    ("Structure issues", "Unresolved ligand atoms"),
    ("Structure issues", "Unresolved pocket atoms"),
    ("Structure issues", "Non-standard pocket residues"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def mkdirs(root: Path) -> None:
    for name in ("scripts", "data", "tables", "qc", "logs"):
        (root / name).mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=columns, low_memory=False)


def require_unique(df: pd.DataFrame, key: str, label: str) -> None:
    if df[key].isna().any():
        raise RuntimeError(f"{label}: {key} has null values")
    dup = int(df[key].duplicated().sum())
    if dup:
        raise RuntimeError(f"{label}: {dup} duplicate {key} values")


def scan_partitioned_parquet(directory: Path, columns: list[str], pair_ids: set[str]) -> pd.DataFrame:
    chunks: list[pd.DataFrame] = []
    for path in sorted(directory.glob("**/*.parquet")):
        schema_names = set(pq.ParquetFile(path).schema_arrow.names)
        if not set(columns).issubset(schema_names):
            continue
        chunk = pq.ParquetFile(path).read(columns=columns).to_pandas()
        chunk = chunk[chunk["pair_id"].isin(pair_ids)]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        return pd.DataFrame(columns=columns)
    return pd.concat(chunks, ignore_index=True)


def normalize_list(values) -> list[str]:
    out: set[str] = set()
    for value in values:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        if isinstance(value, (list, tuple)):
            items = value
        else:
            text = str(value).strip()
            if not text:
                continue
            try:
                parsed = json.loads(text)
                items = parsed if isinstance(parsed, list) else [parsed]
            except Exception:
                items = [x.strip() for x in text.replace(";", ",").split(",")]
        out.update(str(x).strip() for x in items if str(x).strip())
    return sorted(out)


def build_scaffolds(smiles: pd.Series) -> tuple[dict[str, tuple[str | None, str]], dict[str, int]]:
    cache: dict[str, tuple[str | None, str]] = {}
    counts = defaultdict(int)
    for smi in sorted(set(x for x in smiles.dropna().astype(str) if x.strip())):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                cache[smi] = (None, "PARSE_FAILED")
            else:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                out = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)
                cache[smi] = (out or None, "OK" if out else "ACYCLIC_EMPTY")
        except Exception:
            cache[smi] = (None, "SCAFFOLD_FAILED")
        counts[cache[smi][1]] += 1
    return cache, dict(counts)


def metric_row(section: str, prop: str, value, status: str, definition: str,
               source: str, denominator: str, notes: str = "") -> dict:
    return {
        "section": section, "property": prop, "value": value, "status": status,
        "definition": definition, "source": source, "denominator": denominator,
        "notes": notes,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--output-root", required=True, type=Path)
    ap.add_argument("--analysis-timestamp", required=True,
                    help="Fixed ISO-8601 timestamp reused for deterministic reruns")
    ap.add_argument("--git-commit", required=True)
    args = ap.parse_args()
    root, out = args.data_root.resolve(), args.output_root.resolve()
    mkdirs(out)

    paths = {
        "f4_marker": root / "filter_04_crystal_packing_influence" / "database_runs" / FINAL_F4_RUN / "_FROZEN.json",
        "f4_pass": root / "filter_04_crystal_packing_influence" / "database_runs" / FINAL_F4_RUN / "output" / "02_filter4_pass_pairs.tsv.gz",
        "f5_identity": root / "filter_05_equivalent_redocking_case" / "runs" / F5_IDENTITY_RUN / "step01" / "output" / "01_step1_pair_inventory.tsv.gz",
        "p2": root / "processing_2_assembly_ready_structure_preparation" / "runs" / P2_RUN / "output" / "formal_ready_ligand_placements.parquet",
        "f2_placements": root / "filter_2_ligand_qualification_v4" / "runs" / F2_RUN / "output" / "02_retained_assembly_placements.tsv.gz",
        "f2_sources": root / "filter_2_ligand_qualification_v4" / "runs" / F2_RUN / "output" / "01_source_membership.tsv.gz",
        "pb": root / "filter_03_ground_truth_structure_quality_database" / "runs" / F3_STRICT_RUN / "output" / "posebusters_pair_evidence.parquet",
        "f3_v2": root / "filter_03_ground_truth_structure_quality_v2" / "runs" / F3_V2_RUN / "output" / "filter3_pair_quality_v2",
        "prov_summary": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "pair_protein_provenance_summary.parquet",
        "prov_chains": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "pair_receptor_chains.parquet",
        "uniprot_map": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "protein_uniprot_mappings.parquet",
        "uniprot_meta": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "protein_uniprot_reference_metadata.parquet",
        "source_chains": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "protein_source_chains.parquet",
        "function": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "protein_function_annotations.parquet",
        "pocket": root / "protein_provenance_annotation" / "runs" / PROVENANCE_RUN / "output" / "pair_pocket_residues.parquet",
    }
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError("Required frozen inputs missing:\n" + "\n".join(missing))

    marker = json.loads(paths["f4_marker"].read_text(encoding="utf-8"))
    if marker.get("database_role") != "FINAL_DATABASE_CONSTRUCTION_STAGE":
        raise RuntimeError("Filter 4 marker is not the explicit final database-construction stage")

    f4 = read_tsv(paths["f4_pass"])
    require_unique(f4, "pair_id", "Filter 4 PASS")
    if len(f4) != 91860 or set(f4["filter4_decision"]) != {"PASS"}:
        raise RuntimeError(f"Unexpected formal membership: rows={len(f4)}, decisions={set(f4['filter4_decision'])}")
    pair_ids = set(f4["pair_id"])

    identity_cols = ["pair_id", "resolved_ccd_id", "ligand_exact_id", "normalized_ccd_isomeric_smiles"]
    identity = read_tsv(paths["f5_identity"], identity_cols)
    identity = identity[identity["pair_id"].isin(pair_ids)].copy()
    require_unique(identity, "pair_id", "Filter 5 identity projection")
    if len(identity) != len(f4):
        raise RuntimeError("Filter 5 identity projection does not cover formal members")

    p2_cols = ["assembly_ligand_placement_id", "source_ligand_instance_id",
               "missing_heavy_atom_count", "topology_status", "validation_mapping_status"]
    p2 = pq.read_table(paths["p2"], columns=p2_cols).to_pandas()
    p2 = p2[p2["assembly_ligand_placement_id"].isin(set(f4["ligand_assembly_placement_id"]))]
    require_unique(p2, "assembly_ligand_placement_id", "Processing 2 formal placements")

    placement = read_tsv(paths["f2_placements"], ["assembly_ligand_placement_id", "source_ligand_instance_id"])
    placement = placement[placement["assembly_ligand_placement_id"].isin(set(f4["ligand_assembly_placement_id"]))]
    require_unique(placement, "assembly_ligand_placement_id", "Filter 2 retained placements")
    source_cols = ["source_ligand_instance_id", "resolved_ccd_id", "observed_heavy_atom_count",
                   "previous_terminal_route", "previous_reason_code", "terminal_route", "reason_code",
                   "is_suspicious", "contains_metal"]
    source = read_tsv(paths["f2_sources"], source_cols)
    source = source[source["source_ligand_instance_id"].isin(set(placement["source_ligand_instance_id"]))]
    require_unique(source, "source_ligand_instance_id", "Filter 2 source membership")

    pb_cols = ["pair_id", "internal_steric_clash", "posebusters_any_fail", "posebusters_status"]
    pb = pq.read_table(paths["pb"], columns=pb_cols).to_pandas()
    pb = pb[pb["pair_id"].isin(pair_ids)].copy()
    if pb["pair_id"].duplicated().any():
        pb = pb.groupby("pair_id", as_index=False).agg({
            "internal_steric_clash": "all", "posebusters_any_fail": "any", "posebusters_status": "first"
        })

    f3cols = ["pair_id", "pair_status", "pocket_missing_backbone_heavy_atom_count",
              "direct_binding_missing_sidechain_heavy_atom_count",
              "nonbinding_pocket_missing_sidechain_heavy_atom_count"]
    f3 = scan_partitioned_parquet(paths["f3_v2"], f3cols, pair_ids)
    require_unique(f3, "pair_id", "Filter 3 v2 quality evidence")

    chain_cols = ["pair_id", "source_chain_id"]
    chains = pq.read_table(paths["prov_chains"], columns=chain_cols).to_pandas()
    chains = chains[chains["pair_id"].isin(pair_ids)].drop_duplicates()
    if set(chains["pair_id"]) != pair_ids:
        raise RuntimeError("Protein provenance does not cover every formal pair")
    source_chain_ids = set(chains["source_chain_id"].dropna())

    umap = pq.read_table(paths["uniprot_map"], columns=["source_chain_id", "uniprot_accession"]).to_pandas()
    umap = umap[umap["source_chain_id"].isin(source_chain_ids)].dropna().drop_duplicates()
    pair_uniprot = chains.merge(umap, on="source_chain_id", how="left").groupby("pair_id")["uniprot_accession"].apply(normalize_list)

    umeta = pq.read_table(paths["uniprot_meta"], columns=["uniprot_accession", "organism_tax_id", "organism_name"]).to_pandas()
    umeta = umeta.drop_duplicates("uniprot_accession")
    chain_tax = umap.merge(umeta, on="uniprot_accession", how="left")
    src = pq.read_table(paths["source_chains"], columns=["source_chain_id", "source_tax_id", "source_organism_name"]).to_pandas()
    src = src[src["source_chain_id"].isin(source_chain_ids)].drop_duplicates("source_chain_id")
    chain_species = chains.merge(chain_tax, on="source_chain_id", how="left").merge(src, on="source_chain_id", how="left")
    chain_species["species_tax_id"] = chain_species["organism_tax_id"].fillna(chain_species["source_tax_id"])
    chain_species["species_name"] = chain_species["organism_name"].fillna(chain_species["source_organism_name"])
    pair_tax = chain_species.groupby("pair_id")["species_tax_id"].apply(normalize_list)
    pair_species_name = chain_species.groupby("pair_id")["species_name"].apply(normalize_list)

    fun = pq.read_table(paths["function"], columns=["source_chain_id", "resource", "annotation_json"]).to_pandas()
    fun = fun[(fun["source_chain_id"].isin(source_chain_ids)) & (fun["resource"] == "CATH")]
    def cath_id(x):
        try:
            return str(json.loads(x).get("CATH_ID", "")).strip() or None
        except Exception:
            return None
    fun["cath_id"] = fun["annotation_json"].map(cath_id)
    fun = fun.dropna(subset=["cath_id"])[["source_chain_id", "cath_id"]].drop_duplicates()
    pair_cath = chains.merge(fun, on="source_chain_id", how="left").groupby("pair_id")["cath_id"].apply(normalize_list)

    pocket = pq.read_table(paths["pocket"], columns=["pair_id", "pdb_residue_name"]).to_pandas()
    pocket = pocket[pocket["pair_id"].isin(pair_ids)].copy()
    pocket["pdb_residue_name"] = pocket["pdb_residue_name"].fillna("").astype(str).str.upper().str.strip()
    nonstd = pocket[~pocket["pdb_residue_name"].isin(STANDARD_PROTEIN_COMPONENTS)].copy()
    pair_nonstd = nonstd.groupby("pair_id")["pdb_residue_name"].apply(normalize_list)

    scaffold_cache, scaffold_unique_status = build_scaffolds(identity["normalized_ccd_isomeric_smiles"])
    identity["murcko_scaffold"] = identity["normalized_ccd_isomeric_smiles"].map(lambda x: scaffold_cache.get(str(x), (None, "MISSING_SMILES"))[0] if pd.notna(x) else None)
    identity["murcko_status"] = identity["normalized_ccd_isomeric_smiles"].map(lambda x: scaffold_cache.get(str(x), (None, "MISSING_SMILES"))[1] if pd.notna(x) else "MISSING_SMILES")

    entry = f4[["pair_id", "pdb_id", "assembly_id", "model_id", "ligand_assembly_placement_id",
                "component_id", "benchmark_filter3_terminal_status", "filter4_decision", "filter4_reason"]].copy()
    entry = entry.merge(identity, on="pair_id", how="left", validate="one_to_one")
    entry = entry.merge(p2, left_on="ligand_assembly_placement_id", right_on="assembly_ligand_placement_id", how="left", validate="many_to_one")
    entry = entry.merge(placement, on="assembly_ligand_placement_id", how="left", suffixes=("", "_f2"), validate="many_to_one")
    entry["source_ligand_instance_id"] = entry["source_ligand_instance_id"].fillna(entry["source_ligand_instance_id_f2"])
    entry = entry.drop(columns=["source_ligand_instance_id_f2"])
    entry = entry.merge(source, on="source_ligand_instance_id", how="left", suffixes=("", "_source"), validate="many_to_one")
    entry = entry.merge(pb, on="pair_id", how="left", validate="one_to_one")
    entry = entry.merge(f3, on="pair_id", how="left", validate="one_to_one")
    entry["uniprot_accessions"] = entry["pair_id"].map(pair_uniprot).map(lambda x: ";".join(x) if isinstance(x, list) else "")
    entry["species_tax_ids"] = entry["pair_id"].map(pair_tax).map(lambda x: ";".join(x) if isinstance(x, list) else "")
    entry["species_names"] = entry["pair_id"].map(pair_species_name).map(lambda x: ";".join(x) if isinstance(x, list) else "")
    entry["cath_ids"] = entry["pair_id"].map(pair_cath).map(lambda x: ";".join(x) if isinstance(x, list) else "")
    entry["nonstandard_pocket_residue_ids"] = entry["pair_id"].map(pair_nonstd).map(lambda x: ";".join(x) if isinstance(x, list) else "")
    entry["ion_ligand"] = pd.to_numeric(entry["observed_heavy_atom_count"], errors="coerce").eq(1)
    artifact_pattern = r"artifact|additive|solvent|simple_inorganic"
    entry["artifact_ligand"] = (entry["previous_terminal_route"].fillna("").str.contains(artifact_pattern, case=False, regex=True) |
                                entry["previous_reason_code"].fillna("").str.contains(artifact_pattern, case=False, regex=True))
    entry["covalent_ligand"] = ~entry["pair_status"].eq("FINAL_ORDINARY_NONCOVALENT_PAIR")
    entry["steric_overlap"] = entry["internal_steric_clash"].eq(False)
    entry["unresolved_ligand_atom_count"] = pd.to_numeric(entry["missing_heavy_atom_count"], errors="coerce")
    miss_cols = ["pocket_missing_backbone_heavy_atom_count", "direct_binding_missing_sidechain_heavy_atom_count",
                 "nonbinding_pocket_missing_sidechain_heavy_atom_count"]
    for c in miss_cols:
        entry[c] = pd.to_numeric(entry[c], errors="coerce")
    entry["unresolved_pocket_atom_count"] = entry[miss_cols].sum(axis=1, min_count=len(miss_cols))
    entry["unresolved_ligand_atoms"] = entry["unresolved_ligand_atom_count"].gt(0)
    entry["unresolved_pocket_atoms"] = entry["unresolved_pocket_atom_count"].gt(0)
    entry["nonstandard_pocket_residues"] = entry["nonstandard_pocket_residue_ids"].ne("")
    entry["affinity_annotation_available"] = False
    entry = entry.sort_values("pair_id").reset_index(drop=True)
    require_unique(entry, "pair_id", "Final entry properties")
    if len(entry) != 91860:
        raise RuntimeError("Final entry property table changed row count")

    pq.write_table(pa.Table.from_pandas(entry, preserve_index=False), out / "data" / "ours_entry_properties.parquet",
                   compression="zstd", use_dictionary=True, write_statistics=True)

    unique_uniprot = set(x for vals in pair_uniprot for x in vals)
    unique_tax = set(x for vals in pair_tax for x in vals)
    unique_cath = set(x for vals in pair_cath for x in vals)
    metrics = [
        metric_row("Dataset scope", "Total entries", len(entry), "CALCULATED", "Unique Filter 4 PASS pair_id records in the explicit final database-construction release.", "Filter 4 frozen PASS inventory", "91,860 formal entries"),
        metric_row("Dataset scope", "Unique PDB-CCD pairs", entry[["pdb_id", "resolved_ccd_id"]].drop_duplicates().shape[0], "CALCULATED", "Distinct (PDB ID, resolved CCD ID) combinations.", "Filter 4 + Filter 5 identity projection", "91,860 formal entries"),
        metric_row("Dataset scope", "Unique PDB IDs", entry["pdb_id"].nunique(), "CALCULATED", "Distinct lowercase PDB identifiers.", "Filter 4 frozen PASS inventory", "91,860 formal entries"),
        metric_row("Dataset scope", "Unique UniProt IDs", len(unique_uniprot), "CALCULATED", "Distinct UniProt accessions mapped to any participating receptor source chain.", "Frozen protein provenance annotation", "Participating receptor chains; partial mappings retained"),
        metric_row("Dataset scope", "Unique CATH IDs", len(unique_cath), "CALCULATED", "Distinct CATH_ID values mapped to any participating receptor source chain.", "Frozen SIFTS CATH annotations", "Participating receptor chains with CATH annotations"),
        metric_row("Dataset scope", "Unique species", len(unique_tax), "CALCULATED", "Distinct taxonomy IDs; UniProt organism taxon preferred, source-chain taxonomy used as fallback.", "Frozen UniProt metadata + source-chain provenance", "Participating receptor chains with a taxonomy ID"),
        metric_row("Dataset scope", "Affinity annotations", False, "CALCULATED", "Whether a formal integrated Kd/Ki/IC50/EC50 field exists in this frozen database release.", "Frozen release schema inventory", "Formal release", "No formal affinity field is integrated."),
        metric_row("Ligand diversity", "Unique CCD IDs", entry["resolved_ccd_id"].nunique(), "CALCULATED", "Distinct resolved CCD component identifiers.", "Filter 5 frozen identity projection", "91,860 formal entries"),
        metric_row("Ligand diversity", "Unique Murcko scaffolds", entry["murcko_scaffold"].dropna().nunique(), "CALCULATED", "Distinct non-empty canonical isomeric Bemis-Murcko scaffold SMILES; acyclic empty scaffolds excluded.", f"RDKit {rdBase.rdkitVersion}; normalized frozen CCD isomeric SMILES", "Entries with successfully parsed non-empty scaffolds"),
        metric_row("Ligand diversity", "Ion ligands", int(entry["ion_ligand"].sum()), "CALCULATED", "Formal entries whose mapped Filter 2 source ligand has exactly one observed heavy atom.", "Filter 2 frozen source classification", "91,860 formal entries"),
        metric_row("Ligand diversity", "Covalent ligands", int(entry["covalent_ligand"].sum()), "CALCULATED", "Formal entries not labeled FINAL_ORDINARY_NONCOVALENT_PAIR by frozen pair construction.", "Filter 3 v2 frozen pair status", "91,860 formal entries"),
        metric_row("Ligand diversity", "Artifact ligands", int(entry["artifact_ligand"].sum()), "CALCULATED", "Formal entries whose frozen pre-v4 Filter 2 route/reason explicitly contains artifact, additive, solvent, or simple_inorganic.", "Filter 2 frozen provenance fields", "91,860 formal entries", "No ad hoc CCD blacklist was introduced."),
        metric_row("Structure issues", "Missing bonds", None, "NOT_AVAILABLE", "Entry-level missing-bond predicate equivalent to the comparison-table definition.", "No formal frozen entry-level missing-bond field", "Not available", "Topology status is retained in the per-entry table but is not recast as missing bonds."),
        metric_row("Structure issues", "Steric overlaps", int(entry["steric_overlap"].sum()), "CALCULATED", "Entries failing the frozen PoseBusters internal_steric_clash check (False means overlap failure).", "Filter 3 strict PoseBusters evidence", "Formal entries with PoseBusters evidence"),
        metric_row("Structure issues", "Unresolved ligand atoms", int(entry["unresolved_ligand_atoms"].sum()), "CALCULATED", "Entries with Processing 2 missing_heavy_atom_count > 0.", "Processing 2 frozen formal-ready placements", "91,860 formal entries"),
        metric_row("Structure issues", "Unresolved pocket atoms", int(entry["unresolved_pocket_atoms"].sum()), "CALCULATED", "Entries with any frozen Filter 3 v2 missing backbone or side-chain heavy atom in the 6 A pocket.", "Filter 3 v2 frozen quality evidence", "91,860 formal entries"),
        metric_row("Structure issues", "Non-standard pocket residues", int(entry["nonstandard_pocket_residues"].sum()), "CALCULATED", "Entries containing a 6 A pocket residue outside the frozen protein template set (20 canonical residues plus MSE and SEC).", "Frozen pair-pocket residue inventory", "91,860 formal entries"),
    ]
    ours = pd.DataFrame(metrics)
    ours.to_csv(out / "data" / "ours_summary_stats.csv", index=False, lineterminator="\n")

    external_rows = []
    for dataset in ["PDBbind", "HiQBind", "BioLiP2", "PLINDER", "CROWN"]:
        for section, prop in PROPERTIES:
            external_rows.append({
                "dataset": dataset, "section": section, "property": prop, "value": "",
                "status": "NOT_AVAILABLE", "source": "Not available in the current repository",
                "retrieval_date": args.analysis_timestamp[:10],
                "notes": "Left NA rather than inferred from incomparable or unverified public counts.",
            })
    external = pd.DataFrame(external_rows)
    external.to_csv(out / "data" / "external_reference_stats.csv", index=False, lineterminator="\n")

    combined_rows = []
    ours_lookup = {(r.section, r.property): r.value for r in ours.itertuples()}
    ext_lookup = {(r.dataset, r.section, r.property): r.value for r in external.itertuples()}
    datasets = ["Ours", "PDBbind", "HiQBind", "BioLiP2", "PLINDER", "CROWN"]
    for section, prop in PROPERTIES:
        row = {"section": section, "property": prop}
        for dataset in datasets:
            row[dataset] = ours_lookup[(section, prop)] if dataset == "Ours" else ext_lookup[(dataset, section, prop)]
        combined_rows.append(row)
    pd.DataFrame(combined_rows).to_csv(out / "data" / "combined_comparison_stats.csv", index=False, lineterminator="\n")

    # Keep a non-empty field last so Git's whitespace check does not interpret
    # valid empty TSV notes as trailing whitespace.
    ours[["property", "definition", "source", "denominator", "notes", "status"]].to_csv(
        out / "qc" / "metric_definitions.tsv", sep="\t", index=False, lineterminator="\n")

    top_nonstd = (nonstd.groupby("pdb_residue_name").agg(
        residue_rows=("pair_id", "size"), entries=("pair_id", "nunique"))
        .reset_index().rename(columns={"pdb_residue_name": "component_id"})
        .sort_values(["entries", "residue_rows", "component_id"], ascending=[False, False, True]))
    top_nonstd.to_csv(out / "qc" / "top_nonstandard_pocket_residues.tsv", sep="\t", index=False, lineterminator="\n")
    entry["unresolved_ligand_atom_count"].value_counts(dropna=False).sort_index().rename_axis("missing_heavy_atom_count").reset_index(name="entries").to_csv(
        out / "qc" / "unresolved_ligand_atom_distribution.tsv", sep="\t", index=False, lineterminator="\n")
    entry["unresolved_pocket_atom_count"].value_counts(dropna=False).sort_index().rename_axis("missing_pocket_heavy_atom_count").reset_index(name="entries").to_csv(
        out / "qc" / "unresolved_pocket_atom_distribution.tsv", sep="\t", index=False, lineterminator="\n")

    missingness = []
    for col in ["resolved_ccd_id", "normalized_ccd_isomeric_smiles", "murcko_scaffold", "uniprot_accessions",
                "species_tax_ids", "cath_ids", "missing_heavy_atom_count", "internal_steric_clash",
                "unresolved_pocket_atom_count"]:
        s = entry[col]
        unavailable = s.isna() | (s.astype(str) == "")
        missingness.append({"field": col, "available_entries": int((~unavailable).sum()),
                            "missing_entries": int(unavailable.sum()), "total_entries": len(entry)})
    pd.DataFrame(missingness).to_csv(out / "qc" / "missingness_report.tsv", sep="\t", index=False, lineterminator="\n")

    provenance_rows = []
    for name, path in paths.items():
        if path.is_file():
            provenance_rows.append({"input": name, "path": str(path.relative_to(root)), "sha256": sha256(path), "size_bytes": path.stat().st_size})
        else:
            files = sorted(path.glob("**/*.parquet"))
            manifest_hash = hashlib.sha256("".join(f"{sha256(p)}  {p.relative_to(root)}\n" for p in files).encode()).hexdigest()
            provenance_rows.append({"input": name, "path": str(path.relative_to(root)), "sha256": manifest_hash,
                                    "size_bytes": sum(p.stat().st_size for p in files)})
    pd.DataFrame(provenance_rows).to_csv(out / "qc" / "source_provenance.tsv", sep="\t", index=False, lineterminator="\n")

    qc = [
        ("formal_rows", len(entry), 91860, len(entry) == 91860),
        ("unique_pair_id", entry["pair_id"].nunique(), len(entry), entry["pair_id"].nunique() == len(entry)),
        ("identity_join_missing", int(entry["resolved_ccd_id"].isna().sum()), 0, not entry["resolved_ccd_id"].isna().any()),
        ("processing2_join_missing", int(entry["missing_heavy_atom_count"].isna().sum()), 0, not entry["missing_heavy_atom_count"].isna().any()),
        ("filter2_join_missing", int(entry["observed_heavy_atom_count"].isna().sum()), 0, not entry["observed_heavy_atom_count"].isna().any()),
        ("posebusters_join_missing", int(entry["internal_steric_clash"].isna().sum()), 0, not entry["internal_steric_clash"].isna().any()),
        ("filter3_v2_join_missing", int(entry["pair_status"].isna().sum()), 0, not entry["pair_status"].isna().any()),
        ("protein_provenance_pair_missing", len(pair_ids - set(chains["pair_id"])), 0, pair_ids == set(chains["pair_id"])),
    ]
    pd.DataFrame(qc, columns=["check", "observed", "expected", "pass"]).to_csv(
        out / "qc" / "validation_report.tsv", sep="\t", index=False, lineterminator="\n")
    if not all(x[3] for x in qc):
        raise RuntimeError("One or more validation checks failed")

    report = f"""Authoritative input report
==========================
analysis_timestamp: {args.analysis_timestamp}
analysis_git_commit_before_changes: {args.git_commit}
data_root: {root}
release_version: {FINAL_F4_RUN}
database_role: {marker.get('database_role')}
formal_entry_definition: one unique pair_id in frozen Filter 4 PASS inventory
formal_input_path: {paths['f4_pass'].relative_to(root)}
formal_input_sha256: {sha256(paths['f4_pass'])}
formal_rows: {len(entry)}
unique_pair_id: {entry['pair_id'].nunique()}
selection_status: PASS

Authority rationale
-------------------
The selected run has an explicit frozen marker whose database_role is
FINAL_DATABASE_CONSTRUCTION_STAGE, its README defines Filter 4 as the final
database-construction stage, and the PASS inventory validates at 91,860 unique
pair records.  Later Filter 5 and Processing 4 products are benchmark
deduplication/docking-ready derivatives and are not database membership.

Known repository discrepancy
----------------------------
The repository-level manifests/frozen_runs.yaml still describes the older
158,226-case benchmark chain adopted on 2026-08-24.  It is not used as the
authority for this current database-member analysis.  No frozen file was
modified to resolve that documentation lag.
"""
    (out / "qc" / "authoritative_input_report.txt").write_text(report, encoding="utf-8", newline="\n")

    environment = {
        "analysis_timestamp": args.analysis_timestamp,
        "python": sys.version.replace("\n", " "), "platform": platform.platform(),
        "pandas": pd.__version__, "pyarrow": pa.__version__, "rdkit": rdBase.rdkitVersion,
    }
    (out / "qc" / "environment.json").write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    readme = """# Dataset comparison analysis

This directory builds an auditable CROWN-style property table for the current
formal database members.  The scientific inputs are read-only frozen releases;
the analysis does not change membership or scientific rules.

The authoritative population is the 91,860 unique `pair_id` records in the
frozen Filter 4 PASS inventory from
`20260826_filter3_118255_strict_posebusters_01`.  Filter 4 is explicitly marked
as `FINAL_DATABASE_CONSTRUCTION_STAGE`.  Filter 5 and Processing 4 are later
benchmark/docking-ready derivatives, not database membership.

Run `scripts/build_dataset_comparison.py`, then
`scripts/build_latex_table.py`.  Both accept explicit paths; no project path is
hard-coded.  External datasets are left `NA` because no confirmed reference
statistics were present in the repository.  See `qc/metric_definitions.tsv`,
`qc/source_provenance.tsv`, and `qc/authoritative_input_report.txt` before using
the table.

The entry-level Parquet, PDF, and logs are generated artifacts and are excluded
from Git under the repository data-management policy.  Compact scripts, CSV,
TeX, definitions, and QC records are tracked.
"""
    (out / "README.md").write_text(readme, encoding="utf-8", newline="\n")

    log = {
        "status": "PASS", "analysis_timestamp": args.analysis_timestamp,
        "formal_entries": len(entry), "unique_pair_ids": entry["pair_id"].nunique(),
        "scaffold_unique_smiles_status": scaffold_unique_status,
        "output_root": str(out),
    }
    (out / "logs" / "build_dataset_comparison.log").write_text(json.dumps(log, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
