#!/usr/bin/env python3
"""Build six-dataset CATH H-level and ligand-taxonomy comparison annotations.

This script is read-only with respect to all frozen scientific datasets.  It
writes only under analysis/dataset_comparison/{harmonized_cath_v1,
comparison_ligand_taxonomy_v1} and updates compact comparison CSVs.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
import tarfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

DATASETS = ["PDBbind", "HiQBind", "BioLiP2", "PLINDER", "CROWN", "Ours"]
PREPARED_FILES = {
    "PDBbind": "pdbbind_properties_quality.parquet",
    "HiQBind": "hiqbind_properties_quality.parquet",
    "BioLiP2": "biolip2_properties_quality.parquet",
    "PLINDER": "plinder_quality.parquet",
    "CROWN": "crown_quality.parquet",
    "Ours": "ours_properties_harmonized_quality.parquet",
}
VALID_H = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
PLINDER_NON_ION_ELEMENTS = {"C", "H", "*", "N", "O", "P", "S", "Se", "Te"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def parse_artifact_list(path: Path) -> set[str]:
    return {line.strip().upper() for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")}


def classify_smiles(smiles: object) -> dict:
    empty = {
        "taxonomy_graph_status": "UNAVAILABLE", "expected_heavy_atom_count": None,
        "carbon_atom_count": None, "unique_heavy_elements": None,
        "monoatomic_ion_flag": False, "simple_inorganic_flag": False,
        "monoatomic_reason": "GRAPH_UNAVAILABLE", "simple_inorganic_reason": "GRAPH_UNAVAILABLE",
    }
    if not isinstance(smiles, str) or not smiles.strip():
        return empty
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {**empty, "taxonomy_graph_status": "PARSE_FAILED"}
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    elements = sorted({a.GetSymbol() for a in atoms})
    heavy = len(atoms)
    carbon = sum(a.GetAtomicNum() == 6 for a in atoms)
    mono = heavy == 1 and atoms[0].GetSymbol() not in PLINDER_NON_ION_ELEMENTS
    inorganic = heavy >= 2 and carbon == 0
    return {
        "taxonomy_graph_status": "CLASSIFIED", "expected_heavy_atom_count": heavy,
        "carbon_atom_count": carbon, "unique_heavy_elements": ";".join(elements),
        "monoatomic_ion_flag": mono, "simple_inorganic_flag": inorganic,
        "monoatomic_reason": "SINGLE_HEAVY_ATOM_NON_C_H_DUMMY_N_O_P_S_SE_TE" if mono else "NOT_MONOATOMIC_ION",
        "simple_inorganic_reason": "MULTI_HEAVY_ATOM_CARBON_FREE" if inorganic else "NOT_SIMPLE_INORGANIC",
    }


def build_taxonomy(repo: Path, artifact_list: Path) -> None:
    root = repo / "analysis/dataset_comparison/comparison_ligand_taxonomy_v1"
    out, qc = root / "output", root / "qc"
    out.mkdir(parents=True, exist_ok=True); qc.mkdir(parents=True, exist_ok=True)
    prepared = repo / "analysis/dataset_property_comparison/prepared"
    artifacts = parse_artifact_list(artifact_list)
    all_frames, summary, missing, overlaps, top_rows, review_rows = [], [], [], [], [], []
    for dataset in DATASETS:
        p = prepared / PREPARED_FILES[dataset]
        d = pd.read_parquet(p, columns=["source_entry_id", "pdb_id", "ccd_id", "canonical_isomeric_smiles"])
        if len(d) != d.source_entry_id.nunique():
            raise RuntimeError(f"{dataset}: source_entry_id is not unique")
        cache = {s: classify_smiles(s) for s in set(d.canonical_isomeric_smiles.dropna())}
        missing_result = classify_smiles(None)
        ann = pd.DataFrame([cache.get(s, missing_result) for s in d.canonical_isomeric_smiles])
        d = pd.concat([d.reset_index(drop=True), ann], axis=1)
        d["dataset"] = dataset
        d["normalized_ccd_id"] = d.ccd_id.fillna("").astype(str).str.strip().str.upper()
        d["shared_artifact_list_flag"] = d.normalized_ccd_id.isin(artifacts)
        d["shared_artifact_reason"] = d.shared_artifact_list_flag.map(
            {True: "DIRECT_NORMALIZED_CCD_MATCH_PLINDER_V0.2.0_LIST", False: "NO_DIRECT_CCD_MATCH"})
        d["taxonomy_version"] = "comparison_ligand_taxonomy_v1"
        all_frames.append(d)
        n = len(d)
        classified = int((d.taxonomy_graph_status == "CLASSIFIED").sum())
        counts = {
            "monoatomic_ion_entries": int(d.monoatomic_ion_flag.sum()),
            "simple_inorganic_entries": int(d.simple_inorganic_flag.sum()),
            "shared_artifact_list_entries": int(d.shared_artifact_list_flag.sum()),
        }
        for metric, value in counts.items():
            summary.append({"dataset": dataset, "formal_N": n, "metric": metric, "entries": value,
                            "percent_of_formal_N": value / n * 100, "taxonomy_version": "comparison_ligand_taxonomy_v1"})
        missing.append({"dataset": dataset, "formal_N": n, "graph_classified_N": classified,
                        "graph_unavailable_N": n - classified, "graph_coverage_percent": classified / n * 100,
                        "ccd_available_N": int(d.normalized_ccd_id.ne("").sum()),
                        "artifact_list_size": len(artifacts)})
        combo = d.groupby(["monoatomic_ion_flag", "simple_inorganic_flag", "shared_artifact_list_flag"], dropna=False).size()
        for flags, count in combo.items():
            overlaps.append({"dataset": dataset, "monoatomic_ion_flag": flags[0],
                             "simple_inorganic_flag": flags[1], "shared_artifact_list_flag": flags[2],
                             "entries": int(count)})
        for metric, flag in [("monoatomic_ion_entries", "monoatomic_ion_flag"),
                             ("simple_inorganic_entries", "simple_inorganic_flag"),
                             ("shared_artifact_list_entries", "shared_artifact_list_flag")]:
            sub = d[d[flag]].groupby("normalized_ccd_id", dropna=False).agg(
                entries=("source_entry_id", "size"), unique_pdb=("pdb_id", "nunique")).reset_index()
            sub = sub.sort_values(["entries", "normalized_ccd_id"], ascending=[False, True]).head(50)
            for row in sub.itertuples(index=False):
                top_rows.append({"dataset": dataset, "metric": metric, "ccd_id": row.normalized_ccd_id,
                                 "entries": int(row.entries), "unique_pdb": int(row.unique_pdb)})
            candidates = d[d[flag]].copy()
            candidates["sample_key"] = candidates.source_entry_id.map(
                lambda x: hashlib.sha256(f"taxonomy_v1|{dataset}|{metric}|{x}".encode()).hexdigest())
            for row in candidates.sort_values("sample_key").head(25).itertuples(index=False):
                review_rows.append({"dataset": dataset, "metric": metric, "source_entry_id": row.source_entry_id,
                                    "pdb_id": row.pdb_id, "ccd_id": row.normalized_ccd_id,
                                    "canonical_isomeric_smiles": row.canonical_isomeric_smiles,
                                    "heavy_atoms": row.expected_heavy_atom_count,
                                    "carbon_atoms": row.carbon_atom_count,
                                    "elements": row.unique_heavy_elements,
                                    "review_status": "DETERMINISTIC_AUDIT_SAMPLE_NOT_MEMBERSHIP_REVIEW"})
    full = pd.concat(all_frames, ignore_index=True)
    full.to_parquet(out / "six_dataset_entry_taxonomy.parquet", index=False)
    write_tsv(out / "harmonized_ligand_taxonomy_summary.tsv", summary)
    write_tsv(out / "taxonomy_overlap_summary.tsv", overlaps)
    write_tsv(out / "taxonomy_top_ccd.tsv", top_rows)
    write_tsv(qc / "manual_review_sample.tsv", review_rows)
    write_tsv(qc / "taxonomy_coverage_qc.tsv", missing)
    if any(x["formal_N"] != x["graph_classified_N"] + x["graph_unavailable_N"] for x in missing):
        raise RuntimeError("taxonomy coverage does not close")
    if (full.monoatomic_ion_flag & full.simple_inorganic_flag).any():
        raise RuntimeError("monoatomic and simple-inorganic definitions are not mutually exclusive")


def load_cath_domain_map(path: Path):
    domain_to_h, chain_to_domains = {}, defaultdict(list)
    with path.open(encoding="ascii") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.split()
            if len(x) < 5:
                continue
            domain, h = x[0], ".".join(x[1:5])
            if VALID_H.fullmatch(h):
                domain_to_h[domain] = h
                chain_to_domains[(domain[:4].lower(), domain[4])].append(domain)
    return domain_to_h, chain_to_domains


def load_sifts(path: Path, domain_to_h: dict):
    d = pd.read_csv(path, sep="\t", comment="#", compression="gzip", dtype=str)
    d.columns = [x.strip() for x in d.columns]
    d["PDB"] = d.PDB.str.lower(); d["CHAIN"] = d.CHAIN.astype(str)
    d["cath_h_id"] = d.CATH_ID.map(domain_to_h)
    return d


def parse_h_values(value) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple)):
        values = value
    else:
        values = re.split(r"[;,|\s]+", str(value).strip())
    return sorted({str(v) for v in values if VALID_H.fullmatch(str(v))})


def chains_from_binding(value) -> set[str]:
    if not isinstance(value, str):
        return set()
    return set(re.findall(r"(?:^|\s)([^:\s]+):", value))


def pdbbind_protein_chains(external: Path) -> dict[str, set[str]]:
    result = defaultdict(set)
    for archive_name in ["PDBbind_v2020_refined.tar.gz", "PDBbind_v2020_other_PL.tar.gz"]:
        with tarfile.open(external / "pdbbind_v2020" / archive_name, "r:gz") as tf:
            for member in tf:
                if not member.isfile() or not member.name.endswith("_protein.pdb"):
                    continue
                pdb = Path(member.name).parent.name.lower()
                handle = tf.extractfile(member)
                for raw in handle:
                    if raw.startswith(b"ATOM  ") and len(raw) >= 22:
                        result[pdb].add(raw[21:22].decode("ascii", "ignore").strip() or "_")
    return result


def build_cath(repo: Path, external: Path, cath_list: Path, sifts_path: Path) -> None:
    root = repo / "analysis/dataset_comparison/harmonized_cath_v1"
    out, qc = root / "output", root / "qc"
    out.mkdir(parents=True, exist_ok=True); qc.mkdir(parents=True, exist_ok=True)
    prepared = repo / "analysis/dataset_property_comparison/prepared"
    domain_to_h, _ = load_cath_domain_map(cath_list)
    sifts = load_sifts(sifts_path, domain_to_h)
    chain_h = sifts.dropna(subset=["cath_h_id"]).groupby(["PDB", "CHAIN"])["cath_h_id"].agg(lambda x: sorted(set(x))).to_dict()
    rows = []

    def add(dataset, entry, ids, scope_status):
        ids = sorted({x for x in ids if VALID_H.fullmatch(str(x))})
        rows.append({"dataset": dataset, "source_entry_id": str(entry), "cath_h_ids": ids,
                     "cath_h_count": len(ids), "annotation_status": "ANNOTATED" if ids else scope_status})

    # CROWN released per-entry CATH H arrays.
    c = pq.read_table(external / "crown_202606" / "CROWN_metadata.parquet", columns=["basename", "cath_ids"]).to_pydict()
    for entry, ids in zip(c["basename"], c["cath_ids"]): add("CROWN", entry, parse_h_values(ids), "NO_VALID_CATH_H")

    # PLINDER released pocket-level H classification.
    p = pd.read_parquet(external / "plinder_2024-06_v2" / "annotation_table.parquet",
                        columns=["ligand_id", "ligand_is_proper", "system_pocket_CATH"])
    p = p[p.ligand_is_proper]
    for x in p.itertuples(index=False): add("PLINDER", x.ligand_id, parse_h_values(x.system_pocket_CATH), "NO_VALID_CATH_H")

    # Ours formal participating receptor chains.
    ours = pd.read_parquet(prepared / PREPARED_FILES["Ours"], columns=["source_entry_id"])
    ours_ch = pd.read_parquet(
        "/home/linx/data/youcq/autodl-tmp/benchmark_1.0/protein_provenance_annotation/runs/20260826_filter4_pass_01/output/pair_receptor_chains.parquet",
        columns=["pair_id", "pdb_id", "auth_asym_id"])
    ours_map = defaultdict(set)
    for x in ours_ch.itertuples(index=False): ours_map[x.pair_id].update(chain_h.get((str(x.pdb_id).lower(), str(x.auth_asym_id)), []))
    for entry in ours.source_entry_id: add("Ours", entry, ours_map.get(entry, []), "NO_MAPPED_FORMAL_RECEPTOR_CATH_H")

    # BioLiP2 formal binding-site chains from exact raw records.
    b = pd.read_parquet(prepared / PREPARED_FILES["BioLiP2"], columns=["source_entry_id", "pdb_id"])
    bio_chains = {}
    base = pd.read_csv(external / "biolip2_20260626" / "PL_annotation_base_before_20260102.csv", dtype=str)
    for x in base.itertuples(index=False):
        key = "BASE:" + str(getattr(x, "Ligand_file"))
        bio_chains[key] = chains_from_binding(getattr(x, "Binding_residue"))
    for f in sorted((external / "biolip2_20260626/weekly").glob("Q-BioLiP-*.csv")):
        w = pd.read_csv(f, dtype=str)
        for _, x in w.iterrows(): bio_chains["WEEKLY:" + str(x["Ligand Detail"])] = chains_from_binding(x.get("Binding Site"))
    for x in b.itertuples(index=False):
        ids = set()
        for chain in bio_chains.get(x.source_entry_id, set()): ids.update(chain_h.get((str(x.pdb_id).lower(), str(chain)), []))
        add("BioLiP2", x.source_entry_id, ids, "NO_MAPPED_BINDING_SITE_CHAIN_CATH_H")

    # HiQBind receptor chains matched through released protein UniProt accessions and SIFTS.
    h = pd.read_parquet(prepared / PREPARED_FILES["HiQBind"], columns=["source_entry_id", "pdb_id", "source_row_id"])
    raw_h = pd.read_csv(external / "hiqbind_v3/hiqbind_sm_metadata.csv", dtype=str)
    sifts_by_pdb = {k: g for k, g in sifts.groupby("PDB")}
    for x in h.itertuples(index=False):
        val = raw_h.iloc[int(x.source_row_id)]["Protein UniProtID"]
        accessions = {z.strip() for z in str(val).split(",") if z.strip() and z.strip().lower() != "nan"}
        ids = set(); g = sifts_by_pdb.get(str(x.pdb_id).lower())
        if g is not None:
            ids.update(g.loc[g.SP_PRIMARY.isin(accessions), "cath_h_id"].dropna())
        add("HiQBind", x.source_entry_id, ids, "NO_UNIPROT_MATCHED_RECEPTOR_CATH_H")

    # PDBbind released receptor-only protein files; stream archives without extraction.
    pb = pd.read_parquet(prepared / PREPARED_FILES["PDBbind"], columns=["source_entry_id", "pdb_id"])
    pb_chains = pdbbind_protein_chains(external)
    for x in pb.itertuples(index=False):
        ids = set()
        for chain in pb_chains.get(str(x.pdb_id).lower(), set()): ids.update(chain_h.get((str(x.pdb_id).lower(), chain), []))
        add("PDBbind", x.source_entry_id, ids, "NO_MAPPED_RELEASED_RECEPTOR_CHAIN_CATH_H")

    result = pd.DataFrame(rows)
    result.to_parquet(out / "six_dataset_entry_cath_h.parquet", index=False)
    summary = []
    for dataset in DATASETS:
        d = result[result.dataset == dataset]
        ids = {h for values in d.cath_h_ids for h in values}
        expected = len(pd.read_parquet(prepared / PREPARED_FILES[dataset], columns=["source_entry_id"]))
        summary.append({"dataset": dataset, "formal_N": expected, "annotation_rows": len(d),
                        "entries_with_cath_h": int((d.cath_h_count > 0).sum()),
                        "entries_without_cath_h": int((d.cath_h_count == 0).sum()),
                        "entry_coverage_percent": float((d.cath_h_count > 0).mean() * 100),
                        "unique_valid_four_level_cath_h": len(ids), "null_counted_as_id": False})
    write_tsv(out / "harmonized_cath_summary.tsv", summary)
    if next(x for x in summary if x["dataset"] == "CROWN")["unique_valid_four_level_cath_h"] != 2040:
        raise RuntimeError("CROWN harmonized CATH calibration did not reproduce 2040")
    if any(x["formal_N"] != x["annotation_rows"] for x in summary):
        raise RuntimeError("CATH population does not close")
    write_tsv(qc / "cath_population_qc.tsv", [{**x, "status": "PASS"} for x in summary])


def update_tables(repo: Path) -> None:
    cath = pd.read_csv(repo / "analysis/dataset_comparison/harmonized_cath_v1/output/harmonized_cath_summary.tsv", sep="\t")
    tax = pd.read_csv(repo / "analysis/dataset_comparison/comparison_ligand_taxonomy_v1/output/harmonized_ligand_taxonomy_summary.tsv", sep="\t")
    cath_map = dict(zip(cath.dataset, cath.unique_valid_four_level_cath_h))
    wide = tax.pivot(index="metric", columns="dataset", values="entries")
    combined_path = repo / "analysis/dataset_comparison/data/combined_comparison_stats.csv"
    combined = pd.read_csv(combined_path)
    mask = (combined.section == "Dataset scope") & (combined.property == "Unique CATH IDs")
    for ds in DATASETS: combined.loc[mask, ds] = cath_map[ds]
    metric_labels = {
        "monoatomic_ion_entries": "Monoatomic ion entries (harmonized v1)",
        "simple_inorganic_entries": "Simple inorganic entries (harmonized v1)",
        "shared_artifact_list_entries": "Shared artifact-list entries (harmonized v1)",
    }
    combined = combined[~combined.property.isin(metric_labels.values())]
    for metric, label in metric_labels.items():
        row = {"section": "Ligand taxonomy v1", "property": label}
        row.update({ds: int(wide.loc[metric, ds]) for ds in DATASETS})
        combined = pd.concat([combined, pd.DataFrame([row])], ignore_index=True)
    combined.to_csv(combined_path, index=False)

    ours_path = repo / "analysis/dataset_comparison/data/ours_summary_stats.csv"
    ours = pd.read_csv(ours_path)
    mask = (ours.section == "Dataset scope") & (ours.property == "Unique CATH IDs")
    ours.loc[mask, ["value", "status", "definition", "source", "denominator", "notes"]] = [
        cath_map["Ours"], "CALCULATED", "Unique valid four-level CATH H-level classifications; null/unannotated values are not IDs.",
        "Formal receptor chains -> SIFTS domain instance -> official CATH-Plus v4.4.0 mapping",
        "91,860 formal Filter 4 PASS entries", "Harmonized against CROWN metadata expected value 2,040."]
    ours = ours[~ours.property.isin(metric_labels.values())]
    for metric, label in metric_labels.items():
        val = int(wide.loc[metric, "Ours"])
        ours = pd.concat([ours, pd.DataFrame([{"section": "Ligand taxonomy v1", "property": label,
            "value": val, "status": "CALCULATED", "definition": metric,
            "source": "comparison_ligand_taxonomy_v1", "denominator": "91,860 formal entries",
            "notes": "Retrospective harmonized comparison annotation; database membership unchanged."}])], ignore_index=True)
    ours.to_csv(ours_path, index=False)

    reported = [
        {"reporting_source": "CROWN website", "dataset": "CROWN", "metric": "Unique CATH IDs", "value": 2041,
         "status": "CROWN_REPORTED_ONLY", "harmonized_value": 2040,
         "note": "Website value retained for provenance; null/unannotated is excluded from harmonized IDs."},
        {"reporting_source": "CROWN website", "dataset": "PLINDER", "metric": "Ion ligands", "value": 22728,
         "status": "CROWN_REPORTED_ONLY", "harmonized_value": "NA",
         "note": "Not reproduced by trial filters; replaced in formal comparison by monoatomic_ion_entries."},
        {"reporting_source": "CROWN website", "dataset": "PLINDER", "metric": "Artifact ligands", "value": 18626,
         "status": "CROWN_REPORTED_ONLY", "harmonized_value": "NA",
         "note": "Not reproduced by trial filters; replaced by explicit harmonized taxonomy metrics."},
    ]
    write_tsv(repo / "analysis/dataset_comparison/data/crown_reported_statistics_supplementary.tsv", reported)


def write_docs(repo: Path, artifact_list: Path, cath_list: Path) -> None:
    tax = repo / "analysis/dataset_comparison/comparison_ligand_taxonomy_v1"
    cath = repo / "analysis/dataset_comparison/harmonized_cath_v1"
    (tax / "README.md").write_text("""# comparison_ligand_taxonomy_v1

Entry-level retrospective annotations for the six Mode B populations. No database membership changes.

- `monoatomic_ion_entries`: exactly one expected heavy atom whose element is outside C/H/dummy/N/O/P/S/Se/Te, matching the executable PLINDER v0.2.0 single-atom test. Formal charge is not required.
- `simple_inorganic_entries`: at least two expected heavy atoms and zero carbon atoms. This category excludes monoatomic entries.
- `shared_artifact_list_entries`: direct normalized CCD-code match to the frozen PLINDER v0.2.0 curated artifact list. It is an independent overlapping flag; synonym expansion is deliberately not inferred.

The first two labels are mutually exclusive. Artifact-list membership can overlap either and all intersections are reported. Graph-unavailable entries remain in the denominator and are reported in QC.
""", encoding="utf-8")
    (cath / "README.md").write_text("""# harmonized_cath_v1

Formal metric: unique valid four-level CATH H-level classifications. Null or unannotated values never count as an ID. CROWN metadata calibrates exactly at 2,040; the website-reported 2,041 is retained only in the reported-statistics supplementary table.

Receptor scope sources: CROWN released CATH arrays; PLINDER released pocket CATH; Ours formal participating receptor chains; BioLiP binding-site chains; HiQBind released protein UniProt accessions mapped through SIFTS; PDBbind released receptor-only protein files. SIFTS domain instances are converted with official CATH-Plus v4.4.0.
""", encoding="utf-8")
    manifests = [
        {"source": "PLINDER v0.2.0 artifact list", "path": str(artifact_list), "sha256": sha256(artifact_list)},
        {"source": "CATH-Plus v4.4.0 domain list", "path": str(cath_list), "sha256": sha256(cath_list)},
    ]
    write_tsv(tax / "references/source_manifest.tsv", manifests)
    write_tsv(cath / "references/source_manifest.tsv", manifests[1:])


def finalize_hashes(root: Path) -> None:
    files = [p for p in root.rglob("*") if p.is_file() and p.name != "hashes.sha256" and "__pycache__" not in p.parts]
    lines = [f"{sha256(p)}  {p.relative_to(root).as_posix()}" for p in sorted(files)]
    (root / "hashes.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--external", type=Path, required=True)
    args = ap.parse_args()
    artifact = args.repo / "analysis/dataset_comparison/comparable_metric_audit/external/plinder/artifacts_badlist.csv"
    cath_list = args.repo / "analysis/dataset_comparison/comparable_metric_audit/external/cath/cath-domain-list-v4_4_0.txt"
    sifts = Path("/home/linx/data/youcq/autodl-tmp/benchmark_1.0/protein_provenance_annotation/references/sifts_20260826/flatfiles/pdb_chain_cath_uniprot.tsv.gz")
    build_taxonomy(args.repo, artifact)
    build_cath(args.repo, args.external, cath_list, sifts)
    update_tables(args.repo)
    write_docs(args.repo, artifact, cath_list)
    finalize_hashes(args.repo / "analysis/dataset_comparison/comparison_ligand_taxonomy_v1")
    finalize_hashes(args.repo / "analysis/dataset_comparison/harmonized_cath_v1")


if __name__ == "__main__":
    main()
