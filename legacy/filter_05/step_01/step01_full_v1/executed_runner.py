#!/usr/bin/env python3
"""Filter 5 Step 1: exact ligand/receptor identity and candidate blocking."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from rdkit import Chem, rdBase, RDLogger

RDLogger.DisableLog("rdApp.*")


ROOT = Path("/root/autodl-tmp/benchmark_1.0")
FILTER4 = ROOT / "filter_04_crystal_packing_influence/step_05_final_crystal_packing_decision/runs/step05_full_v1/output/02_filter4_pass_pairs.tsv.gz"
P2RUN = ROOT / "processing_2_assembly_ready_structure_preparation/runs/20260810_full_01"
P3RUN = ROOT / "processing_03_direct_contact_qualification/runs/20260811_full_01"
SIFTS = ROOT / "filter_05_equivalent_redocking_case/inputs/sifts_snapshot"
STEPBASE = ROOT / "filter_05_equivalent_redocking_case/step_01_exact_identity_and_candidate_blocking/runs"
EXPECTED = 241_545


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 6, "mtime": 0})


def strip_cif_outer_quote(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    s = str(value).strip()
    # CIF semicolon-delimited text fields are also legal outer quoting.  The
    # frozen SQLite snapshot preserves those delimiters for a small number of
    # long descriptor values.
    if len(s) >= 3 and s.startswith(";") and s.endswith(";") and "\n" in s:
        return s[1:-1].strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in {"'", '"'}:
        return s[1:-1]
    return s


def canon(smiles: str, isomeric: bool) -> tuple[str, str]:
    if not smiles:
        return "", "MISSING_DESCRIPTOR"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "", "RDKIT_PARSE_FAILURE"
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=isomeric), "OK"
    except Exception as exc:
        return "", f"RDKIT_EXCEPTION:{type(exc).__name__}"


def read_sifts(path: Path, needed_pdb: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    df.columns = [c.lower() for c in df.columns]
    df["pdb"] = df["pdb"].str.lower()
    df = df[df["pdb"].isin(needed_pdb)].copy()
    for c in ["res_beg", "res_end", "sp_beg", "sp_end"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


def parse_chain_instance(x: str) -> dict[str, str]:
    p = str(x).split("|")
    return {
        "pdb_id": p[0].lower() if len(p) > 0 else "",
        "assembly_id": p[1] if len(p) > 1 else "",
        "model_id": p[2] if len(p) > 2 else "",
        "label_asym_id": p[3] if len(p) > 3 else "",
        "operator_id": "|".join(p[4:]) if len(p) > 4 else "",
    }


def scan_chain_metadata(required: set[str]) -> pd.DataFrame:
    """Project only five low-cardinality columns from the frozen source-atom dataset."""
    path = P2RUN / "output/prepared_receptor_source_atoms"
    dataset = ds.dataset(path, format="parquet")
    cols = ["pdb_id", "model_id", "entity_id", "label_asym_id", "auth_asym_id"]
    pieces: list[pd.DataFrame] = []
    for i, frag in enumerate(dataset.get_fragments()):
        tab = frag.to_table(columns=cols)
        if tab.num_rows == 0:
            continue
        d = tab.to_pandas().drop_duplicates(cols)
        # Source data includes the benchmark universe; restricting here saves memory.
        d["source_chain_key"] = d["pdb_id"].astype(str) + "|" + d["model_id"].astype(str) + "|" + d["label_asym_id"].astype(str)
        d = d[d["source_chain_key"].isin(required)]
        if not d.empty:
            pieces.append(d)
        if i % 500 == 0:
            print(f"chain metadata fragments {i}", flush=True)
    if not pieces:
        return pd.DataFrame(columns=cols + ["source_chain_key"])
    out = pd.concat(pieces, ignore_index=True).drop_duplicates()
    return out


def frozen_observed_sequence_hashes(missing: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """Conservative Level-2 fallback from frozen source residue identities."""
    if missing.empty:
        return {}
    missing_keys = set(missing["source_chain_key"])
    missing_pdb = set(missing["pdb_id"])
    path = P2RUN / "output/prepared_receptor_source_atoms"
    dataset = ds.dataset(path, format="parquet")
    cols = ["pdb_id", "model_id", "label_asym_id", "label_seq_id", "auth_seq_id", "insertion_code", "label_comp_id"]
    residues: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)
    for i, frag in enumerate(dataset.get_fragments()):
        tab = frag.to_table(columns=cols)
        if tab.num_rows == 0:
            continue
        d = tab.to_pandas().drop_duplicates(cols)
        d = d[d["pdb_id"].isin(missing_pdb)]
        if d.empty:
            continue
        d["source_chain_key"] = d["pdb_id"].astype(str) + "|" + d["model_id"].astype(str) + "|" + d["label_asym_id"].astype(str)
        for r in d[d["source_chain_key"].isin(missing_keys)].itertuples(index=False):
            residues[r.source_chain_key].add((str(r.label_seq_id), str(r.auth_seq_id), str(r.insertion_code), str(r.label_comp_id)))
        if i % 500 == 0:
            print(f"fallback sequence fragments {i}", flush=True)

    def skey(t: tuple[str, str, str, str]):
        try:
            return (0, int(float(t[0])), t[1], t[2])
        except Exception:
            return (1, t[0], t[1], t[2])

    out = {}
    for key, vals in residues.items():
        ordered = sorted(vals, key=skey)
        seq = ";".join("|".join(v) for v in ordered)
        out[key] = (hashlib.sha256(seq.encode()).hexdigest(), seq)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", default="step01_full_v1")
    args = ap.parse_args()
    run = STEPBASE / args.run_id
    out = run / "output"
    valdir = run / "validation"
    logs = run / "logs"
    for d in [out, valdir, logs]:
        d.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    t0 = time.time()

    # Formal membership and frozen Processing 3 receptor-chain interface.
    f4 = pd.read_csv(FILTER4, sep="\t", dtype={"pdb_id": str})
    prov = ds.dataset(P3RUN / "output/provisional_pairs", format="parquet").to_table(
        columns=["pair_id", "ligand_assembly_placement_id", "pdb_id", "component_id", "receptor_chain_instance_ids"]
    ).to_pandas()
    pairs = f4.merge(prov, on=["pair_id", "pdb_id"], how="left", validate="one_to_one", indicator=True)
    if len(pairs) != EXPECTED or pairs["pair_id"].duplicated().any() or (pairs["_merge"] != "both").any():
        raise RuntimeError("formal input membership gate failed")
    pairs = pairs.drop(columns="_merge")

    # Phase B: normalize every component in the frozen active CCD snapshot once.
    db = sqlite3.connect(P2RUN / "input/ccd_active_snapshot.sqlite")
    comp = pd.read_sql_query(
        "SELECT component_id,resolved_ccd_id,name,ccd_type,descriptor_smiles_canonical,canonical_smiles_from_graph "
        "FROM (SELECT component_id,component_id AS resolved_ccd_id,name,ccd_type,descriptor_smiles_canonical,canonical_smiles_from_graph FROM components)", db
    )
    atom_stereo = pd.read_sql_query(
        "SELECT component_id, SUM(CASE WHEN upper(stereo_config) IN ('R','S') THEN 1 ELSE 0 END) AS ccd_atom_rs_count "
        "FROM atoms GROUP BY component_id", db
    )
    db.close()
    comp = comp.merge(atom_stereo, on="component_id", how="left")
    comp["ccd_atom_rs_count"] = comp["ccd_atom_rs_count"].fillna(0).astype(int)
    comp["ccd_descriptor_smiles"] = comp["descriptor_smiles_canonical"].map(strip_cif_outer_quote)
    conn, iso, status, graph_conn, graph_status = [], [], [], [], []
    for raw, graph in zip(comp["ccd_descriptor_smiles"], comp["canonical_smiles_from_graph"]):
        c, s1 = canon(raw, False)
        i, s2 = canon(raw, True)
        g, sg = canon(strip_cif_outer_quote(graph), False)
        conn.append(c); iso.append(i); status.append(s1 if s1 != "OK" else s2)
        graph_conn.append(g); graph_status.append(sg)
    comp["normalized_ccd_connectivity_smiles"] = conn
    comp["normalized_ccd_isomeric_smiles"] = iso
    comp["normalization_status"] = status
    comp["normalized_graph_connectivity_smiles"] = graph_conn
    comp["graph_normalization_status"] = graph_status
    comp["descriptor_graph_connectivity_status"] = [
        "MATCH" if a and a == b else ("DESCRIPTOR_GRAPH_CONNECTIVITY_MISMATCH" if a and b else "NOT_COMPARABLE")
        for a, b in zip(conn, graph_conn)
    ]
    comp["raw_contains_atom_stereo_token"] = comp["ccd_descriptor_smiles"].str.contains("@", regex=False)
    comp["raw_contains_bond_stereo_token"] = comp["ccd_descriptor_smiles"].str.contains(r"[/\\]", regex=True)
    comp["normalized_contains_atom_stereo_token"] = comp["normalized_ccd_isomeric_smiles"].str.contains("@", regex=False)
    comp["normalized_contains_bond_stereo_token"] = comp["normalized_ccd_isomeric_smiles"].str.contains(r"[/\\]", regex=True)
    valid_iso = sorted(x for x in comp.loc[comp.normalization_status == "OK", "normalized_ccd_isomeric_smiles"].unique() if x)
    lex = {s: f"LEX{i:08d}" for i, s in enumerate(valid_iso, 1)}
    comp["ligand_exact_id"] = comp["normalized_ccd_isomeric_smiles"].map(lex).fillna("")
    comp["chemistry_review_flag"] = (
        (comp["normalization_status"] != "OK") |
        (comp["descriptor_graph_connectivity_status"] == "DESCRIPTOR_GRAPH_CONNECTIVITY_MISMATCH") |
        ((comp["ccd_atom_rs_count"] > 0) & ~comp["normalized_contains_atom_stereo_token"])
    )
    ligand_cols = [
        "resolved_ccd_id", "ccd_descriptor_smiles", "canonical_smiles_from_graph",
        "normalized_ccd_connectivity_smiles", "normalized_ccd_isomeric_smiles", "ligand_exact_id",
        "normalization_status", "descriptor_graph_connectivity_status", "ccd_atom_rs_count",
        "raw_contains_atom_stereo_token", "raw_contains_bond_stereo_token",
        "normalized_contains_atom_stereo_token", "normalized_contains_bond_stereo_token", "chemistry_review_flag"
    ]
    write_tsv_gz(comp[ligand_cols], out / "01_ligand_exact_identity_map.tsv.gz")

    # Expand only the frozen receptor chain instances attached to each pair.
    chain_rows = []
    for r in pairs[["pair_id", "receptor_chain_instance_ids"]].itertuples(index=False):
        for cid in str(r.receptor_chain_instance_ids).split(";"):
            d = parse_chain_instance(cid)
            d.update(pair_id=r.pair_id, receptor_chain_instance_id=cid,
                     source_chain_key=f"{d['pdb_id']}|{d['model_id']}|{d['label_asym_id']}")
            chain_rows.append(d)
    chains = pd.DataFrame(chain_rows)
    required_source = set(chains["source_chain_key"])
    metadata = scan_chain_metadata(required_source)
    # A frozen source chain should have a unique author-chain mapping.
    meta_summary = metadata.groupby("source_chain_key", as_index=False).agg(
        entity_id=("entity_id", lambda x: ";".join(sorted(set(map(str, x))))),
        auth_asym_id=("auth_asym_id", lambda x: ";".join(sorted(set(map(str, x)))))
    )
    chains = chains.merge(meta_summary, on="source_chain_key", how="left")
    chains["auth_asym_id"] = chains["auth_asym_id"].fillna("")

    needed_pdb = set(pairs["pdb_id"].str.lower())
    s_all = read_sifts(SIFTS / "pdb_chain_uniprot.tsv.gz", needed_pdb)
    s_obs = read_sifts(SIFTS / "uniprot_segments_observed.tsv.gz", needed_pdb)
    seg_by = defaultdict(list)
    for r in s_all.itertuples(index=False):
        seg_by[(r.pdb, str(r.chain))].append(r)

    fallback_needed = []
    chain_summary = []
    chain_segments = []
    for r in chains.itertuples(index=False):
        auths = [x for x in str(r.auth_asym_id).split(";") if x]
        segs = []
        for auth in auths:
            segs.extend(seg_by.get((r.pdb_id, auth), []))
        # Deduplicate identical bulk SIFTS segment rows.
        uniq = {}
        for s in segs:
            k = (s.sp_primary, s.res_beg, s.res_end, s.sp_beg, s.sp_end)
            uniq[k] = s
        segs = list(uniq.values())
        accs = sorted(set(str(s.sp_primary) for s in segs if str(s.sp_primary)))
        if len(auths) != 1:
            status, method, ident = "SOURCE_CHAIN_MAPPING_AMBIGUOUS", "UNRESOLVED", ""
        elif len(accs) == 1:
            status, method, ident = "SIFTS_UNIPROT_MAPPED", "SIFTS_UNIPROT", accs[0]
        elif len(accs) > 1:
            ordered = sorted(segs, key=lambda s: (int(s.res_beg) if pd.notna(s.res_beg) else 10**9, str(s.sp_primary)))
            ordered_acc = []
            for s in ordered:
                if not ordered_acc or ordered_acc[-1] != str(s.sp_primary):
                    ordered_acc.append(str(s.sp_primary))
            status, method, ident = "RECEPTOR_IDENTITY_COMPLEX_MAPPING", "SIFTS_SEGMENT_AWARE", "CHIMERA[" + ">".join(ordered_acc) + "]"
        else:
            status, method, ident = "SIFTS_MAPPING_MISSING", "FALLBACK_PENDING", ""
            fallback_needed.append(r._asdict())
        for j, s in enumerate(sorted(segs, key=lambda z: (str(z.sp_primary), z.res_beg if pd.notna(z.res_beg) else 10**9)), 1):
            chain_segments.append({
                "pair_id": r.pair_id, "receptor_chain_instance_id": r.receptor_chain_instance_id,
                "pdb_id": r.pdb_id, "label_asym_id": r.label_asym_id, "auth_asym_id": r.auth_asym_id,
                "operator_id": r.operator_id, "mapping_segment_id": f"SEG{j:04d}",
                "uniprot_accession": s.sp_primary, "pdb_label_seq_start": s.res_beg,
                "pdb_label_seq_end": s.res_end, "pdb_auth_seq_start": s.pdb_beg,
                "pdb_auth_seq_end": s.pdb_end, "uniprot_start": s.sp_beg, "uniprot_end": s.sp_end,
                "segment_source": "pdb_chain_uniprot.tsv.gz"
            })
        chain_summary.append({**r._asdict(), "mapping_status": status,
                              "chain_receptor_identity": ident, "receptor_identity_method": method,
                              "uniprot_accessions": ";".join(accs), "mapping_segment_count": len(segs),
                              "frozen_sequence_sha256": ""})
    chain_summary = pd.DataFrame(chain_summary)

    if fallback_needed:
        miss = pd.DataFrame(fallback_needed)
        seqhash = frozen_observed_sequence_hashes(miss)
        for idx, r in chain_summary[chain_summary.mapping_status == "SIFTS_MAPPING_MISSING"].iterrows():
            hseq = seqhash.get(r.source_chain_key)
            if hseq:
                h, _seq = hseq
                chain_summary.loc[idx, ["mapping_status", "chain_receptor_identity", "receptor_identity_method", "frozen_sequence_sha256"]] = [
                    "SIFTS_MISSING_SEQUENCE_FALLBACK", f"SEQ:{h}", "FROZEN_OBSERVED_RESIDUE_SEQUENCE_SHA256", h
                ]
            else:
                chain_summary.loc[idx, ["mapping_status", "receptor_identity_method"]] = ["RECEPTOR_IDENTITY_UNRESOLVED", "UNRESOLVED"]

    chain_out_cols = ["pair_id", "receptor_chain_instance_id", "pdb_id", "label_asym_id", "auth_asym_id", "operator_id",
                      "entity_id", "mapping_status", "uniprot_accessions", "mapping_segment_count",
                      "chain_receptor_identity", "receptor_identity_method", "frozen_sequence_sha256"]
    write_tsv_gz(chain_summary[chain_out_cols], out / "02a_receptor_chain_identity_map.tsv.gz")
    segdf = pd.DataFrame(chain_segments)
    write_tsv_gz(segdf, out / "02b_receptor_sifts_segments.tsv.gz")

    # Pair receptor identity is a sorted multiset: duplicate homomer identities are retained.
    pair_rec = []
    for pid, g in chain_summary.groupby("pair_id", sort=False):
        ids = list(g["chain_receptor_identity"])
        unresolved = any(not x for x in ids)
        complex_map = any(x == "RECEPTOR_IDENTITY_COMPLEX_MAPPING" for x in g["mapping_status"])
        methods = sorted(set(g["receptor_identity_method"]))
        pair_rec.append({
            "pair_id": pid,
            "receptor_identity_key": "|".join(sorted(ids)) if not unresolved else "",
            "receptor_identity_status": "UNRESOLVED" if unresolved else ("COMPLEX_MAPPING" if complex_map else "RESOLVED"),
            "receptor_identity_method": ";".join(methods),
            "receptor_chain_count": len(ids),
            "receptor_identity_review_flag": unresolved or complex_map,
        })
    pair_rec = pd.DataFrame(pair_rec)
    write_tsv_gz(pair_rec, out / "02_receptor_identity_map.tsv.gz")

    # Phase D: exact chemistry + exact biological receptor target blocking.
    ligand_join = comp[["resolved_ccd_id", "ligand_exact_id", "normalized_ccd_isomeric_smiles", "normalization_status", "chemistry_review_flag"]]
    inv = pairs.merge(ligand_join, left_on="component_id", right_on="resolved_ccd_id", how="left", validate="many_to_one")
    inv = inv.merge(pair_rec, on="pair_id", how="left", validate="one_to_one")
    inv["blocking_eligible"] = inv["ligand_exact_id"].fillna("").ne("") & inv["receptor_identity_key"].fillna("").ne("")
    inv["candidate_block_key"] = ""
    ok = inv["blocking_eligible"]
    inv.loc[ok, "candidate_block_key"] = inv.loc[ok, "ligand_exact_id"] + "|" + inv.loc[ok, "receptor_identity_key"]
    block_keys = sorted(inv.loc[ok, "candidate_block_key"].unique())
    bid = {k: f"F5B{i:08d}" for i, k in enumerate(block_keys, 1)}
    inv["candidate_block_id"] = inv["candidate_block_key"].map(bid).fillna("")
    sizes = inv.loc[ok].groupby("candidate_block_id")["pair_id"].size()
    inv["candidate_block_size"] = inv["candidate_block_id"].map(sizes).fillna(0).astype(int)
    inv["step1_pair_status"] = "F5_STEP1_REVIEW"
    inv.loc[ok & (inv["candidate_block_size"] == 1), "step1_pair_status"] = "F5_STEP1_SINGLETON"
    inv.loc[ok & (inv["candidate_block_size"] >= 2), "step1_pair_status"] = "F5_STEP1_SITE_AUDIT_CANDIDATE"
    inv["step1_review_flag"] = (~inv["blocking_eligible"]) | inv["chemistry_review_flag"].fillna(True) | inv["receptor_identity_review_flag"].fillna(True)
    inv["step1_reason"] = ""
    inv.loc[~inv["ligand_exact_id"].fillna("").ne(""), "step1_reason"] += "LIGAND_ID_UNRESOLVED;"
    inv.loc[~inv["receptor_identity_key"].fillna("").ne(""), "step1_reason"] += "RECEPTOR_ID_UNRESOLVED;"
    inv.loc[inv["chemistry_review_flag"].fillna(True), "step1_reason"] += "CHEMISTRY_REVIEW;"
    inv.loc[inv["receptor_identity_review_flag"].fillna(True), "step1_reason"] += "RECEPTOR_MAPPING_REVIEW;"
    inv["step1_reason"] = inv["step1_reason"].str.rstrip(";").replace("", "DETERMINISTIC_BLOCKING")
    inv_cols = ["pair_id", "pdb_id", "assembly_id", "model_id", "ligand_assembly_placement_id", "component_id",
                "resolved_ccd_id", "ligand_exact_id", "normalized_ccd_isomeric_smiles", "receptor_chain_instance_ids",
                "receptor_identity_key", "receptor_identity_status", "receptor_identity_method", "candidate_block_id",
                "candidate_block_size", "step1_pair_status", "step1_review_flag", "step1_reason"]
    write_tsv_gz(inv[inv_cols], out / "03_filter5_step1_pair_inventory.tsv.gz")
    write_tsv_gz(inv.loc[inv.step1_pair_status == "F5_STEP1_SITE_AUDIT_CANDIDATE", inv_cols], out / "04_filter5_step1_candidate_pairs.tsv.gz")
    write_tsv_gz(inv.loc[inv.step1_pair_status == "F5_STEP1_SINGLETON", inv_cols], out / "05_filter5_step1_singletons.tsv.gz")
    write_tsv_gz(inv.loc[inv.step1_review_flag, inv_cols], out / "06_filter5_step1_review.tsv.gz")

    # Summary and validation.
    size_dist = sizes.value_counts().sort_index()
    summary_rows = [
        ("input_pairs", len(inv)), ("input_duplicate_pair_id", int(inv.pair_id.duplicated().sum())),
        ("silent_drop", EXPECTED - len(inv)), ("active_components", len(comp)),
        ("unique_exact_ligand_identities", comp.loc[comp.normalization_status == "OK", "ligand_exact_id"].nunique()),
        ("unique_receptor_identity_keys", inv.loc[inv.receptor_identity_key.ne(""), "receptor_identity_key"].nunique()),
        ("singleton_pairs", int((inv.step1_pair_status == "F5_STEP1_SINGLETON").sum())),
        ("non_singleton_pairs", int((inv.step1_pair_status == "F5_STEP1_SITE_AUDIT_CANDIDATE").sum())),
        ("primary_review_pairs", int((inv.step1_pair_status == "F5_STEP1_REVIEW").sum())),
        ("auxiliary_review_pairs", int(inv.step1_review_flag.sum())),
        ("candidate_blocks_non_singleton", int((sizes >= 2).sum())),
        ("largest_block_size", int(sizes.max() if len(sizes) else 0)),
        ("chemistry_review_components", int(comp.chemistry_review_flag.sum())),
        ("receptor_unresolved_pairs", int((inv.receptor_identity_status == "UNRESOLVED").sum())),
        ("receptor_complex_mapping_pairs", int((inv.receptor_identity_status == "COMPLEX_MAPPING").sum())),
    ]
    for n, count in size_dist.items():
        summary_rows.append((f"block_size_{n}_count", int(count)))
    summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
    summary.to_csv(out / "07_filter5_step1_summary.tsv", sep="\t", index=False)
    largest = inv[inv.candidate_block_size >= 2].groupby(
        ["candidate_block_id", "ligand_exact_id", "receptor_identity_key"], as_index=False
    ).agg(block_size=("pair_id", "size"), pdb_count=("pdb_id", "nunique"), component_ids=("component_id", lambda x: ";".join(sorted(set(x)))))
    largest = largest.sort_values(["block_size", "candidate_block_id"], ascending=[False, True])
    largest.to_csv(out / "08_filter5_step1_largest_blocks.tsv", sep="\t", index=False)

    chem_qc = {
        "active_components": len(comp),
        "descriptor_nonempty": int(comp.ccd_descriptor_smiles.ne("").sum()),
        "descriptor_missing": int(comp.ccd_descriptor_smiles.eq("").sum()),
        "rdkit_parse_success": int((comp.normalization_status == "OK").sum()),
        "rdkit_parse_failure": int((comp.normalization_status != "OK").sum()),
        "unique_normalized_connectivity_smiles": int(comp.normalized_ccd_connectivity_smiles.replace("", pd.NA).nunique()),
        "unique_normalized_isomeric_smiles": int(comp.normalized_ccd_isomeric_smiles.replace("", pd.NA).nunique()),
        "raw_at_token": int(comp.raw_contains_atom_stereo_token.sum()),
        "raw_slash_token": int(comp.raw_contains_bond_stereo_token.sum()),
        "normalized_at_token": int(comp.normalized_contains_atom_stereo_token.sum()),
        "normalized_slash_token": int(comp.normalized_contains_bond_stereo_token.sum()),
        "ccd_atom_rs_components": int((comp.ccd_atom_rs_count > 0).sum()),
        "ccd_atom_rs_normalized_preserves_at": int(((comp.ccd_atom_rs_count > 0) & comp.normalized_contains_atom_stereo_token).sum()),
        "ccd_atom_rs_stereo_lost": int(((comp.ccd_atom_rs_count > 0) & ~comp.normalized_contains_atom_stereo_token).sum()),
        "descriptor_graph_connectivity_mismatch": int((comp.descriptor_graph_connectivity_status == "DESCRIPTOR_GRAPH_CONNECTIVITY_MISMATCH").sum()),
    }
    validation = {
        "status": "PASS",
        "checks": {
            "input_pair_count_241545": len(inv) == EXPECTED,
            "pair_id_unique": not inv.pair_id.duplicated().any(),
            "silent_drop_zero": len(inv) == EXPECTED,
            "primary_partition_reconstructs_input": int((inv.step1_pair_status.isin(["F5_STEP1_SINGLETON", "F5_STEP1_SITE_AUDIT_CANDIDATE", "F5_STEP1_REVIEW"])).sum()) == EXPECTED,
            "candidate_interface_only_non_singletons": bool((inv.loc[inv.step1_pair_status == "F5_STEP1_SITE_AUDIT_CANDIDATE", "candidate_block_size"] >= 2).all()),
            "ligand_map_component_count_45250": len(comp) == 45250,
        },
        "chemistry_qc": chem_qc,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.time() - t0, 3),
        "versions": {"python": sys.version, "rdkit": rdBase.rdkitVersion, "pandas": pd.__version__, "pyarrow": pa.__version__},
    }
    if not all(validation["checks"].values()):
        validation["status"] = "FAIL"
    (valdir / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8")
    report = f"""# Filter 5 Step 1 report\n\nStatus: **{validation['status']}**\n\nThis step performs candidate generation only. It does not delete, collapse, or select representatives. RDKit parses and canonically serializes frozen CCD descriptor chemistry; it does not infer chemistry from coordinates and no standardization, tautomerization, neutralization, reionization, or protonation normalization is applied.\n\n- Formal Filter 4 PASS input: {len(inv):,}\n- Exact ligand identities: {summary.set_index('metric').loc['unique_exact_ligand_identities','value']:,}\n- Receptor identity keys: {summary.set_index('metric').loc['unique_receptor_identity_keys','value']:,}\n- Singleton pairs: {(inv.step1_pair_status == 'F5_STEP1_SINGLETON').sum():,}\n- Non-singleton candidate pairs: {(inv.step1_pair_status == 'F5_STEP1_SITE_AUDIT_CANDIDATE').sum():,}\n- Primary unresolved/review pairs: {(inv.step1_pair_status == 'F5_STEP1_REVIEW').sum():,}\n- Non-singleton blocks: {(sizes >= 2).sum():,}\n- Largest block: {int(sizes.max() if len(sizes) else 0):,}\n- Silent drops: 0\n\nSIFTS snapshot: 2026-08-09, PDB 32.26, UniProt 2026.03. Chain identity uses `pdb_chain_uniprot.tsv.gz`; residue-level mapping for Step 2 uses `uniprot_segments_observed.tsv.gz`. Chimera accessions remain segment-aware and are flagged for review.\n"""
    (out / "09_filter5_step1_report.md").write_text(report, encoding="utf-8")

    # Provenance, schemas, manifest, checksums, and immutable marker only on PASS.
    provenance = {
        "formal_membership_input": str(FILTER4), "formal_membership_sha256": sha256(FILTER4),
        "processing2_run": str(P2RUN), "processing3_run": str(P3RUN),
        "sifts_retrieval_date": "2026-08-17",
        "sifts_sources": [
            {"url": "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz", "sha256": sha256(SIFTS / "pdb_chain_uniprot.tsv.gz")},
            {"url": "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_segments_observed.tsv.gz", "sha256": sha256(SIFTS / "uniprot_segments_observed.tsv.gz")},
        ],
    }
    (run / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    schema = {}
    for p in sorted(out.iterdir()):
        if p.name.endswith(".tsv.gz"):
            schema[p.name] = {c: str(t) for c, t in pd.read_csv(p, sep="\t", nrows=10).dtypes.items()}
        elif p.name.endswith(".tsv"):
            schema[p.name] = {c: str(t) for c, t in pd.read_csv(p, sep="\t", nrows=10).dtypes.items()}
    (run / "output_schema.json").write_text(json.dumps(schema, indent=2), encoding="utf-8")
    manifest_rows = []
    for p in sorted(out.iterdir()):
        if p.is_file():
            rows = ""
            if p.name.endswith(".tsv.gz"):
                with gzip.open(p, "rt") as fh: rows = max(sum(1 for _ in fh) - 1, 0)
            elif p.name.endswith(".tsv"):
                with p.open() as fh: rows = max(sum(1 for _ in fh) - 1, 0)
            manifest_rows.append((p.name, p.stat().st_size, rows, sha256(p)))
    pd.DataFrame(manifest_rows, columns=["file", "bytes", "data_rows", "sha256"]).to_csv(run / "output_manifest.tsv", sep="\t", index=False)
    checksum_files = [p for p in run.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "_FROZEN.json"}]
    with (run / "SHA256SUMS").open("w") as fh:
        for p in sorted(checksum_files):
            fh.write(f"{sha256(p)}  {p.relative_to(run)}\n")
    if validation["status"] == "PASS":
        (run / "_FROZEN.json").write_text(json.dumps({"status": "FROZEN", "frozen_utc": datetime.now(timezone.utc).isoformat(), "validation": "validation/validation.json"}, indent=2), encoding="utf-8")
    print(json.dumps(validation, indent=2), flush=True)
    return 0 if validation["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
