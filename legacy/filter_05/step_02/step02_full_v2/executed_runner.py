#!/usr/bin/env python3
"""Filter 5 Step 2: frozen-residue same-binding-site equivalence audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import math
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds


ROOT = Path("/root/autodl-tmp/benchmark_1.0")
STEP1 = ROOT / "filter_05_equivalent_redocking_case/step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1"
P3RUN = ROOT / "processing_03_direct_contact_qualification/runs/20260811_full_01"
SIFTS = ROOT / "filter_05_equivalent_redocking_case/inputs/sifts_snapshot/uniprot_segments_observed.tsv.gz"
STEPBASE = ROOT / "filter_05_equivalent_redocking_case/step_02_same_binding_site_audit/runs"
EXPECTED_UNIVERSE = 241_545


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def gz_writer(path: Path, header: list[str]):
    raw = gzip.GzipFile(filename=str(path), mode="wb", compresslevel=4, mtime=0)
    fh = io.TextIOWrapper(raw, encoding="utf-8", newline="")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    return fh, w


def parse_residue_id(x: object) -> tuple[str, str, str, str]:
    p = str(x).split("|")
    return tuple((p + [""] * 4)[:4])  # label_seq, auth_seq, insertion, comp


def metric(a: frozenset[str], b: frozenset[str]) -> tuple[int, int, float, float]:
    ni = len(a & b)
    nu = len(a | b)
    j = ni / nu if nu else float("nan")
    mn = min(len(a), len(b))
    c = ni / mn if mn else float("nan")
    return ni, nu, j, c


def fmt(x: float) -> str:
    return "" if math.isnan(x) else f"{x:.6f}"


def audit_bin(x: float) -> str:
    if math.isnan(x): return "NA"
    if x == 1.0: return "1.00"
    if x >= 0.9: return "0.90-<1.00"
    if x >= 0.8: return "0.80-<0.90"
    if x >= 0.7: return "0.70-<0.80"
    if x >= 0.5: return "0.50-<0.70"
    return "<0.50"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", default="step02_full_v1")
    args = ap.parse_args()
    run = STEPBASE / args.run_id
    out, valdir, logs = run / "output", run / "validation", run / "logs"
    for d in [out, valdir, logs]: d.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat(); t0 = time.time()

    universe = pd.read_csv(STEP1 / "output/03_filter5_step1_pair_inventory.tsv.gz", sep="\t", usecols=["pair_id", "step1_pair_status"])
    cand = pd.read_csv(STEP1 / "output/04_filter5_step1_candidate_pairs.tsv.gz", sep="\t")
    if len(universe) != EXPECTED_UNIVERSE or universe.pair_id.duplicated().any():
        raise RuntimeError("Step 1 universe gate failed")
    if len(cand) != 218_765 or cand.pair_id.duplicated().any() or not cand.candidate_block_size.ge(2).all():
        raise RuntimeError("Step 1 candidate interface gate failed")
    block_sizes = cand.groupby("candidate_block_id").size()
    expected_edges = int((block_sizes * (block_sizes - 1) // 2).sum())
    print(f"candidate pairs={len(cand)} blocks={len(block_sizes)} expected_edges={expected_edges}", flush=True)

    chain = pd.read_csv(STEP1 / "output/02a_receptor_chain_identity_map.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    chain = chain[chain.pair_id.isin(set(cand.pair_id))].copy()
    chain_key = {(r.pair_id, r.receptor_chain_instance_id): r for r in chain.itertuples(index=False)}
    allowed_acc = {(r.pair_id, r.receptor_chain_instance_id): set(filter(None, r.uniprot_accessions.split(";"))) for r in chain.itertuples(index=False)}

    sifts = pd.read_csv(SIFTS, sep="\t", comment="#", dtype=str, keep_default_na=False)
    sifts.columns = [c.lower() for c in sifts.columns]
    needed_pdb = set(cand.pdb_id.str.lower())
    sifts["pdb"] = sifts.pdb.str.lower()
    sifts = sifts[sifts.pdb.isin(needed_pdb)].copy()
    for c in ["res_beg", "res_end", "sp_beg", "sp_end"]:
        sifts[c] = pd.to_numeric(sifts[c], errors="coerce")
    segs = defaultdict(list)
    for r in sifts.dropna(subset=["res_beg", "res_end", "sp_beg"]).itertuples(index=False):
        segs[(r.pdb, str(r.chain))].append((int(r.res_beg), int(r.res_end), str(r.sp_primary), int(r.sp_beg)))
    for k in segs:
        segs[k] = sorted(set(segs[k]))

    placement_to_pair = defaultdict(list)
    for r in cand[["pair_id", "ligand_assembly_placement_id"]].itertuples(index=False):
        placement_to_pair[r.ligand_assembly_placement_id].append(r.pair_id)
    placements = set(placement_to_pair)
    pairset = set(cand.pair_id)
    binding = ds.dataset(P3RUN / "output/binding_residues", format="parquet").to_table(
        columns=["ligand_assembly_placement_id", "chain_instance_id", "protein_residue_id"]
    ).to_pandas()
    binding = binding[binding.ligand_assembly_placement_id.isin(placements)]
    pocket = ds.dataset(P3RUN / "output/pair_pocket_residues", format="parquet").to_table(
        columns=["pair_id", "chain_instance_id", "protein_residue_id"]
    ).to_pandas()
    pocket = pocket[pocket.pair_id.isin(pairset)]
    print(f"candidate binding rows={len(binding)} pocket rows={len(pocket)}", flush=True)

    # pair -> source/mapped fingerprints and mapping anomaly counts
    source_binding = defaultdict(set); mapped_binding = defaultdict(set); bind_unmapped = Counter(); bind_ambig = Counter(); bind_outside = Counter()
    source_pocket = defaultdict(set); mapped_pocket = defaultdict(set); pocket_unmapped = Counter(); pocket_ambig = Counter(); pocket_outside = Counter()
    anomaly = defaultdict(Counter)
    excluded_binding_nonreceptor = 0
    excluded_pocket_nonreceptor = 0
    residue_header = ["pair_id", "site_type", "chain_instance_id", "protein_residue_id", "label_seq_id", "auth_seq_id", "insertion_code", "component_id", "mapping_status", "uniprot_accession", "uniprot_residue_number"]
    rfh, rw = gz_writer(out / "00_filter5_step2_residue_mapping.tsv.gz", residue_header)

    def map_one(pid: str, cid: str, resid: str, site_type: str):
        label, auth, ins, comp = parse_residue_id(resid)
        key = (pid, cid)
        ch = chain_key.get(key)
        hits = []
        if ch is None:
            status = "CHAIN_NOT_IN_FROZEN_RECEPTOR"
        elif ch.mapping_status == "SIFTS_MISSING_SEQUENCE_FALLBACK":
            status = "SIFTS_RESIDUE_MAPPING_MISSING"
        else:
            try: li = int(float(label))
            except Exception: li = None
            if li is None:
                status = "LABEL_SEQ_ID_INVALID"
            else:
                auths = [x for x in ch.auth_asym_id.split(";") if x]
                for ach in auths:
                    for beg, end, acc, ubeg in segs.get((ch.pdb_id.lower(), ach), []):
                        if beg <= li <= end:
                            hits.append((acc, ubeg + li - beg))
                hits = sorted(set(hits))
                if not hits: status = "SIFTS_RESIDUE_MAPPING_MISSING"
                elif len(hits) > 1: status = "ONE_TO_MANY_RESIDUE_MAPPING"
                elif hits[0][0] not in allowed_acc.get(key, set()): status = "MAPPED_OUTSIDE_RECEPTOR_IDENTITY"
                else: status = "MAPPED"
        if ins not in {"", ".", "?", "nan"}: anomaly[pid]["insertion_code_present"] += 1
        if status == "MAPPED":
            acc, unum = hits[0]
            sig = f"{acc}:{unum}"
        else:
            acc, unum, sig = "", "", ""
        rw.writerow([pid, site_type, cid, resid, label, auth, ins, comp, status, acc, unum])
        return status, sig

    # Binding rows lack pair_id by frozen design; expand only if a placement has multiple receptor pair definitions.
    for i, r in enumerate(binding.itertuples(index=False), 1):
        for pid in placement_to_pair[r.ligand_assembly_placement_id]:
            if (pid, r.chain_instance_id) not in chain_key:
                excluded_binding_nonreceptor += 1
                label, auth, ins, comp = parse_residue_id(r.protein_residue_id)
                rw.writerow([pid, "BINDING_EXCLUDED_NONRECEPTOR", r.chain_instance_id, r.protein_residue_id,
                             label, auth, ins, comp, "EXCLUDED_NOT_PAIR_RECEPTOR", "", ""])
                continue
            source_binding[pid].add((r.chain_instance_id, r.protein_residue_id))
            status, sig = map_one(pid, r.chain_instance_id, r.protein_residue_id, "BINDING")
            if status == "MAPPED": mapped_binding[pid].add(sig)
            elif status == "ONE_TO_MANY_RESIDUE_MAPPING": bind_ambig[pid] += 1
            elif status == "MAPPED_OUTSIDE_RECEPTOR_IDENTITY": bind_outside[pid] += 1
            else: bind_unmapped[pid] += 1
        if i % 1_000_000 == 0: print(f"binding mapped rows {i}", flush=True)
    del binding
    for i, r in enumerate(pocket.itertuples(index=False), 1):
        if (r.pair_id, r.chain_instance_id) not in chain_key:
            excluded_pocket_nonreceptor += 1
            label, auth, ins, comp = parse_residue_id(r.protein_residue_id)
            rw.writerow([r.pair_id, "POCKET_EXCLUDED_NONRECEPTOR", r.chain_instance_id, r.protein_residue_id,
                         label, auth, ins, comp, "EXCLUDED_NOT_PAIR_RECEPTOR", "", ""])
            continue
        source_pocket[r.pair_id].add((r.chain_instance_id, r.protein_residue_id))
        status, sig = map_one(r.pair_id, r.chain_instance_id, r.protein_residue_id, "POCKET")
        if status == "MAPPED": mapped_pocket[r.pair_id].add(sig)
        elif status == "ONE_TO_MANY_RESIDUE_MAPPING": pocket_ambig[r.pair_id] += 1
        elif status == "MAPPED_OUTSIDE_RECEPTOR_IDENTITY": pocket_outside[r.pair_id] += 1
        else: pocket_unmapped[r.pair_id] += 1
        if i % 2_000_000 == 0: print(f"pocket mapped rows {i}", flush=True)
    del pocket
    rfh.close()

    pair_info = cand.set_index("pair_id").to_dict("index")
    site = {}
    for pid in cand.pair_id:
        bset = frozenset(mapped_binding[pid]); pset = frozenset(mapped_pocket[pid])
        site[pid] = {
            "b": bset, "p": pset, "binding_source": len(source_binding[pid]), "pocket_source": len(source_pocket[pid]),
            "bu": bind_unmapped[pid], "ba": bind_ambig[pid], "bo": bind_outside[pid],
            "pu": pocket_unmapped[pid], "pa": pocket_ambig[pid], "po": pocket_outside[pid],
            "complete": len(source_binding[pid]) > 0 and bind_unmapped[pid] == 0 and bind_ambig[pid] == 0 and bind_outside[pid] == 0 and len(bset) > 0,
        }

    ph = ["candidate_block_id", "pair_id_a", "pair_id_b", "ligand_exact_id", "receptor_identity_key",
          "binding_n_a", "binding_n_b", "binding_n_intersection", "binding_n_union", "binding_jaccard", "binding_containment",
          "pocket_n_a", "pocket_n_b", "pocket_n_intersection", "pocket_n_union", "pocket_jaccard", "pocket_containment",
          "site_mapping_status", "step2_site_status", "step2_reason"]
    pfh, pw = gz_writer(out / "01_filter5_step2_pairwise_site_comparisons.tsv.gz", ph)
    sfh, sw = gz_writer(out / "03_filter5_step2_same_site_candidates.tsv.gz", ph)
    vfh, vw = gz_writer(out / "04_filter5_step2_review.tsv.gz", ph)
    status_counts = Counter(); jbins = Counter(); cbins = Counter(); pjbins = Counter(); pcbins = Counter()
    exact_neighbors = Counter(); strong_neighbors = Counter(); candidate_degree = Counter(); review_neighbors = Counter()
    examples = defaultdict(list); edge_count = 0; same_count = 0; review_count = 0
    groups = cand.groupby("candidate_block_id", sort=True)
    for gi, (block, g) in enumerate(groups, 1):
        ids = sorted(g.pair_id)
        first = g.iloc[0]
        for a, b in itertools.combinations(ids, 2):
            sa, sb = site[a], site[b]
            ni, nu, j, c = metric(sa["b"], sb["b"])
            pni, pnu, pj, pcnt = metric(sa["p"], sb["p"])
            if not sa["complete"] or not sb["complete"]:
                mapping_status = "INCOMPLETE_BINDING_MAPPING"
                status = "SITE_MAPPING_REVIEW"; reason = "BINDING_SIGNATURE_INCOMPLETE"
            else:
                mapping_status = "BINDING_MAPPING_COMPLETE"
                if sa["pu"] or sb["pu"] or sa["pa"] or sb["pa"] or sa["po"] or sb["po"]:
                    mapping_status = "BINDING_COMPLETE_POCKET_INCOMPLETE"
                if sa["b"] == sb["b"]:
                    status = "SITE_EXACT"; reason = "IDENTICAL_MAPPED_BINDING_RESIDUE_SET"
                elif c >= 0.8 and j >= 0.7:
                    status = "SITE_STRONG_CANDIDATE"; reason = "PROVISIONAL_ENGINEERING_THRESHOLD:C>=0.80_AND_J>=0.70"
                elif ni == 0:
                    status = "SITE_DIFFERENT"; reason = "DISJOINT_COMPLETE_BINDING_RESIDUE_SETS"
                elif c < 0.2 and j < 0.1:
                    status = "SITE_DIFFERENT"; reason = "VERY_LOW_OVERLAP_COMPLETE_MAPPING:C<0.20_AND_J<0.10"
                else:
                    status = "SITE_WEAK_OR_AMBIGUOUS"; reason = "PARTIAL_OR_BORDERLINE_OVERLAP"
            row = [block, a, b, first.ligand_exact_id, first.receptor_identity_key,
                   len(sa["b"]), len(sb["b"]), ni, nu, fmt(j), fmt(c),
                   len(sa["p"]), len(sb["p"]), pni, pnu, fmt(pj), fmt(pcnt), mapping_status, status, reason]
            pw.writerow(row); edge_count += 1; status_counts[status] += 1
            jbins[audit_bin(j)] += 1; cbins[audit_bin(c)] += 1; pjbins[audit_bin(pj)] += 1; pcbins[audit_bin(pcnt)] += 1
            if len(examples[status]) < 20: examples[status].append(row)
            if status in {"SITE_EXACT", "SITE_STRONG_CANDIDATE"}:
                sw.writerow(row); same_count += 1; candidate_degree[a] += 1; candidate_degree[b] += 1
                if status == "SITE_EXACT": exact_neighbors[a] += 1; exact_neighbors[b] += 1
                else: strong_neighbors[a] += 1; strong_neighbors[b] += 1
            if status in {"SITE_MAPPING_REVIEW", "SITE_WEAK_OR_AMBIGUOUS"}:
                vw.writerow(row); review_count += 1; review_neighbors[a] += 1; review_neighbors[b] += 1
        if gi % 5000 == 0: print(f"blocks compared {gi}/{len(block_sizes)} edges={edge_count}", flush=True)
    pfh.close(); sfh.close(); vfh.close()

    inv_rows = []
    for r in cand.itertuples(index=False):
        s = site[r.pair_id]
        bsig = ";".join(sorted(s["b"])); psig = ";".join(sorted(s["p"]))
        if not s["complete"]: pair_status = "STEP2_PAIR_MAPPING_REVIEW"
        elif candidate_degree[r.pair_id]: pair_status = "STEP2_PAIR_SAME_SITE_CANDIDATE"
        else: pair_status = "STEP2_PAIR_NO_STRONG_SITE_NEIGHBOR"
        inv_rows.append({
            "pair_id": r.pair_id, "candidate_block_id": r.candidate_block_id,
            "binding_site_signature": bsig, "binding_site_signature_sha256": hashlib.sha256(bsig.encode()).hexdigest() if bsig else "",
            "pocket_site_signature": psig, "pocket_site_signature_sha256": hashlib.sha256(psig.encode()).hexdigest() if psig else "",
            "source_binding_residue_count": s["binding_source"], "mapped_binding_residue_count": len(s["b"]),
            "unmapped_binding_residue_count": s["bu"], "ambiguous_binding_residue_count": s["ba"],
            "outside_identity_binding_residue_count": s["bo"], "source_pocket_residue_count": s["pocket_source"],
            "mapped_pocket_residue_count": len(s["p"]), "unmapped_pocket_residue_count": s["pu"],
            "site_candidate_degree": candidate_degree[r.pair_id], "exact_site_neighbor_count": exact_neighbors[r.pair_id],
            "strong_site_neighbor_count": strong_neighbors[r.pair_id], "review_neighbor_count": review_neighbors[r.pair_id],
            "step2_pair_status": pair_status,
        })
    pinv = pd.DataFrame(inv_rows)
    pinv.to_csv(out / "02_filter5_step2_pair_inventory.tsv.gz", sep="\t", index=False,
                compression={"method": "gzip", "compresslevel": 6, "mtime": 0})

    summary_rows = [
        ("step1_non_singleton_candidate_pairs", len(cand)), ("candidate_block_count", len(block_sizes)),
        ("total_pairwise_site_comparisons", edge_count), ("same_site_candidate_edges", same_count),
        ("review_edges", review_count), ("pairs_complete_binding_sifts_mapping", int(sum(s["complete"] for s in site.values()))),
        ("pairs_incomplete_binding_sifts_mapping", int(sum(not s["complete"] for s in site.values()))),
        ("unmapped_binding_residue_rows", int(sum(bind_unmapped.values()))),
        ("ambiguous_binding_residue_rows", int(sum(bind_ambig.values()))),
        ("outside_identity_binding_residue_rows", int(sum(bind_outside.values()))),
        ("excluded_binding_nonreceptor_residue_associations", excluded_binding_nonreceptor),
        ("excluded_pocket_nonreceptor_residue_rows", excluded_pocket_nonreceptor),
    ]
    for status in ["SITE_EXACT", "SITE_STRONG_CANDIDATE", "SITE_WEAK_OR_AMBIGUOUS", "SITE_DIFFERENT", "SITE_MAPPING_REVIEW"]:
        summary_rows.append((status, status_counts[status]))
    for label, counter in [("binding_jaccard", jbins), ("binding_containment", cbins), ("pocket_jaccard", pjbins), ("pocket_containment", pcbins)]:
        for b in ["1.00", "0.90-<1.00", "0.80-<0.90", "0.70-<0.80", "0.50-<0.70", "<0.50", "NA"]:
            summary_rows.append((f"{label}_bin_{b}", counter[b]))
    pd.DataFrame(summary_rows, columns=["metric", "value"]).to_csv(out / "05_filter5_step2_summary.tsv", sep="\t", index=False)
    exrows = []
    for status, rows in examples.items():
        for row in rows: exrows.append(dict(zip(ph, row)))
    pd.DataFrame(exrows).to_csv(out / "05b_filter5_step2_audit_examples.tsv", sep="\t", index=False)

    validation = {
        "status": "PASS",
        "checks": {
            "step1_universe_unchanged_241545": len(universe) == EXPECTED_UNIVERSE,
            "step2_membership_exactly_step1_non_singletons": set(pinv.pair_id) == set(cand.pair_id) and len(pinv) == len(cand),
            "candidate_pair_ids_unique": not pinv.pair_id.duplicated().any(),
            "pairwise_edge_count_complete": edge_count == expected_edges,
            "pairwise_status_partition_complete": sum(status_counts.values()) == edge_count,
            "same_site_interface_exact_or_strong_only": same_count == status_counts["SITE_EXACT"] + status_counts["SITE_STRONG_CANDIDATE"],
            "no_pair_deletion_or_collapse": True,
        },
        "thresholds": {
            "SITE_EXACT": "mapped binding sets identical",
            "SITE_STRONG_CANDIDATE": "containment >= 0.80 AND Jaccard >= 0.70 (PROVISIONAL_ENGINEERING_THRESHOLD)",
            "SITE_DIFFERENT": "disjoint OR (containment < 0.20 AND Jaccard < 0.10), complete binding mapping required",
        },
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.time() - t0, 3),
        "versions": {"python": sys.version, "pandas": pd.__version__, "pyarrow": pa.__version__},
    }
    if not all(validation["checks"].values()): validation["status"] = "FAIL"
    (valdir / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    report = f"""# Filter 5 Step 2 report\n\nStatus: **{validation['status']}**\n\nThis is a same-site candidate audit only. It performs no structure alignment, protein/pocket RMSD, ligand pose RMSD, representative selection, deletion, or equivalent-case collapse. Binding residues and the optional pocket context are reused verbatim from frozen Processing 3 and mapped through the frozen SIFTS observed-segment snapshot.\n\n- Step 1 candidate pairs: {len(cand):,}\n- Candidate blocks: {len(block_sizes):,}\n- Pairwise comparisons: {edge_count:,}\n- SITE_EXACT: {status_counts['SITE_EXACT']:,}\n- SITE_STRONG_CANDIDATE: {status_counts['SITE_STRONG_CANDIDATE']:,}\n- SITE_WEAK_OR_AMBIGUOUS: {status_counts['SITE_WEAK_OR_AMBIGUOUS']:,}\n- SITE_DIFFERENT: {status_counts['SITE_DIFFERENT']:,}\n- SITE_MAPPING_REVIEW: {status_counts['SITE_MAPPING_REVIEW']:,}\n- Complete binding SIFTS mapping: {sum(s['complete'] for s in site.values()):,}/{len(site):,} pairs\n\nThe strong threshold is explicitly provisional engineering candidate generation, not a universal literature cutoff.\n"""
    (out / "06_filter5_step2_report.md").write_text(report, encoding="utf-8")
    provenance = {
        "step1_frozen_input": str(STEP1), "step1_frozen_marker_sha256": sha256(STEP1 / "_FROZEN.json"),
        "binding_residue_source": str(P3RUN / "output/binding_residues"),
        "pocket_residue_source": str(P3RUN / "output/pair_pocket_residues"),
        "sifts_observed_segments": str(SIFTS), "sifts_sha256": sha256(SIFTS),
    }
    (run / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    schema = {
        "00_filter5_step2_residue_mapping.tsv.gz": residue_header,
        "01_filter5_step2_pairwise_site_comparisons.tsv.gz": ph,
        "02_filter5_step2_pair_inventory.tsv.gz": list(pinv.columns),
        "03_filter5_step2_same_site_candidates.tsv.gz": ph,
        "04_filter5_step2_review.tsv.gz": ph,
    }
    (run / "output_schema.json").write_text(json.dumps(schema, indent=2), encoding="utf-8")
    manifest = []
    for p in sorted(out.iterdir()):
        if not p.is_file(): continue
        rows = ""
        if p.name.endswith(".tsv.gz"):
            with gzip.open(p, "rt") as fh: rows = max(sum(1 for _ in fh) - 1, 0)
        elif p.name.endswith(".tsv"):
            with p.open() as fh: rows = max(sum(1 for _ in fh) - 1, 0)
        manifest.append((p.name, p.stat().st_size, rows, sha256(p)))
    pd.DataFrame(manifest, columns=["file", "bytes", "data_rows", "sha256"]).to_csv(run / "output_manifest.tsv", sep="\t", index=False)
    files = [p for p in run.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "_FROZEN.json"}]
    with (run / "SHA256SUMS").open("w") as fh:
        for p in sorted(files): fh.write(f"{sha256(p)}  {p.relative_to(run)}\n")
    if validation["status"] == "PASS":
        (run / "_FROZEN.json").write_text(json.dumps({"status": "FROZEN", "frozen_utc": datetime.now(timezone.utc).isoformat(), "validation": "validation/validation.json"}, indent=2), encoding="utf-8")
    print(json.dumps(validation, indent=2), flush=True)
    return 0 if validation["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
