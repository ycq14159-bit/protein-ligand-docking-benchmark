#!/usr/bin/env python3
"""Filter 4 Step 5: deterministic final decision, reconciliation, and release."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


VERSION = "filter4_step5_v1.0.0"
SCHEMA_VERSION = "filter4_step5_schema_v1.0.0"

EXPECTED_PATH_SUFFIXES = {
    "step1_run": "step_01_lattice_neighbor_search/runs/step01_full_v3",
    "step2_run": "step_02_biological_assembly_equivalence/runs/step02_full_v3",
    "step3_run": "step_03_direct_ligand_crystal_contact/runs/step03_full_v2",
    "step4_run": "step_04_binding_residue_mediated_crystal_contact/runs/step04_full_v2",
}

MAPPING = {
    "UPSTREAM_NO_NEIGHBOR": ("PASS", "NO_CRYSTALLOGRAPHIC_NEIGHBOR", ""),
    "UPSTREAM_NO_EXTERNAL_NEIGHBOR": ("PASS", "NO_EXTERNAL_CRYSTAL_NEIGHBOR", ""),
    "UPSTREAM_DIRECT_LIGAND_CONTACT_REJECT": ("REJECT", "DIRECT_LIGAND_CRYSTAL_CONTACT", "STEP3"),
    "SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT": ("REJECT", "BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT", "STEP4"),
    "SUCCESS_NO_BINDING_RESIDUE_CRYSTAL_CONTACT": ("PASS", "EXTERNAL_NEIGHBOR_NO_RELEVANT_CONTACT", ""),
    "BA_EQUIVALENCE_REVIEW": ("REVIEW", "BA_EQUIVALENCE_UNRESOLVED", ""),
}

EXPECTED_UPSTREAM = {
    "UPSTREAM_NO_NEIGHBOR": 7663,
    "UPSTREAM_NO_EXTERNAL_NEIGHBOR": 131764,
    "UPSTREAM_DIRECT_LIGAND_CONTACT_REJECT": 57580,
    "SUCCESS_BINDING_RESIDUE_CRYSTAL_CONTACT": 37285,
    "SUCCESS_NO_BINDING_RESIDUE_CRYSTAL_CONTACT": 102118,
    "BA_EQUIVALENCE_REVIEW": 2,
}

EXPECTED_DECISION_REASON = {
    ("PASS", "NO_CRYSTALLOGRAPHIC_NEIGHBOR"): 7663,
    ("PASS", "NO_EXTERNAL_CRYSTAL_NEIGHBOR"): 131764,
    ("PASS", "EXTERNAL_NEIGHBOR_NO_RELEVANT_CONTACT"): 102118,
    ("REJECT", "DIRECT_LIGAND_CRYSTAL_CONTACT"): 57580,
    ("REJECT", "BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT"): 37285,
    ("REVIEW", "BA_EQUIVALENCE_UNRESOLVED"): 2,
}

FINAL_COLUMNS = [
    "pair_id", "pdb_id", "assembly_id", "model_id",
    "step1_status", "step1_has_lattice_neighbor",
    "step2_status", "has_external_crystal_neighbor",
    "step3_status", "direct_ligand_crystal_contact_4A",
    "step4_status", "binding_residue_crystal_contact_4A",
    "n_external_instances", "n_external_ligand_6A", "n_external_pocket_6A",
    "n_direct_contact_instances", "n_direct_contact_units",
    "n_ligand_heavy_atoms_contacted_4A", "fraction_ligand_heavy_atoms_contacted_4A",
    "binding_residue_count", "n_crystal_bridged_binding_residues",
    "fraction_binding_residues_crystal_bridged",
    "n_binding_residue_contacting_external_instances",
    "filter4_decision", "filter4_reason", "reject_stage", "filter4_release_eligible",
]

DESCRIPTIONS = {
    "pair_id": "Frozen Filter 4 pair identifier.",
    "pdb_id": "PDB accession inherited from the frozen upstream inventory.",
    "assembly_id": "Selected assembly identifier parsed from the frozen pair identity.",
    "model_id": "Selected model identifier parsed from the frozen pair identity.",
    "step1_status": "Frozen Step 1 lattice-neighbour status.",
    "step1_has_lattice_neighbor": "Whether Step 1 found at least one neighbour within its frozen 6 A discovery shell.",
    "step2_status": "Frozen Step 2 biological-assembly equivalence status.",
    "has_external_crystal_neighbor": "Frozen Step 2 external crystallographic-neighbour flag.",
    "step3_status": "Frozen Step 3 direct-ligand-contact status.",
    "direct_ligand_crystal_contact_4A": "Frozen Step 3 pair-level direct ligand crystal-contact flag.",
    "step4_status": "Frozen Step 4 binding-residue-mediated-contact status.",
    "binding_residue_crystal_contact_4A": "Frozen Step 4 pair-level binding-residue crystal-contact flag.",
    "n_external_instances": "Frozen Step 2 number of external crystallographic instances.",
    "n_external_ligand_6A": "Frozen Step 2 number of external instances in the ligand 6 A shell.",
    "n_external_pocket_6A": "Frozen Step 2 number of external instances in the pocket 6 A shell.",
    "n_direct_contact_instances": "Frozen Step 3 number of direct-contact external instances.",
    "n_direct_contact_units": "Frozen Step 3 number of direct ligand-contact units.",
    "n_ligand_heavy_atoms_contacted_4A": "Frozen Step 3 count of contacted ligand heavy atoms.",
    "fraction_ligand_heavy_atoms_contacted_4A": "Frozen Step 3 fraction of ligand heavy atoms contacted.",
    "binding_residue_count": "Frozen Step 4 binding-residue count.",
    "n_crystal_bridged_binding_residues": "Frozen Step 4 count of crystal-bridged binding residues.",
    "fraction_binding_residues_crystal_bridged": "Frozen Step 4 fraction of binding residues crystal bridged.",
    "n_binding_residue_contacting_external_instances": "Frozen Step 4 number of external instances contacting binding residues.",
    "filter4_decision": "Deterministic final decision; allowed values PASS, REJECT, REVIEW.",
    "filter4_reason": "Single deterministic reason associated with the final Filter 4 decision.",
    "reject_stage": "STEP3 or STEP4 for rejected pairs; blank otherwise.",
    "filter4_release_eligible": "True only for FILTER4 PASS pairs released to Filter 5.",
}

SOURCE_STAGE = {
    **{x: "PAIR_IDENTITY" for x in ["pair_id", "pdb_id", "assembly_id", "model_id"]},
    **{x: "STEP1" for x in ["step1_status", "step1_has_lattice_neighbor"]},
    **{x: "STEP2" for x in ["step2_status", "has_external_crystal_neighbor", "n_external_instances", "n_external_ligand_6A", "n_external_pocket_6A"]},
    **{x: "STEP3" for x in ["step3_status", "direct_ligand_crystal_contact_4A", "n_direct_contact_instances", "n_direct_contact_units", "n_ligand_heavy_atoms_contacted_4A", "fraction_ligand_heavy_atoms_contacted_4A"]},
    **{x: "STEP4" for x in ["step4_status", "binding_residue_crystal_contact_4A", "binding_residue_count", "n_crystal_bridged_binding_residues", "fraction_binding_residues_crystal_bridged", "n_binding_residue_contacting_external_instances"]},
    **{x: "STEP5" for x in ["filter4_decision", "filter4_reason", "reject_stage", "filter4_release_eligible"]},
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path, block: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block): h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def truth(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().eq("true")


def physical_rows(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        return max(0, sum(1 for _ in fh) - 1)


def write_gzip(frame: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, sep="\t", index=False, compression="gzip", na_rep="", lineterminator="\n")
    os.replace(tmp, path)


def verify_upstream(cfg: dict) -> dict:
    evidence = {}
    for name, suffix in EXPECTED_PATH_SUFFIXES.items():
        path = Path(cfg["input"][name])
        normalized = path.as_posix().rstrip("/")
        if not normalized.endswith(suffix):
            raise RuntimeError(f"forbidden or unexpected {name}: {path}; expected suffix {suffix}")
        if not path.is_dir(): raise RuntimeError(f"missing upstream run: {path}")
        if name == "step1_run":
            validation = json.loads((path / "validation.json").read_text())
            metadata = json.loads((path / "run_metadata.json").read_text())
            if not validation.get("validation_pass") or metadata.get("status") != "VALIDATED":
                raise RuntimeError("Step1 is not validated")
            freeze_status = "VALIDATED_LEGACY_NO_FROZEN_MARKER"
        else:
            marker = json.loads((path / "_FROZEN.json").read_text())
            if marker.get("status") != "FROZEN" or not marker.get("validation_pass"):
                raise RuntimeError(f"{name} is not frozen")
            freeze_status = "FROZEN"
        evidence[name] = {"absolute_path": str(path), "run_name": path.name, "freeze_status": freeze_status,
                          "sha256sums_sha256": sha256(path / "SHA256SUMS")}
    if Path(cfg["input"]["step4_run"]).name == "step04_full_v1": raise RuntimeError("step04_full_v1 is explicitly forbidden")
    if Path(cfg["input"]["step3_run"]).name == "step03_full_v1": raise RuntimeError("step03_full_v1 is explicitly forbidden")
    return evidence


def build_final(source: pd.DataFrame) -> pd.DataFrame:
    unknown = sorted(set(source["step4_status"]) - set(MAPPING))
    if unknown: raise RuntimeError(f"unmapped Step4 statuses: {unknown}")
    parts = source["candidate_pair_id"].str.split("|")
    if not parts.map(len).ge(4).all(): raise RuntimeError("unparseable frozen pair identity")
    mapped = source["step4_status"].map(MAPPING)
    out = pd.DataFrame({
        "pair_id": source["candidate_pair_id"], "pdb_id": source["pdb_id"],
        "assembly_id": parts.str[2], "model_id": parts.str[3],
        "step1_status": source["step1_status"],
        "step1_has_lattice_neighbor": ~source["step1_status"].eq("NO_NEIGHBOR_WITHIN_6A"),
        "step2_status": source["step2_status"], "has_external_crystal_neighbor": truth(source["has_external_crystal_neighbor"]),
        "step3_status": source["step3_status"], "direct_ligand_crystal_contact_4A": truth(source["pair_direct_crystal_contact_4A"]),
        "step4_status": source["step4_status"], "binding_residue_crystal_contact_4A": truth(source["pair_binding_residue_crystal_contact_4A"]),
        "n_external_instances": source["external_crystal_instance_count"],
        "n_external_ligand_6A": source["external_ligand_6A_count"], "n_external_pocket_6A": source["external_pocket_6A_count"],
        "n_direct_contact_instances": source["n_external_instances_direct_contact_4A"], "n_direct_contact_units": source["n_contact_units_4A"],
        "n_ligand_heavy_atoms_contacted_4A": source["n_ligand_heavy_atoms_contacted_4A"],
        "fraction_ligand_heavy_atoms_contacted_4A": source["fraction_ligand_heavy_atoms_contacted_4A"],
        "binding_residue_count": source["binding_residue_count"], "n_crystal_bridged_binding_residues": source["n_crystal_bridged_binding_residues"],
        "fraction_binding_residues_crystal_bridged": source["fraction_binding_residues_crystal_bridged"],
        "n_binding_residue_contacting_external_instances": source["n_external_instances_contacting_binding_residues"],
        "filter4_decision": mapped.map(lambda x: x[0]), "filter4_reason": mapped.map(lambda x: x[1]),
        "reject_stage": mapped.map(lambda x: x[2]),
    })
    out["filter4_release_eligible"] = out["filter4_decision"].eq("PASS")
    return out[FINAL_COLUMNS]


def decision_summary(final: pd.DataFrame) -> pd.DataFrame:
    details = (final.groupby(["filter4_decision", "filter4_reason"], sort=False).size().rename("count").reset_index()
               .rename(columns={"filter4_decision": "decision", "filter4_reason": "reason"}))
    details.insert(0, "row_type", "DETAIL")
    order = {key: i for i, key in enumerate(EXPECTED_DECISION_REASON)}
    details["_order"] = [order.get((r.decision, r.reason), 999) for r in details.itertuples()]
    details = details.sort_values("_order").drop(columns="_order")
    totals = pd.DataFrame([{"row_type": "DECISION_TOTAL", "decision": d, "reason": "", "count": int((final["filter4_decision"] == d).sum())}
                           for d in ["PASS", "REJECT", "REVIEW"]])
    grand = pd.DataFrame([{"row_type": "GRAND_TOTAL", "decision": "TOTAL", "reason": "", "count": len(final)}])
    return pd.concat([details, totals, grand], ignore_index=True)


def make_schema(datasets: dict[str, pd.DataFrame]) -> dict:
    result = {}
    for name, frame in datasets.items():
        cols = []
        for col in frame.columns:
            nullable = bool(frame[col].isna().any() or (frame[col].dtype == object and frame[col].astype(str).eq("").any()))
            cols.append({"column": col, "type": str(frame[col].dtype), "nullable": nullable,
                         "description": DESCRIPTIONS.get(col, "Release summary or inherited audit field."),
                         "source_stage": SOURCE_STAGE.get(col, "STEP5")})
        result[name] = {"columns": cols}
    return {"schema_version": SCHEMA_VERSION, "filter4_decision_allowed_values": ["PASS", "REJECT", "REVIEW"], "datasets": result}


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True, type=Path); ap.add_argument("--run-dir", required=True, type=Path)
    args = ap.parse_args(); cfg = yaml.safe_load(args.config.read_text()); run = args.run_dir
    if run.exists() and any(run.iterdir()): raise SystemExit(f"run directory not empty: {run}")
    output = run / "output"; output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.config, run / "config_snapshot.yaml"); shutil.copy2(Path(__file__), run / "executed_runner.py")
    provenance = verify_upstream(cfg); step4 = Path(cfg["input"]["step4_run"]); step1 = Path(cfg["input"]["step1_run"])
    source = pd.read_csv(step4 / "output/04_pair_binding_residue_contact_inventory.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    upstream_counts = Counter(source["step4_status"])
    upstream_exact = dict(upstream_counts) == EXPECTED_UPSTREAM
    if not upstream_exact: raise RuntimeError(f"Step4 frozen status counts differ: {dict(upstream_counts)}")
    original = pd.read_csv(step1 / "output/01_pair_step1_inventory.tsv.gz", sep="\t", usecols=["candidate_pair_id"], dtype=str, keep_default_na=False)
    final = build_final(source); passed = final[final["filter4_decision"].eq("PASS")].copy()
    rejected = final[final["filter4_decision"].eq("REJECT")].copy(); review = final[final["filter4_decision"].eq("REVIEW")].copy()
    summary = decision_summary(final)
    write_gzip(final, output / "01_filter4_final_pair_inventory.tsv.gz")
    write_gzip(passed, output / "02_filter4_pass_pairs.tsv.gz")
    write_gzip(rejected, output / "03_filter4_rejected_pairs.tsv.gz")
    review.to_csv(output / "04_filter4_review_pairs.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")
    summary.to_csv(output / "05_filter4_decision_summary.tsv", sep="\t", index=False, lineterminator="\n")

    direct = rejected[rejected["filter4_reason"].eq("DIRECT_LIGAND_CRYSTAL_CONTACT")]
    mediated = rejected[rejected["filter4_reason"].eq("BINDING_RESIDUE_MEDIATED_CRYSTAL_CONTACT")]
    p_no_neighbour = passed[passed["filter4_reason"].eq("NO_CRYSTALLOGRAPHIC_NEIGHBOR")]
    p_no_external = passed[passed["filter4_reason"].eq("NO_EXTERNAL_CRYSTAL_NEIGHBOR")]
    p_no_relevant = passed[passed["filter4_reason"].eq("EXTERNAL_NEIGHBOR_NO_RELEVANT_CONTACT")]
    final_reason_counts = Counter(map(tuple, final[["filter4_decision", "filter4_reason"]].to_numpy()))
    final_ids, original_ids = set(final["pair_id"]), set(original["candidate_pair_id"])
    checks = {
        "exact_step4_source_full_v2": step4.name == "step04_full_v2",
        "exact_step3_source_full_v2": Path(cfg["input"]["step3_run"]).name == "step03_full_v2",
        "upstream_status_counts_exact": upstream_exact,
        "final_rows_336412": len(final) == 336412,
        "duplicate_pair_id_zero": not final["pair_id"].duplicated().any(),
        "silent_drop_zero": final_ids == original_ids and len(final) == len(original),
        "final_pair_set_equals_filter4_input": final_ids == original_ids,
        "decision_allowed_values_only": set(final["filter4_decision"]) == {"PASS", "REJECT", "REVIEW"},
        "decision_mapping_exact": all((row.filter4_decision, row.filter4_reason, row.reject_stage) == MAPPING[row.step4_status] for row in final.itertuples()),
        "decision_reason_counts_exact": dict(final_reason_counts) == EXPECTED_DECISION_REASON,
        "pass_241545": len(passed) == 241545,
        "reject_94865": len(rejected) == 94865,
        "review_2": len(review) == 2,
        "decision_partition_336412": len(passed) + len(rejected) + len(review) == 336412,
        "release_only_pass": len(passed) == 241545 and passed["filter4_decision"].eq("PASS").all() and passed["filter4_release_eligible"].all(),
        "reject_release_ineligible": (~rejected["filter4_release_eligible"]).all(),
        "review_release_ineligible": (~review["filter4_release_eligible"]).all(),
        "direct_reject_consistency": len(direct) == 57580 and direct["direct_ligand_crystal_contact_4A"].all() and direct["reject_stage"].eq("STEP3").all(),
        "mediated_reject_consistency": len(mediated) == 37285 and (~mediated["direct_ligand_crystal_contact_4A"]).all() and mediated["binding_residue_crystal_contact_4A"].all() and mediated["reject_stage"].eq("STEP4").all(),
        "reject_reasons_disjoint": set(direct["pair_id"]).isdisjoint(set(mediated["pair_id"])),
        "no_neighbour_pass_consistency": len(p_no_neighbour) == 7663 and (~p_no_neighbour["step1_has_lattice_neighbor"]).all(),
        "no_external_pass_consistency": len(p_no_external) == 131764 and p_no_external["step1_has_lattice_neighbor"].all() and (~p_no_external["has_external_crystal_neighbor"]).all() and pd.to_numeric(p_no_external["n_external_instances"]).eq(0).all(),
        "external_no_relevant_pass_consistency": len(p_no_relevant) == 102118 and p_no_relevant["has_external_crystal_neighbor"].all() and (~p_no_relevant["direct_ligand_crystal_contact_4A"]).all() and (~p_no_relevant["binding_residue_crystal_contact_4A"]).all(),
        "review_consistency": len(review) == 2 and review["step4_status"].eq("BA_EQUIVALENCE_REVIEW").all() and review["filter4_reason"].eq("BA_EQUIVALENCE_UNRESOLVED").all() and set(review["pdb_id"]) == {"3gbn"},
        "new_review_count_zero": len(review) == int(upstream_counts["BA_EQUIVALENCE_REVIEW"]),
    }
    serialized = {}
    for path, expected in [(output/"01_filter4_final_pair_inventory.tsv.gz", len(final)), (output/"02_filter4_pass_pairs.tsv.gz", len(passed)),
                           (output/"03_filter4_rejected_pairs.tsv.gz", len(rejected)), (output/"04_filter4_review_pairs.tsv", len(review)),
                           (output/"05_filter4_decision_summary.tsv", len(summary))]:
        actual = physical_rows(path); serialized[path.name] = {"expected_rows": expected, "physical_rows": actual, "match": actual == expected}
    checks["serialized_physical_row_counts_match"] = all(x["match"] for x in serialized.values())
    validation = {"validated_at": utc(), "validation_pass": all(bool(x) for x in checks.values()), "checks": checks,
                  "upstream_status_counts": dict(upstream_counts), "decision_counts": dict(Counter(final["filter4_decision"])),
                  "decision_reason_counts": [{"decision": k[0], "reason": k[1], "count": v} for k, v in final_reason_counts.items()],
                  "serialization_row_counts": serialized, "input_pair_count": len(original), "output_pair_count": len(final)}
    atomic_json(run / "validation.json", validation)

    report = f"""# Filter 4 Step 5 — Final Crystal-Packing Decision and Release

Status: `{'PASS' if validation['validation_pass'] else 'FAIL'}`  
Decision policy: `filter4_strict_benchmark_policy_v1.0.0`

## Purpose

Filter 4 audits whether the experimentally observed protein–ligand pose is exposed to non-biological crystallographic packing contacts that directly involve the ligand or its frozen binding residues.

## Step hierarchy

- Step 1: local crystallographic lattice-neighbour search
- Step 2: selected Biological Assembly equivalence
- Step 3: direct ligand crystal contact
- Step 4: binding-residue-mediated crystal contact
- Step 5: deterministic final benchmark decision and release

## Frozen decision policy

```text
336,412 Filter 4 input pairs
            |
            v
External crystallographic neighbour?
       |                    |
      NO                   YES
       |                    |
      PASS          Direct ligand contact?
                         |          |
                        YES        NO
                         |          |
                      REJECT   Contact with frozen
                               binding residue?
                                  |       |
                                 YES     NO
                                  |       |
                               REJECT    PASS

BA equivalence unresolved -> REVIEW
```

- PASS: {len(passed):,}
- REJECT: {len(rejected):,}
- REVIEW: {len(review):,}

Direct exclusion (`external crystal -> ligand`) is the most direct crystal-packing evidence. Binding-residue-mediated exclusion (`external crystal -> frozen binding residue -> ligand`) is the shortest one-hop indirect structural constraint. Under the strict benchmark inclusion policy, these structures are excluded because the crystallographic environment directly constrains either the ligand itself or a residue already participating in ligand binding. This does not prove that crystal packing necessarily changed the ligand pose.

No multi-hop residue propagation, full-pocket propagation, vdW secondary criterion, PLIP/Arpeggio interaction classification, energetic scoring, SASA-based severity, PISA energy, or crystal-packing severity score was used. Step 5 performs no structural or contact calculation.

The formal downstream interface for Filter 5 is `02_filter4_pass_pairs.tsv.gz`.
"""
    (output / "06_filter4_final_report.md").write_text(report, encoding="utf-8")
    metadata = {"filter_name": "filter_04_crystal_packing_influence", "filter_version": "1.0.0", "decision_policy_version": "filter4_strict_benchmark_policy_v1.0.0",
                "input_pair_count": len(final), "pass_pair_count": len(passed), "reject_pair_count": len(rejected), "review_pair_count": len(review),
                "step1_source": provenance["step1_run"], "step2_source": provenance["step2_run"], "step3_source": provenance["step3_run"], "step4_source": provenance["step4_run"],
                "created_at": utc(), "validation_pass": validation["validation_pass"],
                "downstream_interface": str(output / "02_filter4_pass_pairs.tsv.gz")}
    atomic_json(run / "filter4_release_metadata.json", metadata); atomic_json(run / "input_provenance.json", provenance)
    atomic_json(run / "output_schema.json", make_schema({"01_filter4_final_pair_inventory": final, "02_filter4_pass_pairs": passed,
                "03_filter4_rejected_pairs": rejected, "04_filter4_review_pairs": review, "05_filter4_decision_summary": summary}))
    files = [p for p in run.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "output_manifest.tsv", "_FROZEN.json"}]
    manifest = [{"relative_path": p.relative_to(run).as_posix(), "size_bytes": p.stat().st_size, "sha256": sha256(p),
                 "schema_version": SCHEMA_VERSION, "generated_by": VERSION} for p in sorted(files)]
    pd.DataFrame(manifest).to_csv(run / "output_manifest.tsv", sep="\t", index=False, lineterminator="\n")
    checksum_files = [p for p in run.rglob("*") if p.is_file() and p.name not in {"SHA256SUMS", "_FROZEN.json"}]
    with (run / "SHA256SUMS").open("w", encoding="utf-8") as fh:
        for p in sorted(checksum_files): fh.write(f"{sha256(p)}  {p.relative_to(run).as_posix()}\n")
    if validation["validation_pass"]:
        atomic_json(run / "_FROZEN.json", {"stage": "filter_04_step_05", "status": "FROZEN", "validation_pass": True,
                    "input_pairs": len(final), "pass_pairs": len(passed), "reject_pairs": len(rejected), "review_pairs": len(review),
                    "decision_policy_version": metadata["decision_policy_version"], "frozen_at": utc(), "sha256sums_sha256": sha256(run / "SHA256SUMS")})
    print(json.dumps(validation, indent=2, default=lambda x: x.item() if hasattr(x, "item") else str(x)))
    if not validation["validation_pass"]: raise SystemExit(2)


if __name__ == "__main__":
    main()
