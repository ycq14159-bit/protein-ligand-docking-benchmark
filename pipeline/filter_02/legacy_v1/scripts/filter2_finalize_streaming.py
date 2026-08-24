from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import sqlite3
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import filter2_pipeline as p


OUT = p.OUT
FULL = OUT / "full"
REPORTS = OUT / "reports"
RELEASE = OUT / "release"
VALIDATION = OUT / "validation"


def utc():
    return datetime.now(timezone.utc).isoformat()


def open_writer(path: Path, fields: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    handle = gzip.open(tmp, "wt", encoding="utf-8", newline="") if path.suffix == ".gz" else tmp.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    return handle, writer, tmp


def finish_writer(item, target: Path):
    handle, _, tmp = item
    handle.close()
    os.replace(tmp, target)


def insert_batches(db, sql: str, rows, size=50000):
    batch = []
    for row in rows:
        batch.append(row)
        if len(batch) >= size:
            db.executemany(sql, batch)
            batch.clear()
    if batch:
        db.executemany(sql, batch)


def historical_file():
    preferred = Path("/root/autodl-tmp/vs_benchmark/data_interaction_refinement_arpeggio_v2/main_benchmark_candidate_v2.csv")
    if preferred.exists():
        return preferred
    ordinary_queue = Path("/root/autodl-tmp/vs_benchmark/data_arpeggio_v2_full/manifests/arpeggio_ordinary_main_queue.tsv")
    if ordinary_queue.exists():
        return ordinary_queue
    candidates = list(Path("/root/autodl-tmp/vs_benchmark").rglob("candidate_pairs_classified_v1.csv"))
    return candidates[0] if candidates else None


def raw_check(row):
    path = Path(row["mmcif_path"])
    if not path.exists():
        return "missing"
    if path.stat().st_size != int(row["mmcif_file_size"]):
        return "size_mismatch"
    return "pass" if p.sha(path) == row["mmcif_checksum"] else "checksum_mismatch"


def main():
    started = utc()
    required = {
        "entries": FULL / "filter_2_entries.tsv.gz",
        "sources": FULL / "filter_2_sources.tsv.gz",
        "assemblies": FULL / "filter_2_assemblies.tsv.gz",
        "conformers": FULL / "filter_2_conformers.tsv.gz",
        "parents": FULL / "filter_2_parents.tsv.gz",
        "connections": FULL / "filter_2_connections.tsv.gz",
        "categories": FULL / "filter_2_categories.tsv.gz",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise SystemExit("Missing merged checkpoint tables: " + ", ".join(missing))

    aliases = {
        "filter_2_entry_inventory.tsv.gz": "filter_2_entries.tsv.gz",
        "filter_2_source_component_instances.tsv.gz": "filter_2_sources.tsv.gz",
        "filter_2_assembly_component_instances.tsv.gz": "filter_2_assemblies.tsv.gz",
        "filter_2_component_conformers.tsv.gz": "filter_2_conformers.tsv.gz",
        "filter_2_parent_mapping.tsv.gz": "filter_2_parents.tsv.gz",
        "filter_2_struct_conn_links.tsv.gz": "filter_2_connections.tsv.gz",
    }
    for dst, src in aliases.items():
        shutil.copy2(FULL / src, FULL / dst)

    used_ids = {row["label_comp_id"] for row in p.iter_tsv(required["sources"])}
    ccd = {row["original_component_id"]: row for row in p.iter_tsv(OUT / "references/ccd_component_cache.tsv.gz") if row["original_component_id"] in used_ids}
    components = [ccd[x] for x in sorted(used_ids) if x in ccd]
    components.extend({
        "original_component_id": x,
        "resolved_ccd_id": "",
        "ccd_identity_status": "ccd_missing",
        "chemical_entity_class": "unknown",
        "artifact_prior": "unknown",
        "classification_reason": "not_in_frozen_ccd",
        "rule_version": p.RULE_VERSION,
    } for x in sorted(used_ids - set(ccd)))
    p.write_tsv(FULL / "filter_2_component_classification.tsv.gz", components, p.COMP_FIELDS, True)
    p.write_tsv(RELEASE / "filter_2_component_classification.tsv.gz", components, p.COMP_FIELDS, True)

    db_path = VALIDATION / "filter_2_keys.sqlite"
    if db_path.exists():
        db_path.unlink()
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA temp_store=FILE")
    db.execute("PRAGMA cache_size=-1048576")
    db.execute("CREATE TABLE source_keys(k TEXT PRIMARY KEY) WITHOUT ROWID")
    db.execute("CREATE TABLE assembly_keys(k TEXT PRIMARY KEY, source_k TEXT, pdb_id TEXT, assembly_id TEXT) WITHOUT ROWID")
    db.execute("CREATE TABLE conformer_keys(k TEXT PRIMARY KEY) WITHOUT ROWID")
    db.execute("CREATE TABLE component_keys(k TEXT PRIMARY KEY) WITHOUT ROWID")
    db.execute("CREATE TABLE qualified_assemblies(pdb_id TEXT, assembly_id TEXT, PRIMARY KEY(pdb_id,assembly_id)) WITHOUT ROWID")
    db.execute("CREATE TABLE route_map(kind TEXT, pdb_id TEXT, component_id TEXT, route TEXT, PRIMARY KEY(kind,pdb_id,component_id,route)) WITHOUT ROWID")
    db.executemany("INSERT INTO component_keys VALUES(?)", ((r["original_component_id"],) for r in components))
    db.executemany("INSERT INTO qualified_assemblies VALUES(?,?)", ((r["pdb_id"], r["assembly_id"]) for r in p.iter_tsv(OUT / "inputs/filter_1_receptor_qualified_assemblies_snapshot.tsv.gz")))
    db.commit()

    full_route_names = {
        "ordinary_small_molecule_candidate": "ordinary_candidates",
        "artifact_review": "artifact_review",
        "polymer_or_modified_residue": "polymer_modified_residues",
        "unresolved_review": "unresolved_review",
    }
    full_writers = {route: open_writer(FULL / f"filter_2_{name}.tsv.gz", p.SOURCE_FIELDS) for route, name in full_route_names.items()}
    full_writers["__special__"] = open_writer(FULL / "filter_2_special_candidates.tsv.gz", p.SOURCE_FIELDS)
    release_writers = {
        "ordinary_small_molecule_candidate": open_writer(RELEASE / "filter_2_ordinary_component_instances.tsv.gz", p.SOURCE_FIELDS),
        "artifact_review": open_writer(RELEASE / "filter_2_artifact_review.tsv.gz", p.SOURCE_FIELDS),
        "unresolved_review": open_writer(RELEASE / "filter_2_unresolved_review.tsv.gz", p.SOURCE_FIELDS),
        "__special__": open_writer(RELEASE / "filter_2_special_component_instances.tsv.gz", p.SOURCE_FIELDS),
    }

    route_stats = defaultdict(lambda: {"source": 0, "assembly": 0, "ccd": set(), "pdb": set(), "qualified_assembly": set()})
    polymer_counts = Counter()
    artifact_counts = Counter()
    covalent_counts = Counter()
    metal_counts = Counter()
    source_total = missing_route = ordinary_count = special_count = artifact_count = polymer_count = unresolved_count = 0
    covalent_linked_count = metal_inorganic_count = 0
    contamination = Counter()
    key_batch = []
    route_batch = []
    for row in p.iter_tsv(required["sources"]):
        source_total += 1
        route = row["filter_2_route"]
        if not route:
            missing_route += 1
        route_stats[route]["source"] += 1
        route_stats[route]["ccd"].add(row["resolved_ccd_id"])
        route_stats[route]["pdb"].add(row["pdb_id"])
        polymer_counts[row["polymer_context"]] += 1
        artifact_counts[row["artifact_prior"]] += 1
        covalent_counts[row["instance_covalent_link_status"]] += 1
        metal_counts[row["instance_metal_link_status"]] += 1
        if row["instance_covalent_link_status"] == "declared_receptor_covalent":
            covalent_linked_count += 1
        if row["chemical_entity_class"] in {"metal_or_inorganic", "organometallic"}:
            metal_inorganic_count += 1
        key_batch.append((row["source_component_instance_id"],))
        route_batch.append(("source", row["pdb_id"], row["label_comp_id"], route))
        if len(key_batch) >= 50000:
            db.executemany("INSERT OR IGNORE INTO source_keys VALUES(?)", key_batch)
            db.executemany("INSERT OR IGNORE INTO route_map VALUES(?,?,?,?)", route_batch)
            db.commit(); key_batch.clear(); route_batch.clear()
        if route in full_writers:
            full_writers[route][1].writerow(row)
        elif route.endswith("_special"):
            full_writers["__special__"][1].writerow(row)
        if route in release_writers:
            release_writers[route][1].writerow(row)
        elif route.endswith("_special"):
            release_writers["__special__"][1].writerow(row)
        if route == "ordinary_small_molecule_candidate":
            ordinary_count += 1
            checks = {
                "ordinary_polymer_residue": row["polymer_context"] != "independent_nonpolymer",
                "ordinary_modified_polymer_residue": row["modified_residue_status"] == "true",
                "ordinary_short_peptide": row["polymer_context"] == "short_peptide",
                "ordinary_DNA_RNA": row["polymer_context"] in {"rna_polymer", "dna_polymer", "hybrid_nucleic_acid"},
                "ordinary_branched_glycan": row["polymer_context"] == "branched_glycan",
                "ordinary_water": row["chemical_entity_class"] == "water",
                "ordinary_metal_inorganic": row["chemical_entity_class"] in {"metal_or_inorganic", "organometallic"},
                "ordinary_unresolved_CCD": row["instance_status"] != "resolved",
                "ordinary_explicit_receptor_covalent": row["instance_covalent_link_status"] == "declared_receptor_covalent",
            }
            contamination.update({k: int(v) for k, v in checks.items()})
        elif route.endswith("_special"):
            special_count += 1
        elif route == "artifact_review": artifact_count += 1
        elif route == "polymer_or_modified_residue": polymer_count += 1
        elif route == "unresolved_review": unresolved_count += 1
    if key_batch:
        db.executemany("INSERT OR IGNORE INTO source_keys VALUES(?)", key_batch)
        db.executemany("INSERT OR IGNORE INTO route_map VALUES(?,?,?,?)", route_batch)
        db.commit()
    for route, writer in full_writers.items():
        finish_writer(writer, FULL / ("filter_2_special_candidates.tsv.gz" if route == "__special__" else f"filter_2_{full_route_names[route]}.tsv.gz"))
    release_targets = {
        "ordinary_small_molecule_candidate": "filter_2_ordinary_component_instances.tsv.gz",
        "artifact_review": "filter_2_artifact_review.tsv.gz",
        "unresolved_review": "filter_2_unresolved_review.tsv.gz",
        "__special__": "filter_2_special_component_instances.tsv.gz",
    }
    for route, writer in release_writers.items():
        finish_writer(writer, RELEASE / release_targets[route])

    ordinary_assembly_writer = open_writer(RELEASE / "filter_2_ordinary_assembly_component_instances.tsv.gz", p.ASSEMBLY_FIELDS)
    assembly_counts = Counter()
    assembly_total = assembly_missing_route = ordinary_assembly_count = 0
    assembly_key_batch = []
    route_batch = []
    for row in p.iter_tsv(required["assemblies"]):
        assembly_total += 1
        route = row["filter_2_route"]
        if not route: assembly_missing_route += 1
        assembly_counts[row["assembly_membership_status"]] += 1
        route_stats[route]["assembly"] += 1
        route_stats[route]["qualified_assembly"].add(row["pdb_id"] + "|" + row["assembly_id"])
        assembly_key_batch.append((row["assembly_component_instance_id"], row["source_component_instance_id"], row["pdb_id"], row["assembly_id"]))
        route_batch.append(("assembly", row["pdb_id"], row["resolved_ccd_id"], route))
        if len(assembly_key_batch) >= 50000:
            db.executemany("INSERT OR IGNORE INTO assembly_keys VALUES(?,?,?,?)", assembly_key_batch)
            db.executemany("INSERT OR IGNORE INTO route_map VALUES(?,?,?,?)", route_batch)
            db.commit(); assembly_key_batch.clear(); route_batch.clear()
        if route == "ordinary_small_molecule_candidate":
            ordinary_assembly_writer[1].writerow(row)
            ordinary_assembly_count += 1
    if assembly_key_batch:
        db.executemany("INSERT OR IGNORE INTO assembly_keys VALUES(?,?,?,?)", assembly_key_batch)
        db.executemany("INSERT OR IGNORE INTO route_map VALUES(?,?,?,?)", route_batch)
        db.commit()
    finish_writer(ordinary_assembly_writer, RELEASE / "filter_2_ordinary_assembly_component_instances.tsv.gz")

    conformer_total = 0
    conf_batch = []
    for row in p.iter_tsv(required["conformers"]):
        conformer_total += 1
        conf_batch.append((row["component_conformer_id"],))
        if len(conf_batch) >= 50000:
            db.executemany("INSERT OR IGNORE INTO conformer_keys VALUES(?)", conf_batch); db.commit(); conf_batch.clear()
    if conf_batch:
        db.executemany("INSERT OR IGNORE INTO conformer_keys VALUES(?)", conf_batch); db.commit()

    entry_total = parse_success = parse_failed = missing_terminal = water_total = 0
    entry_status_counts = Counter()
    failure_counts = Counter()
    water_writer = open_writer(FULL / "filter_2_excluded_water_summary.tsv.gz", ["pdb_id", "water_count"])
    for row in p.iter_tsv(required["entries"]):
        entry_total += 1
        parse_success += row["parse_status"] == "success"
        parse_failed += row["parse_status"] == "failed"
        missing_terminal += not bool(row["entry_status"])
        entry_status_counts[row["entry_status"]] += 1
        if row["entry_status"] != "pass": failure_counts[row["terminal_reason"]] += 1
        count = int(row["water_count"])
        water_total += count
        if count: water_writer[1].writerow({"pdb_id": row["pdb_id"], "water_count": count})
    finish_writer(water_writer, FULL / "filter_2_excluded_water_summary.tsv.gz")

    distinct_source = db.execute("SELECT COUNT(*) FROM source_keys").fetchone()[0]
    distinct_assembly = db.execute("SELECT COUNT(*) FROM assembly_keys").fetchone()[0]
    distinct_conformer = db.execute("SELECT COUNT(*) FROM conformer_keys").fetchone()[0]
    distinct_component = db.execute("SELECT COUNT(*) FROM component_keys").fetchone()[0]
    duplicates = {
        "component": len(components) - distinct_component,
        "source": source_total - distinct_source,
        "assembly": assembly_total - distinct_assembly,
        "conformer": conformer_total - distinct_conformer,
    }
    missing_source_fk = db.execute("SELECT COUNT(*) FROM assembly_keys a LEFT JOIN source_keys s ON a.source_k=s.k WHERE s.k IS NULL").fetchone()[0]
    missing_qualified_assembly_fk = db.execute("SELECT COUNT(*) FROM assembly_keys a LEFT JOIN qualified_assemblies q ON a.pdb_id=q.pdb_id AND a.assembly_id=q.assembly_id WHERE q.pdb_id IS NULL").fetchone()[0]

    route_report = []
    for route, stat in sorted(route_stats.items()):
        route_report.append({
            "filter_2_route": route,
            "source_component_instance_count": stat["source"],
            "assembly_component_instance_count": stat["assembly"],
            "unique_ccd_count": len(stat["ccd"] - {""}),
            "unique_pdb_entry_count": len(stat["pdb"]),
            "unique_qualified_assembly_count": len(stat["qualified_assembly"]),
        })
    p.write_tsv(REPORTS / "filter_2_route_distribution.tsv", route_report, ["filter_2_route","source_component_instance_count","assembly_component_instance_count","unique_ccd_count","unique_pdb_entry_count","unique_qualified_assembly_count"])
    distributions = [
        ("filter_2_entry_flow.tsv", "entry_status", entry_status_counts),
        ("filter_2_artifact_prior_distribution.tsv", "artifact_prior", artifact_counts),
        ("filter_2_polymer_context_distribution.tsv", "polymer_context", polymer_counts),
        ("filter_2_covalent_status_distribution.tsv", "instance_covalent_link_status", covalent_counts),
        ("filter_2_metal_status_distribution.tsv", "instance_metal_link_status", metal_counts),
        ("filter_2_assembly_membership_distribution.tsv", "assembly_membership_status", assembly_counts),
        ("filter_2_failure_reason_distribution.tsv", "terminal_reason", failure_counts),
    ]
    for name, key, counts in distributions:
        p.write_tsv(REPORTS / name, [{key: k, "count": v} for k, v in sorted(counts.items())], [key, "count"])
    component_distributions = [
        ("filter_2_component_class_distribution.tsv", "chemical_entity_class"),
        ("filter_2_ccd_status_distribution.tsv", "ccd_identity_status"),
        ("filter_2_rdkit_status_distribution.tsv", "rdkit_parse_status"),
        ("filter_2_element_distribution.tsv", "element_set"),
        ("filter_2_fragment_distribution.tsv", "fragment_count"),
    ]
    for name, key in component_distributions:
        counts = Counter(row.get(key, "") for row in components)
        p.write_tsv(REPORTS / name, [{key: k, "count": v} for k, v in sorted(counts.items())], [key, "count"])

    cross_count = 0
    cross_summary = Counter()
    cross_writer = open_writer(REPORTS / "filter_2_historical_crosswalk.tsv", ["pdb_id","component_id","historical_candidate_status","new_source_component_route","new_assembly_component_route","route_agreement","route_disagreement","disagreement_reason"])
    hist = historical_file()
    if hist:
        opener = gzip.open if str(hist).endswith(".gz") else open
        delimiter = "," if hist.suffix == ".csv" else "\t"
        with opener(hist, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter=delimiter):
                pid = (row.get("pdb_id") or "").lower()
                cid = (row.get("ligand_id") or row.get("ligand_component_id") or row.get("component_id") or "").upper()
                sr = [x[0] for x in db.execute("SELECT route FROM route_map WHERE kind='source' AND pdb_id=? AND component_id=? ORDER BY route", (pid, cid))]
                ar = [x[0] for x in db.execute("SELECT route FROM route_map WHERE kind='assembly' AND pdb_id=? AND component_id=? ORDER BY route", (pid, cid))]
                mapped = bool(sr)
                historical_status = row.get("preliminary_category") or row.get("refined_category") or row.get("final_action") or "historical_candidate"
                out = {"pdb_id": pid, "component_id": cid, "historical_candidate_status": historical_status, "new_source_component_route": ",".join(sr), "new_assembly_component_route": ",".join(ar), "route_agreement": "mapped" if mapped else "not_directly_comparable", "route_disagreement": "" if mapped else "not_found", "disagreement_reason": "" if mapped else "instance_model_changed_or_not_in_filter1_scope"}
                cross_writer[1].writerow(out); cross_count += 1; cross_summary[out["route_disagreement"] or "mapped"] += 1
    finish_writer(cross_writer, REPORTS / "filter_2_historical_crosswalk.tsv")
    db.commit(); db.close()

    raw_rows = list(p.iter_tsv(OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz"))
    with ThreadPoolExecutor(max_workers=16) as pool:
        raw_counts = Counter(pool.map(raw_check, raw_rows, chunksize=64))
    raw_audit = {"checked": len(raw_rows), **dict(raw_counts)}
    (VALIDATION / "raw_mmcif_checksum_audit.json").write_text(json.dumps(raw_audit, indent=2) + "\n")
    checksum_mismatch = raw_counts["missing"] + raw_counts["size_mismatch"] + raw_counts["checksum_mismatch"]

    processing_modified = p.sha(p.P1 / "release/processing_1_mmcif_index.tsv.gz") != p.sha(OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz")
    filter1_modified = p.sha(p.ENTRY_INPUT) != p.sha(OUT / "inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")
    validation = {
        "input_entries": 248037,
        "input_qualified_assemblies": 360611,
        "entry_inventory_rows": entry_total,
        "parse_success": parse_success,
        "parse_failed": parse_failed,
        "duplicate_component_identity_key": duplicates["component"],
        "duplicate_source_component_instance_key": duplicates["source"],
        "duplicate_assembly_component_instance_key": duplicates["assembly"],
        "duplicate_conformer_key": duplicates["conformer"],
        "assembly_missing_source_mapping": missing_source_fk,
        "assembly_missing_qualified_assembly_mapping": missing_qualified_assembly_fk,
        "missing_component_route": missing_route + assembly_missing_route,
        "missing_terminal_status": missing_terminal,
        "silent_drop": 248037 - entry_total,
        "ordinary_contamination": dict(contamination),
        "processing_1_modified": processing_modified,
        "filter_1_modified": filter1_modified,
        "raw_mmcif_modified": checksum_mismatch > 0,
        "raw_mmcif_checksum_audit": raw_audit,
        "historical_directories_modified": False,
        "assembly_coordinate_materialization_started": False,
        "pair_construction_started": False,
        "distance_calculation_started": False,
        "interaction_annotation_started": False,
        "structure_quality_filtering_started": False,
        "checksum_mismatch": checksum_mismatch,
    }
    validation["release_validation_pass"] = (
        entry_total == 248037 and validation["silent_drop"] == 0 and all(v == 0 for v in duplicates.values())
        and validation["missing_component_route"] == 0 and missing_terminal == 0
        and missing_source_fk == 0 and missing_qualified_assembly_fk == 0
        and all(v == 0 for v in contamination.values())
        and not processing_modified and not filter1_modified and checksum_mismatch == 0
    )
    progress = json.loads((OUT / "checkpoints/progress.json").read_text())
    rdkit_parse_fail = sum(row.get("rdkit_parse_status") != "pass" for row in components)
    rdkit_sanitize_fail = sum(row.get("rdkit_sanitize_status") != "pass" for row in components)
    ccd_invalid = sum(row.get("ccd_identity_status") not in {"exact_ccd_match", "obsolete_id_resolved"} for row in components)
    summary = {
        "finalize_start": started,
        "finalize_end": utc(),
        "full_start": progress["start"],
        "full_end": progress["end"],
        "runtime_seconds": progress["elapsed_seconds"],
        "input_entries": 248037,
        "input_qualified_assemblies": 360611,
        "parse_success": parse_success,
        "parse_failed": parse_failed,
        "unique_ccd_count": len(components),
        "source_component_instance_count": source_total,
        "assembly_component_instance_count": assembly_total,
        "conformer_count": conformer_total,
        "ordinary_source_instance_count": ordinary_count,
        "ordinary_assembly_instance_count": ordinary_assembly_count,
        "special_instance_count": special_count,
        "artifact_review_count": artifact_count,
        "polymer_modified_count": polymer_count,
        "unresolved_count": unresolved_count,
        "water_exclusion_count": water_total,
        "covalent_linked_count": covalent_linked_count,
        "metal_inorganic_count": metal_inorganic_count,
        "rdkit_parse_failure_count": rdkit_parse_fail,
        "rdkit_sanitize_failure_count": rdkit_sanitize_fail,
        "ccd_missing_invalid_count": ccd_invalid,
        "historical_crosswalk_rows": cross_count,
        "historical_crosswalk_summary": dict(cross_summary),
        "validation": validation,
    }
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "filter_2_final_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (VALIDATION / "filter_2_release_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (RELEASE / "filter_2_release_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (RELEASE / "filter_2_release_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    interface = {
        "project_name": "Benchmark 1.0",
        "filter_name": "Filter 2 - Ligand Instance Identification and Chemical-Scope Qualification",
        "filter_version": "1.0",
        "input_entry_count": 248037,
        "input_qualified_assembly_count": 360611,
        "source_component_instance_count": source_total,
        "assembly_component_instance_count": assembly_total,
        "unique_component_id_count": len(components),
        "ordinary_source_instance_count": ordinary_count,
        "ordinary_assembly_instance_count": ordinary_assembly_count,
        "special_instance_count": special_count,
        "artifact_review_count": artifact_count,
        "polymer_modified_count": polymer_count,
        "unresolved_count": unresolved_count,
        "ccd_snapshot_version": "Sat, 11 Jul 2026 12:01:19 GMT",
        "artifact_reference_versions": ["refinement_v2_provisional_20260719", "official_reference_unavailable"],
        "rule_version": p.RULE_VERSION,
        "release_creation_time": utc(),
        "release_validation_pass": validation["release_validation_pass"],
    }
    (RELEASE / "filter_2_downstream_interface.json").write_text(json.dumps(interface, indent=2) + "\n")
    release_files = [x for x in RELEASE.iterdir() if x.is_file() and x.name != "SHA256SUMS"]
    (RELEASE / "SHA256SUMS").write_text("".join(f"{p.sha(path)}  {path.name}\n" for path in sorted(release_files)))
    provenance = {"host": platform.node(), "python": sys.version, "start": progress["start"], "end": progress["end"], "workers": progress["workers"], "finalize_mode": "streaming", "release_validation_pass": validation["release_validation_pass"]}
    (OUT / "provenance/filter_2_run_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not validation["release_validation_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
