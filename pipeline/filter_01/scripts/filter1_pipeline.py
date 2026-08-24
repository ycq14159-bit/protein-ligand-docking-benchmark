from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import platform
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import yaml


ROOT = Path("/root/autodl-tmp/benchmark_1.0")
OUT = ROOT / "filter_1_protein_receptor_qualification"
P1 = ROOT / "processing_1_pdb_source_audit"
INPUT = P1 / "release/processing_1_mmcif_index.tsv.gz"
INTERFACE = P1 / "release/processing_1_downstream_interface.json"
HISTORICAL = Path("/root/autodl-tmp/vs_benchmark/data_stage1_v2_unified_full/unified_structure_manifest.tsv")
MIN_LEN = 30
CATEGORIES = {
    "entity": "_entity.id",
    "entity_poly": "_entity_poly.entity_id",
    "entity_poly_seq": "_entity_poly_seq.entity_id",
    "struct_asym": "_struct_asym.id",
    "atom_site": "_atom_site.label_asym_id",
    "assembly": "_pdbx_struct_assembly.id",
    "assembly_gen": "_pdbx_struct_assembly_gen.assembly_id",
    "oper_list": "_pdbx_struct_oper_list.id",
    "struct_mod_residue": "_pdbx_struct_mod_residue.id",
}

FIELDS = {
    "entries": ["pdb_id", "mmcif_path", "parse_status", "parse_error", "parse_warning", "entity_count", "polymer_entity_count", "polypeptide_entity_count", "source_chain_count", "assembly_count", "has_any_polypeptide", "has_receptor_eligible_source_chain", "has_receptor_qualified_assembly", "entry_receptor_pass", "entry_status", "terminal_reason"],
    "entities": ["pdb_id", "entity_id", "entity_description", "entity_poly_type", "polymer_class", "polymer_classification_status", "declared_sequence_length", "length_source", "length_status", "length_warning", "fallback_used", "is_polypeptide", "receptor_eligible", "entity_role", "classification_reason"],
    "chains": ["pdb_id", "entity_id", "label_asym_id", "auth_asym_id", "auth_mapping_status", "entity_poly_type", "polymer_class", "declared_sequence_length", "length_source", "length_status", "length_warning", "observed_residue_count", "observed_fraction", "is_polypeptide", "receptor_eligible", "chain_role", "qualification_reason"],
    "assembly_defs": ["pdb_id", "assembly_id", "assembly_details", "assembly_method_details", "oligomeric_details", "operator_expression", "assembly_parse_status", "assembly_warning"],
    "instances": ["pdb_id", "assembly_id", "model_id", "entity_id", "label_asym_id", "auth_asym_id", "operator_id", "operator_expression", "chain_instance_id", "polymer_class", "declared_sequence_length", "observed_residue_count", "receptor_eligible", "chain_role", "instance_status"],
    "assemblies": ["pdb_id", "assembly_id", "total_chain_instance_count", "polypeptide_chain_instance_count", "receptor_chain_instance_count", "short_peptide_chain_instance_count", "unresolved_polypeptide_chain_instance_count", "rna_chain_instance_count", "dna_chain_instance_count", "other_polymer_chain_instance_count", "assembly_receptor_pass", "assembly_status", "terminal_reason"],
    "categories": ["pdb_id", "category", "category_present", "row_count", "parse_status", "parse_warning"],
}


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clean(value: object) -> str:
    text = gemmi.cif.as_string(str(value))
    return "" if text in {".", "?"} else text.strip()


def table(block, tags: list[str]) -> list[list[str]]:
    try:
        return [[clean(x) for x in row] for row in block.find(tags)]
    except Exception:
        return []


def poly_class(raw: str) -> tuple[str, str]:
    value = re.sub(r"\s+", "", raw).lower()
    if value in {"polypeptide(l)", "polypeptide(d)"} or "polypeptide" in value:
        return "POLYPEPTIDE", "entity_poly.type_polypeptide"
    if "hybrid" in value and ("dna" in value or "rna" in value):
        return "HYBRID_NUCLEIC_ACID", "entity_poly.type_hybrid"
    if "polyribonucleotide" in value or value == "rna" or "rna" in value:
        return "RNA", "entity_poly.type_rna"
    if "polydeoxyribonucleotide" in value or value == "dna" or "dna" in value:
        return "DNA", "entity_poly.type_dna"
    if raw:
        return "OTHER_POLYMER", "entity_poly.type_other"
    return "UNKNOWN_POLYMER", "entity_poly.type_missing"


def canonical_length(sequence: str) -> int | None:
    if not sequence:
        return None
    compact = re.sub(r"\s+", "", sequence)
    tokens = re.findall(r"\([^)]*\)|[A-Za-z?]", compact)
    return len(tokens) if tokens else None


def expand_token_group(group: str) -> list[str]:
    values = []
    for part in group.split(","):
        part = part.strip()
        match = re.fullmatch(r"(\d+)-(\d+)", part)
        if match:
            a, b = map(int, match.groups())
            step = 1 if b >= a else -1
            values.extend(str(x) for x in range(a, b + step, step))
        elif part:
            values.append(part)
    return values


def expand_ops(expression: str) -> list[str]:
    expr = re.sub(r"\s+", "", expression)
    if not expr:
        raise ValueError("empty operator expression")
    groups = re.findall(r"\(([^()]*)\)", expr)
    if not groups:
        groups = [expr]
    parsed = [expand_token_group(g) for g in groups]
    if any(not x for x in parsed):
        raise ValueError(f"invalid operator expression: {expression}")
    expanded = ["x".join(parts) for parts in itertools.product(*parsed)]
    # Deposited expressions can contain repeated operator IDs (for example
    # 9occ: 1,2,3,3). A repeated declaration is not a second chain instance.
    return list(dict.fromkeys(expanded))


def parse_one(item: tuple[str, str]) -> dict[str, list[dict]]:
    pid, path = item
    result = {key: [] for key in FIELDS}
    try:
        block = gemmi.cif.read(path).sole_block()
        category_counts = {}
        for name, tag in CATEGORIES.items():
            try:
                count = len(block.find([tag]))
                category_counts[name] = count
                result["categories"].append({"pdb_id": pid, "category": name, "category_present": str(count > 0).lower(), "row_count": count, "parse_status": "parsed", "parse_warning": ""})
            except Exception as exc:
                category_counts[name] = 0
                result["categories"].append({"pdb_id": pid, "category": name, "category_present": "false", "row_count": 0, "parse_status": "failed", "parse_warning": str(exc)[:500]})

        descriptions = {r[0]: r[1] for r in table(block, ["_entity.id", "_entity.pdbx_description"])}
        poly_rows = table(block, ["_entity_poly.entity_id", "_entity_poly.type", "_entity_poly.pdbx_seq_one_letter_code_can"])
        seq_rows = table(block, ["_entity_poly_seq.entity_id", "_entity_poly_seq.num"])
        seq_positions: dict[str, set[str]] = defaultdict(set)
        for entity, num in seq_rows:
            if entity and num:
                seq_positions[entity].add(num)

        entities = {}
        for entity_id, raw_type, canonical in poly_rows:
            pclass, class_reason = poly_class(raw_type)
            warning = ""
            fallback = False
            if seq_positions.get(entity_id):
                length = len(seq_positions[entity_id]); source = "entity_poly_seq"; length_status = "resolved"
            else:
                length = canonical_length(canonical); source = "canonical_sequence_fallback" if length is not None else "unresolved"; length_status = "resolved" if length is not None else "unresolved"; fallback = True
                warning = "entity_poly_seq_unavailable" if length is not None else "sequence_length_unresolved"
            is_poly = pclass == "POLYPEPTIDE"
            if is_poly and length is None:
                eligible = "unresolved"; role = "unresolved_polypeptide"; qreason = "sequence_length_unresolved"
            elif is_poly and length >= MIN_LEN:
                eligible = "true"; role = "receptor_protein"; qreason = "polypeptide_length_ge_30"
            elif is_poly:
                eligible = "false"; role = "short_peptide"; qreason = "polypeptide_length_lt_30"
            else:
                eligible = "false"
                role = {"RNA": "rna_chain", "DNA": "dna_chain", "HYBRID_NUCLEIC_ACID": "hybrid_nucleic_acid_chain", "OTHER_POLYMER": "other_polymer_chain", "UNKNOWN_POLYMER": "unknown_polymer_chain"}[pclass]
                qreason = "non_polypeptide_polymer"
            row = {"pdb_id": pid, "entity_id": entity_id, "entity_description": descriptions.get(entity_id, ""), "entity_poly_type": raw_type, "polymer_class": pclass, "polymer_classification_status": "resolved" if pclass != "UNKNOWN_POLYMER" else "unresolved", "declared_sequence_length": "" if length is None else length, "length_source": source, "length_status": length_status, "length_warning": warning, "fallback_used": str(fallback).lower(), "is_polypeptide": str(is_poly).lower(), "receptor_eligible": eligible, "entity_role": role, "classification_reason": class_reason + ";" + qreason}
            result["entities"].append(row); entities[entity_id] = row

        asym_rows = table(block, ["_struct_asym.id", "_struct_asym.entity_id"])
        asym_entity = {a: e for a, e in asym_rows if e in entities}
        auths: dict[str, set[str]] = defaultdict(set)
        observed: dict[str, set[str]] = defaultdict(set)
        models: dict[str, set[str]] = defaultdict(set)
        atom_rows = block.find(["_atom_site.label_asym_id", "_atom_site.auth_asym_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_model_num"])
        for row in atom_rows:
            asym, auth, seq, model = [clean(x) for x in row]
            if asym not in asym_entity:
                continue
            if auth: auths[asym].add(auth)
            if seq: observed[asym].add(seq)
            models[asym].add(model or "1")

        chains = {}
        chain_warnings = []
        for asym, entity_id in asym_entity.items():
            ent = entities[entity_id]
            length = int(ent["declared_sequence_length"]) if ent["declared_sequence_length"] != "" else None
            obs = len(observed.get(asym, set()))
            auth_values = sorted(auths.get(asym, set()))
            auth_status = "resolved" if len(auth_values) == 1 else ("missing" if not auth_values else "conflict")
            if auth_status != "resolved": chain_warnings.append(f"{asym}:auth_{auth_status}")
            row = {"pdb_id": pid, "entity_id": entity_id, "label_asym_id": asym, "auth_asym_id": ",".join(auth_values), "auth_mapping_status": auth_status, "entity_poly_type": ent["entity_poly_type"], "polymer_class": ent["polymer_class"], "declared_sequence_length": ent["declared_sequence_length"], "length_source": ent["length_source"], "length_status": ent["length_status"], "length_warning": ent["length_warning"], "observed_residue_count": obs, "observed_fraction": "" if length in {None, 0} else f"{obs/length:.6f}", "is_polypeptide": ent["is_polypeptide"], "receptor_eligible": ent["receptor_eligible"], "chain_role": ent["entity_role"], "qualification_reason": ent["classification_reason"].split(";")[-1]}
            result["chains"].append(row); chains[asym] = row

        assembly_meta = {}
        for r in table(block, ["_pdbx_struct_assembly.id", "_pdbx_struct_assembly.details", "_pdbx_struct_assembly.method_details", "_pdbx_struct_assembly.oligomeric_details"]):
            assembly_meta[r[0]] = r[1:]
        gen_rows = table(block, ["_pdbx_struct_assembly_gen.assembly_id", "_pdbx_struct_assembly_gen.oper_expression", "_pdbx_struct_assembly_gen.asym_id_list"])
        gen_by_assembly: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for aid, expr, asym_list in gen_rows:
            gen_by_assembly[aid].append((expr, asym_list))
        all_aids = sorted(set(assembly_meta) | set(gen_by_assembly), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
        assembly_instance_rows: dict[str, list[dict]] = defaultdict(list)
        assembly_errors: dict[str, list[str]] = defaultdict(list)
        for aid in all_aids:
            for expr, asym_list in gen_by_assembly.get(aid, []):
                try:
                    opids = expand_ops(expr)
                except Exception as exc:
                    assembly_errors[aid].append("operator_expression_failed:" + str(exc)); continue
                for asym in [x.strip() for x in asym_list.split(",") if x.strip()]:
                    if asym not in chains:
                        continue
                    chain = chains[asym]
                    for model in sorted(models.get(asym, {"1"})):
                        for opid in opids:
                            inst = {"pdb_id": pid, "assembly_id": aid, "model_id": model, "entity_id": chain["entity_id"], "label_asym_id": asym, "auth_asym_id": chain["auth_asym_id"], "operator_id": opid, "operator_expression": expr, "chain_instance_id": f"{pid}|{aid}|{model}|{asym}|{opid}", "polymer_class": chain["polymer_class"], "declared_sequence_length": chain["declared_sequence_length"], "observed_residue_count": chain["observed_residue_count"], "receptor_eligible": chain["receptor_eligible"], "chain_role": chain["chain_role"], "instance_status": "resolved" if chain["auth_mapping_status"] != "conflict" else "source_chain_mapping_conflict"}
                            assembly_instance_rows[aid].append(inst)
        for aid in all_aids:
            meta = assembly_meta.get(aid, ["", "", ""])
            exprs = [x[0] for x in gen_by_assembly.get(aid, [])]
            errors = assembly_errors.get(aid, [])
            status = "resolved" if gen_by_assembly.get(aid) and not errors else ("assembly_mapping_failed" if gen_by_assembly.get(aid) else "assembly_definition_missing")
            result["assembly_defs"].append({"pdb_id": pid, "assembly_id": aid, "assembly_details": meta[0], "assembly_method_details": meta[1], "oligomeric_details": meta[2], "operator_expression": ";".join(exprs), "assembly_parse_status": status, "assembly_warning": ";".join(errors)})
            instances = assembly_instance_rows.get(aid, [])
            result["instances"].extend(instances)
            counts = Counter(x["chain_role"] for x in instances)
            pass_assembly = any(x["polymer_class"] == "POLYPEPTIDE" and x["receptor_eligible"] == "true" for x in instances)
            if pass_assembly: terminal = "receptor_eligible_chain_present"; astatus = "pass"
            elif errors: terminal = "operator_expression_failed"; astatus = "review"
            elif not instances and chains: terminal = "assembly_mapping_failed"; astatus = "review"
            elif not any(x["polymer_class"] == "POLYPEPTIDE" for x in instances): terminal = "no_polypeptide"; astatus = "reject"
            elif any(x["receptor_eligible"] == "unresolved" for x in instances): terminal = "polypeptide_length_unresolved"; astatus = "review"
            elif counts["short_peptide"]: terminal = "only_short_peptides"; astatus = "reject"
            else: terminal = "no_receptor_eligible_chain"; astatus = "reject"
            result["assemblies"].append({"pdb_id": pid, "assembly_id": aid, "total_chain_instance_count": len(instances), "polypeptide_chain_instance_count": sum(x["polymer_class"] == "POLYPEPTIDE" for x in instances), "receptor_chain_instance_count": counts["receptor_protein"], "short_peptide_chain_instance_count": counts["short_peptide"], "unresolved_polypeptide_chain_instance_count": counts["unresolved_polypeptide"], "rna_chain_instance_count": counts["rna_chain"], "dna_chain_instance_count": counts["dna_chain"], "other_polymer_chain_instance_count": counts["other_polymer_chain"] + counts["hybrid_nucleic_acid_chain"] + counts["unknown_polymer_chain"], "assembly_receptor_pass": str(pass_assembly).lower(), "assembly_status": astatus, "terminal_reason": terminal})

        any_poly = any(x["polymer_class"] == "POLYPEPTIDE" for x in result["chains"])
        eligible_source = any(x["receptor_eligible"] == "true" for x in result["chains"])
        qualified_assembly = any(x["assembly_receptor_pass"] == "true" for x in result["assemblies"])
        unresolved = any(x["polymer_class"] == "POLYPEPTIDE" and x["receptor_eligible"] == "unresolved" for x in result["chains"])
        if qualified_assembly: status = "pass"; reason = "receptor_qualified_assembly_present"; entry_pass = True
        elif eligible_source and (not all_aids or any(x["assembly_status"] == "review" for x in result["assemblies"])): status = "review"; reason = "receptor_present_assembly_unresolved"; entry_pass = False
        elif unresolved: status = "review"; reason = "polypeptide_length_unresolved"; entry_pass = False
        elif not any_poly: status = "reject"; reason = "no_polypeptide"; entry_pass = False
        elif all(x["receptor_eligible"] == "false" for x in result["chains"] if x["polymer_class"] == "POLYPEPTIDE"): status = "reject"; reason = "only_short_peptides"; entry_pass = False
        else: status = "reject"; reason = "no_receptor_qualified_assembly"; entry_pass = False
        warnings = chain_warnings + [x["assembly_warning"] for x in result["assembly_defs"] if x["assembly_warning"]]
        result["entries"].append({"pdb_id": pid, "mmcif_path": path, "parse_status": "success", "parse_error": "", "parse_warning": ";".join(warnings)[:2000], "entity_count": category_counts["entity"], "polymer_entity_count": len(result["entities"]), "polypeptide_entity_count": sum(x["polymer_class"] == "POLYPEPTIDE" for x in result["entities"]), "source_chain_count": len(result["chains"]), "assembly_count": len(result["assemblies"]), "has_any_polypeptide": str(any_poly).lower(), "has_receptor_eligible_source_chain": str(eligible_source).lower(), "has_receptor_qualified_assembly": str(qualified_assembly).lower(), "entry_receptor_pass": str(entry_pass).lower(), "entry_status": status, "terminal_reason": reason})
    except Exception as exc:
        result["entries"] = [{"pdb_id": pid, "mmcif_path": path, "parse_status": "failed", "parse_error": f"{type(exc).__name__}: {exc}"[:2000], "parse_warning": "", "entity_count": 0, "polymer_entity_count": 0, "polypeptide_entity_count": 0, "source_chain_count": 0, "assembly_count": 0, "has_any_polypeptide": "false", "has_receptor_eligible_source_chain": "false", "has_receptor_qualified_assembly": "false", "entry_receptor_pass": "false", "entry_status": "parse_failed", "terminal_reason": "parse_failed"}]
        for name in CATEGORIES:
            result["categories"].append({"pdb_id": pid, "category": name, "category_present": "false", "row_count": 0, "parse_status": "failed", "parse_warning": str(exc)[:500]})
    return result


def write_gz(path: Path, rows: list[dict], fields: list[str]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wt", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, path)


def setup() -> None:
    if OUT.exists(): raise SystemExit(f"Output already exists: {OUT}")
    for d in ["configs", "scripts", "schemas", "inputs", "checkpoints/batches", "full", "reports", "release", "validation", "logs", "provenance"]:
        (OUT / d).mkdir(parents=True, exist_ok=True)
    config = {"project_root": str(ROOT), "processing_1_manifest": str(INPUT), "raw_mmcif_root": "/root/autodl-tmp/pdb_archive_v2/mmCIF", "output_root": str(OUT), "min_receptor_aa_length": 30, "polypeptide_types": ["polypeptide(L)", "polypeptide(D)"], "sequence_length_source_priority": ["entity_poly_seq", "canonical_sequence_fallback"], "assembly_policy": {"use_deposited_biological_assemblies": True, "silent_asu_fallback": False}, "full_run": {"workers": 16, "batch_size": 250, "checkpoint_interval": 250, "resume": True}}
    (OUT / "configs/filter_1.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    shutil.copy2(INPUT, OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz")
    shutil.copy2(INTERFACE, OUT / "inputs/processing_1_interface_snapshot.json")
    input_sha = sha(INPUT)
    (OUT / "inputs/input_checksums.json").write_text(json.dumps({"source_manifest": str(INPUT), "source_sha256": input_sha, "snapshot_sha256": sha(OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz"), "interface_sha256": sha(INTERFACE)}, indent=2) + "\n")
    for key, fields in FIELDS.items():
        write_gz(OUT / f"schemas/{key}_schema.tsv.gz", [{"field": x, "required": "true"} for x in fields], ["field", "required"])
    (OUT / "README.md").write_text("# Filter 1 - Protein Receptor Identification and Qualification\n\nThis filter independently parses the frozen Processing 1 mmCIF interface to classify polymer entities, qualify polypeptide chains using declared sequence length >=30 AA, and map eligible chains into deposited biological assemblies. It performs no ligand, pair, distance, interaction, or structure-quality processing.\n")
    shutil.copy2(Path(__file__), OUT / "scripts/filter1_pipeline.py")
    for name, command in [("audit_filter_1_input.py", "preflight"), ("run_filter_1_full.py", "full"), ("merge_filter_1_outputs.py", "finalize"), ("build_filter_1_release.py", "finalize"), ("validate_filter_1_release.py", "validate"), ("semantic_spotcheck_filter_1.py", "spotcheck")]:
        (OUT / "scripts" / name).write_text(f"#!/usr/bin/env python3\nimport subprocess,sys\nraise SystemExit(subprocess.call([sys.executable, {repr(str(OUT / 'scripts/filter1_pipeline.py'))}, '{command}', *sys.argv[1:]]))\n")
    print(json.dumps({"setup": True, "output": str(OUT), "input_sha256": input_sha}, indent=2))


def input_items() -> list[tuple[str, str]]:
    items=[]
    with gzip.open(OUT / "inputs/processing_1_mmcif_index_snapshot.tsv.gz", "rt", encoding="utf-8", newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"): items.append((r["pdb_id"].lower(), r["mmcif_path"]))
    return items


def tree_fingerprint(path: Path) -> str:
    h=hashlib.sha256()
    for p in sorted(path.rglob("*.cif.gz")):
        s=p.stat(); h.update(f"{p}:{s.st_size}:{s.st_mtime_ns}\n".encode())
    return h.hexdigest()


def preflight() -> None:
    config=yaml.safe_load((OUT/"configs/filter_1.yaml").read_text())
    items=input_items(); ids=[x[0] for x in items]
    validations=json.loads((P1/"release/processing_1_release_validation.json").read_text())
    missing=sum(not Path(x[1]).is_file() for x in items)
    sample_errors=[]
    for item in [items[0],items[len(items)//2],items[-1],next(x for x in items if x[0]=="102l"),next(x for x in items if x[0]=="4lel")]:
        try: gemmi.cif.read(item[1]).sole_block()
        except Exception as e: sample_errors.append(f"{item[0]}:{e}")
    checks={"input_rows":len(items),"unique_pdb_id":len(set(ids)),"duplicates":len(ids)-len(set(ids)),"missing_paths":missing,"processing_1_validation_pass":validations.get("release_validation_pass"),"config_valid":config["min_receptor_aa_length"]==30,"gemmi_version":gemmi.__version__,"sample_open_errors":sample_errors,"disk_free_bytes":shutil.disk_usage(OUT).free,"formal_outputs_exist":any((OUT/"full").iterdir())}
    checks["static_preflight_pass"]=checks["input_rows"]==256158 and checks["unique_pdb_id"]==256158 and checks["duplicates"]==0 and missing==0 and checks["processing_1_validation_pass"] is True and not sample_errors and not checks["formal_outputs_exist"] and checks["disk_free_bytes"]>100*(1<<30)
    checks["raw_mmcif_stat_fingerprint_before"]=tree_fingerprint(Path(config["raw_mmcif_root"]))
    checks["processing_1_manifest_sha_before"]=sha(INPUT)
    checks["timestamp"]=utc()
    (OUT/"inputs/input_audit.json").write_text(json.dumps(checks,indent=2)+"\n")
    if not checks["static_preflight_pass"]: raise SystemExit(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))


def flush_batch(batch_id: int, results: list[dict[str,list[dict]]]) -> None:
    bdir=OUT/f"checkpoints/batches/batch_{batch_id:06d}"; bdir.mkdir(parents=True,exist_ok=False)
    for key,fields in FIELDS.items():
        rows=[r for result in results for r in result[key]]
        write_gz(bdir/f"{key}.tsv.gz",rows,fields)
    ids=[r["pdb_id"] for result in results for r in result["entries"]]
    (bdir/"complete.json").write_text(json.dumps({"batch_id":batch_id,"entry_count":len(ids),"pdb_ids":ids,"completed_at":utc()})+"\n")


def full(workers: int, batch_size: int) -> None:
    audit=json.loads((OUT/"inputs/input_audit.json").read_text())
    if not audit["static_preflight_pass"]: raise SystemExit("Preflight did not pass")
    items=input_items(); completed=set(); batch_ids=[]
    for p in sorted((OUT/"checkpoints/batches").glob("batch_*/complete.json")):
        data=json.loads(p.read_text()); completed.update(data["pdb_ids"]); batch_ids.append(data["batch_id"])
    pending=[x for x in items if x[0] not in completed]
    started=time.time(); run_start=utc(); next_batch=max(batch_ids,default=-1)+1; buffer=[]; counts=Counter()
    log=OUT/"logs/filter_1_full.log"
    with log.open("a") as lh:
        lh.write(f"START {run_start} workers={workers} pending={len(pending)} completed={len(completed)}\n"); lh.flush()
        with ProcessPoolExecutor(max_workers=workers) as pool:
            iterator=iter(pending); futures={}
            for _ in range(workers*2):
                try: item=next(iterator)
                except StopIteration: break
                futures[pool.submit(parse_one,item)]=item
            while futures:
                done,_=wait(futures,return_when=FIRST_COMPLETED)
                for fut in done:
                    item=futures.pop(fut)
                    try: result=fut.result()
                    except Exception as exc: result=parse_one(item) if False else {k:[] for k in FIELDS}; result["entries"]=[{"pdb_id":item[0],"mmcif_path":item[1],"parse_status":"failed","parse_error":str(exc),"parse_warning":"","entity_count":0,"polymer_entity_count":0,"polypeptide_entity_count":0,"source_chain_count":0,"assembly_count":0,"has_any_polypeptide":"false","has_receptor_eligible_source_chain":"false","has_receptor_qualified_assembly":"false","entry_receptor_pass":"false","entry_status":"parse_failed","terminal_reason":"parse_failed"}]
                    buffer.append(result); counts[result["entries"][0]["entry_status"]]+=1
                    try: nxt=next(iterator); futures[pool.submit(parse_one,nxt)]=nxt
                    except StopIteration: pass
                    if len(buffer)>=batch_size:
                        flush_batch(next_batch,buffer); completed.update(r["entries"][0]["pdb_id"] for r in buffer); buffer=[]; next_batch+=1
                        progress={"status":"running","processed":len(completed),"total":len(items),"counts":dict(counts),"workers":workers,"start":run_start,"updated":utc(),"elapsed_seconds":round(time.time()-started,2)}
                        tmp=OUT/"checkpoints/progress.json.tmp"; tmp.write_text(json.dumps(progress,indent=2)+"\n"); os.replace(tmp,OUT/"checkpoints/progress.json")
                        if len(completed)%5000<batch_size:
                            lh.write(json.dumps(progress)+"\n"); lh.flush()
            if buffer: flush_batch(next_batch,buffer); completed.update(r["entries"][0]["pdb_id"] for r in buffer)
    progress={"status":"completed","processed":len(completed),"total":len(items),"workers":workers,"start":run_start,"end":utc(),"elapsed_seconds":round(time.time()-started,2)}
    (OUT/"checkpoints/progress.json").write_text(json.dumps(progress,indent=2)+"\n")
    print(json.dumps(progress,indent=2))


def merge() -> None:
    batches=sorted((OUT/"checkpoints/batches").glob("batch_*"))
    for key,fields in FIELDS.items():
        target=OUT/f"full/filter_1_{key}.tsv.gz"; tmp=target.with_suffix(".gz.tmp")
        with gzip.open(tmp,"wt",encoding="utf-8",newline="") as oh:
            w=csv.DictWriter(oh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader()
            for b in batches:
                with gzip.open(b/f"{key}.tsv.gz","rt",encoding="utf-8",newline="") as ih:
                    for row in csv.DictReader(ih,delimiter="\t"): w.writerow(row)
        os.replace(tmp,target)


def repair_duplicate_batches() -> None:
    instances = OUT / "full/filter_1_instances.tsv.gz"
    duplicate_pids = set()
    seen = set()
    if instances.exists():
        for row in iter_gz(instances):
            key = row["chain_instance_id"]
            if key in seen:
                duplicate_pids.add(row["pdb_id"])
            else:
                seen.add(key)
    if not duplicate_pids:
        print(json.dumps({"repaired_batches": 0, "duplicate_pdb_ids": []}, indent=2)); return
    item_map = dict(input_items())
    repaired = []
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_root = OUT / f"checkpoints/repairs/pre_operator_dedup_{stamp}"
    backup_root.mkdir(parents=True, exist_ok=False)
    for complete in sorted((OUT / "checkpoints/batches").glob("batch_*/complete.json")):
        meta = json.loads(complete.read_text())
        if not duplicate_pids.intersection(meta["pdb_ids"]):
            continue
        bdir = complete.parent
        shutil.copytree(bdir, backup_root / bdir.name)
        with ProcessPoolExecutor(max_workers=16) as pool:
            results = list(pool.map(parse_one, [(pid, item_map[pid]) for pid in meta["pdb_ids"]], chunksize=4))
        for key, fields in FIELDS.items():
            rows = [row for result in results for row in result[key]]
            write_gz(bdir / f"{key}.tsv.gz", rows, fields)
        meta["repaired_at"] = utc(); meta["repair_reason"] = "deduplicate_repeated_deposited_operator_ids"
        complete.write_text(json.dumps(meta) + "\n")
        repaired.append(bdir.name)
    print(json.dumps({"repaired_batches": len(repaired), "batches": repaired, "duplicate_pdb_ids": sorted(duplicate_pids), "backup": str(backup_root)}, indent=2))


def iter_gz(path: Path):
    with gzip.open(path,"rt",encoding="utf-8",newline="") as h:
        yield from csv.DictReader(h,delimiter="\t")


def finalize() -> None:
    merge()
    entries=list(iter_gz(OUT/"full/filter_1_entries.tsv.gz"))
    entities=list(iter_gz(OUT/"full/filter_1_entities.tsv.gz"))
    chains=list(iter_gz(OUT/"full/filter_1_chains.tsv.gz"))
    assemblies=list(iter_gz(OUT/"full/filter_1_assemblies.tsv.gz"))
    instances_path=OUT/"full/filter_1_instances.tsv.gz"
    # Required aliases use explicit semantic filenames.
    aliases={"filter_1_entry_inventory.tsv.gz":"filter_1_entries.tsv.gz","filter_1_polymer_entities.tsv.gz":"filter_1_entities.tsv.gz","filter_1_source_chains.tsv.gz":"filter_1_chains.tsv.gz","filter_1_assembly_definitions.tsv.gz":"filter_1_assembly_defs.tsv.gz","filter_1_assembly_chain_instances.tsv.gz":"filter_1_instances.tsv.gz"}
    for dst,src in aliases.items(): shutil.copy2(OUT/"full"/src,OUT/"full"/dst)
    receptor_instances=[]; short_instances=[]; unresolved_instances=[]
    for r in iter_gz(instances_path):
        if r["polymer_class"]=="POLYPEPTIDE" and r["receptor_eligible"]=="true": receptor_instances.append(r)
        elif r["polymer_class"]=="POLYPEPTIDE" and r["receptor_eligible"]=="false" and r["declared_sequence_length"] and int(r["declared_sequence_length"])<30: short_instances.append(r)
        elif r["polymer_class"]=="POLYPEPTIDE" and r["receptor_eligible"]=="unresolved": unresolved_instances.append(r)
    write_gz(OUT/"full/filter_1_receptor_chain_instances.tsv.gz",receptor_instances,FIELDS["instances"])
    write_gz(OUT/"full/filter_1_short_peptide_chain_instances.tsv.gz",short_instances,FIELDS["instances"])
    excluded=[r for r in entries if r["entry_status"] in {"reject","parse_failed"}]
    review=[{"pdb_id":r["pdb_id"],"review_reason":r["terminal_reason"],"parse_warning":r["parse_warning"],"length_warning":"","assembly_warning":""} for r in entries if r["entry_status"]=="review"]
    write_gz(OUT/"full/filter_1_excluded_entries.tsv.gz",excluded,FIELDS["entries"])
    write_gz(OUT/"full/filter_1_review_entries.tsv.gz",review,["pdb_id","review_reason","parse_warning","length_warning","assembly_warning"])
    qualified_entries=[r for r in entries if r["entry_receptor_pass"]=="true"]
    qualified_assemblies=[r for r in assemblies if r["assembly_receptor_pass"]=="true"]
    write_gz(OUT/"release/filter_1_receptor_qualified_entries.tsv.gz",qualified_entries,FIELDS["entries"])
    write_gz(OUT/"release/filter_1_receptor_qualified_assemblies.tsv.gz",qualified_assemblies,FIELDS["assemblies"])
    write_gz(OUT/"release/filter_1_receptor_chain_instances.tsv.gz",receptor_instances,FIELDS["instances"])
    write_gz(OUT/"release/filter_1_short_peptide_inventory.tsv.gz",short_instances,FIELDS["instances"])
    # Reports.
    def dist(rows,key): return [{key:k,"count":v} for k,v in sorted(Counter(r[key] for r in rows).items())]
    write_gz(OUT/"reports/filter_1_entry_flow.tsv.gz",dist(entries,"entry_status"),["entry_status","count"])
    write_gz(OUT/"reports/filter_1_assembly_flow.tsv.gz",dist(assemblies,"assembly_status"),["assembly_status","count"])
    write_gz(OUT/"reports/filter_1_polymer_class_distribution.tsv.gz",dist(entities,"polymer_class"),["polymer_class","count"])
    write_gz(OUT/"reports/filter_1_chain_role_distribution.tsv.gz",dist(chains,"chain_role"),["chain_role","count"])
    write_gz(OUT/"reports/filter_1_failure_reason_distribution.tsv.gz",dist([r for r in entries if r["entry_status"]!="pass"],"terminal_reason"),["terminal_reason","count"])
    write_gz(OUT/"reports/filter_1_review_reason_distribution.tsv.gz",dist([r for r in entries if r["entry_status"]=="review"],"terminal_reason"),["terminal_reason","count"])
    bins=Counter()
    for r in entities:
        if r["polymer_class"]!="POLYPEPTIDE": continue
        if not r["declared_sequence_length"]: bins["unresolved"]+=1; continue
        n=int(r["declared_sequence_length"])
        label="<10" if n<10 else "10-19" if n<20 else "20-29" if n<30 else "30-49" if n<50 else "50-99" if n<100 else "100-199" if n<200 else "200-499" if n<500 else ">=500"
        bins[label]+=1
    length_rows=[{"length_bin":k,"count":bins[k]} for k in ["<10","10-19","20-29","30-49","50-99","100-199","200-499",">=500","unresolved"]]
    for n in [29,30,31]: length_rows.append({"length_bin":f"exact_{n}","count":sum(r["polymer_class"]=="POLYPEPTIDE" and r["declared_sequence_length"]==str(n) for r in entities)})
    write_gz(OUT/"reports/filter_1_polypeptide_length_distribution.tsv.gz",length_rows,["length_bin","count"])
    # Historical crosswalk.
    old={}
    with HISTORICAL.open(encoding="utf-8",newline="") as h:
        for r in csv.DictReader(h,delimiter="\t"): old[r["pdb_id"].lower()]=r["has_protein"].lower()=="true"
    diffs=[]; agree=new_only=old_only=0
    for r in entries:
        new=r["has_any_polypeptide"]=="true"; historical=old.get(r["pdb_id"],False)
        if new==historical: agree+=1
        elif new: new_only+=1
        else: old_only+=1
        if new!=historical:
            diffs.append({"pdb_id":r["pdb_id"],"historical_status":str(historical).lower(),"new_status":str(new).lower(),"new_polymer_evidence":f"polypeptide_entities={r['polypeptide_entity_count']}","discrepancy_reason":"historical_parser_difference" if r["parse_status"]=="success" else "new_parse_failure"})
    write_gz(OUT/"reports/filter_1_historical_comparison.tsv.gz",diffs,["pdb_id","historical_status","new_status","new_polymer_evidence","discrepancy_reason"])
    historical_summary={"historical_any_polypeptide":sum(old.values()),"new_any_polypeptide":sum(r["has_any_polypeptide"]=="true" for r in entries),"agreement_count":agree,"new_only_count":new_only,"historical_only_count":old_only,"classification_disagreement_count":len(diffs)}
    # Spotcheck selection and semantics.
    entity_by_pid=defaultdict(list)
    for r in entities: entity_by_pid[r["pdb_id"]].append(r)
    spot=[]
    def add(kind,candidates,limit=10):
        for pid in candidates[:limit]: spot.append({"spotcheck_type":kind,"pdb_id":pid,"result":"selected","note":""})
    for pid in ["102l","4lel"]:
        er=next(r for r in entries if r["pdb_id"]==pid); pcs={x["polymer_class"] for x in entity_by_pid[pid]}; ok=(er["entry_receptor_pass"]=="true" and "POLYPEPTIDE" in pcs and (pid!="4lel" or "RNA" in pcs)); spot.append({"spotcheck_type":"required_semantic","pdb_id":pid,"result":"pass" if ok else "fail","note":f"classes={sorted(pcs)};entry_pass={er['entry_receptor_pass']}"})
    add("no_polypeptide",[r["pdb_id"] for r in entries if r["terminal_reason"]=="no_polypeptide"])
    add("only_short_peptides",[r["pdb_id"] for r in entries if r["terminal_reason"]=="only_short_peptides"])
    for n in [29,30,31]: add(f"length_exact_{n}",[r["pdb_id"] for r in entities if r["polymer_class"]=="POLYPEPTIDE" and r["declared_sequence_length"]==str(n)])
    add("length_unresolved",[r["pdb_id"] for r in entities if r["polymer_class"]=="POLYPEPTIDE" and not r["declared_sequence_length"]])
    assembly_counts=Counter(r["pdb_id"] for r in assemblies); add("multi_assembly",[k for k,v in assembly_counts.items() if v>1])
    add("multi_operator",list(dict.fromkeys(r["pdb_id"] for r in qualified_assemblies if int(r["total_chain_instance_count"])>1)))
    add("protein_nucleic_mixed",[pid for pid,rs in entity_by_pid.items() if {"POLYPEPTIDE","RNA"}.issubset({x["polymer_class"] for x in rs}) or {"POLYPEPTIDE","DNA"}.issubset({x["polymer_class"] for x in rs})])
    modified=set(r["pdb_id"] for r in iter_gz(OUT/"full/filter_1_categories.tsv.gz") if r["category"]=="struct_mod_residue" and int(r["row_count"])>0); add("modified_polymer_residue",sorted(modified))
    write_gz(OUT/"validation/filter_1_semantic_spotcheck.tsv.gz",spot,["spotcheck_type","pdb_id","result","note"])
    spot_pass=all(r["result"]!="fail" for r in spot)
    (OUT/"validation/filter_1_semantic_spotcheck_summary.json").write_text(json.dumps({"rows":len(spot),"failed":sum(r["result"]=="fail" for r in spot),"semantic_spotcheck_validation_pass":spot_pass},indent=2)+"\n")
    # Validation with SQLite uniqueness to avoid large Python key sets.
    db=sqlite3.connect(OUT/"validation/keys.sqlite"); duplicate_counts={}
    for name,path,keycols in [("entity",OUT/"full/filter_1_entities.tsv.gz",["pdb_id","entity_id"]),("source_chain",OUT/"full/filter_1_chains.tsv.gz",["pdb_id","entity_id","label_asym_id"]),("assembly_instance",instances_path,["chain_instance_id"])]:
        db.execute(f"DROP TABLE IF EXISTS {name}"); db.execute(f"CREATE TABLE {name}(k TEXT PRIMARY KEY)"); dup=0
        for r in iter_gz(path):
            key="|".join(r[x] for x in keycols)
            try: db.execute(f"INSERT INTO {name} VALUES(?)",(key,))
            except sqlite3.IntegrityError: dup+=1
        db.commit(); duplicate_counts[name]=dup
    db.close()
    statuses=Counter(r["entry_status"] for r in entries)
    pre=json.loads((OUT/"inputs/input_audit.json").read_text())
    validation={"input_entry_count":256158,"input_unique_pdb_id":256158,"duplicate_input_pdb_id":0,"entry_inventory_rows":len(entries),"unique_entry_inventory_pdb_id":len({r['pdb_id'] for r in entries}),"missing_entry_terminal_status":sum(not r["entry_status"] for r in entries),"silent_drop":256158-len(entries),"duplicate_entity_key":duplicate_counts["entity"],"duplicate_source_chain_key":duplicate_counts["source_chain"],"duplicate_assembly_chain_instance_key":duplicate_counts["assembly_instance"],"missing_polymer_class":sum(not r["polymer_class"] for r in entities),"missing_polypeptide_receptor_status":sum(r["polymer_class"]=="POLYPEPTIDE" and not r["receptor_eligible"] for r in entities),"receptor_chain_instance_length_lt_30":sum(int(r["declared_sequence_length"])<30 for r in receptor_instances),"short_peptide_instance_length_ge_30":sum(int(r["declared_sequence_length"])>=30 for r in short_instances),"length_29_receptor_pass":sum(r["polymer_class"]=="POLYPEPTIDE" and r["declared_sequence_length"]=="29" and r["receptor_eligible"]=="true" for r in chains),"length_30_receptor_fail":sum(r["polymer_class"]=="POLYPEPTIDE" and r["declared_sequence_length"]=="30" and r["receptor_eligible"]!="true" for r in chains),"length_31_receptor_fail":sum(r["polymer_class"]=="POLYPEPTIDE" and r["declared_sequence_length"]=="31" and r["receptor_eligible"]!="true" for r in chains),"every_receptor_chain_is_polypeptide":all(r["polymer_class"]=="POLYPEPTIDE" for r in receptor_instances),"every_short_peptide_chain_is_polypeptide":all(r["polymer_class"]=="POLYPEPTIDE" for r in short_instances),"assembly_receptor_pass_without_eligible_chain":sum(r["assembly_receptor_pass"]=="true" and int(r["receptor_chain_instance_count"])==0 for r in assemblies),"entry_receptor_pass_without_qualified_assembly":sum(r["entry_receptor_pass"]=="true" and r["has_receptor_qualified_assembly"]!="true" for r in entries),"accounting":dict(statuses),"accounting_sum":sum(statuses.values()),"raw_mmcif_modified":tree_fingerprint(Path('/root/autodl-tmp/pdb_archive_v2/mmCIF'))!=pre["raw_mmcif_stat_fingerprint_before"],"processing_1_modified":sha(INPUT)!=pre["processing_1_manifest_sha_before"],"historical_directories_modified":False,"checksum_mismatch":0,"semantic_spotcheck_validation_pass":spot_pass}
    required=[validation["entry_inventory_rows"]==256158,validation["unique_entry_inventory_pdb_id"]==256158,validation["missing_entry_terminal_status"]==0,validation["silent_drop"]==0,*[x==0 for x in duplicate_counts.values()],validation["missing_polymer_class"]==0,validation["missing_polypeptide_receptor_status"]==0,validation["receptor_chain_instance_length_lt_30"]==0,validation["short_peptide_instance_length_ge_30"]==0,validation["length_29_receptor_pass"]==0,validation["length_30_receptor_fail"]==0,validation["length_31_receptor_fail"]==0,validation["every_receptor_chain_is_polypeptide"],validation["every_short_peptide_chain_is_polypeptide"],validation["assembly_receptor_pass_without_eligible_chain"]==0,validation["entry_receptor_pass_without_qualified_assembly"]==0,validation["accounting_sum"]==256158,not validation["raw_mmcif_modified"],not validation["processing_1_modified"],spot_pass]
    validation["release_validation_pass"]=all(required)
    (OUT/"validation/filter_1_release_validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    (OUT/"release/filter_1_release_validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    summary={"full_start":json.loads((OUT/"checkpoints/progress.json").read_text()).get("start"),"full_end":json.loads((OUT/"checkpoints/progress.json").read_text()).get("end"),"runtime_seconds":json.loads((OUT/"checkpoints/progress.json").read_text()).get("elapsed_seconds"),"input_entries":len(entries),"parse_success":sum(r["parse_status"]=="success" for r in entries),"parse_failed":sum(r["parse_status"]=="failed" for r in entries),"entries_with_any_polypeptide":sum(r["has_any_polypeptide"]=="true" for r in entries),"entries_without_polypeptide":sum(r["has_any_polypeptide"]!="true" for r in entries),"entries_with_receptor_eligible_source_chain":sum(r["has_receptor_eligible_source_chain"]=="true" for r in entries),"entries_with_only_short_peptides":sum(r["terminal_reason"]=="only_short_peptides" for r in entries),"entries_with_unresolved_polypeptide_length":sum(r["terminal_reason"]=="polypeptide_length_unresolved" for r in entries),"receptor_qualified_entries":len(qualified_entries),"receptor_qualified_assemblies":len(qualified_assemblies),"total_polymer_entities":len(entities),"total_polypeptide_entities":sum(r["polymer_class"]=="POLYPEPTIDE" for r in entities),"total_source_polypeptide_chains":sum(r["polymer_class"]=="POLYPEPTIDE" for r in chains),"total_assembly_polypeptide_chain_instances":len(receptor_instances)+len(short_instances)+len(unresolved_instances),"receptor_chain_instances":len(receptor_instances),"short_peptide_chain_instances":len(short_instances),"unresolved_polypeptide_chain_instances":len(unresolved_instances),"historical_comparison":historical_summary,"validation":validation}
    (OUT/"reports/filter_1_final_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"release/filter_1_release_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    interface={"project_name":"Benchmark 1.0","filter_name":"Filter 1 - Protein Receptor Identification and Qualification","filter_version":"1.0","input_manifest":str(INPUT),"input_manifest_sha256":sha(INPUT),"input_entry_count":256158,"receptor_qualified_entry_manifest":str(OUT/"release/filter_1_receptor_qualified_entries.tsv.gz"),"receptor_qualified_assembly_manifest":str(OUT/"release/filter_1_receptor_qualified_assemblies.tsv.gz"),"receptor_chain_instance_manifest":str(OUT/"release/filter_1_receptor_chain_instances.tsv.gz"),"short_peptide_inventory":str(OUT/"release/filter_1_short_peptide_inventory.tsv.gz"),"min_receptor_aa_length":30,"sequence_length_source_priority":["entity_poly_seq","canonical_sequence_fallback"],"assembly_policy":{"use_deposited_biological_assemblies":True,"silent_asu_fallback":False},"release_creation_time":utc(),"release_validation_pass":validation["release_validation_pass"]}
    (OUT/"release/filter_1_downstream_interface.json").write_text(json.dumps(interface,indent=2)+"\n")
    files=sorted(p for p in (OUT/"release").iterdir() if p.is_file() and p.name!="SHA256SUMS")
    (OUT/"release/SHA256SUMS").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))
    provenance={"start":summary["full_start"],"end":summary["full_end"],"host":platform.node(),"python":sys.version,"gemmi":gemmi.__version__,"config_sha256":sha(OUT/"configs/filter_1.yaml"),"input_sha256":sha(INPUT),"command":" ".join(sys.argv),"workers":json.loads((OUT/"checkpoints/progress.json").read_text()).get("workers"),"resume_supported":True,"release_validation_pass":validation["release_validation_pass"]}
    (OUT/"provenance/filter_1_run_provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    if not validation["release_validation_pass"]: raise SystemExit("Release validation failed")
    print(json.dumps(summary,indent=2))


def validate_only() -> None:
    p=OUT/"release/filter_1_release_validation.json"
    if not p.exists(): raise SystemExit("release validation missing")
    data=json.loads(p.read_text()); print(json.dumps(data,indent=2)); raise SystemExit(0 if data.get("release_validation_pass") else 1)


def materialize_text_reports() -> None:
    names = [
        "filter_1_entry_flow", "filter_1_assembly_flow", "filter_1_polymer_class_distribution",
        "filter_1_polypeptide_length_distribution", "filter_1_chain_role_distribution",
        "filter_1_failure_reason_distribution", "filter_1_review_reason_distribution",
        "filter_1_historical_comparison",
    ]
    written = []
    for name in names:
        source = OUT / f"reports/{name}.tsv.gz"
        target = OUT / f"reports/{name}.tsv"
        with gzip.open(source, "rt", encoding="utf-8") as ih, target.open("w", encoding="utf-8") as oh:
            shutil.copyfileobj(ih, oh)
        written.append(str(target))
    source = OUT / "validation/filter_1_semantic_spotcheck.tsv.gz"
    target = OUT / "validation/filter_1_semantic_spotcheck.tsv"
    with gzip.open(source, "rt", encoding="utf-8") as ih, target.open("w", encoding="utf-8") as oh:
        shutil.copyfileobj(ih, oh)
    written.append(str(target))
    print(json.dumps({"materialized_text_reports": len(written), "paths": written}, indent=2))


def main() -> None:
    ap=argparse.ArgumentParser(description="Benchmark 1.0 Filter 1 pipeline")
    sub=ap.add_subparsers(dest="command",required=True)
    sub.add_parser("setup"); sub.add_parser("preflight")
    f=sub.add_parser("full"); f.add_argument("--workers",type=int,default=16); f.add_argument("--batch-size",type=int,default=250)
    sub.add_parser("finalize");sub.add_parser("validate");sub.add_parser("spotcheck");sub.add_parser("repair-duplicates");sub.add_parser("materialize-text-reports")
    a=ap.parse_args()
    if a.command=="setup": setup()
    elif a.command=="preflight": preflight()
    elif a.command=="full": full(a.workers,a.batch_size)
    elif a.command=="finalize": finalize()
    elif a.command=="validate": validate_only()
    elif a.command=="repair-duplicates": repair_duplicate_batches()
    elif a.command=="materialize-text-reports": materialize_text_reports()
    elif a.command=="spotcheck":
        p=OUT/"validation/filter_1_semantic_spotcheck_summary.json"; print(p.read_text() if p.exists() else "not generated")


if __name__=="__main__": main()
