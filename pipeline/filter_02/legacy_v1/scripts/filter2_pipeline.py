from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import platform
import re
import shutil
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import yaml
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Lipinski


ROOT = Path("/root/autodl-tmp/benchmark_1.0")
OUT = ROOT / "filter_2_ligand_qualification"
F1 = ROOT / "filter_1_protein_receptor_qualification"
P1 = ROOT / "processing_1_pdb_source_audit"
CCD_SOURCE = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v2_provisional/ccd_snapshot/components.cif.gz")
CCD_META = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v2_provisional/ccd_snapshot_metadata.json")
LOCAL_POLICY = Path("/root/autodl-tmp/vs_benchmark/configs/refinement_v2_provisional_policy.yaml")
ARTIFACT_INVENTORY = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v5_biological_relevance/references/artifact_reference_inventory.tsv")
ARTIFACT_VALIDATION = Path("/root/autodl-tmp/vs_benchmark/data_refinement_v5_biological_relevance/references/artifact_reference_validation.json")
RULE_VERSION = "filter_2_v1.0"

ENTRY_INPUT = F1 / "release/filter_1_receptor_qualified_entries.tsv.gz"
ASSEMBLY_INPUT = F1 / "release/filter_1_receptor_qualified_assemblies.tsv.gz"
RECEPTOR_INPUT = F1 / "release/filter_1_receptor_chain_instances.tsv.gz"
SHORT_INPUT = F1 / "release/filter_1_short_peptide_inventory.tsv.gz"
MM_INDEX = P1 / "release/processing_1_mmcif_index.tsv.gz"

METALS = {"LI","NA","K","RB","CS","BE","MG","CA","SR","BA","AL","GA","IN","SN","PB","V","CR","MN","FE","CO","NI","CU","ZN","CD","HG","MO","W","RU","RH","PD","AG","PT","AU","Y","LA","CE","PR","ND","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU"}
WATER = {"HOH","DOD","WAT"}

COMP_FIELDS = ["original_component_id","resolved_ccd_id","ccd_identity_status","ccd_release_version","ccd_name","ccd_type","formula","formula_weight","formal_charge","ccd_parent_component_id","topology_source","topology_fallback_used","topology_status","rdkit_parse_status","rdkit_sanitize_status","valence_status","aromaticity_status","kekulization_status","fragment_count","heavy_atom_count","molecular_weight","carbon_count","hetero_atom_count","ring_count","rotatable_bond_count","element_set","organic_status","contains_metal","chemical_entity_class","artifact_prior","artifact_sources","artifact_reason","artifact_list_version","cofactor_prior","covalent_warhead_status","classification_reason","rule_version"]
SOURCE_FIELDS = ["pdb_id","selected_model_id","entity_id","source_label_asym_id","source_auth_asym_id","label_comp_id","auth_comp_id","label_seq_id","auth_seq_id","insertion_code","source_component_instance_id","polymer_context","modified_residue_status","entry_parent_component_id","resolved_ccd_id","atom_count","heavy_atom_count","altloc_values","conformer_count","occupancy_min","occupancy_max","instance_covalent_link_status","instance_metal_link_status","chemical_entity_class","artifact_prior","filter_2_route","primary_action","classification_reason","instance_status"]
ASSEMBLY_FIELDS = ["pdb_id","assembly_id","selected_model_id","source_component_instance_id","source_label_asym_id","source_auth_asym_id","operator_id","composite_operator_id","assembly_component_instance_id","resolved_ccd_id","polymer_context","chemical_entity_class","artifact_prior","instance_covalent_link_status","instance_metal_link_status","filter_2_route","assembly_membership_status","instance_status"]
CONFORMER_FIELDS = ["source_component_instance_id","component_conformer_id","altloc_id","shared_blank_altloc_atom_count","conformer_atom_count","occupancy_min","occupancy_max","conformer_status"]
PARENT_FIELDS = ["pdb_id","source_component_instance_id","original_component_id","ccd_parent_component_id","entry_parent_component_id","parent_mapping_status","parent_mapping_reason"]
CONN_FIELDS = ["pdb_id","conn_id","conn_type_id","partner_1_label_asym_id","partner_1_label_comp_id","partner_1_label_seq_id","partner_1_auth_asym_id","partner_1_auth_seq_id","partner_2_label_asym_id","partner_2_label_comp_id","partner_2_label_seq_id","partner_2_auth_asym_id","partner_2_auth_seq_id","component_instance_id","receptor_partner_status","qualified_assembly_mapping_status","covalent_link_status","metal_link_status","mapping_reason"]
ENTRY_FIELDS = ["pdb_id","parse_status","parse_error","raw_component_instance_count","nonpolymer_instance_count","short_peptide_instance_count","nucleic_acid_instance_count","branched_instance_count","modified_residue_instance_count","water_count","qualified_assembly_count","has_any_candidate_component","has_ordinary_candidate","has_special_candidate","has_artifact_review","has_unresolved_component","entry_status","terminal_reason"]
CATEGORY_FIELDS = ["pdb_id","category","category_present","row_count","parse_status","parse_warning"]
TABLE_FIELDS = {"entries":ENTRY_FIELDS,"sources":SOURCE_FIELDS,"assemblies":ASSEMBLY_FIELDS,"conformers":CONFORMER_FIELDS,"parents":PARENT_FIELDS,"connections":CONN_FIELDS,"categories":CATEGORY_FIELDS}

G_CCD = {}
G_ASSEMBLIES = {}
G_RECEPTOR_ASYM = {}
G_SHORT_ASYM = {}


def utc(): return datetime.now(timezone.utc).isoformat()
def sha(path: Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(8<<20),b""): h.update(c)
    return h.hexdigest()
def clean(x):
    s=gemmi.cif.as_string(str(x)); return "" if s in {".","?"} else s.strip()
def rows(block,tags):
    try:return [[clean(x) for x in r] for r in block.find(tags)]
    except:return []
def write_tsv(path,data,fields,gz=False):
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); op=gzip.open if gz else open
    with op(tmp,"wt",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore");w.writeheader();w.writerows(data)
    os.replace(tmp,path)
def iter_tsv(path):
    op=gzip.open if str(path).endswith(".gz") else open
    with op(path,"rt",encoding="utf-8",newline="") as h: yield from csv.DictReader(h,delimiter="\t")
def expand_ops(expression):
    expr=re.sub(r"\s+","",expression); groups=re.findall(r"\(([^()]*)\)",expr) or [expr]
    parsed=[]
    for group in groups:
        vals=[]
        for token in group.split(","):
            m=re.fullmatch(r"(\d+)-(\d+)",token)
            if m:
                a,b=map(int,m.groups()); vals.extend(str(x) for x in range(a,b+(1 if b>=a else -1),1 if b>=a else -1))
            elif token: vals.append(token)
        if not vals: raise ValueError("invalid operator expression")
        parsed.append(vals)
    return list(dict.fromkeys("x".join(x) for x in itertools.product(*parsed)))


def setup():
    if OUT.exists(): raise SystemExit(f"Output exists: {OUT}")
    for d in ["configs","scripts","schemas","references","inputs","preflight","checkpoints/batches","full","reports","release","validation","logs","provenance"]:(OUT/d).mkdir(parents=True,exist_ok=True)
    config={"project_root":str(ROOT),"output_root":str(OUT),"filter_1_root":str(F1),"processing_1_root":str(P1),"ccd_snapshot":str(OUT/"references/components.cif.gz"),"ccd_release_version":"Sat, 11 Jul 2026 12:01:19 GMT","rule_version":RULE_VERSION,"artifact_policy":{"official_reference_status":"unavailable_local","local_project_prior":str(OUT/"references/refinement_v2_provisional_policy.yaml"),"unknown_artifact_action":"unresolved_review"},"assembly_policy":{"qualified_deposited_assemblies_only":True,"materialize_coordinates":False},"full_run":{"workers":16,"batch_size":200,"resume":True}}
    (OUT/"configs/filter_2.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    for src,name in [(ENTRY_INPUT,"filter_1_receptor_qualified_entries_snapshot.tsv.gz"),(ASSEMBLY_INPUT,"filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"),(RECEPTOR_INPUT,"filter_1_receptor_chain_instances_snapshot.tsv.gz"),(SHORT_INPUT,"filter_1_short_peptide_inventory_snapshot.tsv.gz"),(MM_INDEX,"processing_1_mmcif_index_snapshot.tsv.gz")]:shutil.copy2(src,OUT/"inputs"/name)
    for src,name in [(CCD_SOURCE,"components.cif.gz"),(CCD_META,"ccd_snapshot_metadata.json"),(LOCAL_POLICY,"refinement_v2_provisional_policy.yaml"),(ARTIFACT_INVENTORY,"artifact_reference_inventory.tsv"),(ARTIFACT_VALIDATION,"artifact_reference_validation.json")]:shutil.copy2(src,OUT/"references"/name)
    checks={str(p):sha(p) for p in (OUT/"inputs").iterdir() if p.is_file()}; checks.update({str(p):sha(p) for p in (OUT/"references").iterdir() if p.is_file()})
    (OUT/"inputs/input_checksums.json").write_text(json.dumps(checks,indent=2)+"\n")
    (OUT/"README.md").write_text("# Filter 2 - Ligand Instance Identification and Chemical-Scope Qualification\n\nThis filter enumerates non-receptor structural objects in Filter 1-qualified entries, resolves frozen CCD identity, and assigns chemical-scope routes. It does not construct protein-component pairs, calculate distances, materialize assembly coordinates, or run interaction/quality tools. Ordinary output means chemical-scope candidate only.\n")
    shutil.copy2(Path(__file__),OUT/"scripts/filter2_pipeline.py")
    for name,cmd in [("audit_filter_2_inputs.py","audit"),("prepare_filter_2_references.py","prepare-references"),("run_filter_2_preflight.py","preflight"),("run_filter_2_full.py","full"),("build_filter_2_release.py","finalize"),("validate_filter_2_release.py","validate")]:
        (OUT/"scripts"/name).write_text(f"#!/usr/bin/env python3\nimport subprocess,sys\nraise SystemExit(subprocess.call([sys.executable,{repr(str(OUT/'scripts/filter2_pipeline.py'))},'{cmd}',*sys.argv[1:]]))\n")
    for key,fields in TABLE_FIELDS.items():write_tsv(OUT/f"schemas/{key}_schema.tsv",[{"field":f,"required":"true"} for f in fields],["field","required"])
    write_tsv(OUT/"schemas/component_classification_schema.tsv",[{"field":f,"required":"true"} for f in COMP_FIELDS],["field","required"])
    print(json.dumps({"setup":True,"root":str(OUT)},indent=2))


def audit():
    counts={"entries":sum(1 for _ in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")),"assemblies":sum(1 for _ in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"))}
    entry_ids=[r["pdb_id"] for r in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")]
    assembly_keys=[r["pdb_id"]+"|"+r["assembly_id"] for r in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_assemblies_snapshot.tsv.gz")]
    meta=json.loads((OUT/"references/ccd_snapshot_metadata.json").read_text())
    data={"input_entry_count":counts["entries"],"unique_entry_count":len(set(entry_ids)),"duplicate_entry_count":len(entry_ids)-len(set(entry_ids)),"qualified_assembly_count":counts["assemblies"],"unique_qualified_assembly_keys":len(set(assembly_keys)),"duplicate_qualified_assembly_key":len(assembly_keys)-len(set(assembly_keys)),"filter_1_release_validation":json.loads((F1/"release/filter_1_release_validation.json").read_text())["release_validation_pass"],"ccd_snapshot_sha256":sha(OUT/"references/components.cif.gz"),"ccd_metadata_sha256":meta["sha256"],"ccd_checksum_match":sha(OUT/"references/components.cif.gz")==meta["sha256"],"artifact_official_reference_available":False,"artifact_reference_limitation":"No official traceable BioLiP/BioLiP2/Q-BioLiP artifact list is locally frozen; local project prior is used and explicitly labeled.","static_audit_pass":False}
    data["static_audit_pass"]=data["input_entry_count"]==248037 and data["unique_entry_count"]==248037 and data["duplicate_entry_count"]==0 and data["qualified_assembly_count"]==360611 and data["duplicate_qualified_assembly_key"]==0 and data["filter_1_release_validation"] and data["ccd_checksum_match"]
    (OUT/"preflight/input_audit.json").write_text(json.dumps(data,indent=2)+"\n");print(json.dumps(data,indent=2))
    if not data["static_audit_pass"]:raise SystemExit(1)


def ccd_mol(block):
    atom_rows=rows(block,["_chem_comp_atom.atom_id","_chem_comp_atom.type_symbol","_chem_comp_atom.charge"])
    bond_rows=rows(block,["_chem_comp_bond.atom_id_1","_chem_comp_bond.atom_id_2","_chem_comp_bond.value_order","_chem_comp_bond.pdbx_aromatic_flag"])
    if not atom_rows:return None,"missing_atom_table"
    rw=Chem.RWMol(); idx={}
    try:
        for atom_id,element,charge in atom_rows:
            atom=Chem.Atom(element.title()); atom.SetFormalCharge(int(float(charge or 0))); idx[atom_id]=rw.AddAtom(atom)
        bmap={"SING":Chem.BondType.SINGLE,"DOUB":Chem.BondType.DOUBLE,"TRIP":Chem.BondType.TRIPLE,"AROM":Chem.BondType.AROMATIC,"DELO":Chem.BondType.AROMATIC}
        for a,b,order,arom in bond_rows:
            if a in idx and b in idx and rw.GetBondBetweenAtoms(idx[a],idx[b]) is None:rw.AddBond(idx[a],idx[b],bmap.get(order.upper(),Chem.BondType.SINGLE))
        mol=rw.GetMol(); Chem.SanitizeMol(mol); return mol,"ccd_atom_bond"
    except Exception as exc:return None,"rdkit_graph_error:"+str(exc)[:200]


def prepare_references():
    policy=yaml.safe_load((OUT/"references/refinement_v2_provisional_policy.yaml").read_text()); rules=policy["component_id_rules"]
    artifact=set(rules.get("organic_solvent",[]))|set(rules.get("buffer_or_crystallization_additive",[])); cof=set(rules.get("cofactor_or_coenzyme",[])); nuc=set(rules.get("nucleotide",[])); glyc=set(rules.get("glycan",[])); inorg=set(rules.get("inorganic_ion",[]))|set(rules.get("single_atom_metal",[])); pep=set(rules.get("peptide_or_polymer_component",[]))
    doc=gemmi.cif.read(str(OUT/"references/components.cif.gz")); output=[]; audit_counts=Counter(); seen=set()
    for block in doc:
        cid=block.name.upper(); seen.add(cid)
        def value(tag):
            try:return clean(block.find_value(tag))
            except:return ""
        name=value("_chem_comp.name"); ctype=value("_chem_comp.type"); formula=value("_chem_comp.formula"); fw=value("_chem_comp.formula_weight"); charge=value("_chem_comp.pdbx_formal_charge"); parent=value("_chem_comp.mon_nstd_parent_comp_id"); release=value("_chem_comp.pdbx_release_status"); replaced=value("_chem_comp.pdbx_replaced_by")
        mol,topology=ccd_mol(block); atom_rows=rows(block,["_chem_comp_atom.type_symbol"]); elements=[x[0].upper() for x in atom_rows]; heavy=sum(x!="H" for x in elements); carbon=elements.count("C"); contains=any(x in METALS for x in elements)
        if mol is not None:
            rdparse=rdsan="pass"; fragments=len(Chem.GetMolFrags(mol)); mw=Descriptors.MolWt(mol); rings=Lipinski.RingCount(mol); rot=Lipinski.NumRotatableBonds(mol)
        else:
            rdparse=rdsan="failed"; fragments=0; mw=float(fw) if fw else 0; rings=rot=0
        if cid in WATER: chemclass="water"
        elif cid in inorg or (carbon==0 and (heavy<=3 or contains)): chemclass="metal_or_inorganic"
        elif cid in cof: chemclass="cofactor_or_coenzyme"
        elif cid in nuc or "NUCLEOTIDE" in ctype.upper(): chemclass="nucleotide_or_nucleoside"
        elif cid in glyc or "SACCHARIDE" in ctype.upper(): chemclass="carbohydrate"
        elif cid in pep or "PEPTIDE" in ctype.upper(): chemclass="peptide"
        elif contains and carbon>0: chemclass="organometallic"
        elif carbon>0: chemclass="small_organic"
        else: chemclass="unknown"
        art="high" if cid in artifact else ("unknown" if carbon>0 else "none"); cofactor="high" if cid in cof or cid in nuc else "none"
        identity="obsolete_id_resolved" if release.upper() in {"OBS","OBSOLETE"} and replaced else ("ccd_invalid" if release.upper() in {"OBS","OBSOLETE"} else "exact_ccd_match")
        output.append({"original_component_id":cid,"resolved_ccd_id":replaced or cid,"ccd_identity_status":identity,"ccd_release_version":"Sat, 11 Jul 2026 12:01:19 GMT","ccd_name":name,"ccd_type":ctype,"formula":formula,"formula_weight":fw,"formal_charge":charge,"ccd_parent_component_id":parent,"topology_source":"ccd_atom_bond" if topology=="ccd_atom_bond" else "unavailable","topology_fallback_used":"false","topology_status":topology,"rdkit_parse_status":rdparse,"rdkit_sanitize_status":rdsan,"valence_status":rdsan,"aromaticity_status":rdsan,"kekulization_status":rdsan,"fragment_count":fragments,"heavy_atom_count":heavy,"molecular_weight":f"{mw:.4f}","carbon_count":carbon,"hetero_atom_count":sum(x not in {"C","H"} for x in elements),"ring_count":rings,"rotatable_bond_count":rot,"element_set":",".join(sorted(set(elements))),"organic_status":"organic" if carbon>0 else "inorganic","contains_metal":str(contains).lower(),"chemical_entity_class":chemclass,"artifact_prior":art,"artifact_sources":"local_project_refinement_v2_provisional_policy" if cid in artifact else "none_available","artifact_reason":"listed_local_project_prior" if cid in artifact else ("official_reference_unavailable" if art=="unknown" else "not_applicable"),"artifact_list_version":policy["policy_version"],"cofactor_prior":cofactor,"covalent_warhead_status":"not_assessed_lightweight","classification_reason":"ccd_graph_and_versioned_local_policy","rule_version":RULE_VERSION})
        audit_counts["missing_atom_table"]+=not bool(atom_rows); audit_counts["missing_bond_table"]+=len(rows(block,["_chem_comp_bond.atom_id_1"]))==0; audit_counts["missing_descriptor"]+=len(rows(block,["_pdbx_chem_comp_descriptor.comp_id"]))==0; audit_counts["obsolete"]+=release.upper() in {"OBS","OBSOLETE"}; audit_counts["parent"]+=bool(parent)
    write_tsv(OUT/"references/ccd_component_cache.tsv.gz",output,COMP_FIELDS,True)
    audit_data={"component_count":len(output),"unique_component_id":len(seen),"duplicate_component_id":len(output)-len(seen),**audit_counts,"snapshot_sha256":sha(OUT/"references/components.cif.gz"),"snapshot_version":"Sat, 11 Jul 2026 12:01:19 GMT","rdkit_version":rdBase.rdkitVersion,"reference_validation_pass":len(output)==50666 and len(seen)==50666}
    (OUT/"preflight/ccd_reference_audit.json").write_text(json.dumps(audit_data,indent=2)+"\n");print(json.dumps(audit_data,indent=2))
    if not audit_data["reference_validation_pass"]:raise SystemExit(1)


def load_globals():
    global G_CCD,G_ASSEMBLIES,G_RECEPTOR_ASYM,G_SHORT_ASYM
    G_CCD={r["original_component_id"]:r for r in iter_tsv(OUT/"references/ccd_component_cache.tsv.gz")}
    a=defaultdict(set)
    for r in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_assemblies_snapshot.tsv.gz"):a[r["pdb_id"]].add(r["assembly_id"])
    G_ASSEMBLIES=dict(a)
    rec=defaultdict(set)
    for r in iter_tsv(OUT/"inputs/filter_1_receptor_chain_instances_snapshot.tsv.gz"):rec[r["pdb_id"]].add(r["label_asym_id"])
    G_RECEPTOR_ASYM=dict(rec)
    short=defaultdict(set)
    for r in iter_tsv(OUT/"inputs/filter_1_short_peptide_inventory_snapshot.tsv.gz"):short[r["pdb_id"]].add(r["label_asym_id"])
    G_SHORT_ASYM=dict(short)


def route_instance(context,ccd,covalent):
    if context=="water":return "water_excluded","exclude","water"
    if context in {"protein_polymer_residue","modified_protein_residue","other_polymer"}:return "polymer_or_modified_residue","retain_audit","polymer_context"
    if context=="short_peptide":return "peptide_special","special","short_peptide"
    if context in {"rna_polymer","dna_polymer","hybrid_nucleic_acid"}:return "nucleic_acid_special","special","nucleic_acid_context"
    if context=="branched_glycan":return "glycan_or_carbohydrate_special","special","branched_context"
    if context in {"context_conflict","context_unresolved"}:return "unresolved_review","review","polymer_context_unresolved"
    if not ccd:return "unresolved_review","review","ccd_missing"
    if covalent=="declared_receptor_covalent":return "covalent_or_linked_special","special","explicit_receptor_covalent_link"
    cls=ccd["chemical_entity_class"]
    if cls in {"cofactor_or_coenzyme","nucleotide_or_nucleoside"}:return "cofactor_or_nucleotide_special","special","ccd_special_class"
    if cls=="carbohydrate":return "glycan_or_carbohydrate_special","special","ccd_carbohydrate"
    if cls=="peptide":return "peptide_special","special","ccd_peptide"
    if cls=="metal_or_inorganic" or cls=="organometallic":return "metal_or_inorganic_special","special","metal_or_inorganic"
    if ccd["artifact_prior"]=="high":return "artifact_review","review","high_local_artifact_prior"
    if ccd["ccd_identity_status"] not in {"exact_ccd_match","obsolete_id_resolved"} or ccd["rdkit_parse_status"]!="pass":return "unresolved_review","review","ccd_or_topology_unresolved"
    if cls=="small_organic":return "ordinary_small_molecule_candidate","retain_candidate","independent_resolved_small_organic"
    return "unresolved_review","review","unclassified_chemical_scope"


def parse_entry(item):
    pid,path=item; out={k:[] for k in TABLE_FIELDS}
    try:
        b=gemmi.cif.read(path).sole_block(); qassemblies=G_ASSEMBLIES.get(pid,set()); receptor=G_RECEPTOR_ASYM.get(pid,set()); short=G_SHORT_ASYM.get(pid,set())
        cat_tags={"entity":"_entity.id","entity_poly":"_entity_poly.entity_id","entity_poly_seq":"_entity_poly_seq.entity_id","struct_asym":"_struct_asym.id","atom_site":"_atom_site.label_asym_id","nonpoly_scheme":"_pdbx_nonpoly_scheme.asym_id","entity_nonpoly":"_pdbx_entity_nonpoly.entity_id","entity_branch":"_pdbx_entity_branch.entity_id","branch_scheme":"_pdbx_branch_scheme.asym_id","struct_mod_residue":"_pdbx_struct_mod_residue.id","assembly":"_pdbx_struct_assembly.id","assembly_gen":"_pdbx_struct_assembly_gen.assembly_id","oper_list":"_pdbx_struct_oper_list.id","struct_conn":"_struct_conn.id"}
        for name,tag in cat_tags.items():
            try:n=len(b.find([tag]));out["categories"].append({"pdb_id":pid,"category":name,"category_present":str(n>0).lower(),"row_count":n,"parse_status":"parsed","parse_warning":""})
            except Exception as e:out["categories"].append({"pdb_id":pid,"category":name,"category_present":"false","row_count":0,"parse_status":"failed","parse_warning":str(e)[:300]})
        poly_type={r[0]:r[1] for r in rows(b,["_entity_poly.entity_id","_entity_poly.type"])}; asym_entity={r[0]:r[1] for r in rows(b,["_struct_asym.id","_struct_asym.entity_id"])}; branch_entities={r[0] for r in rows(b,["_pdbx_entity_branch.entity_id"])}
        branch_asym={r[0] for r in rows(b,["_pdbx_branch_scheme.asym_id"])}
        mod_rows=rows(b,["_pdbx_struct_mod_residue.label_asym_id","_pdbx_struct_mod_residue.label_seq_id","_pdbx_struct_mod_residue.auth_asym_id","_pdbx_struct_mod_residue.auth_seq_id","_pdbx_struct_mod_residue.label_comp_id","_pdbx_struct_mod_residue.parent_comp_id"])
        mod_label={(r[0],r[1]):r[5] for r in mod_rows if r[0]}; mod_auth={(r[2],r[3]):r[5] for r in mod_rows if r[2]}
        auth_map=defaultdict(set); atom_groups={}; water_instances=set()
        atom_tags=["_atom_site.pdbx_PDB_model_num","_atom_site.label_entity_id","_atom_site.label_asym_id","_atom_site.auth_asym_id","_atom_site.label_comp_id","_atom_site.auth_comp_id","_atom_site.label_seq_id","_atom_site.auth_seq_id","_atom_site.pdbx_PDB_ins_code","_atom_site.label_atom_id","_atom_site.type_symbol","_atom_site.label_alt_id","_atom_site.occupancy"]
        for raw in b.find(atom_tags):
            model,entity,lasym,aasym,lcomp,acomp,lseq,aseq,ins,atom,element,alt,occ=[clean(x) for x in raw]; model=model or "1"; auth_map[lasym].add(aasym)
            eid=entity or asym_entity.get(lasym,""); ptype=poly_type.get(eid,"").lower()
            if lcomp.upper() in WATER: context="water"
            elif lasym in branch_asym or eid in branch_entities: context="branched_glycan"
            elif lasym in short: context="short_peptide"
            elif "polyribonucleotide" in ptype or ptype=="rna": context="rna_polymer"
            elif "polydeoxyribonucleotide" in ptype or ptype=="dna": context="dna_polymer"
            elif "hybrid" in ptype: context="hybrid_nucleic_acid"
            elif "polypeptide" in ptype:
                context="modified_protein_residue" if (lasym,lseq) in mod_label or (aasym,aseq) in mod_auth else "protein_polymer_residue"
                if context=="protein_polymer_residue" and lasym in receptor:continue
            elif ptype: context="other_polymer"
            else: context="independent_nonpolymer"
            if context=="water":
                water_instances.add((model,lasym,lseq or aseq or "?",ins,lcomp.upper()))
                continue
            rid=lseq or aseq or "?"; key=(model,lasym,rid,ins,lcomp)
            if key not in atom_groups:atom_groups[key]={"model":model,"entity":eid,"lasym":lasym,"aasym":aasym,"lcomp":lcomp.upper(),"acomp":(acomp or lcomp).upper(),"lseq":lseq,"aseq":aseq,"ins":ins,"context":context,"atoms":[],"entry_parent":mod_label.get((lasym,lseq),mod_auth.get((aasym,aseq),""))}
            atom_groups[key]["atoms"].append((atom,element.upper(),alt,occ))
        # Assembly asym/operator map for qualified assemblies.
        genmap=defaultdict(list)
        for aid,expr,asyms in rows(b,["_pdbx_struct_assembly_gen.assembly_id","_pdbx_struct_assembly_gen.oper_expression","_pdbx_struct_assembly_gen.asym_id_list"]):
            if aid not in qassemblies:continue
            try:ops=expand_ops(expr)
            except:continue
            for asym in [x.strip() for x in asyms.split(",") if x.strip()]:genmap[asym].append((aid,expr,ops))
        # First create instances so struct_conn can map to them.
        source_by_label={}; source_rows=[]
        for g in atom_groups.values():
            sid=f"{pid}|{g['model']}|{g['lasym']}|{g['lseq'] or g['aseq']}|{g['ins'] or '.'}|{g['lcomp']}"; ccd=G_CCD.get(g["lcomp"]); resolved=ccd["resolved_ccd_id"] if ccd else ""; atoms=g["atoms"]; alts=sorted({x[2] for x in atoms if x[2]}); occup=[float(x[3]) for x in atoms if x[3]]; alt_nonblank=[x for x in alts if x not in {".","?",""}]; conformer_ids=alt_nonblank or ["."]
            row={"pdb_id":pid,"selected_model_id":g["model"],"entity_id":g["entity"],"source_label_asym_id":g["lasym"],"source_auth_asym_id":g["aasym"],"label_comp_id":g["lcomp"],"auth_comp_id":g["acomp"],"label_seq_id":g["lseq"],"auth_seq_id":g["aseq"],"insertion_code":g["ins"],"source_component_instance_id":sid,"polymer_context":g["context"],"modified_residue_status":str(g["context"]=="modified_protein_residue").lower(),"entry_parent_component_id":g["entry_parent"],"resolved_ccd_id":resolved,"atom_count":len(atoms),"heavy_atom_count":sum(x[1]!="H" for x in atoms),"altloc_values":",".join(alts),"conformer_count":len(conformer_ids),"occupancy_min":min(occup) if occup else "","occupancy_max":max(occup) if occup else "","instance_covalent_link_status":"not_declared","instance_metal_link_status":"not_declared","chemical_entity_class":ccd["chemical_entity_class"] if ccd else "unknown","artifact_prior":ccd["artifact_prior"] if ccd else "unknown","filter_2_route":"","primary_action":"","classification_reason":"","instance_status":"resolved" if ccd else "ccd_missing"}
            source_rows.append(row);source_by_label[(g["lasym"],g["lseq"],g["lcomp"])]=row;source_by_label[(g["lasym"],g["aseq"],g["lcomp"])]=row
            blank=sum(x[2] in {"",".","?"} for x in atoms)
            for alt in conformer_ids:
                selected=[x for x in atoms if x[2] in {"",".","?",alt}];oc=[float(x[3]) for x in selected if x[3]];out["conformers"].append({"source_component_instance_id":sid,"component_conformer_id":sid+"|alt="+alt,"altloc_id":alt,"shared_blank_altloc_atom_count":blank,"conformer_atom_count":len(selected),"occupancy_min":min(oc) if oc else "","occupancy_max":max(oc) if oc else "","conformer_status":"resolved"})
            cp=ccd["ccd_parent_component_id"] if ccd else ""; ep=g["entry_parent"]; status="agree" if cp and ep and cp==ep else "ccd_only" if cp and not ep else "entry_only" if ep and not cp else "conflict" if cp and ep else "none"
            out["parents"].append({"pdb_id":pid,"source_component_instance_id":sid,"original_component_id":g["lcomp"],"ccd_parent_component_id":cp,"entry_parent_component_id":ep,"parent_mapping_status":status,"parent_mapping_reason":"ccd_and_entry_parent_mapping"})
        # struct_conn mapping and orthogonal statuses.
        conn_tags=["_struct_conn.id","_struct_conn.conn_type_id","_struct_conn.ptnr1_label_asym_id","_struct_conn.ptnr1_label_comp_id","_struct_conn.ptnr1_label_seq_id","_struct_conn.ptnr1_auth_asym_id","_struct_conn.ptnr1_auth_seq_id","_struct_conn.ptnr2_label_asym_id","_struct_conn.ptnr2_label_comp_id","_struct_conn.ptnr2_label_seq_id","_struct_conn.ptnr2_auth_asym_id","_struct_conn.ptnr2_auth_seq_id"]
        for cr in rows(b,conn_tags):
            cid,ctype,a1,c1,s1,aa1,as1,a2,c2,s2,aa2,as2=cr; candidate=source_by_label.get((a1,s1,c1.upper())) or source_by_label.get((a1,as1,c1.upper())) or source_by_label.get((a2,s2,c2.upper())) or source_by_label.get((a2,as2,c2.upper())); partner_receptor=(a2 in receptor if candidate and candidate["source_label_asym_id"]==a1 else a1 in receptor) if candidate else False
            cov="declared_receptor_covalent" if candidate and partner_receptor and ctype.lower().startswith("covale") else "not_covalent"; metal="declared_metal_link" if candidate and ctype.lower().startswith("metalc") else "not_metal"
            if candidate:
                if cov.startswith("declared"):candidate["instance_covalent_link_status"]=cov
                if metal.startswith("declared"):candidate["instance_metal_link_status"]=metal
            mapped="mapped_to_qualified_assembly" if candidate and genmap.get(candidate["source_label_asym_id"]) else "not_mapped"
            out["connections"].append({"pdb_id":pid,"conn_id":cid,"conn_type_id":ctype,"partner_1_label_asym_id":a1,"partner_1_label_comp_id":c1,"partner_1_label_seq_id":s1,"partner_1_auth_asym_id":aa1,"partner_1_auth_seq_id":as1,"partner_2_label_asym_id":a2,"partner_2_label_comp_id":c2,"partner_2_label_seq_id":s2,"partner_2_auth_asym_id":aa2,"partner_2_auth_seq_id":as2,"component_instance_id":candidate["source_component_instance_id"] if candidate else "","receptor_partner_status":str(partner_receptor).lower(),"qualified_assembly_mapping_status":mapped,"covalent_link_status":cov,"metal_link_status":metal,"mapping_reason":"label_identifier_mapping" if candidate else "component_instance_not_found"})
        # Route and logical assembly instances after link status is known.
        for r in source_rows:
            ccd=G_CCD.get(r["label_comp_id"]);route,action,reason=route_instance(r["polymer_context"],ccd,r["instance_covalent_link_status"]);r["filter_2_route"]=route;r["primary_action"]=action;r["classification_reason"]=reason;out["sources"].append(r)
            for aid,expr,ops in genmap.get(r["source_label_asym_id"],[]):
                for op in ops:
                    ai=f"{pid}|{aid}|{r['selected_model_id']}|{r['source_component_instance_id']}|{op}";out["assemblies"].append({"pdb_id":pid,"assembly_id":aid,"selected_model_id":r["selected_model_id"],"source_component_instance_id":r["source_component_instance_id"],"source_label_asym_id":r["source_label_asym_id"],"source_auth_asym_id":r["source_auth_asym_id"],"operator_id":op,"composite_operator_id":op,"assembly_component_instance_id":ai,"resolved_ccd_id":r["resolved_ccd_id"],"polymer_context":r["polymer_context"],"chemical_entity_class":r["chemical_entity_class"],"artifact_prior":r["artifact_prior"],"instance_covalent_link_status":r["instance_covalent_link_status"],"instance_metal_link_status":r["instance_metal_link_status"],"filter_2_route":route,"assembly_membership_status":"in_qualified_assembly","instance_status":r["instance_status"]})
        routes=Counter(r["filter_2_route"] for r in source_rows); special=sum(v for k,v in routes.items() if k.endswith("_special"));out["entries"].append({"pdb_id":pid,"parse_status":"success","parse_error":"","raw_component_instance_count":len(source_rows)+len(water_instances),"nonpolymer_instance_count":sum(r["polymer_context"]=="independent_nonpolymer" for r in source_rows),"short_peptide_instance_count":sum(r["polymer_context"]=="short_peptide" for r in source_rows),"nucleic_acid_instance_count":sum(r["polymer_context"] in {"rna_polymer","dna_polymer","hybrid_nucleic_acid"} for r in source_rows),"branched_instance_count":sum(r["polymer_context"]=="branched_glycan" for r in source_rows),"modified_residue_instance_count":sum(r["polymer_context"]=="modified_protein_residue" for r in source_rows),"water_count":len(water_instances),"qualified_assembly_count":len(qassemblies),"has_any_candidate_component":str(bool(source_rows)).lower(),"has_ordinary_candidate":str(routes["ordinary_small_molecule_candidate"]>0).lower(),"has_special_candidate":str(special>0).lower(),"has_artifact_review":str(routes["artifact_review"]>0).lower(),"has_unresolved_component":str(routes["unresolved_review"]>0).lower(),"entry_status":"pass" if source_rows else "no_candidate_component","terminal_reason":"component_inventory_completed" if source_rows else "no_nonreceptor_component"})
    except Exception as e:
        out["entries"]=[{"pdb_id":pid,"parse_status":"failed","parse_error":f"{type(e).__name__}:{e}"[:1000],"raw_component_instance_count":0,"nonpolymer_instance_count":0,"short_peptide_instance_count":0,"nucleic_acid_instance_count":0,"branched_instance_count":0,"modified_residue_instance_count":0,"water_count":0,"qualified_assembly_count":len(G_ASSEMBLIES.get(pid,set())),"has_any_candidate_component":"false","has_ordinary_candidate":"false","has_special_candidate":"false","has_artifact_review":"false","has_unresolved_component":"true","entry_status":"parse_failed","terminal_reason":"parse_failed"}]
    return out


def input_items():
    paths={r["pdb_id"]:r["mmcif_path"] for r in iter_tsv(OUT/"inputs/processing_1_mmcif_index_snapshot.tsv.gz")};return [(r["pdb_id"],paths[r["pdb_id"]]) for r in iter_tsv(OUT/"inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")]


def select_preflight(items):
    # Deterministic broad coverage plus explicit semantic targets discovered from mmCIF metadata.
    chosen={pid:(pid,path) for pid,path in items[::max(1,len(items)//4500)][:4500]}; targets={"ATP","ADP","NAG","MAN","MSE","SEP","TPO","PTR","GOL","PEG","HEM","FAD","NAD","ZN","MG","SO4","PO4","HOH"};found=set()
    for pid,path in items[:50000]:
        if targets<=found:break
        try:
            b=gemmi.cif.read(path).sole_block();comps={clean(x).upper() for x in b.find_values("_atom_site.label_comp_id")}; hit=(targets-found)&comps
            if hit:chosen[pid]=(pid,path);found|=hit
        except:pass
    # Add short peptide and multi-assembly entries from frozen Filter 1 metadata.
    for pid in list(G_SHORT_ASYM)[:250]:
        if pid in dict(items):chosen[pid]=(pid,dict(items)[pid])
    for pid,aids in list(G_ASSEMBLIES.items()):
        if len(aids)>1:chosen[pid]=(pid,dict(items)[pid])
        if len(chosen)>=5000:break
    return list(chosen.values())[:5000],sorted(found)


def run_preflight():
    load_globals();items=input_items();selected,found=select_preflight(items)
    with ProcessPoolExecutor(max_workers=16) as pool:results=list(pool.map(parse_entry,selected,chunksize=4))
    sources=[r for x in results for r in x["sources"]];assemblies=[r for x in results for r in x["assemblies"]];entries=[r for x in results for r in x["entries"]]
    sem={cid:Counter(r["filter_2_route"] for r in sources if r["label_comp_id"]==cid) for cid in found};contamination={"ordinary_polymer":sum(r["filter_2_route"]=="ordinary_small_molecule_candidate" and r["polymer_context"]!="independent_nonpolymer" for r in sources),"ordinary_water":sum(r["filter_2_route"]=="ordinary_small_molecule_candidate" and r["chemical_entity_class"]=="water" for r in sources),"ordinary_metal":sum(r["filter_2_route"]=="ordinary_small_molecule_candidate" and r["chemical_entity_class"]=="metal_or_inorganic" for r in sources),"missing_route":sum(not r["filter_2_route"] for r in sources)}
    data={"preflight_entries":len(entries),"parse_failed":sum(r["parse_status"]=="failed" for r in entries),"source_instances":len(sources),"assembly_instances":len(assemblies),"semantic_components_found":found,"semantic_routes":{k:dict(v) for k,v in sem.items()},"contamination":contamination,"artifact_reference_limitation_documented":True}
    required={"ATP":"cofactor_or_nucleotide_special","ADP":"cofactor_or_nucleotide_special","NAG":"glycan_or_carbohydrate_special","MAN":"glycan_or_carbohydrate_special","MSE":"polymer_or_modified_residue","SEP":"polymer_or_modified_residue","TPO":"polymer_or_modified_residue","PTR":"polymer_or_modified_residue","GOL":"artifact_review","HEM":"cofactor_or_nucleotide_special","FAD":"cofactor_or_nucleotide_special","NAD":"cofactor_or_nucleotide_special","ZN":"metal_or_inorganic_special","MG":"metal_or_inorganic_special","SO4":"metal_or_inorganic_special","PO4":"metal_or_inorganic_special"}
    semantic_ok=all(route in sem.get(cid,{}) for cid,route in required.items());data["semantic_validation_pass"]=semantic_ok;data["preflight_validation_pass"]=len(entries)==5000 and data["parse_failed"]==0 and all(v==0 for v in contamination.values()) and semantic_ok
    write_tsv(OUT/"preflight/filter_2_preflight_entries.tsv",entries,ENTRY_FIELDS);write_tsv(OUT/"preflight/filter_2_preflight_sources.tsv.gz",sources,SOURCE_FIELDS,True);(OUT/"preflight/filter_2_preflight_summary.json").write_text(json.dumps(data,indent=2)+"\n");print(json.dumps(data,indent=2))
    if not data["preflight_validation_pass"]:raise SystemExit(1)


def flush_batch(bid,results):
    d=OUT/f"checkpoints/batches/batch_{bid:06d}";d.mkdir(parents=True,exist_ok=False)
    for key,fields in TABLE_FIELDS.items():write_tsv(d/f"{key}.tsv.gz",[r for x in results for r in x[key]],fields,True)
    ids=[x["entries"][0]["pdb_id"] for x in results];(d/"complete.json").write_text(json.dumps({"batch_id":bid,"pdb_ids":ids,"entry_count":len(ids),"completed_at":utc()})+"\n")


def run_full(workers,batch_size):
    pre=json.loads((OUT/"preflight/filter_2_preflight_summary.json").read_text());
    if not pre["preflight_validation_pass"]:raise SystemExit("preflight failed")
    load_globals();items=input_items();completed=set();bids=[]
    for p in (OUT/"checkpoints/batches").glob("batch_*/complete.json"):
        x=json.loads(p.read_text());completed.update(x["pdb_ids"]);bids.append(x["batch_id"])
    pending=[x for x in items if x[0] not in completed];start=time.time();start_iso=utc();bid=max(bids,default=-1)+1;buf=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        it=iter(pending);futs={}
        for _ in range(workers*2):
            try:x=next(it);futs[pool.submit(parse_entry,x)]=x
            except StopIteration:break
        while futs:
            done,_=wait(futs,return_when=FIRST_COMPLETED)
            for f in done:
                item=futs.pop(f)
                try:buf.append(f.result())
                except Exception:buf.append(parse_entry(item))
                try:x=next(it);futs[pool.submit(parse_entry,x)]=x
                except StopIteration:pass
                if len(buf)>=batch_size:
                    flush_batch(bid,buf);completed.update(x["entries"][0]["pdb_id"] for x in buf);buf=[];bid+=1;tmp=OUT/"checkpoints/progress.json.tmp";tmp.write_text(json.dumps({"status":"running","processed":len(completed),"total":len(items),"start":start_iso,"workers":workers,"updated":utc(),"elapsed_seconds":round(time.time()-start,2)},indent=2)+"\n");os.replace(tmp,OUT/"checkpoints/progress.json")
        if buf:flush_batch(bid,buf);completed.update(x["entries"][0]["pdb_id"] for x in buf)
    progress={"status":"completed","processed":len(completed),"total":len(items),"start":start_iso,"end":utc(),"workers":workers,"elapsed_seconds":round(time.time()-start,2)};(OUT/"checkpoints/progress.json").write_text(json.dumps(progress,indent=2)+"\n");print(json.dumps(progress,indent=2))


def merge():
    batches=sorted((OUT/"checkpoints/batches").glob("batch_*"))
    for key,fields in TABLE_FIELDS.items():
        target=OUT/f"full/filter_2_{key}.tsv.gz";tmp=target.with_suffix(".gz.tmp")
        with gzip.open(tmp,"wt",encoding="utf-8",newline="") as oh:
            w=csv.DictWriter(oh,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader()
            for b in batches:
                for r in iter_tsv(b/f"{key}.tsv.gz"):w.writerow(r)
        os.replace(tmp,target)


def finalize():
    merge();entries=list(iter_tsv(OUT/"full/filter_2_entries.tsv.gz"));sources=list(iter_tsv(OUT/"full/filter_2_sources.tsv.gz"));assemblies=list(iter_tsv(OUT/"full/filter_2_assemblies.tsv.gz"));conformers=list(iter_tsv(OUT/"full/filter_2_conformers.tsv.gz"));parents_path=OUT/"full/filter_2_parents.tsv.gz";connections_path=OUT/"full/filter_2_connections.tsv.gz"
    aliases={"filter_2_entry_inventory.tsv.gz":"filter_2_entries.tsv.gz","filter_2_source_component_instances.tsv.gz":"filter_2_sources.tsv.gz","filter_2_assembly_component_instances.tsv.gz":"filter_2_assemblies.tsv.gz","filter_2_component_conformers.tsv.gz":"filter_2_conformers.tsv.gz","filter_2_parent_mapping.tsv.gz":"filter_2_parents.tsv.gz","filter_2_struct_conn_links.tsv.gz":"filter_2_connections.tsv.gz"}
    for dst,src in aliases.items():shutil.copy2(OUT/"full"/src,OUT/"full"/dst)
    used_ids=sorted({r["label_comp_id"] for r in sources});ccd={r["original_component_id"]:r for r in iter_tsv(OUT/"references/ccd_component_cache.tsv.gz")};components=[ccd[x] for x in used_ids if x in ccd];missing=[{"original_component_id":x,"resolved_ccd_id":"","ccd_identity_status":"ccd_missing","chemical_entity_class":"unknown","artifact_prior":"unknown","classification_reason":"not_in_frozen_ccd","rule_version":RULE_VERSION} for x in used_ids if x not in ccd];components+=missing;write_tsv(OUT/"full/filter_2_component_classification.tsv.gz",components,COMP_FIELDS,True)
    ordinary=[r for r in sources if r["filter_2_route"]=="ordinary_small_molecule_candidate"];special=[r for r in sources if r["filter_2_route"].endswith("_special")];artifact=[r for r in sources if r["filter_2_route"]=="artifact_review"];poly=[r for r in sources if r["filter_2_route"]=="polymer_or_modified_residue"];unresolved=[r for r in sources if r["filter_2_route"]=="unresolved_review"]
    for name,data in [("ordinary_candidates",ordinary),("special_candidates",special),("artifact_review",artifact),("polymer_modified_residues",poly),("unresolved_review",unresolved)]:write_tsv(OUT/f"full/filter_2_{name}.tsv.gz",data,SOURCE_FIELDS,True)
    water=[{"pdb_id":r["pdb_id"],"water_count":r["water_count"]} for r in entries if int(r["water_count"])>0];write_tsv(OUT/"full/filter_2_excluded_water_summary.tsv.gz",water,["pdb_id","water_count"],True)
    ord_assembly=[r for r in assemblies if r["filter_2_route"]=="ordinary_small_molecule_candidate"]
    for name,data,fields in [("filter_2_ordinary_component_instances.tsv.gz",ordinary,SOURCE_FIELDS),("filter_2_ordinary_assembly_component_instances.tsv.gz",ord_assembly,ASSEMBLY_FIELDS),("filter_2_special_component_instances.tsv.gz",special,SOURCE_FIELDS),("filter_2_artifact_review.tsv.gz",artifact,SOURCE_FIELDS),("filter_2_unresolved_review.tsv.gz",unresolved,SOURCE_FIELDS),("filter_2_component_classification.tsv.gz",components,COMP_FIELDS)]:write_tsv(OUT/"release"/name,data,fields,True)
    # Reports with route-level distinct units.
    def report(name,rows_,key):
        d=defaultdict(lambda:{"source_component_instance_count":0,"assembly_component_instance_count":0,"ccd":set(),"pdb":set(),"assembly":set()})
        for r in rows_:
            k=r[key];d[k]["source_component_instance_count"]+=1;d[k]["ccd"].add(r.get("resolved_ccd_id",""));d[k]["pdb"].add(r["pdb_id"])
        for r in assemblies:
            k=r[key] if key in r else r["filter_2_route"];d[k]["assembly_component_instance_count"]+=1;d[k]["assembly"].add(r["pdb_id"]+"|"+r["assembly_id"])
        out=[]
        for k,v in sorted(d.items()):out.append({key:k,"source_component_instance_count":v["source_component_instance_count"],"assembly_component_instance_count":v["assembly_component_instance_count"],"unique_ccd_count":len(v["ccd"]-{""}),"unique_pdb_entry_count":len(v["pdb"]),"unique_qualified_assembly_count":len(v["assembly"])})
        write_tsv(OUT/"reports"/name,out,[key,"source_component_instance_count","assembly_component_instance_count","unique_ccd_count","unique_pdb_entry_count","unique_qualified_assembly_count"])
    report("filter_2_route_distribution.tsv",sources,"filter_2_route")
    simple=[("filter_2_entry_flow.tsv",entries,"entry_status"),("filter_2_artifact_prior_distribution.tsv",sources,"artifact_prior"),("filter_2_polymer_context_distribution.tsv",sources,"polymer_context"),("filter_2_covalent_status_distribution.tsv",sources,"instance_covalent_link_status"),("filter_2_metal_status_distribution.tsv",sources,"instance_metal_link_status"),("filter_2_assembly_membership_distribution.tsv",assemblies,"assembly_membership_status")]
    for name,data,key in simple:write_tsv(OUT/"reports"/name,[{key:k,"count":v} for k,v in sorted(Counter(r[key] for r in data).items())],[key,"count"])
    for name,key in [("filter_2_component_class_distribution.tsv","chemical_entity_class"),("filter_2_ccd_status_distribution.tsv","ccd_identity_status"),("filter_2_rdkit_status_distribution.tsv","rdkit_parse_status"),("filter_2_element_distribution.tsv","element_set"),("filter_2_fragment_distribution.tsv","fragment_count")]:write_tsv(OUT/"reports"/name,[{key:k,"count":v} for k,v in sorted(Counter(r.get(key,"") for r in components).items())],[key,"count"])
    write_tsv(OUT/"reports/filter_2_failure_reason_distribution.tsv",[{"terminal_reason":k,"count":v} for k,v in sorted(Counter(r["terminal_reason"] for r in entries if r["entry_status"]!="pass").items())],["terminal_reason","count"])
    # Historical crosswalk by pdb/component/chain/residue, audit only.
    historical_path=Path('/root/autodl-tmp/vs_benchmark/data_interaction_refinement_arpeggio_v2/main_benchmark_candidate_v2.csv')
    if not historical_path.exists():
        candidates=list(Path('/root/autodl-tmp/vs_benchmark').rglob('candidate_pairs_classified_v1.csv'));historical_path=candidates[0] if candidates else Path('/nonexistent')
    hist_rows=[]
    if historical_path.exists():
        op=gzip.open if str(historical_path).endswith('.gz') else open
        with op(historical_path,'rt',encoding='utf-8',newline='') as h:
            for r in csv.DictReader(h,delimiter=',' if historical_path.suffix=='.csv' else '\t'):hist_rows.append(r)
    new_by=defaultdict(list)
    for r in sources:new_by[(r['pdb_id'],r['label_comp_id'])].append(r)
    assembly_routes_by=defaultdict(set)
    for r in assemblies:
        assembly_routes_by[(r['pdb_id'],r['resolved_ccd_id'])].add(r['filter_2_route'])
    cross=[]
    for h in hist_rows:
        pid=(h.get('pdb_id') or '').lower();cid=(h.get('ligand_id') or h.get('component_id') or '').upper();matches=new_by.get((pid,cid),[]);routes=sorted({r['filter_2_route'] for r in matches});cross.append({"pdb_id":pid,"component_id":cid,"historical_candidate_status":h.get('preliminary_category','historical_candidate'),"new_source_component_route":','.join(routes),"new_assembly_component_route":','.join(sorted(assembly_routes_by.get((pid,cid),set()))),"route_agreement":"not_directly_comparable" if not routes else "mapped","route_disagreement":"" if routes else "not_found","disagreement_reason":"instance_model_changed_or_not_in_filter1_scope" if not routes else ""})
    write_tsv(OUT/"reports/filter_2_historical_crosswalk.tsv",cross,["pdb_id","component_id","historical_candidate_status","new_source_component_route","new_assembly_component_route","route_agreement","route_disagreement","disagreement_reason"])
    # SQLite duplicate validation.
    db=sqlite3.connect(OUT/"validation/filter_2_keys.sqlite");dups={}
    for name,data,key in [("component",components,"original_component_id"),("source",sources,"source_component_instance_id"),("assembly",assemblies,"assembly_component_instance_id"),("conformer",conformers,"component_conformer_id")]:
        db.execute(f"DROP TABLE IF EXISTS {name}");db.execute(f"CREATE TABLE {name}(k TEXT PRIMARY KEY)");dup=0
        for r in data:
            try:db.execute(f"INSERT INTO {name} VALUES(?)",(r[key],))
            except sqlite3.IntegrityError:dup+=1
        db.commit();dups[name]=dup
    db.close()
    contamination={"ordinary_polymer_residue":sum(r["polymer_context"]!="independent_nonpolymer" for r in ordinary),"ordinary_modified_polymer_residue":sum(r["modified_residue_status"]=="true" for r in ordinary),"ordinary_short_peptide":sum(r["polymer_context"]=="short_peptide" for r in ordinary),"ordinary_DNA_RNA":sum(r["polymer_context"] in {"rna_polymer","dna_polymer","hybrid_nucleic_acid"} for r in ordinary),"ordinary_branched_glycan":sum(r["polymer_context"]=="branched_glycan" for r in ordinary),"ordinary_water":sum(r["chemical_entity_class"]=="water" for r in ordinary),"ordinary_metal_inorganic":sum(r["chemical_entity_class"] in {"metal_or_inorganic","organometallic"} for r in ordinary),"ordinary_unresolved_CCD":sum(r["instance_status"]!="resolved" for r in ordinary),"ordinary_explicit_receptor_covalent":sum(r["instance_covalent_link_status"]=="declared_receptor_covalent" for r in ordinary)}
    def verify_raw(row):
        p=Path(row["mmcif_path"])
        if not p.exists():return "missing"
        if p.stat().st_size!=int(row["mmcif_file_size"]):return "size_mismatch"
        return "pass" if sha(p)==row["mmcif_checksum"] else "checksum_mismatch"
    raw_rows=list(iter_tsv(OUT/"inputs/processing_1_mmcif_index_snapshot.tsv.gz"))
    with ThreadPoolExecutor(max_workers=16) as pool:raw_checks=Counter(pool.map(verify_raw,raw_rows,chunksize=64))
    raw_audit={"checked":len(raw_rows),**dict(raw_checks)}
    (OUT/"validation/raw_mmcif_checksum_audit.json").write_text(json.dumps(raw_audit,indent=2)+"\n")
    raw_mismatch=raw_checks["missing"]+raw_checks["size_mismatch"]+raw_checks["checksum_mismatch"]
    progress=json.loads((OUT/"checkpoints/progress.json").read_text());validation={"input_entries":248037,"input_qualified_assemblies":360611,"entry_inventory_rows":len(entries),"parse_success":sum(r["parse_status"]=="success" for r in entries),"parse_failed":sum(r["parse_status"]=="failed" for r in entries),"duplicate_component_identity_key":dups["component"],"duplicate_source_component_instance_key":dups["source"],"duplicate_assembly_component_instance_key":dups["assembly"],"duplicate_conformer_key":dups["conformer"],"missing_component_route":sum(not r["filter_2_route"] for r in sources),"missing_terminal_status":sum(not r["entry_status"] for r in entries),"silent_drop":248037-len(entries),"ordinary_contamination":contamination,"processing_1_modified":sha(P1/"release/processing_1_mmcif_index.tsv.gz")!=json.loads((OUT/"inputs/input_checksums.json").read_text())[str(OUT/"inputs/processing_1_mmcif_index_snapshot.tsv.gz")],"filter_1_modified":sha(ENTRY_INPUT)!=json.loads((OUT/"inputs/input_checksums.json").read_text())[str(OUT/"inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")],"raw_mmcif_modified":raw_mismatch>0,"raw_mmcif_checksum_audit":raw_audit,"historical_directories_modified":False,"assembly_coordinate_materialization_started":False,"pair_construction_started":False,"distance_calculation_started":False,"interaction_annotation_started":False,"structure_quality_filtering_started":False,"checksum_mismatch":raw_mismatch}
    # Correct source-vs-snapshot comparisons (snapshots are exact copies).
    validation["processing_1_modified"]=sha(P1/"release/processing_1_mmcif_index.tsv.gz")!=sha(OUT/"inputs/processing_1_mmcif_index_snapshot.tsv.gz");validation["filter_1_modified"]=sha(ENTRY_INPUT)!=sha(OUT/"inputs/filter_1_receptor_qualified_entries_snapshot.tsv.gz")
    validation["release_validation_pass"]=len(entries)==248037 and validation["silent_drop"]==0 and all(v==0 for v in dups.values()) and validation["missing_component_route"]==0 and validation["missing_terminal_status"]==0 and all(v==0 for v in contamination.values()) and not validation["processing_1_modified"] and not validation["filter_1_modified"] and not validation["raw_mmcif_modified"]
    summary={"full_start":progress["start"],"full_end":progress["end"],"runtime_seconds":progress["elapsed_seconds"],"input_entries":248037,"input_qualified_assemblies":360611,"parse_success":validation["parse_success"],"parse_failed":validation["parse_failed"],"unique_ccd_count":len(components),"source_component_instance_count":len(sources),"assembly_component_instance_count":len(assemblies),"conformer_count":len(conformers),"ordinary_source_instance_count":len(ordinary),"ordinary_assembly_instance_count":len(ord_assembly),"special_instance_count":len(special),"artifact_review_count":len(artifact),"polymer_modified_count":len(poly),"unresolved_count":len(unresolved),"water_exclusion_count":sum(int(r["water_count"]) for r in entries),"covalent_linked_count":sum(r["instance_covalent_link_status"]=="declared_receptor_covalent" for r in sources),"metal_inorganic_count":sum(r["chemical_entity_class"] in {"metal_or_inorganic","organometallic"} for r in sources),"rdkit_parse_failure_count":sum(r.get("rdkit_parse_status")!="pass" for r in components),"ccd_missing_invalid_count":sum(r.get("ccd_identity_status") not in {"exact_ccd_match","obsolete_id_resolved"} for r in components),"historical_crosswalk_rows":len(cross),"validation":validation}
    (OUT/"reports/filter_2_final_summary.json").write_text(json.dumps(summary,indent=2)+"\n");(OUT/"validation/filter_2_release_validation.json").write_text(json.dumps(validation,indent=2)+"\n");(OUT/"release/filter_2_release_validation.json").write_text(json.dumps(validation,indent=2)+"\n");(OUT/"release/filter_2_release_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    interface={"project_name":"Benchmark 1.0","filter_name":"Filter 2 - Ligand Instance Identification and Chemical-Scope Qualification","filter_version":"1.0","input_entry_count":248037,"input_qualified_assembly_count":360611,"source_component_instance_count":len(sources),"assembly_component_instance_count":len(assemblies),"unique_component_id_count":len(components),"ordinary_source_instance_count":len(ordinary),"ordinary_assembly_instance_count":len(ord_assembly),"special_instance_count":len(special),"artifact_review_count":len(artifact),"polymer_modified_count":len(poly),"unresolved_count":len(unresolved),"ccd_snapshot_version":"Sat, 11 Jul 2026 12:01:19 GMT","artifact_reference_versions":["refinement_v2_provisional_20260719","official_reference_unavailable"],"rule_version":RULE_VERSION,"release_creation_time":utc(),"release_validation_pass":validation["release_validation_pass"]};(OUT/"release/filter_2_downstream_interface.json").write_text(json.dumps(interface,indent=2)+"\n")
    files=[p for p in (OUT/"release").iterdir() if p.is_file() and p.name!="SHA256SUMS"];(OUT/"release/SHA256SUMS").write_text("".join(f"{sha(p)}  {p.name}\n" for p in sorted(files)))
    (OUT/"provenance/filter_2_run_provenance.json").write_text(json.dumps({"host":platform.node(),"python":sys.version,"gemmi":gemmi.__version__,"rdkit":rdBase.rdkitVersion,"start":progress["start"],"end":progress["end"],"workers":progress["workers"],"config_sha256":sha(OUT/"configs/filter_2.yaml"),"release_validation_pass":validation["release_validation_pass"]},indent=2)+"\n")
    print(json.dumps(summary,indent=2));
    if not validation["release_validation_pass"]:raise SystemExit(1)


def validate():
    p=OUT/"release/filter_2_release_validation.json";data=json.loads(p.read_text());print(json.dumps(data,indent=2));raise SystemExit(0 if data["release_validation_pass"] else 1)


def main():
    ap=argparse.ArgumentParser(description="Benchmark 1.0 Filter 2 pipeline");sub=ap.add_subparsers(dest="cmd",required=True)
    for c in ["setup","audit","prepare-references","preflight","finalize","validate"]:sub.add_parser(c)
    f=sub.add_parser("full");f.add_argument("--workers",type=int,default=16);f.add_argument("--batch-size",type=int,default=200)
    a=ap.parse_args()
    if a.cmd=="setup":setup()
    elif a.cmd=="audit":audit()
    elif a.cmd=="prepare-references":prepare_references()
    elif a.cmd=="preflight":run_preflight()
    elif a.cmd=="full":run_full(a.workers,a.batch_size)
    elif a.cmd=="finalize":finalize()
    elif a.cmd=="validate":validate()
if __name__=="__main__":main()
