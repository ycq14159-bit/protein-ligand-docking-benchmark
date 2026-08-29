from __future__ import annotations
import argparse,hashlib,json,shutil
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def main(run:Path):
 out=run/'output';rel=run/'release';rel.mkdir(exist_ok=True)
 m=pd.read_csv(out/'01_source_membership.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
 s=pd.read_csv(out/'04_formal_candidate_sources.tsv.gz',sep='\t',usecols=['source_ligand_instance_id'],dtype=str)
 p=pd.read_csv(out/'02_retained_assembly_placements.tsv.gz',sep='\t',usecols=['source_ligand_instance_id','assembly_ligand_placement_id'],dtype=str)
 u=pd.read_csv(out/'03_retained_source_without_assembly_mapping.tsv.gz',sep='\t',dtype=str)
 sid=set(s.source_ligand_instance_id);ps=set(p.source_ligand_instance_id)
 validation={
  'membership_rows':len(m),'membership_source_pk_unique':m.source_ligand_instance_id.nunique()==len(m),
  'formal_source_rows':len(s),'formal_source_pk_unique':s.source_ligand_instance_id.nunique()==len(s),
  'formal_placement_rows':len(p),'formal_placement_pk_unique':p.assembly_ligand_placement_id.nunique()==len(p),
  'placement_source_fk_missing':len(ps-sid),'formal_source_without_placement':len(sid-ps),
  'additional_assembly_copies':len(p)-len(s),'formal_total_identity':len(p)==len(s)+(len(p)-len(s)),
  'membership_pass_before_mapping':int(m.dual_source_pass.eq('True').sum()),'assembly_mapping_review_sources':len(u),
  'suspicious_rows':int(m.is_suspicious.eq('True').sum()),
  'strict_confirmed_suspicious':int((m.is_suspicious.eq('True')&m.dual_source_status.eq('LIGAND_RELEVANCE_CONFIRMED')).sum()),
 }
 expected={'membership_rows':533610,'formal_source_rows':183904,'formal_placement_rows':236383,'additional_assembly_copies':52479,'membership_pass_before_mapping':183909,'assembly_mapping_review_sources':5,'suspicious_rows':354325,'strict_confirmed_suspicious':4624}
 validation['expected_counts']=expected;validation['counts_match']=all(validation[k]==v for k,v in expected.items());validation['pass']=validation['counts_match'] and validation['membership_source_pk_unique'] and validation['formal_source_pk_unique'] and validation['formal_placement_pk_unique'] and validation['placement_source_fk_missing']==0 and validation['formal_source_without_placement']==0
 if not validation['pass']:raise RuntimeError(validation)
 (rel/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
 shutil.copy2(out/'summary.json',rel/'release_summary.json');shutil.copy2(out/'status_counts.tsv',rel/'status_counts.tsv');shutil.copy2(run/'reconstruction/recovered_scope_validation.json',rel/'scope_recovery_validation.json')
 files=[out/'01_source_membership.tsv.gz',out/'02_retained_assembly_placements.tsv.gz',out/'03_retained_source_without_assembly_mapping.tsv.gz',out/'04_formal_candidate_sources.tsv.gz',rel/'validation.json',rel/'release_summary.json',rel/'status_counts.tsv',rel/'scope_recovery_validation.json']
 rows=[]
 for f in files:rows.append({'relative_path':str(f.relative_to(run)),'size_bytes':f.stat().st_size,'sha256':sha(f),'role':'scientific_output' if 'output/' in str(f.relative_to(run)) else 'release_control'})
 pd.DataFrame(rows).to_csv(rel/'output_manifest.tsv',sep='\t',index=False)
 sums=''.join(f"{r['sha256']}  {r['relative_path']}\n" for r in rows);(rel/'SHA256SUMS').write_text(sums)
 marker={'status':'FROZEN','run_id':run.name,'stage':'filter_2_ligand_qualification_v4_dual_source_strict','frozen_at_utc':datetime.now(timezone.utc).isoformat(),'scientific_rule':'SUSPICIOUS only passes when BioLiP EXACT and Q-BioLiP EXACT_ASSEMBLY_CCD_CHAIN and Relevant=yes; non-suspicious sources bypass; assembly-unmapped retained sources are REVIEW','formal_output':{'candidate_source_ligands':183904,'additional_assembly_ligands':52479,'total_ligand_records':236383},'membership_pass_before_mapping':183909,'assembly_mapping_review':5,'sha256sums_sha256':sha(rel/'SHA256SUMS'),'validation_sha256':sha(rel/'validation.json'),'scope_recovery_caveat':'Original P1 v2 entry parquet was absent after migration. Source-containing scope was reconstructed from frozen mmCIF plus wwPDB validation XML and required exact closure to historical 106842/533610/718007 anchors. Current XML returned 403 for 6o8v,6xok,6xrv; those three use frozen mmCIF pass plus historical closure. 8ymg was explicitly excluded because XML metrics are missing.'}
 (run/'_FROZEN.json').write_text(json.dumps(marker,indent=2)+'\n')
 current=run.parents[1]/'CURRENT_RUN.json';current.write_text(json.dumps({'current_run':run.name,'path':str(run),'status':'FROZEN','updated_at_utc':marker['frozen_at_utc']},indent=2)+'\n')
 print(json.dumps({'validation':validation,'marker':marker},indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);x=a.parse_args();main(x.run)
