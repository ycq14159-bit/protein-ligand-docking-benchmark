from __future__ import annotations
import argparse, hashlib, json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''): h.update(b)
    return h.hexdigest()

def classify(r):
    b=r.overlap_class; ql=r.q_match_level; qr=r.q_relevance_class
    q_exact=ql=='EXACT_ASSEMBLY_CCD_CHAIN'
    if b=='EXACT' and q_exact and qr=='RELEVANT': return 'LIGAND_RELEVANCE_CONFIRMED',True
    if b in {'EXACT','AMBIGUOUS'} and q_exact and qr=='IRRELEVANT': return 'RELEVANCE_CONFLICT',False
    if b=='AMBIGUOUS' and q_exact and qr=='RELEVANT': return 'BIOLIP_MAPPING_AMBIGUOUS',False
    if b=='NOT_RETAINED' and q_exact and qr=='RELEVANT': return 'SINGLE_SOURCE_RELEVANT_Q',False
    if b=='NOT_ASSESSABLE' and q_exact and qr=='RELEVANT': return 'SINGLE_SOURCE_RELEVANT_Q_BIOLIP_UNAVAILABLE',False
    if b=='NOT_RETAINED' and q_exact and qr=='IRRELEVANT': return 'EXTERNAL_IRRELEVANCE_SUPPORTED',False
    if b=='NOT_ASSESSABLE' and q_exact and qr=='IRRELEVANT': return 'Q_IRRELEVANCE_SUPPORTED_BIOLIP_UNAVAILABLE',False
    if qr=='MIXED_RELEVANT_IRRELEVANT' or not q_exact: return 'RELEVANCE_AMBIGUOUS',False
    return 'RELEVANCE_NOT_ASSESSABLE',False

def main(root:Path,audit:Path,out:Path,entry_scope:Path|None):
    out.mkdir(parents=True,exist_ok=True)
    p1=root/'processing_01_source_and_entry_quality_v2/runs/20260820_full_01/output/entry_validation_and_qualification.parquet'
    f1=root/'filter_1_protein_receptor_qualification/full/filter_1_entries.tsv.gz'
    f2=root/'filter_2_ligand_qualification_v3/runs/20260804_full_01/output'
    af=audit/'06_combined/two_audit_source_comparison.tsv.gz'
    required=[f2/'provisional_source_ligands.tsv.gz',f2/'ligand_assembly_logical_placements.tsv.gz',af]
    if entry_scope is None: required += [p1,f1]
    else: required += [entry_scope]
    for p in required:
        if not p.exists(): raise FileNotFoundError(p)
    if entry_scope is None:
        e=pq.read_table(p1,columns=['pdb_id','entry_qualification_status']).to_pandas()
        eligible=set(e.loc[e.entry_qualification_status.eq('ENTRY_XRAY_ELIGIBLE'),'pdb_id'].astype(str).str.lower())
        f=pd.read_csv(f1,sep='\t',usecols=['pdb_id','entry_receptor_pass'])
        f.pdb_id=f.pdb_id.astype(str).str.lower()
        entries=set(f.loc[f.entry_receptor_pass.astype(bool)&f.pdb_id.isin(eligible),'pdb_id'])
        if len(entries)!=142049: raise RuntimeError(f'entry scope {len(entries)} != 142049')
    else:
        es=pd.read_parquet(entry_scope);entries=set(es.pdb_id.astype(str).str.lower())
        if len(entries)!=106842: raise RuntimeError(f'source-containing recovered scope {len(entries)} != 106842')

    src=pd.read_csv(f2/'provisional_source_ligands.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    src.pdb_id=src.pdb_id.str.lower(); src=src[src.pdb_id.isin(entries)].copy()
    if len(src)!=533610 or src.source_ligand_instance_id.nunique()!=533610: raise RuntimeError(f'source scope mismatch {len(src)}')
    a=pd.read_csv(af,sep='\t',dtype=str,keep_default_na=False)
    if len(a)!=487982 or a.source_ligand_instance_id.nunique()!=487982: raise RuntimeError('audit universe mismatch')
    decision=[classify(r) for r in a.itertuples(index=False)]
    a['dual_source_status']=[x[0] for x in decision]; a['dual_source_pass']=[x[1] for x in decision]
    slim=a[['source_ligand_instance_id','overlap_class','q_match_level','q_relevance_class','dual_source_status','dual_source_pass']]
    m=src.merge(slim,on='source_ligand_instance_id',how='left',validate='one_to_one',indicator=True)
    m['is_suspicious']=m._merge.eq('both')
    m.loc[~m.is_suspicious,'dual_source_status']='NOT_SUSPICIOUS_BYPASS'
    m.loc[~m.is_suspicious,'dual_source_pass']=True
    m['dual_source_pass']=m.dual_source_pass.astype(bool)
    m.drop(columns=['_merge'],inplace=True)

    pl=pd.read_csv(f2/'ligand_assembly_logical_placements.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    pl.pdb_id=pl.pdb_id.str.lower(); pl=pl[pl.pdb_id.isin(entries)].copy()
    if len(pl)!=718007 or pl.assembly_ligand_placement_id.nunique()!=718007: raise RuntimeError(f'placement scope mismatch {len(pl)}')
    retained=m[m.dual_source_pass].copy(); keep=set(retained.source_ligand_instance_id)
    rpl=pl[pl.source_ligand_instance_id.isin(keep)].copy()
    mapped_counts=rpl.groupby('source_ligand_instance_id').size()
    mapped_sources=set(mapped_counts.index)
    unmapped=retained[~retained.source_ligand_instance_id.isin(mapped_sources)].copy()
    formal_sources=retained[retained.source_ligand_instance_id.isin(mapped_sources)].copy()
    additional=int((mapped_counts-1).clip(lower=0).sum())
    formal_total=len(formal_sources)+additional
    if formal_total != len(rpl): raise RuntimeError('formal total identity failed')

    m.to_csv(out/'01_source_membership.tsv.gz',sep='\t',index=False,compression='gzip')
    rpl.to_csv(out/'02_retained_assembly_placements.tsv.gz',sep='\t',index=False,compression='gzip')
    unmapped.to_csv(out/'03_retained_source_without_assembly_mapping.tsv.gz',sep='\t',index=False,compression='gzip')
    formal_sources.to_csv(out/'04_formal_candidate_sources.tsv.gz',sep='\t',index=False,compression='gzip')
    counts=Counter(m.dual_source_status)
    summary={
      'run_type':'VERSIONED_READ_ONLY_MEMBERSHIP_PROJECTION_NO_FROZEN_SCIENTIFIC_INPUT_MODIFIED',
      'created_utc':datetime.now(timezone.utc).isoformat(),
      'formal_entry_scope':142049,'source_containing_entry_scope':len(entries),'pre_rule_candidate_source_ligands':len(src),'pre_rule_assembly_placement_rows':len(pl),
      'suspicious_in_formal_scope':int(m.is_suspicious.sum()),
      'suspicious_strict_confirmed_in_formal_scope':int((m.is_suspicious&m.dual_source_pass).sum()),
      'non_suspicious_bypass_sources':int((~m.is_suspicious).sum()),
      'dual_source_membership_pass_sources_before_assembly_mapping':len(retained),
      'formal_output_candidate_source_ligands':len(formal_sources),
      'formal_output_additional_assembly_ligands':additional,
      'formal_output_total_ligand_records':formal_total,
      'retained_assembly_placement_rows':len(rpl),'retained_sources_with_mapping':len(mapped_sources),
      'retained_sources_without_mapping':len(unmapped),
      'stopped_suspicious_sources':int((m.is_suspicious&~m.dual_source_pass).sum()),
      'status_counts':dict(counts),
      'validation':{
        'formal_entry_scope_142049':True,'source_containing_entry_scope_106842':len(entries)==106842 if entry_scope is not None else True,'source_scope_533610':len(src)==533610,
        'placement_scope_718007':len(pl)==718007,'source_pk_unique':m.source_ligand_instance_id.nunique()==len(m),
        'placement_pk_unique':pl.assembly_ligand_placement_id.nunique()==len(pl),
        'retained_placement_fk_missing':int((~rpl.source_ligand_instance_id.isin(keep)).sum()),
        'source_decision_closure':int(m.dual_source_status.value_counts().sum())==len(m),
        'formal_total_identity':formal_total==len(formal_sources)+additional==len(rpl)
      },
      'input_sha256':{str(p.relative_to(root)):sha256(p) for p in ([f2/'provisional_source_ligands.tsv.gz',f2/'ligand_assembly_logical_placements.tsv.gz'] + ([p1,f1] if entry_scope is None else [entry_scope]))},
      'audit_input_sha256':sha256(af)
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    pd.DataFrame(sorted(counts.items()),columns=['dual_source_status','source_ligands']).to_csv(out/'status_counts.tsv',sep='\t',index=False)
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--entry-scope',type=Path);x=ap.parse_args();main(x.root,x.audit,x.out,x.entry_scope)
