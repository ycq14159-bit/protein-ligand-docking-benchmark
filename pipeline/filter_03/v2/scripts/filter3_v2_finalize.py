#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, os, shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

def utc(): return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''): h.update(chunk)
    return h.hexdigest()
def atomic_json(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)
def link(source,destination):
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists(): return
    try: os.link(source,destination)
    except OSError: shutil.copy2(source,destination)
def split_codes(series):
    count=Counter()
    for value in series.fillna('').astype(str):
        count.update(code for code in value.split(';') if code)
    return count

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); args=ap.parse_args(); run=Path(args.run_dir); output=run/'output'; release=run/'release'; output.mkdir(parents=True,exist_ok=True); release.mkdir(parents=True,exist_ok=True)
    sources={
      'filter3_pair_quality_v2':(run/'work/final_batches','pair_quality_v2.parquet'),
      'pocket_completeness_v2':(run/'work/preclassification_batches','pocket_completeness.parquet'),
      'binding_residue_quality_v2':(run/'work/preclassification_batches','binding_residue_quality_v2.parquet'),
      'structural_gap_audit_v2':(run/'work/preclassification_batches','structural_gap_audit_v2.parquet'),
      'posebusters_new_results':(run/'work/posebusters_new_batches','posebusters_results.parquet'),
    }
    for role,(root,name) in sources.items():
        for bucket in range(256):
            source=root/f'bucket_id={bucket:03d}'/name
            if source.exists(): link(source,output/role/f'bucket_id={bucket:03d}'/'part-000000.parquet')
    pair_files=sorted((output/'filter3_pair_quality_v2').rglob('*.parquet'))
    if len(pair_files)!=256: raise RuntimeError(f'expected 256 pair partitions, found {len(pair_files)}')
    status=Counter(); reasons=Counter(); warnings=Counter(); pair_count=0; pdbs=set(); preview=[]
    for path in pair_files:
        frame=pq.ParquetFile(path).read().to_pandas(); pair_count+=len(frame); pdbs.update(frame['pdb_id'].astype(str)); status.update(frame['filter3_v2_terminal_status']); reasons.update(split_codes(frame['reason_codes'])); warnings.update(split_codes(frame['warning_codes']))
        preview.append(frame.head(5))
    preview_frame=pd.concat(preview,ignore_index=True).sort_values('pair_id').head(1000)
    preview_frame.to_csv(release/'filter3_v2_pair_quality_preview.tsv',sep='\t',index=False)
    schema={}
    for role in sources:
        files=sorted((output/role).rglob('*.parquet'))
        if files:
            arrow=pq.ParquetFile(files[0]).schema_arrow
            schema[role]={'schema_version':'filter3_v2_schema_2.0.0','columns':[{'name':field.name,'type':str(field.type),'nullable':field.nullable} for field in arrow]}
    atomic_json(release/'output_schema.json',schema)
    crosswalk=[]
    for path in pair_files:
        frame=pq.ParquetFile(path).read(columns=['pair_id','filter3_v2_terminal_status']).to_pandas()
        old_path=Path('/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/runs/20260813_release_correction_02/output/filter3_pair_quality')/path.parent.name/'part-000000.parquet'
        old=pq.ParquetFile(old_path).read(columns=['pair_id','terminal_status']).to_pandas() if old_path.exists() else pd.DataFrame()
        if not old.empty:
            merged=frame.merge(old,on='pair_id',how='left'); crosswalk.append(merged.groupby(['terminal_status','filter3_v2_terminal_status']).size().reset_index(name='count'))
    if crosswalk:
        cw=pd.concat(crosswalk).groupby(['terminal_status','filter3_v2_terminal_status'],as_index=False)['count'].sum(); cw.to_csv(release/'v1_vs_v2_terminal_crosswalk.tsv',sep='\t',index=False)
    summary={'stage':'Filter 3 v2 - Protein-Ligand Ground-Truth Structure Quality Qualification','run_id':run.name,'completed_at':utc(),'input_pair_count':pair_count,'unique_pdb_count':len(pdbs),'terminal_status_counts':dict(status),'reason_code_counts':dict(reasons),'warning_code_counts':dict(warnings),'retained_pair_count':status['FILTER3_HIGH_QUALITY']+status['FILTER3_GOOD_QUALITY'],'rsrz_used_for_membership':False,'occupancy_is_hard_gate':False,'ligand_completeness_inherited_from_processing2':True}
    atomic_json(release/'filter3_v2_release_summary.json',summary)
    rows=[]
    generated=utc()
    for path in sorted(output.rglob('*.parquet')):
        pf=pq.ParquetFile(path); rows.append({'relative_path':str(path.relative_to(run)),'file_role':path.relative_to(output).parts[0],'file_format':'parquet','row_count':pf.metadata.num_rows,'column_count':len(pf.schema_arrow.names),'size_bytes':path.stat().st_size,'sha256':sha(path),'schema_version':'filter3_v2_schema_2.0.0','created_at':generated,'generated_by':'filter3_v2_finalize.py'})
    pd.DataFrame(rows).to_csv(release/'output_manifest.tsv',sep='\t',index=False)
    interface={'source_run_id':run.name,'status':'VALIDATED_PENDING_FREEZE','formal_pair_dataset':str(output/'filter3_pair_quality_v2'),'primary_key':'pair_id','terminal_status_field':'filter3_v2_terminal_status','retain_statuses':['FILTER3_HIGH_QUALITY','FILTER3_GOOD_QUALITY'],'created_at':utc()}
    atomic_json(release/'filter3_v2_downstream_interface.json',interface)
    with open(release/'SHA256SUMS','w') as handle:
        for path in sorted([p for p in release.iterdir() if p.is_file() and p.name!='SHA256SUMS']): handle.write(f'{sha(path)}  {path.name}\n')
if __name__=='__main__': main()
