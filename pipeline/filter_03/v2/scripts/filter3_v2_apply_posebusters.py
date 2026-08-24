#!/usr/bin/env python3
from __future__ import annotations
import argparse, concurrent.futures, json, os, time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import yaml

def clean(v):
    if v is None: return ''
    t=str(v).strip(); return '' if t.lower() in {'','.','?','none','false','nan'} else t
def bv(v): return False if v is None or pd.isna(v) else bool(v)
def atomic_json(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)
def write(frame,path):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+'.tmp'); pq.write_table(pa.Table.from_pandas(frame,preserve_index=False),tmp,compression='zstd'); os.replace(tmp,path)
def read_bucket(root,bucket):
    paths=[str(p) for p in sorted((Path(root)/f'bucket_id={bucket:03d}').glob('*.parquet'))]
    return pd.DataFrame() if not paths else ds.dataset(paths,format='parquet').to_table().to_pandas(split_blocks=True,self_destruct=True)
def utc(): return datetime.now(timezone.utc).isoformat()

def process(bucket,config,run_dir):
    run_dir=Path(run_dir); out=run_dir/'work/final_batches'/f'bucket_id={bucket:03d}'; marker=out/'_COMPLETE.json'
    if marker.exists(): return json.loads(marker.read_text())
    pre=pq.ParquetFile(run_dir/'work/preclassification_batches'/f'bucket_id={bucket:03d}'/'pair_quality_pre_posebusters.parquet').read().to_pandas()
    old=read_bucket(Path(config['input']['reusable_evidence_output'])/'posebusters_raw_geometry',bucket)
    new=read_bucket(run_dir/'work/posebusters_new_batches',bucket)
    combined=pd.concat([old,new],ignore_index=True,sort=False) if not old.empty or not new.empty else pd.DataFrame()
    pb={r['source_ligand_instance_id']:r for r in combined.to_dict('records')} if not combined.empty else {}
    ligand=read_bucket(Path(config['input']['reusable_evidence_output'])/'ligand_validation_mapping',bucket)
    source_by_placement=dict(zip(ligand['ligand_assembly_placement_id'],ligand['source_ligand_instance_id']))
    rows=[]
    for row in pre.to_dict('records'):
        if row['filter3_v2_terminal_status']!='FILTER3_POSEBUSTERS_PENDING': rows.append(row); continue
        source=source_by_placement.get(row['ligand_assembly_placement_id'],''); evidence=pb.get(source)
        warnings=[x for x in clean(row.get('warning_codes')).split(';') if x]
        if evidence is None or clean(evidence.get('posebusters_status'))!='COMPLETED':
            row['filter3_v2_terminal_status']='FILTER3_TECHNICAL_FAILURE'; row['decision']='FAIL'; row['destination']='technical_failure'; row['reason_codes']='POSEBUSTERS_EXECUTION_FAILURE'
        else:
            fatal_chem=any(not bv(evidence.get(name)) for name in ('sanitization','all_atoms_connected','no_radicals'))
            fatal_clash=not bv(evidence.get('internal_steric_clash'))
            pbwarn=any(not bv(evidence.get(name)) for name in config['posebusters']['warning_checks'])
            reasons=[]
            if fatal_chem: reasons.append('POSEBUSTERS_FATAL_CHEMISTRY')
            if fatal_clash: reasons.append('POSEBUSTERS_FATAL_INTERNAL_CLASH')
            if pbwarn: warnings.append('POSEBUSTERS_NONFATAL_WARNING')
            if reasons:
                row['filter3_v2_terminal_status']='FILTER3_REJECT'; row['decision']='REJECT'; row['destination']='excluded'; row['reason_codes']=';'.join(reasons)
            else:
                high=(float(row['entry_resolution'])<=2.5 and row.get('ligand_mean_occupancy') is not None and float(row['ligand_mean_occupancy'])>=0.8 and row.get('pocket_mean_occupancy') is not None and float(row['pocket_mean_occupancy'])>=0.8 and not warnings)
                row['filter3_v2_terminal_status']='FILTER3_HIGH_QUALITY' if high else 'FILTER3_GOOD_QUALITY'; row['decision']='PASS'; row['destination']='ordinary_ground_truth'; row['reason_codes']=''
            row['warning_codes']=';'.join(sorted(set(warnings)))
        rows.append(row)
    frame=pd.DataFrame(rows); write(frame,out/'pair_quality_v2.parquet')
    result={'status':'COMPLETED','bucket_id':bucket,'pair_count':len(frame),'terminal_status_counts':dict(Counter(frame['filter3_v2_terminal_status'])),'finished_at':utc()}; atomic_json(marker,result); return result

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); ap.add_argument('--run-dir',required=True); ap.add_argument('--workers',type=int,default=8); args=ap.parse_args(); config=yaml.safe_load(Path(args.config).read_text()); started=time.time(); results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(process,b,config,args.run_dir) for b in range(256)]
        for f in concurrent.futures.as_completed(futures):
            r=f.result(); results.append(r); counts=Counter(); [counts.update(x['terminal_status_counts']) for x in results]
            progress={'status':'RUNNING','phase':'FILTER3_V2_FINAL_CLASSIFICATION','bucket_completed':len(results),'bucket_total':256,'pair_count':sum(x['pair_count'] for x in results),'terminal_status_counts':dict(counts),'runtime_seconds':time.time()-started,'updated_at':utc()}; atomic_json(Path(args.run_dir)/'final_status.json',progress); print(json.dumps(progress),flush=True)
    progress['status']='COMPLETED'; progress['phase']='FILTER3_V2_FINAL_CLASSIFICATION_COMPLETE'; atomic_json(Path(args.run_dir)/'final_status.json',progress)
if __name__=='__main__': main()
