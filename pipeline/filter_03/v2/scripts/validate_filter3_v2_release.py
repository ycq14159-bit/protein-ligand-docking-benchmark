#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

ALLOWED={'FILTER3_HIGH_QUALITY','FILTER3_GOOD_QUALITY','FILTER3_REJECT','FILTER3_VALIDATION_DATA_UNAVAILABLE','FILTER3_NON_XRAY_PROTOCOL_PENDING','FILTER3_TECHNICAL_FAILURE'}
def sha(path):
 h=hashlib.sha256();
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''): h.update(chunk)
 return h.hexdigest()
def atomic(path,value):
 path=Path(path); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)
def write_release_hashes(release):
 paths=sorted(p for p in release.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
 lines=[f'{sha(path)}  {path.relative_to(release).as_posix()}' for path in paths]
 target=release/'SHA256SUMS'; tmp=target.with_suffix('.tmp'); tmp.write_text('\n'.join(lines)+'\n'); os.replace(tmp,target)
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',required=True); args=ap.parse_args(); run=Path(args.run_dir); files=sorted((run/'output/filter3_pair_quality_v2').rglob('*.parquet'))
 ids=set(); duplicate=0; statuses=Counter(); rows=0; missing_status=0; rsrz_true=0; p3_ids=set()
 for path in files:
  frame=pq.ParquetFile(path).read(columns=['pair_id','filter3_v2_terminal_status','rsrz_used_for_membership']).to_pandas(); rows+=len(frame); duplicate+=int(frame['pair_id'].duplicated().sum())+sum(x in ids for x in frame['pair_id']); ids.update(frame['pair_id']); statuses.update(frame['filter3_v2_terminal_status']); missing_status+=int(frame['filter3_v2_terminal_status'].isna().sum()); rsrz_true+=int(frame['rsrz_used_for_membership'].fillna(False).sum())
 p3=Path('/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/runs/20260811_full_01/output/provisional_pairs')
 for path in sorted(p3.rglob('*.parquet')): p3_ids.update(pq.ParquetFile(path).read(columns=['pair_id']).column('pair_id').to_pylist())
 baseline=json.loads((run/'audit/immutability_baseline.json').read_text()); immutable={item['path']:sha(item['path'])==item['sha256'] for item in baseline['files']}
 checks={'input_count_equals_output_count':rows==744580,'sample_id_unique':duplicate==0 and len(ids)==rows,'pair_set_equals_processing3':ids==p3_ids,'no_missing_terminal_status':missing_status==0,'terminal_status_allowed':set(statuses)<=ALLOWED,'no_posebusters_pending':'FILTER3_POSEBUSTERS_PENDING' not in statuses,'accounting_closed':sum(statuses.values())==rows,'rsrz_not_used_for_membership':rsrz_true==0,'upstream_and_v1_immutable':all(immutable.values()),'output_partition_count':len(files)==256}
 result={'run_id':run.name,'validated_at':datetime.now(timezone.utc).isoformat(),'input_count':len(p3_ids),'output_count':rows,'unique_pair_id_count':len(ids),'duplicate_pair_id_count':duplicate,'missing_pair_count':len(p3_ids-ids),'unexpected_pair_count':len(ids-p3_ids),'terminal_status_counts':dict(statuses),'checks':checks,'immutability_checks':immutable,'validation_pass':all(checks.values())}
 atomic(run/'release/filter3_v2_release_validation.json',result)
 if not result['validation_pass']: raise SystemExit(2)
 freeze={'status':'FROZEN','run_id':run.name,'stage':'filter_03_ground_truth_structure_quality_v2','frozen_at':datetime.now(timezone.utc).isoformat(),'accounting_pass':True,'schema_pass':True,'validation_pass':True,'manifest_sha256':sha(run/'release/output_manifest.tsv'),'code_version_reference':'scripts_manifest_sha256'}; atomic(run/'_FROZEN.json',freeze)
 interface=json.loads((run/'release/filter3_v2_downstream_interface.json').read_text()); interface['status']='FROZEN'; atomic(run/'release/filter3_v2_downstream_interface.json',interface)
 write_release_hashes(run/'release')
 current={'current_run_id':run.name,'status':'FROZEN','relative_path':f'runs/{run.name}','manifest_sha256':freeze['manifest_sha256'],'updated_at':freeze['frozen_at']}; atomic(run.parent.parent/'CURRENT_RUN.json',current)
if __name__=='__main__': main()
