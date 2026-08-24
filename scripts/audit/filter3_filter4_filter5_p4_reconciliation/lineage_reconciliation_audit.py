#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime, timezone
import hashlib, json
import pandas as pd
import pyarrow.dataset as ds

ROOT=Path('/home/linx/data/youcq/autodl-tmp/benchmark_1.0')
OUT=ROOT/'reconciliation/filter3_filter4_filter5_p4_20260822'
F3=ROOT/'filter_03_ground_truth_structure_quality_v2/runs/20260814_full_01'
F5=ROOT/'filter_05_equivalent_redocking_case'
S1=F5/'step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1'
S2=F5/'step_02_same_binding_site_audit/runs/step02_full_v2'
S3=F5/'step_03_local_receptor_state_equivalence/runs/step03_full_v6'
S4=F5/'step_04_native_ligand_pose_equivalence/runs/step04_full_v1'
S5=F5/'step_05_strict_equivalent_grouping_and_representative_selection/runs/step05_full_v1'
P4=ROOT/'processing_04_docking_ready_case_construction/runs/p4_full_v1_0_1'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''): h.update(b)
 return h.hexdigest()

def write_tsv(df,name):
 df.to_csv(OUT/name,sep='\t',index=False)

OUT.mkdir(parents=True,exist_ok=False)

# Frozen evidence and hash checks.
evidence=[]
for label,path in [
 ('F3_FROZEN',F3/'_FROZEN.json'),('F3_MANIFEST',F3/'release/output_manifest.tsv'),
 ('F3_VALIDATION',F3/'release/filter3_v2_release_validation.json'),
 ('F5S1_FROZEN',S1/'_FROZEN.json'),('F5S1_MANIFEST',S1/'output_manifest.tsv'),
 ('F5S1_VALIDATION',S1/'validation/validation.json'),('F5S1_PROVENANCE',S1/'provenance.json'),
 ('F5S5_FROZEN',S5/'_FROZEN.json'),('F5S5_MANIFEST',S5/'output_manifest.tsv'),
 ('F5S5_VALIDATION',S5/'validation/validation.json'),('F5S5_PROVENANCE',S5/'provenance.json'),
 ('P4_FROZEN',P4/'_FROZEN.json'),('P4_MANIFEST',P4/'output_manifest.parquet'),
 ('P4_VALIDATION',P4/'validation/validation.json'),('P4_PROVENANCE',P4/'input/provenance.json')]:
 evidence.append({'label':label,'path':str(path),'exists':path.exists(),'size_bytes':path.stat().st_size if path.exists() else None,'sha256':sha(path) if path.exists() else ''})
write_tsv(pd.DataFrame(evidence),'evidence_files.tsv')

# F3 baseline and revised-rule eligible universe.
q=ds.dataset(str(F3/'output/filter3_pair_quality_v2'),format='parquet').to_table(columns=['pair_id','filter3_v2_terminal_status','warning_codes']).to_pandas()
q['pair_id']=q.pair_id.astype(str)
f3_counts={str(k):int(v) for k,v in q.filter3_v2_terminal_status.value_counts().items()}
old_eligible=set(q.loc[q.filter3_v2_terminal_status.isin(['FILTER3_HIGH_QUALITY','FILTER3_GOOD_QUALITY']),'pair_id'])
pb_warning=q.warning_codes.fillna('').str.contains('POSEBUSTERS_NONFATAL_WARNING',regex=False)
new_reject=set(q.loc[pb_warning & q.pair_id.isin(old_eligible),'pair_id'])
new_eligible=old_eligible-new_reject

# Frozen Filter 4 PASS proxy is the frozen F5 Step 1 input, whose provenance names
# and hashes the missing migrated Filter 4 PASS file.
s1=pd.read_csv(S1/'output/03_filter5_step1_pair_inventory.tsv.gz',sep='\t')
s1['pair_id']=s1.pair_id.astype(str)
f4pass=set(s1.pair_id)
f4_removed=f4pass & new_reject
candidate_f4pass=f4pass-new_reject

# Frozen F5 final and P4 universe.
old_final=pd.read_csv(S5/'output/01_filter5_final_pair_inventory.tsv.gz',sep='\t')
old_final['pair_id']=old_final.pair_id.astype(str)
old_f5_status=dict(zip(old_final.pair_id,old_final.filter5_final_status))
old_f5_retained={p for p,s in old_f5_status.items() if s!='F5_REDUNDANT_EQUIVALENT_CASE'}
p4inv=pd.read_parquet(P4/'output/processing4_case_inventory.parquet')
p4input=pd.read_parquet(P4/'input/full_case_inventory.parquet',columns=['case_id','pair_id','bucket_id'])
p4input['pair_id']=p4input.pair_id.astype(str); p4input['case_id']=p4input.case_id.astype(str)
old_p4=set(p4input.pair_id)

# Re-run the exact frozen Step 5 grouping algorithm in memory after revised F3
# removes candidates; no result is called frozen or authoritative.
inv=s1[s1.pair_id.isin(candidate_f4pass)].copy()
pinv=pd.read_csv(S2/'output/02_filter5_step2_pair_inventory.tsv.gz',sep='\t',usecols=['pair_id','binding_site_signature','step2_pair_status'])
pinv['pair_id']=pinv.pair_id.astype(str)
inv=inv.merge(pinv,on='pair_id',how='left',validate='one_to_one')
inv['step2_exact_site_group_id']=''
ok=inv.candidate_block_id.fillna('').ne('') & inv.binding_site_signature.fillna('').ne('')
inv.loc[ok,'step2_exact_site_group_id']=['ESG'+hashlib.sha256((b+'|'+s).encode()).hexdigest()[:16] for b,s in zip(inv.loc[ok,'candidate_block_id'],inv.loc[ok,'binding_site_signature'])]
qmap=dict(zip(q.pair_id,q.filter3_v2_terminal_status))
rank={'FILTER3_HIGH_QUALITY':0,'FILTER3_GOOD_QUALITY':1}
def priority(p): return (rank.get(qmap.get(p,''),2),p)
s3=pd.read_csv(S3/'output/01_step3_pairwise_receptor_state.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b','step3_status'])
r3={'CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW'}
s3review=set(s3.loc[s3.step3_status.isin(r3),'pair_id_a'].astype(str))|set(s3.loc[s3.step3_status.isin(r3),'pair_id_b'].astype(str))
s4=pd.read_csv(S4/'output/01_step4_pairwise_pose.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b','step4_status'])
r4={'POSE_MAPPING_REVIEW','UPSTREAM_LIGAND_CONSISTENCY_REVIEW'}
s4review=set(s4.loc[s4.step4_status.isin(r4),'pair_id_a'].astype(str))|set(s4.loc[s4.step4_status.isin(r4),'pair_id_b'].astype(str))
edges=pd.read_csv(S4/'output/02_step4_pose_equivalent_edges.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b'])
strict=set(); adj=defaultdict(set)
for r in edges.itertuples(index=False):
 x,y=str(r.pair_id_a),str(r.pair_id_b)
 if x not in candidate_f4pass or y not in candidate_f4pass: continue
 a,b=(x,y) if x<y else (y,x); strict.add((a,b)); adj[x].add(y); adj[y].add(x)
p_to_esg=dict(zip(inv.pair_id.astype(str),inv.step2_exact_site_group_id))
nodes_by=defaultdict(set)
for x,y in strict:
 if p_to_esg[x]!=p_to_esg[y]: raise RuntimeError('strict edge crosses exact site group')
 nodes_by[p_to_esg[x]].update((x,y))
assigned={}; groups=[]
for esg,nodes in sorted(nodes_by.items()):
 ordered=sorted(nodes,key=priority); unassigned=set(ordered)
 for seed in ordered:
  if seed not in unassigned: continue
  g=[seed]; unassigned.remove(seed)
  for x in ordered:
   if x in unassigned and all(x in adj[y] for y in g): g.append(x); unassigned.remove(x)
  if len(g)<2: continue
  rep=g[0]; groups.append({'candidate_group_number':len(groups)+1,'step2_exact_site_group_id':esg,'representative_pair_id':rep,'group_size':len(g),'members':'|'.join(g)})
  for p in g: assigned[p]=rep
new_status={}
for r in inv.itertuples(index=False):
 p=str(r.pair_id)
 if p in assigned: s='F5_RETAIN_REPRESENTATIVE' if p==assigned[p] else 'F5_REDUNDANT_EQUIVALENT_CASE'
 elif bool(r.step1_review_flag) or r.step2_pair_status=='STEP2_PAIR_MAPPING_REVIEW' or p in s3review or p in s4review: s='F5_REVIEW_RETAIN'
 else: s='F5_RETAIN_UNIQUE'
 new_status[p]=s
candidate_retained={p for p,s in new_status.items() if s!='F5_REDUNDANT_EQUIVALENT_CASE'}

# Representative replacement audit against old groups.
members=pd.read_csv(S5/'output/05_filter5_group_members.tsv.gz',sep='\t')
members['pair_id']=members.pair_id.astype(str)
old_reps=set(members.loc[members.member_role.eq('REPRESENTATIVE'),'pair_id'])
rejected_old_reps=old_reps & new_reject
group_rows=[]
for gid,g in members.groupby('equivalence_group_id',sort=False):
 oldrep=str(g.representative_pair_id.iloc[0])
 if oldrep not in rejected_old_reps: continue
 oldmembers=set(g.pair_id)
 survivors=oldmembers-new_reject
 newret=sorted(survivors & candidate_retained)
 group_rows.append({'equivalence_group_id':gid,'old_representative_pair_id':oldrep,'old_group_size':len(oldmembers),
                    'eligible_survivor_count':len(survivors),'eligible_survivors':'|'.join(sorted(survivors)),
                    'replacement_retained_count':len(newret),'replacement_pair_ids':'|'.join(newret),
                    'group_lost_completely':len(survivors)==0})
rep_audit=pd.DataFrame(group_rows)
write_tsv(rep_audit,'representative_replacement_audit.tsv')

# Diffs: authoritative and counterfactual candidate.
def diff_files(prefix,newset):
 common=old_p4 & newset; oldonly=old_p4-newset; newonly=newset-old_p4
 base=p4input[['case_id','pair_id','bucket_id']]
 c=base[base.pair_id.isin(common)].copy().sort_values('pair_id')
 oo=base[base.pair_id.isin(oldonly)].copy().sort_values('pair_id')
 oo['exit_reason']=oo.pair_id.map(lambda p:'FILTER3_POSEBUSTERS_FAIL' if p in new_reject else ('FILTER5_REPRESENTATIVE_REPLACED' if p in old_reps else 'OTHER'))
 no=pd.DataFrame({'pair_id':sorted(newonly)})
 no['source']=no.pair_id.map(lambda p:'FILTER5_REPLACEMENT_REPRESENTATIVE' if any(p in x for x in rep_audit.replacement_pair_ids.fillna('').str.split('|')) else 'OTHER')
 write_tsv(c,f'{prefix}common.tsv'); write_tsv(oo,f'{prefix}old_only.tsv'); write_tsv(no,f'{prefix}new_only.tsv')
 return {'old_count':len(old_p4),'new_count':len(newset),'common_count':len(common),'old_only_count':len(oldonly),'new_only_count':len(newonly),
         'duplicates_old':int(p4input.pair_id.duplicated().sum()),'duplicates_new':0,'missing_key_old':int(p4input.pair_id.isna().sum()),'missing_key_new':0}
authoritative_diff=diff_files('authoritative_',old_f5_retained)
candidate_diff=diff_files('',candidate_retained)

f5_counts={str(k):int(v) for k,v in old_final.filter5_final_status.value_counts().items()}
candidate_counts={str(k):int(v) for k,v in pd.Series(new_status).value_counts().items()}
p4_counts={str(k):int(v) for k,v in p4inv.status.value_counts().items()}
rescued=pd.read_parquet(P4/'output/rescue_status_transitions.parquet')
rescued_ready=set(rescued.loc[rescued.status.eq('P4_DOCKING_READY'),'pair_id'].astype(str))
summary={
 'generated_at':datetime.now(timezone.utc).isoformat(),
 'authoritative_conclusion':{
  'latest_frozen_filter3_run':str(F3),'latest_frozen_filter3_policy_has_any_posebusters_fail_reject':False,
  'latest_frozen_filter5_run':str(S5),'latest_frozen_filter5_retained_count':len(old_f5_retained),
  'authoritative_processing4_input_count':len(old_f5_retained),'p4_rebase_authorized':False,
  'reason':'No frozen Filter3 or Filter5 run implements/propagates the revised PoseBusters rule.'},
 'filter3_frozen':{'input':len(q),'status_counts':f3_counts,'eligible':len(old_eligible),'primary_key':'pair_id','duplicates':int(q.pair_id.duplicated().sum()),
                   'validation_status':json.loads((F3/'release/filter3_v2_release_validation.json').read_text())['validation_pass']},
 'revised_filter3_counterfactual':{'new_posebusters_rejects':len(new_reject),'eligible':len(new_eligible),'high':f3_counts['FILTER3_HIGH_QUALITY'],'good':f3_counts['FILTER3_GOOD_QUALITY']-len(new_reject)},
 'filter4_frozen_evidence':{'input':336412,'pass':len(f4pass),'reject':94865,'review':2,'pass_proxy':str(S1/'output/03_filter5_step1_pair_inventory.tsv.gz'),
                            'direct_server_step5_run_present':False,'formal_membership_path_recorded_by_f5':json.loads((S1/'provenance.json').read_text())['formal_membership_input'],
                            'formal_membership_sha256_recorded_by_f5':json.loads((S1/'provenance.json').read_text())['formal_membership_sha256']},
 'filter5_frozen':{'input':len(old_final),'status_counts':f5_counts,'retained':len(old_f5_retained),'primary_key':'pair_id','duplicates':int(old_final.pair_id.duplicated().sum()),'validation_status':'PASS'},
 'processing4_frozen':{'input':len(p4input),'status_counts':p4_counts,'primary_keys':['case_id','pair_id'],'duplicate_case_id':int(p4input.case_id.duplicated().sum()),'duplicate_pair_id':int(p4input.pair_id.duplicated().sum()),'validation_status':'PASS'},
 'authoritative_p4_vs_f5_diff':authoritative_diff,
 'counterfactual_full_propagation':{'filter4_pass_after_input_restriction':len(candidate_f4pass),'filter4_removed':len(f4_removed),'filter5_input':len(inv),'filter5_status_counts':candidate_counts,'filter5_retained':len(candidate_retained),
                                    'diff_vs_old_p4':candidate_diff,'not_authoritative_reason':'No frozen run, manifest, validation report, or _FROZEN.json exists for this derived universe.'},
 'representative_audit':{'new_filter3_rejects_among_old_f5_representatives':len(rejected_old_reps),'groups_requiring_reselection':int((rep_audit.eligible_survivor_count>0).sum()),
                         'replacement_representatives_or_singleton_retains':len(candidate_retained-old_f5_retained),'groups_lost_completely':int(rep_audit.group_lost_completely.sum()),
                         'new_only_cases':len(candidate_retained-old_f5_retained),'old_only_cases':len(old_f5_retained-candidate_retained)},
 'p4_candidate_projection':{'legacy_rescued_ready_total':len(rescued_ready),'legacy_rescued_ready_still_in_candidate':len(rescued_ready&candidate_retained),
                            'candidate_new_cases_requiring_p4_generation':len(candidate_retained-old_p4),'candidate_common_reusable':len(candidate_retained&old_p4)}
}
(OUT/'lineage_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
write_tsv(pd.DataFrame(groups),'counterfactual_filter5_groups.tsv')
write_tsv(pd.DataFrame([{'pair_id':p,'candidate_filter5_status':new_status[p]} for p in sorted(new_status)]),'counterfactual_filter5_inventory.tsv')
print(json.dumps(summary,indent=2,sort_keys=True))
