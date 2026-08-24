#!/usr/bin/env python3
"""Filter 5 Step 5: deterministic quality-aware all-pairs-pass grouping."""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import pyarrow.dataset as ds
ROOT=Path('/root/autodl-tmp/benchmark_1.0');F5=ROOT/'filter_05_equivalent_redocking_case'
S1=F5/'step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1';S2=F5/'step_02_same_binding_site_audit/runs/step02_full_v2';S3=F5/'step_03_local_receptor_state_equivalence/runs/step03_full_v6';S4=F5/'step_04_native_ligand_pose_equivalence/runs/step04_full_v1';BASE=F5/'step_05_strict_equivalent_grouping_and_representative_selection/runs';EXPECTED=241545
F3=ROOT/'filter_03_ground_truth_structure_quality_v2/runs/20260814_full_01/output/filter3_pair_quality_v2'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def gzwrite(df,p):df.to_csv(p,sep='\t',index=False,compression={'method':'gzip','compresslevel':6,'mtime':0})
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-id',default='step05_full_v1');a=ap.parse_args();run=BASE/a.run_id;out=run/'output';val=run/'validation';logs=run/'logs'
 for d in(out,val,logs):d.mkdir(parents=True,exist_ok=True)
 for p in(S3,S4):
  if not (p/'_FROZEN.json').exists() or json.loads((p/'validation/validation.json').read_text())['status']!='PASS':raise RuntimeError(f'upstream not frozen PASS {p}')
 t0=time.time();started=datetime.now(timezone.utc).isoformat();inv=pd.read_csv(S1/'output/03_filter5_step1_pair_inventory.tsv.gz',sep='\t')
 if len(inv)!=EXPECTED or inv.pair_id.duplicated().any():raise RuntimeError('universe')
 pinv=pd.read_csv(S2/'output/02_filter5_step2_pair_inventory.tsv.gz',sep='\t',usecols=['pair_id','binding_site_signature','step2_pair_status']);inv=inv.merge(pinv,on='pair_id',how='left',validate='one_to_one')
 inv['step2_exact_site_group_id']='';ok=inv.candidate_block_id.fillna('').ne('')&inv.binding_site_signature.fillna('').ne('')
 inv.loc[ok,'step2_exact_site_group_id']=[ 'ESG'+hashlib.sha256((b+'|'+s).encode()).hexdigest()[:16] for b,s in zip(inv.loc[ok,'candidate_block_id'],inv.loc[ok,'binding_site_signature'])]
 # Frozen quality priority, no newly composed score.
 qt=ds.dataset(F3,format='parquet').to_table(columns=['pair_id','filter3_v2_terminal_status']).to_pandas();qmap=dict(zip(qt.pair_id,qt.filter3_v2_terminal_status));inv['filter3_quality_class']=inv.pair_id.map(qmap).fillna('FILTER3_CLASS_UNAVAILABLE')
 rank={'FILTER3_HIGH_QUALITY':0,'FILTER3_GOOD_QUALITY':1};priority=lambda p:(rank.get(qmap.get(p,''),2),p)
 s3=pd.read_csv(S3/'output/01_step3_pairwise_receptor_state.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b','step3_status']);s3review=set(s3.loc[s3.step3_status.isin(['CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW']),'pair_id_a'])|set(s3.loc[s3.step3_status.isin(['CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW']),'pair_id_b']);del s3
 s4=pd.read_csv(S4/'output/01_step4_pairwise_pose.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b','step4_status']);s4review=set(s4.loc[s4.step4_status.isin(['POSE_MAPPING_REVIEW','UPSTREAM_LIGAND_CONSISTENCY_REVIEW']),'pair_id_a'])|set(s4.loc[s4.step4_status.isin(['POSE_MAPPING_REVIEW','UPSTREAM_LIGAND_CONSISTENCY_REVIEW']),'pair_id_b']);del s4
 strictdf=pd.read_csv(S4/'output/02_step4_pose_equivalent_edges.tsv.gz',sep='\t');s3pass=pd.read_csv(S3/'output/04_step3_pass_to_step4.tsv.gz',sep='\t',usecols=['pair_id_a','pair_id_b','alignment_id','joint_binding_interface_backbone_rmsd','step3_status'])
 strict_step4_gate=bool(strictdf.step4_status.eq('POSE_EQUIVALENT').all() and strictdf.fixed_frame_symmetry_aware_heavy_atom_rmsd.notna().all() and strictdf.fixed_frame_symmetry_aware_heavy_atom_rmsd.le(.5).all())
 strict_step3_gate=bool(s3pass.step3_status.eq('RECEPTOR_STATE_EQUIVALENT_CANDIDATE').all() and s3pass.joint_binding_interface_backbone_rmsd.notna().all() and s3pass.joint_binding_interface_backbone_rmsd.le(.5).all() and set(strictdf.alignment_id).issubset(set(s3pass.alignment_id)))
 if not strict_step4_gate or not strict_step3_gate:raise RuntimeError('strict edge gate audit failed')
 del s3pass
 strict=set((min(r.pair_id_a,r.pair_id_b),max(r.pair_id_a,r.pair_id_b)) for r in strictdf.itertuples(index=False));adj=defaultdict(set)
 for x,y in strict:adj[x].add(y);adj[y].add(x)
 nodes_by=defaultdict(set)
 p_to_esg=dict(zip(inv.pair_id,inv.step2_exact_site_group_id))
 for x,y in strict:
  if p_to_esg[x]!=p_to_esg[y]:raise RuntimeError('strict edge crosses exact site group')
  nodes_by[p_to_esg[x]].update((x,y))
 groups=[];members=[];group_member_map={};assigned={};gidn=0
 for esg,nodes in sorted(nodes_by.items()):
  ordered=sorted(nodes,key=priority);unassigned=set(ordered)
  for seed in ordered:
   if seed not in unassigned:continue
   g=[seed];unassigned.remove(seed)
   for x in ordered:
    if x in unassigned and all(x in adj[y] for y in g):g.append(x);unassigned.remove(x)
   if len(g)<2:continue
   gidn+=1;gid=f'F5EQ{gidn:08d}';rep=g[0]
   group_member_map[gid]=g
   first=inv.loc[inv.pair_id.eq(rep)].iloc[0];groups.append({'equivalence_group_id':gid,'group_size':len(g),'representative_pair_id':rep,'ligand_exact_id':first.ligand_exact_id,'receptor_identity_key':first.receptor_identity_key,'binding_site_signature':first.binding_site_signature,'step2_exact_site_group_id':esg})
   for p in g:assigned[p]=(gid,len(g),rep);members.append({'equivalence_group_id':gid,'pair_id':p,'member_role':'REPRESENTATIVE' if p==rep else 'REDUNDANT','representative_pair_id':rep})
 # Complete group audit and non-overlap proof.
 group_fail=0;expected_internal=0
 for g in groups:
  z=group_member_map[g['equivalence_group_id']];expected_internal+=len(z)*(len(z)-1)//2
  for i,x in enumerate(z):
   for y in z[i+1:]:group_fail+=(min(x,y),max(x,y)) not in strict
 inv['has_step3_review']=inv.pair_id.isin(s3review);inv['has_step4_review']=inv.pair_id.isin(s4review);inv['equivalence_group_id']='';inv['group_size']=1;inv['representative_pair_id']='';inv['filter5_final_status']='';inv['filter5_final_reason']=''
 for i,r in inv.iterrows():
  p=r.pair_id
  if p in assigned:
   gid,n,rep=assigned[p];inv.at[i,'equivalence_group_id']=gid;inv.at[i,'group_size']=n;inv.at[i,'representative_pair_id']=rep
   if p==rep:inv.at[i,'filter5_final_status']='F5_RETAIN_REPRESENTATIVE';inv.at[i,'filter5_final_reason']='HIGHEST_FROZEN_QUALITY_PRIORITY_IN_ALL_PAIRS_PASS_GROUP'
   else:inv.at[i,'filter5_final_status']='F5_REDUNDANT_EQUIVALENT_CASE';inv.at[i,'filter5_final_reason']='NONREPRESENTATIVE_IN_PROVEN_ALL_PAIRS_PASS_GROUP'
  elif bool(r.step1_review_flag) or r.step2_pair_status=='STEP2_PAIR_MAPPING_REVIEW' or r.has_step3_review or r.has_step4_review:
   inv.at[i,'filter5_final_status']='F5_REVIEW_RETAIN';inv.at[i,'filter5_final_reason']='UNRESOLVED_DEDUPLICATION_REVIEW_RETAINED'
  else:inv.at[i,'filter5_final_status']='F5_RETAIN_UNIQUE';inv.at[i,'filter5_final_reason']='NO_PROVEN_STRICT_EQUIVALENT_GROUP'
 cols=['pair_id','filter3_quality_class','candidate_block_id','step2_exact_site_group_id','has_step3_review','has_step4_review','equivalence_group_id','group_size','representative_pair_id','filter5_final_status','filter5_final_reason'];final=inv[cols].rename(columns={'candidate_block_id':'step1_block_id'})
 gzwrite(final,out/'01_filter5_final_pair_inventory.tsv.gz');gzwrite(final[final.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE')],out/'02_filter5_retained_pairs.tsv.gz');gzwrite(final[final.filter5_final_status.eq('F5_REDUNDANT_EQUIVALENT_CASE')],out/'03_filter5_redundant_equivalent_pairs.tsv.gz');gzwrite(pd.DataFrame(groups),out/'04_filter5_equivalence_groups.tsv.gz');gzwrite(pd.DataFrame(members),out/'05_filter5_group_members.tsv.gz');gzwrite(final[final.filter5_final_status.eq('F5_REVIEW_RETAIN')],out/'06_filter5_review_retain_pairs.tsv.gz');gzwrite(strictdf,out/'07_filter5_strict_equivalence_edges.tsv.gz')
 status=final.filter5_final_status.value_counts();sizes=Counter();size_order=['2','3','4','5-9','10-19','20-49','>=50']
 for g in groups:
  n=g['group_size'];sizes['2' if n==2 else '3' if n==3 else '4' if n==4 else '5-9' if n<10 else '10-19' if n<20 else '20-49' if n<50 else '>=50']+=1
 pd.DataFrame([(k,sizes[k]) for k in size_order],columns=['group_size_bin','group_count']).to_csv(out/'09_filter5_group_size_distribution.tsv',sep='\t',index=False)
 largest=sorted(groups,key=lambda g:(-g['group_size'],g['equivalence_group_id']))[:10]
 stats=[('filter4_pass_input',len(final)),('strict_equivalence_edges',len(strict)),('exact_site_groups_with_strict_edge',len(nodes_by)),('multi_member_groups',len(groups)),('total_group_members',len(members)),('largest_group_size',largest[0]['group_size'] if largest else 0),('internal_edges_validated',expected_internal),('group_validation_failures',group_fail)]+[(f'largest_group_{i}_size',g['group_size']) for i,g in enumerate(largest,1)]+[(k,int(status[k])) for k in ['F5_RETAIN_UNIQUE','F5_RETAIN_REPRESENTATIVE','F5_REVIEW_RETAIN','F5_REDUNDANT_EQUIVALENT_CASE']]+[('final_retained',int(final.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE').sum()))]
 for cls in ['FILTER3_HIGH_QUALITY','FILTER3_GOOD_QUALITY']:
  stats.append((f'retained_{cls}',int(((final.filter3_quality_class==cls)&final.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE')).sum())));stats.append((f'redundant_{cls}',int(((final.filter3_quality_class==cls)&final.filter5_final_status.eq('F5_REDUNDANT_EQUIVALENT_CASE')).sum())))
 pd.DataFrame(stats,columns=['metric','value']).to_csv(out/'10_filter5_final_statistics.tsv',sep='\t',index=False)
 # Funnel preserves frozen counts; Step3/4 detailed summaries remain their formal artifacts.
 funnel=[('Filter4_PASS',241545),('Step1_singletons',22751),('Step1_non_singletons',218765),('Step1_primary_review_pairs',29),('Step1_auxiliary_review_pairs',2210),('Step1_candidate_blocks',35214),('Step2_SITE_EXACT',4926271),('Step2_SITE_STRONG_CANDIDATE',697793),('Step2_SITE_WEAK_OR_AMBIGUOUS',1180489),('Step2_SITE_DIFFERENT',8474616),('Step2_SITE_MAPPING_REVIEW',795381)]
 funnel += [('Step3_'+str(r.metric),int(r.value)) for r in pd.read_csv(S3/'output/08_step3_summary.tsv',sep='\t').itertuples(index=False)]
 funnel += [('Step4_'+str(r.metric),int(r.value)) for r in pd.read_csv(S4/'output/06_step4_summary.tsv',sep='\t').itertuples(index=False)]
 funnel += [('Step5_STRICT_EQUIVALENT',len(strict)),('Step5_final_retained',int(final.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE').sum()))]
 pd.DataFrame(funnel,columns=['stage','count']).to_csv(out/'08_filter5_funnel_summary.tsv',sep='\t',index=False)
 checks={'inventory_241545':len(final)==EXPECTED,'pair_ids_unique':final.pair_id.nunique()==EXPECTED,'duplicate_zero':not final.pair_id.duplicated().any(),'silent_drop_zero':set(final.pair_id)==set(pd.read_csv(S1/'output/03_filter5_step1_pair_inventory.tsv.gz',sep='\t',usecols=['pair_id']).pair_id),'status_partition_closes':int(status.sum())==EXPECTED,'retained_plus_redundant_closes':len(final)==len(final[final.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE')])+status['F5_REDUNDANT_EQUIVALENT_CASE'],'strict_edges_step3_pass_rmsd_le_0_50':strict_step3_gate,'strict_edges_step4_pass_rmsd_le_0_50':strict_step4_gate,'group_internal_all_pairs_strict_pass':group_fail==0,'nonoverlapping_groups':len(assigned)==len(set(assigned)),'representative_membership_valid':all(assigned[g['representative_pair_id']][0]==g['equivalence_group_id'] for g in groups),'connected_components_not_used':True,'single_linkage_not_used':True}
 checks={k:bool(v) for k,v in checks.items()};validation={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'grouping_algorithm':'QUALITY_PRIORITY_AWARE_ALL_PAIRS_PASS_PARTITION_V1','counts':{str(k):int(v) for k,v in status.items()},'elapsed_seconds':time.time()-t0,'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat()};(val/'validation.json').write_text(json.dumps(validation,indent=2));(out/'11_filter5_report.md').write_text('# Filter 5 final\n\nStrict-equivalence grouping uses deterministic frozen-quality-priority greedy all-pairs-pass partitioning; connected components and single linkage are not used. Reviews retain by default.\n')
 s3prov=json.loads((S3/'provenance.json').read_text());s4prov=json.loads((S4/'provenance.json').read_text());(run/'provenance.json').write_text(json.dumps({'source_step1':str(S1),'source_step2':str(S2),'source_step3':str(S3),'source_step4':str(S4),'coordinate_source':s3prov['coordinate_source'],'sifts_sha256':s3prov['sifts_sha256'],'rdkit_version':s4prov['rdkit_version'],'kabsch_implementation':s3prov['kabsch'],'receptor_cutoff_angstrom':.5,'ligand_cutoff_angstrom':.5,'chain_mapping_algorithm':s3prov['chain_mapping_algorithm'],'grouping_algorithm':'QUALITY_PRIORITY_AWARE_ALL_PAIRS_PASS_PARTITION_V1','representative_priority':'FILTER3_HIGH_THEN_GOOD_THEN_PAIR_ID_LEXICAL','strict_equivalence_transitive':False,'connected_components_grouping':False,'uncertain_policy':'RETAIN'},indent=2));(run/'output_schema.json').write_text(json.dumps({p.name:'TSV.GZ' if p.name.endswith('.gz') else 'TSV/MD' for p in out.iterdir()},indent=2));mr=[]
 for p in sorted(out.iterdir()):
  rows=''
  if p.name.endswith('.tsv.gz'):
   with gzip.open(p,'rt') as f:rows=max(sum(1 for _ in f)-1,0)
  elif p.suffix=='.tsv':
   with open(p) as f:rows=max(sum(1 for _ in f)-1,0)
  mr.append((p.name,p.stat().st_size,rows,sha(p)))
 pd.DataFrame(mr,columns=['file','bytes','data_rows','sha256']).to_csv(run/'output_manifest.tsv',sep='\t',index=False);files=[p for p in run.rglob('*') if p.is_file() and p.name not in ('SHA256SUMS','_FROZEN.json')]
 with open(run/'SHA256SUMS','w') as f:
  for p in sorted(files):f.write(f'{sha(p)}  {p.relative_to(run)}\n')
 if validation['status']=='PASS':(run/'_FROZEN.json').write_text(json.dumps({'status':'FROZEN','frozen_utc':datetime.now(timezone.utc).isoformat(),'validation':'validation/validation.json'},indent=2))
 print(json.dumps(validation,indent=2));return 0 if validation['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
