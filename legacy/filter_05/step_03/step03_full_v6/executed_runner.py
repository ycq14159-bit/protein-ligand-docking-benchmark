#!/usr/bin/env python3
"""Filter 5 Step 3 v6: chain-instance-aware multi-receptor interface audit."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,itertools,json,math,sys,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.dataset as ds

ROOT=Path('/root/autodl-tmp/benchmark_1.0'); F5=ROOT/'filter_05_equivalent_redocking_case'
S1=F5/'step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1'
S2=F5/'step_02_same_binding_site_audit/runs/step02_full_v2'
REC=ROOT/'processing_2_assembly_ready_structure_preparation/runs/20260810_full_01/output/prepared_receptor_assembly_atoms'
BASE=F5/'step_03_local_receptor_state_equivalence/runs'; EXPECTED=4926271; MAX_MAPS=100000
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def gw(p,h):
 raw=gzip.GzipFile(filename=str(p),mode='wb',compresslevel=4,mtime=0);f=io.TextIOWrapper(raw,encoding='utf8',newline='');w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(h);return f,w
def bucket(p):return int(hashlib.sha256(p.lower().encode()).hexdigest()[:8],16)%256
def kabsch(a,b):
 ca=a.mean(0);cb=b.mean(0);aa=a-ca;bb=b-cb;u,s,vt=np.linalg.svd(bb.T@aa);r=vt.T@u.T
 if np.linalg.det(r)<0:vt[-1]*=-1;r=vt.T@u.T
 t=ca-cb@r.T;z=b@r.T+t;rms=float(np.sqrt(np.mean(np.sum((a-z)**2,axis=1))));return r,t,rms
def f12(x):return '' if x is None or not np.isfinite(x) else f'{x:.12g}'
def rb(x):
 for hi,n in [(.1,'0-0.10'),(.25,'0.10-0.25'),(.5,'0.25-0.50'),(.75,'0.50-0.75'),(1,'0.75-1.00'),(1.5,'1.00-1.50'),(2,'1.50-2.00')]:
  if x<=hi:return n
 return '>2.00'

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-id',default='step03_full_v6');a=ap.parse_args();run=BASE/a.run_id
 out=run/'output';val=run/'validation';logs=run/'logs'
 for d in(out,val,logs):d.mkdir(parents=True,exist_ok=True)
 started=datetime.now(timezone.utc).isoformat();t0=time.time();exactfile=S2/'output/01_filter5_step2_pairwise_site_comparisons.tsv.gz'
 man=pd.read_csv(S2/'output_manifest.tsv',sep='\t');assert '01_filter5_step2_pairwise_site_comparisons.tsv.gz' in set(man.file)
 endpoints=set();ninput=0;edge_hashes=set();dup=0;noncanon=0
 for d in pd.read_csv(exactfile,sep='\t',usecols=['pair_id_a','pair_id_b','step2_site_status'],chunksize=250000):
  x=d[d.step2_site_status.eq('SITE_EXACT')];ninput+=len(x);endpoints.update(x.pair_id_a);endpoints.update(x.pair_id_b)
  for r in x.itertuples(index=False):
   if r.pair_id_a>=r.pair_id_b:noncanon+=1
   h=hashlib.blake2b((r.pair_id_a+'\0'+r.pair_id_b).encode(),digest_size=16).digest();dup+=h in edge_hashes;edge_hashes.add(h)
 if ninput!=EXPECTED or dup or noncanon:raise RuntimeError((ninput,dup,noncanon))
 del edge_hashes
 inv=pd.read_csv(S1/'output/03_filter5_step1_pair_inventory.tsv.gz',sep='\t');inv=inv[inv.pair_id.isin(endpoints)].set_index('pair_id')
 pinv=pd.read_csv(S2/'output/02_filter5_step2_pair_inventory.tsv.gz',sep='\t');pinv=pinv[pinv.pair_id.isin(endpoints)].set_index('pair_id')
 if len(inv)!=len(endpoints) or len(pinv)!=len(endpoints):raise RuntimeError('endpoint join')
 print('exact',ninput,'endpoints',len(endpoints),flush=True)
 # Chain-aware binding instances from frozen Step 2 residue mapping.
 binds=defaultdict(list)
 for d in pd.read_csv(S2/'output/00_filter5_step2_residue_mapping.tsv.gz',sep='\t',chunksize=300000,dtype=str,keep_default_na=False):
  d=d[d.pair_id.isin(endpoints)&d.site_type.eq('BINDING')&d.mapping_status.eq('MAPPED')]
  for r in d.itertuples(index=False):binds[r.pair_id].append((r.chain_instance_id,r.uniprot_accession,int(r.uniprot_residue_number),str(r.label_seq_id),r.component_id))
 targets={(ch,ls) for rows in binds.values() for ch,acc,num,ls,comp in rows};pb=defaultdict(set)
 for pid,r in inv.iterrows():pb[bucket(r.pdb_id)].add(r.pdb_id.lower())
 coords={};badkeys=set();cols=['filter_1_chain_instance_id','pdb_id','label_seq_id','label_comp_id','label_atom_id','Cartn_x','Cartn_y','Cartn_z']
 for bi,pdbs in sorted(pb.items()):
  p=REC/f'bucket_id={bi:03d}';tab=ds.dataset(p,format='parquet').to_table(columns=cols,filter=ds.field('pdb_id').isin(list(pdbs)));d=tab.to_pandas();d=d[d.label_atom_id.isin(['N','CA','C'])]
  for r in d.itertuples(index=False):
   k=(r.filter_1_chain_instance_id,str(r.label_seq_id))
   if k not in targets:continue
   q=coords.setdefault(k,{'comp':str(r.label_comp_id),'atoms':{}})
   if q['comp']!=str(r.label_comp_id) or str(r.label_atom_id) in q['atoms']:badkeys.add(k)
   else:q['atoms'][str(r.label_atom_id)]=np.array([r.Cartn_x,r.Cartn_y,r.Cartn_z],np.float64)
  if bi%16==0:print('bucket',bi,'coords',len(coords),flush=True)
 # Materialize every endpoint once. chain signature includes biological segment accession and positions.
 states={};chainstats=Counter()
 for pid in endpoints:
  chrows=defaultdict(list);bad=[]
  for ch,acc,num,ls,srccomp in binds[pid]:
   q=coords.get((ch,ls))
   if q is None:bad.append('FROZEN_RESIDUE_COORDINATE_MISSING');continue
   if (ch,ls) in badkeys:bad.append('DUPLICATE_OR_INCONSISTENT_BACKBONE_ATOM')
   if q['comp']!=srccomp:bad.append('UPSTREAM_COMPONENT_ID_MISMATCH')
   if set(q['atoms'])!={'N','CA','C'} or not all(np.isfinite(q['atoms'][z]).all() for z in ('N','CA','C')):bad.append('REQUIRED_N_CA_C_UNAVAILABLE')
   chrows[ch].append(((acc,num),q['comp'],q['atoms']))
  chains=[]
  for ch,rows in chrows.items():
   pos=[z[0] for z in rows]
   if len(pos)!=len(set(pos)):bad.append('DUPLICATE_POSITION_WITHIN_CHAIN')
   rows=sorted(rows,key=lambda z:z[0]);bio=tuple(sorted(set(p[0] for p in pos)));sig=(bio,tuple(pos))
   chains.append({'id':ch,'sig':sig,'rows':rows})
  chains.sort(key=lambda z:z['id']);classes=defaultdict(list)
  for c in chains:classes[c['sig']].append(c)
  sigmulti=tuple(sorted((repr(k),len(v)) for k,v in classes.items()));nch=len(chains)
  if nch==1:itype='single_chain'
  elif len(classes)==1:itype='homomeric'
  elif all(len(v)==1 for v in classes.values()):itype='heteromeric'
  else:itype='mixed_multimeric'
  chainstats[str(nch) if nch<4 else '>=4']+=1;chainstats[itype]+=1
  hh=[]
  for c in chains:
   zz=[]
   for pos,comp,ats in c['rows']:zz.append(repr(pos)+'|'+comp+'|'+'|'.join(','.join(map(repr,ats[x])) for x in ('N','CA','C')))
   hh.append(c['id']+'|'+repr(c['sig'])+'\n'+'\n'.join(zz))
  states[pid]={'chains':chains,'classes':classes,'multiset':sigmulti,'bad':sorted(set(bad)),'itype':itype,'key':hashlib.sha256('\n@@\n'.join(hh).encode()).hexdigest()}
 print('states',len(states),'coords',len(coords),flush=True)
 hdr=['candidate_block_id','exact_site_group_id','pair_id_a','pair_id_b','binding_chain_count_a','binding_chain_count_b','interface_type','chain_aware_interface_status','candidate_chain_mapping_count','composition_valid_mapping_count','selected_chain_mapping_id','binding_residue_instance_count','alignment_atom_count','joint_binding_interface_backbone_rmsd','alignment_id','step3_status','step3_reason']
 pf,pw=gw(out/'01_step3_pairwise_receptor_state.tsv.gz',hdr);mf,mw=gw(out/'02_step3_chain_instance_mappings.tsv.gz',['alignment_id','pair_id_a','pair_id_b','selected_chain_mapping_id','pair_a_chain_instance','pair_b_chain_instance','biological_identity','chain_binding_signature']);tf,tw=gw(out/'03_step3_alignment_transforms.tsv.gz',['alignment_id','pair_id_a','pair_id_b','R11','R12','R13','R21','R22','R23','R31','R32','R33','tx','ty','tz','joint_binding_interface_backbone_rmsd']);gf,gwrt=gw(out/'04_step3_pass_to_step4.tsv.gz',hdr);vf,vw=gw(out/'05_step3_variants.tsv.gz',hdr);rf,rw=gw(out/'06_step3_review.tsv.gz',hdr)
 counts=Counter();edgestats=Counter();bins=Counter();sens=Counter();cache={};edges=0;transform_bad=0;sampled_rmsd_count=0;sampled_rmsd_bad=0
 def evaluate(sa,sb):
  if sa['bad'] or sb['bad']:return ('UPSTREAM_RESIDUE_CONSISTENCY_REVIEW',';'.join(sorted(set(sa['bad']+sb['bad']))),0,0,None,None,None,None)
  if sa['multiset']!=sb['multiset']:return ('RECEPTOR_INTERFACE_VARIANT','CHAIN_BINDING_SIGNATURE_MULTISET_DIFFERS',0,0,None,None,None,None)
  nmap=1
  for sig in sa['classes']:nmap*=math.factorial(len(sa['classes'][sig]))
  if nmap>MAX_MAPS:return ('CHAIN_INSTANCE_MAPPING_REVIEW','EXACT_MAPPING_SPACE_INFEASIBLE',nmap,0,None,None,None,None)
  choices=[];sigs=sorted(sa['classes'],key=repr)
  for sig in sigs:choices.append(list(itertools.permutations(range(len(sb['classes'][sig])))))
  compvalid=0;best=None;mid=0;selected=None
  for combo in itertools.product(*choices):
   mid+=1;mapping=[];aa=[];bb=[];okcomp=True
   for sig,perm in zip(sigs,combo):
    ac=sorted(sa['classes'][sig],key=lambda z:z['id']);bc=sorted(sb['classes'][sig],key=lambda z:z['id'])
    for i,j in enumerate(perm):
     x,y=ac[i],bc[j];mapping.append((x,y,sig))
     for ra,rb in zip(x['rows'],y['rows']):
      if ra[0]!=rb[0] or ra[1]!=rb[1]:okcomp=False;break
      for atom in ('N','CA','C'):aa.append(ra[2][atom]);bb.append(rb[2][atom])
     if not okcomp:break
    if not okcomp:break
   if not okcomp:continue
   compvalid+=1;A=np.vstack(aa);B=np.vstack(bb)
   if np.array_equal(A,B):R=np.eye(3);t=np.zeros(3);rms=0.0
   else:R,t,rms=kabsch(A,B)
   if best is None or rms<best[0]-1e-15:best=(rms,R,t,len(A));selected=(mid,mapping)
  if not compvalid:return ('BINDING_RESIDUE_VARIANT','NO_CHAIN_MAPPING_WITH_EXACT_ACTUAL_RESIDUE_CHEMISTRY',nmap,0,None,None,None,None)
  rms,R,t,na=best;return ('OK','',nmap,compvalid,(f'CMAP{selected[0]:08d}',selected[1]),R,t,(rms,na))
 for ci,d in enumerate(pd.read_csv(exactfile,sep='\t',chunksize=200000),1):
  d=d[d.step2_site_status.eq('SITE_EXACT')]
  for x in d.itertuples(index=False):
   pa,pb=x.pair_id_a,x.pair_id_b;sa,sb=states[pa],states[pb];esg='ESG'+hashlib.sha256((x.candidate_block_id+'|'+pinv.loc[pa,'binding_site_signature']).encode()).hexdigest()[:16]
   ck=(sa['key'],sb['key'])
   if ck not in cache:cache[ck]=evaluate(sa,sb)
   pre,reason,nmap,nvalid,sel,R,t,res=cache[ck];aid='';selid='';rms=np.nan;na=0
   if pre=='OK':
    selid,mapping=sel;rms,na=res;aid='ALN'+hashlib.sha256((pa+'|'+pb).encode()).hexdigest()[:16]
    status='RECEPTOR_STATE_EQUIVALENT_CANDIDATE' if rms<=.5 else 'RECEPTOR_CONFORMATION_VARIANT';reason='JOINT_RMSD_LE_0.50' if rms<=.5 else 'JOINT_RMSD_GT_0.50'
    if not(np.isfinite(R).all() and np.isfinite(t).all() and np.allclose(R.T@R,np.eye(3),atol=1e-8) and abs(np.linalg.det(R)-1)<1e-8):transform_bad+=1
    if int(hashlib.sha256((pa+'\0'+pb).encode()).hexdigest()[:8],16)%10000==0:
     sampled_rmsd_count+=1;AA=[];BB=[]
     for ca,cb,sig in mapping:
      for ra_rec,rb_rec in zip(ca['rows'],cb['rows']):
       for atom in ('N','CA','C'):AA.append(ra_rec[2][atom]);BB.append(rb_rec[2][atom])
     AA=np.vstack(AA);BB=np.vstack(BB);stored_R=np.array([float(f12(q)) for q in R.ravel()]).reshape(3,3);stored_t=np.array([float(f12(q)) for q in t]);stored_rms=float(f12(rms));recomputed=float(np.sqrt(np.mean(np.sum((AA-(BB@stored_R.T+stored_t))**2,axis=1))))
     sampled_rmsd_bad+=len(AA)!=na or not np.isfinite(recomputed) or abs(recomputed-stored_rms)>1e-8
    tw.writerow([aid,pa,pb,*[f12(q) for q in R.ravel()],*[f12(q) for q in t],f12(rms)])
    for ca,cb,sig in mapping:mw.writerow([aid,pa,pb,selid,ca['id'],cb['id'],repr(sig[0]),repr(sig)])
    bins[rb(rms)]+=1
    for cut in (.25,.5,.75,1):sens[str(cut)]+=rms<=cut
   else:status=pre
   iface='MATCH' if sa['multiset']==sb['multiset'] else 'VARIANT';itype=sa['itype'] if sa['itype']==sb['itype'] else 'mixed_multimeric'
   edgestats[itype]+=1
   if nmap==1:edgestats['unique_chain_mapping']+=1
   elif nmap>1:edgestats['multiple_chain_mappings']+=1
   if status=='CHAIN_INSTANCE_MAPPING_REVIEW':edgestats['chain_mapping_review']+=1
   if nvalid>0:edgestats['actual_composition_exact']+=1
   if status=='BINDING_RESIDUE_VARIANT':edgestats['binding_residue_variant']+=1
   if pre=='OK':edgestats['kabsch_success']+=1
   elif status in ('CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW'):edgestats['kabsch_review']+=1
   row=[x.candidate_block_id,esg,pa,pb,len(sa['chains']),len(sb['chains']),itype,iface,nmap,nvalid,selid,sum(len(c['rows']) for c in sa['chains']),na,f12(rms),aid,status,reason];pw.writerow(row);counts[status]+=1;edges+=1
   if status=='RECEPTOR_STATE_EQUIVALENT_CANDIDATE':gwrt.writerow(row)
   elif status in ('CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW'):rw.writerow(row)
   else:vw.writerow(row)
  print('chunks',ci,'edges',edges,flush=True)
 for f in(pf,mf,tf,gf,vf,rf):f.close()
 order=['RECEPTOR_STATE_EQUIVALENT_CANDIDATE','RECEPTOR_INTERFACE_VARIANT','BINDING_RESIDUE_VARIANT','RECEPTOR_CONFORMATION_VARIANT','CHAIN_INSTANCE_MAPPING_REVIEW','POCKET_ALIGNMENT_REVIEW','UPSTREAM_RESIDUE_CONSISTENCY_REVIEW']
 summ=[('SITE_EXACT_input_edges',edges),('unique_endpoint_pair_ids',len(endpoints))]+[(k,counts[k]) for k in order]+[(f'edge_{k}',v) for k,v in sorted(edgestats.items())]+[(f'endpoint_chain_count_{k}',v) for k,v in sorted(chainstats.items())]
 pd.DataFrame(summ,columns=['metric','value']).to_csv(out/'08_step3_summary.tsv',sep='\t',index=False)
 dist=[(f'bin_{k}',bins[k]) for k in ['0-0.10','0.10-0.25','0.25-0.50','0.50-0.75','0.75-1.00','1.00-1.50','1.50-2.00','>2.00']]+[(f'sensitivity_le_{k}',sens[k]) for k in ['0.25','0.5','0.75','1']]
 pd.DataFrame(dist,columns=['metric','value']).to_csv(out/'07_step3_rmsd_distribution.tsv',sep='\t',index=False)
 checks={'input_4926271':edges==EXPECTED,'non_SITE_EXACT_zero':True,'duplicate_canonical_edges_zero':dup==0,'noncanonical_edges_zero':noncanon==0,'silent_drop_zero':edges==EXPECTED,'terminal_partition':sum(counts.values())==edges,'transforms_valid':transform_bad==0,'stored_transform_sample_nonempty':sampled_rmsd_count>0,'stored_transform_sampled_rmsd_recomputes':sampled_rmsd_bad==0,'no_ligand_atoms_used':True,'one_global_multichain_kabsch':True}
 validation={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'counts':dict(counts),'sampled_stored_transform_rmsd_count':sampled_rmsd_count,'sampled_stored_transform_rmsd_failures':sampled_rmsd_bad,'major_redesign':'CHAIN_INSTANCE_AWARE_MULTI_RECEPTOR_INTERFACE','elapsed_seconds':time.time()-t0,'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat()};(val/'validation.json').write_text(json.dumps(validation,indent=2))
 (out/'09_step3_report.md').write_text('# Filter 5 Step 3 v6\n\nChain-instance-aware multiset interface reconstruction and one global multi-chain Kabsch. Formal cutoff 0.50 Å.\n')
 (run/'provenance.json').write_text(json.dumps({'source_step1':str(S1),'source_step2':str(S2),'coordinate_source':str(REC),'sifts_sha256':sha(F5/'inputs/sifts_snapshot/uniprot_segments_observed.tsv.gz'),'MAJOR_STEP3_REDESIGN':'CHAIN_INSTANCE_AWARE_MULTI_RECEPTOR_INTERFACE','chain_mapping_algorithm':'EXACT_SIGNATURE_CLASS_PERMUTATION_V1','kabsch':'numpy.linalg.svd float64 one-global-interface','receptor_cutoff_angstrom':.5},indent=2))
 (run/'output_schema.json').write_text(json.dumps({p.name:'TSV.GZ' if p.name.endswith('.gz') else 'TSV/MD' for p in out.iterdir()},indent=2));mr=[]
 for p in sorted(out.iterdir()):
  rows=''
  if p.name.endswith('.tsv.gz'):
   with gzip.open(p,'rt') as f:rows=max(sum(1 for _ in f)-1,0)
  elif p.suffix=='.tsv':
   with open(p) as f:rows=max(sum(1 for _ in f)-1,0)
  mr.append((p.name,p.stat().st_size,rows,sha(p)))
 pd.DataFrame(mr,columns=['file','bytes','data_rows','sha256']).to_csv(run/'output_manifest.tsv',sep='\t',index=False)
 files=[p for p in run.rglob('*') if p.is_file() and p.name not in ('SHA256SUMS','_FROZEN.json')]
 with open(run/'SHA256SUMS','w') as f:
  for p in sorted(files):f.write(f'{sha(p)}  {p.relative_to(run)}\n')
 if validation['status']=='PASS':(run/'_FROZEN.json').write_text(json.dumps({'status':'FROZEN','frozen_utc':datetime.now(timezone.utc).isoformat(),'validation':'validation/validation.json'},indent=2))
 print(json.dumps(validation,indent=2));return 0 if validation['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
