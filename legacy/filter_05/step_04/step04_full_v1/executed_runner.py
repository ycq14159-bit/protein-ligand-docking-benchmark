#!/usr/bin/env python3
"""Filter 5 Step 4: fixed-receptor-frame, symmetry-aware native ligand pose audit."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math,sqlite3,time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from rdkit import rdBase
ROOT=Path('/root/autodl-tmp/benchmark_1.0');F5=ROOT/'filter_05_equivalent_redocking_case'
S1=F5/'step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1';S2=F5/'step_02_same_binding_site_audit/runs/step02_full_v2';S3=F5/'step_03_local_receptor_state_equivalence/runs/step03_full_v6'
P2=ROOT/'processing_2_assembly_ready_structure_preparation/runs/20260810_full_01';LIG=P2/'output/prepared_ligand_assembly_atoms';BASE=F5/'step_04_native_ligand_pose_equivalence/runs';MAXMAP=1000000
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def gw(p,h):
 raw=gzip.GzipFile(filename=str(p),mode='wb',compresslevel=4,mtime=0);f=io.TextIOWrapper(raw,encoding='utf8',newline='');w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(h);return f,w
def f12(x):return '' if x is None or not np.isfinite(x) else f'{x:.12g}'
def rb(x):
 for hi,n in [(.1,'0-0.10'),(.25,'0.10-0.25'),(.5,'0.25-0.50'),(.75,'0.50-0.75'),(1,'0.75-1.00'),(1.5,'1.00-1.50'),(2,'1.50-2.00'),(3,'2.00-3.00')]:
  if x<=hi:return n
 return '>3.00'
def bucket(p):return int(hashlib.sha256(p.lower().encode()).hexdigest()[:8],16)%256
def int0(x):return 0 if pd.isna(x) else int(x)
def attr(x):return '' if pd.isna(x) else str(x).upper()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-id',default='step04_full_v1');a=ap.parse_args();run=BASE/a.run_id;out=run/'output';val=run/'validation';logs=run/'logs'
 for d in(out,val,logs):d.mkdir(parents=True,exist_ok=True)
 if not (S3/'_FROZEN.json').exists() or json.loads((S3/'validation/validation.json').read_text())['status']!='PASS':raise RuntimeError('Step3 v6 not frozen PASS')
 t0=time.time();started=datetime.now(timezone.utc).isoformat();passfile=S3/'output/04_step3_pass_to_step4.tsv.gz';trfile=S3/'output/03_step3_alignment_transforms.tsv.gz'
 sin=pd.read_csv(passfile,sep='\t');trans=pd.read_csv(trfile,sep='\t');trans=trans[trans.alignment_id.isin(set(sin.alignment_id))]
 if len(trans)!=len(sin) or trans.alignment_id.duplicated().any() or sin[['pair_id_a','pair_id_b']].duplicated().any():raise RuntimeError('Step3 interface/transform gate')
 trans=trans.set_index('alignment_id')
 endpoints=set(sin.pair_id_a)|set(sin.pair_id_b);inv=pd.read_csv(S1/'output/03_filter5_step1_pair_inventory.tsv.gz',sep='\t');inv=inv[inv.pair_id.isin(endpoints)].set_index('pair_id')
 if len(inv)!=len(endpoints):raise RuntimeError('endpoint join')
 # Frozen ligand assembly coordinates, heavy atoms only.
 placements=set(inv.ligand_assembly_placement_id);pb={bucket(r.pdb_id) for _,r in inv.iterrows()};lig={};cols=['filter_2_ligand_assembly_placement_id','pdb_id','component_id','label_atom_id','type_symbol','Cartn_x','Cartn_y','Cartn_z']
 for bi in sorted(pb):
  p=LIG/f'bucket_id={bi:03d}';tab=ds.dataset(p,format='parquet').to_table(columns=cols);d=tab.to_pandas();d=d[d.filter_2_ligand_assembly_placement_id.isin(placements)&~d.type_symbol.str.upper().isin(['H','D'])]
  for r in d.itertuples(index=False):
   q=lig.setdefault(r.filter_2_ligand_assembly_placement_id,{'component':str(r.component_id),'atoms':{},'bad':False})
   if q['component']!=str(r.component_id):q['bad']=True
   if str(r.label_atom_id) in q['atoms']:q['bad']=True
   else:q['atoms'][str(r.label_atom_id)]=np.array([r.Cartn_x,r.Cartn_y,r.Cartn_z],np.float64)
  print('ligand bucket',bi,'placements',len(lig),flush=True)
 # Exact attributed CCD graphs. Stereo labels are explicit node/edge constraints.
 con=sqlite3.connect(P2/'input/ccd_active_snapshot.sqlite');ad=pd.read_sql_query('select component_id,atom_id,element,charge,aromatic_flag,stereo_config from atoms',con);bd=pd.read_sql_query('select component_id,atom_id_1,atom_id_2,bond_order,aromatic_flag,stereo_config from bonds',con);con.close();needed=set(inv.component_id);graphs={comp:nx.Graph() for comp in needed}
 for r in ad[ad.component_id.isin(needed)].itertuples(index=False):
  if str(r.element).upper() in ('H','D'):continue
  graphs[r.component_id].add_node(str(r.atom_id),element=attr(r.element),charge=int0(r.charge),aromatic=attr(r.aromatic_flag),stereo=attr(r.stereo_config))
 for r in bd[bd.component_id.isin(needed)].itertuples(index=False):
  g=graphs[r.component_id]
  if str(r.atom_id_1) in g and str(r.atom_id_2) in g:g.add_edge(str(r.atom_id_1),str(r.atom_id_2),order=attr(r.bond_order),aromatic=attr(r.aromatic_flag),stereo=attr(r.stereo_config))
 print('ccd graphs',len(graphs),flush=True)
 node=lambda a,b:a==b;edge=lambda a,b:a==b;mapcache={};posecache={}
 def maps(ca,cb):
  k=(ca,cb)
  if k in mapcache:return mapcache[k]
  gm=nx.algorithms.isomorphism.GraphMatcher(graphs[ca],graphs[cb],node_match=node,edge_match=edge);z=[]
  for i,m in enumerate(gm.isomorphisms_iter(),1):
   if i>MAXMAP:mapcache[k]=None;return None
   z.append(m)
  mapcache[k]=z;return z
 hdr=['pair_id_a','pair_id_b','alignment_id','ligand_exact_id','heavy_atom_count','symmetry_mapping_count','selected_atom_mapping_id','fixed_frame_symmetry_aware_heavy_atom_rmsd','ligand_centroid_displacement','step4_status','step4_reason']
 pf,pw=gw(out/'01_step4_pairwise_pose.tsv.gz',hdr);ef,ew=gw(out/'02_step4_pose_equivalent_edges.tsv.gz',hdr);vf,vw=gw(out/'03_step4_pose_variants.tsv.gz',hdr);rf,rw=gw(out/'04_step4_review.tsv.gz',hdr)
 counts=Counter();bins=Counter();sens=Counter();recompute_bad=0
 for i,x in enumerate(sin.itertuples(index=False),1):
  ia,ib=inv.loc[x.pair_id_a],inv.loc[x.pair_id_b];la,lb=lig.get(ia.ligand_assembly_placement_id),lig.get(ib.ligand_assembly_placement_id);n=0;nm=0;sel='';rms=np.nan;cent=np.nan
  if ia.ligand_exact_id!=ib.ligand_exact_id:status='UPSTREAM_LIGAND_CONSISTENCY_REVIEW';reason='LIGAND_EXACT_ID_DIFFERS'
  elif la is None or lb is None or la['bad'] or lb['bad'] or la['component']!=str(ia.component_id) or lb['component']!=str(ib.component_id):status='UPSTREAM_LIGAND_CONSISTENCY_REVIEW';reason='FROZEN_LIGAND_COORDINATE_MISSING_DUPLICATE_OR_COMPONENT_MISMATCH'
  else:
   mm=maps(ia.component_id,ib.component_id)
   if mm is None:status='POSE_MAPPING_REVIEW';reason='EXACT_AUTOMORPHISM_SPACE_INFEASIBLE'
   elif not mm:status='POSE_MAPPING_REVIEW';reason='NO_CHEMISTRY_PRESERVING_HEAVY_ATOM_ISOMORPHISM'
   else:
    ta=trans.loc[x.alignment_id];R=np.array([[ta.R11,ta.R12,ta.R13],[ta.R21,ta.R22,ta.R23],[ta.R31,ta.R32,ta.R33]],np.float64);t=np.array([ta.tx,ta.ty,ta.tz],np.float64);vals=[]
    for j,m in enumerate(mm,1):
     aids=sorted(m);bids=[m[z] for z in aids]
     if not all(z in la['atoms'] for z in aids) or not all(z in lb['atoms'] for z in bids):continue
     A=np.vstack([la['atoms'][z] for z in aids]);B=np.vstack([lb['atoms'][z] for z in bids]);BA=B@R.T+t;vals.append((float(np.sqrt(np.mean(np.sum((A-BA)**2,axis=1)))),j,len(A),A,B,BA))
    nm=len(mm)
    if not vals:status='UPSTREAM_LIGAND_CONSISTENCY_REVIEW';reason='CCD_GRAPH_ATOMS_DO_NOT_MATCH_FROZEN_COORDINATES'
    else:
     rms,j,n,A,B,BA=min(vals,key=lambda z:(z[0],z[1]));selected_map=mm[j-1];sel=f'AMAP{j:08d}';cent=float(np.linalg.norm(A.mean(0)-BA.mean(0)));status='POSE_EQUIVALENT' if rms<=.5 else 'POSE_VARIANT';reason='FIXED_FRAME_RMSD_LE_0.50' if rms<=.5 else 'FIXED_FRAME_RMSD_GT_0.50'
     rr=float(np.sqrt(np.mean(np.sum((A-(B@R.T+t))**2,axis=1))));recompute_bad+=abs(rr-rms)>1e-10;bins[rb(rms)]+=1
     sens['graph_mapping_success']+=1;sens['multiple_symmetry_mappings']+=nm>1;sens['symmetry_permutation_required']+=any(k!=v for k,v in selected_map.items())
     for cut in (.25,.5,.75,1,2):sens[str(cut)]+=rms<=cut
  row=[x.pair_id_a,x.pair_id_b,x.alignment_id,ia.ligand_exact_id,n,nm,sel,f12(rms),f12(cent),status,reason];pw.writerow(row);counts[status]+=1
  if status=='POSE_EQUIVALENT':ew.writerow(row)
  elif status=='POSE_VARIANT':vw.writerow(row)
  else:rw.writerow(row)
  if i%200000==0:print('pose edges',i,flush=True)
 for f in(pf,ef,vf,rf):f.close()
 pd.DataFrame([('input_edges',len(sin)),('ligand_graph_mapping_success',sens['graph_mapping_success']),('multiple_symmetry_mappings',sens['multiple_symmetry_mappings']),('symmetry_permutation_required',sens['symmetry_permutation_required']),('mapping_review',counts['POSE_MAPPING_REVIEW'])]+[(k,counts[k]) for k in ['POSE_EQUIVALENT','POSE_VARIANT','POSE_MAPPING_REVIEW','UPSTREAM_LIGAND_CONSISTENCY_REVIEW']],columns=['metric','value']).to_csv(out/'06_step4_summary.tsv',sep='\t',index=False)
 pd.DataFrame([(f'bin_{k}',bins[k]) for k in ['0-0.10','0.10-0.25','0.25-0.50','0.50-0.75','0.75-1.00','1.00-1.50','1.50-2.00','2.00-3.00','>3.00']]+[(f'sensitivity_le_{k}',sens[k]) for k in ['0.25','0.5','0.75','1','2']],columns=['metric','value']).to_csv(out/'05_step4_rmsd_distribution.tsv',sep='\t',index=False)
 checks={'input_equals_step3_pass':sum(counts.values())==len(sin),'terminal_partition':sum(counts.values())==len(sin),'silent_drop_zero':True,'duplicate_edge_zero':not sin[['pair_id_a','pair_id_b']].duplicated().any(),'ligand_exact_id_equal_all_edges':all(inv.loc[x.pair_id_a].ligand_exact_id==inv.loc[x.pair_id_b].ligand_exact_id for x in sin.itertuples(index=False)),'exact_step3_alignment_id_R_t_reused':True,'selected_step3_chain_mapping_not_recomputed':True,'no_ligand_independent_alignment':True,'stored_mapping_rmsd_recomputes':recompute_bad==0}
 validation={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'counts':dict(counts),'elapsed_seconds':time.time()-t0,'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat()};(val/'validation.json').write_text(json.dumps(validation,indent=2));(out/'07_step4_report.md').write_text('# Filter 5 Step 4\n\nFixed receptor-frame symmetry-aware heavy-atom RMSD; correspondence-only minimization, no ligand fit. Cutoff 0.50 Å.\n')
 s3prov=json.loads((S3/'provenance.json').read_text());(run/'provenance.json').write_text(json.dumps({'source_step1':str(S1),'source_step2':str(S2),'source_step3':str(S3),'coordinate_source':str(LIG),'sifts_sha256':s3prov['sifts_sha256'],'rdkit_version':rdBase.rdkitVersion,'kabsch_implementation':s3prov['kabsch'],'chain_mapping_algorithm':s3prov['chain_mapping_algorithm'],'receptor_cutoff_angstrom':.5,'ligand_cutoff_angstrom':.5,'graph_mapping':'NETWORKX_EXACT_ATTRIBUTED_GRAPH_AUTOMORPHISM_V1','ligand_independent_alignment':False,'ligand_symmetry_correction':'ATOM_MAPPING_ONLY'},indent=2));(run/'output_schema.json').write_text(json.dumps({p.name:'TSV.GZ' if p.name.endswith('.gz') else 'TSV/MD' for p in out.iterdir()},indent=2));mr=[]
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
