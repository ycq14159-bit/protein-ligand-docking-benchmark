#!/usr/bin/env python3
"""Append official-source cross-benchmark scout and finalize audit manifests."""

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

ROOT=Path('/home/linx/data/youcq/autodl-tmp/benchmark_1.0')
OUT=ROOT/'audits/benchmark_figure_scout/20260823_draft_01'
SOURCE_LOG=Path('/tmp/cross_benchmark_source_log.tsv')
RESEARCH_NOTE=Path('/tmp/cross_benchmark_research.md')
FIELD_NOTE=Path('/tmp/field_scout_findings.md')
COL={'blue':'#2563EB','orange':'#F59E0B','green':'#16A34A','red':'#DC2626','purple':'#7C3AED','gray':'#64748B','light':'#E2E8F0'}

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8<<20),b''): h.update(b)
    return h.hexdigest()

def utc(): return datetime.now(timezone.utc).isoformat()

def save(fig,stem):
    png=OUT/'cross_benchmark'/f'{stem}.png'; svg=OUT/'cross_benchmark'/f'{stem}.svg'
    fig.savefig(png,dpi=220,bbox_inches='tight',facecolor='white'); fig.savefig(svg,bbox_inches='tight',facecolor='white'); plt.close(fig)
    return str(png.relative_to(OUT))

def append(path,text):
    old=path.read_text() if path.exists() else ''
    path.write_text(old.rstrip()+"\n\n"+text.strip()+"\n")

def main():
    qc=json.loads((OUT/'qc/qc_report.json').read_text())
    if not qc.get('validation_pass'): raise RuntimeError('base Figure Scout QC did not pass')
    if not SOURCE_LOG.exists() or not RESEARCH_NOTE.exists(): raise RuntimeError('cross-benchmark research files missing')
    shutil.copy2(SOURCE_LOG,OUT/'research/cross_benchmark_source_log.tsv')
    shutil.copy2(RESEARCH_NOTE,OUT/'research/cross_benchmark_research.md')
    if FIELD_NOTE.exists(): shutil.copy2(FIELD_NOTE,OUT/'research/frozen_field_scout_findings.md')
    data=pd.read_csv(SOURCE_LOG,sep='\t',dtype=str).fillna('')
    final=pd.read_parquet(OUT/'data_derived/final_retained_characterization.parquet',columns=['pdb_id'])
    ours_pdb=final.pdb_id.nunique()

    # Scale panel uses facets and labels because statistical units differ.
    structural=['CROWN','PLINDER','PDBbind','HiQBind','BioLiP2 / Q-BioLiP']
    evaluation=['This benchmark','PoseBusters Benchmark','PoseX Self-Docking','Astex Diverse Set','CASF-2016','DockGen']
    counts={r.benchmark:int(r.reported_records_or_cases.replace('>','')) for r in data.itertuples()}
    counts['This benchmark']=158226
    fig,axes=plt.subplots(1,3,figsize=(16,6.2),gridspec_kw={'width_ratios':[1.15,1.15,1.25]})
    for ax,names,title in [(axes[0],structural,'Structural resources\n(heterogeneous record units)'),(axes[1],evaluation,'Evaluation benchmarks\n(case/complex units)')]:
        vals=[counts[n] for n in names]; colors=[COL['purple'] if n=='This benchmark' else COL['blue'] for n in names]
        display={'PoseBusters Benchmark':'PoseBusters','PoseX Self-Docking':'PoseX','Astex Diverse Set':'Astex'}
        y=np.arange(len(names)); ax.barh(y,vals,color=colors); ax.set_yticks(y,[display.get(n,n) for n in names]); ax.invert_yaxis(); ax.set_xscale('log'); ax.set_xlabel('Reported records/cases (log scale)'); ax.set_title(title,loc='left',fontweight='bold')
        ax.spines[['top','right']].set_visible(False); ax.grid(axis='x',color=COL['light'],lw=.7); ax.set_axisbelow(True)
        for yi,(n,v) in enumerate(zip(names,vals)):
            prefix='>' if n=='PLINDER' else ''
            if n in {'PLINDER','BioLiP2 / Q-BioLiP'}:
                ax.text(v*.96,yi,f'{prefix}{v:,}',va='center',ha='right',fontsize=8,color='white')
            else:
                ax.text(v*1.08,yi,f'{prefix}{v:,}',va='center',fontsize=8)
    levels=['overlap_entry','overlap_chemical','overlap_instance','overlap_system']
    level_names=['Entry','Chemical','Ligand instance','Binding system']
    mapping={'YES':3,'HIGH_PARTIAL':2.5,'PARTIAL':2,'NEED_DATA':0,'':0}
    arr=np.array([[mapping.get(str(row[c]),1) for c in levels] for _,row in data.iterrows()])
    cmap=ListedColormap(['#E2E8F0','#FDE68A','#FDBA74','#86EFAC'])
    im=axes[2].imshow(arr,aspect='auto',vmin=0,vmax=3,cmap=cmap)
    axes[2].set_xticks(range(4),level_names,rotation=30,ha='right'); axes[2].set_yticks(range(len(data)),data.benchmark); axes[2].set_title('Overlap feasibility from\nlightweight official metadata',loc='left',fontweight='bold')
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            axes[2].text(j,i,data.iloc[i][levels[j]].replace('_',' '),ha='center',va='center',fontsize=6.5)
    fig.suptitle('Cross-benchmark landscape — units and overlap levels kept explicit',x=.035,y=.98,ha='left',fontsize=15,fontweight='bold')
    fig.subplots_adjust(left=.07,right=.98,bottom=.17,top=.76,wspace=.38)
    output=save(fig,'fig07_cross_benchmark_landscape')

    overlap_def=pd.DataFrame([
        ['entry','normalized lowercase four-character PDB ID','Current lightweight scope'],
        ['chemical_full','full Standard InChIKey','Preserves stereo/protonation distinctions'],
        ['chemical_connectivity','first InChIKey block','Connectivity-only auxiliary view'],
        ['ligand_instance','PDB + model + chain/asym + residue + insertion + AltLoc','PARTIAL / NEED_DATA for most resources'],
        ['binding_system','assembly + operator/placement + receptor chain set + ligand instance','NEED_DATA for whole-suite comparison'],
    ],columns=['overlap_level','minimum_key','status_or_note'])
    overlap_def.to_csv(OUT/'cross_benchmark/overlap_definition.tsv',sep='\t',index=False)

    method_cols=['density_or_local_quality_qc','ligand_completeness','pocket_completeness','crystal_contact_handling','biological_assembly','strict_equivalence_dedup','evaluation_tasks']
    method=pd.DataFrame({'benchmark':['This benchmark']+data.benchmark.tolist()})
    for c in method_cols: method[c]='TO_VERIFY_FROM_METHODS'
    method.loc[0,method_cols[:-1]]='YES_FROZEN'; method.loc[0,'evaluation_tasks']='FUTURE_NO_PREDICTIONS'
    method.to_csv(OUT/'cross_benchmark/methodology_matrix_scout.tsv',sep='\t',index=False)

    inv=pd.read_csv(OUT/'figure_inventory.csv')
    m=inv.figure_id.eq('F07')
    inv.loc[m,['status','output_file','recommendation']]=['READY_PARTIAL',output,'KEEP']
    extra={'figure_id':'F07B','title':'Cross-benchmark methodological comparison','question':'Which curation and evaluation capabilities differ?','population':'Official methods documentation','metrics':', '.join(method_cols),'status':'NEED_METHODS_VERIFICATION','source':'official-source research log','output_file':'cross_benchmark/methodology_matrix_scout.tsv','recommendation':'NEED_DATA'}
    inv=pd.concat([inv,pd.DataFrame([extra])],ignore_index=True)
    inv.to_csv(OUT/'figure_inventory.csv',index=False)

    prov=pd.read_csv(OUT/'figure_provenance.csv')
    pm=prov.figure_id.eq('F07')
    prov.loc[pm,'status']='READY_PARTIAL'; prov.loc[pm,'output_file']=output; prov.loc[pm,'recommendation']='KEEP'; prov.loc[pm,'n_total']=len(data)+1; prov.loc[pm,'n_metric_available']=len(data)+1; prov.loc[pm,'n_missing']=0; prov.loc[pm,'denominator']='Each named official version; heterogeneous units shown in facets'; prov.loc[pm,'source_files']='research/cross_benchmark_source_log.tsv'; prov.loc[pm,'generated_timestamp']=utc()
    prov.to_csv(OUT/'figure_provenance.csv',index=False)

    append(OUT/'figure_review.md',f'''## Cross-benchmark review

F07 scale comparison is `KEEP`, but only in facets. It must never become a single ranked bar chart because structural resources and evaluation suites count different units. The overlap-feasibility matrix is more honest than a premature UpSet. Entry overlap is currently feasible for most resources; whole-suite ligand-instance and binding-system overlap remain `PARTIAL / NEED_DATA`.

The most serious public-version traps are PoseBusters 428 vs formal 308, HiQBind 32,275 total vs 31,572 current small-molecule rows, dynamic BioLiP2/Q-BioLiP snapshots, and DockGen 189 vs the later PoseBench 91. CASF-2016 should be described as the PDBbind v2016 core set. This benchmark has {ours_pdb:,} unique PDB entries among 158,226 frozen retained pair cases.''')
    append(OUT/'figure_review.md','''## Unexpected internal results worth following up

- Final-ligand concentration is the clearest new composition issue: SO4/GOL/EDO alone account for 35.0%, and the top six CCD identities (adding PO4/ACT/PEG) for 42.44%. Keep F05C and open a separate cognate-vs-additive audit; do not conflate it with the heavy-atom question.
- HIGH cases are removed more often than GOOD cases at both F4 and F5. Keep this as a descriptive panel and avoid causal language.
- Mixed HIGH/GOOD equivalence groups correctly select HIGH, but representatives are not better on continuous resolution/RSCC medians. Describe the actual rule, not an inferred optimization objective.
- F4 severity panels are scientifically useful only with their Step 3/4 eligibility denominators.''')
    append(OUT/'data_gap_report.md','''## Cross-benchmark gap details

Official lightweight metadata supports scale and entry-level work now. CROWN, PLINDER, HiQBind and PoseX also expose useful chemical fields, but a common normalized InChIKey freeze is still required. Most resources do not expose model/AltLoc/assembly/operator locators, so instance/system overlap cannot yet be claimed. The methods-comparison matrix is left as a verification scaffold rather than filled by assumption.''')
    append(OUT/'qc_report.md','''## Additional provenance findings

- Filter 1 has release validation and checksums but no single unified `_FROZEN.json`.
- Filter 4 Step 1 is a validated legacy run without its own frozen marker; Step 5 pins its inputs and hashes.
- Filter 4/5 lack top-level `CURRENT_RUN.json`, so the Figure Scout fixes exact run paths.
- A failed Filter 2 freeze-attempt marker remains under an audit directory; automatic discovery must not glob every `*FROZEN*` file.
- Cross-benchmark scale values were not treated as a common denominator.''')

    qc['cross_benchmark']={'official_source_rows':len(data),'ours_unique_pdb_in_frozen_final':int(ours_pdb),'scale_unit_mixing_prevented':True,'instance_or_system_overlap_claimed':False,'methodology_matrix_status':'NEED_METHODS_VERIFICATION'}
    (OUT/'qc/qc_report.json').write_text(json.dumps(qc,indent=2,default=str)+'\n')
    shutil.copy2(Path(__file__),OUT/'scripts'/Path(__file__).name)

    checks={'base_qc_pass':qc['validation_pass'],'figure_inventory_unique':inv.figure_id.is_unique,'all_ready_outputs_exist':True,'cross_source_rows_10':len(data)==10,'frozen_final_158226':len(final)==158226}
    for row in inv.itertuples():
        if str(row.status).startswith('READY') and str(row.output_file): checks['all_ready_outputs_exist'] &= (OUT/str(row.output_file)).exists()
    final_qc={'generated_at':utc(),'validation_pass':bool(all(checks.values())),'checks':{k:bool(v) for k,v in checks.items()},'figure_count':len(inv)}
    (OUT/'qc/finalization_validation.json').write_text(json.dumps(final_qc,indent=2)+'\n')
    if not final_qc['validation_pass']: raise SystemExit(2)

    rows=[]
    for p in sorted(x for x in OUT.rglob('*') if x.is_file() and x.name not in {'output_manifest.tsv','SHA256SUMS'}):
        rows.append({'relative_path':str(p.relative_to(OUT)),'size_bytes':p.stat().st_size,'sha256':sha256(p)})
    pd.DataFrame(rows).to_csv(OUT/'output_manifest.tsv',sep='\t',index=False)
    with (OUT/'SHA256SUMS').open('w') as h:
        for p in sorted(x for x in OUT.rglob('*') if x.is_file() and x.name!='SHA256SUMS'):
            h.write(f'{sha256(p)}  {p.relative_to(OUT)}\n')
    print(json.dumps({'output':str(OUT),'validation_pass':True,'figures':len(inv),'ours_unique_pdb':int(ours_pdb)},indent=2))

if __name__=='__main__': main()
