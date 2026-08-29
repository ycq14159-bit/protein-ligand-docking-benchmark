from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

def main(root:Path,run:Path):
    d=pd.read_parquet(run/'reconstruction/reconstructed_entry_gate_for_f2_source_pdbs.parquet')
    base=set(d.loc[d.reconstructed_status.eq('ELIGIBLE'),'pdb_id']); reason={p:'MMCIF_ELIGIBLE_XML_NOT_REQUIRED_OR_CONFIRMED' for p in base}
    for sub in ['validation_recovery','validation_recovery_reject']:
        z=pd.read_csv(run/sub/'validation_xml_recovery.tsv.gz',sep='\t')
        for r in z.itertuples(index=False):
            if bool(r.xml_gate_pass): base.add(r.pdb_id);reason[r.pdb_id]='XML_GATE_PASS_RECOVERY'
    for sub in ['validation_recovery_pass_boundary','validation_recovery_pass_remaining']:
        z=pd.read_csv(run/sub/'validation_xml_recovery.tsv.gz',sep='\t')
        for r in z.itertuples(index=False):
            if bool(r.xml_gate_pass): base.add(r.pdb_id);reason[r.pdb_id]='XML_GATE_PASS_CONFIRMED'
            elif r.xml_status=='FAILED': base.add(r.pdb_id);reason[r.pdb_id]='MMCIF_GATE_PASS_REMOTE_XML_403_HISTORICAL_CLOSURE'
            else: base.discard(r.pdb_id);reason.pop(r.pdb_id,None)
    rows=pd.DataFrame({'pdb_id':sorted(base)});rows['scope_recovery_reason']=rows.pdb_id.map(reason)
    f2=root/'filter_2_ligand_qualification_v3/runs/20260804_full_01/output'
    s=pd.read_csv(f2/'provisional_source_ligands.tsv.gz',sep='\t',usecols=['pdb_id','source_ligand_instance_id'],dtype=str);s.pdb_id=s.pdb_id.str.lower();ss=s[s.pdb_id.isin(base)]
    p=pd.read_csv(f2/'ligand_assembly_logical_placements.tsv.gz',sep='\t',usecols=['pdb_id','assembly_ligand_placement_id'],dtype=str);p.pdb_id=p.pdb_id.str.lower();pp=p[p.pdb_id.isin(base)]
    validation={'source_containing_entries':len(rows),'sources':len(ss),'source_unique':ss.source_ligand_instance_id.nunique(),'placements':len(pp),'placement_unique':pp.assembly_ligand_placement_id.nunique(),'anchors':{'entries_106842':len(rows)==106842,'sources_533610':len(ss)==533610,'placements_718007':len(pp)==718007},'remote_xml_403_fallback_ids':sorted(rows.loc[rows.scope_recovery_reason.eq('MMCIF_GATE_PASS_REMOTE_XML_403_HISTORICAL_CLOSURE'),'pdb_id']),'explicit_xml_metrics_missing_excluded':['8ymg']}
    if not all(validation['anchors'].values()):raise RuntimeError(validation)
    rows.to_parquet(run/'reconstruction/recovered_source_containing_entry_scope_106842.parquet',index=False)
    (run/'reconstruction/recovered_scope_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',type=Path,required=True);a.add_argument('--run',type=Path,required=True);x=a.parse_args();main(x.root,x.run)
