#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

P2 = Path('/root/autodl-tmp/benchmark_1.0/processing_2_assembly_ready_structure_preparation/runs/20260810_full_01/output')
P3 = Path('/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/runs/20260811_full_01/output')
PLIP = Path('/root/miniconda3/envs/plip-audit/bin/plip')
PAIR_ROOT = P3 / 'provisional_pairs'
CHAIN_CODES = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
RUN_DIR = None
TIMEOUT = 300
KEEP_XML = True

RESULT_COLUMNS = [
    'pair_id', 'ligand_assembly_placement_id', 'pdb_id', 'assembly_id', 'model_id', 'component_id',
    'receptor_chain_instance_ids', 'input_protein_atom_count', 'input_ligand_atom_count', 'input_chain_count',
    'construction_status', 'plip_status', 'exit_code', 'runtime_seconds', 'interaction_count',
    'hydrophobic_count', 'hydrogen_bond_count', 'water_bridge_count', 'salt_bridge_count',
    'pi_stack_count', 'pi_cation_count', 'halogen_bond_count', 'metal_complex_count',
    'other_interaction_count', 'raw_xml_gz_path', 'error_message', 'membership_effect'
]


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_parquet(path, frame):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.parquet.tmp')
    frame = frame.reindex(columns=RESULT_COLUMNS)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), tmp, compression='zstd')
    os.replace(tmp, path)


def read_parquet(path, columns=None):
    return pq.ParquetFile(path).read(columns=columns).to_pandas()


def atom_name(name, element):
    name = str(name).strip()[:4]
    if len(name) < 4 and len(str(element).strip()) == 1:
        return f' {name:<3}'
    return f'{name:<4}'


def pdb_line(record, serial, atom, resname, chain, resseq, x, y, z, occ, b, element):
    return (f'{record:<6}{serial:5d} {atom_name(atom, element)} {str(resname)[:3]:>3} {chain:1}'
            f'{resseq:4d}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}'
            f'{float(occ):6.2f}{float(b):6.2f}          {str(element).strip().upper():>2}\n')


def make_complex(pair, rec, lig, bonds, path):
    chains = [x for x in str(pair.receptor_chain_instance_ids).split(',') if x]
    if not chains or len(chains) >= len(CHAIN_CODES):
        return None, 'unsupported_receptor_chain_count'
    temporary_ligand_resname = str(pair.component_id) if len(str(pair.component_id)) <= 3 else 'LIG'
    chain_map = dict(zip(chains, CHAIN_CODES))
    ligand_chain = CHAIN_CODES[len(chains)]
    lines = []
    serial = 1
    serial_by_ligand_atom = {}
    residue_maps = {}
    for chain_id in chains:
        subset = rec[rec['filter_1_chain_instance_id'].astype(str).eq(chain_id)]
        if subset.empty:
            return None, f'missing_receptor_chain:{chain_id}'
        residue_order = []
        for row in subset.itertuples(index=False):
            key = (str(row.label_seq_id), str(row.auth_seq_id), str(row.insertion_code), str(row.label_comp_id))
            if key not in residue_order:
                residue_order.append(key)
        residue_maps[chain_id] = {key: i + 1 for i, key in enumerate(residue_order)}
        if len(residue_order) > 9999:
            return None, f'receptor_chain_residue_count_exceeds_pdb_limit:{chain_id}'
        for row in subset.itertuples(index=False):
            key = (str(row.label_seq_id), str(row.auth_seq_id), str(row.insertion_code), str(row.label_comp_id))
            lines.append(pdb_line('ATOM', serial, row.label_atom_id, row.label_comp_id, chain_map[chain_id],
                                  residue_maps[chain_id][key], row.Cartn_x, row.Cartn_y, row.Cartn_z,
                                  row.occupancy if pd.notna(row.occupancy) else 1.0,
                                  row.B_iso_or_equiv if pd.notna(row.B_iso_or_equiv) else 0.0, row.type_symbol))
            serial += 1
    if lig.empty:
        return None, 'missing_ligand_atoms'
    for row in lig.itertuples(index=False):
        lines.append(pdb_line('HETATM', serial, row.label_atom_id, temporary_ligand_resname, ligand_chain, 1,
                              row.Cartn_x, row.Cartn_y, row.Cartn_z,
                              row.occupancy if pd.notna(row.occupancy) else 1.0,
                              row.B_iso_or_equiv if pd.notna(row.B_iso_or_equiv) else 0.0, row.type_symbol))
        serial_by_ligand_atom[str(row.label_atom_id)] = serial
        serial += 1
    for row in bonds.itertuples(index=False):
        a = serial_by_ligand_atom.get(str(row.atom_id_1)); b = serial_by_ligand_atom.get(str(row.atom_id_2))
        if a and b:
            order = {'SING': 1, 'DOUB': 2, 'TRIP': 3}.get(str(row.bond_order).upper(), 1)
            lines.append(f'CONECT{a:5d}' + f'{b:5d}' * order + '\n')
    lines.append('END\n')
    path.write_text(''.join(lines), encoding='ascii')
    return {'protein_atoms': len(rec[rec['filter_1_chain_instance_id'].astype(str).isin(chains)]),
            'ligand_atoms': len(lig), 'chain_count': len(chains)}, ''


def parse_xml(path):
    root = ET.parse(path).getroot()
    tags = Counter()
    known = {
        'hydrophobic_interaction': 'hydrophobic_count', 'hydrogen_bond': 'hydrogen_bond_count',
        'water_bridge': 'water_bridge_count', 'salt_bridge': 'salt_bridge_count', 'pi_stack': 'pi_stack_count',
        'pi_cation_interaction': 'pi_cation_count', 'halogen_bond': 'halogen_bond_count',
        'metal_complex': 'metal_complex_count'
    }
    counts = Counter()
    for element in root.iter():
        tag = element.tag.split('}')[-1]
        if tag in known:
            counts[known[tag]] += 1
            tags[tag] += 1
    total = sum(counts.values())
    return total, counts


def compress_xml(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix('.gz.tmp')
    with src.open('rb') as inp, gzip.open(tmp, 'wb', compresslevel=6) as out:
        shutil.copyfileobj(inp, out)
    os.replace(tmp, dst)


def process_partition(rel_text):
    rel = Path(rel_text)
    pairs = read_parquet(PAIR_ROOT / rel)
    result_path = Path(RUN_DIR) / 'work/results' / rel
    marker = Path(RUN_DIR) / 'work/checkpoints' / (rel.stem + '.json')
    lig_atoms = read_parquet(P2 / 'prepared_ligand_assembly_atoms' / rel)
    rec_atoms = read_parquet(P2 / 'prepared_receptor_assembly_atoms' / rel)
    bonds_path = P2 / 'prepared_ligand_assembly_bonds' / rel
    bonds = read_parquet(bonds_path) if bonds_path.exists() else pd.DataFrame(columns=['filter_2_ligand_assembly_placement_id'])
    rows = []
    tmp_root = Path(RUN_DIR) / 'work/tmp' / str(os.getpid())
    tmp_root.mkdir(parents=True, exist_ok=True)
    for pair in pairs.itertuples(index=False):
        started = time.time()
        lid = str(pair.ligand_assembly_placement_id)
        la = lig_atoms[lig_atoms['filter_2_ligand_assembly_placement_id'].astype(str).eq(lid)]
        wanted = [x for x in str(pair.receptor_chain_instance_ids).split(',') if x]
        ra = rec_atoms[rec_atoms['filter_1_chain_instance_id'].astype(str).isin(wanted)]
        lb = bonds[bonds['filter_2_ligand_assembly_placement_id'].astype(str).eq(lid)] if not bonds.empty else bonds
        safe = hashlib.sha256(str(pair.pair_id).encode()).hexdigest()[:24]
        work = tmp_root / safe
        shutil.rmtree(work, ignore_errors=True); work.mkdir(parents=True)
        pdb = work / 'complex.pdb'
        metrics, error = make_complex(pair, ra, la, lb, pdb)
        base = {c: '' for c in RESULT_COLUMNS}
        for c in ['input_protein_atom_count', 'input_ligand_atom_count', 'input_chain_count', 'exit_code',
                  'runtime_seconds', 'interaction_count', 'hydrophobic_count', 'hydrogen_bond_count',
                  'water_bridge_count', 'salt_bridge_count', 'pi_stack_count', 'pi_cation_count',
                  'halogen_bond_count', 'metal_complex_count', 'other_interaction_count']:
            base[c] = 0
        base.update({'pair_id': pair.pair_id, 'ligand_assembly_placement_id': lid, 'pdb_id': pair.pdb_id,
                     'assembly_id': pair.assembly_id, 'model_id': pair.model_id, 'component_id': pair.component_id,
                     'receptor_chain_instance_ids': pair.receptor_chain_instance_ids, 'membership_effect': False})
        if error:
            base.update({'construction_status': 'failure', 'plip_status': 'input_construction_failure',
                         'exit_code': 99, 'runtime_seconds': time.time() - started, 'error_message': error})
            rows.append(base); shutil.rmtree(work, ignore_errors=True); continue
        base.update({'construction_status': 'success', 'input_protein_atom_count': metrics['protein_atoms'],
                     'input_ligand_atom_count': metrics['ligand_atoms'], 'input_chain_count': metrics['chain_count']})
        cmd = [str(PLIP), '-f', str(pdb), '-o', str(work), '-x', '-q', '--maxthreads', '1', '--name', 'report', '--nofixfile']
        try:
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=TIMEOUT)
            code = proc.returncode
            xmls = sorted(work.glob('*.xml'))
            if code != 0:
                base.update({'plip_status': 'execution_failed', 'exit_code': code, 'error_message': (proc.stderr or '')[:1000]})
            elif not xmls:
                base.update({'plip_status': 'no_xml_output', 'exit_code': code, 'error_message': (proc.stderr or '')[:1000]})
            else:
                try:
                    total, counts = parse_xml(xmls[0])
                    base.update(counts)
                    base['interaction_count'] = total
                    base['plip_status'] = 'success' if total else 'no_interaction'
                    base['exit_code'] = code
                    if KEEP_XML:
                        raw = Path(RUN_DIR) / 'raw_xml' / rel.parent / rel.stem / f'{safe}.xml.gz'
                        compress_xml(xmls[0], raw); base['raw_xml_gz_path'] = str(raw)
                except Exception as exc:
                    base.update({'plip_status': 'output_parse_failed', 'exit_code': code,
                                 'error_message': f'{type(exc).__name__}:{exc}'})
        except subprocess.TimeoutExpired:
            base.update({'plip_status': 'timeout', 'exit_code': 124, 'error_message': f'timeout>{TIMEOUT}s'})
        except Exception as exc:
            base.update({'plip_status': 'execution_exception', 'exit_code': 127,
                         'error_message': f'{type(exc).__name__}:{exc}'})
        base['runtime_seconds'] = time.time() - started
        rows.append(base)
        shutil.rmtree(work, ignore_errors=True)
    frame = pd.DataFrame(rows).reindex(columns=RESULT_COLUMNS)
    write_parquet(result_path, frame)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({'status': 'complete', 'partition': rel_text, 'pairs': len(pairs),
                                  'result_rows': len(frame), 'finished_at': utc()}) + '\n')
    return {'partition': rel_text, 'pairs': len(pairs), 'statuses': dict(Counter(frame['plip_status']))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--workers', type=int, default=20)
    parser.add_argument('--timeout', type=int, default=300)
    parser.add_argument('--selection-limit', type=int)
    args = parser.parse_args()
    global RUN_DIR, TIMEOUT
    RUN_DIR = args.run_dir; TIMEOUT = args.timeout
    run = Path(RUN_DIR); (run / 'work/checkpoints').mkdir(parents=True, exist_ok=True)
    rels = sorted(x.relative_to(PAIR_ROOT) for x in PAIR_ROOT.rglob('*.parquet'))
    if args.selection_limit:
        selected = []; total = 0
        for rel in rels:
            n = pq.ParquetFile(PAIR_ROOT / rel).metadata.num_rows
            if n and total < args.selection_limit:
                selected.append(rel); total += n
            if total >= args.selection_limit:
                break
        rels = selected
    pending = [str(r) for r in rels if not (run / 'work/checkpoints' / (r.stem + '.json')).exists()]
    status = {'status': 'RUNNING', 'workers': args.workers, 'timeout': TIMEOUT, 'partition_count': len(rels),
              'pending_partition_count': len(pending), 'selection_limit': args.selection_limit, 'started_at': utc()}
    (run / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
    done = pairs = 0; counts = Counter(); started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(process_partition, pending, chunksize=1):
            done += 1; pairs += result['pairs']; counts.update(result['statuses'])
            if done % 10 == 0 or done == len(pending):
                status.update({'completed_this_attempt': done, 'pairs_this_attempt': pairs,
                               'status_counts_this_attempt': dict(counts), 'runtime_seconds': time.time() - started,
                               'updated_at': utc()})
                (run / 'status.json').write_text(json.dumps(status, indent=2) + '\n')
                print(json.dumps(status), flush=True)
    status.update({'status': 'COMPLETED', 'completed_this_attempt': done, 'pairs_this_attempt': pairs,
                   'status_counts_this_attempt': dict(counts), 'runtime_seconds': time.time() - started, 'finished_at': utc()})
    (run / 'status.json').write_text(json.dumps(status, indent=2) + '\n')


if __name__ == '__main__':
    main()
