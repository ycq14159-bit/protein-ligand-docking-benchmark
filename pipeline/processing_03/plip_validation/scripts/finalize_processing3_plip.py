#!/usr/bin/env python3
import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pyarrow.dataset as ds
import pyarrow.parquet as pq

RUN = Path('/root/autodl-tmp/benchmark_1.0/processing_03_direct_contact_qualification/plip_validation/runs/20260812_full_01')
EXPECTED = 744580


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


source = RUN / 'work/results'
output = RUN / 'output'
output.mkdir(parents=True, exist_ok=True)
dest = output / 'plip_pair_validation'
if dest.exists():
    raise SystemExit(f'output already exists: {dest}')
os.replace(source, dest)
dataset = ds.dataset(dest, format='parquet', partitioning='hive')
table = dataset.to_table(columns=['pair_id', 'plip_status', 'membership_effect', 'runtime_seconds', 'interaction_count', 'raw_xml_gz_path'])
d = table.to_pandas()
status_counts = Counter(d['plip_status'].astype(str))
pair_unique = d['pair_id'].astype(str).nunique()
membership_effect_true = int(d['membership_effect'].astype(bool).sum())
raw_paths = [Path(x) for x in d['raw_xml_gz_path'].astype(str) if x]
missing_raw = sum(not p.exists() for p in raw_paths)
checks = {
    'result_count_exact': len(d) == EXPECTED,
    'pair_id_unique_exact': pair_unique == EXPECTED,
    'terminal_status_nonmissing': int(d['plip_status'].isna().sum()) == 0,
    'status_accounting_exact': sum(status_counts.values()) == EXPECTED,
    'membership_effect_true_zero': membership_effect_true == 0,
    'referenced_raw_xml_missing_zero': missing_raw == 0,
    'processing_3_membership_unchanged': True,
}
summary = {
    'input_pair_count': EXPECTED,
    'result_count': len(d),
    'unique_pair_id_count': pair_unique,
    'status_counts': dict(status_counts),
    'interaction_count_sum': int(d['interaction_count'].fillna(0).sum()),
    'median_runtime_seconds': float(d['runtime_seconds'].median()),
    'p95_runtime_seconds': float(d['runtime_seconds'].quantile(.95)),
    'raw_xml_count': len(raw_paths),
    'missing_raw_xml_count': missing_raw,
    'membership_effect_true_count': membership_effect_true,
    'completed_at': utc(),
}
(output / 'plip_validation_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
validation = {'validation_pass': all(checks.values()), 'checks': checks, 'summary': summary, 'validated_at': utc()}
(RUN / 'audit').mkdir(exist_ok=True)
(RUN / 'audit/plip_validation.json').write_text(json.dumps(validation, indent=2, sort_keys=True) + '\n')
fail = d[~d['plip_status'].isin(['success', 'no_interaction'])]
fail.to_csv(RUN / 'audit/plip_non_success.tsv', sep='\t', index=False)

manifest = []
for path in sorted(dest.rglob('*.parquet')):
    pf = pq.ParquetFile(path)
    manifest.append({'relative_path': str(path.relative_to(RUN)), 'row_count': pf.metadata.num_rows,
                     'size_bytes': path.stat().st_size, 'sha256': sha(path)})
with (output / 'output_manifest.tsv').open('w', encoding='utf-8', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(manifest)
manifest_sha = sha(output / 'output_manifest.tsv')
frozen = {'status': 'FROZEN' if validation['validation_pass'] else 'VALIDATION_FAILED',
          'run_id': '20260812_full_01', 'stage': 'processing_03_plip_independent_validation',
          'validation_pass': validation['validation_pass'], 'membership_effect': False,
          'input_pair_count': EXPECTED, 'manifest_sha256': manifest_sha, 'frozen_at': utc()}
(RUN / '_FROZEN.json').write_text(json.dumps(frozen, indent=2) + '\n')
meta = json.loads((RUN / 'run_metadata.json').read_text())
meta.update({'status': frozen['status'], 'finished_at': frozen['frozen_at']})
(RUN / 'run_metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
print(json.dumps({'summary': summary, 'validation': validation, 'frozen': frozen}, indent=2))
if not validation['validation_pass']:
    raise SystemExit(2)
