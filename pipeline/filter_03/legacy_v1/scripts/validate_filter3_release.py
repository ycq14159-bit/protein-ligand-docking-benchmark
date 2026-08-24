#!/usr/bin/env python3
import json
import os
from collections import Counter
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.dataset as ds

run_id = os.environ.get('FILTER3_RUN_ID', '20260812_full_01')
run = Path('/root/autodl-tmp/benchmark_1.0/filter_03_ground_truth_structure_quality/runs') / run_id
dataset = ds.dataset(run / 'output/filter3_pair_quality', format='parquet', partitioning='hive')
table = dataset.to_table(columns=['pair_id','ligand_assembly_placement_id','experimental_method_class','validation_mapping_status','terminal_status','decision','destination'])
rows = table.to_pylist()
status = Counter(row['terminal_status'] for row in rows)
pair_ids = [row['pair_id'] for row in rows]
placements = [row['ligand_assembly_placement_id'] for row in rows]
nonxray_bad = [row for row in rows if row['experimental_method_class'] != 'xray' and row['terminal_status'] not in {'FILTER3_NON_XRAY_PROTOCOL_PENDING','FILTER3_TECHNICAL_FAILURE'}]
mapping_quality_reject = [row for row in rows if row['validation_mapping_status'] != 'VALIDATION_MAPPING_OK' and row['terminal_status'] == 'FILTER3_REJECT']
checks = {
 'row_count': len(rows),
 'pair_id_unique': len(pair_ids) == len(set(pair_ids)),
 'placement_id_unique': len(placements) == len(set(placements)),
 'status_counts': dict(status),
 'status_sum': sum(status.values()),
 'nonxray_incorrect_terminal_count': len(nonxray_bad),
 'mapping_failure_as_quality_reject_count': len(mapping_quality_reject),
 'frozen_exists': (run / '_FROZEN.json').exists(),
 'release_validation': json.loads((run / 'release/filter3_release_validation.json').read_text()),
 'output_manifest_rows': sum(1 for _ in (run / 'release/output_manifest.tsv').open()) - 1,
}
checks['validation_pass'] = (
 checks['row_count'] == 744580 and checks['pair_id_unique'] and checks['placement_id_unique']
 and checks['status_sum'] == 744580 and checks['nonxray_incorrect_terminal_count'] == 0
 and checks['mapping_failure_as_quality_reject_count'] == 0 and checks['frozen_exists']
 and checks['release_validation']['validation_pass'] is True
)
print(json.dumps(checks, indent=2, sort_keys=True))
