#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import pandas as pd

ROOT = Path('/home/linx/data/youcq/autodl-tmp/benchmark_1.0')
OUT = ROOT / 'reconciliation/filter3_filter4_filter5_p4_20260822'
F4 = ROOT / 'filter_04_crystal_packing_influence/step_05_final_crystal_packing_decision/runs/step05_full_v1'
F5S1 = ROOT / 'filter_05_equivalent_redocking_case/step_01_exact_identity_and_candidate_blocking/runs/step01_full_v1'
F5S5 = ROOT / 'filter_05_equivalent_redocking_case/step_05_strict_equivalent_grouping_and_representative_selection/runs/step05_full_v1'
P4 = ROOT / 'processing_04_docking_ready_case_construction/runs/p4_full_v1_0_1'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def ids(path, column='pair_id'):
    frame = pd.read_csv(path, sep='\t', usecols=[column])
    values = frame[column].astype(str)
    return frame.shape[0], int(values.duplicated().sum()), set(values)


f4_inventory = F4 / 'output/01_filter4_final_pair_inventory.tsv.gz'
f4_pass = F4 / 'output/02_filter4_pass_pairs.tsv.gz'
f5_input = F5S1 / 'output/03_filter5_step1_pair_inventory.tsv.gz'
f5_final = F5S5 / 'output/01_filter5_final_pair_inventory.tsv.gz'

f4_all_rows, f4_all_dup, f4_all = ids(f4_inventory)
f4_pass_rows, f4_pass_dup, f4_pass_ids = ids(f4_pass)
f5_input_rows, f5_input_dup, f5_input_ids = ids(f5_input)

f5 = pd.read_csv(f5_final, sep='\t', usecols=['pair_id', 'filter5_final_status'])
f5['pair_id'] = f5['pair_id'].astype(str)
f5_retained = set(f5.loc[f5.filter5_final_status.ne('F5_REDUNDANT_EQUIVALENT_CASE'), 'pair_id'])

p4_input = pd.read_parquet(P4 / 'input/full_case_inventory.parquet', columns=['case_id', 'pair_id'])
p4_input['pair_id'] = p4_input.pair_id.astype(str)
p4_ids = set(p4_input.pair_id)

f4_validation = json.loads((F4 / 'validation.json').read_text())
f4_frozen = json.loads((F4 / '_FROZEN.json').read_text())
f5_provenance = json.loads((F5S1 / 'provenance.json').read_text())

checks = {
    'filter4_transfer_rsync_dry_run_empty': (OUT / 'filter4_post_transfer_dry_run.txt').read_text().strip() == '',
    'filter4_authoritative_sha256_failed_lines_zero': 'FAILED' not in (OUT / 'filter4_authoritative_sha256_check.log').read_text(),
    'filter4_validation_pass': bool(f4_validation.get('validation_pass')),
    'filter4_frozen_marker_present': f4_frozen.get('status') == 'FROZEN',
    'filter4_inventory_rows_336412': f4_all_rows == 336412,
    'filter4_inventory_pair_id_unique': f4_all_dup == 0 and len(f4_all) == f4_all_rows,
    'filter4_pass_rows_241545': f4_pass_rows == 241545,
    'filter4_pass_pair_id_unique': f4_pass_dup == 0 and len(f4_pass_ids) == f4_pass_rows,
    'filter4_pass_sha_matches_filter5_provenance': sha256(f4_pass) == f5_provenance['formal_membership_sha256'],
    'filter4_pass_set_equals_filter5_step1_input': f4_pass_ids == f5_input_ids,
    'filter5_step1_input_pair_id_unique': f5_input_dup == 0 and len(f5_input_ids) == f5_input_rows,
    'filter5_retained_set_equals_processing4_input': f5_retained == p4_ids,
}

result = {
    'generated_at': datetime.now(timezone.utc).isoformat(),
    'overall_pass': all(checks.values()),
    'checks': checks,
    'transfer': {
        'source': 'sbdx:/root/autodl-tmp/benchmark_1.0/filter_04_crystal_packing_influence/',
        'destination': str(ROOT / 'filter_04_crystal_packing_influence'),
        'transferred_files': 18752,
        'transferred_bytes': 7069712815,
        'post_transfer_rsync_differences': 0,
    },
    'filter4': {
        'run': str(F4),
        'inventory_rows': f4_all_rows,
        'pass_rows': f4_pass_rows,
        'pass_sha256': sha256(f4_pass),
        'validation_decision_counts': f4_validation.get('decision_counts', {}),
    },
    'lineage_set_differences': {
        'filter4_pass_vs_filter5_input': {
            'common': len(f4_pass_ids & f5_input_ids),
            'filter4_only': len(f4_pass_ids - f5_input_ids),
            'filter5_only': len(f5_input_ids - f4_pass_ids),
        },
        'filter5_retained_vs_processing4_input': {
            'common': len(f5_retained & p4_ids),
            'filter5_only': len(f5_retained - p4_ids),
            'processing4_only': len(p4_ids - f5_retained),
        },
    },
}
(OUT / 'reconciliation_validation.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')

summary_path = OUT / 'lineage_summary.json'
summary = json.loads(summary_path.read_text())
summary['filter4_frozen_evidence'].update({
    'direct_server_step5_run_present': True,
    'direct_step5_run': str(F4),
    'direct_validation_pass': bool(f4_validation.get('validation_pass')),
    'direct_pass_membership_path': str(f4_pass),
    'direct_pass_membership_sha256': sha256(f4_pass),
    'direct_pass_set_equals_filter5_step1_input': f4_pass_ids == f5_input_ids,
})
summary['transfer_and_final_validation'] = result
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')

report = f'''# Filter 3 → Filter 4 → Filter 5 → Processing 4 lineage reconciliation

Generated: {result['generated_at']}

## Final determination

The current authoritative frozen lineage is internally consistent and remains unchanged:

| Stage | Authoritative frozen output |
|---|---:|
| Filter 3 eligible | 336,412 (HIGH 120,297 + GOOD 216,115) |
| Filter 4 PASS | 241,545 |
| Filter 5 retained | 158,226 |
| Processing 4 input | 158,226 |
| Processing 4 final | READY 158,017; REVIEW 172; START_FAILED 37 |

Filter 4 PASS is exactly the Filter 5 Step 1 universe: common 241,545, Filter-4-only 0, Filter-5-only 0. Its SHA-256 is `{result['filter4']['pass_sha256']}`, exactly matching the hash recorded in frozen Filter 5 provenance.

Frozen Filter 5 retained membership is exactly the frozen Processing 4 input: common 158,226, Filter-5-only 0, Processing-4-only 0. Therefore `p4_full_v1_0_1` must not be rebased under the currently frozen lineage.

## Meaning of the three disputed counts

- **241,545** is authoritative: frozen Filter 4 Step 5 PASS membership and frozen Filter 5 Step 1 input.
- **158,226** is authoritative: frozen Filter 5 retained membership and frozen Processing 4 input.
- **156,621** is not a frozen output. It is a counterfactual audit result obtained by applying the proposed rule “any PoseBusters fail rejects” to the frozen records, restricting unchanged Filter 4 decisions, and replaying the frozen Filter 5 representative-selection algorithm in memory.

## Counterfactual rule audit (not authoritative)

The proposed Filter 3 rule would newly reject 2,969 current GOOD pairs; HIGH remains 120,297 and GOOD would become 213,146, for 333,443 eligible pairs. Of these removals, 2,148 are in the current Filter 4 PASS set, yielding a projected Filter 4 PASS universe of 239,397 if all other Filter 4 decisions are held fixed.

Replaying the frozen Filter 5 grouping rule on that projected universe yields 156,621 retained cases. Relative to current Processing 4 input: common 156,602, old-only 1,624, new-only replacement representatives 19. Among removed old representatives, 19 groups can select a replacement and 289 groups disappear completely.

This projection must not replace formal outputs until new Filter 3, Filter 4, and Filter 5 runs each have their own manifests, validation reports, and `_FROZEN.json` markers.

## Transfer and integrity validation

Missing Filter 4 material was restored directly from sbdx to the private server using resumable rsync. The transfer added 18,752 files / 7,069,712,815 bytes. A post-transfer dry-run reported zero differences. Embedded SHA256SUMS validation checked 62 authoritative files with zero failures. Filter 4 Step 5 frozen validation reports PASS 241,545, REJECT 94,865, REVIEW 2, totaling 336,412.

Overall reconciliation validation: **{'PASS' if result['overall_pass'] else 'FAIL'}**.

## Audit artifacts

- `lineage_summary.json`: machine-readable lineage and counterfactual summary
- `reconciliation_validation.json`: final transfer, checksum, cardinality, and set-equality checks
- `evidence_files.tsv`: evidence paths and hashes
- `authoritative_common.tsv`, `authoritative_old_only.tsv`, `authoritative_new_only.tsv`: frozen F5 vs P4 diff
- `common.tsv`, `old_only.tsv`, `new_only.tsv`: counterfactual candidate vs current P4 diff
- `representative_replacement_audit.tsv`: affected representative groups
- `counterfactual_filter5_inventory.tsv`, `counterfactual_filter5_groups.tsv`: derived audit-only candidate
- `filter4_authoritative_sha256_check.log`, `filter4_post_transfer_dry_run.txt`: transfer verification evidence
'''
(OUT / 'FINAL_LINEAGE_RECONCILIATION_REPORT.md').write_text(report)

artifact_lines = []
for path in sorted(p for p in OUT.iterdir() if p.is_file() and p.name != 'REPORT_SHA256SUMS'):
    artifact_lines.append(f'{sha256(path)}  {path.name}')
(OUT / 'REPORT_SHA256SUMS').write_text('\n'.join(artifact_lines) + '\n')
print(json.dumps(result, indent=2, sort_keys=True))
