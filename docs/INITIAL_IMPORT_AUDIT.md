# Initial import safety audit

Audit target: public GitHub repository. Result: **PASS**.

- Candidate files before Git metadata: 190
- Candidate bytes before Git metadata: 2,480,573
- Files over 20/50/100 MB: 0/0/0
- Blocked scientific-data/credential files: 0
- High-confidence secret matches: 0
- Non-placeholder credential-assignment matches: 0
- Binary files: 0
- Preserved legacy files containing historical absolute data paths: 94

Legacy absolute paths are provenance/portability findings, not credentials. They remain unchanged because this import does not rewrite pre-Git scientific code. New metadata uses `BENCHMARK_DATA_ROOT`.

## Twenty largest candidate files

| Path | Bytes |
|---|---:|
| `pipeline/filter_02/legacy_v2/scripts/filter2_v2_pipeline.py` | 73,231 |
| `scripts/figures/benchmark_figure_scout/build_benchmark_figure_scout.py` | 71,414 |
| `pipeline/filter_02/legacy_v1/scripts/filter2_pipeline.py` | 56,461 |
| `pipeline/filter_04/step_04/scripts/filter4_step4_pipeline.py` | 55,204 |
| `legacy/filter_04/step_04/step04_pilot_v2/executed_runner.py` | 55,204 |
| `legacy/filter_04/step_04/step04_full_v2/executed_runner.py` | 55,204 |
| `legacy/filter_04/step_04/step04_pilot_v1/executed_runner.py` | 52,552 |
| `legacy/filter_04/step_04/step04_full_v1/executed_runner.py` | 52,552 |
| `legacy/filter_04/step_04/step04_smoke_v1/executed_runner.py` | 51,703 |
| `scripts/data_engineering/entry_work_packages/build_entry_work_packages.py` | 51,219 |
| `pipeline/processing_02/scripts/processing2_pipeline.py` | 50,602 |
| `pipeline/filter_02/v3/scripts/filter2_v3_pipeline.py` | 48,987 |
| `pipeline/filter_01/scripts/filter1_pipeline.py` | 48,756 |
| `pipeline/filter_04/step_02/scripts/filter4_step2_pipeline.py` | 40,874 |
| `legacy/filter_04/step_02/step02_pilot_v3/executed_runner.py` | 40,874 |
| `legacy/filter_04/step_02/step02_full_v3/executed_runner.py` | 40,874 |
| `pipeline/filter_04/step_03/scripts/filter4_step3_pipeline.py` | 40,761 |
| `legacy/filter_04/step_03/step03_pilot_v2/executed_runner.py` | 40,761 |
| `legacy/filter_04/step_03/step03_full_v2/executed_runner.py` | 40,761 |
| `pipeline/processing_03/scripts/processing3_pipeline.py` | 40,539 |
