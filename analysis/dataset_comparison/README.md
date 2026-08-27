# Dataset comparison analysis

This directory builds an auditable CROWN-style property table for the current
formal database members.  The scientific inputs are read-only frozen releases;
the analysis does not change membership or scientific rules.

The authoritative population is the 91,860 unique `pair_id` records in the
frozen Filter 4 PASS inventory from
`20260826_filter3_118255_strict_posebusters_01`.  Filter 4 is explicitly marked
as `FINAL_DATABASE_CONSTRUCTION_STAGE`.  Filter 5 and Processing 4 are later
benchmark/docking-ready derivatives, not database membership.

Run `scripts/build_dataset_comparison.py`, then
`scripts/build_latex_table.py`.  Both accept explicit paths; no project path is
hard-coded.  External datasets are left `NA` because no confirmed reference
statistics were present in the repository.  See `qc/metric_definitions.tsv`,
`qc/source_provenance.tsv`, and `qc/authoritative_input_report.txt` before using
the table.

The entry-level Parquet, PDF, and logs are generated artifacts and are excluded
from Git under the repository data-management policy.  Compact scripts, CSV,
TeX, definitions, and QC records are tracked.
