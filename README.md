
# Protein–small molecule docking benchmark

This repository is the code-and-provenance track for an offline structural-bioinformatics benchmark. The scientific data track remains outside Git under `BENCHMARK_DATA_ROOT`.

Current authoritative frozen lineage:

`Processing 1: 144,408 eligible entries` → `Filter 1: 142,049 receptor-qualified entries` → `Filter 2: 236,383 ligand records` → `Processing 2: 206,218 ligand placements` → `Processing 3: 176,900 pairs` → `Filter 3: 118,255` → `Filter 4 PASS: 91,860` → `Filter 5 retained: 65,162` → `Processing 4 docking-ready: 64,100`

The current benchmark branch applies strict PoseBusters rejection in Filter 3,
retains only exact-equivalent deduplication in Filter 5, and treats Processing 4
as packaging with one explicit scope gate: ligand heavy-atom count must be greater
than 3. The 1,062 Filter 5 cases outside that scope remain recorded rather than
being silently dropped.

Git adoption began on **2026-08-24**. Files named `v1`, `v2`, `v3`, `full_v*`, or `executed_runner` that predate this date are legacy snapshots; they are not reconstructed Git history.

## Repository scope

Git manages source code, configuration, rules, validation/audit/analysis scripts, small provenance manifests, and documentation. Raw structures, derived coordinates, Parquet datasets, run outputs, caches, environments, and docking-ready case files remain in the data track and are excluded by policy.

Set `BENCHMARK_DATA_ROOT` before running code against a local data installation. Some preserved legacy scripts contain historical absolute paths; do not treat them as portable entry points without an explicit versioned port.

Historical implementations remain for provenance, but current entry points and
frozen-run anchors are identified in `manifests/frozen_runs.yaml`. No software
or data license is granted by this repository unless a license is added explicitly.
