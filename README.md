
# Protein–small molecule docking benchmark

This repository is the code-and-provenance track for an offline structural-bioinformatics benchmark. The scientific data track remains outside Git under `BENCHMARK_DATA_ROOT`.

Current authoritative frozen lineage:

`Filter 3: 336,412` → `Filter 4 PASS: 241,545` → `Filter 5 retained: 158,226` → `Processing 4 docking-ready: 158,017`

Git adoption began on **2026-08-24**. Files named `v1`, `v2`, `v3`, `full_v*`, or `executed_runner` that predate this date are legacy snapshots; they are not reconstructed Git history.

## Repository scope

Git manages source code, configuration, rules, validation/audit/analysis scripts, small provenance manifests, and documentation. Raw structures, derived coordinates, Parquet datasets, run outputs, caches, environments, and docking-ready case files remain in the data track and are excluded by policy.

Set `BENCHMARK_DATA_ROOT` before running code against a local data installation. Some preserved legacy scripts contain historical absolute paths; do not treat them as portable entry points without an explicit versioned port.

This initial import preserves legacy scientific logic. It does not rerun or modify any Filter/Processing result. No software or data license is granted by this repository unless a license is added explicitly.
