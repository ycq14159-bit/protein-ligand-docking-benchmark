# CATH / ion / artifact comparable-metric audit

This directory audits three formerly `NA` cross-dataset metrics for the frozen Filter 4 PASS population (91,860 unique protein-ligand pairs). It performs no membership filtering and changes no frozen scientific output.

## Outcome

All three formal values remain `NA` because the required external calibration did not pass:

- CROWN metadata has 2,040 valid four-level CATH H IDs, while the maintenance comparison reports 2,041. The reported value is obtained only by treating missing CATH as an additional category, which violates the locked definition.
- PLINDER 2024-06/v2 has 1,357,906 ligand rows, 990,260 systems and 484,562 systems containing a proper ligand; none is the reported 649,915 population. Neither row-level nor system-level ion/artifact counts reproduce 22,728/18,626. The 2024-04/v1 diagnostic also fails.
- The executable PLINDER v0.2.0 element set and size boundaries differ from the simplified prose in the task, so implementation was not guessed.

The STOP rule is enforced by scripts 05 and 06: they record `NOT_RUN_BLOCKED` and do not create Ours annotation parquet files. Script 08 records that the formal CSV/LaTeX/Overleaf products remain unchanged.

## Run

Place the frozen external objects under `external/{cath,crown,plinder}` using the manifests, then run:

```bash
python scripts/run_all.py
```

Large external data and derived parquet files are intentionally ignored by Git. Compact manifests, scripts, calibration tables and hashes are versioned.

## Scientific scope

PLINDER/CROWN-compatible ion and artifact labels are retrospective comparison annotations and do not alter database membership. The old values 58,346, 28 and 3,526 are retained only as `INTERNAL_PROXY_NOT_COMPARABLE`. The separate 317 non-binding pocket side-chain warnings are not reprocessed here.

