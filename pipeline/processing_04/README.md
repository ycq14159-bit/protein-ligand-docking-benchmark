# Processing 4: Docking-Ready Case Construction

The current frozen release `p4_benchmark_f5v3_65162_gt3_v3_1_0_01` converts
65,162 Filter 5 v3 retained pairs into docking-ready case files. Processing 4
is a preparation stage, not a new scientific quality filter. Its only scope
gate is ligand heavy-atom count greater than 3.

Formal closure:

- input: 4,751 HIGH + 60,411 GOOD = 65,162 pairs;
- READY: 4,313 HIGH + 59,787 GOOD = 64,100 cases;
- out of scope (`heavy_atom_count <= 3`): 438 HIGH + 624 GOOD = 1,062;
- preparation review: 0;
- ligand-start generation failure: 0.

Each READY case contains `receptor.cif`, `receptor.pdb`,
`ligand_reference.sdf`, `ligand.smi`, `ligand_start.sdf`, `site.json`, and
`metadata.json`. The ligand start is independent of the native pose: wwPDB CCD
ideal coordinates are preferred (63,949 cases); deterministic ETKDGv3 + UFF is
used for 149 cases, and ETKDGv3 without UFF parameters for 2 cases. Native
ligand coordinates are forbidden as a fallback.

The canonical site uses the native ligand heavy-atom bounding-box center, a
5 Å margin on each side, and a 22.5 Å minimum dimension. Method-specific
docking configuration remains outside Processing 4.

Current implementation and executed configuration:

```text
scripts/processing4_benchmark_f5v3_exact_only.py
config/processing4_benchmark_f5v3_exact_only.json
```

Set `BENCHMARK_DATA_ROOT` (or the stage-specific
`PROCESSING4_BENCHMARK_ROOT`) to the external scientific-data root before
running. The repository itself never contains case coordinates or Parquet
outputs.

Older scripts are retained as historical implementations and must not be
mistaken for the current frozen entry point.
