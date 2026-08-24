# Processing 4: Docking-Ready Case Construction

Processing 4 converts the 158,226 Filter 5 retained pairs into self-contained
docking cases. It is preparation, not a new scientific filter. It never reruns
assembly identification, contact detection, PLIP, crystal-packing analysis,
quality filtering, or equivalent-case grouping.

Each successful case contains `receptor.cif`, `receptor.pdb`,
`ligand_reference.sdf`, `ligand.smi`, `ligand_start.sdf`, `site.json`, and
`metadata.json`. Native ligand coordinates are used only for the reference and
site definition. The start conformer is generated from a coordinate-free copy
of the frozen graph with ETKDGv3 and UFF.

The runner is resumable at case and bucket level. It writes a case atomically,
then records `_SUCCESS.json`. Review/failure cases remain in the inventory and
are never silently dropped.

Typical commands:

```bash
python scripts/processing4_pipeline.py prepare-run --run-dir runs/p4_smoke100_v1 --mode smoke --limit 100
python scripts/processing4_pipeline.py run --run-dir runs/p4_smoke100_v1 --workers 2
python scripts/processing4_pipeline.py validate --run-dir runs/p4_smoke100_v1
```

The pilot site policy in `config/processing4_v1.json` must be frozen after the
smoke validation before the full run is launched.
