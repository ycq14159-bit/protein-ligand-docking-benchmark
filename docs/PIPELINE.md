
# Pipeline

The benchmark follows two unit-preserving lanes before pair construction: receptor qualification and ligand qualification/placement. Processing 2 prepares frozen assembly coordinates and ligand topology; Processing 3 constructs ordinary non-covalent protein–ligand pairs; Filter 3 applies local experimental-quality rules (including direct rejection on any evaluated PoseBusters failure); Filter 4 audits crystallographic neighbours; Filter 5 removes strict exact-equivalent redocking cases; Processing 4 constructs docking-ready case files.

Authoritative counts and exact run paths are maintained in `manifests/frozen_runs.yaml`. Do not infer authority from the lexically latest directory name.

Current frozen benchmark lineage:

`176,900 Processing 3 pairs` → `118,255 Filter 3 HIGH+GOOD` → `91,860 Filter 4 PASS` → `65,162 Filter 5 retained` → `64,100 Processing 4 READY`

Filter 5 v3 no longer performs the former near-similarity grouping based on receptor similarity, ligand Tanimoto similarity, or 6 Å pocket-environment similarity. Processing 4 prefers wwPDB CCD ideal coordinates for the independent ligand start and uses deterministic RDKit ETKDGv3 only when the official ideal conformer cannot be used. Native ligand coordinates are never a ligand-start fallback.
