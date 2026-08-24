
# Pipeline

The benchmark follows two unit-preserving lanes before pair construction: receptor qualification and ligand qualification/placement. Processing 2 prepares frozen assembly coordinates and ligand topology; Processing 3 constructs ordinary non-covalent protein–ligand pairs; Filter 3 applies local experimental-quality rules; Filter 4 audits crystallographic neighbours; Filter 5 applies strict equivalence grouping; Processing 4 constructs docking-ready case files.

Authoritative counts and exact run paths are maintained in `manifests/frozen_runs.yaml`. Do not infer authority from the lexically latest directory name.

The currently frozen Filter 3 does not implement the later approved “any PoseBusters fail rejects” policy. The derived 156,621 population is counterfactual until a complete versioned rerun is frozen and validated.
