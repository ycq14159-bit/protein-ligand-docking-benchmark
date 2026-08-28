# Comparable metric audit — harmonized v2 decision

The formal CATH metric is `unique valid four-level CATH H-level classifications`; null/unannotated values are excluded. Current CROWN metadata reproduces the harmonized expected value 2,040 exactly. The website value 2,041 is retained only as `CROWN_REPORTED_ONLY`.

The CROWN-reported PLINDER values Ion=22,728 and Artifact=18,626 are no longer targets for trial-filter reconstruction. They are retained only in `data/crown_reported_statistics_supplementary.tsv`.

Formal six-dataset ligand comparison uses `comparison_ligand_taxonomy_v1`: `monoatomic_ion_entries`, `simple_inorganic_entries`, and `shared_artifact_list_entries`. These are retrospective entry-level annotations and do not alter database membership.

See sibling directories `harmonized_cath_v1/` and `comparison_ligand_taxonomy_v1/` for row-level annotations, definitions, QC, and reproducible code.
