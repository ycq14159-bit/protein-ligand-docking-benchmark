# Dataset comparison analysis

This directory contains auditable comparisons for the frozen 91,860-member Filter 4 PASS database and the six fixed Mode B populations. These analyses are retrospective annotations only; they do not change database membership or Filter 1–4 scientific rules.

## Formal harmonized metrics

- `Unique CATH IDs` is now the number of unique valid four-level CATH H-level classifications. Null/unannotated values are excluded. CROWN calibrates at 2,040; website 2,041 is reported-only.
- `comparison_ligand_taxonomy_v1` applies the same entry-level graph/list definitions to PDBbind, HiQBind, BioLiP2, PLINDER, CROWN and Ours:
  - monoatomic ion entries;
  - simple inorganic entries;
  - shared artifact-list entries.
- PLINDER Ion=22,728 and Artifact=18,626 from the CROWN website are retained only in `data/crown_reported_statistics_supplementary.tsv`; they are not reconstructed with trial filters.

Row-level derived Parquet files remain on the private server and are excluded from Git. Compact scripts, definitions, summaries, LaTeX tables, QC records and SHA256 manifests are tracked. See `harmonized_cath_v1/`, `comparison_ligand_taxonomy_v1/`, and `comparable_metric_audit/` before using the numbers.
