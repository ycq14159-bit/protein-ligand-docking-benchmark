# Harmonized six-dataset property comparison

Mode B is **GO**. Historical CROWN replication remains a separate **HOLD** workstream and does not block the paper-primary harmonized comparison.

Locked Mode B datasets:

- PDBbind v2020: 19,443 official general-set protein-ligand complexes.
- HiQBind Figshare v3 small-molecule metadata: 31,572 rows.
- BioLiP2/Q-BioLiP: official base plus weekly annotations through 2026-06-26.
- PLINDER 2024-06/v2: full `index/annotation_table.parquet`.
- CROWN 2026-06: 141,261-row metadata release.
- Ours: frozen 20260826 Filter 4 PASS, 91,860 pairs.

Formal Mode B plotting populations are PDBbind 19,443; HiQBind 31,572;
BioLiP2/Q-BioLiP 548,615; PLINDER 616,723; CROWN 141,261; and Ours
91,860. BioLiP2 uses 521,146 eligible official-base interactions plus 27,469
`Relevant=yes` weekly interactions after ligand-dictionary membership and explicit
removal of nucleic-acid, k-mer, metal and single-atom-ion labels. The complete
raw annotation ledger is preserved outside Git.

Raw external data live outside Git under `E:/todo/external_datasets/` and `/home/linx/data/youcq/autodl-tmp/external_datasets/`. PDBbind raw and derived row-level material remain private and are never committed.

## Execution checkpoints

1. Acquisition: archive/file integrity, SHA256, row count and version lock.
2. Adapters: descriptor coverage, ligand-instance mapping and wwPDB validation coverage.
3. Figure: common RDKit 2025.09.5 definitions, Gaussian KDE `bw_adjust=0.85`, fixed axes and full-population missingness denominators.

The official PDBbind ligand SDF is preferred for descriptors; the accompanying
official MOL2 is used only when the SDF cannot be sanitized. Ligand validation
is accepted only after unique deposited-coordinate mapping. HiQBind uses exact
PDB+CCD+chain+residue mapping. BioLiP2 lacks residue numbers in the released
interaction key, so only unique PDB+CCD+chain matches are accepted and its high
RSR/RSCC missingness is shown rather than imputed.

The historical 19,449/31,573/86,458/649,915/153,005 values are retained only for future CROWN replication provenance.
