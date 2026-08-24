
# Phase 0 engineering inventory

Read-only filesystem audit on 2026-08-24 found **2,628,665 files** totaling **234.379 GiB** under the data root. An extension-only selection was unsafe because it identified 598,728 JSON files, most of which are generated per-case records. The import therefore uses semantic source/config directories and explicit legacy runner snapshots.

Imported source/config/document files: **166** (2.30 MiB), including **69** explicitly marked pre-Git legacy files.

Imported file types: `.java` 7, `.json` 6, `.md` 4, `.pre-portable-path-20260821t095107z` 1, `.py` 109, `.sh` 21, `.xml` 1, `.yaml` 17.

Pipeline/import distribution: `configs` 1, `filter_01` 8, `filter_02` 25, `filter_03` 25, `filter_04` 27, `legacy` 33, `processing_01` 14, `processing_02` 3, `processing_03` 6, `processing_04` 9, `scripts` 14, `shared` 1.

## Deliberately excluded

- 787,804 compressed files, 255,966 Parquet files, 202,582 PDB files, 202,531 CIF files, 370,800 SDF files, and 186,062 SMILES files.
- Run `input/`, `output/`, `work/`, `logs/`, `tmp/`, cache, checkpoint, and environment trees.
- SQLite databases, binary JARs, generated figures, derived tables, and docking-ready case directories.
- Vendored CROWN reference code under Processing 2, pending an explicit public-repository license review; benchmark-authored Processing 2 code is included.
- Duplicate notebook checkpoints and a malformed duplicate PLIP directory containing a carriage-return path component.

The detailed one-to-one source import is recorded in `manifests/source_import.tsv`.
