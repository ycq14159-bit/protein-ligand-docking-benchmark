
# Data management

Git and data are separate tracks. Git contains code and compact metadata. `${BENCHMARK_DATA_ROOT}` contains public-source snapshots, large tabular outputs, structures, intermediate work products, and frozen releases.

No raw mmCIF/PDB, SDF/MOL2, Parquet, SQLite, compressed run table, generated case directory, cache, log, or environment belongs in Git. This initial import does not use Git LFS or DVC.

Only small marker, validation, summary, and configuration files are hash-anchored in the manifests. A full hash of hundreds of gigabytes was intentionally not performed.
