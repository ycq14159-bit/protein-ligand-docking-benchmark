
# Versioning

Git history starts on 2026-08-24. Pre-existing `v1`/`v2`/`v3`/`full_v*` scripts are preserved as legacy source snapshots rather than rewritten into artificial historical commits.

Future changes should use normal commits and versioned runs. A frozen scientific run must name its code commit, configuration, input manifests, output manifest, validation report, and immutable data-relative path. Do not create a frozen-stage Git tag unless the repository contents are proven to match the code that produced that run.

The initial import tag is `initial-import-20260824`. It marks Git adoption only; it is not a scientific release tag.
