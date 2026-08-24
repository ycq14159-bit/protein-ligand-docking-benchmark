
# Reproducibility

1. Check out a named Git commit or tag.
2. Set `BENCHMARK_DATA_ROOT` to the separately provisioned data root.
3. Resolve the intended stage/run through `manifests/frozen_runs.yaml` rather than directory discovery.
4. Verify the listed compact anchors before analysis.
5. Use the environment documentation/configuration associated with the stage.
6. Never overwrite a frozen run; create a new versioned run directory and validation record.

Legacy code is preserved verbatim and may contain non-portable historical paths. Porting must be a new reviewed commit, not an in-place rewrite of the imported snapshot.
