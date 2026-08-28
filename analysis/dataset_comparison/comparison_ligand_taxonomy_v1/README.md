# comparison_ligand_taxonomy_v1

Entry-level retrospective annotations for the six Mode B populations. No database membership changes.

- `monoatomic_ion_entries`: exactly one expected heavy atom whose element is outside C/H/dummy/N/O/P/S/Se/Te, matching the executable PLINDER v0.2.0 single-atom test. Formal charge is not required.
- `simple_inorganic_entries`: at least two expected heavy atoms and zero carbon atoms. This category excludes monoatomic entries.
- `shared_artifact_list_entries`: direct normalized CCD-code match to the frozen PLINDER v0.2.0 curated artifact list. It is an independent overlapping flag; synonym expansion is deliberately not inferred.

The first two labels are mutually exclusive. Artifact-list membership can overlap either and all intersections are reported. Graph-unavailable entries remain in the denominator and are reported in QC.
