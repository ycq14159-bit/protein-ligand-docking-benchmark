# Filter 5: Strict Exact-Equivalent Redocking-Case Deduplication

The current frozen release `20260828_exact_only_v3_01` accepts 7,674 HIGH +
84,186 GOOD = 91,860 Filter 4 PASS pairs and retains 4,751 HIGH + 60,411 GOOD =
65,162 benchmark pairs.

Filter 5 v3 groups cases only when strict ligand identity and strict receptor
binding-site identity establish exact redocking-case equivalence. Ligand
identity preserves stereochemistry, formal charge, and bond order. Receptor
identity uses SIFTS/UniProt mapping plus direct-binding residue positions and
residue names. The representative is selected by resolution, then R-free,
absolute R-free/R-work gap, R-work, and pair ID.

Formal closure:

```text
12,417 exact representatives
+ 34,911 exact-unique retained
+ 17,834 review-retained
+ 26,698 exact redundant
= 91,860 input pairs
```

Review-retained cases are retained because exact redundancy was not proven;
they must not be described as confirmed unique. The former near-similarity
step using receptor similarity, ligand Tanimoto similarity, and 6 Å pocket
similarity has been removed from the membership policy.

Current implementation and policy snapshot:

```text
scripts/benchmark_filter5_v3_exact_only_finalize.py
config/filter5_v3_exact_only_policy.json
```
