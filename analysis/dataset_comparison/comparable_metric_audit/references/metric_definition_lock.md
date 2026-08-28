# Comparable metric definition lock

Access date: 2026-08-28 (Asia/Hong_Kong).

## Unique CATH IDs

- External metric: `Unique CATH IDs`.
- Locked operational definition: unique valid four-level `C.A.T.H` homologous-superfamily classification IDs over formal receptor chains.
- Domain-instance source for Ours: SIFTS; classification mapping source: official CATH-Plus v4.4.0.
- External target: current CROWN maintenance population, 141,261 rows and 2,041 reported IDs.
- Calibration finding: current CROWN metadata contains 2,040 valid four-level IDs. Adding one synthetic missing-CATH category produces 2,041, but a missing category is not a CATH H-level classification. Status is therefore `BLOCKED_CATH_DEFINITION_MISMATCH`.

## Ion ligands

- External metric: `Ion ligands`.
- Locked operational definition: the actual PLINDER release classifier, not formal charge and not observed one-heavy-atom count.
- PLINDER v0.2.0 code classifies a molecule as ion when it has one heavy atom and that atom is not matched by `[#6,#1,#0,#7,#8,#15,#16,#34,#52]` (C, H, dummy, N, O, P, S, Se, Te).
- External target: CROWN comparison reports PLINDER total 649,915 and ion 22,728.
- Calibration finding: neither official 2024-06/v2 annotation nor the diagnostic 2024-04/v1 annotation reproduces the population/counting unit. Status is `BLOCKED_PLINDER_CLASSIFICATION_MISMATCH`.

## Artifact ligands

- External metric: `Artifact ligands`.
- Locked operational definition: actual PLINDER `identify_artifact` implementation and release artifact list, not historical Filter 2 provenance.
- PLINDER v0.2.0 precedence is ion first (never artifact), then CCD/synonym artifact-list match, CCD dummy atoms, and `is_excluded_mol`.
- Actual v0.2.0 `is_excluded_mol` fails when heavy atoms are **less than 5** or carbon atoms are **less than 2**, absolute formal charge is over 2, or longest unbranched carbon linker is over 12. This means the executable boundary is HA >= 5 and C >= 2, not the prose `>5` and `>2` supplied in the task. The discrepancy is frozen rather than silently reconciled.
- External target: CROWN comparison reports PLINDER artifact 18,626.
- Calibration finding: official annotation objects do not reproduce the target. Status is `BLOCKED_PLINDER_CLASSIFICATION_MISMATCH`.

## Non-membership statement

These labels are retrospective cross-dataset comparison annotations. They do not alter database membership, scientific Filter 1–4 rules, or pair IDs. Historical proxies remain available as `INTERNAL_PROXY_NOT_COMPARABLE`: 58,346 SIFTS domain instances, 28 observed one-heavy-atom entries, and 3,526 historical simple-inorganic provenance entries.

