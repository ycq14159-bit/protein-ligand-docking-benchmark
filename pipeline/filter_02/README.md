# Filter 2: Ligand Recognition and Qualification

The current frozen release is
`filter_2_ligand_qualification_v4/runs/20260825_dual_source_strict_01`.
It accepts 142,049 receptor-qualified PDB entries, enumerates 533,610 source
ligand instances and 718,007 assembly placements before the new relevance
rule, and releases 183,904 source ligands plus 52,479 additional assembly
placements = 236,383 ligand records.

For `SUSPICIOUS_LIGAND` cases, membership requires all three conditions:

```text
BioLiP = EXACT_RETAINED
AND Q-BioLiP = EXACT_ASSEMBLY_CCD_CHAIN
AND Q-BioLiP Relevant = yes
```

Non-suspicious sources bypass this dual-source rule. All stopped suspicious
rows remain in the membership ledger with explicit provenance states. The
current v4 projection/finalization/freeze scripts are under `scripts/v4/`;
legacy v1-v3 implementations remain for provenance.
