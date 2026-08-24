# Filter 3 → Filter 4 → Filter 5 → Processing 4 lineage reconciliation

Generated: 2026-08-22T11:10:10.472886+00:00

## Final determination

The current authoritative frozen lineage is internally consistent and remains unchanged:

| Stage | Authoritative frozen output |
|---|---:|
| Filter 3 eligible | 336,412 (HIGH 120,297 + GOOD 216,115) |
| Filter 4 PASS | 241,545 |
| Filter 5 retained | 158,226 |
| Processing 4 input | 158,226 |
| Processing 4 final | READY 158,017; REVIEW 172; START_FAILED 37 |

Filter 4 PASS is exactly the Filter 5 Step 1 universe: common 241,545, Filter-4-only 0, Filter-5-only 0. Its SHA-256 is `e2e6db2012c4766c719c5fc5f9028b76604a0c98dc319941e20a90b7f243d3a0`, exactly matching the hash recorded in frozen Filter 5 provenance.

Frozen Filter 5 retained membership is exactly the frozen Processing 4 input: common 158,226, Filter-5-only 0, Processing-4-only 0. Therefore `p4_full_v1_0_1` must not be rebased under the currently frozen lineage.

## Meaning of the three disputed counts

- **241,545** is authoritative: frozen Filter 4 Step 5 PASS membership and frozen Filter 5 Step 1 input.
- **158,226** is authoritative: frozen Filter 5 retained membership and frozen Processing 4 input.
- **156,621** is not a frozen output. It is a counterfactual audit result obtained by applying the proposed rule “any PoseBusters fail rejects” to the frozen records, restricting unchanged Filter 4 decisions, and replaying the frozen Filter 5 representative-selection algorithm in memory.

## Counterfactual rule audit (not authoritative)

The proposed Filter 3 rule would newly reject 2,969 current GOOD pairs; HIGH remains 120,297 and GOOD would become 213,146, for 333,443 eligible pairs. Of these removals, 2,148 are in the current Filter 4 PASS set, yielding a projected Filter 4 PASS universe of 239,397 if all other Filter 4 decisions are held fixed.

Replaying the frozen Filter 5 grouping rule on that projected universe yields 156,621 retained cases. Relative to current Processing 4 input: common 156,602, old-only 1,624, new-only replacement representatives 19. Among removed old representatives, 19 groups can select a replacement and 289 groups disappear completely.

This projection must not replace formal outputs until new Filter 3, Filter 4, and Filter 5 runs each have their own manifests, validation reports, and `_FROZEN.json` markers.

## Transfer and integrity validation

Missing Filter 4 material was restored directly from sbdx to the private server using resumable rsync. The transfer added 18,752 files / 7,069,712,815 bytes. A post-transfer dry-run reported zero differences. Embedded SHA256SUMS validation checked 62 authoritative files with zero failures. Filter 4 Step 5 frozen validation reports PASS 241,545, REJECT 94,865, REVIEW 2, totaling 336,412.

Overall reconciliation validation: **PASS**.

## Audit artifacts

- `lineage_summary.json`: machine-readable lineage and counterfactual summary
- `reconciliation_validation.json`: final transfer, checksum, cardinality, and set-equality checks
- `evidence_files.tsv`: evidence paths and hashes
- `authoritative_common.tsv`, `authoritative_old_only.tsv`, `authoritative_new_only.tsv`: frozen F5 vs P4 diff
- `common.tsv`, `old_only.tsv`, `new_only.tsv`: counterfactual candidate vs current P4 diff
- `representative_replacement_audit.tsv`: affected representative groups
- `counterfactual_filter5_inventory.tsv`, `counterfactual_filter5_groups.tsv`: derived audit-only candidate
- `filter4_authoritative_sha256_check.log`, `filter4_post_transfer_dry_run.txt`: transfer verification evidence
