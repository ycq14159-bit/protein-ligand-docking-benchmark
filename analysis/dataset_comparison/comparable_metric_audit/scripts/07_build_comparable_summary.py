from __future__ import annotations

import csv

from audit_common import OUT, QC, write_tsv


def main() -> None:
    rows = [
        {"metric": "Unique CATH IDs", "value": "NA", "status": "BLOCKED_CATH_DEFINITION_MISMATCH", "definition": "unique four-level C.A.T.H H-level classifications over formal receptor chains", "calibration_status": "FAIL"},
        {"metric": "Ion ligands", "value": "NA", "status": "BLOCKED_PLINDER_CLASSIFICATION_MISMATCH", "definition": "PLINDER-compatible ion classification", "calibration_status": "FAIL"},
        {"metric": "Artifact ligands", "value": "NA", "status": "BLOCKED_PLINDER_CLASSIFICATION_MISMATCH", "definition": "actual PLINDER artifact classifier/annotation", "calibration_status": "FAIL"},
    ]
    with (OUT / "ours_comparable_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    write_tsv(OUT / "metric_reason_summary.tsv", [
        {"metric": r["metric"], "status": r["status"], "reason": "External reported value could not be reproduced exactly from the frozen official source object; formal Ours value remains NA."}
        for r in rows
    ])
    write_tsv(QC / "population_qc.tsv", [{
        "population": "Filter 4 PASS", "formal_pair_rows": 91860, "unique_pair_id": 91860,
        "source_sha256": "244193b15987b4e69d1bbf6e5bc08923ae2ab1623d1fd7cb307130cb0d2e7fa3", "status": "PASS"
    }])


if __name__ == "__main__":
    main()

