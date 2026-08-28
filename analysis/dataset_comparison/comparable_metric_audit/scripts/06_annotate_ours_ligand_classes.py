from __future__ import annotations

from audit_common import OUT, QC, write_tsv


def main() -> None:
    write_tsv(QC / "classification_partition_qc.tsv", [{
        "formal_pairs": 91860, "classification_rows": "NA", "ion": "NA", "artifact": "NA", "proper": "NA",
        "status": "NOT_RUN_BLOCKED", "reason": "BLOCKED_PLINDER_CLASSIFICATION_MISMATCH"
    }])
    write_tsv(QC / "manual_review_sample.tsv", [{
        "status": "NOT_RUN_BLOCKED", "reason": "No Ours classification is permitted before exact external calibration."
    }])
    path = OUT / "ours_ligand_classification.parquet"
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite potentially misleading blocked output: {path}")


if __name__ == "__main__":
    main()

