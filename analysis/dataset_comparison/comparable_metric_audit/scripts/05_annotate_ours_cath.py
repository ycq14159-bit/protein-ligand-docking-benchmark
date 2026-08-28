from __future__ import annotations

from audit_common import OUT, QC, write_tsv


def main() -> None:
    write_tsv(QC / "ours_cath_annotation_status.tsv", [{
        "status": "NOT_RUN_BLOCKED", "reason": "BLOCKED_CATH_DEFINITION_MISMATCH",
        "output_created": False,
        "note": "External calibration reconstructed 2040 valid H-level IDs, not the reported 2041; Ours annotation is prohibited by the task STOP rule."
    }])
    for path in (OUT / "ours_cath_mapping.parquet", OUT / "ours_comparable_annotations.parquet"):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite potentially misleading blocked output: {path}")


if __name__ == "__main__":
    main()

