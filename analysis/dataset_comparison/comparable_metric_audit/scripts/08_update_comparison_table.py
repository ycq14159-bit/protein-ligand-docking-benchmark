from __future__ import annotations

from audit_common import QC, write_tsv


def main() -> None:
    write_tsv(QC / "table_update_status.tsv", [
        {"metric": "Unique CATH IDs", "action": "UNCHANGED_NA", "reason": "CATH calibration failed"},
        {"metric": "Ion ligands", "action": "UNCHANGED_NA", "reason": "PLINDER calibration failed"},
        {"metric": "Artifact ligands", "action": "UNCHANGED_NA", "reason": "PLINDER calibration failed"},
        {"metric": "LaTeX/Overleaf", "action": "NOT_REGENERATED", "reason": "No metric passed calibration"},
    ])


if __name__ == "__main__":
    main()

