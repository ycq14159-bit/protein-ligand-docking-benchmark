from __future__ import annotations

from audit_common import QC, write_tsv


def main() -> None:
    write_tsv(QC / "ligand_class_external_calibration.tsv", [
        {"metric": "PLINDER Total entries", "reported": 649915, "status": "CROWN_REPORTED_ONLY",
         "formal_harmonized_replacement": "Mode B proper-ligand population N=616723"},
        {"metric": "PLINDER Ion ligands", "reported": 22728, "status": "CROWN_REPORTED_ONLY",
         "formal_harmonized_replacement": "comparison_ligand_taxonomy_v1/monoatomic_ion_entries"},
        {"metric": "PLINDER Artifact ligands", "reported": 18626, "status": "CROWN_REPORTED_ONLY",
         "formal_harmonized_replacement": "comparison_ligand_taxonomy_v1 explicit metrics"},
    ])


if __name__ == "__main__":
    main()
