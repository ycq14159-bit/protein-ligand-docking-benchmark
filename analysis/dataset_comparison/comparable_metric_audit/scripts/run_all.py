from __future__ import annotations

import runpy
from pathlib import Path


def main() -> None:
    here = Path(__file__).resolve().parent
    for name in [
        "01_inspect_external_definitions.py", "02_prepare_cath_mapping.py", "03_calibrate_cath.py",
        "04_calibrate_plinder_ligand_classes.py", "05_annotate_ours_cath.py",
        "06_annotate_ours_ligand_classes.py", "07_build_comparable_summary.py",
        "08_update_comparison_table.py",
    ]:
        runpy.run_path(str(here / name), run_name="__main__")


if __name__ == "__main__":
    main()

