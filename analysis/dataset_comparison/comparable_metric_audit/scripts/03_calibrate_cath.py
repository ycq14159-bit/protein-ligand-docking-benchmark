from __future__ import annotations

import re
import pyarrow.parquet as pq

from audit_common import CROWN_REPORTED_CATH, EXT, QC, write_tsv


def main() -> None:
    table = pq.read_table(EXT / "crown" / "CROWN_metadata_2026-06.parquet", columns=["cath_ids"])
    valid = set()
    missing_rows = 0
    assignments = 0
    pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
    for value in table.column("cath_ids").to_pylist():
        items = list(value or [])
        if not items:
            missing_rows += 1
        assignments += len(items)
        valid.update(str(item) for item in items if pattern.fullmatch(str(item)))
    reproduced = len(valid)
    status = "CALIBRATION_PASS" if reproduced == CROWN_REPORTED_CATH else "BLOCKED_CATH_DEFINITION_MISMATCH"
    rows = [
        {"metric": "CROWN population rows", "reported": 141261, "reproduced": table.num_rows,
         "difference": table.num_rows - 141261, "counting_unit": "CROWN metadata rows", "status": "PASS" if table.num_rows == 141261 else "FAIL"},
        {"metric": "Unique CATH IDs", "reported": CROWN_REPORTED_CATH, "reproduced": reproduced,
         "difference": reproduced - CROWN_REPORTED_CATH, "counting_unit": "unique valid four-level cath_h_id", "status": status},
        {"metric": "Diagnostic including missing category", "reported": CROWN_REPORTED_CATH, "reproduced": reproduced + int(missing_rows > 0),
         "difference": reproduced + int(missing_rows > 0) - CROWN_REPORTED_CATH,
         "counting_unit": "unique cath_h_id plus one synthetic missing category", "status": "DIAGNOSTIC_ONLY"},
        {"metric": "CATH assignments", "reported": "", "reproduced": assignments, "difference": "",
         "counting_unit": "array elements", "status": "INFORMATIONAL"},
        {"metric": "Rows without CATH", "reported": "", "reproduced": missing_rows, "difference": "",
         "counting_unit": "CROWN metadata rows", "status": "INFORMATIONAL"},
    ]
    write_tsv(QC / "cath_external_calibration.tsv", rows)


if __name__ == "__main__":
    main()

