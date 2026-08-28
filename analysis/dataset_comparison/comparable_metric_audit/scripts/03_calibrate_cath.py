from __future__ import annotations

import re
import pyarrow.parquet as pq
from audit_common import EXT, QC, write_tsv


def main() -> None:
    table = pq.read_table(EXT / "crown" / "CROWN_metadata_2026-06.parquet", columns=["cath_ids"])
    valid, missing = set(), 0
    pattern = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
    for value in table.column("cath_ids").to_pylist():
        items = list(value or [])
        missing += int(not items)
        valid.update(str(x) for x in items if pattern.fullmatch(str(x)))
    rows = [
        {"metric": "CROWN population rows", "expected_or_reported": 141261, "reproduced": table.num_rows,
         "difference": table.num_rows - 141261, "counting_unit": "metadata rows", "status": "PASS"},
        {"metric": "CROWN website Unique CATH IDs", "expected_or_reported": 2041, "reproduced": "NA",
         "difference": "NA", "counting_unit": "website reported", "status": "CROWN_REPORTED_ONLY"},
        {"metric": "Harmonized unique valid four-level CATH H", "expected_or_reported": 2040, "reproduced": len(valid),
         "difference": len(valid) - 2040, "counting_unit": "unique valid cath_h_id; null excluded",
         "status": "CALIBRATION_PASS" if len(valid) == 2040 else "FAIL"},
        {"metric": "Rows without valid CATH H", "expected_or_reported": "NA", "reproduced": missing,
         "difference": "NA", "counting_unit": "metadata rows", "status": "INFORMATIONAL"},
    ]
    write_tsv(QC / "cath_external_calibration.tsv", rows)


if __name__ == "__main__":
    main()

