from __future__ import annotations

import pandas as pd

from audit_common import EXT, QC, write_tsv


def main() -> None:
    source = EXT / "cath" / "cath-domain-list-v4_4_0.txt"
    records = []
    with source.open(encoding="ascii") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            domain, c, a, t, h = fields[:5]
            records.append((domain, f"{c}.{a}.{t}.{h}"))
    frame = pd.DataFrame(records, columns=["cath_domain_instance", "cath_h_id"])
    frame.to_parquet(EXT / "cath" / "cath_domain_to_h.parquet", index=False)
    write_tsv(QC / "cath_mapping_qc.tsv", [{
        "release": "CATH-Plus v4.4.0", "mapping_rows": len(frame),
        "unique_domain_instances": frame.cath_domain_instance.nunique(),
        "unique_cath_h_id": frame.cath_h_id.nunique(),
        "duplicate_domain_instances": int(frame.cath_domain_instance.duplicated().sum()),
        "status": "PASS" if not frame.cath_domain_instance.duplicated().any() else "FAIL",
    }])


if __name__ == "__main__":
    main()

