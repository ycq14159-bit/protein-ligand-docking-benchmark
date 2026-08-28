from __future__ import annotations

import gc
import pyarrow.compute as pc
import pyarrow.parquet as pq

from audit_common import (EXT, PLINDER_REPORTED_ARTIFACT, PLINDER_REPORTED_ION,
                          PLINDER_REPORTED_TOTAL, QC, write_tsv)


def row(metric, release, unit, reproduced, reported):
    return {"metric": metric, "release": release, "counting_unit": unit, "reported": reported,
            "reproduced": int(reproduced), "difference": int(reproduced) - int(reported),
            "status": "CALIBRATION_PASS" if int(reproduced) == int(reported) else "BLOCKED_PLINDER_CLASSIFICATION_MISMATCH"}


def distinct_where(table, value_column, flag_column):
    return pc.count_distinct(pc.filter(table[value_column], table[flag_column])).as_py()


def true_count(table, column):
    return pc.sum(pc.cast(table[column], "int64")).as_py()


def main() -> None:
    p2 = EXT / "plinder" / "annotation_table_2024-06_v2.parquet"
    c2 = ["entry_pdb_id", "system_id", "system_id_no_biounit", "ligand_ccd_code",
          "ligand_is_ion", "ligand_is_artifact", "ligand_is_proper",
          "system_ligand_has_ion", "system_ligand_has_artifact", "system_type"]
    d2 = pq.read_table(p2, columns=c2)
    system_count_v2 = pc.count_distinct(d2["system_id"]).as_py()
    no_biounit_count_v2 = pc.count_distinct(d2["system_id_no_biounit"]).as_py()
    proper_system_count_v2 = distinct_where(d2, "system_id", "ligand_is_proper")
    ion_system_count_v2 = distinct_where(d2, "system_id", "system_ligand_has_ion")
    artifact_system_count_v2 = distinct_where(d2, "system_id", "system_ligand_has_artifact")
    rows = [
        row("Total", "2024-06/v2", "ligand annotation rows", d2.num_rows, PLINDER_REPORTED_TOTAL),
        row("Total", "2024-06/v2", "unique system_id", system_count_v2, PLINDER_REPORTED_TOTAL),
        row("Total", "2024-06/v2", "unique system_id_no_biounit", no_biounit_count_v2, PLINDER_REPORTED_TOTAL),
        row("Total", "2024-06/v2", "proper-ligand rows", true_count(d2, "ligand_is_proper"), PLINDER_REPORTED_TOTAL),
        row("Total", "2024-06/v2", "systems containing a proper ligand", proper_system_count_v2, PLINDER_REPORTED_TOTAL),
        row("Ion", "2024-06/v2", "ligand annotation rows", true_count(d2, "ligand_is_ion"), PLINDER_REPORTED_ION),
        row("Ion", "2024-06/v2", "unique systems with system_ligand_has_ion", ion_system_count_v2, PLINDER_REPORTED_ION),
        row("Artifact", "2024-06/v2", "ligand annotation rows", true_count(d2, "ligand_is_artifact"), PLINDER_REPORTED_ARTIFACT),
        row("Artifact", "2024-06/v2", "unique systems with system_ligand_has_artifact", artifact_system_count_v2, PLINDER_REPORTED_ARTIFACT),
    ]
    del d2
    gc.collect()
    p1 = EXT / "plinder" / "annotation_table_2024-04_v1.parquet"
    c1 = ["entry_pdb_id", "system_id", "system_type", "ligand_ccd_code", "ligand_is_ion", "ligand_is_artifact"]
    d1 = pq.read_table(p1, columns=c1)
    system_count_v1 = pc.count_distinct(d1["system_id"]).as_py()
    ion_system_count_v1 = distinct_where(d1, "system_id", "ligand_is_ion")
    artifact_system_count_v1 = distinct_where(d1, "system_id", "ligand_is_artifact")
    proper_mask = pc.and_(pc.invert(d1["ligand_is_ion"]), pc.invert(d1["ligand_is_artifact"]))
    proper_system_count_v1 = pc.count_distinct(pc.filter(d1["system_id"], proper_mask)).as_py()
    rows.extend([
        row("Total", "2024-04/v1", "ligand annotation rows", d1.num_rows, PLINDER_REPORTED_TOTAL),
        row("Total", "2024-04/v1", "unique system_id", system_count_v1, PLINDER_REPORTED_TOTAL),
        row("Total", "2024-04/v1", "non-ion/non-artifact systems", proper_system_count_v1, PLINDER_REPORTED_TOTAL),
        row("Ion", "2024-04/v1", "ligand annotation rows", true_count(d1, "ligand_is_ion"), PLINDER_REPORTED_ION),
        row("Ion", "2024-04/v1", "unique systems containing ion", ion_system_count_v1, PLINDER_REPORTED_ION),
        row("Artifact", "2024-04/v1", "ligand annotation rows", true_count(d1, "ligand_is_artifact"), PLINDER_REPORTED_ARTIFACT),
        row("Artifact", "2024-04/v1", "unique systems containing artifact", artifact_system_count_v1, PLINDER_REPORTED_ARTIFACT),
    ])
    write_tsv(QC / "ligand_class_external_calibration.tsv", rows)
    write_tsv(QC / "ion_audit.tsv", [r for r in rows if r["metric"] == "Ion"])
    write_tsv(QC / "artifact_audit.tsv", [r for r in rows if r["metric"] == "Artifact"])


if __name__ == "__main__":
    main()
