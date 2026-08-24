from filter3_v2_pipeline import classify_terminal


def base():
    return {
        "method": "xray", "technical_failure": False, "unavailable": [],
        "reject": [], "warnings": [], "resolution": 2.0,
        "ligand_occupancy": 1.0, "pocket_occupancy": 1.0,
        "posebusters_warning": False,
    }


def test_low_occupancy_is_good_not_reject():
    row = base()
    row["ligand_occupancy"] = 0.5
    status, _, warnings = classify_terminal(row)
    assert status == "FILTER3_GOOD_QUALITY"
    assert "LIGAND_OCCUPANCY_WARNING" in warnings


def test_missing_occupancy_is_good_not_unavailable():
    row = base()
    row["pocket_occupancy"] = None
    status, _, warnings = classify_terminal(row)
    assert status == "FILTER3_GOOD_QUALITY"
    assert "POCKET_OCCUPANCY_UNAVAILABLE" in warnings


def test_hard_failure_rejects():
    row = base()
    row["reject"] = ["LIGAND_DENSITY_QUALITY_FAIL"]
    assert classify_terminal(row)[0] == "FILTER3_REJECT"


def test_missing_hard_metric_is_unavailable():
    row = base()
    row["unavailable"] = ["POCKET_DENSITY_METRIC_UNAVAILABLE"]
    assert classify_terminal(row)[0] == "FILTER3_VALIDATION_DATA_UNAVAILABLE"


def test_non_xray_is_pending_not_reject():
    row = base()
    row["method"] = "cryo_em"
    assert classify_terminal(row)[0] == "FILTER3_NON_XRAY_PROTOCOL_PENDING"

