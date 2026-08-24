from pathlib import Path


def test_frozen_thresholds_are_present():
    text = Path(__file__).with_name("filter3_config.yaml").read_text()
    for token in ("3.0", "0.40", "0.45", "0.05", "0.80", "0.30"):
        assert token in text


def test_plip_and_filter4_are_out_of_scope():
    text = Path(__file__).with_name("filter3_config.yaml").read_text()
    assert "plip: true" in text
    assert "crystal_packing: true" in text
