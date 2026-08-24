import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).with_name("processing2_pipeline.py")
spec = importlib.util.spec_from_file_location("processing2_pipeline", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_blank_altloc_and_first_tie_winner():
    module.selftest(type("Args", (), {})())


def test_strict_completeness_policy_is_not_crown_missing_o_exception():
    assert "missing_CCD_heavy_atoms_equals_zero" in SCRIPT.read_text()
    assert "crown_missing_O_exception_enabled: false" not in SCRIPT.read_text()

