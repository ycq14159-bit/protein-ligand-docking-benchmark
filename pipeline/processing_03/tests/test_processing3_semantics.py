import importlib.util
from pathlib import Path

script = Path(__file__).with_name("processing3_pipeline.py")
spec = importlib.util.spec_from_file_location("processing3_pipeline", script)
p3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p3)


def test_clash_boundary():
    assert p3.severe_steric_clash(1.999999)
    assert not p3.severe_steric_clash(2.0)


def test_vdw_contact_bounds():
    vdw = {"C": 1.7, "N": 1.55}
    upper = 1.7 + 1.55 + 0.5
    assert not p3.qualifying_direct_contact(1.999, "C", "N", vdw)
    assert p3.qualifying_direct_contact(2.0, "C", "N", vdw)
    assert p3.qualifying_direct_contact(upper, "C", "N", vdw)
    assert not p3.qualifying_direct_contact(upper + 1e-6, "C", "N", vdw)
    assert not p3.qualifying_direct_contact(3.0, "X", "N", vdw)


def test_altloc_not_in_terminal_key_schema():
    assert "alt_id" not in p3.OUTPUT_COLUMNS["placement_terminal_status"]
    assert "ligand_assembly_placement_id" in p3.OUTPUT_COLUMNS["placement_terminal_status"]


def test_nonprotein_covalent_partner_is_not_unresolved():
    import pandas as pd
    from types import SimpleNamespace
    meta = SimpleNamespace(component_id="LIG", label_asym_id="L", auth_asym_id="A", auth_seq_id="10", operator_path="1")
    receptor = pd.DataFrame([{"label_asym_id": "P", "auth_asym_id": "A", "operator_path": "1", "chain_instance_id": "protein-1"}])
    conn = {"conn_id": "c1", "conn_type_id": "covale",
            "p1": {"label_comp_id": "LIG", "auth_comp_id": "LIG", "label_asym_id": "L", "auth_asym_id": "A", "auth_seq_id": "10", "symmetry": "1_555"},
            "p2": {"label_comp_id": "OTH", "auth_comp_id": "OTH", "label_asym_id": "X", "auth_asym_id": "A", "auth_seq_id": "11", "symmetry": "1_555"}}
    mapped, unresolved = p3.ligand_declared_covalent(meta, receptor, [conn])
    assert mapped == []
    assert unresolved == []
