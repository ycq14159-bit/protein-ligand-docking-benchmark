import importlib.util
from pathlib import Path

import numpy as np
from rdkit import Chem


SCRIPT = Path(__file__).parents[1] / "scripts/processing4_pipeline.py"
SPEC = importlib.util.spec_from_file_location("processing4_pipeline", SCRIPT)
P4 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P4)


def test_stable_bucket_known_value():
    assert P4.stable_bucket("13zy") == 0


def test_case_id_is_safe_and_deterministic():
    pair = "P3|13zy|1|example"
    assert P4.case_id(pair, "13zy") == P4.case_id(pair, "13zy")
    assert "|" not in P4.case_id(pair, "13zy")


def test_start_generation_does_not_need_native_coordinates():
    mol = Chem.MolFromSmiles("N[C@@H](C)C(=O)O")
    cfg = {"etkdg_random_seed": 24301, "etkdg_enforce_chirality": True,
           "uff_max_iterations": 500}
    start, _ = P4.independent_start(mol, cfg)
    assert start.GetNumConformers() == 1
    assert P4.graph_key(start) == P4.graph_key(mol)


def test_unspecified_coordinate_perceived_stereo_is_normalized_to_frozen_source():
    source = Chem.MolFromSmiles("COP(=O)([O-])O[P@@](=O)(O)OC")
    perceived = Chem.MolFromSmiles("CO[P@@](=O)([O-])O[P@@](=O)(O)OC")
    assert P4.graph_key(perceived, source) == P4.graph_key(source)


def test_site_box_covers_native_ligand():
    import pandas as pd
    atoms = pd.DataFrame({"type_symbol": ["C", "N"], "Cartn_x": [-2.0, 3.0],
                          "Cartn_y": [1.0, 4.0], "Cartn_z": [0.0, 2.0]})
    cfg = {"site_margin_angstrom_per_side": 5.0,
           "site_minimum_dimension_angstrom": 20.0,
           "site_center": "native_ligand_heavy_atom_bounding_box_center",
           "site_definition_version": "test"}
    site = P4.native_site(atoms, cfg)
    assert site["site_center"] == {"x": 0.5, "y": 2.5, "z": 1.0}
    assert all(site["search_box"][k] >= 20.0 for k in ("size_x", "size_y", "size_z"))
