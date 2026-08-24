
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_required_metadata_exists():
    for relative in [
        "manifests/data_locations.yaml",
        "manifests/frozen_runs.yaml",
        "docs/DATA_MANAGEMENT.md",
        "docs/VERSIONING.md",
    ]:
        assert (ROOT / relative).is_file()


def test_no_bulk_scientific_data_extensions():
    blocked = {".parquet", ".cif", ".mmcif", ".pdb", ".sdf", ".mol2", ".sqlite", ".gz"}
    bad = [str(path.relative_to(ROOT)) for path in ROOT.rglob("*") if path.is_file() and path.suffix.lower() in blocked]
    assert bad == []
