"""Inline tests for canonical MPAS grid fixtures.

Validates the declarative ``.esm`` grid instances and the family schema
committed under ``discretizations/grids/mpas/``. Fixtures are small
declarative configs per the 2026-04-20 mayor scope correction (bead
``dsc-3nw``); per-binding runtimes derive geometry on demand from the
referenced mesh file, so these tests cover structure, the GRIDS_API §8
schema shape, and the closed-sphere Euler identities — not accessor output.
Accessor conformance is covered separately in ``test_mpas.py`` and the
cross-binding tests under ``tests/conformance/grids/``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_GRID_DIR = _REPO_ROOT / "discretizations" / "grids"
_FIXTURE_DIR = _GRID_DIR / "mpas"
_SCHEMA_PATH = _GRID_DIR / "mpas.schema.json"

_FIXTURE_NAMES = ("x1.642", "x1.2562", "x1.10242")


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


# --- Schema -------------------------------------------------------------------


def test_schema_file_exists():
    assert _SCHEMA_PATH.is_file()


def test_schema_matches_grids_api_section_8_shape():
    schema = _load(_SCHEMA_PATH)
    assert schema["family"] == "mpas"
    assert schema["version"] == "1.0.0"
    assert "loader" in schema["required"]

    opts = schema["options"]
    for key in ("loader", "R", "ghosts", "dtype"):
        assert key in opts, f"schema must declare option {key!r}"
    assert opts["R"]["units"] == "m"
    assert opts["R"]["default"] == 6371000.0
    assert opts["dtype"]["default"] == "float64"
    assert opts["ghosts"]["default"] == 0

    tols = schema["tolerances"]
    assert tols["default_rel"] == 1e-14
    # Cell-area / arc-length fields relax to 1e-12 per `docs/GRIDS_API.md` §4.2.
    for field in ("area", "dc_edge", "dv_edge"):
        assert field in tols["per_field"]

    for binding in ("julia", "python", "rust", "ts"):
        assert binding in schema["ecosystem_hints"]


# --- Fixture set --------------------------------------------------------------


def test_exactly_three_canonical_fixtures_committed():
    assert _FIXTURE_DIR.is_dir()
    found = sorted(p.name for p in _FIXTURE_DIR.glob("*.esm"))
    assert found == sorted(f"{n}.esm" for n in _FIXTURE_NAMES)


@pytest.mark.parametrize("name", _FIXTURE_NAMES)
def test_declarative_shape(name):
    esm = _load(_FIXTURE_DIR / f"{name}.esm")
    assert esm["family"] == "mpas"
    assert esm["version"] == "1.0.0"
    assert esm["schema_version"] == "1.0.0"
    assert esm["topology"] == "unstructured"

    opts = esm["options"]
    assert opts["R"] == 6371000.0
    assert opts["dtype"] == "float64"
    assert opts["ghosts"] == 0

    loader = opts["loader"]
    assert loader is not None
    assert loader["reader"] == "mpas_mesh"
    assert loader["check"] == "strict"
    assert loader["path"].startswith("meshes/mpas/")
    assert loader["path"].endswith(".grid.nc")


@pytest.mark.parametrize("name", _FIXTURE_NAMES)
def test_no_inline_geometry(name):
    # Per the 2026-04-20 mayor correction on bead dsc-3nw: fixtures are
    # declarative configs, not serialized geometry blobs.
    esm = _load(_FIXTURE_DIR / f"{name}.esm")
    forbidden = (
        "cells",
        "edges",
        "vertices",
        "lon_cell",
        "lat_cell",
        "area_cell",
        "cells_on_cell",
        "cells_on_edge",
        "edges_on_cell",
    )
    for key in forbidden:
        assert key not in esm
        assert key not in esm["options"]


@pytest.mark.parametrize("name", _FIXTURE_NAMES)
def test_dimensions_obey_icosahedral_euler(name):
    # Quasi-uniform MPAS Voronoi meshes dual to an icosahedral triangulation
    # of the sphere satisfy:
    #   n_edges    = 3 * (n_cells - 2)
    #   n_vertices = 2 * (n_cells - 2)
    # from V - E + F = 2 with F = 2*(n_cells - 2) Delaunay triangles and
    # V = n_cells generators, and max_edges = 6 (12 pentagons, rest hexagons).
    esm = _load(_FIXTURE_DIR / f"{name}.esm")
    dims = esm["dimensions"]
    nc = dims["n_cells"]
    assert dims["n_edges"] == 3 * (nc - 2)
    assert dims["n_vertices"] == 2 * (nc - 2)
    assert dims["max_edges"] == 6


@pytest.mark.parametrize("name", _FIXTURE_NAMES)
def test_mesh_id_matches_filename_and_cell_count(name):
    esm = _load(_FIXTURE_DIR / f"{name}.esm")
    prov = esm["provenance"]
    assert prov["family"] == "mpas"
    assert prov["mesh_id"] == name
    expected_n = int(name.split(".")[1])  # "x1.<n_cells>"
    assert esm["dimensions"]["n_cells"] == expected_n
    assert name in esm["options"]["loader"]["path"]


def test_fixtures_are_json_round_trippable():
    for name in _FIXTURE_NAMES:
        esm = _load(_FIXTURE_DIR / f"{name}.esm")
        # Byte-for-byte determinism of the declarative payload.
        assert json.loads(json.dumps(esm)) == esm


def test_fixtures_cover_a_spread_of_resolutions():
    # Bead acceptance: at least 2-3 canonical resolutions per family.
    cell_counts = sorted(
        _load(_FIXTURE_DIR / f"{n}.esm")["dimensions"]["n_cells"] for n in _FIXTURE_NAMES
    )
    assert len(cell_counts) >= 3
    # Resolutions must not collapse onto a single scale.
    assert cell_counts[-1] >= 4 * cell_counts[0]


# --- Cross-fixture invariants -------------------------------------------------


def test_all_fixtures_share_schema_version_and_radius():
    radii = set()
    schema_versions = set()
    for name in _FIXTURE_NAMES:
        esm = _load(_FIXTURE_DIR / f"{name}.esm")
        radii.add(esm["options"]["R"])
        schema_versions.add(esm["schema_version"])
    assert schema_versions == {"1.0.0"}
    assert radii == {6371000.0}
