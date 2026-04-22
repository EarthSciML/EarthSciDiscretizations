"""Tests for the DUO icosahedral grid accessor runtime.

Covers the GRIDS_API.md §2.4 / §3.2 / §7 / §10 contract: topology consistency
(Euler characteristic, subdivision counts, neighbor reciprocity), metric
accuracy (total area = 4 pi R^2, sphere-radius invariance), loader contract
(§10: builtin path scheme, DuoLoader / mapping acceptance), and declarative
.esm round-trip shape.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import earthsci_toolkit
from earthsci_toolkit import grids
from earthsci_toolkit.grids import DuoGrid, DuoLoader, duo


def _builtin(level: int) -> DuoLoader:
    return DuoLoader(
        path=f"builtin://icosahedral/{level}",
        reader="builtin_icosahedral",
        check="strict",
    )


# ------------------------------------------------------------------ API --


def test_exported_from_grids_namespace():
    assert grids.duo is duo
    assert grids.DuoGrid is DuoGrid
    assert grids.DuoLoader is DuoLoader


def test_loader_is_keyword_only():
    with pytest.raises(TypeError):
        duo(_builtin(0))  # type: ignore[misc] - positional forbidden per §2.4


def test_missing_loader_raises_typeerror():
    with pytest.raises(TypeError):
        duo()  # type: ignore[call-arg]


def test_default_R_and_dtype():
    g = duo(loader=_builtin(0))
    assert g.R == 6.371e6
    assert g.dtype == "float64"
    assert g.ghosts == 0


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        duo(loader=_builtin(0), dtype="float16")


def test_invalid_ghosts_rejected():
    with pytest.raises(ValueError):
        duo(loader=_builtin(0), ghosts=-1)
    with pytest.raises(TypeError):
        duo(loader=_builtin(0), ghosts=1.5)  # type: ignore[arg-type]


def test_invalid_R_rejected():
    with pytest.raises(ValueError):
        duo(loader=_builtin(0), R=0.0)
    with pytest.raises(ValueError):
        duo(loader=_builtin(0), R=-1.0)


# ------------------------------------------------------------ topology --


def test_level0_is_bare_icosahedron():
    g = duo(loader=_builtin(0), R=1.0)
    assert g.n_cells == 20
    assert g.n_vertices == 12
    assert g.n_edges == 30
    # Euler characteristic of a sphere: V - E + F = 2
    assert g.n_vertices - g.n_edges + g.n_cells == 2
    assert g.family == "duo"
    assert g.topology == "unstructured"


@pytest.mark.parametrize("r", [0, 1, 2, 3])
def test_subdivision_counts(r):
    g = duo(loader=_builtin(r), R=1.0)
    assert g.n_cells == 20 * 4**r
    assert g.n_vertices == 10 * 4**r + 2
    assert g.n_edges == 30 * 4**r
    assert g.n_vertices - g.n_edges + g.n_cells == 2


def test_cell_neighbors_symmetric_and_valid():
    g = duo(loader=_builtin(2), R=1.0)
    Nc = g.n_cells
    # Closed mesh: every cell has exactly 3 neighbors (no sentinels).
    assert np.all(g.cell_neighbors >= 0)
    assert np.all(g.cell_neighbors < Nc)
    # Reciprocity: if c -> n via some edge, then n -> c via some edge.
    for c in range(Nc):
        for k in range(3):
            nb = int(g.cell_neighbors[k, c])
            back = (
                int(g.cell_neighbors[0, nb]),
                int(g.cell_neighbors[1, nb]),
                int(g.cell_neighbors[2, nb]),
            )
            assert c in back


def test_neighbors_accessor_matches_cell_neighbors():
    g = duo(loader=_builtin(1), R=1.0)
    for c in range(g.n_cells):
        assert g.neighbors(c) == (
            int(g.cell_neighbors[0, c]),
            int(g.cell_neighbors[1, c]),
            int(g.cell_neighbors[2, c]),
        )


def test_vertex_faces_cover_all_cells():
    g = duo(loader=_builtin(1), R=1.0)
    seen: set[int] = set()
    for faces_at_v in g.vertex_faces:
        seen.update(faces_at_v)
    assert seen == set(range(g.n_cells))


# --------------------------------------------------------------- area --


@pytest.mark.parametrize("r", [0, 1, 2, 3])
def test_total_area_unit_sphere(r):
    g = duo(loader=_builtin(r), R=1.0)
    assert abs(g.total_area() - 4 * math.pi) < 1e-10 * (4 * math.pi)


def test_total_area_earth_radius():
    R = 6.371e6
    g = duo(loader=_builtin(2), R=R)
    expected = 4 * math.pi * R * R
    assert abs(g.total_area() - expected) / expected < 1e-10


def test_all_areas_positive():
    g = duo(loader=_builtin(2), R=1.0)
    assert np.all(g.area > 0)


# ------------------------------------------------------------- centers --


def test_cell_centers_shape_and_range():
    g = duo(loader=_builtin(1), R=1.0)
    lon, lat = g.cell_centers()
    assert lon.shape == (g.n_cells,)
    assert lat.shape == (g.n_cells,)
    assert np.all(lon >= -math.pi - 1e-12)
    assert np.all(lon <= math.pi + 1e-12)
    assert np.all(lat >= -math.pi / 2 - 1e-12)
    assert np.all(lat <= math.pi / 2 + 1e-12)


def test_cell_center_scalar_matches_bulk():
    g = duo(loader=_builtin(1), R=1.0)
    lon, lat = g.cell_centers()
    for c in (0, 5, g.n_cells - 1):
        s_lon, s_lat = g.cell_centers(c)
        assert math.isclose(s_lon, float(lon[c]), abs_tol=1e-12)
        assert math.isclose(s_lat, float(lat[c]), abs_tol=1e-12)


def test_cell_centers_lie_on_sphere_of_radius_R():
    R = 6.371e6
    g = duo(loader=_builtin(2), R=R)
    radii = np.sqrt(
        g.cell_cart[0, :] ** 2 + g.cell_cart[1, :] ** 2 + g.cell_cart[2, :] ** 2
    )
    assert np.allclose(radii, R, rtol=1e-12)


# ------------------------------------------------------------ metrics --


def test_metric_eval_delivers_area_lon_lat_xyz():
    R = 6.371e6
    g = duo(loader=_builtin(1), R=R)
    for c in (0, 5, g.n_cells - 1):
        assert g.metric_eval("area", c) == float(g.area[c])
        assert g.metric_eval("lon", c) == float(g.lon[c])
        assert g.metric_eval("lat", c) == float(g.lat[c])
        x = g.metric_eval("x", c)
        y = g.metric_eval("y", c)
        z = g.metric_eval("z", c)
        assert math.isclose(math.sqrt(x * x + y * y + z * z), R, rel_tol=1e-12)


def test_metric_eval_rejects_unknown_name():
    g = duo(loader=_builtin(0), R=1.0)
    with pytest.raises(ValueError):
        g.metric_eval("nonsense", 0)


def test_metric_eval_rejects_out_of_range_cell():
    g = duo(loader=_builtin(0), R=1.0)
    with pytest.raises(ValueError):
        g.metric_eval("area", g.n_cells)
    with pytest.raises(ValueError):
        g.metric_eval("area", -1)


def test_neighbors_rejects_out_of_range_cell():
    g = duo(loader=_builtin(0), R=1.0)
    with pytest.raises(ValueError):
        g.neighbors(g.n_cells)


# ---------------------------------------------------------- loader ----


def test_loader_accepts_duoloader_struct_and_mapping():
    ldr_struct = DuoLoader(path="builtin://icosahedral/1", reader="builtin_icosahedral")
    ldr_map = {
        "path": "builtin://icosahedral/1",
        "reader": "builtin_icosahedral",
        "check": "strict",
    }
    g1 = duo(loader=ldr_struct, R=1.0)
    g2 = duo(loader=ldr_map, R=1.0)
    assert g1.n_cells == g2.n_cells
    np.testing.assert_array_equal(g1.lon, g2.lon)
    np.testing.assert_array_equal(g1.lat, g2.lat)


def test_loader_missing_path_rejected():
    with pytest.raises(TypeError):
        duo(loader={"reader": "builtin_icosahedral"})


def test_loader_duo_mesh_reader_not_yet_implemented():
    with pytest.raises(ValueError):
        duo(loader={"path": "/tmp/does_not_exist.duo", "reader": "duo_mesh"})


def test_loader_unknown_scheme_rejected():
    with pytest.raises(ValueError):
        duo(loader={"path": "bogus://x", "reader": "unknown_reader"})


def test_loader_negative_level_rejected():
    with pytest.raises(ValueError):
        duo(loader={"path": "builtin://icosahedral/-1", "reader": "builtin_icosahedral"})


def test_loader_unparseable_level_rejected():
    with pytest.raises(ValueError):
        duo(
            loader={
                "path": "builtin://icosahedral/notanint",
                "reader": "builtin_icosahedral",
            }
        )


def test_loader_type_rejected():
    with pytest.raises(TypeError):
        duo(loader=42)  # type: ignore[arg-type]


# --------------------------------------------------------------- to_esm --


def test_to_esm_is_declarative():
    g = duo(loader=_builtin(2), R=6.371e6)
    doc = g.to_esm()
    assert doc["family"] == "duo"
    assert doc["topology"] == "unstructured"
    assert doc["dtype"] == "float64"
    assert doc["n_cells"] == g.n_cells
    assert doc["n_vertices"] == g.n_vertices
    assert doc["n_edges"] == g.n_edges
    assert doc["options"]["level"] == 2
    assert doc["options"]["R"] == 6.371e6
    assert doc["options"]["loader"]["path"] == "builtin://icosahedral/2"
    assert doc["options"]["loader"]["reader"] == "builtin_icosahedral"
    assert doc["schema_version"] == "1.0.0"
    assert "provenance" in doc

    # Mayor's correction: no inline geometry arrays in .esm.
    for forbidden in ("cells", "vertices_xyz", "edges_list", "faces"):
        assert forbidden not in doc

    import json

    wire = json.dumps(doc)
    # Declarative config must stay small - no per-cell payload.
    assert len(wire) < 2_000


def test_to_esm_json_roundtrip():
    import json

    g = duo(loader=_builtin(1), R=1.0)
    doc = g.to_esm()
    reparsed = json.loads(json.dumps(doc))
    assert reparsed["options"]["level"] == 1
    assert reparsed["options"]["R"] == 1.0


def test_provenance_carries_binding_identity():
    g = duo(loader=_builtin(0), R=1.0)
    prov = g.to_esm()["provenance"]
    assert prov["binding"] == "python"
    assert prov["binding_version"] == earthsci_toolkit.__version__
    assert prov["family"] == "duo"
    assert prov["level"] == 0
    assert prov["reader"] == "builtin_icosahedral"


# -------------------------------------------------------------- dtypes --


def test_float32_dtype_propagates_to_arrays():
    g = duo(loader=_builtin(1), R=1.0, dtype="float32")
    assert g.dtype == "float32"
    assert g.lon.dtype == np.float32
    assert g.lat.dtype == np.float32
    assert g.area.dtype == np.float32
    assert g.cell_cart.dtype == np.float32
    assert g.vertices.dtype == np.float32
    assert g.to_esm()["dtype"] == "float32"
    # Float32 total-area tolerance loosens but must still be sane.
    assert abs(g.total_area() - 4 * math.pi) < 1e-4 * (4 * math.pi)


def test_numpy_dtype_accepted():
    g = duo(loader=_builtin(0), R=1.0, dtype=np.dtype("float64"))
    assert g.dtype == "float64"


# ---------------------------------------------------------- determinism --


def test_two_calls_produce_identical_grids():
    g1 = duo(loader=_builtin(2), R=1.0)
    g2 = duo(loader=_builtin(2), R=1.0)
    np.testing.assert_array_equal(g1.lon, g2.lon)
    np.testing.assert_array_equal(g1.lat, g2.lat)
    np.testing.assert_array_equal(g1.area, g2.area)
    np.testing.assert_array_equal(g1.faces, g2.faces)
    np.testing.assert_array_equal(g1.cell_neighbors, g2.cell_neighbors)
