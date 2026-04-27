"""Tests for the MPAS unstructured Voronoi grid accessor runtime.

Mirrors ``typescript/tests/mpas.test.ts`` so cross-binding conformance can
compare accessor outputs at pinned query points (GRIDS_API.md §3.5 / §4,
mayor's 2026-04-20 correction on bead dsc-3nw).

Covers the GRIDS_API.md §2.4 / §3.2 / §7 / §10 contract: topology consistency
(neighbor reciprocity, edge-cell refs), metric accuracy (total area =
4 pi R^2 on the synthetic tetra mesh, great-circle edge arcs), loader
contract (§10), and declarative .esm round-trip shape.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from earthsci_discretizations import grids
from earthsci_discretizations.grids import (
    MpasGrid,
    MpasLoader,
    MpasMeshData,
    check_mesh,
    mpas,
    mpas_mesh_data,
)

# --- Synthetic tetrahedral Voronoi mesh --------------------------------------


def _tetra_mesh(R: float = 1.0) -> MpasMeshData:
    """Tiny synthetic global MPAS-like mesh: 4 triangular Voronoi cells at
    tetrahedral vertex directions. Each cell has 3 neighbors (the other
    three cells); total area = 4 pi R^2; per-cell area = pi R^2.

    Matches ``typescript/tests/mpas.test.ts::tetraMesh`` so pinned query
    points agree cross-binding.
    """
    verts = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    n_cells = 4
    max_edges = 3
    lon_cell = np.empty(n_cells)
    lat_cell = np.empty(n_cells)
    x_cell = np.empty(n_cells)
    y_cell = np.empty(n_cells)
    z_cell = np.empty(n_cells)
    for c, (x, y, z) in enumerate(verts):
        n = math.sqrt(x * x + y * y + z * z)
        xh, yh, zh = x / n, y / n, z / n
        lon_cell[c] = math.atan2(yh, xh)
        lat_cell[c] = math.asin(zh)
        x_cell[c] = R * xh
        y_cell[c] = R * yh
        z_cell[c] = R * zh

    area_cell = np.full(n_cells, math.pi * R * R)
    n_edges_on_cell = np.full(n_cells, 3, dtype=np.int32)

    cells_on_cell = np.full((n_cells, max_edges), -1, dtype=np.int32)
    for c in range(n_cells):
        j = 0
        for k in range(n_cells):
            if k == c:
                continue
            cells_on_cell[c, j] = k
            j += 1

    n_edges = 6
    cells_on_edge = np.full((n_edges, 2), -1, dtype=np.int32)
    edge_lookup: dict[tuple[int, int], int] = {}
    e = 0
    for i in range(n_cells):
        for j in range(i + 1, n_cells):
            cells_on_edge[e, 0] = i
            cells_on_edge[e, 1] = j
            edge_lookup[(i, j)] = e
            e += 1

    edges_on_cell = np.full((n_cells, max_edges), -1, dtype=np.int32)
    for c in range(n_cells):
        for j in range(max_edges):
            nb = int(cells_on_cell[c, j])
            key = (c, nb) if c < nb else (nb, c)
            edges_on_cell[c, j] = edge_lookup[key]

    lon_edge = np.empty(n_edges)
    lat_edge = np.empty(n_edges)
    dc_edge = np.empty(n_edges)
    for ei in range(n_edges):
        c1 = int(cells_on_edge[ei, 0])
        c2 = int(cells_on_edge[ei, 1])
        v1 = (x_cell[c1] / R, y_cell[c1] / R, z_cell[c1] / R)
        v2 = (x_cell[c2] / R, y_cell[c2] / R, z_cell[c2] / R)
        d = max(-1.0, min(1.0, v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]))
        dc_edge[ei] = R * math.acos(d)
        mx = 0.5 * (v1[0] + v2[0])
        my = 0.5 * (v1[1] + v2[1])
        mz = 0.5 * (v1[2] + v2[2])
        n = math.sqrt(mx * mx + my * my + mz * mz)
        lon_edge[ei] = math.atan2(my / n, mx / n)
        lat_edge[ei] = math.asin(mz / n)

    dv_edge = dc_edge.copy()

    return mpas_mesh_data(
        lon_cell=lon_cell,
        lat_cell=lat_cell,
        area_cell=area_cell,
        n_edges_on_cell=n_edges_on_cell,
        cells_on_cell=cells_on_cell,
        edges_on_cell=edges_on_cell,
        lon_edge=lon_edge,
        lat_edge=lat_edge,
        cells_on_edge=cells_on_edge,
        dc_edge=dc_edge,
        dv_edge=dv_edge,
        max_edges=max_edges,
        x_cell=x_cell,
        y_cell=y_cell,
        z_cell=z_cell,
        n_vertices=n_cells,
        R=R,
    )


# --- API ----------------------------------------------------------------------


def test_exported_from_grids_namespace():
    assert grids.mpas is mpas
    assert grids.MpasGrid is MpasGrid
    assert grids.MpasLoader is MpasLoader
    assert grids.MpasMeshData is MpasMeshData
    assert grids.mpas_mesh_data is mpas_mesh_data


def test_mesh_is_keyword_only():
    with pytest.raises(TypeError):
        mpas(_tetra_mesh())  # type: ignore[misc] - positional forbidden per §2.4


def test_missing_mesh_and_loader_raises_typeerror():
    with pytest.raises(TypeError):
        mpas()  # type: ignore[call-arg]


def test_mesh_type_checked():
    with pytest.raises(TypeError):
        mpas(mesh={"not": "a mesh"})  # type: ignore[arg-type]


def test_constructs_from_in_memory_mesh():
    g = mpas(mesh=_tetra_mesh(), R=1.0)
    assert g.family == "mpas"
    assert g.topology == "unstructured"
    assert g.dtype == "float64"
    assert g.n_cells == 4
    assert g.n_edges == 6
    assert g.max_edges == 3
    assert g.ghosts == 0
    assert g.loader is None


def test_invalid_dtype_rejected():
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), dtype="float16")


def test_invalid_ghosts_rejected():
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), ghosts=2)
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), ghosts=-1)
    with pytest.raises(TypeError):
        mpas(mesh=_tetra_mesh(), ghosts=1.5)  # type: ignore[arg-type]


def test_invalid_R_rejected():
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), R=0.0)
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), R=-1.0)
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), R=float("nan"))


# --- Topology -----------------------------------------------------------------


def test_neighbor_symmetry_and_edge_refs():
    g = mpas(mesh=_tetra_mesh())
    for c in range(g.n_cells):
        nbs = g.neighbors(c)
        assert len(nbs) == 3
        for nb in nbs:
            assert 0 <= nb < g.n_cells
            assert c in g.neighbors(nb)
    for e in range(g.n_edges):
        c1 = int(g.mesh.cells_on_edge[e, 0])
        c2 = int(g.mesh.cells_on_edge[e, 1])
        assert 0 <= c1 < g.n_cells
        assert 0 <= c2 < g.n_cells
        assert c1 != c2


def test_strict_mode_catches_broken_reciprocity():
    mesh = _tetra_mesh()
    # Corrupt: cell 0 forgets its first neighbor; cell 1 still claims 0.
    mesh.cells_on_cell[0, 0] = -1
    with pytest.raises(AssertionError, match="neighbor symmetry"):
        check_mesh(mesh, True)
    check_mesh(mesh, False)  # lenient mode passes


def test_loader_strict_propagates_to_check_mesh():
    mesh = _tetra_mesh()
    mesh.cells_on_cell[0, 0] = -1
    with pytest.raises(AssertionError, match="neighbor symmetry"):
        mpas(
            mesh=mesh,
            loader=MpasLoader(path="/tmp/fake.nc", check="strict"),
        )
    # Lenient loader skips the reciprocity check.
    g = mpas(mesh=mesh, loader=MpasLoader(path="/tmp/fake.nc", check="lenient"))
    assert g.loader is not None
    assert g.loader.check == "lenient"


def test_strict_mode_catches_out_of_range_neighbor():
    mesh = _tetra_mesh()
    mesh.cells_on_cell[0, 0] = 999
    with pytest.raises(AssertionError):
        check_mesh(mesh, True)


# --- Metric accuracy ----------------------------------------------------------


def test_total_area_matches_4_pi_R2():
    R = 2.5
    g = mpas(mesh=_tetra_mesh(R), R=R)
    expected = 4.0 * math.pi * R * R
    assert abs(g.total_area() - expected) / expected < 1e-12
    acc = sum(g.cell_area(c) for c in range(g.n_cells))
    assert abs(acc - expected) / expected < 1e-12


def test_edge_arc_is_great_circle():
    R = 2.5
    g = mpas(mesh=_tetra_mesh(R), R=R)
    for e in range(g.n_edges):
        c1 = int(g.mesh.cells_on_edge[e, 0])
        c2 = int(g.mesh.cells_on_edge[e, 1])
        x1, y1, z1 = g.cell_center_cart(c1)
        x2, y2, z2 = g.cell_center_cart(c2)
        d = (x1 * x2 + y1 * y2 + z1 * z2) / (R * R)
        d = max(-1.0, min(1.0, d))
        expected = R * math.acos(d)
        assert abs(g.edge_length(e) - expected) / R < 1e-12
        assert abs(g.cell_distance(e) - expected) / R < 1e-12


def test_cell_centers_lie_on_sphere():
    R = 6.371e6
    g = mpas(mesh=_tetra_mesh(R), R=R)
    for c in range(g.n_cells):
        x, y, z = g.cell_center_cart(c)
        r = math.sqrt(x * x + y * y + z * z)
        assert abs(r - R) / R < 1e-12


# --- Pinned accessor outputs (cross-binding reference) -----------------------


def test_pinned_cell0_on_tetra_mesh():
    """Vertex 0 = (1,1,1)/sqrt(3); matches TypeScript test's pinned query point."""
    g = mpas(mesh=_tetra_mesh(1.0), R=1.0)
    lon, lat = g.cell_centers(0)
    assert abs(lon - math.pi / 4.0) < 1e-14
    assert abs(lat - math.asin(1.0 / math.sqrt(3.0))) < 1e-14
    x, y, z = g.cell_center_cart(0)
    assert abs(math.sqrt(x * x + y * y + z * z) - 1.0) < 1e-14
    assert g.metric_eval("lon", 0) == lon
    assert g.metric_eval("area", 0) == g.cell_area(0)
    assert g.metric_eval("dc_edge", 0) == g.cell_distance(0)


def test_metric_eval_cell_and_edge():
    g = mpas(mesh=_tetra_mesh(1.0), R=1.0)
    assert g.metric_eval("x", 0) == g.cell_center_cart(0)[0]
    assert g.metric_eval("y", 0) == g.cell_center_cart(0)[1]
    assert g.metric_eval("z", 0) == g.cell_center_cart(0)[2]
    assert g.metric_eval("n_edges_on_cell", 0) == 3
    assert g.metric_eval("lat_edge", 0) == float(g.mesh.lat_edge[0])
    assert g.metric_eval("lon_edge", 0) == float(g.mesh.lon_edge[0])
    assert g.metric_eval("dv_edge", 0) == float(g.mesh.dv_edge[0])
    with pytest.raises(ValueError):
        g.metric_eval("nonsense", 0)


def test_cell_centers_arrays_when_no_arg():
    g = mpas(mesh=_tetra_mesh())
    lon, lat = g.cell_centers()
    assert isinstance(lon, np.ndarray)
    assert isinstance(lat, np.ndarray)
    assert lon.shape == (g.n_cells,)
    assert lat.shape == (g.n_cells,)


def test_rejects_invalid_indices():
    g = mpas(mesh=_tetra_mesh())
    with pytest.raises(IndexError):
        g.cell_centers(-1)
    with pytest.raises(IndexError):
        g.cell_centers(g.n_cells)
    with pytest.raises(TypeError):
        g.neighbors(3.5)  # type: ignore[arg-type]
    with pytest.raises(IndexError):
        g.cell_area(-1)
    with pytest.raises(IndexError):
        g.edge_length(-1)
    with pytest.raises(IndexError):
        g.edge_length(g.n_edges)
    with pytest.raises(IndexError):
        g.metric_eval("lon_edge", g.n_edges)


# --- Loader (§10) -------------------------------------------------------------


def test_loader_mapping_is_coerced():
    g = mpas(
        mesh=_tetra_mesh(),
        loader={"path": "/tmp/x.nc"},  # defaults: reader=auto, check=strict
    )
    assert g.loader == MpasLoader(path="/tmp/x.nc", reader="auto", check="strict")


def test_loader_rejects_bad_reader_and_check():
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), loader=MpasLoader(path="x", reader="bogus"))
    with pytest.raises(ValueError):
        mpas(mesh=_tetra_mesh(), loader=MpasLoader(path="x", check="bogus"))
    with pytest.raises(TypeError):
        mpas(mesh=_tetra_mesh(), loader=MpasLoader(path=""))
    with pytest.raises(TypeError):
        mpas(mesh=_tetra_mesh(), loader={"reader": "auto"})  # missing path


def test_path_based_loading_without_reader_fn_raises():
    with pytest.raises(TypeError, match="reader_fn"):
        mpas(loader=MpasLoader(path="/tmp/fake.nc", reader="nc4"))


def test_path_based_loading_uses_reader_fn():
    mesh = _tetra_mesh()
    captured = {}

    def reader_fn(path: str) -> MpasMeshData:
        captured["path"] = path
        return mesh

    g = mpas(
        loader=MpasLoader(path="/tmp/fake.nc", reader="mpas_mesh", check="lenient"),
        reader_fn=reader_fn,
    )
    assert captured["path"] == "/tmp/fake.nc"
    assert g.loader == MpasLoader(path="/tmp/fake.nc", reader="mpas_mesh", check="lenient")
    esm = g.to_esm()
    loader_blk = esm["options"]["loader"]
    assert loader_blk == {
        "path": "/tmp/fake.nc",
        "reader": "mpas_mesh",
        "check": "lenient",
    }


def test_reader_fn_must_return_mpas_mesh_data():
    with pytest.raises(TypeError, match="reader_fn"):
        mpas(
            loader=MpasLoader(path="/tmp/fake.nc"),
            reader_fn=lambda _p: {"not": "a mesh"},  # type: ignore[arg-type,return-value]
        )


# --- Declarative .esm lowering ------------------------------------------------


def test_to_esm_shape_matches_ts_binding():
    g = mpas(
        mesh=_tetra_mesh(1.0),
        loader=MpasLoader(path="/tmp/fake.nc", reader="mpas_mesh", check="lenient"),
        R=1.0,
    )
    esm = g.to_esm()
    assert esm["family"] == "mpas"
    assert esm["version"] == "1.0.0"
    assert esm["dtype"] == "float64"
    assert esm["topology"] == "unstructured"
    assert esm["ghosts"] == 0
    assert esm["n_cells"] == 4
    assert esm["n_edges"] == 6
    assert esm["max_edges"] == 3
    assert esm["options"]["R"] == 1.0
    loader_blk = esm["options"]["loader"]
    assert loader_blk == {
        "path": "/tmp/fake.nc",
        "reader": "mpas_mesh",
        "check": "lenient",
    }
    prov = esm["provenance"]
    assert prov["binding"] == "python"
    assert prov["source_sha"] == "dsc-uct"
    assert prov["reader_version"] == "0.1.0"


def test_to_esm_excludes_inline_geometry():
    """2026-04-20 scope correction: .esm is a declarative config, not a
    serialized geometry blob (§6 / mayor bead dsc-3nw correction)."""
    g = mpas(mesh=_tetra_mesh(1.0))
    esm = g.to_esm()
    for forbidden in (
        "cells",
        "edges",
        "vertices",
        "lon_cell",
        "lat_cell",
        "area_cell",
        "cells_on_cell",
        "cells_on_edge",
    ):
        assert forbidden not in esm


def test_to_esm_json_round_trips():
    g = mpas(mesh=_tetra_mesh(1.0))
    esm = g.to_esm()
    assert json.loads(json.dumps(esm)) == esm


def test_to_esm_omits_loader_when_mesh_only():
    g = mpas(mesh=_tetra_mesh(), R=1.0)
    esm = g.to_esm()
    assert esm["options"]["loader"] is None
    assert g.loader is None


# --- Miscellaneous ------------------------------------------------------------


def test_dtype_propagation_float32():
    g = mpas(mesh=_tetra_mesh(), dtype="float32")
    assert g.dtype == "float32"
    assert g.to_esm()["dtype"] == "float32"


def test_determinism_same_mesh_same_lowering():
    m1 = _tetra_mesh(1.0)
    m2 = _tetra_mesh(1.0)
    g1 = mpas(mesh=m1, R=1.0)
    g2 = mpas(mesh=m2, R=1.0)
    # provenance is identical byte-for-byte (no time / hostname leakage).
    assert g1.to_esm() == g2.to_esm()
    for c in range(g1.n_cells):
        assert g1.cell_centers(c) == g2.cell_centers(c)
        assert g1.neighbors(c) == g2.neighbors(c)
        assert g1.cell_area(c) == g2.cell_area(c)


def test_repr_contains_key_dimensions():
    g = mpas(mesh=_tetra_mesh(), R=1.5)
    r = repr(g)
    assert "MpasGrid" in r
    assert "n_cells=4" in r
    assert "n_edges=6" in r
