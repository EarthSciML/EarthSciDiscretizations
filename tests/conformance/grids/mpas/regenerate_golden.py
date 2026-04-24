#!/usr/bin/env python3
"""Regenerate fixtures.json and golden accessor outputs for the mpas
cross-language conformance harness.

Two reasons Python is the current reference binding instead of Julia (per
docs/GRIDS_API.md §4.3):

1. The Julia mpas runtime is being landed in bead dsc-7j0; at the time this
   harness was authored src/grids/mpas.jl did not yet exist. When dsc-7j0
   lands, a follow-up should port this regenerator to Julia and re-emit
   golden/*.json from the Julia binding (the tolerance already covers the
   libm-cross-binding drift; regenerating should not change byte sums).
2. The mpas accessor runtime does not use transcendental-heavy math on
   pinned query points — cell_centers / neighbors / area / edge metrics
   read from stored arrays. Cross-binding drift is therefore driven almost
   entirely by float JSON round-trip, not libm families, so Python-vs-Julia
   reference choice is effectively neutral here.

Usage:
    python3 tests/conformance/grids/mpas/regenerate_golden.py [--fixtures]

Without --fixtures, only golden/*.json is regenerated (fixtures.json stays
committed). With --fixtures, fixtures.json is also regenerated from the
tetra-mesh builder below — use this only after a deliberate fixture change.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PY_SRC = REPO / "python" / "src"

sys.path.insert(0, str(PY_SRC))

import numpy as np  # noqa: E402

from earthsci_toolkit.grids.mpas import (  # noqa: E402
    MpasMeshData,
    mpas,
    mpas_mesh_data,
)


def _tetra_mesh_arrays(R: float) -> dict:
    """Build the synthetic tetrahedral Voronoi mesh as plain lists.

    Matches ``python/tests/test_mpas.py::_tetra_mesh`` and
    ``typescript/tests/mpas.test.ts::tetraMesh`` so pinned query outputs
    agree cross-binding.
    """
    verts = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    n_cells = 4
    max_edges = 3
    lon_cell = [0.0] * n_cells
    lat_cell = [0.0] * n_cells
    x_cell = [0.0] * n_cells
    y_cell = [0.0] * n_cells
    z_cell = [0.0] * n_cells
    for c, (x, y, z) in enumerate(verts):
        n = math.sqrt(x * x + y * y + z * z)
        xh, yh, zh = x / n, y / n, z / n
        lon_cell[c] = math.atan2(yh, xh)
        lat_cell[c] = math.asin(zh)
        x_cell[c] = R * xh
        y_cell[c] = R * yh
        z_cell[c] = R * zh

    area_cell = [math.pi * R * R] * n_cells
    n_edges_on_cell = [3] * n_cells

    cells_on_cell = [[-1] * max_edges for _ in range(n_cells)]
    for c in range(n_cells):
        j = 0
        for k in range(n_cells):
            if k == c:
                continue
            cells_on_cell[c][j] = k
            j += 1

    n_edges = 6
    cells_on_edge = [[-1, -1] for _ in range(n_edges)]
    edge_lookup: dict[tuple[int, int], int] = {}
    e = 0
    for i in range(n_cells):
        for j in range(i + 1, n_cells):
            cells_on_edge[e] = [i, j]
            edge_lookup[(i, j)] = e
            edge_lookup[(j, i)] = e
            e += 1

    edges_on_cell = [[-1] * max_edges for _ in range(n_cells)]
    for c in range(n_cells):
        for j in range(max_edges):
            nb = cells_on_cell[c][j]
            key = (c, nb) if c < nb else (nb, c)
            edges_on_cell[c][j] = edge_lookup[key]

    lon_edge = [0.0] * n_edges
    lat_edge = [0.0] * n_edges
    dc_edge = [0.0] * n_edges
    for ei in range(n_edges):
        c1, c2 = cells_on_edge[ei]
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

    dv_edge = list(dc_edge)

    return {
        "max_edges": max_edges,
        "n_vertices": n_cells,
        "lon_cell": lon_cell,
        "lat_cell": lat_cell,
        "x_cell": x_cell,
        "y_cell": y_cell,
        "z_cell": z_cell,
        "area_cell": area_cell,
        "n_edges_on_cell": n_edges_on_cell,
        "cells_on_cell": _flatten(cells_on_cell),
        "edges_on_cell": _flatten(edges_on_cell),
        "lon_edge": lon_edge,
        "lat_edge": lat_edge,
        "cells_on_edge": _flatten(cells_on_edge),
        "dc_edge": dc_edge,
        "dv_edge": dv_edge,
    }


def _flatten(rows: list[list[int]]) -> list[int]:
    out: list[int] = []
    for r in rows:
        out.extend(int(v) for v in r)
    return out


def _build_fixtures_spec() -> dict:
    return {
        "version": "1.0.0",
        "family": "mpas",
        "tolerance": {
            "relative": 1e-12,
            "notes": (
                "Per docs/GRIDS_API.md §4.2. 1e-12 is tighter than the libm "
                "drift observed cross-binding on these accessors (no heavy "
                "transcendentals on the hot path — mesh arrays are stored) "
                "and loose enough to absorb JSON float round-trip."
            ),
        },
        "reference_binding": "python",
        "reference_binding_notes": (
            "Interim reference until Julia mpas (bead dsc-7j0) lands. Per "
            "docs/GRIDS_API.md §4.3, Julia is the canonical reference; when "
            "the Julia mpas runtime is available a follow-up should "
            "regenerate golden/*.json via a Julia port of this script."
        ),
        "fixtures": [
            {
                "name": "tetra_small",
                "opts": {"R": 1.0, "ghosts": 0, "dtype": "float64"},
                "mesh": _tetra_mesh_arrays(1.0),
                "query_cells": [0, 1, 2, 3],
                "query_edges": [0, 1, 2, 3, 4, 5],
            },
            {
                "name": "tetra_realistic",
                "opts": {"R": 6371000.0, "ghosts": 0, "dtype": "float64"},
                "mesh": _tetra_mesh_arrays(6371000.0),
                "query_cells": [0, 1, 2, 3],
                "query_edges": [0, 1, 2, 3, 4, 5],
            },
        ],
    }


def _mesh_from_fixture(fixture: dict) -> MpasMeshData:
    m = fixture["mesh"]
    R = float(fixture["opts"]["R"])
    return mpas_mesh_data(
        lon_cell=np.asarray(m["lon_cell"], dtype=np.float64),
        lat_cell=np.asarray(m["lat_cell"], dtype=np.float64),
        area_cell=np.asarray(m["area_cell"], dtype=np.float64),
        n_edges_on_cell=np.asarray(m["n_edges_on_cell"], dtype=np.int32),
        cells_on_cell=np.asarray(m["cells_on_cell"], dtype=np.int32),
        edges_on_cell=np.asarray(m["edges_on_cell"], dtype=np.int32),
        lon_edge=np.asarray(m["lon_edge"], dtype=np.float64),
        lat_edge=np.asarray(m["lat_edge"], dtype=np.float64),
        cells_on_edge=np.asarray(m["cells_on_edge"], dtype=np.int32),
        dc_edge=np.asarray(m["dc_edge"], dtype=np.float64),
        dv_edge=np.asarray(m["dv_edge"], dtype=np.float64),
        max_edges=int(m["max_edges"]),
        x_cell=np.asarray(m["x_cell"], dtype=np.float64),
        y_cell=np.asarray(m["y_cell"], dtype=np.float64),
        z_cell=np.asarray(m["z_cell"], dtype=np.float64),
        n_vertices=int(m["n_vertices"]),
        R=R,
    )


_CELL_METRICS = ("lon", "lat", "area", "x", "y", "z", "n_edges_on_cell")
_EDGE_METRICS = ("lon_edge", "lat_edge", "dc_edge", "dv_edge")


def _build_golden(fixture: dict) -> dict:
    opts = fixture["opts"]
    mesh = _mesh_from_fixture(fixture)
    grid = mpas(
        mesh=mesh,
        R=float(opts["R"]),
        dtype=str(opts["dtype"]),
        ghosts=int(opts["ghosts"]),
    )
    cells = [int(c) for c in fixture["query_cells"]]
    edges = [int(e) for e in fixture["query_edges"]]

    centers = []
    carts = []
    neighbors_list = []
    areas = []
    cell_metrics = {name: [] for name in _CELL_METRICS}
    for c in cells:
        lon, lat = grid.cell_centers(c)
        centers.append({"lon": float(lon), "lat": float(lat)})
        x, y, z = grid.cell_center_cart(c)
        carts.append({"x": float(x), "y": float(y), "z": float(z)})
        neighbors_list.append([int(n) for n in grid.neighbors(c)])
        areas.append(float(grid.cell_area(c)))
        for name in _CELL_METRICS:
            cell_metrics[name].append(float(grid.metric_eval(name, c)))

    edge_lengths = []
    cell_distances = []
    edge_metrics = {name: [] for name in _EDGE_METRICS}
    for e in edges:
        edge_lengths.append(float(grid.edge_length(e)))
        cell_distances.append(float(grid.cell_distance(e)))
        for name in _EDGE_METRICS:
            edge_metrics[name].append(float(grid.metric_eval(name, e)))

    return {
        "fixture": fixture["name"],
        "family": "mpas",
        "opts": opts,
        "n_cells": int(grid.n_cells),
        "n_edges": int(grid.n_edges),
        "max_edges": int(grid.max_edges),
        "query_cells": cells,
        "query_edges": edges,
        "cell_centers": centers,
        "cell_centers_cart": carts,
        "neighbors": neighbors_list,
        "cell_area": areas,
        "cell_metrics": cell_metrics,
        "edge_length": edge_lengths,
        "cell_distance": cell_distances,
        "edge_metrics": edge_metrics,
    }


def _dump_json(obj: object, path: Path) -> None:
    text = json.dumps(obj, indent=2, ensure_ascii=False)
    path.write_text(text + "\n")
    print(f"wrote {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--fixtures",
        action="store_true",
        help="Also regenerate fixtures.json (not just golden/*.json).",
    )
    args = ap.parse_args()

    fixtures_path = HERE / "fixtures.json"
    golden_dir = HERE / "golden"
    golden_dir.mkdir(parents=True, exist_ok=True)

    if args.fixtures or not fixtures_path.exists():
        spec = _build_fixtures_spec()
        _dump_json(spec, fixtures_path)
    else:
        spec = json.loads(fixtures_path.read_text())

    for fx in spec["fixtures"]:
        out = _build_golden(fx)
        _dump_json(out, golden_dir / f"{fx['name']}.json")


if __name__ == "__main__":
    main()
