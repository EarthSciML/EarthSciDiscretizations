#!/usr/bin/env python3
"""Regenerate golden accessor outputs for the latlon conformance harness.

The Julia latlon binding is not yet on main (tracked by bead dsc-1ts); per
`docs/GRIDS_API.md` §4.3 Julia is the nominal reference binding for ULP ties,
but for the initial corpus the Python binding is used. Once the Julia binding
lands, this script should be retired (or kept as a cross-check) and the golden
regenerated from Julia.

The script activates an editable install of the Python binding so that it runs
against whatever is on the current branch. Run from the repo root:

    python tests/conformance/grids/latlon/regenerate_golden.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]
PYTHON_SRC = REPO_ROOT / "python" / "src"
sys.path.insert(0, str(PYTHON_SRC))

from earthsci_toolkit import grids  # noqa: E402

FIXTURES_PATH = HERE / "fixtures.json"
GOLDEN_DIR = HERE / "golden"

DIR_KEYS = ("W", "E", "S", "N")
METRIC_NAMES = (
    "J",
    "g_lonlon",
    "g_latlat",
    "g_lonlat",
    "ginv_lonlon",
    "ginv_latlat",
    "ginv_lonlat",
)


def _build_grid(opts: dict):
    variant = opts["variant"]
    kwargs = {
        "variant": variant,
        "R": float(opts["R"]),
        "ghosts": int(opts["ghosts"]),
        "dtype": opts["dtype"],
        "pole_policy": opts["pole_policy"],
    }
    if variant == "regular":
        kwargs["nlon"] = int(opts["nlon"])
        kwargs["nlat"] = int(opts["nlat"])
    else:
        kwargs["nlon_per_row"] = [int(x) for x in opts["nlon_per_row"]]
        if "lat_edges" in opts:
            kwargs["lat_edges"] = [float(x) for x in opts["lat_edges"]]
    return grids.lat_lon(**kwargs)


def _neighbor_entry(nbr):
    if nbr is None:
        return None
    j, i = nbr
    return [int(j), int(i)]


def build_output(fixture: dict) -> dict:
    opts = fixture["opts"]
    grid = _build_grid(opts)
    qps = fixture["query_points"]

    centers = []
    neighbors = {k: [] for k in DIR_KEYS}
    metrics = {name: [] for name in METRIC_NAMES}
    areas = []

    for qp in qps:
        j, i = int(qp[0]), int(qp[1])
        lon, lat = grid.cell_centers(j=j, i=i)
        centers.append({"lon": float(lon), "lat": float(lat)})

        nbrs = grid.neighbors(j, i)
        for dkey in DIR_KEYS:
            neighbors[dkey].append(_neighbor_entry(nbrs[dkey]))

        for name in METRIC_NAMES:
            metrics[name].append(float(grid.metric_eval(name, j, i)))
        areas.append(float(grid.metric_eval("area", j, i)))

    return {
        "fixture": fixture["name"],
        "family": "lat_lon",
        "variant": opts["variant"],
        "opts": opts,
        "n_cells": int(grid.n_cells),
        "query_points": qps,
        "cell_centers": centers,
        "neighbors_W": neighbors["W"],
        "neighbors_E": neighbors["E"],
        "neighbors_S": neighbors["S"],
        "neighbors_N": neighbors["N"],
        "metric_J": metrics["J"],
        "metric_g_lonlon": metrics["g_lonlon"],
        "metric_g_latlat": metrics["g_latlat"],
        "metric_g_lonlat": metrics["g_lonlat"],
        "metric_ginv_lonlon": metrics["ginv_lonlon"],
        "metric_ginv_latlat": metrics["ginv_latlat"],
        "metric_ginv_lonlat": metrics["ginv_lonlat"],
        "area": areas,
    }


def main() -> None:
    spec = json.loads(FIXTURES_PATH.read_text())
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    for fx in spec["fixtures"]:
        out = build_output(fx)
        path = GOLDEN_DIR / f"{fx['name']}.json"
        with path.open("w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
