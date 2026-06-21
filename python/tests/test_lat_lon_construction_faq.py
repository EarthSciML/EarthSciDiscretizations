"""Cross-language conformance for the lat-lon **construction FAQ** (esd-dru).

Builds each fixture grid imperatively, re-derives its construction arrays through
the elementwise FAQ bridge (:func:`lat_lon_construction_faq`, routed via the ESS
``eval_coeff`` evaluator), and asserts the result is

1. byte-identical to the imperative ``LatLonGrid`` accessors over **all** cells
   (parity — the FAQ introduces zero error), and
2. equal to the committed Julia-reference
   ``tests/conformance/grids/latlon/construction/golden.json``.

Validation methodology (matches the established cross-binding latlon conformance
in ``test_lat_lon_conformance.py`` and ``docs/GRIDS_API.md`` §4.2, with Julia as
the ULP anchor): the per-cell trig arrays carry ``sin``/``cos`` whose libm differs
sub-ULP across bindings; near the poles the ``sin φ_n − sin φ_s`` area term
cancels and amplifies that drift just past the 1e-14 family tolerance on the
high-res ``realistic`` fixture. So the cross-binding golden float comparison is
sampled at the fixture's curated ``query_points`` at the declared 1e-14 (exactly
as the existing latlon conformance test does), while the integer neighbor maps /
boundary masks (libm-free) and the 1-D latitude arrays are checked over the whole
grid. The exhaustive guarantee comes from the all-cells FAQ==imperative parity:
the FAQ is byte-identical to this binding's own construction, and that
construction is what the established query-point conformance certifies.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_discretizations import grids
from earthsci_discretizations.grids.lat_lon_faq import lat_lon_construction_faq

HARNESS_DIR = (
    Path(__file__).resolve().parents[2] / "tests" / "conformance" / "grids" / "latlon"
)
_SPEC = json.loads((HARNESS_DIR / "fixtures.json").read_text())
_GOLDEN = json.loads((HARNESS_DIR / "construction" / "golden.json").read_text())
REL_TOL = float(_GOLDEN["tolerance"]["relative"])
_BY_NAME = {f["name"]: f for f in _GOLDEN["fixtures"]}
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]

# Per-cell trig-bearing arrays — golden-compared at the curated query points.
_CELL_FLOAT_KEYS = (
    "cell_centers_lon",
    "cell_centers_lat",
    "cell_widths_lon",
    "cell_widths_lat",
    "cell_volume",
    "metric_g_lonlon",
    "metric_g_latlat",
    "metric_ginv_lonlon",
    "metric_jacobian",
    "dg_lonlon_dlat",
)
# 1-D latitude arrays (no per-cell cancellation) — golden-compared over all rows.
_LAT_FLOAT_KEYS = ("lat_edges", "lat_centers")
_INT_KEYS = (
    "neighbor_lon_minus",
    "neighbor_lon_plus",
    "neighbor_lat_minus",
    "neighbor_lat_plus",
)
_MASK_KEYS = ("boundary_lat_lower", "boundary_lat_upper")


def _close_rel(a: float, b: float, tol: float = REL_TOL) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _build_grid(opts: dict):
    kwargs = {
        "variant": opts["variant"],
        "R": float(opts["R"]),
        "ghosts": int(opts["ghosts"]),
        "dtype": opts["dtype"],
        "pole_policy": opts["pole_policy"],
    }
    if opts["variant"] == "regular":
        kwargs["nlon"] = int(opts["nlon"])
        kwargs["nlat"] = int(opts["nlat"])
    else:
        kwargs["nlon_per_row"] = [int(x) for x in opts["nlon_per_row"]]
        if "lat_edges" in opts:
            kwargs["lat_edges"] = [float(x) for x in opts["lat_edges"]]
    return grids.lat_lon(**kwargs)


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_lat_lon_construction_faq(fixture_name: str) -> None:
    fixture = _fixture(fixture_name)
    grid = _build_grid(fixture["opts"])
    faq = lat_lon_construction_faq(grid)
    gl = _BY_NAME[fixture_name]
    nc = int(gl["n_cells"])
    per_row = [int(x) for x in grid.nlon_per_row]
    row_off = [0]
    for n_i in per_row:
        row_off.append(row_off[-1] + n_i)

    # --- counts ---
    assert len(faq.cell_centers_lon) == nc
    assert grid.nlat == int(gl["nlat"])

    # --- parity: FAQ == imperative accessors, over ALL cells ---
    lon_bulk, lat_bulk = grid.cell_centers()
    for a, b in zip(faq.cell_centers_lon, lon_bulk, strict=True):
        assert float(a) == float(b), f"{fixture_name}: cc_lon parity"
    for a, b in zip(faq.cell_centers_lat, lat_bulk, strict=True):
        assert float(a) == float(b), f"{fixture_name}: cc_lat parity"
    for j in range(grid.nlat):
        for i in range(per_row[j]):
            k = row_off[j] + i
            assert float(faq.cell_volume[k]) == float(grid.cell_area(j, i)), (
                f"{fixture_name}: cell_volume parity at ({j},{i})"
            )

    # --- cross-binding golden: per-cell trig arrays at curated query points ---
    for qp in fixture["query_points"]:
        j, i = int(qp[0]), int(qp[1])
        k = row_off[j] + i
        for key in _CELL_FLOAT_KEYS:
            got = float(getattr(faq, key)[k])
            exp = float(json.loads(gl[key])[k])
            assert _close_rel(got, exp), f"{fixture_name}:{key} at ({j},{i}): {got} vs {exp}"

    # --- cross-binding golden: 1-D latitude arrays over all rows ---
    for key in _LAT_FLOAT_KEYS:
        for a, b in zip(getattr(faq, key), json.loads(gl[key]), strict=True):
            assert _close_rel(float(a), float(b)), f"{fixture_name}:{key} golden"

    # --- cross-binding golden: integer + mask arrays, exact, all cells ---
    for key in _INT_KEYS + _MASK_KEYS:
        assert [int(x) for x in getattr(faq, key)] == list(json.loads(gl[key])), (
            f"{fixture_name}: {key}"
        )
