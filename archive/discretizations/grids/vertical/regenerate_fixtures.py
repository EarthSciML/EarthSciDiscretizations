"""Regenerate canonical .esm fixtures for the vertical grid family.

The `.esm` files in this directory are declarative configs produced by calling
`earthsci_toolkit.grids.vertical(...).to_esm()` for a small set of pinned
resolutions across the six coordinate kinds. Fixtures are written with
`json.dumps(..., sort_keys=True, indent=2)` so diffs stay minimal.

Usage
-----
From the repo root::

    PYTHONPATH=python/src python3 discretizations/grids/vertical/regenerate_fixtures.py

The script is idempotent: running it on a clean binding regenerates the same
bytes. When the reference (Python) binding's `to_esm` changes, rerun this and
commit the diff; the Julia inline tests in `test/test_vertical_fixtures.jl`
will verify the resulting `.esm` files stay schema-valid.
"""

from __future__ import annotations

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


def _load_vertical():
    try:
        from earthsci_toolkit.grids import vertical
    except ModuleNotFoundError:
        repo = HERE.parents[2]
        sys.path.insert(0, str(repo / "python" / "src"))
        from earthsci_toolkit.grids import vertical
    return vertical


def build_fixtures(vertical):
    # CAM-like 12-level hybrid eta ak/bk profile, increasing in altitude
    # (interface 0 = surface, interface nz = top). Values chosen to be a
    # plausible representative eta profile (ak in Pa, bk dimensionless) and
    # to satisfy the strict-decrease requirement on synthesized sigma.
    eta_ak = [
        0.0,
        219.4067,
        489.5209,
        988.2500,
        1805.2010,
        2983.7240,
        4462.3340,
        6160.5880,
        7851.2430,
        7731.2710,
        5652.6050,
        2609.2010,
        0.0,
    ]
    eta_bk = [
        1.0,
        0.9851122,
        0.9634060,
        0.9242153,
        0.8558366,
        0.7447181,
        0.5830062,
        0.3762173,
        0.1803898,
        0.0691445,
        0.0164329,
        0.0019929,
        0.0,
    ]
    return {
        # Uniform sigma at two canonical resolutions. sigma is the coarsest
        # coordinate family: pure geometric subdivision of [0, 1].
        "sigma_uniform_n16": vertical(coordinate="sigma", nz=16),
        "sigma_uniform_n64": vertical(coordinate="sigma", nz=64),
        # Geometric altitude profile from surface to ~32 km, 32 layers of
        # increasing thickness. Models a troposphere + stratosphere column
        # without assuming a hybrid coordinate.
        "z_troposphere_l32": vertical(
            coordinate="z",
            levels=[
                0.0,
                100.0,
                200.0,
                300.0,
                500.0,
                750.0,
                1000.0,
                1500.0,
                2000.0,
                2500.0,
                3000.0,
                3750.0,
                4500.0,
                5500.0,
                6500.0,
                7500.0,
                8500.0,
                9500.0,
                10500.0,
                11500.0,
                12500.0,
                14000.0,
                15500.0,
                17000.0,
                18500.0,
                20000.0,
                22000.0,
                24000.0,
                26000.0,
                28000.0,
                30000.0,
                32000.0,
                34000.0,
            ],
        ),
        # Hybrid eta coordinate, 12 layers. ak/bk are a representative
        # CAM-style profile; synthesized sigma must be strictly decreasing.
        "eta_hybrid_l12": vertical(
            coordinate="eta",
            ak=eta_ak,
            bk=eta_bk,
            p0=1.0e5,
        ),
        # Isentropic potential-temperature surfaces from 280 K to 380 K.
        "theta_isentropic_l10": vertical(
            coordinate="theta",
            levels=[280.0, 290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0,
                    360.0, 370.0, 380.0],
        ),
    }


def main() -> None:
    vertical = _load_vertical()
    fixtures = build_fixtures(vertical)
    for name, grid in fixtures.items():
        doc = grid.to_esm()
        path = HERE / f"{name}.esm"
        path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n")
        print(f"wrote {path.relative_to(HERE.parents[2])}")


if __name__ == "__main__":
    main()
