"""Vertical grid family accessor runtime (1D column).

Conforms to the cross-binding contract in ``docs/GRIDS_API.md`` §2.4, §3.2, §7.

Supported coordinate kinds (rhymes with the Julia sibling, `src/grids/vertical.jl`):

======================  ============================  =================================
Kind                    Value domain                  Required options
======================  ============================  =================================
``"sigma"``             [0, 1]; 1 = surface, 0 = top  ``nz`` (uniform) or ``levels``
``"eta"``               hybrid sigma-pressure (NCAR)  ``ak``, ``bk`` (length ``nz+1``)
``"z"``                 geometric altitude (m)        ``levels`` (strictly increasing)
``"theta"``             potential temperature (K)     ``levels`` (strictly increasing)
``"hybrid_sigma_theta"`` blended sigma→theta          ``levels`` + optional ``transition``
``"z_star"``            generalized height            ``levels`` (strictly increasing)
======================  ============================  =================================

Per the 2026-04-20 scope correction, the ``.esm`` lowering is a small declarative
config (family + coordinate kind + interface levels, plus hybrid coefficients
when applicable), NOT a serialized geometry blob. Centers and widths are derived
on demand from ``levels`` via pure arithmetic.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np

__all__ = ["vertical", "VerticalGrid"]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}

_COORDINATES = frozenset(
    {"sigma", "eta", "z", "theta", "hybrid_sigma_theta", "z_star"}
)

_FAMILY_VERSION = "1.0.0"

_METRIC_NAMES = frozenset({"dz", "z", "sigma", "pressure", "ak", "bk"})


def _coerce_dtype(dtype: Any) -> str:
    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            return "float64"
        if dtype == np.float32:
            return "float32"
        raise ValueError(f"unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(f"dtype must be 'float64' or 'float32', got {dtype!r}")
    return dtype


def _uniform_sigma_levels(nz: int, np_dtype: type) -> np.ndarray:
    if nz < 1:
        raise ValueError(f"vertical: nz must be >= 1; got {nz}")
    # Index 0 = surface (sigma=1), index nz = top (sigma=0).
    return np.asarray(
        [1.0 - k / nz for k in range(nz + 1)], dtype=np_dtype
    )


def _coerce_levels(
    levels: Any,
    np_dtype: type,
    *,
    must_decrease: bool,
    domain: tuple[float, float] | None = None,
    label: str = "levels",
) -> np.ndarray:
    if not isinstance(levels, (list, tuple, np.ndarray)) and not isinstance(
        levels, Iterable
    ):
        raise TypeError(
            f"vertical: `{label}` must be a sequence; got {type(levels).__name__}"
        )
    arr = np.asarray(list(levels), dtype=np_dtype)
    if arr.ndim != 1:
        raise ValueError(f"vertical: `{label}` must be 1-dimensional")
    if arr.size < 2:
        raise ValueError(
            f"vertical: `{label}` must have >= 2 entries; got {arr.size}"
        )
    if domain is not None:
        lo, hi = float(domain[0]), float(domain[1])
        if np.any(arr < lo) or np.any(arr > hi):
            raise ValueError(
                f"vertical: `{label}` entries must lie in [{lo}, {hi}]"
            )
    diffs = np.diff(arr)
    if must_decrease:
        if not np.all(diffs < 0):
            raise ValueError(
                f"vertical: `{label}` must be strictly decreasing"
            )
    else:
        if not np.all(diffs > 0):
            raise ValueError(
                f"vertical: `{label}` must be strictly increasing"
            )
    return arr


def _coerce_hybrid(
    coeffs: Any, np_dtype: type, expected_len: int, label: str
) -> np.ndarray:
    if coeffs is None:
        raise ValueError(f"vertical: `{label}` is required")
    arr = np.asarray(list(coeffs), dtype=np_dtype)
    if arr.ndim != 1:
        raise ValueError(f"vertical: `{label}` must be 1-dimensional")
    if arr.size != expected_len:
        raise ValueError(
            f"vertical: `{label}` must have length nz+1 = {expected_len}; "
            f"got {arr.size}"
        )
    return arr


def _centers_and_widths(levels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = 0.5 * (levels[:-1] + levels[1:])
    widths = np.abs(np.diff(levels))
    return centers, widths


class VerticalGrid:
    """1D vertical column.

    The grid is defined by ``nz + 1`` interface ``levels`` and derives cell
    ``centers`` and ``widths`` arithmetically. The interpretation of the
    interface values depends on ``coordinate``:

    * ``"sigma"`` / ``"hybrid_sigma_theta"``: dimensionless sigma ∈ [0, 1]
      (surface first → top last; ``levels`` is strictly decreasing).
    * ``"z"`` / ``"z_star"``: metres (strictly increasing).
    * ``"theta"``: potential temperature in Kelvin (strictly increasing).
    * ``"eta"``: synthesized sigma = ``ak/p0 + bk`` (strictly decreasing).

    See ``docs/GRIDS_API.md`` §2.4, §3.2, §7.
    """

    family = "vertical"
    topology = "column"

    def __init__(
        self,
        *,
        coordinate: str,
        levels: np.ndarray,
        ak: np.ndarray,
        bk: np.ndarray,
        p0: float,
        dtype: str,
        ghosts: int,
    ) -> None:
        self.coordinate = coordinate
        self.levels = levels
        self.nz = int(levels.size - 1)
        self.centers, self.widths = _centers_and_widths(levels)
        self.ak = ak
        self.bk = bk
        self.p0 = float(p0)
        self.dtype = dtype
        self.ghosts = int(ghosts)

    @property
    def n_cells(self) -> int:
        return self.nz

    @property
    def n_vertices(self) -> int:
        return self.nz + 1

    @property
    def n_edges(self) -> int:
        return self.nz

    @property
    def provenance(self) -> dict:
        import earthsci_toolkit

        return {
            "binding": "python",
            "binding_version": earthsci_toolkit.__version__,
            "family": self.family,
            "version": _FAMILY_VERSION,
            "coordinate": self.coordinate,
            "dtype": self.dtype,
        }

    def cell_centers(self, k: int | None = None):
        """Mid-layer coordinate values (native units).

        With no argument, returns the full centers array of shape ``(nz,)``.
        With an integer ``k``, returns the scalar center at that layer.
        """
        if k is None:
            return self.centers
        self._check_layer(k)
        return float(self.centers[k])

    def cell_widths(self, k: int | None = None):
        """Layer thicknesses (native units, always positive).

        With no argument, returns the full widths array of shape ``(nz,)``.
        With an integer ``k``, returns the scalar width at that layer.
        """
        if k is None:
            return self.widths
        self._check_layer(k)
        return float(self.widths[k])

    def neighbors(self, k: int) -> dict[str, int]:
        """Axis-aligned vertical neighbours of layer ``k``.

        Returns a dict with ``"down"`` → ``k-1`` and/or ``"up"`` → ``k+1``.
        Top and bottom layers drop the out-of-range side.
        """
        self._check_layer(k)
        out: dict[str, int] = {}
        if k > 0:
            out["down"] = k - 1
        if k < self.nz - 1:
            out["up"] = k + 1
        return out

    def metric_eval(self, name: str, k: int) -> float:
        """Evaluate a named metric at layer ``k``.

        Supported names:

        * ``"dz"`` - layer thickness (native units).
        * ``"z"`` - cell-centre value (native units).
        * ``"sigma"`` - sigma at cell centre; valid only for sigma-like
          coordinates (``sigma``, ``hybrid_sigma_theta``, ``eta``).
        * ``"pressure"`` - reference pressure at cell centre
          (``p = ak + bk * p0`` averaged across the layer's two interfaces).
          Requires hybrid coefficients.
        * ``"ak"``, ``"bk"`` - hybrid coefficients averaged across the
          layer's two interfaces. Requires hybrid coefficients.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"vertical: unknown metric name: {name!r}")
        self._check_layer(k)
        if name == "dz":
            return float(self.widths[k])
        if name == "z":
            return float(self.centers[k])
        if name == "sigma":
            if self.coordinate in ("sigma", "hybrid_sigma_theta", "eta"):
                return float(self.centers[k])
            raise ValueError(
                f"vertical: 'sigma' undefined for coordinate {self.coordinate!r}"
            )
        if name == "pressure":
            if self.ak.size == 0 or self.bk.size == 0:
                raise ValueError(
                    "vertical: 'pressure' requires hybrid ak/bk "
                    f"(coordinate {self.coordinate!r} has none)"
                )
            p_lo = float(self.ak[k]) + float(self.bk[k]) * self.p0
            p_hi = float(self.ak[k + 1]) + float(self.bk[k + 1]) * self.p0
            return 0.5 * (p_lo + p_hi)
        if name == "ak":
            if self.ak.size == 0:
                raise ValueError(
                    "vertical: 'ak' unavailable (no hybrid coefficients)"
                )
            return 0.5 * (float(self.ak[k]) + float(self.ak[k + 1]))
        # name == "bk"
        if self.bk.size == 0:
            raise ValueError(
                "vertical: 'bk' unavailable (no hybrid coefficients)"
            )
        return 0.5 * (float(self.bk[k]) + float(self.bk[k + 1]))

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering.

        Returns a §6-schema-shaped config: family + coordinate kind +
        interface levels (+ optional hybrid coefficients) + provenance.
        No derived arrays (centres, widths) appear in the wire form.
        """
        options: dict[str, Any] = {
            "coordinate": self.coordinate,
            "nz": int(self.nz),
            "levels": [float(x) for x in self.levels],
        }
        if self.ak.size > 0:
            options["ak"] = [float(x) for x in self.ak]
        if self.bk.size > 0:
            options["bk"] = [float(x) for x in self.bk]
        if (
            self.coordinate in ("eta", "hybrid_sigma_theta")
            or self.ak.size > 0
            or self.bk.size > 0
        ):
            options["p0"] = float(self.p0)
        return {
            "family": self.family,
            "topology": self.topology,
            "dtype": self.dtype,
            "ndim": 1,
            "ghosts": int(self.ghosts),
            "n_cells": int(self.nz),
            "n_vertices": int(self.nz + 1),
            "n_edges": int(self.nz),
            "options": options,
            "provenance": self.provenance,
            "schema_version": _FAMILY_VERSION,
        }

    def __repr__(self) -> str:
        return (
            f"VerticalGrid(coordinate={self.coordinate!r}, nz={self.nz}, "
            f"dtype={self.dtype!r}, ghosts={self.ghosts})"
        )

    def _check_layer(self, k: int) -> None:
        if not isinstance(k, (int, np.integer)) or isinstance(k, bool):
            raise TypeError(f"vertical: k must be int, got {type(k).__name__}")
        if not 0 <= int(k) < self.nz:
            raise ValueError(
                f"vertical: k out of range [0, {self.nz}): {k!r}"
            )


def vertical(
    *,
    coordinate: str,
    nz: int | None = None,
    levels: Any = None,
    ak: Any = None,
    bk: Any = None,
    p0: float = 1.0e5,
    transition: float | None = None,
    dtype: Any = "float64",
    ghosts: int = 0,
) -> VerticalGrid:
    """Generate a 1D vertical grid.

    See ``docs/GRIDS_API.md`` §2.4 for the cross-binding signature contract.
    All options are keyword-only; ``coordinate`` is required.

    Parameters
    ----------
    coordinate : str
        One of ``"sigma"``, ``"eta"``, ``"z"``, ``"theta"``,
        ``"hybrid_sigma_theta"``, ``"z_star"``.
    nz : int, optional
        Number of layers. Required for uniform ``"sigma"`` (or
        ``"hybrid_sigma_theta"``) when ``levels`` is omitted. Must be >= 1.
    levels : sequence of float, optional
        Explicit interface values (length ``nz + 1``). Required for ``"z"``,
        ``"theta"``, ``"z_star"``. Optional for ``"sigma"`` /
        ``"hybrid_sigma_theta"`` (overrides ``nz`` if both given).
    ak, bk : sequence of float, optional
        Hybrid A/B coefficients at interfaces (length ``nz + 1``). Required
        for ``"eta"``. Optional for ``"hybrid_sigma_theta"``.
    p0 : float, optional
        Reference surface pressure (Pa). Default 1.0e5.
    transition : float, optional
        ``hybrid_sigma_theta``-only: sigma value in ``(0, 1)`` at which the
        coordinate transitions from pure-sigma to pure-theta.
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision. Default ``"float64"``.
    ghosts : int, optional
        Halo layer width. Default 0.
    """
    if coordinate is None:
        raise TypeError("vertical: required option `coordinate` is missing")
    if not isinstance(coordinate, str):
        raise TypeError(
            f"vertical: `coordinate` must be str, got {type(coordinate).__name__}"
        )
    if coordinate not in _COORDINATES:
        raise ValueError(
            f"vertical: unknown coordinate {coordinate!r}; expected one of "
            f"{sorted(_COORDINATES)}"
        )

    if isinstance(ghosts, bool) or not isinstance(ghosts, (int, np.integer)):
        raise TypeError(f"vertical: ghosts must be int, got {type(ghosts).__name__}")
    ghosts = int(ghosts)
    if ghosts < 0:
        raise ValueError(f"vertical: ghosts must be >= 0, got {ghosts}")

    dtype_str = _coerce_dtype(dtype)
    np_dtype = _DTYPE_MAP[dtype_str]

    if nz is not None:
        if isinstance(nz, bool) or not isinstance(nz, (int, np.integer)):
            raise TypeError(f"vertical: nz must be int, got {type(nz).__name__}")
        nz = int(nz)
        if nz < 1:
            raise ValueError(f"vertical: nz must be >= 1; got {nz}")

    p0 = float(p0)
    if not (p0 > 0 and math.isfinite(p0)):
        raise ValueError(f"vertical: p0 must be positive and finite, got {p0}")

    if coordinate == "sigma":
        if levels is not None:
            lv = _coerce_levels(
                levels, np_dtype, must_decrease=True, domain=(0.0, 1.0)
            )
            if nz is not None and lv.size != nz + 1:
                raise ValueError(
                    f"vertical: sigma: nz={nz} inconsistent with levels length {lv.size}"
                )
        else:
            if nz is None:
                raise TypeError(
                    "vertical: sigma requires `nz` or `levels`"
                )
            lv = _uniform_sigma_levels(nz, np_dtype)
        ak_arr = np.zeros(0, dtype=np_dtype)
        bk_arr = np.zeros(0, dtype=np_dtype)

    elif coordinate in ("z", "theta", "z_star"):
        if levels is None:
            raise TypeError(
                f"vertical: {coordinate!r} requires explicit `levels`"
            )
        lv = _coerce_levels(levels, np_dtype, must_decrease=False)
        if nz is not None and lv.size != nz + 1:
            raise ValueError(
                f"vertical: {coordinate}: nz={nz} inconsistent with levels length {lv.size}"
            )
        ak_arr = np.zeros(0, dtype=np_dtype)
        bk_arr = np.zeros(0, dtype=np_dtype)

    elif coordinate == "eta":
        if ak is None:
            raise TypeError("vertical: 'eta' requires `ak` (length nz+1)")
        if bk is None:
            raise TypeError("vertical: 'eta' requires `bk` (length nz+1)")
        ak_probe = np.asarray(list(ak), dtype=np_dtype)
        bk_probe = np.asarray(list(bk), dtype=np_dtype)
        if ak_probe.size != bk_probe.size:
            raise ValueError(
                "vertical: 'eta' ak/bk must have equal length; "
                f"got {ak_probe.size} vs {bk_probe.size}"
            )
        nz_eff = nz if nz is not None else (ak_probe.size - 1)
        if nz_eff < 1:
            raise ValueError(f"vertical: nz must be >= 1; got {nz_eff}")
        ak_arr = _coerce_hybrid(ak, np_dtype, nz_eff + 1, "ak")
        bk_arr = _coerce_hybrid(bk, np_dtype, nz_eff + 1, "bk")
        sigma = ak_arr / np_dtype(p0) + bk_arr
        if not np.all(np.diff(sigma) < 0):
            raise ValueError(
                "vertical: 'eta' synthesized sigma (ak/p0 + bk) must be "
                "strictly decreasing"
            )
        lv = sigma.astype(np_dtype, copy=False)

    elif coordinate == "hybrid_sigma_theta":
        if levels is None and nz is None:
            raise TypeError(
                "vertical: 'hybrid_sigma_theta' requires `nz` or `levels`"
            )
        if levels is not None:
            lv = _coerce_levels(
                levels, np_dtype, must_decrease=True, domain=(0.0, 1.0)
            )
            if nz is not None and lv.size != nz + 1:
                raise ValueError(
                    f"vertical: hybrid_sigma_theta: nz={nz} inconsistent with "
                    f"levels length {lv.size}"
                )
        else:
            lv = _uniform_sigma_levels(nz, np_dtype)
        if transition is not None:
            if not (0.0 < float(transition) < 1.0):
                raise ValueError(
                    f"vertical: hybrid_sigma_theta `transition` must be in (0, 1); "
                    f"got {transition}"
                )
        ak_arr = (
            _coerce_hybrid(ak, np_dtype, lv.size, "ak")
            if ak is not None
            else np.zeros(0, dtype=np_dtype)
        )
        bk_arr = (
            _coerce_hybrid(bk, np_dtype, lv.size, "bk")
            if bk is not None
            else np.zeros(0, dtype=np_dtype)
        )

    else:  # pragma: no cover - guarded above
        raise AssertionError(f"unreachable coordinate: {coordinate!r}")

    return VerticalGrid(
        coordinate=coordinate,
        levels=lv,
        ak=ak_arr,
        bk=bk_arr,
        p0=p0,
        dtype=dtype_str,
        ghosts=ghosts,
    )
