"""Arakawa staggering accessor runtime.

Conforms to the cross-binding contract in ``docs/GRIDS_API.md`` §2.4, §3.2, §7.

Arakawa staggering is a *transform* over an underlying base grid. Given a base
grid (Cartesian or lat-lon) plus a stagger label (A/B/C/D/E), the runtime
provides on-demand accessors for cell centers, u-face, v-face, and corner
locations per the staggering convention.

Per the 2026-04-20 scope correction, the ``.esm`` lowering is a small
declarative config (family + base-grid ref + stagger + dimensions + extents).
Geometry is derived from that config by the accessors here, not serialized as
an inline blob.

Staggering conventions (2-D, horizontal):

* A: h, u, v colocated at cell centers.
* B: h at cell centers; u, v colocated at corners.
* C: h at cell centers; u at u-faces (east-west); v at v-faces (north-south).
* D: h at cell centers; u at v-faces; v at u-faces (swapped from C).
* E: rotated B-grid; topologically equivalent to B with a 45° rotation flag
  in the lowered config.

Phase 1 (``cartesian``) and Phase 3 (``lat_lon``) base-grid modules have not
landed on ``main`` yet; a minimal :class:`CartesianBase` is provided here so
the accessors are exercisable today. Future per-family base grids supply the
same primitives via the :class:`BaseGrid` protocol and plug in without
touching this module.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Protocol, runtime_checkable

import numpy as np

__all__ = [
    "arakawa",
    "ArakawaGrid",
    "BaseGrid",
    "CartesianBase",
]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}

_STAGGERS = frozenset({"A", "B", "C", "D", "E"})

# Canonical variable-location table (§2.4 + Julia/Rust parity).
_VARIABLE_LOCATIONS: Mapping[str, tuple[str, str, str]] = {
    "A": ("cell_center", "cell_center", "cell_center"),
    "B": ("cell_center", "corner", "corner"),
    "C": ("cell_center", "u_edge", "v_edge"),
    "D": ("cell_center", "v_edge", "u_edge"),
    # E is topologically B; the 45° rotation is carried separately.
    "E": ("cell_center", "corner", "corner"),
}

_LOCATIONS = frozenset({"cell_center", "u_edge", "v_edge", "corner"})

_METRIC_NAMES = frozenset({"dx", "dy", "area"})


def _location_shape(loc: str, nx: int, ny: int) -> tuple[int, int]:
    """Shape of ``loc`` on a base grid of ``(nx, ny)`` interior cells."""
    if loc == "cell_center":
        return (nx, ny)
    if loc == "u_edge":
        return (nx + 1, ny)
    if loc == "v_edge":
        return (nx, ny + 1)
    if loc == "corner":
        return (nx + 1, ny + 1)
    raise ValueError(f"arakawa: unknown location: {loc!r}")


@runtime_checkable
class BaseGrid(Protocol):
    """Primitives every Arakawa base grid must expose.

    Phase 1 (``cartesian``) and Phase 3 (``lat_lon``) grid types will subclass
    or satisfy this Protocol. Until they land, :class:`CartesianBase` provides
    a minimal in-module implementation.
    """

    nx: int
    ny: int

    def cell_center(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of cell center ``(i, j)`` (0-based)."""

    def x_edge(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of the u-face location at ``(i, j)``."""

    def y_edge(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of the v-face location at ``(i, j)``."""

    def corner(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of the cell corner at ``(i, j)``."""

    def dx(self) -> float:
        """Uniform cell spacing in x."""

    def dy(self) -> float:
        """Uniform cell spacing in y."""

    def to_esm(self) -> dict:
        """Declarative ``.esm``-shaped summary of the base grid."""


class CartesianBase:
    """Minimal Cartesian base grid: uniform rectangular mesh over
    ``[xlo, xhi] × [ylo, yhi]`` with ``(nx, ny)`` interior cells.

    Matches the Julia and Rust ``CartesianBase`` helpers used as a temporary
    stand-in until the Phase 1 Cartesian grid module lands in Python.
    """

    def __init__(
        self,
        *,
        xlo: float,
        xhi: float,
        ylo: float,
        yhi: float,
        nx: int,
        ny: int,
    ) -> None:
        if isinstance(nx, bool) or not isinstance(nx, (int, np.integer)):
            raise TypeError(f"CartesianBase: nx must be int, got {type(nx).__name__}")
        if isinstance(ny, bool) or not isinstance(ny, (int, np.integer)):
            raise TypeError(f"CartesianBase: ny must be int, got {type(ny).__name__}")
        nx = int(nx)
        ny = int(ny)
        if nx < 1:
            raise ValueError(f"CartesianBase: nx must be >= 1, got {nx}")
        if ny < 1:
            raise ValueError(f"CartesianBase: ny must be >= 1, got {ny}")
        xlo = float(xlo)
        xhi = float(xhi)
        ylo = float(ylo)
        yhi = float(yhi)
        for label, val in (("xlo", xlo), ("xhi", xhi), ("ylo", ylo), ("yhi", yhi)):
            if not math.isfinite(val):
                raise ValueError(f"CartesianBase: {label} must be finite, got {val!r}")
        if xhi <= xlo:
            raise ValueError(
                f"CartesianBase: xhi must be > xlo, got xlo={xlo}, xhi={xhi}"
            )
        if yhi <= ylo:
            raise ValueError(
                f"CartesianBase: yhi must be > ylo, got ylo={ylo}, yhi={yhi}"
            )
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.nx = nx
        self.ny = ny

    def dx(self) -> float:
        return (self.xhi - self.xlo) / float(self.nx)

    def dy(self) -> float:
        return (self.yhi - self.ylo) / float(self.ny)

    def cell_center(self, i: int, j: int) -> tuple[float, float]:
        dx = self.dx()
        dy = self.dy()
        return (self.xlo + (i + 0.5) * dx, self.ylo + (j + 0.5) * dy)

    def x_edge(self, i: int, j: int) -> tuple[float, float]:
        dx = self.dx()
        dy = self.dy()
        return (self.xlo + i * dx, self.ylo + (j + 0.5) * dy)

    def y_edge(self, i: int, j: int) -> tuple[float, float]:
        dx = self.dx()
        dy = self.dy()
        return (self.xlo + (i + 0.5) * dx, self.ylo + j * dy)

    def corner(self, i: int, j: int) -> tuple[float, float]:
        dx = self.dx()
        dy = self.dy()
        return (self.xlo + i * dx, self.ylo + j * dy)

    def to_esm(self) -> dict:
        return {
            "family": "cartesian",
            "nx": int(self.nx),
            "ny": int(self.ny),
            "extent": [[self.xlo, self.ylo], [self.xhi, self.yhi]],
        }

    def __repr__(self) -> str:
        return (
            f"CartesianBase(xlo={self.xlo}, xhi={self.xhi}, "
            f"ylo={self.ylo}, yhi={self.yhi}, nx={self.nx}, ny={self.ny})"
        )


class ArakawaGrid:
    """Staggered grid: a :class:`BaseGrid` plus a stagger label.

    See module docstring. Indexing is 0-based, consistent with the Python
    binding's other families.

    Attributes
    ----------
    base : BaseGrid
        Underlying base grid supplying geometric primitives.
    stagger : str
        One of ``"A"``, ``"B"``, ``"C"``, ``"D"``, ``"E"``.
    dtype : str
        ``"float64"`` or ``"float32"``.
    ghosts : int
        Halo cell width.
    """

    family = "arakawa"
    topology = "block_structured"

    def __init__(
        self,
        *,
        base: BaseGrid,
        stagger: str,
        dtype: str,
        ghosts: int,
    ) -> None:
        self.base = base
        self.stagger = stagger
        self.dtype = dtype
        self.ghosts = ghosts

    @property
    def nx(self) -> int:
        return int(self.base.nx)

    @property
    def ny(self) -> int:
        return int(self.base.ny)

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny

    @property
    def rotated(self) -> bool:
        return self.stagger == "E"

    def variable_location(self, var: str) -> str:
        """Location where named variable ``var`` (``"h"``, ``"u"``, ``"v"``)
        lives under this grid's stagger.
        """
        h, u, v = _VARIABLE_LOCATIONS[self.stagger]
        if var == "h":
            return h
        if var == "u":
            return u
        if var == "v":
            return v
        raise ValueError(f"arakawa: unknown variable {var!r}; expected 'h', 'u', or 'v'")

    def location_shape(self, loc: str) -> tuple[int, int]:
        """Shape ``(ni, nj)`` of ``loc`` on this grid."""
        if loc not in _LOCATIONS:
            raise ValueError(f"arakawa: unknown location {loc!r}")
        return _location_shape(loc, self.nx, self.ny)

    def variable_shape(self, var: str) -> tuple[int, int]:
        """Shape of the location that ``var`` lives on under this grid's stagger."""
        return self.location_shape(self.variable_location(var))

    def cell_centers(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of cell center ``(i, j)`` (0-based)."""
        self._check_bounds("cell_center", i, j)
        return self.base.cell_center(int(i), int(j))

    def u_face(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of the u-variable location at ``(i, j)`` under this stagger."""
        loc = self.variable_location("u")
        self._check_bounds(loc, i, j)
        return self._coord_at(loc, int(i), int(j))

    def v_face(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of the v-variable location at ``(i, j)`` under this stagger."""
        loc = self.variable_location("v")
        self._check_bounds(loc, i, j)
        return self._coord_at(loc, int(i), int(j))

    def corners(self, i: int, j: int) -> tuple[float, float]:
        """``(x, y)`` of cell corner ``(i, j)`` (0-based; ``0..nx`` × ``0..ny``)."""
        self._check_bounds("corner", i, j)
        return self.base.corner(int(i), int(j))

    def coord(self, loc: str, i: int, j: int) -> tuple[float, float]:
        """Generic coordinate accessor for any location ``loc``."""
        if loc not in _LOCATIONS:
            raise ValueError(f"arakawa: unknown location {loc!r}")
        self._check_bounds(loc, i, j)
        return self._coord_at(loc, int(i), int(j))

    def neighbors(
        self, loc: str, i: int, j: int
    ) -> dict[str, tuple[int, int] | None]:
        """Axial neighbors of ``(loc, i, j)`` as ``{dir: (i', j') | None}``.

        Neighbor is ``None`` at a domain boundary. All neighbors stay on the
        same location (no cross-location stepping — that's the caller's job
        per stagger semantics).
        """
        if loc not in _LOCATIONS:
            raise ValueError(f"arakawa: unknown location {loc!r}")
        self._check_bounds(loc, i, j)
        ni, nj = self.location_shape(loc)
        ii = int(i)
        jj = int(j)
        return {
            "W": (ii - 1, jj) if ii > 0 else None,
            "E": (ii + 1, jj) if ii + 1 < ni else None,
            "S": (ii, jj - 1) if jj > 0 else None,
            "N": (ii, jj + 1) if jj + 1 < nj else None,
        }

    def cell_neighbors(self, i: int, j: int) -> dict[str, tuple[int, int] | None]:
        """Shortcut for ``neighbors("cell_center", i, j)``."""
        return self.neighbors("cell_center", i, j)

    def metric_eval(self, name: str, i: int, j: int) -> float:
        """Evaluate a named metric at cell ``(i, j)``.

        Valid ``name`` values:

        * ``"dx"`` — zonal spacing
        * ``"dy"`` — meridional spacing
        * ``"area"`` — cell area / volume

        For :class:`CartesianBase`, metrics are uniform and ``(i, j)`` are
        only used to validate bounds.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"arakawa: unknown metric {name!r}")
        self._check_bounds("cell_center", i, j)
        dx = float(self.base.dx())
        dy = float(self.base.dy())
        if name == "dx":
            return dx
        if name == "dy":
            return dy
        return dx * dy  # area

    @property
    def provenance(self) -> dict:
        import earthsci_toolkit

        return {
            "binding": "python",
            "binding_version": earthsci_toolkit.__version__,
            "source": "earthsci_toolkit.grids.arakawa",
            "stagger": self.stagger,
        }

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering per the 2026-04-20 scope correction.

        Returns a §6-schema-shaped config: family, stagger, base reference,
        dtype, ghosts, n_cells, topology. No inline geometry arrays.
        """
        return {
            "family": self.family,
            "version": "1.0.0",
            "dtype": self.dtype,
            "topology": self.topology,
            "ghosts": int(self.ghosts),
            "n_cells": int(self.n_cells),
            "stagger": self.stagger,
            "rotated": self.rotated,
            "base": self.base.to_esm(),
            "provenance": self.provenance,
        }

    def __repr__(self) -> str:
        return (
            f"ArakawaGrid(stagger={self.stagger!r}, nx={self.nx}, ny={self.ny}, "
            f"dtype={self.dtype!r}, ghosts={self.ghosts})"
        )

    # internals -----------------------------------------------------------

    def _coord_at(self, loc: str, i: int, j: int) -> tuple[float, float]:
        if loc == "cell_center":
            return self.base.cell_center(i, j)
        if loc == "u_edge":
            return self.base.x_edge(i, j)
        if loc == "v_edge":
            return self.base.y_edge(i, j)
        return self.base.corner(i, j)

    def _check_bounds(self, loc: str, i, j) -> None:
        if i is None or j is None:
            raise TypeError("arakawa: indices (i, j) required")
        if isinstance(i, bool) or not isinstance(i, (int, np.integer)):
            raise TypeError(f"arakawa: i must be int, got {type(i).__name__}")
        if isinstance(j, bool) or not isinstance(j, (int, np.integer)):
            raise TypeError(f"arakawa: j must be int, got {type(j).__name__}")
        ni, nj = self.location_shape(loc)
        ii = int(i)
        jj = int(j)
        if not 0 <= ii < ni:
            raise ValueError(f"arakawa: i out of range [0, {ni}) for {loc!r}: {i!r}")
        if not 0 <= jj < nj:
            raise ValueError(f"arakawa: j out of range [0, {nj}) for {loc!r}: {j!r}")


def arakawa(
    *,
    base: BaseGrid,
    stagger: str,
    dtype: str = "float64",
    ghosts: int = 0,
) -> ArakawaGrid:
    """Generate an Arakawa-staggered grid as a transform over a base grid.

    See ``docs/GRIDS_API.md`` §2.4 for the cross-binding signature contract.
    All options are keyword-only; ``base`` and ``stagger`` are required.

    Parameters
    ----------
    base : BaseGrid
        Underlying base grid providing cell_center / x_edge / y_edge / corner
        primitives. :class:`CartesianBase` supplies a minimal implementation
        until Phase 1 / Phase 3 base grids land in Python.
    stagger : {"A", "B", "C", "D", "E"}
        Arakawa staggering label.
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision declared in the ``.esm`` lowering. Default
        ``"float64"``.
    ghosts : int, optional
        Halo cell width. Default ``0``.
    """
    if base is None:
        raise TypeError("arakawa: base is required")
    if not all(
        hasattr(base, attr)
        for attr in ("nx", "ny", "cell_center", "x_edge", "y_edge", "corner", "dx", "dy", "to_esm")
    ):
        raise TypeError(
            "arakawa: base must satisfy the BaseGrid protocol "
            "(nx, ny, cell_center, x_edge, y_edge, corner, dx, dy, to_esm)"
        )

    if not isinstance(stagger, str):
        raise TypeError(
            f"arakawa: stagger must be str, got {type(stagger).__name__}"
        )
    if stagger not in _STAGGERS:
        raise ValueError(
            f"arakawa: stagger must be one of {sorted(_STAGGERS)}, got {stagger!r}"
        )

    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            dtype = "float64"
        elif dtype == np.float32:
            dtype = "float32"
        else:
            raise ValueError(f"arakawa: unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(
            f"arakawa: dtype must be 'float64' or 'float32', got {dtype!r}"
        )

    if isinstance(ghosts, bool) or not isinstance(ghosts, (int, np.integer)):
        raise TypeError(f"arakawa: ghosts must be int, got {type(ghosts).__name__}")
    ghosts = int(ghosts)
    if ghosts < 0:
        raise ValueError(f"arakawa: ghosts must be >= 0, got {ghosts}")

    return ArakawaGrid(base=base, stagger=stagger, dtype=dtype, ghosts=ghosts)
