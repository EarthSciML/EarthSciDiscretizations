"""Cubed-sphere grid accessor runtime (FV3 gnomonic convention).

Conforms to the cross-binding contract in ``docs/GRIDS_API.md`` §2.4, §3.2, §7.

Per the 2026-04-20 scope correction, the ``.esm`` lowering is a small
declarative config (family, dimensions, extents, generator reference), not a
serialized geometry blob. Geometry is derived on demand by the accessors:

* :meth:`CubedSphereGrid.cell_centers` - cell-center (lon, lat)
* :meth:`CubedSphereGrid.neighbors` - 4-way face connectivity, including
  cross-panel edges
* :meth:`CubedSphereGrid.metric_eval` - metric-tensor components at cell
  centers (covariant, inverse, and derived quantities like face area)

Cross-binding conformance (``docs/GRIDS_API.md`` §4) compares these accessor
outputs at pinned query points rather than serialized bytes.
"""

from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

__all__ = ["cubed_sphere", "CubedSphereGrid"]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}

_DIRECTIONS = ("W", "E", "S", "N")

# Standard cubed-sphere panel connectivity (0-indexed).
# Mirrors src/grids/panel_connectivity.jl (Julia's 1..6 shifted to Python's 0..5).
# Entry: (neighbor_panel, neighbor_edge, reverse_along_edge)
_PANEL_CONNECTIVITY: Mapping[int, Mapping[str, tuple[int, str, bool]]] = {
    0: {"W": (4, "E", False), "E": (1, "W", False),
        "S": (5, "N", False), "N": (2, "S", False)},
    1: {"W": (0, "E", False), "E": (3, "W", False),
        "S": (5, "E", True),  "N": (2, "E", False)},
    2: {"W": (4, "N", True),  "E": (1, "N", False),
        "S": (0, "N", False), "N": (3, "N", True)},
    3: {"W": (1, "E", False), "E": (4, "W", False),
        "S": (5, "S", True),  "N": (2, "N", True)},
    4: {"W": (3, "E", False), "E": (0, "W", False),
        "S": (5, "W", False), "N": (2, "W", True)},
    5: {"W": (4, "S", False), "E": (1, "S", True),
        "S": (3, "S", True),  "N": (0, "S", False)},
}


def _gnomonic_to_cart(xi: float, eta: float, panel: int) -> tuple[float, float, float]:
    """Unit-sphere Cartesian (x, y, z) for computational coords (xi, eta) on ``panel``."""
    X = math.tan(xi)
    Y = math.tan(eta)
    r = math.sqrt(1.0 + X * X + Y * Y)
    if panel == 0:
        return (1.0 / r, X / r, Y / r)
    if panel == 1:
        return (-X / r, 1.0 / r, Y / r)
    if panel == 2:
        return (-Y / r, X / r, 1.0 / r)
    if panel == 3:
        return (-1.0 / r, -X / r, Y / r)
    if panel == 4:
        return (X / r, -1.0 / r, Y / r)
    if panel == 5:
        return (Y / r, X / r, -1.0 / r)
    raise ValueError(f"invalid panel: {panel!r}")


def _gnomonic_to_lonlat(xi: float, eta: float, panel: int) -> tuple[float, float]:
    x, y, z = _gnomonic_to_cart(xi, eta, panel)
    lon = math.atan2(y, x)
    lat = math.asin(max(-1.0, min(1.0, z)))
    return lon, lat


def _gnomonic_metric(xi: float, eta: float, R: float) -> tuple[float, float, float, float]:
    """Covariant metric ``(J, g_xixi, g_etaeta, g_xieta)`` at (xi, eta) for radius R.

    The gnomonic metric is panel-independent: each panel is a rigid rotation of
    the sphere, so the (xi, eta)->sphere metric is the same on all six panels.
    """
    X = math.tan(xi)
    Y = math.tan(eta)
    D2 = 1.0 + X * X + Y * Y
    sx = 1.0 + X * X
    sy = 1.0 + Y * Y
    D3 = D2 * math.sqrt(D2)
    D4 = D2 * D2
    R2 = R * R
    J = R2 * sx * sy / D3
    gxx = R2 * sx * sx * sy / D4
    gyy = R2 * sy * sy * sx / D4
    gxy = -R2 * X * Y * sx * sy / D4
    return J, gxx, gyy, gxy


def _spherical_quad_area(corners: list[tuple[float, float, float]], R: float) -> float:
    """Area of a spherical quadrilateral given 4 unit-sphere corner vectors and radius.

    Sums interior angles via the planar-projection approach in
    ``src/grids/metric_tensors.jl::compute_cell_area`` - stable for the small
    non-planar quads at cube corners.
    """
    n = 4
    total = 0.0
    for k in range(n):
        v_prev = np.asarray(corners[(k - 1) % n], dtype=np.float64)
        v_curr = np.asarray(corners[k], dtype=np.float64)
        v_next = np.asarray(corners[(k + 1) % n], dtype=np.float64)
        t1 = np.cross(v_curr, np.cross(v_prev - v_curr, v_curr))
        t2 = np.cross(v_curr, np.cross(v_next - v_curr, v_curr))
        n1 = float(np.linalg.norm(t1))
        n2 = float(np.linalg.norm(t2))
        if n1 < 1e-15 or n2 < 1e-15:
            total += math.pi / 2
        else:
            cos_ang = float(np.dot(t1, t2)) / (n1 * n2)
            total += math.acos(max(-1.0, min(1.0, cos_ang)))
    return R * R * (total - (n - 2) * math.pi)


_METRIC_NAMES = frozenset({
    "J", "g_xixi", "g_etaeta", "g_xieta",
    "ginv_xixi", "ginv_etaeta", "ginv_xieta",
    "area",
})


class CubedSphereGrid:
    """Gnomonic cubed-sphere grid (6 panels x Nc x Nc cells).

    See ``docs/GRIDS_API.md`` §2.4, §3.2, §7. The grid stores only the
    declarative parameters plus pre-tabulated 1D axis arrays. All geometric
    quantities (lon, lat, metric, area) are derived on demand from the
    gnomonic_c6 generator rather than materialized eagerly; bulk properties
    cache arrays on first access.
    """

    family = "cubed_sphere"
    topology = "block_structured"

    def __init__(
        self,
        *,
        Nc: int,
        R: float,
        dtype: str,
        ghosts: int,
        xi_edges: np.ndarray,
        eta_edges: np.ndarray,
        xi_centers: np.ndarray,
        eta_centers: np.ndarray,
    ) -> None:
        self.Nc = Nc
        self.R = R
        self.dtype = dtype
        self.ghosts = ghosts
        self.xi_edges = xi_edges
        self.eta_edges = eta_edges
        self.xi_centers = xi_centers
        self.eta_centers = eta_centers
        self._lon_cache: np.ndarray | None = None
        self._lat_cache: np.ndarray | None = None
        self._area_cache: np.ndarray | None = None

    @property
    def n_cells(self) -> int:
        return 6 * self.Nc * self.Nc

    @property
    def provenance(self) -> dict:
        import earthsci_discretizations

        return {
            "binding": "python",
            "binding_version": earthsci_discretizations.__version__,
            "source": "earthsci_discretizations.grids.cubed_sphere",
            "generator": "gnomonic_c6",
        }

    def cell_centers(
        self, p: int | None = None, i: int | None = None, j: int | None = None
    ) -> tuple[float, float] | tuple[np.ndarray, np.ndarray]:
        """Cell-center geographic coords.

        With no arguments, returns ``(lon, lat)`` arrays of shape ``(6, Nc, Nc)``.
        With ``(p, i, j)``, returns the scalar ``(lon, lat)`` pair for that cell.
        """
        if p is None and i is None and j is None:
            return self.lon, self.lat
        self._check_cell(p, i, j)
        return _gnomonic_to_lonlat(
            float(self.xi_centers[i]), float(self.eta_centers[j]), int(p)
        )

    def neighbors(self, p: int, i: int, j: int) -> dict[str, tuple[int, int, int]]:
        """Face neighbors of cell ``(p, i, j)`` as ``{dir: (p', i', j')}``.

        Interior cells resolve within the same panel; boundary cells cross to
        the neighbor panel via :data:`_PANEL_CONNECTIVITY`, with the along-edge
        index reversed when the edges are oriented opposite.
        """
        self._check_cell(p, i, j)
        Nc = self.Nc
        out: dict[str, tuple[int, int, int]] = {}
        out["W"] = (p, i - 1, j) if i > 0 else self._cross_panel_neighbor(p, "W", i, j)
        out["E"] = (
            (p, i + 1, j) if i < Nc - 1 else self._cross_panel_neighbor(p, "E", i, j)
        )
        out["S"] = (p, i, j - 1) if j > 0 else self._cross_panel_neighbor(p, "S", i, j)
        out["N"] = (
            (p, i, j + 1) if j < Nc - 1 else self._cross_panel_neighbor(p, "N", i, j)
        )
        return out

    def metric_eval(self, name: str, p: int, i: int, j: int) -> float:
        """Evaluate a metric field at cell center ``(p, i, j)``.

        Valid ``name`` values:

        * ``"J"`` - Jacobian determinant
        * ``"g_xixi"``, ``"g_etaeta"``, ``"g_xieta"`` - covariant metric
        * ``"ginv_xixi"``, ``"ginv_etaeta"``, ``"ginv_xieta"`` - inverse metric
        * ``"area"`` - cell face area (uses panel-specific corner geometry)

        The covariant/inverse metric is panel-independent; ``area`` varies
        only at cube-corner cells where the quadrilateral is non-planar.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"unknown metric name: {name!r}")
        self._check_cell(p, i, j)
        if name == "area":
            return self._cell_area(int(p), int(i), int(j))
        xi = float(self.xi_centers[i])
        eta = float(self.eta_centers[j])
        J, gxx, gyy, gxy = _gnomonic_metric(xi, eta, self.R)
        if name == "J":
            return J
        if name == "g_xixi":
            return gxx
        if name == "g_etaeta":
            return gyy
        if name == "g_xieta":
            return gxy
        det = gxx * gyy - gxy * gxy
        if name == "ginv_xixi":
            return gyy / det
        if name == "ginv_etaeta":
            return gxx / det
        return -gxy / det  # ginv_xieta

    @property
    def lon(self) -> np.ndarray:
        if self._lon_cache is None:
            self._materialize_centers()
        assert self._lon_cache is not None
        return self._lon_cache

    @property
    def lat(self) -> np.ndarray:
        if self._lat_cache is None:
            self._materialize_centers()
        assert self._lat_cache is not None
        return self._lat_cache

    @property
    def area(self) -> np.ndarray:
        if self._area_cache is None:
            np_dtype = _DTYPE_MAP[self.dtype]
            arr = np.zeros((6, self.Nc, self.Nc), dtype=np_dtype)
            for p in range(6):
                for i in range(self.Nc):
                    for j in range(self.Nc):
                        arr[p, i, j] = self._cell_area(p, i, j)
            self._area_cache = arr
        return self._area_cache

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering per the 2026-04-20 scope correction.

        Returns a §6-schema-shaped config: family, dimensions, extents, and
        a named generator reference. No inline geometry arrays.
        """
        return {
            "family": self.family,
            "version": "1.0.0",
            "dtype": self.dtype,
            "topology": self.topology,
            "generator": "gnomonic_c6",
            "params": {
                "Nc": int(self.Nc),
                "R": float(self.R),
                "ghosts": int(self.ghosts),
            },
            "provenance": self.provenance,
        }

    def __repr__(self) -> str:
        return (
            f"CubedSphereGrid(Nc={self.Nc}, R={self.R}, "
            f"dtype={self.dtype!r}, ghosts={self.ghosts})"
        )

    # internals -----------------------------------------------------------

    def _check_cell(self, p, i, j) -> None:
        if p is None or i is None or j is None:
            raise TypeError("cell index (p, i, j) required")
        if not 0 <= int(p) < 6:
            raise ValueError(f"panel out of range [0, 6): {p!r}")
        if not 0 <= int(i) < self.Nc:
            raise ValueError(f"i out of range [0, {self.Nc}): {i!r}")
        if not 0 <= int(j) < self.Nc:
            raise ValueError(f"j out of range [0, {self.Nc}): {j!r}")

    def _cross_panel_neighbor(
        self, p: int, direction: str, i: int, j: int
    ) -> tuple[int, int, int]:
        nb_panel, nb_edge, reverse = _PANEL_CONNECTIVITY[p][direction]
        a = j if direction in ("W", "E") else i
        if reverse:
            a = self.Nc - 1 - a
        last = self.Nc - 1
        if nb_edge == "W":
            return (nb_panel, 0, a)
        if nb_edge == "E":
            return (nb_panel, last, a)
        if nb_edge == "S":
            return (nb_panel, a, 0)
        return (nb_panel, a, last)

    def _cell_area(self, p: int, i: int, j: int) -> float:
        xi_w = float(self.xi_edges[i])
        xi_e = float(self.xi_edges[i + 1])
        eta_s = float(self.eta_edges[j])
        eta_n = float(self.eta_edges[j + 1])
        corners = [
            _gnomonic_to_cart(xi_w, eta_s, p),
            _gnomonic_to_cart(xi_e, eta_s, p),
            _gnomonic_to_cart(xi_e, eta_n, p),
            _gnomonic_to_cart(xi_w, eta_n, p),
        ]
        return _spherical_quad_area(corners, self.R)

    def _materialize_centers(self) -> None:
        Nc = self.Nc
        np_dtype = _DTYPE_MAP[self.dtype]
        lon = np.zeros((6, Nc, Nc), dtype=np_dtype)
        lat = np.zeros((6, Nc, Nc), dtype=np_dtype)
        xi_c = np.asarray(self.xi_centers, dtype=np.float64)
        eta_c = np.asarray(self.eta_centers, dtype=np.float64)
        Xi, Eta = np.meshgrid(xi_c, eta_c, indexing="ij")
        X = np.tan(Xi)
        Y = np.tan(Eta)
        r = np.sqrt(1.0 + X * X + Y * Y)
        # Per-panel Cartesian (x, y, z) on unit sphere, vectorized.
        panel_xyz = [
            (1.0 / r, X / r, Y / r),
            (-X / r, np.ones_like(r) / r, Y / r),
            (-Y / r, X / r, np.ones_like(r) / r),
            (-np.ones_like(r) / r, -X / r, Y / r),
            (X / r, -np.ones_like(r) / r, Y / r),
            (Y / r, X / r, -np.ones_like(r) / r),
        ]
        for p, (xv, yv, zv) in enumerate(panel_xyz):
            lon[p] = np.arctan2(yv, xv).astype(np_dtype, copy=False)
            lat[p] = np.arcsin(np.clip(zv, -1.0, 1.0)).astype(np_dtype, copy=False)
        self._lon_cache = lon
        self._lat_cache = lat


def cubed_sphere(
    *,
    Nc,
    R: float = 6.371e6,
    dtype="float64",
    ghosts: int = 0,
) -> CubedSphereGrid:
    """Generate a gnomonic cubed-sphere grid (6 panels x Nc x Nc cells).

    See ``docs/GRIDS_API.md`` §2.4 for the cross-binding signature contract.
    All options are keyword-only; ``Nc`` is required. The returned grid
    lowers to a small declarative ``.esm`` via :meth:`CubedSphereGrid.to_esm`;
    geometry is derived on demand.

    Parameters
    ----------
    Nc : int
        Cells per panel edge. Must be >= 1.
    R : float, optional
        Sphere radius in meters. Default 6.371e6 (Earth).
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision. Default "float64".
    ghosts : int, optional
        Halo cell width. Default 0.
    """
    if isinstance(Nc, bool) or not isinstance(Nc, (int, np.integer)):
        raise TypeError(f"Nc must be int, got {type(Nc).__name__}")
    Nc = int(Nc)
    if Nc < 1:
        raise ValueError(f"Nc must be >= 1, got {Nc}")

    if isinstance(ghosts, bool) or not isinstance(ghosts, (int, np.integer)):
        raise TypeError(f"ghosts must be int, got {type(ghosts).__name__}")
    ghosts = int(ghosts)
    if ghosts < 0:
        raise ValueError(f"ghosts must be >= 0, got {ghosts}")

    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            dtype = "float64"
        elif dtype == np.float32:
            dtype = "float32"
        else:
            raise ValueError(f"unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(f"dtype must be 'float64' or 'float32', got {dtype!r}")

    R = float(R)
    if not (R > 0 and math.isfinite(R)):
        raise ValueError(f"R must be a positive finite number, got {R!r}")

    np_dtype = _DTYPE_MAP[dtype]
    dxi = (math.pi / 2) / Nc
    xi_edges = np.asarray(
        [-math.pi / 4 + k * dxi for k in range(Nc + 1)], dtype=np_dtype
    )
    eta_edges = xi_edges.copy()
    xi_centers = 0.5 * (xi_edges[:-1] + xi_edges[1:])
    eta_centers = xi_centers.copy()

    return CubedSphereGrid(
        Nc=Nc,
        R=R,
        dtype=dtype,
        ghosts=ghosts,
        xi_edges=xi_edges,
        eta_edges=eta_edges,
        xi_centers=xi_centers,
        eta_centers=eta_centers,
    )
