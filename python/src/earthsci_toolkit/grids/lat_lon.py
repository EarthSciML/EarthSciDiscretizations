"""Lat-lon grid accessor runtime (regular + reduced-Gaussian variants).

Conforms to the cross-binding contract in ``docs/GRIDS_API.md`` §2.4, §3.2, §7.

Per the 2026-04-20 scope correction, the ``.esm`` lowering is a small
declarative config (family, dimensions, generator reference) - not a serialized
geometry blob. Geometry is derived on demand by the accessors:

* :meth:`LatLonGrid.cell_center` - cell-center ``(lon, lat)`` in radians
* :meth:`LatLonGrid.neighbors` - 4-way connectivity (longitudinally periodic,
  bounded at the poles under the default ``pole_policy="none"``)
* :meth:`LatLonGrid.metric_eval` - metric-tensor components at cell centers

Two variants are supported:

* ``variant="regular"`` - strictly uniform in lon and lat.
* ``variant="reduced_gaussian"`` - per-row ``nlon_per_row`` with optional
  user-supplied ``lat_edges`` for genuine Gaussian quadrature latitudes.
  When ``lat_edges`` is omitted it defaults to equal-angle latitudes, which
  is a test-only convenience.

Polar-singularity handling is declared via ``pole_policy`` but only
``"none"`` is implemented in this phase. Non-``"none"`` policies raise
``ValueError`` at construction so downstream callers see the
not-yet-implemented surface rather than silently incorrect geometry.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence

import numpy as np

__all__ = ["lat_lon", "LatLonGrid"]


_DTYPE_MAP: Mapping[str, type] = {"float64": np.float64, "float32": np.float32}

_DIRECTIONS = ("W", "E", "S", "N")

_VARIANTS = frozenset({"regular", "reduced_gaussian"})
_POLE_POLICIES = frozenset({"none", "average", "fold"})

_METRIC_NAMES = frozenset({
    "J",
    "g_lonlon", "g_latlat", "g_lonlat",
    "ginv_lonlon", "ginv_latlat", "ginv_lonlat",
    "area",
})


def _map_i(i: int, from_n: int, to_n: int) -> int:
    """Map column ``i`` in a row of width ``from_n`` to nearest-center column in width ``to_n``."""
    if from_n == to_n:
        return i
    frac = (i + 0.5) / from_n
    k = int(math.floor(frac * to_n))
    return min(k, to_n - 1)


class LatLonGrid:
    """Regular or reduced-Gaussian lat-lon grid.

    See ``docs/GRIDS_API.md`` §2.4, §3.2, §7. The grid stores only the
    declarative parameters plus the pre-tabulated 1-D latitude arrays; all
    per-cell geometric quantities are derived on demand, with bulk
    ``(lon, lat)`` and area arrays cached on first access.
    """

    family = "lat_lon"
    topology = "rectilinear"

    def __init__(
        self,
        *,
        variant: str,
        nlat: int,
        nlon_per_row: np.ndarray,
        R: float,
        dtype: str,
        ghosts: int,
        pole_policy: str,
        lon_start: float,
        lat_edges: np.ndarray,
        lat_centers: np.ndarray,
    ) -> None:
        self.variant = variant
        self.nlat = nlat
        self.nlon_per_row = nlon_per_row
        self.R = R
        self.dtype = dtype
        self.ghosts = ghosts
        self.pole_policy = pole_policy
        self.lon_start = lon_start
        self.lat_edges = lat_edges
        self.lat_centers = lat_centers
        self._lon_bulk_cache: np.ndarray | None = None
        self._lat_bulk_cache: np.ndarray | None = None
        self._area_bulk_cache: np.ndarray | None = None

    @property
    def n_cells(self) -> int:
        return int(self.nlon_per_row.sum())

    @property
    def nlon_uniform(self) -> int | None:
        """Uniform ``nlon`` for the regular variant; ``None`` otherwise."""
        if self.variant == "regular":
            return int(self.nlon_per_row[0])
        return None

    def nlon(self, j: int) -> int:
        """Number of longitude cells in row ``j``."""
        if not 0 <= int(j) < self.nlat:
            raise ValueError(f"j out of range [0, {self.nlat}): {j!r}")
        return int(self.nlon_per_row[j])

    @property
    def provenance(self) -> dict:
        import earthsci_toolkit

        return {
            "binding": "python",
            "binding_version": earthsci_toolkit.__version__,
            "source": "earthsci_toolkit.grids.lat_lon",
            "generator": self._generator_name(),
        }

    def lon_edges(self, j: int) -> np.ndarray:
        """Longitude edges for row ``j`` (length ``nlon_per_row[j] + 1``)."""
        n = self.nlon(j)
        dlon = 2.0 * math.pi / n
        np_dtype = _DTYPE_MAP[self.dtype]
        return np.asarray(
            [self.lon_start + k * dlon for k in range(n + 1)], dtype=np_dtype
        )

    def lon_centers(self, j: int) -> np.ndarray:
        """Longitude centers for row ``j`` (length ``nlon_per_row[j]``)."""
        n = self.nlon(j)
        dlon = 2.0 * math.pi / n
        np_dtype = _DTYPE_MAP[self.dtype]
        return np.asarray(
            [self.lon_start + (k + 0.5) * dlon for k in range(n)], dtype=np_dtype
        )

    def cell_center(self, j: int, i: int) -> tuple[float, float]:
        """Cell-center ``(lon, lat)`` for cell ``(j, i)``, both in radians."""
        self._check_cell(j, i)
        n = int(self.nlon_per_row[j])
        dlon = 2.0 * math.pi / n
        lon = self.lon_start + (int(i) + 0.5) * dlon
        return lon, float(self.lat_centers[int(j)])

    def cell_centers(
        self, j: int | None = None, i: int | None = None
    ) -> tuple[float, float] | tuple[np.ndarray, np.ndarray]:
        """Cell-center geographic coords.

        With no arguments, returns flat ragged-row-major ``(lon, lat)`` arrays
        of shape ``(n_cells,)``. With ``(j, i)``, returns the scalar
        ``(lon, lat)`` pair for that cell.
        """
        if j is None and i is None:
            return self._lon_bulk(), self._lat_bulk()
        if j is None or i is None:
            raise TypeError("cell index (j, i) required")
        return self.cell_center(int(j), int(i))

    def row_offset(self, j: int) -> int:
        """Starting flat index of row ``j`` in the ragged bulk layout."""
        if not 0 <= int(j) <= self.nlat:
            raise ValueError(f"row offset query out of range [0, {self.nlat}]: {j!r}")
        return int(self.nlon_per_row[: int(j)].sum())

    def neighbors(self, j: int, i: int) -> dict[str, tuple[int, int] | None]:
        """Face neighbors of cell ``(j, i)`` as ``{dir: (j', i')}`` (or ``None`` at a pole).

        Longitude wraps periodically. Under the default ``pole_policy="none"``
        the S-neighbor of the first row and N-neighbor of the last row are
        ``None``. For reduced-Gaussian grids the N/S neighbor is the
        nearest-center cell in the adjacent row (accounting for differing
        ``nlon``).
        """
        self._check_cell(j, i)
        jj = int(j)
        ii = int(i)
        n_i = int(self.nlon_per_row[jj])
        out: dict[str, tuple[int, int] | None] = {}
        out["W"] = (jj, n_i - 1 if ii == 0 else ii - 1)
        out["E"] = (jj, 0 if ii + 1 == n_i else ii + 1)
        if jj == 0:
            out["S"] = self._pole_neighbor()
        else:
            n_s = int(self.nlon_per_row[jj - 1])
            out["S"] = (jj - 1, _map_i(ii, n_i, n_s))
        if jj + 1 == self.nlat:
            out["N"] = self._pole_neighbor()
        else:
            n_n = int(self.nlon_per_row[jj + 1])
            out["N"] = (jj + 1, _map_i(ii, n_i, n_n))
        return out

    def cell_area(self, j: int, i: int) -> float:
        """Spherical-rectangle area for cell ``(j, i)``:
        ``R^2 * dlon * (sin(lat_n) - sin(lat_s))``.
        """
        self._check_cell(j, i)
        jj = int(j)
        n = int(self.nlon_per_row[jj])
        dlon = 2.0 * math.pi / n
        lat_s = float(self.lat_edges[jj])
        lat_n = float(self.lat_edges[jj + 1])
        return self.R * self.R * dlon * (math.sin(lat_n) - math.sin(lat_s))

    def metric_eval(self, name: str, j: int, i: int) -> float:
        """Evaluate a metric field at cell center ``(j, i)``.

        Valid ``name`` values:

        * ``"J"`` - Jacobian determinant ``R^2 |cos(lat)|``
        * ``"g_lonlon"``, ``"g_latlat"``, ``"g_lonlat"`` - covariant metric
        * ``"ginv_lonlon"``, ``"ginv_latlat"``, ``"ginv_lonlat"`` - inverse
        * ``"area"`` - spherical-rectangle cell area

        The lat-lon metric is longitudinally independent, so ``i`` is unused
        for the non-area metrics; it is still validated.
        """
        if name not in _METRIC_NAMES:
            raise ValueError(f"unknown metric name: {name!r}")
        self._check_cell(j, i)
        if name == "area":
            return self.cell_area(j, i)
        lat = float(self.lat_centers[int(j)])
        cos_lat = math.cos(lat)
        r2 = self.R * self.R
        g_ll = r2 * cos_lat * cos_lat
        g_pp = r2
        if name == "J":
            return r2 * abs(cos_lat)
        if name == "g_lonlon":
            return g_ll
        if name == "g_latlat":
            return g_pp
        if name == "g_lonlat":
            return 0.0
        if name == "ginv_lonlon":
            return 1.0 / g_ll if g_ll > 0.0 else math.inf
        if name == "ginv_latlat":
            return 1.0 / g_pp
        return 0.0  # ginv_lonlat

    @property
    def lon(self) -> np.ndarray:
        """Flat ragged-row-major array of cell-center longitudes, length ``n_cells``."""
        return self._lon_bulk()

    @property
    def lat(self) -> np.ndarray:
        """Flat ragged-row-major array of cell-center latitudes, length ``n_cells``."""
        return self._lat_bulk()

    @property
    def area(self) -> np.ndarray:
        """Flat ragged-row-major array of cell areas, length ``n_cells``."""
        if self._area_bulk_cache is None:
            np_dtype = _DTYPE_MAP[self.dtype]
            arr = np.empty(self.n_cells, dtype=np_dtype)
            k = 0
            for jj in range(self.nlat):
                n = int(self.nlon_per_row[jj])
                dlon = 2.0 * math.pi / n
                lat_s = float(self.lat_edges[jj])
                lat_n = float(self.lat_edges[jj + 1])
                a = self.R * self.R * dlon * (math.sin(lat_n) - math.sin(lat_s))
                arr[k : k + n] = a
                k += n
            self._area_bulk_cache = arr
        return self._area_bulk_cache

    def to_esm(self) -> dict:
        """Declarative ``.esm`` lowering per the 2026-04-20 scope correction.

        For the regular variant, parameters are scalar (``nlon``, ``nlat``,
        ``R``, ...). For the reduced-Gaussian variant the declaration carries
        the row-width schedule (``nlon_per_row``) and the latitude edges,
        because those two arrays *are* the reduction: declarative inputs,
        not derived geometry.
        """
        if self.variant == "regular":
            params: dict = {
                "nlon": int(self.nlon_per_row[0]),
                "nlat": int(self.nlat),
                "R": float(self.R),
                "ghosts": int(self.ghosts),
                "pole_policy": self.pole_policy,
                "lon_start": float(self.lon_start),
            }
        else:
            params = {
                "nlat": int(self.nlat),
                "nlon_per_row": [int(x) for x in self.nlon_per_row],
                "lat_edges": [float(x) for x in self.lat_edges],
                "R": float(self.R),
                "ghosts": int(self.ghosts),
                "pole_policy": self.pole_policy,
                "lon_start": float(self.lon_start),
            }
        return {
            "family": self.family,
            "version": "1.0.0",
            "dtype": self.dtype,
            "topology": self.topology,
            "variant": self.variant,
            "generator": self._generator_name(),
            "params": params,
            "provenance": self.provenance,
        }

    def __repr__(self) -> str:
        if self.variant == "regular":
            return (
                f"LatLonGrid(variant='regular', nlon={int(self.nlon_per_row[0])}, "
                f"nlat={self.nlat}, R={self.R}, dtype={self.dtype!r}, "
                f"ghosts={self.ghosts})"
            )
        return (
            f"LatLonGrid(variant='reduced_gaussian', nlat={self.nlat}, "
            f"R={self.R}, dtype={self.dtype!r}, ghosts={self.ghosts})"
        )

    # internals -----------------------------------------------------------

    def _generator_name(self) -> str:
        return (
            "lat_lon_regular"
            if self.variant == "regular"
            else "lat_lon_reduced_gaussian"
        )

    def _pole_neighbor(self) -> tuple[int, int] | None:
        # Only ``"none"`` reaches here; other policies are rejected at construction.
        return None

    def _check_cell(self, j, i) -> None:
        if j is None or i is None:
            raise TypeError("cell index (j, i) required")
        jj = int(j)
        if not 0 <= jj < self.nlat:
            raise ValueError(f"j out of range [0, {self.nlat}): {j!r}")
        n = int(self.nlon_per_row[jj])
        ii = int(i)
        if not 0 <= ii < n:
            raise ValueError(f"i out of range [0, {n}) for row {jj}: {i!r}")

    def _materialize_bulk(self) -> None:
        np_dtype = _DTYPE_MAP[self.dtype]
        lon = np.empty(self.n_cells, dtype=np_dtype)
        lat = np.empty(self.n_cells, dtype=np_dtype)
        k = 0
        for jj in range(self.nlat):
            n = int(self.nlon_per_row[jj])
            dlon = 2.0 * math.pi / n
            lat_c = float(self.lat_centers[jj])
            for ii in range(n):
                lon[k] = self.lon_start + (ii + 0.5) * dlon
                lat[k] = lat_c
                k += 1
        self._lon_bulk_cache = lon
        self._lat_bulk_cache = lat

    def _lon_bulk(self) -> np.ndarray:
        if self._lon_bulk_cache is None:
            self._materialize_bulk()
        assert self._lon_bulk_cache is not None
        return self._lon_bulk_cache

    def _lat_bulk(self) -> np.ndarray:
        if self._lat_bulk_cache is None:
            self._materialize_bulk()
        assert self._lat_bulk_cache is not None
        return self._lat_bulk_cache


def _coerce_int(name: str, value) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be int, got {type(value).__name__}")
    return int(value)


def _coerce_dtype(dtype) -> str:
    if isinstance(dtype, np.dtype):
        if dtype == np.float64:
            return "float64"
        if dtype == np.float32:
            return "float32"
        raise ValueError(f"unsupported dtype: {dtype}")
    if dtype not in _DTYPE_MAP:
        raise ValueError(f"dtype must be 'float64' or 'float32', got {dtype!r}")
    return dtype


def _coerce_finite_float(name: str, value) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a finite number, got {value!r}") from exc
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return v


def _validate_lat_edges(edges: np.ndarray, nlat: int) -> None:
    if edges.shape != (nlat + 1,):
        raise ValueError(
            f"lat_edges length {edges.shape[0]} does not match nlat+1={nlat + 1}"
        )
    if not np.all(np.isfinite(edges)):
        raise ValueError("lat_edges must be finite and strictly increasing")
    if not np.all(np.diff(edges) > 0):
        raise ValueError("lat_edges must be finite and strictly increasing")
    if edges[0] < -math.pi / 2 - 1e-12 or edges[-1] > math.pi / 2 + 1e-12:
        raise ValueError("lat_edges must lie in [-pi/2, pi/2]")


def _validate_lat_centers(centers: np.ndarray, edges: np.ndarray, nlat: int) -> None:
    if centers.shape != (nlat,):
        raise ValueError(
            f"lat_centers length {centers.shape[0]} does not match nlat={nlat}"
        )
    if not np.all(np.isfinite(centers)):
        raise ValueError("lat_centers must all be finite")
    for k in range(nlat):
        if not (edges[k] <= centers[k] <= edges[k + 1]):
            raise ValueError(
                f"lat_centers[{k}]={centers[k]!r} outside enclosing edges"
                f" [{edges[k]!r}, {edges[k + 1]!r}]"
            )


def _coerce_1d_int_array(name: str, value: Iterable) -> np.ndarray:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a sequence of ints, got {type(value).__name__}")
    try:
        arr = np.asarray(list(value), dtype=np.int64)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a sequence of ints") from exc
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1-D sequence, got shape {arr.shape}")
    return arr


def _coerce_1d_float_array(name: str, value: Iterable) -> np.ndarray:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a sequence of floats, got {type(value).__name__}")
    try:
        arr = np.asarray(list(value), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a sequence of floats") from exc
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1-D sequence, got shape {arr.shape}")
    return arr


def lat_lon(
    *,
    variant: str = "regular",
    nlon: int | None = None,
    nlat: int | None = None,
    nlon_per_row: Sequence[int] | None = None,
    lat_edges: Sequence[float] | None = None,
    lat_centers: Sequence[float] | None = None,
    R: float = 6.371e6,
    dtype="float64",
    ghosts: int = 0,
    pole_policy: str = "none",
    lon_start: float | None = None,
) -> LatLonGrid:
    """Generate a lat-lon grid (regular or reduced-Gaussian variant).

    See ``docs/GRIDS_API.md`` §2.4 for the cross-binding signature contract.
    All options are keyword-only.

    Parameters
    ----------
    variant : {"regular", "reduced_gaussian"}, optional
        Family variant. Default ``"regular"``.
    nlon : int, optional
        Cells in longitude. Required for ``"regular"``; forbidden for
        ``"reduced_gaussian"`` (use ``nlon_per_row`` instead).
    nlat : int, optional
        Cells in latitude. Required for ``"regular"``; defaults to
        ``len(nlon_per_row)`` for ``"reduced_gaussian"``.
    nlon_per_row : sequence of int, optional
        Per-row longitude cell counts. Required for ``"reduced_gaussian"``;
        forbidden for ``"regular"``.
    lat_edges : sequence of float, optional
        Explicit latitude edges (length ``nlat + 1``), strictly increasing,
        within ``[-pi/2, pi/2]``. Defaults to equal-angle edges when omitted.
    lat_centers : sequence of float, optional
        Explicit latitude centers (length ``nlat``). Defaults to the
        midpoints of the latitude edges.
    R : float, optional
        Sphere radius in meters. Default 6.371e6 (Earth).
    dtype : {"float64", "float32"} or numpy.dtype, optional
        Element precision. Default ``"float64"``.
    ghosts : int, optional
        Halo cell width. Default 0.
    pole_policy : {"none", "average", "fold"}, optional
        Polar-singularity handling. Only ``"none"`` is implemented; other
        values raise ``ValueError``. Default ``"none"``.
    lon_start : float, optional
        Starting longitude edge in radians. Default ``-pi``.
    """
    if variant not in _VARIANTS:
        raise ValueError(
            f"variant must be 'regular' or 'reduced_gaussian', got {variant!r}"
        )

    if pole_policy not in _POLE_POLICIES:
        raise ValueError(
            f"pole_policy must be one of {{'none', 'average', 'fold'}}, "
            f"got {pole_policy!r}"
        )
    if pole_policy != "none":
        raise ValueError(
            "non-'none' pole policies (average, fold) are declared but not"
            " implemented in this phase"
        )

    dtype_s = _coerce_dtype(dtype)

    R_val = _coerce_finite_float("R", R)
    if R_val <= 0.0:
        raise ValueError(f"R must be a positive finite number, got {R!r}")

    ghosts_i = _coerce_int("ghosts", ghosts)
    if ghosts_i < 0:
        raise ValueError(f"ghosts must be >= 0, got {ghosts_i}")

    lon_start_val = -math.pi if lon_start is None else _coerce_finite_float(
        "lon_start", lon_start
    )

    if variant == "regular":
        if nlon_per_row is not None:
            raise ValueError("nlon_per_row is not allowed for variant='regular'")
        if nlon is None:
            raise TypeError("missing required keyword argument: 'nlon'")
        if nlat is None:
            raise TypeError("missing required keyword argument: 'nlat'")
        nlon_i = _coerce_int("nlon", nlon)
        nlat_i = _coerce_int("nlat", nlat)
        if nlon_i < 1:
            raise ValueError(f"nlon must be >= 1, got {nlon_i}")
        if nlat_i < 1:
            raise ValueError(f"nlat must be >= 1, got {nlat_i}")
        per_row = np.full(nlat_i, nlon_i, dtype=np.int64)
    else:
        if nlon is not None:
            raise ValueError(
                "nlon is not allowed for variant='reduced_gaussian';"
                " use nlon_per_row"
            )
        if nlon_per_row is None:
            raise TypeError(
                "missing required keyword argument: 'nlon_per_row'"
            )
        per_row = _coerce_1d_int_array("nlon_per_row", nlon_per_row)
        if nlat is None:
            nlat_i = int(per_row.shape[0])
        else:
            nlat_i = _coerce_int("nlat", nlat)
        if nlat_i < 1:
            raise ValueError("nlat must be >= 1")
        if per_row.shape[0] != nlat_i:
            raise ValueError(
                f"nlon_per_row length {per_row.shape[0]} does not match"
                f" nlat={nlat_i}"
            )
        for j, n in enumerate(per_row):
            if int(n) < 1:
                raise ValueError(
                    f"nlon_per_row[{j}]={int(n)} must be >= 1"
                )

    if lat_edges is None:
        dlat = math.pi / nlat_i
        edges = np.asarray(
            [-math.pi / 2 + k * dlat for k in range(nlat_i + 1)], dtype=np.float64
        )
    else:
        edges = _coerce_1d_float_array("lat_edges", lat_edges)
        _validate_lat_edges(edges, nlat_i)

    if lat_centers is None:
        centers = 0.5 * (edges[:-1] + edges[1:])
    else:
        centers = _coerce_1d_float_array("lat_centers", lat_centers)
        _validate_lat_centers(centers, edges, nlat_i)

    return LatLonGrid(
        variant=variant,
        nlat=nlat_i,
        nlon_per_row=per_row,
        R=R_val,
        dtype=dtype_s,
        ghosts=ghosts_i,
        pole_policy=pole_policy,
        lon_start=lon_start_val,
        lat_edges=edges,
        lat_centers=centers,
    )
