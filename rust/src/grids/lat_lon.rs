//! Lat-lon grid accessor runtime (regular + reduced-Gaussian variants).
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.5, §3.3, §7.
//!
//! Per the 2026-04-20 scope correction, the `.esm` lowering is a small
//! declarative config (family, dimensions, generator reference) — not a
//! serialized geometry blob. Geometry is derived on demand by accessors:
//!
//! * [`LatLonGrid::cell_center`] — cell-center (lon, lat) in radians
//! * [`LatLonGrid::neighbors`] — 4-way connectivity (longitudinally periodic,
//!   bounded at the poles under the default [`PolePolicy::None`])
//! * [`LatLonGrid::metric_eval`] — metric-tensor components at cell centers
//!
//! Two variants are supported:
//! * [`LatLonVariant::Regular`] — strictly uniform in lon and lat.
//! * [`LatLonVariant::ReducedGaussian`] — per-row `nlon_per_row` with
//!   optional user-supplied `lat_edges` for genuine Gaussian quadrature
//!   latitudes. When `lat_edges` is omitted it defaults to equal-angle
//!   latitudes, which is a test-only convenience — production reduced-
//!   Gaussian grids must supply the quadrature nodes explicitly.
//!
//! Polar-singularity handling is declared via [`PolePolicy`] but only
//! [`PolePolicy::None`] is implemented in this phase. Non-`None` policies
//! fail at build time with a clear error so downstream callers see the
//! not-yet-implemented surface rather than silently incorrect geometry.

use std::f64::consts::{FRAC_PI_2, PI};

use serde_json::{json, Value};

use crate::{Dtype, Grid, GridError, Result};

/// Cardinal directions used in accessor return keys (matches cubed_sphere).
pub const DIRECTIONS: [&str; 4] = ["W", "E", "S", "N"];

/// Cardinal direction. Indices into the `neighbors(...)` return.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Direction {
    W,
    E,
    S,
    N,
}

impl Direction {
    pub fn idx(self) -> usize {
        match self {
            Direction::W => 0,
            Direction::E => 1,
            Direction::S => 2,
            Direction::N => 3,
        }
    }
}

/// Lat-lon family variant.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LatLonVariant {
    /// Strictly uniform in longitude and latitude.
    Regular,
    /// Per-row `nlon_per_row`; latitudes may be non-uniform (Gaussian
    /// quadrature nodes in production use).
    ReducedGaussian,
}

/// Polar-singularity handling policy.
///
/// Declared across the API surface per the phase-2 bead, but only
/// [`PolePolicy::None`] is implemented. Non-`None` variants are reserved for
/// the pole-handling follow-up bead and currently rejected at build time.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PolePolicy {
    /// Boundary cells at the poles have no north/south neighbor. This is the
    /// default and the only policy implemented in this phase.
    None,
    /// Declared: average scalar fields across the pole. Not implemented.
    Average,
    /// Declared: fold cells through the pole. Not implemented.
    Fold,
}

/// Metric field names accepted by [`LatLonGrid::metric_eval`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MetricName {
    /// Jacobian determinant of the (lon, lat) -> sphere-surface map: `R^2 cos(lat)`.
    J,
    /// Covariant metric component `g_{lonlon} = R^2 cos^2(lat)`.
    GLonLon,
    /// Covariant metric component `g_{latlat} = R^2`.
    GLatLat,
    /// Covariant metric component `g_{lonlat}` (zero for this family).
    GLonLat,
    /// Inverse-metric component `g^{lonlon}`.
    GInvLonLon,
    /// Inverse-metric component `g^{latlat}`.
    GInvLatLat,
    /// Inverse-metric component `g^{lonlat}` (zero for this family).
    GInvLonLat,
    /// Spherical-rectangle cell area (closed form for equal-lon rows).
    Area,
}

impl MetricName {
    /// Parse the canonical wire-form name (e.g. `"g_lonlon"`).
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "J" => MetricName::J,
            "g_lonlon" => MetricName::GLonLon,
            "g_latlat" => MetricName::GLatLat,
            "g_lonlat" => MetricName::GLonLat,
            "ginv_lonlon" => MetricName::GInvLonLon,
            "ginv_latlat" => MetricName::GInvLatLat,
            "ginv_lonlat" => MetricName::GInvLonLat,
            "area" => MetricName::Area,
            _ => return None,
        })
    }
}

/// Neighbor cell indexed as `(j, i)` (row, column within row).
///
/// `None` is returned at the poles under [`PolePolicy::None`] to signal a
/// genuine boundary rather than a wrap.
pub type NeighborCell = Option<(usize, usize)>;

/// Regular or reduced-Gaussian lat-lon grid.
///
/// See `docs/GRIDS_API.md` §2.5, §3.3, §7. The struct stores only the
/// declarative parameters plus the pre-tabulated 1-D latitude arrays; all
/// per-cell geometric quantities are derived on demand.
#[derive(Clone, Debug)]
pub struct LatLonGrid {
    variant: LatLonVariant,
    nlat: usize,
    /// Length `nlat`. For `Regular`, every entry equals the scalar `nlon`.
    nlon_per_row: Vec<usize>,
    r: f64,
    dtype: Dtype,
    ghosts: usize,
    pole_policy: PolePolicy,
    /// Starting longitude edge in radians. Default `-π`.
    lon_start: f64,
    /// Length `nlat + 1`; strictly increasing; clamped to `[-π/2, π/2]`.
    lat_edges: Vec<f64>,
    /// Length `nlat`; cell-center latitudes.
    lat_centers: Vec<f64>,
}

impl LatLonGrid {
    /// Variant of this grid.
    #[inline]
    pub fn variant(&self) -> LatLonVariant {
        self.variant
    }

    /// Number of latitude rows.
    #[inline]
    pub fn nlat(&self) -> usize {
        self.nlat
    }

    /// Number of longitude cells in row `j`, or `None` if `j` is out of range.
    #[inline]
    pub fn nlon(&self, j: usize) -> Option<usize> {
        self.nlon_per_row.get(j).copied()
    }

    /// Uniform `nlon` for the `Regular` variant, `None` otherwise.
    #[inline]
    pub fn nlon_uniform(&self) -> Option<usize> {
        match self.variant {
            LatLonVariant::Regular => self.nlon_per_row.first().copied(),
            LatLonVariant::ReducedGaussian => None,
        }
    }

    /// Sphere radius (meters).
    #[inline]
    pub fn r(&self) -> f64 {
        self.r
    }

    /// Halo cell width.
    #[inline]
    pub fn ghosts(&self) -> usize {
        self.ghosts
    }

    /// Pole policy (always [`PolePolicy::None`] in this phase).
    #[inline]
    pub fn pole_policy(&self) -> PolePolicy {
        self.pole_policy
    }

    /// Topology string (always `"rectilinear"` for this family per §7).
    #[inline]
    pub fn topology(&self) -> &'static str {
        "rectilinear"
    }

    /// Total cell count across all rows.
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nlon_per_row.iter().sum()
    }

    /// Read-only view of the latitude edges (length `nlat + 1`).
    #[inline]
    pub fn lat_edges(&self) -> &[f64] {
        &self.lat_edges
    }

    /// Read-only view of the latitude centers (length `nlat`).
    #[inline]
    pub fn lat_centers(&self) -> &[f64] {
        &self.lat_centers
    }

    /// Read-only view of the per-row longitude cell counts (length `nlat`).
    #[inline]
    pub fn nlon_per_row(&self) -> &[usize] {
        &self.nlon_per_row
    }

    /// Starting longitude edge in radians.
    #[inline]
    pub fn lon_start(&self) -> f64 {
        self.lon_start
    }

    /// Longitude edges for row `j` (length `nlon_per_row[j] + 1`).
    pub fn lon_edges(&self, j: usize) -> Result<Vec<f64>> {
        let n = self.row_nlon(j)?;
        let dlon = 2.0 * PI / n as f64;
        Ok((0..=n).map(|k| self.lon_start + k as f64 * dlon).collect())
    }

    /// Longitude centers for row `j` (length `nlon_per_row[j]`).
    pub fn lon_centers(&self, j: usize) -> Result<Vec<f64>> {
        let n = self.row_nlon(j)?;
        let dlon = 2.0 * PI / n as f64;
        Ok((0..n)
            .map(|k| self.lon_start + (k as f64 + 0.5) * dlon)
            .collect())
    }

    /// Cell-center `(lon, lat)` for cell `(j, i)`, both in radians.
    pub fn cell_center(&self, j: usize, i: usize) -> Result<(f64, f64)> {
        self.check_cell(j, i)?;
        let n = self.nlon_per_row[j];
        let dlon = 2.0 * PI / n as f64;
        let lon = self.lon_start + (i as f64 + 0.5) * dlon;
        Ok((lon, self.lat_centers[j]))
    }

    /// Bulk cell-center `(lon, lat)` arrays, row-major: `[row 0 | row 1 | …]`.
    ///
    /// For the `Regular` variant all rows have the same length and the output
    /// is equivalent to a flattened `(nlat, nlon)` array; for
    /// `ReducedGaussian` the output is a flat ragged concatenation whose row
    /// starts are `sum(nlon_per_row[..j])`.
    pub fn cell_centers_bulk(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.n_cells();
        let mut lon = Vec::with_capacity(n);
        let mut lat = Vec::with_capacity(n);
        for j in 0..self.nlat {
            let nlon = self.nlon_per_row[j];
            let dlon = 2.0 * PI / nlon as f64;
            let lat_c = self.lat_centers[j];
            for i in 0..nlon {
                lon.push(self.lon_start + (i as f64 + 0.5) * dlon);
                lat.push(lat_c);
            }
        }
        (lon, lat)
    }

    /// Starting flat index of row `j` in the ragged bulk layout.
    pub fn row_offset(&self, j: usize) -> Result<usize> {
        if j > self.nlat {
            return Err(GridError::InvalidOption(
                "j",
                format!("row offset query out of range [0, {}]: {j}", self.nlat),
            ));
        }
        Ok(self.nlon_per_row[..j].iter().sum())
    }

    /// Face neighbors of cell `(j, i)` as `[W, E, S, N]`.
    ///
    /// Longitude wraps periodically. Under the default [`PolePolicy::None`]
    /// the S-neighbor of the first row and the N-neighbor of the last row
    /// are `None`. For reduced-Gaussian grids the N/S neighbor is the
    /// nearest-center cell in the adjacent row (accounting for differing
    /// `nlon`).
    pub fn neighbors(&self, j: usize, i: usize) -> Result<[NeighborCell; 4]> {
        self.check_cell(j, i)?;
        let n_i = self.nlon_per_row[j];
        let w = Some((j, if i == 0 { n_i - 1 } else { i - 1 }));
        let e = Some((j, if i + 1 == n_i { 0 } else { i + 1 }));
        let s = if j == 0 {
            self.pole_neighbor()
        } else {
            let n_s = self.nlon_per_row[j - 1];
            Some((j - 1, map_i(i, n_i, n_s)))
        };
        let n = if j + 1 == self.nlat {
            self.pole_neighbor()
        } else {
            let n_n = self.nlon_per_row[j + 1];
            Some((j + 1, map_i(i, n_i, n_n)))
        };
        Ok([w, e, s, n])
    }

    /// Return a single neighbor in direction `dir`.
    pub fn neighbor(&self, j: usize, i: usize, dir: Direction) -> Result<NeighborCell> {
        Ok(self.neighbors(j, i)?[dir.idx()])
    }

    /// Spherical-rectangle area for cell `(j, i)`:
    /// `R² · Δlon · (sin(lat_n) − sin(lat_s))`.
    pub fn cell_area(&self, j: usize, i: usize) -> Result<f64> {
        self.check_cell(j, i)?;
        let n = self.nlon_per_row[j];
        let dlon = 2.0 * PI / n as f64;
        let lat_s = self.lat_edges[j];
        let lat_n = self.lat_edges[j + 1];
        Ok(self.r * self.r * dlon * (lat_n.sin() - lat_s.sin()))
    }

    /// Cell areas for every cell, in the same ragged row-major order as
    /// [`Self::cell_centers_bulk`].
    pub fn area_bulk(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.n_cells());
        for j in 0..self.nlat {
            let n = self.nlon_per_row[j];
            let dlon = 2.0 * PI / n as f64;
            let lat_s = self.lat_edges[j];
            let lat_n = self.lat_edges[j + 1];
            let area = self.r * self.r * dlon * (lat_n.sin() - lat_s.sin());
            for _ in 0..n {
                out.push(area);
            }
        }
        out
    }

    /// Evaluate a metric field at cell center `(j, i)`.
    ///
    /// The lat-lon metric is longitudinally independent, so `i` is unused
    /// for the non-area metrics; it is still validated.
    pub fn metric_eval(&self, name: MetricName, j: usize, i: usize) -> Result<f64> {
        self.check_cell(j, i)?;
        if name == MetricName::Area {
            return self.cell_area(j, i);
        }
        let lat = self.lat_centers[j];
        let cos_lat = lat.cos();
        let r2 = self.r * self.r;
        let g_ll = r2 * cos_lat * cos_lat;
        let g_pp = r2;
        let jac = r2 * cos_lat.abs();
        Ok(match name {
            MetricName::J => jac,
            MetricName::GLonLon => g_ll,
            MetricName::GLatLat => g_pp,
            MetricName::GLonLat => 0.0,
            MetricName::GInvLonLon => {
                if g_ll > 0.0 {
                    1.0 / g_ll
                } else {
                    f64::INFINITY
                }
            }
            MetricName::GInvLatLat => 1.0 / g_pp,
            MetricName::GInvLonLat => 0.0,
            MetricName::Area => unreachable!("handled above"),
        })
    }

    /// Evaluate a metric by its canonical wire-form name (`"g_lonlon"` etc.).
    pub fn metric_eval_by_name(&self, name: &str, j: usize, i: usize) -> Result<f64> {
        let m = MetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, j, i)
    }

    /// Provenance block per `docs/GRIDS_API.md` §6.4.
    pub fn provenance(&self) -> Value {
        json!({
            "binding": "rust",
            "binding_version": env!("CARGO_PKG_VERSION"),
            "source": "earthsci_grids::grids::lat_lon",
            "generator": self.generator_name(),
        })
    }

    // internals -----------------------------------------------------------

    fn generator_name(&self) -> &'static str {
        match self.variant {
            LatLonVariant::Regular => "lat_lon_regular",
            LatLonVariant::ReducedGaussian => "lat_lon_reduced_gaussian",
        }
    }

    fn variant_str(&self) -> &'static str {
        match self.variant {
            LatLonVariant::Regular => "regular",
            LatLonVariant::ReducedGaussian => "reduced_gaussian",
        }
    }

    fn pole_policy_str(&self) -> &'static str {
        match self.pole_policy {
            PolePolicy::None => "none",
            PolePolicy::Average => "average",
            PolePolicy::Fold => "fold",
        }
    }

    fn row_nlon(&self, j: usize) -> Result<usize> {
        self.nlon_per_row.get(j).copied().ok_or_else(|| {
            GridError::InvalidOption("j", format!("j out of range [0, {}): {j}", self.nlat))
        })
    }

    fn check_cell(&self, j: usize, i: usize) -> Result<()> {
        if j >= self.nlat {
            return Err(GridError::InvalidOption(
                "j",
                format!("j out of range [0, {}): {j}", self.nlat),
            ));
        }
        let n = self.nlon_per_row[j];
        if i >= n {
            return Err(GridError::InvalidOption(
                "i",
                format!("i out of range [0, {n}) for row {j}: {i}"),
            ));
        }
        Ok(())
    }

    fn pole_neighbor(&self) -> NeighborCell {
        // Only `None` reaches here; non-`None` policies are rejected at build time.
        None
    }
}

/// Map a column index `i` in a row of width `from` to the nearest-center
/// column in a row of width `to`.
#[inline]
fn map_i(i: usize, from: usize, to: usize) -> usize {
    if from == to {
        return i;
    }
    let frac = (i as f64 + 0.5) / from as f64;
    let k = (frac * to as f64).floor() as usize;
    k.min(to - 1)
}

impl Grid for LatLonGrid {
    fn family(&self) -> &'static str {
        "lat_lon"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative `.esm` lowering per the 2026-04-20 scope correction.
    ///
    /// For `Regular`, parameters are scalar (`nlon`, `nlat`, `R`, …). For
    /// `ReducedGaussian`, the declaration carries the row-width schedule
    /// (`nlon_per_row`) and the latitude edges, because those two arrays
    /// *are* the reduction — they are declarative inputs, not derived
    /// geometry.
    fn to_esm(&self) -> Value {
        let params = match self.variant {
            LatLonVariant::Regular => json!({
                "nlon": self.nlon_per_row[0],
                "nlat": self.nlat,
                "R": self.r,
                "ghosts": self.ghosts,
                "pole_policy": self.pole_policy_str(),
                "lon_start": self.lon_start,
            }),
            LatLonVariant::ReducedGaussian => json!({
                "nlat": self.nlat,
                "nlon_per_row": self.nlon_per_row,
                "lat_edges": self.lat_edges,
                "R": self.r,
                "ghosts": self.ghosts,
                "pole_policy": self.pole_policy_str(),
                "lon_start": self.lon_start,
            }),
        };
        json!({
            "family": self.family(),
            "version": "1.0.0",
            "dtype": dtype_str(self.dtype),
            "topology": self.topology(),
            "variant": self.variant_str(),
            "generator": self.generator_name(),
            "params": params,
            "provenance": self.provenance(),
        })
    }
}

fn dtype_str(d: Dtype) -> &'static str {
    match d {
        Dtype::F64 => "float64",
        Dtype::F32 => "float32",
    }
}

// --- builder ---------------------------------------------------------------

/// Create a new [`LatLonBuilder`] with default options.
///
/// ```
/// use earthsci_grids::{lat_lon, Dtype, Grid};
/// let g = lat_lon::builder().nlon(72).nlat(36).build().unwrap();
/// assert_eq!(g.family(), "lat_lon");
/// assert_eq!(g.n_cells(), 72 * 36);
/// assert_eq!(g.dtype(), Dtype::F64);
/// ```
pub fn builder() -> LatLonBuilder {
    LatLonBuilder::default()
}

/// Builder for [`LatLonGrid`]. Required options are validated at
/// [`LatLonBuilder::build`] time so partial builders can be reused in tests.
#[must_use]
#[derive(Clone, Debug, Default)]
pub struct LatLonBuilder {
    variant: Option<LatLonVariant>,
    nlon: Option<usize>,
    nlat: Option<usize>,
    nlon_per_row: Option<Vec<usize>>,
    lat_edges: Option<Vec<f64>>,
    lat_centers: Option<Vec<f64>>,
    r: Option<f64>,
    dtype: Option<Dtype>,
    ghosts: Option<usize>,
    pole_policy: Option<PolePolicy>,
    lon_start: Option<f64>,
}

impl LatLonBuilder {
    /// Select the grid variant. Defaults to [`LatLonVariant::Regular`].
    pub fn variant(mut self, v: LatLonVariant) -> Self {
        self.variant = Some(v);
        self
    }

    /// Cells in longitude. Required for `Regular`; forbidden for
    /// `ReducedGaussian` (use [`Self::nlon_per_row`] instead).
    pub fn nlon(mut self, nlon: usize) -> Self {
        self.nlon = Some(nlon);
        self
    }

    /// Cells in latitude. Required for `Regular`; defaults to
    /// `nlon_per_row.len()` for `ReducedGaussian`.
    pub fn nlat(mut self, nlat: usize) -> Self {
        self.nlat = Some(nlat);
        self
    }

    /// Per-row longitude cell counts. Required for `ReducedGaussian`;
    /// forbidden for `Regular`.
    pub fn nlon_per_row(mut self, v: Vec<usize>) -> Self {
        self.nlon_per_row = Some(v);
        self
    }

    /// Explicit latitude edges (length `nlat + 1`), strictly increasing,
    /// within `[-π/2, π/2]`. Defaults to equal-angle edges when omitted.
    pub fn lat_edges(mut self, v: Vec<f64>) -> Self {
        self.lat_edges = Some(v);
        self
    }

    /// Explicit latitude centers (length `nlat`). Defaults to the midpoints
    /// of the latitude edges.
    pub fn lat_centers(mut self, v: Vec<f64>) -> Self {
        self.lat_centers = Some(v);
        self
    }

    /// Sphere radius in meters. Defaults to `6.371e6` (Earth).
    pub fn r(mut self, r: f64) -> Self {
        self.r = Some(r);
        self
    }

    /// Element precision. Defaults to [`Dtype::F64`].
    pub fn dtype(mut self, dtype: Dtype) -> Self {
        self.dtype = Some(dtype);
        self
    }

    /// Halo cell width. Defaults to 0.
    pub fn ghosts(mut self, ghosts: usize) -> Self {
        self.ghosts = Some(ghosts);
        self
    }

    /// Pole-singularity policy. Only [`PolePolicy::None`] is currently
    /// implemented; others fail at `build()`.
    pub fn pole_policy(mut self, p: PolePolicy) -> Self {
        self.pole_policy = Some(p);
        self
    }

    /// Starting longitude edge in radians. Defaults to `-π`.
    pub fn lon_start(mut self, s: f64) -> Self {
        self.lon_start = Some(s);
        self
    }

    /// Validate options and construct the grid.
    pub fn build(self) -> Result<LatLonGrid> {
        let variant = self.variant.unwrap_or(LatLonVariant::Regular);

        let r = self.r.unwrap_or(6.371e6);
        if !(r > 0.0 && r.is_finite()) {
            return Err(GridError::InvalidOption(
                "R",
                format!("must be a positive finite number, got {r}"),
            ));
        }

        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);

        let pole_policy = self.pole_policy.unwrap_or(PolePolicy::None);
        if !matches!(pole_policy, PolePolicy::None) {
            return Err(GridError::InvalidOption(
                "pole_policy",
                "non-None pole policies (average, fold) are declared but not implemented"
                    .to_string(),
            ));
        }

        let lon_start = self.lon_start.unwrap_or(-PI);
        if !lon_start.is_finite() {
            return Err(GridError::InvalidOption(
                "lon_start",
                format!("must be finite, got {lon_start}"),
            ));
        }

        let (nlat, nlon_per_row) = match variant {
            LatLonVariant::Regular => {
                if self.nlon_per_row.is_some() {
                    return Err(GridError::InvalidOption(
                        "nlon_per_row",
                        "not allowed for variant=regular".to_string(),
                    ));
                }
                let nlon = self.nlon.ok_or(GridError::MissingOption("nlon"))?;
                let nlat = self.nlat.ok_or(GridError::MissingOption("nlat"))?;
                if nlon < 1 {
                    return Err(GridError::InvalidOption(
                        "nlon",
                        format!("must be >= 1, got {nlon}"),
                    ));
                }
                if nlat < 1 {
                    return Err(GridError::InvalidOption(
                        "nlat",
                        format!("must be >= 1, got {nlat}"),
                    ));
                }
                (nlat, vec![nlon; nlat])
            }
            LatLonVariant::ReducedGaussian => {
                if self.nlon.is_some() {
                    return Err(GridError::InvalidOption(
                        "nlon",
                        "not allowed for variant=reduced_gaussian; use nlon_per_row".to_string(),
                    ));
                }
                let nlon_per_row = self
                    .nlon_per_row
                    .ok_or(GridError::MissingOption("nlon_per_row"))?;
                let nlat = self.nlat.unwrap_or(nlon_per_row.len());
                if nlat == 0 {
                    return Err(GridError::InvalidOption("nlat", "must be >= 1".to_string()));
                }
                if nlon_per_row.len() != nlat {
                    return Err(GridError::InvalidOption(
                        "nlon_per_row",
                        format!("length {} does not match nlat={nlat}", nlon_per_row.len()),
                    ));
                }
                for (j, &n) in nlon_per_row.iter().enumerate() {
                    if n < 1 {
                        return Err(GridError::InvalidOption(
                            "nlon_per_row",
                            format!("row {j} has nlon={n}; must be >= 1"),
                        ));
                    }
                }
                (nlat, nlon_per_row)
            }
        };

        let lat_edges = if let Some(edges) = self.lat_edges {
            validate_lat_edges(&edges, nlat)?;
            edges
        } else {
            let dlat = PI / nlat as f64;
            (0..=nlat).map(|k| -FRAC_PI_2 + k as f64 * dlat).collect()
        };

        let lat_centers = if let Some(centers) = self.lat_centers {
            validate_lat_centers(&centers, &lat_edges, nlat)?;
            centers
        } else {
            (0..nlat)
                .map(|k| 0.5 * (lat_edges[k] + lat_edges[k + 1]))
                .collect()
        };

        Ok(LatLonGrid {
            variant,
            nlat,
            nlon_per_row,
            r,
            dtype,
            ghosts,
            pole_policy,
            lon_start,
            lat_edges,
            lat_centers,
        })
    }
}

fn validate_lat_edges(edges: &[f64], nlat: usize) -> Result<()> {
    if edges.len() != nlat + 1 {
        return Err(GridError::InvalidOption(
            "lat_edges",
            format!("length {} does not match nlat+1={}", edges.len(), nlat + 1),
        ));
    }
    for k in 0..nlat {
        if edges[k] >= edges[k + 1] || !edges[k].is_finite() {
            return Err(GridError::InvalidOption(
                "lat_edges",
                "must be finite and strictly increasing".to_string(),
            ));
        }
    }
    if !edges[nlat].is_finite() {
        return Err(GridError::InvalidOption(
            "lat_edges",
            "must be finite and strictly increasing".to_string(),
        ));
    }
    if edges[0] < -FRAC_PI_2 - 1e-12 || edges[nlat] > FRAC_PI_2 + 1e-12 {
        return Err(GridError::InvalidOption(
            "lat_edges",
            "must lie in [-pi/2, pi/2]".to_string(),
        ));
    }
    Ok(())
}

fn validate_lat_centers(centers: &[f64], edges: &[f64], nlat: usize) -> Result<()> {
    if centers.len() != nlat {
        return Err(GridError::InvalidOption(
            "lat_centers",
            format!("length {} does not match nlat={nlat}", centers.len()),
        ));
    }
    for k in 0..nlat {
        if !(centers[k].is_finite() && edges[k] <= centers[k] && centers[k] <= edges[k + 1]) {
            return Err(GridError::InvalidOption(
                "lat_centers",
                format!(
                    "center {k}={} outside enclosing edges [{}, {}]",
                    centers[k],
                    edges[k],
                    edges[k + 1]
                ),
            ));
        }
    }
    Ok(())
}

// --- tests -----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn g_regular(nlon: usize, nlat: usize) -> LatLonGrid {
        builder().nlon(nlon).nlat(nlat).build().unwrap()
    }

    // --- API ---

    #[test]
    fn builder_requires_nlon() {
        let err = builder().nlat(4).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("nlon")));
    }

    #[test]
    fn builder_requires_nlat() {
        let err = builder().nlon(4).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("nlat")));
    }

    #[test]
    fn builder_rejects_zero_nlon() {
        let err = builder().nlon(0).nlat(4).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nlon", _)));
    }

    #[test]
    fn builder_rejects_zero_nlat() {
        let err = builder().nlon(4).nlat(0).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nlat", _)));
    }

    #[test]
    fn builder_rejects_nonfinite_r() {
        let err = builder().nlon(4).nlat(4).r(f64::NAN).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));
        let err = builder().nlon(4).nlat(4).r(-1.0).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));
    }

    #[test]
    fn builder_rejects_unimplemented_pole_policy() {
        let err = builder()
            .nlon(4)
            .nlat(4)
            .pole_policy(PolePolicy::Average)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("pole_policy", _)));
        let err = builder()
            .nlon(4)
            .nlat(4)
            .pole_policy(PolePolicy::Fold)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("pole_policy", _)));
    }

    #[test]
    fn builder_rejects_nlon_per_row_in_regular() {
        let err = builder()
            .nlon(4)
            .nlat(2)
            .nlon_per_row(vec![4, 4])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nlon_per_row", _)));
    }

    #[test]
    fn builder_rejects_nlon_in_reduced_gaussian() {
        let err = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlon(4)
            .nlon_per_row(vec![4, 8])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nlon", _)));
    }

    #[test]
    fn reduced_gaussian_requires_nlon_per_row() {
        let err = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlat(4)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::MissingOption("nlon_per_row")));
    }

    #[test]
    fn reduced_gaussian_length_mismatch_errors() {
        let err = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlat(3)
            .nlon_per_row(vec![4, 8])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nlon_per_row", _)));
    }

    #[test]
    fn defaults_match_contract() {
        let g = builder().nlon(8).nlat(4).build().unwrap();
        assert_eq!(g.variant(), LatLonVariant::Regular);
        assert_eq!(g.nlat(), 4);
        assert_eq!(g.nlon_uniform(), Some(8));
        assert_eq!(g.r(), 6.371e6);
        assert_eq!(g.ghosts(), 0);
        assert_eq!(g.dtype(), Dtype::F64);
        assert_eq!(g.family(), "lat_lon");
        assert_eq!(g.topology(), "rectilinear");
        assert_eq!(g.pole_policy(), PolePolicy::None);
        assert_eq!(g.n_cells(), 8 * 4);
    }

    #[test]
    fn lat_edges_default_span_pole_to_pole() {
        let g = g_regular(4, 6);
        let e = g.lat_edges();
        assert_eq!(e.len(), 7);
        assert_relative_eq!(e[0], -FRAC_PI_2, epsilon = 1e-15);
        assert_relative_eq!(e[6], FRAC_PI_2, epsilon = 1e-15);
    }

    #[test]
    fn invalid_lat_edges_rejected() {
        // not increasing
        let err = builder()
            .nlon(2)
            .nlat(2)
            .lat_edges(vec![0.0, 0.0, 0.5])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("lat_edges", _)));
        // out of range
        let err = builder()
            .nlon(2)
            .nlat(1)
            .lat_edges(vec![-2.0, 2.0])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("lat_edges", _)));
    }

    // --- topology ---

    #[test]
    fn interior_neighbors_local() {
        let g = g_regular(8, 6);
        let ns = g.neighbors(2, 3).unwrap();
        assert_eq!(ns[Direction::W.idx()], Some((2, 2)));
        assert_eq!(ns[Direction::E.idx()], Some((2, 4)));
        assert_eq!(ns[Direction::S.idx()], Some((1, 3)));
        assert_eq!(ns[Direction::N.idx()], Some((3, 3)));
    }

    #[test]
    fn longitude_wraps_periodically() {
        let g = g_regular(8, 4);
        let ns = g.neighbors(2, 0).unwrap();
        assert_eq!(ns[Direction::W.idx()], Some((2, 7)));
        let ns = g.neighbors(2, 7).unwrap();
        assert_eq!(ns[Direction::E.idx()], Some((2, 0)));
    }

    #[test]
    fn south_pole_has_no_s_neighbor() {
        let g = g_regular(8, 4);
        let ns = g.neighbors(0, 3).unwrap();
        assert_eq!(ns[Direction::S.idx()], None);
        assert_eq!(ns[Direction::N.idx()], Some((1, 3)));
    }

    #[test]
    fn north_pole_has_no_n_neighbor() {
        let g = g_regular(8, 4);
        let ns = g.neighbors(3, 3).unwrap();
        assert_eq!(ns[Direction::N.idx()], None);
        assert_eq!(ns[Direction::S.idx()], Some((2, 3)));
    }

    #[test]
    fn neighbor_cells_in_range() {
        let g = g_regular(5, 4);
        for j in 0..g.nlat() {
            let nlon = g.nlon(j).unwrap();
            for i in 0..nlon {
                for (jj, ii) in g.neighbors(j, i).unwrap().into_iter().flatten() {
                    assert!(jj < g.nlat());
                    assert!(ii < g.nlon(jj).unwrap());
                }
            }
        }
    }

    #[test]
    fn reduced_gaussian_neighbor_rounding() {
        let g = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlon_per_row(vec![4, 8, 4])
            .build()
            .unwrap();
        // Row 1 cell (1, 5) -> center frac ((5 + 0.5)/8) = 0.6875 -> row 0 (width 4)
        // column floor(0.6875 * 4) = 2.
        let ns = g.neighbors(1, 5).unwrap();
        assert_eq!(ns[Direction::S.idx()], Some((0, 2)));
        // (row 1 cell 0) center frac = 0.0625 -> row 0 col 0.
        let ns = g.neighbors(1, 0).unwrap();
        assert_eq!(ns[Direction::S.idx()], Some((0, 0)));
    }

    #[test]
    fn reduced_gaussian_cell_counts() {
        let g = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlon_per_row(vec![4, 8, 4])
            .build()
            .unwrap();
        assert_eq!(g.n_cells(), 16);
        assert_eq!(g.nlat(), 3);
        assert_eq!(g.row_offset(0).unwrap(), 0);
        assert_eq!(g.row_offset(1).unwrap(), 4);
        assert_eq!(g.row_offset(2).unwrap(), 12);
    }

    // --- centers ---

    #[test]
    fn cell_center_at_equator_equals_row_center_lat() {
        let g = g_regular(4, 2);
        let (_, lat) = g.cell_center(0, 0).unwrap();
        // nlat=2 -> edges at [-pi/2, 0, pi/2]; row 0 center = -pi/4.
        assert_relative_eq!(lat, -std::f64::consts::FRAC_PI_4, epsilon = 1e-15);
        let (_, lat) = g.cell_center(1, 0).unwrap();
        assert_relative_eq!(lat, std::f64::consts::FRAC_PI_4, epsilon = 1e-15);
    }

    #[test]
    fn cell_centers_bulk_matches_scalar() {
        let g = g_regular(6, 4);
        let (lon, lat) = g.cell_centers_bulk();
        for j in [0_usize, 1, 3] {
            for i in [0_usize, 2, 5] {
                let (s_lon, s_lat) = g.cell_center(j, i).unwrap();
                let idx = g.row_offset(j).unwrap() + i;
                assert_relative_eq!(s_lon, lon[idx], epsilon = 1e-14);
                assert_relative_eq!(s_lat, lat[idx], epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn cell_centers_lon_in_range() {
        let g = g_regular(12, 6);
        let (lon, lat) = g.cell_centers_bulk();
        let lon_range = (-PI - 1e-12)..=(PI + 1e-12);
        let lat_range = (-FRAC_PI_2 - 1e-12)..=(FRAC_PI_2 + 1e-12);
        for &v in &lon {
            assert!(lon_range.contains(&v));
        }
        for &v in &lat {
            assert!(lat_range.contains(&v));
        }
    }

    // --- area ---

    #[test]
    fn total_sphere_area_matches_4_pi_r_squared() {
        let r_val = 6.371e6_f64;
        let g = builder().nlon(72).nlat(36).r(r_val).build().unwrap();
        let total: f64 = g.area_bulk().iter().sum();
        let expected = 4.0 * PI * r_val * r_val;
        assert!((total - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn total_unit_sphere_area_is_4pi() {
        let g = builder().nlon(10).nlat(8).r(1.0).build().unwrap();
        let total: f64 = g.area_bulk().iter().sum();
        assert!((total - 4.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn reduced_gaussian_total_area_is_4pi_r_squared() {
        let g = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlon_per_row(vec![4, 8, 12, 12, 8, 4])
            .r(1.0)
            .build()
            .unwrap();
        let total: f64 = g.area_bulk().iter().sum();
        assert!((total - 4.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn area_positive_everywhere() {
        let g = g_regular(6, 4);
        for a in g.area_bulk() {
            assert!(a > 0.0);
        }
    }

    #[test]
    fn area_via_metric_matches_bulk() {
        let g = g_regular(6, 4);
        let bulk = g.area_bulk();
        for j in [0_usize, 1, 3] {
            for i in [0_usize, 2, 5] {
                let m = g.metric_eval(MetricName::Area, j, i).unwrap();
                let idx = g.row_offset(j).unwrap() + i;
                assert_relative_eq!(m, bulk[idx], max_relative = 1e-12);
            }
        }
    }

    // --- metric ---

    #[test]
    fn metric_at_equator_row_is_r_squared() {
        // Row near lat=0 -> cos(lat) ~= 1, so g_lonlon ~= g_latlat = R^2.
        // Large nlat pushes the equator row arbitrarily close to lat=0.
        let r = 3.0_f64;
        let g = builder().nlon(4).nlat(400).r(r).build().unwrap();
        let j_eq = 200; // lat center close to 0
        let g_ll = g.metric_eval(MetricName::GLonLon, j_eq, 0).unwrap();
        let g_pp = g.metric_eval(MetricName::GLatLat, j_eq, 0).unwrap();
        assert!((g_ll - r * r).abs() / (r * r) < 1e-4);
        assert_relative_eq!(g_pp, r * r, epsilon = 1e-15);
    }

    #[test]
    fn metric_g_lonlat_is_zero() {
        let g = g_regular(4, 4);
        for j in 0..g.nlat() {
            assert_eq!(g.metric_eval(MetricName::GLonLat, j, 0).unwrap(), 0.0);
            assert_eq!(g.metric_eval(MetricName::GInvLonLat, j, 0).unwrap(), 0.0);
        }
    }

    #[test]
    fn inverse_metric_is_inverse() {
        let g = builder().nlon(6).nlat(6).r(1.0).build().unwrap();
        for j in [0_usize, 2, 5] {
            let g_ll = g.metric_eval(MetricName::GLonLon, j, 0).unwrap();
            let g_pp = g.metric_eval(MetricName::GLatLat, j, 0).unwrap();
            let i_ll = g.metric_eval(MetricName::GInvLonLon, j, 0).unwrap();
            let i_pp = g.metric_eval(MetricName::GInvLatLat, j, 0).unwrap();
            assert_relative_eq!(g_ll * i_ll, 1.0, epsilon = 1e-12);
            assert_relative_eq!(g_pp * i_pp, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn metric_eval_rejects_out_of_range_cell() {
        let g = g_regular(4, 4);
        assert!(g.metric_eval(MetricName::J, 4, 0).is_err());
        assert!(g.metric_eval(MetricName::J, 0, 4).is_err());
    }

    #[test]
    fn metric_eval_by_name_rejects_unknown() {
        let g = g_regular(4, 4);
        let err = g.metric_eval_by_name("not_a_metric", 0, 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn j_positive_everywhere() {
        let g = g_regular(6, 6);
        for j in 0..g.nlat() {
            for i in 0..g.nlon(j).unwrap() {
                assert!(g.metric_eval(MetricName::J, j, i).unwrap() > 0.0);
            }
        }
    }

    // --- to_esm ---

    #[test]
    fn to_esm_regular_is_declarative() {
        let g = builder()
            .nlon(72)
            .nlat(36)
            .r(6.371e6)
            .ghosts(2)
            .build()
            .unwrap();
        let doc = g.to_esm();
        assert_eq!(doc["family"], "lat_lon");
        assert_eq!(doc["version"], "1.0.0");
        assert_eq!(doc["dtype"], "float64");
        assert_eq!(doc["topology"], "rectilinear");
        assert_eq!(doc["variant"], "regular");
        assert_eq!(doc["generator"], "lat_lon_regular");
        assert_eq!(doc["params"]["nlon"], 72);
        assert_eq!(doc["params"]["nlat"], 36);
        assert_eq!(doc["params"]["R"], 6.371e6);
        assert_eq!(doc["params"]["ghosts"], 2);
        assert_eq!(doc["params"]["pole_policy"], "none");
        assert!(doc.get("provenance").is_some());

        let wire = serde_json::to_string(&doc).unwrap();
        assert!(!wire.contains("cells"));
        assert!(!wire.contains("lon_array"));
        assert!(wire.len() < 2_000, "wire too large: {} bytes", wire.len());
    }

    #[test]
    fn to_esm_reduced_gaussian_carries_schedule() {
        let g = builder()
            .variant(LatLonVariant::ReducedGaussian)
            .nlon_per_row(vec![4, 8, 4])
            .r(1.0)
            .build()
            .unwrap();
        let doc = g.to_esm();
        assert_eq!(doc["variant"], "reduced_gaussian");
        assert_eq!(doc["generator"], "lat_lon_reduced_gaussian");
        assert_eq!(doc["params"]["nlat"], 3);
        assert_eq!(doc["params"]["nlon_per_row"], json!([4, 8, 4]));
        assert!(doc["params"]["lat_edges"].is_array());
    }

    #[test]
    fn to_esm_roundtrips_via_json() {
        let g = builder().nlon(4).nlat(3).r(1.0).build().unwrap();
        let doc = g.to_esm();
        let text = serde_json::to_string(&doc).unwrap();
        let reparsed: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(reparsed["params"]["nlon"], 4);
        assert_eq!(reparsed["params"]["nlat"], 3);
        assert_eq!(reparsed["params"]["R"], 1.0);
    }

    #[test]
    fn provenance_identifies_rust_binding() {
        let g = g_regular(4, 4);
        let doc = g.to_esm();
        let prov = &doc["provenance"];
        assert_eq!(prov["binding"], "rust");
        assert_eq!(prov["binding_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(prov["generator"], "lat_lon_regular");
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let g = builder().nlon(4).nlat(4).dtype(Dtype::F32).build().unwrap();
        assert_eq!(g.dtype(), Dtype::F32);
        assert_eq!(g.to_esm()["dtype"], "float32");
    }

    // --- map_i helper ---

    #[test]
    fn map_i_identity_when_widths_match() {
        assert_eq!(map_i(0, 8, 8), 0);
        assert_eq!(map_i(7, 8, 8), 7);
    }

    #[test]
    fn map_i_rounds_to_nearest() {
        assert_eq!(map_i(0, 4, 8), 1);
        assert_eq!(map_i(3, 4, 8), 7);
        assert_eq!(map_i(0, 8, 4), 0);
        assert_eq!(map_i(7, 8, 4), 3);
    }
}
