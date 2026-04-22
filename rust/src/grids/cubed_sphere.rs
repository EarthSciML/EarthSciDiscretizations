//! Cubed-sphere grid accessor runtime (FV3 gnomonic convention).
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.5, §3.3, §7.
//!
//! Per the 2026-04-20 scope correction, the `.esm` lowering is a small
//! declarative config (family, dimensions, generator reference) — not a
//! serialized geometry blob. Geometry is derived on demand by accessors:
//!
//! * [`CubedSphereGrid::cell_centers`] — cell-center (lon, lat)
//! * [`CubedSphereGrid::neighbors`] — 4-way face connectivity (cross-panel aware)
//! * [`CubedSphereGrid::metric_eval`] — metric-tensor components at cell centers
//!
//! Cross-binding conformance (`docs/GRIDS_API.md` §4) compares accessor outputs
//! at pinned query points rather than serialized bytes.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use serde_json::{json, Value};

use crate::{Dtype, Grid, GridError, Result};

/// Number of panels on a cubed-sphere grid.
pub const N_PANELS: usize = 6;

/// Cardinal directions used in accessor return keys.
pub const DIRECTIONS: [&str; 4] = ["W", "E", "S", "N"];

/// `(neighbor_panel, neighbor_edge, reverse_along_edge)` for each (panel, edge).
///
/// Mirrors `python/.../cubed_sphere.py::_PANEL_CONNECTIVITY` and Julia's
/// `src/grids/panel_connectivity.jl` (0-indexed here).
const PANEL_CONNECTIVITY: [[(usize, Edge, bool); 4]; N_PANELS] = [
    // panel 0
    [
        /* W */ (4, Edge::E, false),
        /* E */ (1, Edge::W, false),
        /* S */ (5, Edge::N, false),
        /* N */ (2, Edge::S, false),
    ],
    // panel 1
    [
        (0, Edge::E, false),
        (3, Edge::W, false),
        (5, Edge::E, true),
        (2, Edge::E, false),
    ],
    // panel 2
    [
        (4, Edge::N, true),
        (1, Edge::N, false),
        (0, Edge::N, false),
        (3, Edge::N, true),
    ],
    // panel 3
    [
        (1, Edge::E, false),
        (4, Edge::W, false),
        (5, Edge::S, true),
        (2, Edge::N, true),
    ],
    // panel 4
    [
        (3, Edge::E, false),
        (0, Edge::W, false),
        (5, Edge::W, false),
        (2, Edge::W, true),
    ],
    // panel 5
    [
        (4, Edge::S, false),
        (1, Edge::S, true),
        (3, Edge::S, true),
        (0, Edge::S, false),
    ],
];

/// Panel edge identifier. Used in cross-panel neighbor lookup.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Edge {
    W,
    E,
    S,
    N,
}

/// Cardinal direction. Indices into the `neighbors(...)` return.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Direction {
    W,
    E,
    S,
    N,
}

impl Direction {
    fn idx(self) -> usize {
        match self {
            Direction::W => 0,
            Direction::E => 1,
            Direction::S => 2,
            Direction::N => 3,
        }
    }
}

/// Metric field names accepted by [`CubedSphereGrid::metric_eval`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MetricName {
    /// Jacobian determinant of the (xi, eta) → sphere-surface map.
    J,
    /// Covariant metric component g_{ξξ}.
    GXiXi,
    /// Covariant metric component g_{ηη}.
    GEtaEta,
    /// Covariant metric component g_{ξη}.
    GXiEta,
    /// Inverse-metric component g^{ξξ}.
    GInvXiXi,
    /// Inverse-metric component g^{ηη}.
    GInvEtaEta,
    /// Inverse-metric component g^{ξη}.
    GInvXiEta,
    /// Spherical-quad face area of the cell (uses panel geometry at corners).
    Area,
}

impl MetricName {
    /// Parse the canonical wire-form name (e.g. `"g_xixi"`).
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "J" => MetricName::J,
            "g_xixi" => MetricName::GXiXi,
            "g_etaeta" => MetricName::GEtaEta,
            "g_xieta" => MetricName::GXiEta,
            "ginv_xixi" => MetricName::GInvXiXi,
            "ginv_etaeta" => MetricName::GInvEtaEta,
            "ginv_xieta" => MetricName::GInvXiEta,
            "area" => MetricName::Area,
            _ => return None,
        })
    }
}

// --- gnomonic math ---------------------------------------------------------

/// Unit-sphere Cartesian `(x, y, z)` for computational coords `(xi, eta)` on `panel`.
#[inline]
fn gnomonic_to_cart(xi: f64, eta: f64, panel: usize) -> (f64, f64, f64) {
    let x_ = xi.tan();
    let y_ = eta.tan();
    let r = (1.0 + x_ * x_ + y_ * y_).sqrt();
    match panel {
        0 => (1.0 / r, x_ / r, y_ / r),
        1 => (-x_ / r, 1.0 / r, y_ / r),
        2 => (-y_ / r, x_ / r, 1.0 / r),
        3 => (-1.0 / r, -x_ / r, y_ / r),
        4 => (x_ / r, -1.0 / r, y_ / r),
        5 => (y_ / r, x_ / r, -1.0 / r),
        _ => unreachable!("panel index validated by caller"),
    }
}

/// Gnomonic (lon, lat) at (xi, eta) on `panel`.
#[inline]
fn gnomonic_to_lonlat(xi: f64, eta: f64, panel: usize) -> (f64, f64) {
    let (x, y, z) = gnomonic_to_cart(xi, eta, panel);
    let lon = y.atan2(x);
    let lat = z.clamp(-1.0, 1.0).asin();
    (lon, lat)
}

/// Covariant metric `(J, g_xixi, g_etaeta, g_xieta)` at (xi, eta) for radius R.
///
/// Gnomonic metric is panel-independent: each panel is a rigid rotation of the
/// sphere, so the (xi, eta) → sphere metric is identical across panels.
#[inline]
fn gnomonic_metric(xi: f64, eta: f64, r: f64) -> (f64, f64, f64, f64) {
    let x_ = xi.tan();
    let y_ = eta.tan();
    let d2 = 1.0 + x_ * x_ + y_ * y_;
    let sx = 1.0 + x_ * x_;
    let sy = 1.0 + y_ * y_;
    let d3 = d2 * d2.sqrt();
    let d4 = d2 * d2;
    let r2 = r * r;
    let j = r2 * sx * sy / d3;
    let gxx = r2 * sx * sx * sy / d4;
    let gyy = r2 * sy * sy * sx / d4;
    let gxy = -r2 * x_ * y_ * sx * sy / d4;
    (j, gxx, gyy, gxy)
}

/// Area of a spherical quadrilateral given 4 unit-sphere corners and radius.
///
/// Uses the planar-projection interior-angle sum (mirrors
/// `python/.../cubed_sphere.py::_spherical_quad_area`). Stable for the small
/// non-planar quads at cube corners.
fn spherical_quad_area(corners: &[(f64, f64, f64); 4], r: f64) -> f64 {
    let mut total = 0.0_f64;
    for k in 0..4 {
        let prev = corners[(k + 3) % 4];
        let curr = corners[k];
        let next = corners[(k + 1) % 4];
        let t1 = cross(curr, cross(sub(prev, curr), curr));
        let t2 = cross(curr, cross(sub(next, curr), curr));
        let n1 = norm(t1);
        let n2 = norm(t2);
        total += if n1 < 1e-15 || n2 < 1e-15 {
            FRAC_PI_2
        } else {
            (dot(t1, t2) / (n1 * n2)).clamp(-1.0, 1.0).acos()
        };
    }
    r * r * (total - 2.0 * PI)
}

#[inline]
fn cross(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

#[inline]
fn sub(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

#[inline]
fn dot(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

#[inline]
fn norm(a: (f64, f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1 + a.2 * a.2).sqrt()
}

// --- public grid object ----------------------------------------------------

/// Gnomonic cubed-sphere grid: 6 panels × `nc` × `nc` cells.
///
/// See `docs/GRIDS_API.md` §2.5, §3.3, §7. The struct stores only the
/// declarative parameters plus pre-tabulated 1D (xi, eta) axis arrays; all
/// geometric quantities are derived on demand from the `gnomonic_c6` generator.
#[derive(Clone, Debug)]
pub struct CubedSphereGrid {
    nc: usize,
    r: f64,
    dtype: Dtype,
    ghosts: usize,
    xi_edges: Vec<f64>,
    eta_edges: Vec<f64>,
    xi_centers: Vec<f64>,
    eta_centers: Vec<f64>,
}

impl CubedSphereGrid {
    /// Cells per panel edge.
    #[inline]
    pub fn nc(&self) -> usize {
        self.nc
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

    /// Topology string (always `"block_structured"` for this family).
    #[inline]
    pub fn topology(&self) -> &'static str {
        "block_structured"
    }

    /// Total cell count, `6 * nc * nc`.
    #[inline]
    pub fn n_cells(&self) -> usize {
        N_PANELS * self.nc * self.nc
    }

    /// Read-only view of the xi cell-edge coordinates (length `nc + 1`).
    #[inline]
    pub fn xi_edges(&self) -> &[f64] {
        &self.xi_edges
    }

    /// Read-only view of the eta cell-edge coordinates (length `nc + 1`).
    #[inline]
    pub fn eta_edges(&self) -> &[f64] {
        &self.eta_edges
    }

    /// Read-only view of the xi cell-center coordinates (length `nc`).
    #[inline]
    pub fn xi_centers(&self) -> &[f64] {
        &self.xi_centers
    }

    /// Read-only view of the eta cell-center coordinates (length `nc`).
    #[inline]
    pub fn eta_centers(&self) -> &[f64] {
        &self.eta_centers
    }

    /// Cell-center `(lon, lat)` for cell `(p, i, j)`.
    pub fn cell_centers(&self, p: usize, i: usize, j: usize) -> Result<(f64, f64)> {
        self.check_cell(p, i, j)?;
        Ok(gnomonic_to_lonlat(
            self.xi_centers[i],
            self.eta_centers[j],
            p,
        ))
    }

    /// Bulk materialization of cell-center `(lon, lat)` arrays of shape `(6, nc, nc)`.
    ///
    /// Each array is flat row-major with stride `nc*nc` per panel and stride `nc`
    /// per `i` index, matching the Python `(6, Nc, Nc)` layout.
    pub fn cell_centers_bulk(&self) -> (Vec<f64>, Vec<f64>) {
        let nc = self.nc;
        let mut lon = vec![0.0; N_PANELS * nc * nc];
        let mut lat = vec![0.0; N_PANELS * nc * nc];
        for p in 0..N_PANELS {
            for i in 0..nc {
                for j in 0..nc {
                    let (lo, la) = gnomonic_to_lonlat(self.xi_centers[i], self.eta_centers[j], p);
                    let idx = p * nc * nc + i * nc + j;
                    lon[idx] = lo;
                    lat[idx] = la;
                }
            }
        }
        (lon, lat)
    }

    /// Face neighbors of cell `(p, i, j)` as `[W, E, S, N] => (p', i', j')`.
    ///
    /// Interior cells resolve within the same panel; boundary cells cross to the
    /// neighbor panel via [`PANEL_CONNECTIVITY`], with the along-edge index
    /// reversed when the edges are oriented opposite.
    pub fn neighbors(&self, p: usize, i: usize, j: usize) -> Result<[(usize, usize, usize); 4]> {
        self.check_cell(p, i, j)?;
        let nc = self.nc;
        let w = if i > 0 {
            (p, i - 1, j)
        } else {
            self.cross_panel_neighbor(p, Direction::W, i, j)
        };
        let e = if i + 1 < nc {
            (p, i + 1, j)
        } else {
            self.cross_panel_neighbor(p, Direction::E, i, j)
        };
        let s = if j > 0 {
            (p, i, j - 1)
        } else {
            self.cross_panel_neighbor(p, Direction::S, i, j)
        };
        let n = if j + 1 < nc {
            (p, i, j + 1)
        } else {
            self.cross_panel_neighbor(p, Direction::N, i, j)
        };
        Ok([w, e, s, n])
    }

    /// Return a single neighbor in direction `dir`.
    pub fn neighbor(
        &self,
        p: usize,
        i: usize,
        j: usize,
        dir: Direction,
    ) -> Result<(usize, usize, usize)> {
        Ok(self.neighbors(p, i, j)?[dir.idx()])
    }

    /// Evaluate a metric field at cell center `(p, i, j)`.
    pub fn metric_eval(&self, name: MetricName, p: usize, i: usize, j: usize) -> Result<f64> {
        self.check_cell(p, i, j)?;
        if name == MetricName::Area {
            return Ok(self.cell_area(p, i, j));
        }
        let xi = self.xi_centers[i];
        let eta = self.eta_centers[j];
        let (jac, gxx, gyy, gxy) = gnomonic_metric(xi, eta, self.r);
        Ok(match name {
            MetricName::J => jac,
            MetricName::GXiXi => gxx,
            MetricName::GEtaEta => gyy,
            MetricName::GXiEta => gxy,
            MetricName::GInvXiXi | MetricName::GInvEtaEta | MetricName::GInvXiEta => {
                let det = gxx * gyy - gxy * gxy;
                match name {
                    MetricName::GInvXiXi => gyy / det,
                    MetricName::GInvEtaEta => gxx / det,
                    MetricName::GInvXiEta => -gxy / det,
                    _ => unreachable!(),
                }
            }
            MetricName::Area => unreachable!("handled above"),
        })
    }

    /// Evaluate a metric by its canonical wire-form name (`"g_xixi"` etc.).
    pub fn metric_eval_by_name(&self, name: &str, p: usize, i: usize, j: usize) -> Result<f64> {
        let m = MetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, p, i, j)
    }

    /// Cell face area for every cell as a flat row-major `(6, nc, nc)` array.
    pub fn area_bulk(&self) -> Vec<f64> {
        let nc = self.nc;
        let mut out = vec![0.0; N_PANELS * nc * nc];
        for p in 0..N_PANELS {
            for i in 0..nc {
                for j in 0..nc {
                    out[p * nc * nc + i * nc + j] = self.cell_area(p, i, j);
                }
            }
        }
        out
    }

    /// Provenance block per `docs/GRIDS_API.md` §6.4.
    pub fn provenance(&self) -> Value {
        json!({
            "binding": "rust",
            "binding_version": env!("CARGO_PKG_VERSION"),
            "source": "earthsci_grids::grids::cubed_sphere",
            "generator": "gnomonic_c6",
        })
    }

    // internals -----------------------------------------------------------

    fn check_cell(&self, p: usize, i: usize, j: usize) -> Result<()> {
        if p >= N_PANELS {
            return Err(GridError::InvalidOption(
                "p",
                format!("panel out of range [0, 6): {p}"),
            ));
        }
        if i >= self.nc {
            return Err(GridError::InvalidOption(
                "i",
                format!("i out of range [0, {}): {i}", self.nc),
            ));
        }
        if j >= self.nc {
            return Err(GridError::InvalidOption(
                "j",
                format!("j out of range [0, {}): {j}", self.nc),
            ));
        }
        Ok(())
    }

    fn cross_panel_neighbor(
        &self,
        p: usize,
        dir: Direction,
        i: usize,
        j: usize,
    ) -> (usize, usize, usize) {
        let (nb_panel, nb_edge, reverse) = PANEL_CONNECTIVITY[p][dir.idx()];
        let mut a = match dir {
            Direction::W | Direction::E => j,
            Direction::S | Direction::N => i,
        };
        if reverse {
            a = self.nc - 1 - a;
        }
        let last = self.nc - 1;
        match nb_edge {
            Edge::W => (nb_panel, 0, a),
            Edge::E => (nb_panel, last, a),
            Edge::S => (nb_panel, a, 0),
            Edge::N => (nb_panel, a, last),
        }
    }

    fn cell_area(&self, p: usize, i: usize, j: usize) -> f64 {
        let xi_w = self.xi_edges[i];
        let xi_e = self.xi_edges[i + 1];
        let eta_s = self.eta_edges[j];
        let eta_n = self.eta_edges[j + 1];
        let corners = [
            gnomonic_to_cart(xi_w, eta_s, p),
            gnomonic_to_cart(xi_e, eta_s, p),
            gnomonic_to_cart(xi_e, eta_n, p),
            gnomonic_to_cart(xi_w, eta_n, p),
        ];
        spherical_quad_area(&corners, self.r)
    }
}

impl Grid for CubedSphereGrid {
    fn family(&self) -> &'static str {
        "cubed_sphere"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative `.esm` lowering per the 2026-04-20 scope correction.
    ///
    /// Emits a §6-schema-shaped config — family, dtype, topology, generator
    /// reference, and params — rather than inline geometry arrays.
    fn to_esm(&self) -> Value {
        json!({
            "family": self.family(),
            "version": "1.0.0",
            "dtype": dtype_str(self.dtype),
            "topology": self.topology(),
            "generator": "gnomonic_c6",
            "params": {
                "Nc": self.nc,
                "R": self.r,
                "ghosts": self.ghosts,
            },
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

/// Create a new [`CubedSphereBuilder`] with default options.
///
/// ```
/// use earthsci_grids::{cubed_sphere, Dtype, Grid};
/// let g = cubed_sphere::builder().nc(48).r(6.371e6).build().unwrap();
/// assert_eq!(g.family(), "cubed_sphere");
/// assert_eq!(g.n_cells(), 6 * 48 * 48);
/// assert_eq!(g.dtype(), Dtype::F64);
/// ```
pub fn builder() -> CubedSphereBuilder {
    CubedSphereBuilder::default()
}

/// Builder for [`CubedSphereGrid`]. Required options are validated at
/// [`CubedSphereBuilder::build`] time so partial builders can be reused in tests.
#[must_use]
#[derive(Clone, Debug, Default)]
pub struct CubedSphereBuilder {
    nc: Option<usize>,
    r: Option<f64>,
    dtype: Option<Dtype>,
    ghosts: Option<usize>,
}

impl CubedSphereBuilder {
    /// Cells per panel edge. Required. Must be ≥ 1.
    pub fn nc(mut self, nc: usize) -> Self {
        self.nc = Some(nc);
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

    /// Validate options and construct the grid.
    pub fn build(self) -> Result<CubedSphereGrid> {
        let nc = self.nc.ok_or(GridError::MissingOption("Nc"))?;
        if nc < 1 {
            return Err(GridError::InvalidOption(
                "Nc",
                format!("must be >= 1, got {nc}"),
            ));
        }

        let r = self.r.unwrap_or(6.371e6);
        if !(r > 0.0 && r.is_finite()) {
            return Err(GridError::InvalidOption(
                "R",
                format!("must be a positive finite number, got {r}"),
            ));
        }

        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);

        let dxi = FRAC_PI_2 / (nc as f64);
        let mut xi_edges = Vec::with_capacity(nc + 1);
        for k in 0..=nc {
            xi_edges.push(-FRAC_PI_4 + (k as f64) * dxi);
        }
        let eta_edges = xi_edges.clone();
        let mut xi_centers = Vec::with_capacity(nc);
        for k in 0..nc {
            xi_centers.push(0.5 * (xi_edges[k] + xi_edges[k + 1]));
        }
        let eta_centers = xi_centers.clone();

        Ok(CubedSphereGrid {
            nc,
            r,
            dtype,
            ghosts,
            xi_edges,
            eta_edges,
            xi_centers,
            eta_centers,
        })
    }
}

// --- tests -----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn g(nc: usize) -> CubedSphereGrid {
        builder().nc(nc).build().unwrap()
    }

    // --- API ---

    #[test]
    fn builder_requires_nc() {
        let err = builder().build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("Nc")));
    }

    #[test]
    fn builder_rejects_zero_nc() {
        let err = builder().nc(0).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("Nc", _)));
    }

    #[test]
    fn builder_rejects_nonfinite_r() {
        let err = builder().nc(4).r(f64::NAN).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));
        let err = builder().nc(4).r(-1.0).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));
    }

    #[test]
    fn defaults_match_contract() {
        let g = builder().nc(4).build().unwrap();
        assert_eq!(g.nc(), 4);
        assert_eq!(g.r(), 6.371e6);
        assert_eq!(g.ghosts(), 0);
        assert_eq!(g.dtype(), Dtype::F64);
        assert_eq!(g.family(), "cubed_sphere");
        assert_eq!(g.topology(), "block_structured");
        assert_eq!(g.n_cells(), 6 * 4 * 4);
    }

    // --- topology ---

    #[test]
    fn panel_connectivity_symmetric() {
        for (p, edges) in PANEL_CONNECTIVITY.iter().enumerate() {
            for (d, &(nb_panel, nb_edge, _)) in edges.iter().enumerate() {
                let back = match nb_edge {
                    Edge::W => PANEL_CONNECTIVITY[nb_panel][Direction::W.idx()],
                    Edge::E => PANEL_CONNECTIVITY[nb_panel][Direction::E.idx()],
                    Edge::S => PANEL_CONNECTIVITY[nb_panel][Direction::S.idx()],
                    Edge::N => PANEL_CONNECTIVITY[nb_panel][Direction::N.idx()],
                };
                assert_eq!(back.0, p, "panel {p} dir {d} does not round-trip");
            }
        }
    }

    #[test]
    fn interior_neighbors_local() {
        let gr = g(8);
        let ns = gr.neighbors(0, 3, 4).unwrap();
        assert_eq!(ns[Direction::W.idx()], (0, 2, 4));
        assert_eq!(ns[Direction::E.idx()], (0, 4, 4));
        assert_eq!(ns[Direction::S.idx()], (0, 3, 3));
        assert_eq!(ns[Direction::N.idx()], (0, 3, 5));
    }

    #[test]
    fn neighbor_cells_in_range() {
        let gr = g(5);
        for p in 0..N_PANELS {
            for i in 0..gr.nc() {
                for j in 0..gr.nc() {
                    for (pp, ii, jj) in gr.neighbors(p, i, j).unwrap() {
                        assert!(pp < N_PANELS);
                        assert!(ii < gr.nc());
                        assert!(jj < gr.nc());
                    }
                }
            }
        }
    }

    #[test]
    fn cross_panel_neighbor_returns_to_source_panel() {
        let gr = g(6);
        for (p, panel_edges) in PANEL_CONNECTIVITY.iter().enumerate() {
            for &d in &[Direction::W, Direction::E, Direction::S, Direction::N] {
                let (i, j) = match d {
                    Direction::W => (0, 2),
                    Direction::E => (gr.nc() - 1, 2),
                    Direction::S => (2, 0),
                    Direction::N => (2, gr.nc() - 1),
                };
                let (nbp, nbi, nbj) = gr.neighbor(p, i, j, d).unwrap();
                // stepping the same direction from the neighbor back across the
                // panel boundary should land on the source panel.
                let back_dir = match panel_edges[d.idx()].1 {
                    Edge::W => Direction::W,
                    Edge::E => Direction::E,
                    Edge::S => Direction::S,
                    Edge::N => Direction::N,
                };
                let (back_p, _, _) = gr.neighbor(nbp, nbi, nbj, back_dir).unwrap();
                assert_eq!(back_p, p, "{p}/{d:?} did not return to {p}");
            }
        }
    }

    // --- centers ---

    #[test]
    fn gnomonic_origin_maps_to_cardinals() {
        let (lon0, lat0) = gnomonic_to_lonlat(0.0, 0.0, 0);
        assert_relative_eq!(lon0, 0.0, epsilon = 1e-15);
        assert_relative_eq!(lat0, 0.0, epsilon = 1e-15);
        let (_, lat2) = gnomonic_to_lonlat(0.0, 0.0, 2);
        assert_relative_eq!(lat2, FRAC_PI_2, epsilon = 1e-15);
        let (_, lat5) = gnomonic_to_lonlat(0.0, 0.0, 5);
        assert_relative_eq!(lat5, -FRAC_PI_2, epsilon = 1e-15);
    }

    #[test]
    fn cart_on_unit_sphere() {
        for &xi in &[-0.3_f64, 0.0, 0.2] {
            for &eta in &[-0.25_f64, 0.0, 0.35] {
                for p in 0..6 {
                    let (x, y, z) = gnomonic_to_cart(xi, eta, p);
                    assert!((x * x + y * y + z * z - 1.0).abs() < 1e-14);
                }
            }
        }
    }

    #[test]
    fn cell_centers_bulk_matches_scalar() {
        let gr = g(4);
        let (lon, lat) = gr.cell_centers_bulk();
        for p in [0usize, 2, 5] {
            for &(i, j) in &[(0_usize, 0_usize), (1, 2), (3, 3)] {
                let (s_lon, s_lat) = gr.cell_centers(p, i, j).unwrap();
                let idx = p * 4 * 4 + i * 4 + j;
                assert_relative_eq!(s_lon, lon[idx], epsilon = 1e-12);
                assert_relative_eq!(s_lat, lat[idx], epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn cell_centers_range() {
        let gr = g(8);
        let (lon, lat) = gr.cell_centers_bulk();
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
        let gr = builder().nc(24).r(r_val).build().unwrap();
        let total: f64 = gr.area_bulk().iter().sum();
        let expected = 4.0 * PI * r_val * r_val;
        assert!((total - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn total_unit_sphere_area_is_4pi() {
        let gr = builder().nc(16).r(1.0).build().unwrap();
        let total: f64 = gr.area_bulk().iter().sum();
        assert!((total - 4.0 * PI).abs() < 1e-9);
    }

    #[test]
    fn area_positive_everywhere() {
        let gr = g(8);
        for a in gr.area_bulk() {
            assert!(a > 0.0);
        }
    }

    #[test]
    fn area_via_metric_matches_bulk() {
        let gr = g(4);
        let bulk = gr.area_bulk();
        for p in [0_usize, 1, 3] {
            for i in [0_usize, 2, 3] {
                for j in [0_usize, 1, 3] {
                    let m = gr.metric_eval(MetricName::Area, p, i, j).unwrap();
                    let idx = p * 4 * 4 + i * 4 + j;
                    assert_relative_eq!(m, bulk[idx], max_relative = 1e-12);
                }
            }
        }
    }

    // --- metric ---

    #[test]
    fn metric_at_origin_is_r_squared_identity() {
        let r = 3.0_f64;
        let (j, gxx, gyy, gxy) = gnomonic_metric(0.0, 0.0, r);
        assert_relative_eq!(j, r * r, epsilon = 1e-14);
        assert_relative_eq!(gxx, r * r, epsilon = 1e-14);
        assert_relative_eq!(gyy, r * r, epsilon = 1e-14);
        assert!(gxy.abs() < 1e-14);
    }

    #[test]
    fn metric_jacobian_equals_sqrt_det() {
        let r = 2.5_f64;
        for &(xi, eta) in &[(0.0, 0.0), (0.1, -0.2), (-0.3, 0.25), (0.5, 0.4)] {
            let (j, gxx, gyy, gxy) = gnomonic_metric(xi, eta, r);
            let det = gxx * gyy - gxy * gxy;
            assert!(det > 0.0);
            assert_relative_eq!(j, det.sqrt(), max_relative = 1e-12);
        }
    }

    #[test]
    fn inverse_metric_is_inverse() {
        let gr = builder().nc(6).r(1.0).build().unwrap();
        for &(i, j) in &[(0_usize, 0_usize), (2, 3), (5, 5)] {
            let gxx = gr.metric_eval(MetricName::GXiXi, 0, i, j).unwrap();
            let gyy = gr.metric_eval(MetricName::GEtaEta, 0, i, j).unwrap();
            let gxy = gr.metric_eval(MetricName::GXiEta, 0, i, j).unwrap();
            let ixx = gr.metric_eval(MetricName::GInvXiXi, 0, i, j).unwrap();
            let iyy = gr.metric_eval(MetricName::GInvEtaEta, 0, i, j).unwrap();
            let ixy = gr.metric_eval(MetricName::GInvXiEta, 0, i, j).unwrap();
            let m00 = gxx * ixx + gxy * ixy;
            let m01 = gxx * ixy + gxy * iyy;
            let m10 = gxy * ixx + gyy * ixy;
            let m11 = gxy * ixy + gyy * iyy;
            assert_relative_eq!(m00, 1.0, epsilon = 1e-12);
            assert_relative_eq!(m11, 1.0, epsilon = 1e-12);
            assert!(m01.abs() < 1e-12);
            assert!(m10.abs() < 1e-12);
        }
    }

    #[test]
    fn metric_eval_rejects_out_of_range_cell() {
        let gr = g(4);
        assert!(gr.metric_eval(MetricName::J, 6, 0, 0).is_err());
        assert!(gr.metric_eval(MetricName::J, 0, 4, 0).is_err());
    }

    #[test]
    fn metric_eval_by_name_rejects_unknown() {
        let gr = g(4);
        let err = gr.metric_eval_by_name("not_a_metric", 0, 0, 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn j_positive_everywhere() {
        let gr = g(6);
        for p in 0..N_PANELS {
            for i in 0..gr.nc() {
                for j in 0..gr.nc() {
                    assert!(gr.metric_eval(MetricName::J, p, i, j).unwrap() > 0.0);
                }
            }
        }
    }

    // --- to_esm ---

    #[test]
    fn to_esm_is_declarative() {
        let gr = builder().nc(48).r(6.371e6).ghosts(3).build().unwrap();
        let doc = gr.to_esm();
        assert_eq!(doc["family"], "cubed_sphere");
        assert_eq!(doc["version"], "1.0.0");
        assert_eq!(doc["dtype"], "float64");
        assert_eq!(doc["topology"], "block_structured");
        assert_eq!(doc["generator"], "gnomonic_c6");
        assert_eq!(doc["params"]["Nc"], 48);
        assert_eq!(doc["params"]["R"], 6.371e6);
        assert_eq!(doc["params"]["ghosts"], 3);
        assert!(doc.get("provenance").is_some());

        // No inline geometry arrays anywhere in the serialized config.
        let wire = serde_json::to_string(&doc).unwrap();
        assert!(!wire.contains("cells"));
        assert!(!wire.contains("lon_array"));
        // And the wire form must stay small (declarative, not a blob).
        assert!(wire.len() < 2_000, "wire too large: {} bytes", wire.len());
    }

    #[test]
    fn to_esm_roundtrips_via_json() {
        let gr = builder().nc(4).r(1.0).build().unwrap();
        let doc = gr.to_esm();
        let text = serde_json::to_string(&doc).unwrap();
        let reparsed: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(reparsed["params"]["Nc"], 4);
        assert_eq!(reparsed["params"]["R"], 1.0);
    }

    #[test]
    fn provenance_identifies_rust_binding() {
        let gr = g(4);
        let doc = gr.to_esm();
        let prov = &doc["provenance"];
        assert_eq!(prov["binding"], "rust");
        assert_eq!(prov["binding_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(prov["generator"], "gnomonic_c6");
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let gr = builder().nc(4).dtype(Dtype::F32).build().unwrap();
        assert_eq!(gr.dtype(), Dtype::F32);
        assert_eq!(gr.to_esm()["dtype"], "float32");
    }
}
