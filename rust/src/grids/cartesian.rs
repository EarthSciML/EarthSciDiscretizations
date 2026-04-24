//! Cartesian grid family (1D / 2D / 3D, uniform + non-uniform) — accessor runtime.
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.5, §3.3, §7.
//!
//! Per the 2026-04-20 scope correction (mayor), the `.esm` lowering is a small
//! declarative config (family + dimensions + extents + optional non-uniform
//! edges), NOT a serialized geometry blob. Geometry is derived on demand via:
//!
//! * [`CartesianGrid::cell_center`] — cell-center coordinates
//! * [`CartesianGrid::cell_widths`] — per-axis widths at a cell
//! * [`CartesianGrid::cell_volume`] — length (1D), area (2D), or volume (3D)
//! * [`CartesianGrid::neighbors`] — axis-aligned face neighbors
//! * [`CartesianGrid::metric_eval`] — scalar metric fields at a cell
//!
//! Mirrors the Julia reference binding in `src/grids/cartesian.jl`. Indices
//! are 0-based in Rust (Julia is 1-based); semantics and accessor outputs are
//! otherwise identical.
//!
//! Two construction modes are supported through [`CartesianBuilder`]:
//! * **Uniform** — `nx` (+ `ny`, `nz`) with `extent` per axis; edges and
//!   centers are derived by pure math.
//! * **Non-uniform** — per-axis `edges` arrays (strictly increasing, ≥ 2
//!   entries) fully specify the grid. Uniform and non-uniform modes are
//!   mutually exclusive.

use serde_json::{json, Value};

use crate::{Dtype, Grid, GridError, Result};

const API_VERSION: &str = "1.0.0";

/// Axis-aligned face-neighbor reference.
///
/// `axis` is a 0-based axis index (`0 = x`, `1 = y`, `2 = z`). `side` is
/// `-1` (low face) or `+1` (high face). `index` is the neighbor cell's
/// multi-index with the same arity as the grid's [`CartesianGrid::ndim`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CartesianNeighbor {
    pub axis: usize,
    pub side: i8,
    pub index: Vec<usize>,
}

/// Scalar metric fields on a cartesian grid. See [`CartesianGrid::metric_eval`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MetricName {
    /// Cell measure: length (1D), area (2D), volume (3D).
    Volume,
    /// Metric-tensor determinant square root. `1` for cartesian.
    Jacobian,
    /// Width along axis 0. Valid for any `ndim >= 1`.
    Dx,
    /// Width along axis 1. Requires `ndim >= 2`.
    Dy,
    /// Width along axis 2. Requires `ndim == 3`.
    Dz,
    /// Face area normal to axis 0. `1` (1D), `dy` (2D), `dy*dz` (3D).
    FaceAreaX,
    /// Face area normal to axis 1. Requires `ndim >= 2`.
    FaceAreaY,
    /// Face area normal to axis 2. Requires `ndim == 3`.
    FaceAreaZ,
}

impl MetricName {
    /// Parse the canonical wire-form name.
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "volume" => MetricName::Volume,
            "jacobian" => MetricName::Jacobian,
            "dx" => MetricName::Dx,
            "dy" => MetricName::Dy,
            "dz" => MetricName::Dz,
            "face_area_x" => MetricName::FaceAreaX,
            "face_area_y" => MetricName::FaceAreaY,
            "face_area_z" => MetricName::FaceAreaZ,
            _ => return None,
        })
    }
}

/// Uniform or non-uniform cartesian grid in 1D / 2D / 3D.
///
/// The struct stores only the declarative parameters plus the pre-tabulated
/// 1-D edge / center / width arrays per axis; per-cell geometric quantities
/// are derived on demand.
#[derive(Clone, Debug)]
pub struct CartesianGrid {
    ndim: usize,
    n: Vec<usize>,
    extent: Vec<(f64, f64)>,
    edges: Vec<Vec<f64>>,
    centers: Vec<Vec<f64>>,
    widths: Vec<Vec<f64>>,
    uniform: Vec<bool>,
    ghosts: usize,
    dtype: Dtype,
}

impl CartesianGrid {
    /// Number of axes (1, 2, or 3).
    #[inline]
    pub fn ndim(&self) -> usize {
        self.ndim
    }

    /// Per-axis cell counts (length `ndim`).
    #[inline]
    pub fn n(&self) -> &[usize] {
        &self.n
    }

    /// Per-axis `(lo, hi)` extents (length `ndim`).
    #[inline]
    pub fn extent(&self) -> &[(f64, f64)] {
        &self.extent
    }

    /// Per-axis edges. `edges()[d]` has length `n()[d] + 1`.
    #[inline]
    pub fn edges(&self) -> &[Vec<f64>] {
        &self.edges
    }

    /// Per-axis cell centers. `centers()[d]` has length `n()[d]`.
    #[inline]
    pub fn centers(&self) -> &[Vec<f64>] {
        &self.centers
    }

    /// Per-axis cell widths. `widths()[d]` has length `n()[d]`.
    #[inline]
    pub fn widths(&self) -> &[Vec<f64>] {
        &self.widths
    }

    /// Per-axis uniform flags. `uniform()[d]` is `true` iff all widths on
    /// axis `d` are equal to within `8 * eps * |w0|` tolerance.
    #[inline]
    pub fn uniform(&self) -> &[bool] {
        &self.uniform
    }

    /// Halo cell width.
    #[inline]
    pub fn ghosts(&self) -> usize {
        self.ghosts
    }

    /// Topology string (always `"rectilinear"` for this family per §7).
    #[inline]
    pub fn topology(&self) -> &'static str {
        "rectilinear"
    }

    /// Total cell count across all axes.
    pub fn n_cells(&self) -> usize {
        self.n.iter().product()
    }

    /// Cell-center coordinates for cell `idx` as a `Vec` of length `ndim`.
    pub fn cell_center(&self, idx: &[usize]) -> Result<Vec<f64>> {
        self.check_cell(idx)?;
        Ok((0..self.ndim).map(|d| self.centers[d][idx[d]]).collect())
    }

    /// Per-axis cell widths at cell `idx` as a `Vec` of length `ndim`.
    pub fn cell_widths(&self, idx: &[usize]) -> Result<Vec<f64>> {
        self.check_cell(idx)?;
        Ok((0..self.ndim).map(|d| self.widths[d][idx[d]]).collect())
    }

    /// Cell measure: length (1D), area (2D), or volume (3D).
    pub fn cell_volume(&self, idx: &[usize]) -> Result<f64> {
        self.check_cell(idx)?;
        Ok((0..self.ndim).map(|d| self.widths[d][idx[d]]).product())
    }

    /// Axis-aligned face neighbors of cell `idx`. Boundary faces are omitted.
    ///
    /// Entries come in deterministic order: axis-major, low side before high
    /// side. A cell at a boundary on axis `d` / side `s` produces no entry
    /// for `(d, s)` rather than a sentinel value.
    pub fn neighbors(&self, idx: &[usize]) -> Result<Vec<CartesianNeighbor>> {
        self.check_cell(idx)?;
        let mut out = Vec::with_capacity(2 * self.ndim);
        for d in 0..self.ndim {
            if idx[d] > 0 {
                let mut n_idx = idx.to_vec();
                n_idx[d] -= 1;
                out.push(CartesianNeighbor {
                    axis: d,
                    side: -1,
                    index: n_idx,
                });
            }
            if idx[d] + 1 < self.n[d] {
                let mut n_idx = idx.to_vec();
                n_idx[d] += 1;
                out.push(CartesianNeighbor {
                    axis: d,
                    side: 1,
                    index: n_idx,
                });
            }
        }
        Ok(out)
    }

    /// Evaluate a scalar metric field at cell `idx`.
    pub fn metric_eval(&self, name: MetricName, idx: &[usize]) -> Result<f64> {
        self.check_cell(idx)?;
        match name {
            MetricName::Volume => Ok((0..self.ndim).map(|d| self.widths[d][idx[d]]).product()),
            MetricName::Jacobian => Ok(1.0),
            MetricName::Dx => Ok(self.widths[0][idx[0]]),
            MetricName::Dy => {
                self.require_axis(1, "dy")?;
                Ok(self.widths[1][idx[1]])
            }
            MetricName::Dz => {
                self.require_axis(2, "dz")?;
                Ok(self.widths[2][idx[2]])
            }
            MetricName::FaceAreaX => Ok(self.face_area(0, idx)),
            MetricName::FaceAreaY => {
                self.require_axis(1, "face_area_y")?;
                Ok(self.face_area(1, idx))
            }
            MetricName::FaceAreaZ => {
                self.require_axis(2, "face_area_z")?;
                Ok(self.face_area(2, idx))
            }
        }
    }

    /// Evaluate a scalar metric by its canonical wire-form name.
    pub fn metric_eval_by_name(&self, name: &str, idx: &[usize]) -> Result<f64> {
        let m = MetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, idx)
    }

    /// Metric tensor at cell `idx`. Identity matrix (`ndim x ndim`) for the
    /// cartesian family; `idx` is still validated.
    pub fn metric_g(&self, idx: &[usize]) -> Result<Vec<Vec<f64>>> {
        self.check_cell(idx)?;
        let n = self.ndim;
        Ok((0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect())
    }

    /// Provenance block per `docs/GRIDS_API.md` §6.4.
    pub fn provenance(&self) -> Value {
        json!({
            "binding": "rust",
            "binding_version": env!("CARGO_PKG_VERSION"),
            "api_version": API_VERSION,
            "source": "earthsci_grids::grids::cartesian",
            "generator": "cartesian",
        })
    }

    // internals -----------------------------------------------------------

    fn check_cell(&self, idx: &[usize]) -> Result<()> {
        if idx.len() != self.ndim {
            return Err(GridError::InvalidOption(
                "idx",
                format!(
                    "index arity {} does not match grid ndim {}",
                    idx.len(),
                    self.ndim
                ),
            ));
        }
        for (d, &ix) in idx.iter().enumerate().take(self.ndim) {
            if ix >= self.n[d] {
                return Err(GridError::InvalidOption(
                    "idx",
                    format!("axis {d} index {ix} out of range [0, {})", self.n[d]),
                ));
            }
        }
        Ok(())
    }

    fn require_axis(&self, axis: usize, name: &'static str) -> Result<()> {
        if axis >= self.ndim {
            return Err(GridError::InvalidOption(
                "metric_name",
                format!("metric {name:?} requires ndim > {axis}, got {}", self.ndim),
            ));
        }
        Ok(())
    }

    fn face_area(&self, axis: usize, idx: &[usize]) -> f64 {
        // Product of widths on all axes except `axis`. Empty product for 1D.
        let mut a = 1.0;
        for (d, &ix) in idx.iter().enumerate().take(self.ndim) {
            if d == axis {
                continue;
            }
            a *= self.widths[d][ix];
        }
        a
    }
}

impl Grid for CartesianGrid {
    fn family(&self) -> &'static str {
        "cartesian"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative `.esm` lowering per the 2026-04-20 scope correction.
    ///
    /// Emits `n`, `extent`, and `uniform` for every axis. Non-uniform axes
    /// additionally carry their full edge arrays under `edges` (all-uniform
    /// grids omit the `edges` key to keep the wire form compact).
    fn to_esm(&self) -> Value {
        let extent: Vec<Value> = self.extent.iter().map(|(lo, hi)| json!([lo, hi])).collect();

        let mut out = json!({
            "family": self.family(),
            "version": API_VERSION,
            "dtype": dtype_str(self.dtype),
            "topology": self.topology(),
            "ndim": self.ndim,
            "ghosts": self.ghosts,
            "n_cells": self.n_cells(),
            "n": self.n,
            "extent": extent,
            "uniform": self.uniform,
            "provenance": self.provenance(),
        });

        if self.uniform.iter().any(|&u| !u) {
            let edges: Vec<Value> = (0..self.ndim)
                .map(|d| {
                    if self.uniform[d] {
                        json!([])
                    } else {
                        json!(self.edges[d])
                    }
                })
                .collect();
            out.as_object_mut()
                .expect("to_esm root is an object")
                .insert("edges".to_string(), Value::Array(edges));
        }

        out
    }
}

fn dtype_str(d: Dtype) -> &'static str {
    match d {
        Dtype::F64 => "float64",
        Dtype::F32 => "float32",
    }
}

// --- builder ---------------------------------------------------------------

/// Create a new [`CartesianBuilder`] with default options.
///
/// ```
/// use earthsci_grids::{cartesian, Dtype, Grid};
/// let g = cartesian::builder()
///     .nx(8)
///     .ny(4)
///     .extent(vec![(0.0, 1.0), (0.0, 0.5)])
///     .build()
///     .unwrap();
/// assert_eq!(g.family(), "cartesian");
/// assert_eq!(g.n_cells(), 8 * 4);
/// assert_eq!(g.dtype(), Dtype::F64);
/// ```
pub fn builder() -> CartesianBuilder {
    CartesianBuilder::default()
}

/// Builder for [`CartesianGrid`]. Required options are validated at
/// [`CartesianBuilder::build`] time so partial builders can be reused in tests.
#[must_use]
#[derive(Clone, Debug, Default)]
pub struct CartesianBuilder {
    nx: Option<usize>,
    ny: Option<usize>,
    nz: Option<usize>,
    extent: Option<Vec<(f64, f64)>>,
    edges: Option<Vec<Vec<f64>>>,
    dtype: Option<Dtype>,
    ghosts: Option<usize>,
}

impl CartesianBuilder {
    /// Cell count along axis 0 (x). Required in uniform mode.
    pub fn nx(mut self, nx: usize) -> Self {
        self.nx = Some(nx);
        self
    }

    /// Cell count along axis 1 (y). Presence selects 2D+ in uniform mode.
    pub fn ny(mut self, ny: usize) -> Self {
        self.ny = Some(ny);
        self
    }

    /// Cell count along axis 2 (z). Presence selects 3D in uniform mode.
    pub fn nz(mut self, nz: usize) -> Self {
        self.nz = Some(nz);
        self
    }

    /// Per-axis `(lo, hi)` extents. Length must match the selected ndim.
    /// Required in uniform mode; forbidden in non-uniform mode.
    pub fn extent(mut self, extent: Vec<(f64, f64)>) -> Self {
        self.extent = Some(extent);
        self
    }

    /// Per-axis edge arrays. Each inner array must be strictly increasing
    /// with ≥ 2 entries. Presence selects non-uniform mode; mutually
    /// exclusive with `nx`/`ny`/`nz` and `extent`.
    pub fn edges(mut self, edges: Vec<Vec<f64>>) -> Self {
        self.edges = Some(edges);
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
    pub fn build(self) -> Result<CartesianGrid> {
        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);

        if let Some(edges) = self.edges {
            // Non-uniform path: edges fully specify geometry.
            if self.nx.is_some() || self.ny.is_some() || self.nz.is_some() {
                return Err(GridError::InvalidOption(
                    "edges",
                    "non-uniform `edges` is mutually exclusive with `nx`/`ny`/`nz`".to_string(),
                ));
            }
            if self.extent.is_some() {
                return Err(GridError::InvalidOption(
                    "edges",
                    "non-uniform `edges` is mutually exclusive with `extent`".to_string(),
                ));
            }
            if edges.is_empty() {
                return Err(GridError::InvalidOption(
                    "edges",
                    "must contain at least one axis".to_string(),
                ));
            }
            if edges.len() > 3 {
                return Err(GridError::InvalidOption(
                    "edges",
                    format!("ndim must be ≤ 3; got {}", edges.len()),
                ));
            }
            let ndim = edges.len();
            let mut n = Vec::with_capacity(ndim);
            let mut centers = Vec::with_capacity(ndim);
            let mut widths = Vec::with_capacity(ndim);
            let mut uniform = Vec::with_capacity(ndim);
            let mut extent = Vec::with_capacity(ndim);
            let mut edges_out = Vec::with_capacity(ndim);
            for (d, e_in) in edges.into_iter().enumerate() {
                if e_in.len() < 2 {
                    return Err(GridError::InvalidOption(
                        "edges",
                        format!("axis {d} has {} entries; must be ≥ 2", e_in.len()),
                    ));
                }
                for k in 0..(e_in.len() - 1) {
                    if !e_in[k].is_finite() || !e_in[k + 1].is_finite() {
                        return Err(GridError::InvalidOption(
                            "edges",
                            format!("axis {d} contains non-finite edge"),
                        ));
                    }
                    if e_in[k + 1] <= e_in[k] {
                        return Err(GridError::InvalidOption(
                            "edges",
                            format!("axis {d} edges must be strictly increasing"),
                        ));
                    }
                }
                let n_d = e_in.len() - 1;
                let c_d: Vec<f64> = (0..n_d).map(|k| 0.5 * (e_in[k] + e_in[k + 1])).collect();
                let w_d: Vec<f64> = (0..n_d).map(|k| e_in[k + 1] - e_in[k]).collect();
                let u_d = is_uniform(&w_d);
                n.push(n_d);
                extent.push((e_in[0], e_in[n_d]));
                centers.push(c_d);
                widths.push(w_d);
                uniform.push(u_d);
                edges_out.push(e_in);
            }
            return Ok(CartesianGrid {
                ndim,
                n,
                extent,
                edges: edges_out,
                centers,
                widths,
                uniform,
                ghosts,
                dtype,
            });
        }

        // Uniform path: nx (+ optional ny/nz) + extent.
        let nx = self.nx.ok_or(GridError::MissingOption("nx"))?;
        let ndim = match (self.ny, self.nz) {
            (None, None) => 1,
            (Some(_), None) => 2,
            (Some(_), Some(_)) => 3,
            (None, Some(_)) => {
                return Err(GridError::InvalidOption(
                    "nz",
                    "`nz` requires `ny` to also be set (2D before 3D)".to_string(),
                ));
            }
        };
        let ns = [Some(nx), self.ny, self.nz];
        for (d, opt_n) in ns.iter().take(ndim).enumerate() {
            let v = opt_n.expect("ndim selection guarantees presence");
            if v < 1 {
                let name = axis_count_name(d);
                return Err(GridError::InvalidOption(
                    name,
                    format!("must be ≥ 1, got {v}"),
                ));
            }
        }
        let extent = self.extent.ok_or(GridError::MissingOption("extent"))?;
        if extent.len() != ndim {
            return Err(GridError::InvalidOption(
                "extent",
                format!("length {} does not match ndim {}", extent.len(), ndim),
            ));
        }

        let mut n = Vec::with_capacity(ndim);
        let mut edges_out = Vec::with_capacity(ndim);
        let mut centers = Vec::with_capacity(ndim);
        let mut widths = Vec::with_capacity(ndim);
        for d in 0..ndim {
            let (lo, hi) = extent[d];
            if !lo.is_finite() || !hi.is_finite() {
                return Err(GridError::InvalidOption(
                    "extent",
                    format!("axis {d} extent must be finite"),
                ));
            }
            if hi <= lo {
                return Err(GridError::InvalidOption(
                    "extent",
                    format!("axis {d} must satisfy hi > lo; got ({lo}, {hi})"),
                ));
            }
            let n_d = ns[d].expect("ndim guarantees presence");
            let dx = (hi - lo) / n_d as f64;
            let e_d: Vec<f64> = (0..=n_d).map(|k| lo + k as f64 * dx).collect();
            let c_d: Vec<f64> = (0..n_d).map(|k| 0.5 * (e_d[k] + e_d[k + 1])).collect();
            let w_d: Vec<f64> = (0..n_d).map(|k| e_d[k + 1] - e_d[k]).collect();
            n.push(n_d);
            edges_out.push(e_d);
            centers.push(c_d);
            widths.push(w_d);
        }

        Ok(CartesianGrid {
            ndim,
            n,
            extent,
            edges: edges_out,
            centers,
            widths,
            uniform: vec![true; ndim],
            ghosts,
            dtype,
        })
    }
}

fn axis_count_name(d: usize) -> &'static str {
    match d {
        0 => "nx",
        1 => "ny",
        _ => "nz",
    }
}

fn is_uniform(widths: &[f64]) -> bool {
    if widths.len() <= 1 {
        return true;
    }
    let w0 = widths[0];
    let tol = (f64::EPSILON * w0.abs()).max(f64::EPSILON) * 8.0;
    widths.iter().all(|w| (w - w0).abs() <= tol)
}

// --- tests -----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn g_1d(nx: usize) -> CartesianGrid {
        builder().nx(nx).extent(vec![(0.0, 1.0)]).build().unwrap()
    }

    fn g_2d(nx: usize, ny: usize) -> CartesianGrid {
        builder()
            .nx(nx)
            .ny(ny)
            .extent(vec![(0.0, 1.0), (0.0, 0.5)])
            .build()
            .unwrap()
    }

    fn g_3d(nx: usize, ny: usize, nz: usize) -> CartesianGrid {
        builder()
            .nx(nx)
            .ny(ny)
            .nz(nz)
            .extent(vec![(0.0, 1.0), (0.0, 0.5), (0.0, 0.25)])
            .build()
            .unwrap()
    }

    // --- builder validation ---

    #[test]
    fn builder_requires_nx() {
        let err = builder().extent(vec![(0.0, 1.0)]).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("nx")));
    }

    #[test]
    fn builder_requires_extent_in_uniform_mode() {
        let err = builder().nx(4).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("extent")));
    }

    #[test]
    fn builder_rejects_zero_nx() {
        let err = builder()
            .nx(0)
            .extent(vec![(0.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nx", _)));
    }

    #[test]
    fn builder_rejects_zero_ny() {
        let err = builder()
            .nx(4)
            .ny(0)
            .extent(vec![(0.0, 1.0), (0.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("ny", _)));
    }

    #[test]
    fn builder_rejects_nz_without_ny() {
        let err = builder()
            .nx(4)
            .nz(2)
            .extent(vec![(0.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nz", _)));
    }

    #[test]
    fn builder_rejects_extent_arity_mismatch() {
        let err = builder()
            .nx(4)
            .ny(2)
            .extent(vec![(0.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
    }

    #[test]
    fn builder_rejects_hi_not_greater_than_lo() {
        let err = builder()
            .nx(4)
            .extent(vec![(1.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
        let err = builder()
            .nx(4)
            .extent(vec![(2.0, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
    }

    #[test]
    fn builder_rejects_nonfinite_extent() {
        let err = builder()
            .nx(4)
            .extent(vec![(f64::NAN, 1.0)])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
    }

    #[test]
    fn builder_rejects_mixing_edges_and_nx() {
        let err = builder()
            .nx(4)
            .edges(vec![vec![0.0, 0.5, 1.0]])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
    }

    #[test]
    fn builder_rejects_mixing_edges_and_extent() {
        let err = builder()
            .extent(vec![(0.0, 1.0)])
            .edges(vec![vec![0.0, 0.5, 1.0]])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
    }

    #[test]
    fn builder_rejects_short_edges() {
        let err = builder().edges(vec![vec![0.0]]).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
    }

    #[test]
    fn builder_rejects_nonincreasing_edges() {
        let err = builder()
            .edges(vec![vec![0.0, 0.5, 0.5, 1.0]])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
        let err = builder()
            .edges(vec![vec![0.0, 0.5, 0.4, 1.0]])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
    }

    #[test]
    fn builder_rejects_ndim_gt_3() {
        let err = builder()
            .edges(vec![
                vec![0.0, 1.0],
                vec![0.0, 1.0],
                vec![0.0, 1.0],
                vec![0.0, 1.0],
            ])
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("edges", _)));
    }

    #[test]
    fn defaults_match_contract() {
        let g = g_1d(4);
        assert_eq!(g.ndim(), 1);
        assert_eq!(g.n(), &[4]);
        assert_eq!(g.ghosts(), 0);
        assert_eq!(g.dtype(), Dtype::F64);
        assert_eq!(g.family(), "cartesian");
        assert_eq!(g.topology(), "rectilinear");
        assert_eq!(g.uniform(), &[true]);
        assert_eq!(g.n_cells(), 4);
    }

    // --- uniform geometry ---

    #[test]
    fn uniform_edges_span_extent() {
        let g = g_1d(4);
        let e = &g.edges()[0];
        assert_eq!(e.len(), 5);
        assert_relative_eq!(e[0], 0.0, epsilon = 1e-15);
        assert_relative_eq!(e[4], 1.0, epsilon = 1e-15);
    }

    #[test]
    fn uniform_cell_widths_constant() {
        let g = g_1d(4);
        for i in 0..4 {
            let w = g.cell_widths(&[i]).unwrap();
            assert_relative_eq!(w[0], 0.25, epsilon = 1e-15);
        }
    }

    #[test]
    fn uniform_cell_centers_midpoints() {
        let g = g_1d(4);
        for i in 0..4 {
            let c = g.cell_center(&[i]).unwrap();
            assert_relative_eq!(c[0], 0.125 + i as f64 * 0.25, epsilon = 1e-15);
        }
    }

    #[test]
    fn cell_volume_is_width_product() {
        let g = g_2d(4, 2);
        // dx = 0.25, dy = 0.25 -> area = 0.0625
        let v = g.cell_volume(&[2, 1]).unwrap();
        assert_relative_eq!(v, 0.0625, epsilon = 1e-15);
    }

    #[test]
    fn cell_volume_3d_is_dx_dy_dz() {
        let g = g_3d(4, 2, 5);
        // dx=0.25, dy=0.25, dz=0.05 -> vol=0.003125
        let v = g.cell_volume(&[0, 0, 0]).unwrap();
        assert_relative_eq!(v, 0.25 * 0.25 * 0.05, epsilon = 1e-15);
    }

    // --- non-uniform geometry ---

    #[test]
    fn nonuniform_reports_uniform_false_only_when_nonuniform() {
        let g = builder()
            .edges(vec![vec![0.0, 0.1, 0.5, 1.0]])
            .build()
            .unwrap();
        assert_eq!(g.uniform(), &[false]);
        let g = builder()
            .edges(vec![vec![0.0, 0.25, 0.5, 0.75, 1.0]])
            .build()
            .unwrap();
        assert_eq!(g.uniform(), &[true]);
    }

    #[test]
    fn nonuniform_widths_match_edge_diffs() {
        let g = builder()
            .edges(vec![vec![0.0, 0.1, 0.5, 1.0]])
            .build()
            .unwrap();
        let w = g.widths();
        assert_relative_eq!(w[0][0], 0.1, epsilon = 1e-15);
        assert_relative_eq!(w[0][1], 0.4, epsilon = 1e-15);
        assert_relative_eq!(w[0][2], 0.5, epsilon = 1e-15);
    }

    #[test]
    fn nonuniform_extent_is_first_and_last_edge() {
        let g = builder()
            .edges(vec![vec![-2.0, -1.0, 0.5, 3.0]])
            .build()
            .unwrap();
        assert_eq!(g.extent(), &[(-2.0, 3.0)]);
    }

    // --- neighbors ---

    #[test]
    fn neighbors_interior_1d() {
        let g = g_1d(4);
        let ns = g.neighbors(&[1]).unwrap();
        assert_eq!(
            ns,
            vec![
                CartesianNeighbor {
                    axis: 0,
                    side: -1,
                    index: vec![0]
                },
                CartesianNeighbor {
                    axis: 0,
                    side: 1,
                    index: vec![2]
                },
            ]
        );
    }

    #[test]
    fn neighbors_boundary_omits_entries() {
        let g = g_1d(4);
        let ns = g.neighbors(&[0]).unwrap();
        assert_eq!(ns.len(), 1);
        assert_eq!(ns[0].side, 1);
        let ns = g.neighbors(&[3]).unwrap();
        assert_eq!(ns.len(), 1);
        assert_eq!(ns[0].side, -1);
    }

    #[test]
    fn neighbors_3d_corner_has_three() {
        let g = g_3d(4, 2, 5);
        let ns = g.neighbors(&[0, 0, 0]).unwrap();
        assert_eq!(ns.len(), 3);
        for n in ns {
            assert_eq!(n.side, 1);
        }
    }

    #[test]
    fn neighbors_3d_interior_has_six() {
        let g = g_3d(4, 3, 5);
        let ns = g.neighbors(&[1, 1, 2]).unwrap();
        assert_eq!(ns.len(), 6);
    }

    #[test]
    fn neighbors_rejects_out_of_range() {
        let g = g_2d(4, 2);
        assert!(g.neighbors(&[4, 0]).is_err());
        assert!(g.neighbors(&[0, 2]).is_err());
    }

    #[test]
    fn neighbors_rejects_wrong_arity() {
        let g = g_2d(4, 2);
        assert!(g.neighbors(&[0]).is_err());
        assert!(g.neighbors(&[0, 0, 0]).is_err());
    }

    // --- metrics ---

    #[test]
    fn metric_volume_equals_cell_volume() {
        let g = g_2d(4, 2);
        for i in 0..4 {
            for j in 0..2 {
                let m = g.metric_eval(MetricName::Volume, &[i, j]).unwrap();
                let v = g.cell_volume(&[i, j]).unwrap();
                assert_relative_eq!(m, v, epsilon = 1e-15);
            }
        }
    }

    #[test]
    fn metric_jacobian_is_one() {
        let g = g_3d(4, 2, 5);
        let j = g.metric_eval(MetricName::Jacobian, &[0, 0, 0]).unwrap();
        assert_relative_eq!(j, 1.0, epsilon = 1e-15);
    }

    #[test]
    fn metric_dx_dy_dz_match_widths() {
        let g = g_3d(4, 2, 5);
        let dx = g.metric_eval(MetricName::Dx, &[2, 1, 3]).unwrap();
        let dy = g.metric_eval(MetricName::Dy, &[2, 1, 3]).unwrap();
        let dz = g.metric_eval(MetricName::Dz, &[2, 1, 3]).unwrap();
        assert_relative_eq!(dx, 0.25, epsilon = 1e-15);
        assert_relative_eq!(dy, 0.25, epsilon = 1e-15);
        assert_relative_eq!(dz, 0.05, epsilon = 1e-15);
    }

    #[test]
    fn metric_dy_rejected_in_1d() {
        let g = g_1d(4);
        assert!(g.metric_eval(MetricName::Dy, &[0]).is_err());
    }

    #[test]
    fn metric_dz_rejected_in_2d() {
        let g = g_2d(4, 2);
        assert!(g.metric_eval(MetricName::Dz, &[0, 0]).is_err());
    }

    #[test]
    fn face_area_1d_is_unit() {
        let g = g_1d(4);
        let a = g.metric_eval(MetricName::FaceAreaX, &[0]).unwrap();
        assert_relative_eq!(a, 1.0, epsilon = 1e-15);
    }

    #[test]
    fn face_area_2d_orthogonal_width() {
        let g = g_2d(4, 2);
        // dx = 0.25, dy = 0.25
        let ax = g.metric_eval(MetricName::FaceAreaX, &[0, 0]).unwrap();
        let ay = g.metric_eval(MetricName::FaceAreaY, &[0, 0]).unwrap();
        assert_relative_eq!(ax, 0.25, epsilon = 1e-15); // dy
        assert_relative_eq!(ay, 0.25, epsilon = 1e-15); // dx
    }

    #[test]
    fn face_area_3d_is_orthogonal_width_product() {
        let g = g_3d(4, 2, 5);
        // dx=0.25, dy=0.25, dz=0.05
        let ax = g.metric_eval(MetricName::FaceAreaX, &[0, 0, 0]).unwrap();
        let ay = g.metric_eval(MetricName::FaceAreaY, &[0, 0, 0]).unwrap();
        let az = g.metric_eval(MetricName::FaceAreaZ, &[0, 0, 0]).unwrap();
        assert_relative_eq!(ax, 0.25 * 0.05, epsilon = 1e-15);
        assert_relative_eq!(ay, 0.25 * 0.05, epsilon = 1e-15);
        assert_relative_eq!(az, 0.25 * 0.25, epsilon = 1e-15);
    }

    #[test]
    fn metric_eval_by_name_parses_wire_form() {
        let g = g_2d(4, 2);
        let v = g.metric_eval_by_name("volume", &[0, 0]).unwrap();
        assert_relative_eq!(v, 0.0625, epsilon = 1e-15);
    }

    #[test]
    fn metric_eval_by_name_rejects_unknown() {
        let g = g_2d(4, 2);
        let err = g.metric_eval_by_name("not_a_metric", &[0, 0]).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn metric_eval_rejects_out_of_range() {
        let g = g_2d(4, 2);
        assert!(g.metric_eval(MetricName::Volume, &[4, 0]).is_err());
        assert!(g.metric_eval(MetricName::Volume, &[0, 2]).is_err());
    }

    #[test]
    fn metric_g_is_identity() {
        let g = g_3d(2, 2, 2);
        let m = g.metric_g(&[0, 0, 0]).unwrap();
        for (i, row) in m.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(v, expected, epsilon = 1e-15);
            }
        }
    }

    // --- to_esm ---

    #[test]
    fn to_esm_uniform_is_declarative() {
        let g = builder()
            .nx(8)
            .ny(4)
            .nz(2)
            .extent(vec![(0.0, 1.0), (0.0, 0.5), (0.0, 0.25)])
            .ghosts(1)
            .build()
            .unwrap();
        let doc = g.to_esm();
        assert_eq!(doc["family"], "cartesian");
        assert_eq!(doc["version"], "1.0.0");
        assert_eq!(doc["dtype"], "float64");
        assert_eq!(doc["topology"], "rectilinear");
        assert_eq!(doc["ndim"], 3);
        assert_eq!(doc["ghosts"], 1);
        assert_eq!(doc["n_cells"], 8 * 4 * 2);
        assert_eq!(doc["n"], json!([8, 4, 2]));
        assert_eq!(doc["uniform"], json!([true, true, true]));
        // Uniform grids omit edges to keep the wire form small.
        assert!(doc.get("edges").is_none());
        assert!(doc.get("provenance").is_some());

        let wire = serde_json::to_string(&doc).unwrap();
        // No serialized geometry payload: just the n_cells scalar, no
        // enumerated per-cell data.
        assert!(!wire.contains("\"cell_centers\""));
        assert!(!wire.contains("\"widths\""));
        assert!(wire.len() < 2_000, "wire too large: {} bytes", wire.len());
    }

    #[test]
    fn to_esm_nonuniform_carries_edges() {
        let g = builder()
            .edges(vec![vec![0.0, 0.1, 0.5, 1.0]])
            .build()
            .unwrap();
        let doc = g.to_esm();
        assert_eq!(doc["ndim"], 1);
        assert_eq!(doc["uniform"], json!([false]));
        assert_eq!(doc["edges"], json!([[0.0, 0.1, 0.5, 1.0]]));
        assert_eq!(doc["extent"], json!([[0.0, 1.0]]));
    }

    #[test]
    fn to_esm_mixed_uniformity_edges_are_sparse() {
        // x non-uniform, y uniform -> edges array carries only the
        // non-uniform axis's edges; the uniform axis gets an empty array.
        let g = builder()
            .edges(vec![vec![0.0, 0.1, 1.0], vec![0.0, 0.25, 0.5, 0.75, 1.0]])
            .build()
            .unwrap();
        let doc = g.to_esm();
        assert_eq!(doc["uniform"], json!([false, true]));
        assert_eq!(doc["edges"][0], json!([0.0, 0.1, 1.0]));
        assert_eq!(doc["edges"][1], json!([]));
    }

    #[test]
    fn to_esm_roundtrips_via_json() {
        let g = builder()
            .nx(4)
            .ny(2)
            .extent(vec![(0.0, 1.0), (0.0, 0.5)])
            .build()
            .unwrap();
        let doc = g.to_esm();
        let text = serde_json::to_string(&doc).unwrap();
        let reparsed: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(reparsed["ndim"], 2);
        assert_eq!(reparsed["n"], json!([4, 2]));
    }

    #[test]
    fn provenance_identifies_rust_binding() {
        let g = g_1d(4);
        let doc = g.to_esm();
        let prov = &doc["provenance"];
        assert_eq!(prov["binding"], "rust");
        assert_eq!(prov["binding_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(prov["generator"], "cartesian");
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let g = builder()
            .nx(4)
            .extent(vec![(0.0, 1.0)])
            .dtype(Dtype::F32)
            .build()
            .unwrap();
        assert_eq!(g.dtype(), Dtype::F32);
        assert_eq!(g.to_esm()["dtype"], "float32");
    }

    // --- is_uniform helper ---

    #[test]
    fn is_uniform_singleton() {
        assert!(is_uniform(&[1.0]));
        assert!(is_uniform(&[]));
    }

    #[test]
    fn is_uniform_detects_equal_widths() {
        assert!(is_uniform(&[0.25, 0.25, 0.25, 0.25]));
        assert!(!is_uniform(&[0.1, 0.4, 0.5]));
    }
}
