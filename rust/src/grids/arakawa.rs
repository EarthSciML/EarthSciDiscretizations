//! Arakawa staggering accessor runtime.
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.5, §3.3, §7.
//!
//! Arakawa staggering is a *transform* over an underlying base grid. Given a
//! base grid (Cartesian or lat-lon) plus a stagger label (A/B/C/D/E), the
//! runtime provides on-demand accessors for cell centers, u-face, v-face, and
//! corner locations per the staggering convention.
//!
//! Per the 2026-04-20 scope correction, the `.esm` lowering is a small
//! declarative config (family + base-grid ref + stagger + dimensions +
//! extents); geometry is derived from that config by the accessors, not
//! serialized as an inline blob.
//!
//! Stagger conventions (2-D, horizontal):
//!
//! * A: h, u, v colocated at cell centers.
//! * B: h at cell centers; u, v colocated at corners.
//! * C: h at cell centers; u at u-faces (east-west); v at v-faces (north-south).
//! * D: h at cell centers; u at v-faces; v at u-faces (swapped from C).
//! * E: rotated B-grid; topologically equivalent to B with a 45° rotation flag
//!   in the lowered config.

use serde_json::{json, Value};

use crate::{Dtype, Grid, GridError, Result};

/// Arakawa staggering label.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Stagger {
    A,
    B,
    C,
    D,
    E,
}

impl Stagger {
    /// Canonical wire-form name (`"A"` .. `"E"`).
    pub fn as_str(self) -> &'static str {
        match self {
            Stagger::A => "A",
            Stagger::B => "B",
            Stagger::C => "C",
            Stagger::D => "D",
            Stagger::E => "E",
        }
    }

    /// Parse the canonical wire-form name.
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "A" => Stagger::A,
            "B" => Stagger::B,
            "C" => Stagger::C,
            "D" => Stagger::D,
            "E" => Stagger::E,
            _ => return None,
        })
    }

    /// Per-stagger `(h_loc, u_loc, v_loc)` variable location triple.
    pub fn variable_locations(self) -> (Location, Location, Location) {
        match self {
            Stagger::A => (
                Location::CellCenter,
                Location::CellCenter,
                Location::CellCenter,
            ),
            Stagger::B => (Location::CellCenter, Location::Corner, Location::Corner),
            Stagger::C => (Location::CellCenter, Location::UEdge, Location::VEdge),
            Stagger::D => (Location::CellCenter, Location::VEdge, Location::UEdge),
            // E: topologically like B, 45° rotation carried separately.
            Stagger::E => (Location::CellCenter, Location::Corner, Location::Corner),
        }
    }
}

/// Staggered-variable location on a base grid of `(nx, ny)` interior cells.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Location {
    /// Cell centers — shape `(nx, ny)`.
    CellCenter,
    /// East-west cell faces — shape `(nx + 1, ny)`.
    UEdge,
    /// North-south cell faces — shape `(nx, ny + 1)`.
    VEdge,
    /// Cell corners — shape `(nx + 1, ny + 1)`.
    Corner,
}

impl Location {
    /// Array shape of `loc` on a base grid of `(nx, ny)` interior cells.
    pub fn shape(self, nx: usize, ny: usize) -> (usize, usize) {
        match self {
            Location::CellCenter => (nx, ny),
            Location::UEdge => (nx + 1, ny),
            Location::VEdge => (nx, ny + 1),
            Location::Corner => (nx + 1, ny + 1),
        }
    }
}

/// Named variable addressable on an [`ArakawaGrid`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Variable {
    H,
    U,
    V,
}

// ---------------------------------------------------------------------------
// Base grid abstraction.
//
// Phase 1 (cartesian) and Phase 2/3 (lat-lon) grid types are not yet on main.
// Arakawa is parameterised over a `BaseGrid` trait; a minimal `CartesianBase`
// is provided here so the accessors are exercisable today. Future base-grid
// families implement the same trait and plug in without touching this file.
// ---------------------------------------------------------------------------

/// Primitives every arakawa base grid must expose.
pub trait BaseGrid: std::fmt::Debug {
    /// Interior cells along x.
    fn nx(&self) -> usize;
    /// Interior cells along y.
    fn ny(&self) -> usize;
    /// `(x, y)` of cell center `(i, j)` (0-based).
    fn cell_center(&self, i: usize, j: usize) -> (f64, f64);
    /// `(x, y)` of the u-face location at `(i, j)`.
    fn x_edge(&self, i: usize, j: usize) -> (f64, f64);
    /// `(x, y)` of the v-face location at `(i, j)`.
    fn y_edge(&self, i: usize, j: usize) -> (f64, f64);
    /// `(x, y)` of the cell corner at `(i, j)`.
    fn corner(&self, i: usize, j: usize) -> (f64, f64);
    /// Uniform cell spacing in x (for `metric_eval(:dx, ...)`).
    fn dx(&self) -> f64;
    /// Uniform cell spacing in y (for `metric_eval(:dy, ...)`).
    fn dy(&self) -> f64;
    /// Declarative `.esm`-shaped summary of the base grid.
    fn to_esm(&self) -> Value;
}

/// Minimal Cartesian base grid: uniform rectangular mesh over `[xlo, xhi] × [ylo, yhi]`.
#[derive(Clone, Debug, PartialEq)]
pub struct CartesianBase {
    xlo: f64,
    xhi: f64,
    ylo: f64,
    yhi: f64,
    nx: usize,
    ny: usize,
}

impl CartesianBase {
    /// Build a validated `CartesianBase`.
    pub fn new(xlo: f64, xhi: f64, ylo: f64, yhi: f64, nx: usize, ny: usize) -> Result<Self> {
        if nx == 0 {
            return Err(GridError::InvalidOption("nx", "must be >= 1".into()));
        }
        if ny == 0 {
            return Err(GridError::InvalidOption("ny", "must be >= 1".into()));
        }
        if !(xlo.is_finite() && xhi.is_finite() && ylo.is_finite() && yhi.is_finite()) {
            return Err(GridError::InvalidOption(
                "extent",
                "xlo/xhi/ylo/yhi must be finite".into(),
            ));
        }
        if xhi <= xlo {
            return Err(GridError::InvalidOption(
                "extent",
                format!("xhi must be > xlo, got xlo={xlo}, xhi={xhi}"),
            ));
        }
        if yhi <= ylo {
            return Err(GridError::InvalidOption(
                "extent",
                format!("yhi must be > ylo, got ylo={ylo}, yhi={yhi}"),
            ));
        }
        Ok(Self {
            xlo,
            xhi,
            ylo,
            yhi,
            nx,
            ny,
        })
    }

    pub fn xlo(&self) -> f64 {
        self.xlo
    }
    pub fn xhi(&self) -> f64 {
        self.xhi
    }
    pub fn ylo(&self) -> f64 {
        self.ylo
    }
    pub fn yhi(&self) -> f64 {
        self.yhi
    }
}

impl BaseGrid for CartesianBase {
    fn nx(&self) -> usize {
        self.nx
    }
    fn ny(&self) -> usize {
        self.ny
    }

    fn cell_center(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        let dy = self.dy();
        (
            self.xlo + (i as f64 + 0.5) * dx,
            self.ylo + (j as f64 + 0.5) * dy,
        )
    }

    fn x_edge(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        let dy = self.dy();
        (self.xlo + (i as f64) * dx, self.ylo + (j as f64 + 0.5) * dy)
    }

    fn y_edge(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        let dy = self.dy();
        (self.xlo + (i as f64 + 0.5) * dx, self.ylo + (j as f64) * dy)
    }

    fn corner(&self, i: usize, j: usize) -> (f64, f64) {
        let dx = self.dx();
        let dy = self.dy();
        (self.xlo + (i as f64) * dx, self.ylo + (j as f64) * dy)
    }

    fn dx(&self) -> f64 {
        (self.xhi - self.xlo) / (self.nx as f64)
    }

    fn dy(&self) -> f64 {
        (self.yhi - self.ylo) / (self.ny as f64)
    }

    fn to_esm(&self) -> Value {
        json!({
            "family": "cartesian",
            "nx": self.nx,
            "ny": self.ny,
            "extent": [[self.xlo, self.ylo], [self.xhi, self.yhi]],
        })
    }
}

// ---------------------------------------------------------------------------
// Metric field names.
// ---------------------------------------------------------------------------

/// Metric fields accepted by [`ArakawaGrid::metric_eval`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MetricName {
    /// Zonal spacing at the cell.
    Dx,
    /// Meridional spacing at the cell.
    Dy,
    /// Cell area / volume.
    Area,
}

impl MetricName {
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "dx" => MetricName::Dx,
            "dy" => MetricName::Dy,
            "area" => MetricName::Area,
            _ => return None,
        })
    }
}

// ---------------------------------------------------------------------------
// Arakawa grid: base grid + stagger.
// ---------------------------------------------------------------------------

/// Staggered grid: a [`BaseGrid`] (boxed) plus a stagger label.
///
/// The base is stored as a `Box<dyn BaseGrid>` so new base families can be
/// plugged in without altering this struct's type signature.
#[derive(Debug)]
pub struct ArakawaGrid {
    base: Box<dyn BaseGrid + Send + Sync>,
    stagger: Stagger,
    dtype: Dtype,
    ghosts: usize,
}

impl ArakawaGrid {
    /// Interior cells along x.
    pub fn nx(&self) -> usize {
        self.base.nx()
    }

    /// Interior cells along y.
    pub fn ny(&self) -> usize {
        self.base.ny()
    }

    /// Staggering label.
    pub fn stagger(&self) -> Stagger {
        self.stagger
    }

    /// Halo cell width.
    pub fn ghosts(&self) -> usize {
        self.ghosts
    }

    /// Topology string.
    pub fn topology(&self) -> &'static str {
        "block_structured"
    }

    /// Total interior cell count.
    pub fn n_cells(&self) -> usize {
        self.base.nx() * self.base.ny()
    }

    /// Shape of `loc` on this grid: `(ni, nj)`.
    pub fn location_shape(&self, loc: Location) -> (usize, usize) {
        loc.shape(self.base.nx(), self.base.ny())
    }

    /// Shape of the location that `var` lives on under this grid's stagger.
    pub fn variable_shape(&self, var: Variable) -> (usize, usize) {
        self.location_shape(self.variable_location(var))
    }

    /// Location of `var` under this grid's stagger.
    pub fn variable_location(&self, var: Variable) -> Location {
        let (h, u, v) = self.stagger.variable_locations();
        match var {
            Variable::H => h,
            Variable::U => u,
            Variable::V => v,
        }
    }

    /// `(x, y)` of cell center `(i, j)` (0-based).
    pub fn cell_centers(&self, i: usize, j: usize) -> Result<(f64, f64)> {
        self.check_bounds(Location::CellCenter, i, j)?;
        Ok(self.base.cell_center(i, j))
    }

    /// `(x, y)` of the u-variable location at `(i, j)` under this stagger.
    pub fn u_face(&self, i: usize, j: usize) -> Result<(f64, f64)> {
        let loc = self.variable_location(Variable::U);
        self.check_bounds(loc, i, j)?;
        Ok(self.coord_at(loc, i, j))
    }

    /// `(x, y)` of the v-variable location at `(i, j)` under this stagger.
    pub fn v_face(&self, i: usize, j: usize) -> Result<(f64, f64)> {
        let loc = self.variable_location(Variable::V);
        self.check_bounds(loc, i, j)?;
        Ok(self.coord_at(loc, i, j))
    }

    /// `(x, y)` of cell corner `(i, j)` (0-based, `0..nx+1 × 0..ny+1`).
    pub fn corners(&self, i: usize, j: usize) -> Result<(f64, f64)> {
        self.check_bounds(Location::Corner, i, j)?;
        Ok(self.base.corner(i, j))
    }

    /// `(x, y)` of `(loc, i, j)` — generic accessor.
    pub fn coord(&self, loc: Location, i: usize, j: usize) -> Result<(f64, f64)> {
        self.check_bounds(loc, i, j)?;
        Ok(self.coord_at(loc, i, j))
    }

    /// Face neighbors of `(loc, i, j)` as `[W, E, S, N]`.
    ///
    /// Neighbor is `Some((i', j'))` when in-range, `None` at a domain
    /// boundary. All neighbors stay on the same location (no cross-location
    /// stepping — that's the caller's job per stagger semantics).
    pub fn neighbors(
        &self,
        loc: Location,
        i: usize,
        j: usize,
    ) -> Result<[Option<(usize, usize)>; 4]> {
        self.check_bounds(loc, i, j)?;
        let (ni, nj) = self.location_shape(loc);
        let w = if i > 0 { Some((i - 1, j)) } else { None };
        let e = if i + 1 < ni { Some((i + 1, j)) } else { None };
        let s = if j > 0 { Some((i, j - 1)) } else { None };
        let n = if j + 1 < nj { Some((i, j + 1)) } else { None };
        Ok([w, e, s, n])
    }

    /// Cell-center neighbors shortcut: `neighbors(CellCenter, i, j)`.
    pub fn cell_neighbors(&self, i: usize, j: usize) -> Result<[Option<(usize, usize)>; 4]> {
        self.neighbors(Location::CellCenter, i, j)
    }

    /// Evaluate a named metric at cell `(i, j)`.
    pub fn metric_eval(&self, name: MetricName, i: usize, j: usize) -> Result<f64> {
        self.check_bounds(Location::CellCenter, i, j)?;
        let dx = self.base.dx();
        let dy = self.base.dy();
        Ok(match name {
            MetricName::Dx => dx,
            MetricName::Dy => dy,
            MetricName::Area => dx * dy,
        })
    }

    /// Evaluate a metric by its canonical wire-form name (`"dx"`, `"dy"`, `"area"`).
    pub fn metric_eval_by_name(&self, name: &str, i: usize, j: usize) -> Result<f64> {
        let m = MetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, i, j)
    }

    /// Provenance block per `docs/GRIDS_API.md` §6.4.
    pub fn provenance(&self) -> Value {
        json!({
            "binding": "rust",
            "binding_version": env!("CARGO_PKG_VERSION"),
            "source": "earthsci_grids::grids::arakawa",
            "stagger": self.stagger.as_str(),
        })
    }

    // internals -----------------------------------------------------------

    fn coord_at(&self, loc: Location, i: usize, j: usize) -> (f64, f64) {
        match loc {
            Location::CellCenter => self.base.cell_center(i, j),
            Location::UEdge => self.base.x_edge(i, j),
            Location::VEdge => self.base.y_edge(i, j),
            Location::Corner => self.base.corner(i, j),
        }
    }

    fn check_bounds(&self, loc: Location, i: usize, j: usize) -> Result<()> {
        let (ni, nj) = self.location_shape(loc);
        if i >= ni {
            return Err(GridError::InvalidOption(
                "i",
                format!("i out of range [0, {ni}) for {loc:?}: {i}"),
            ));
        }
        if j >= nj {
            return Err(GridError::InvalidOption(
                "j",
                format!("j out of range [0, {nj}) for {loc:?}: {j}"),
            ));
        }
        Ok(())
    }
}

impl Grid for ArakawaGrid {
    fn family(&self) -> &'static str {
        "arakawa"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative `.esm` lowering per the 2026-04-20 scope correction.
    fn to_esm(&self) -> Value {
        json!({
            "family": self.family(),
            "version": "1.0.0",
            "dtype": dtype_str(self.dtype),
            "topology": self.topology(),
            "ghosts": self.ghosts,
            "n_cells": self.n_cells(),
            "stagger": self.stagger.as_str(),
            "rotated": self.stagger == Stagger::E,
            "base": self.base.to_esm(),
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

// ---------------------------------------------------------------------------
// Builder.
// ---------------------------------------------------------------------------

/// Create a new [`ArakawaBuilder`] with default options.
///
/// ```
/// use earthsci_grids::arakawa::{self, CartesianBase, Stagger};
/// use earthsci_grids::Grid;
/// let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
/// let g = arakawa::builder().base(base).stagger(Stagger::C).build().unwrap();
/// assert_eq!(g.family(), "arakawa");
/// ```
pub fn builder() -> ArakawaBuilder {
    ArakawaBuilder::default()
}

/// Builder for [`ArakawaGrid`]. Required options (`base`, `stagger`) are
/// enforced at `.build()` time.
#[must_use]
#[derive(Default)]
pub struct ArakawaBuilder {
    base: Option<Box<dyn BaseGrid + Send + Sync>>,
    stagger: Option<Stagger>,
    dtype: Option<Dtype>,
    ghosts: Option<usize>,
}

impl ArakawaBuilder {
    /// Supply the underlying base grid. Required.
    pub fn base<B: BaseGrid + Send + Sync + 'static>(mut self, base: B) -> Self {
        self.base = Some(Box::new(base));
        self
    }

    /// Supply a boxed base grid (for erased callers). Required.
    pub fn base_boxed(mut self, base: Box<dyn BaseGrid + Send + Sync>) -> Self {
        self.base = Some(base);
        self
    }

    /// Staggering label. Required.
    pub fn stagger(mut self, stagger: Stagger) -> Self {
        self.stagger = Some(stagger);
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
    pub fn build(self) -> Result<ArakawaGrid> {
        let base = self.base.ok_or(GridError::MissingOption("base"))?;
        let stagger = self.stagger.ok_or(GridError::MissingOption("stagger"))?;
        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);
        Ok(ArakawaGrid {
            base,
            stagger,
            dtype,
            ghosts,
        })
    }
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn cartesian(nx: usize, ny: usize) -> CartesianBase {
        CartesianBase::new(0.0, 1.0, 0.0, 1.0, nx, ny).unwrap()
    }

    fn g(stagger: Stagger, nx: usize, ny: usize) -> ArakawaGrid {
        builder()
            .base(cartesian(nx, ny))
            .stagger(stagger)
            .build()
            .unwrap()
    }

    // --- stagger → variable location table ---

    #[test]
    fn stagger_variable_locations_match_spec() {
        assert_eq!(
            Stagger::A.variable_locations(),
            (
                Location::CellCenter,
                Location::CellCenter,
                Location::CellCenter
            )
        );
        assert_eq!(
            Stagger::B.variable_locations(),
            (Location::CellCenter, Location::Corner, Location::Corner)
        );
        assert_eq!(
            Stagger::C.variable_locations(),
            (Location::CellCenter, Location::UEdge, Location::VEdge)
        );
        assert_eq!(
            Stagger::D.variable_locations(),
            (Location::CellCenter, Location::VEdge, Location::UEdge)
        );
        // E: topologically like B.
        assert_eq!(
            Stagger::E.variable_locations(),
            (Location::CellCenter, Location::Corner, Location::Corner)
        );
    }

    #[test]
    fn location_shape_table() {
        assert_eq!(Location::CellCenter.shape(10, 20), (10, 20));
        assert_eq!(Location::UEdge.shape(10, 20), (11, 20));
        assert_eq!(Location::VEdge.shape(10, 20), (10, 21));
        assert_eq!(Location::Corner.shape(10, 20), (11, 21));
    }

    // --- builder / error contract ---

    #[test]
    fn cartesian_base_rejects_zero_nx() {
        let err = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 0, 4).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("nx", _)));
    }

    #[test]
    fn cartesian_base_rejects_reversed_extent() {
        let err = CartesianBase::new(1.0, 0.0, 0.0, 1.0, 4, 4).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
    }

    #[test]
    fn cartesian_base_rejects_nonfinite() {
        let err = CartesianBase::new(f64::NAN, 1.0, 0.0, 1.0, 4, 4).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("extent", _)));
    }

    #[test]
    fn builder_requires_base() {
        let err = builder().stagger(Stagger::C).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("base")));
    }

    #[test]
    fn builder_requires_stagger() {
        let err = builder().base(cartesian(4, 4)).build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("stagger")));
    }

    #[test]
    fn defaults_match_contract() {
        let gr = g(Stagger::C, 4, 4);
        assert_eq!(gr.family(), "arakawa");
        assert_eq!(gr.topology(), "block_structured");
        assert_eq!(gr.dtype(), Dtype::F64);
        assert_eq!(gr.ghosts(), 0);
        assert_eq!(gr.nx(), 4);
        assert_eq!(gr.ny(), 4);
        assert_eq!(gr.n_cells(), 16);
        assert_eq!(gr.stagger(), Stagger::C);
    }

    // --- C-grid accessor recovers Cartesian coordinates ---

    #[test]
    fn c_grid_accessors_recover_cartesian_coordinates() {
        let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 10, 10).unwrap();
        let gr = builder().base(base).stagger(Stagger::C).build().unwrap();

        let (cx, cy) = gr.cell_centers(0, 0).unwrap();
        assert_relative_eq!(cx, 0.05, epsilon = 1e-14);
        assert_relative_eq!(cy, 0.05, epsilon = 1e-14);

        let (cx, cy) = gr.cell_centers(9, 9).unwrap();
        assert_relative_eq!(cx, 0.95, epsilon = 1e-14);
        assert_relative_eq!(cy, 0.95, epsilon = 1e-14);

        // C: u at u-faces. u(0, 0) sits on western boundary at mid-y of row 0.
        let (ux, uy) = gr.u_face(0, 0).unwrap();
        assert_relative_eq!(ux, 0.0, epsilon = 1e-14);
        assert_relative_eq!(uy, 0.05, epsilon = 1e-14);
        // u(10, 0) is the eastern boundary of the last column.
        let (ux, _) = gr.u_face(10, 0).unwrap();
        assert_relative_eq!(ux, 1.0, epsilon = 1e-14);

        // C: v at v-faces.
        let (vx, vy) = gr.v_face(0, 0).unwrap();
        assert_relative_eq!(vx, 0.05, epsilon = 1e-14);
        assert_relative_eq!(vy, 0.0, epsilon = 1e-14);
        let (_, vy) = gr.v_face(0, 10).unwrap();
        assert_relative_eq!(vy, 1.0, epsilon = 1e-14);

        // Corners.
        assert_eq!(gr.corners(0, 0).unwrap(), (0.0, 0.0));
        assert_eq!(gr.corners(10, 10).unwrap(), (1.0, 1.0));
    }

    // --- A-grid: colocation ---

    #[test]
    fn a_grid_colocates_u_v_at_centers() {
        let gr = g(Stagger::A, 4, 4);
        assert_eq!(gr.variable_shape(Variable::H), (4, 4));
        assert_eq!(gr.variable_shape(Variable::U), (4, 4));
        assert_eq!(gr.variable_shape(Variable::V), (4, 4));
        assert_eq!(gr.u_face(1, 2).unwrap(), gr.cell_centers(1, 2).unwrap());
        assert_eq!(gr.v_face(1, 2).unwrap(), gr.cell_centers(1, 2).unwrap());
    }

    // --- B-grid: u,v at corners ---

    #[test]
    fn b_grid_puts_uv_at_corners() {
        let gr = g(Stagger::B, 4, 4);
        assert_eq!(gr.variable_shape(Variable::H), (4, 4));
        assert_eq!(gr.variable_shape(Variable::U), (5, 5));
        assert_eq!(gr.variable_shape(Variable::V), (5, 5));
        assert_eq!(gr.u_face(0, 0).unwrap(), gr.corners(0, 0).unwrap());
        assert_eq!(gr.v_face(2, 1).unwrap(), gr.corners(2, 1).unwrap());
    }

    // --- D-grid: swaps C's u/v faces ---

    #[test]
    fn d_grid_swaps_c_uv_faces() {
        let base_c = CartesianBase::new(0.0, 2.0, 0.0, 1.0, 4, 4).unwrap();
        let base_d = base_c.clone();
        let c = builder().base(base_c).stagger(Stagger::C).build().unwrap();
        let d = builder().base(base_d).stagger(Stagger::D).build().unwrap();
        assert_eq!(d.variable_shape(Variable::U), c.variable_shape(Variable::V));
        assert_eq!(d.variable_shape(Variable::V), c.variable_shape(Variable::U));
        assert_eq!(d.u_face(0, 0).unwrap(), c.v_face(0, 0).unwrap());
        assert_eq!(d.v_face(0, 0).unwrap(), c.u_face(0, 0).unwrap());
    }

    // --- Neighbors + bounds ---

    #[test]
    fn neighbors_interior_cell() {
        let gr = g(Stagger::C, 4, 4);
        let [w, e, s, n] = gr.cell_neighbors(1, 1).unwrap();
        assert_eq!(w, Some((0, 1)));
        assert_eq!(e, Some((2, 1)));
        assert_eq!(s, Some((1, 0)));
        assert_eq!(n, Some((1, 2)));
    }

    #[test]
    fn neighbors_corner_cell_has_none_on_boundary() {
        let gr = g(Stagger::C, 4, 4);
        let [w, e, s, n] = gr.cell_neighbors(0, 0).unwrap();
        assert_eq!(w, None);
        assert_eq!(s, None);
        assert_eq!(e, Some((1, 0)));
        assert_eq!(n, Some((0, 1)));
    }

    #[test]
    fn neighbors_bounds_check() {
        let gr = g(Stagger::C, 4, 4);
        assert!(gr.cell_neighbors(4, 0).is_err());
        assert!(gr.cell_neighbors(0, 4).is_err());
    }

    #[test]
    fn neighbors_uedge_last_column() {
        let gr = g(Stagger::C, 4, 4);
        // UEdge shape is (5, 4); (4, 3) is at the top-right of the UEdge grid.
        let [w, e, s, n] = gr.neighbors(Location::UEdge, 4, 3).unwrap();
        assert_eq!(e, None);
        assert_eq!(n, None);
        assert_eq!(w, Some((3, 3)));
        assert_eq!(s, Some((4, 2)));
    }

    // --- Metric eval ---

    #[test]
    fn metric_eval_cartesian_uniform() {
        let base = CartesianBase::new(0.0, 2.0, 0.0, 6.0, 4, 3).unwrap();
        let gr = builder().base(base).stagger(Stagger::C).build().unwrap();
        assert_relative_eq!(
            gr.metric_eval(MetricName::Dx, 0, 0).unwrap(),
            0.5,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            gr.metric_eval(MetricName::Dy, 0, 0).unwrap(),
            2.0,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            gr.metric_eval(MetricName::Area, 3, 2).unwrap(),
            1.0,
            epsilon = 1e-14
        );
    }

    #[test]
    fn metric_eval_rejects_unknown_name() {
        let gr = g(Stagger::C, 4, 4);
        let err = gr.metric_eval_by_name("not_a_metric", 0, 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    #[test]
    fn metric_eval_rejects_out_of_range_cell() {
        let gr = g(Stagger::C, 4, 4);
        assert!(gr.metric_eval(MetricName::Dx, 4, 0).is_err());
    }

    // --- to_esm: declarative config ---

    #[test]
    fn to_esm_is_declarative() {
        let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 5, 5).unwrap();
        let gr = builder()
            .base(base)
            .stagger(Stagger::C)
            .ghosts(2)
            .build()
            .unwrap();
        let doc = gr.to_esm();
        assert_eq!(doc["family"], "arakawa");
        assert_eq!(doc["version"], "1.0.0");
        assert_eq!(doc["dtype"], "float64");
        assert_eq!(doc["topology"], "block_structured");
        assert_eq!(doc["stagger"], "C");
        assert_eq!(doc["ghosts"], 2);
        assert_eq!(doc["n_cells"], 25);
        assert_eq!(doc["rotated"], false);

        let base_doc = &doc["base"];
        assert_eq!(base_doc["family"], "cartesian");
        assert_eq!(base_doc["nx"], 5);
        assert_eq!(base_doc["ny"], 5);

        // Per §0 correction: no inline geometry arrays.
        let wire = serde_json::to_string(&doc).unwrap();
        assert!(!wire.contains("cells\":["));
        assert!(!wire.contains("edges\":["));
        assert!(!wire.contains("vertices\":["));
        assert!(wire.len() < 2_000, "wire too large: {} bytes", wire.len());
    }

    #[test]
    fn to_esm_rotated_flag_for_e() {
        let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
        let base2 = base.clone();
        let e = builder().base(base).stagger(Stagger::E).build().unwrap();
        let b = builder().base(base2).stagger(Stagger::B).build().unwrap();
        assert_eq!(e.to_esm()["rotated"], true);
        assert_eq!(b.to_esm()["rotated"], false);
        // E mirrors B topologically.
        assert_eq!(e.variable_shape(Variable::H), b.variable_shape(Variable::H));
        assert_eq!(e.variable_shape(Variable::U), b.variable_shape(Variable::U));
        assert_eq!(e.variable_shape(Variable::V), b.variable_shape(Variable::V));
    }

    #[test]
    fn provenance_identifies_rust_binding() {
        let gr = g(Stagger::C, 4, 4);
        let doc = gr.to_esm();
        let prov = &doc["provenance"];
        assert_eq!(prov["binding"], "rust");
        assert_eq!(prov["binding_version"], env!("CARGO_PKG_VERSION"));
        assert_eq!(prov["stagger"], "C");
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
        let gr = builder()
            .base(base)
            .stagger(Stagger::C)
            .dtype(Dtype::F32)
            .build()
            .unwrap();
        assert_eq!(gr.dtype(), Dtype::F32);
        assert_eq!(gr.to_esm()["dtype"], "float32");
    }

    #[test]
    fn stagger_from_name_roundtrips() {
        for s in [Stagger::A, Stagger::B, Stagger::C, Stagger::D, Stagger::E] {
            assert_eq!(Stagger::from_name(s.as_str()), Some(s));
        }
        assert_eq!(Stagger::from_name("Q"), None);
    }
}
