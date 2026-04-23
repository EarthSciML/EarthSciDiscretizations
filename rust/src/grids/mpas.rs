//! MPAS unstructured Voronoi grid family (loader-backed).
//!
//! Conforms to `docs/GRIDS_API.md` §2.5 / §10 and the mayor's 2026-04-20
//! scope correction (bead dsc-3ut / dsc-3nw): the `.esm` lowering is a small
//! declarative config (family, dimensions, loader ref), not a serialized
//! geometry blob. Geometry is derived on demand by the accessors
//! ([`MpasGrid::cell_centers`], [`MpasGrid::neighbors`],
//! [`MpasGrid::cell_area`], [`MpasGrid::edge_length`],
//! [`MpasGrid::metric_eval`]).
//!
//! No new runtime dependencies are introduced. Path-based loading requires
//! the caller to supply a reader via [`Builder::reader_fn`] that translates
//! a filesystem path into an [`MpasMeshData`]; this keeps NetCDF I/O out of
//! the default dependency graph while satisfying the §10 loader contract.
//! In-memory construction via [`MpasMeshData::from_input`] is the primary
//! path for tests and host-built meshes.
//!
//! Index convention: cells, edges, and adjacency slots are 0-based (Rust /
//! TypeScript convention). The sentinel for "no neighbor" in adjacency
//! arrays is `-1` ([`NO_NEIGHBOR`]). This differs from the Julia reference
//! (1-based with `0`-sentinel per NetCDF convention) but does not affect
//! cross-binding conformance, which compares accessor outputs at pinned
//! query points rather than serialized adjacency arrays.

use serde_json::{json, Map, Value};

use crate::{Dtype, Grid, GridError, Result};

const MPAS_FAMILY_VERSION: &str = "1.0.0";
const DEFAULT_EARTH_RADIUS_M: f64 = 6.371e6;

/// Adjacency sentinel for "no neighbor" / "no edge" / "external boundary".
pub const NO_NEIGHBOR: i32 = -1;

/// Loader struct per GRIDS_API §10.
#[derive(Clone, Debug, PartialEq)]
pub struct MpasLoader {
    /// Filesystem path (may be sandboxed). For callers using a pure
    /// in-memory mesh the loader may be omitted entirely.
    pub path: String,
    /// `"auto"` | `"nc4"` | `"mpas_mesh"`. Default: `"auto"`.
    pub reader: String,
    /// `"strict"` | `"lenient"`. Default: `"strict"`.
    pub check: String,
}

impl MpasLoader {
    /// New loader with reader=`auto`, check=`strict` — the contract default.
    pub fn new<P: Into<String>>(path: P) -> Self {
        Self {
            path: path.into(),
            reader: "auto".to_string(),
            check: "strict".to_string(),
        }
    }

    /// Validate the enum-ish string fields.
    fn validate(&self) -> Result<()> {
        if self.path.is_empty() {
            return Err(GridError::InvalidOption(
                "loader.path",
                "must be a non-empty string".to_string(),
            ));
        }
        match self.reader.as_str() {
            "auto" | "nc4" | "mpas_mesh" => {}
            other => {
                return Err(GridError::InvalidOption(
                    "loader.reader",
                    format!("must be 'auto' | 'nc4' | 'mpas_mesh' (got '{other}')"),
                ));
            }
        }
        match self.check.as_str() {
            "strict" | "lenient" => {}
            other => {
                return Err(GridError::InvalidOption(
                    "loader.check",
                    format!("must be 'strict' | 'lenient' (got '{other}')"),
                ));
            }
        }
        Ok(())
    }

    fn to_json(&self) -> Value {
        json!({
            "path": self.path,
            "reader": self.reader,
            "check": self.check,
        })
    }
}

/// In-memory input for building an [`MpasMeshData`]. Array lengths are
/// enforced in [`MpasMeshData::from_input`]; `x_cell`/`y_cell`/`z_cell`
/// default to the spherical projection of (lon_cell, lat_cell, R).
#[derive(Clone, Debug)]
pub struct MpasMeshInput {
    pub lon_cell: Vec<f64>,
    pub lat_cell: Vec<f64>,
    pub area_cell: Vec<f64>,
    pub n_edges_on_cell: Vec<i32>,
    /// Row-major `(n_cells × max_edges)`. Entry `-1` denotes "no neighbor".
    pub cells_on_cell: Vec<i32>,
    /// Row-major `(n_cells × max_edges)`. Entry `-1` denotes "no edge".
    pub edges_on_cell: Vec<i32>,
    pub lon_edge: Vec<f64>,
    pub lat_edge: Vec<f64>,
    /// Row-major `(n_edges × 2)`. Entry `-1` denotes "external boundary".
    pub cells_on_edge: Vec<i32>,
    pub dc_edge: Vec<f64>,
    pub dv_edge: Vec<f64>,
    pub max_edges: usize,
    pub x_cell: Option<Vec<f64>>,
    pub y_cell: Option<Vec<f64>>,
    pub z_cell: Option<Vec<f64>>,
    pub n_vertices: Option<usize>,
    pub r: Option<f64>,
}

/// Fully-materialized MPAS Voronoi mesh. Returned by [`MpasMeshData::from_input`]
/// and by user-supplied `reader_fn` callbacks.
#[derive(Clone, Debug)]
pub struct MpasMeshData {
    pub n_cells: usize,
    pub n_edges: usize,
    pub n_vertices: usize,
    pub max_edges: usize,
    pub lon_cell: Vec<f64>,
    pub lat_cell: Vec<f64>,
    pub x_cell: Vec<f64>,
    pub y_cell: Vec<f64>,
    pub z_cell: Vec<f64>,
    pub area_cell: Vec<f64>,
    pub n_edges_on_cell: Vec<i32>,
    pub cells_on_cell: Vec<i32>,
    pub edges_on_cell: Vec<i32>,
    pub lon_edge: Vec<f64>,
    pub lat_edge: Vec<f64>,
    pub cells_on_edge: Vec<i32>,
    pub dc_edge: Vec<f64>,
    pub dv_edge: Vec<f64>,
}

fn check_len<T>(v: &[T], label: &'static str, expected: usize) -> Result<()> {
    if v.len() != expected {
        return Err(GridError::InvalidOption(
            label,
            format!("length {} != expected {expected}", v.len()),
        ));
    }
    Ok(())
}

impl MpasMeshData {
    /// Validated constructor. Derives cartesian cell coordinates from
    /// `(lon, lat, R)` when not supplied.
    pub fn from_input(input: MpasMeshInput) -> Result<Self> {
        let n_cells = input.lon_cell.len();
        let n_edges = input.lon_edge.len();
        let max_edges = input.max_edges;
        if max_edges == 0 {
            return Err(GridError::InvalidOption(
                "max_edges",
                "must be a positive integer".to_string(),
            ));
        }
        let r = input.r.unwrap_or(DEFAULT_EARTH_RADIUS_M);

        check_len(&input.lat_cell, "lat_cell", n_cells)?;
        check_len(&input.area_cell, "area_cell", n_cells)?;
        check_len(&input.n_edges_on_cell, "n_edges_on_cell", n_cells)?;
        check_len(&input.cells_on_cell, "cells_on_cell", n_cells * max_edges)?;
        check_len(&input.edges_on_cell, "edges_on_cell", n_cells * max_edges)?;
        check_len(&input.lat_edge, "lat_edge", n_edges)?;
        check_len(&input.cells_on_edge, "cells_on_edge", n_edges * 2)?;
        check_len(&input.dc_edge, "dc_edge", n_edges)?;
        check_len(&input.dv_edge, "dv_edge", n_edges)?;

        let x_cell = match input.x_cell {
            Some(v) => {
                check_len(&v, "x_cell", n_cells)?;
                v
            }
            None => input
                .lon_cell
                .iter()
                .zip(&input.lat_cell)
                .map(|(lon, lat)| r * lat.cos() * lon.cos())
                .collect(),
        };
        let y_cell = match input.y_cell {
            Some(v) => {
                check_len(&v, "y_cell", n_cells)?;
                v
            }
            None => input
                .lon_cell
                .iter()
                .zip(&input.lat_cell)
                .map(|(lon, lat)| r * lat.cos() * lon.sin())
                .collect(),
        };
        let z_cell = match input.z_cell {
            Some(v) => {
                check_len(&v, "z_cell", n_cells)?;
                v
            }
            None => input.lat_cell.iter().map(|lat| r * lat.sin()).collect(),
        };

        Ok(MpasMeshData {
            n_cells,
            n_edges,
            n_vertices: input.n_vertices.unwrap_or(0),
            max_edges,
            lon_cell: input.lon_cell,
            lat_cell: input.lat_cell,
            x_cell,
            y_cell,
            z_cell,
            area_cell: input.area_cell,
            n_edges_on_cell: input.n_edges_on_cell,
            cells_on_cell: input.cells_on_cell,
            edges_on_cell: input.edges_on_cell,
            lon_edge: input.lon_edge,
            lat_edge: input.lat_edge,
            cells_on_edge: input.cells_on_edge,
            dc_edge: input.dc_edge,
            dv_edge: input.dv_edge,
        })
    }
}

/// Strict-mode validation: bounds-check adjacency arrays and enforce
/// neighbor-link reciprocity. Lenient mode skips the reciprocity check.
pub fn check_mesh(m: &MpasMeshData, strict: bool) -> Result<()> {
    let n_cells = m.n_cells as i32;
    let n_edges = m.n_edges as i32;
    let max_edges = m.max_edges;
    for c in 0..m.n_cells {
        let k = m.n_edges_on_cell[c];
        if k < 0 || (k as usize) > max_edges {
            return Err(GridError::SchemaViolation(format!(
                "n_edges_on_cell[{c}]={k} out of [0, {max_edges}]"
            )));
        }
        for j in 0..(k as usize) {
            let nb = m.cells_on_cell[c * max_edges + j];
            if nb < NO_NEIGHBOR || nb >= n_cells {
                return Err(GridError::SchemaViolation(format!(
                    "cells_on_cell[{c},{j}]={nb} out of [-1, {}]",
                    n_cells - 1
                )));
            }
            let e = m.edges_on_cell[c * max_edges + j];
            if e < NO_NEIGHBOR || e >= n_edges {
                return Err(GridError::SchemaViolation(format!(
                    "edges_on_cell[{c},{j}]={e} out of [-1, {}]",
                    n_edges - 1
                )));
            }
        }
    }
    for e in 0..m.n_edges {
        for s in 0..2 {
            let c = m.cells_on_edge[e * 2 + s];
            if c < NO_NEIGHBOR || c >= n_cells {
                return Err(GridError::SchemaViolation(format!(
                    "cells_on_edge[{e},{s}]={c} out of [-1, {}]",
                    n_cells - 1
                )));
            }
        }
    }
    if !strict {
        return Ok(());
    }
    for c in 0..m.n_cells {
        let k = m.n_edges_on_cell[c] as usize;
        for j in 0..k {
            let nb = m.cells_on_cell[c * max_edges + j];
            if nb < 0 {
                continue;
            }
            let nb = nb as usize;
            let kb = m.n_edges_on_cell[nb] as usize;
            let mut found = false;
            for jj in 0..kb {
                if m.cells_on_cell[nb * max_edges + jj] == c as i32 {
                    found = true;
                    break;
                }
            }
            if !found {
                return Err(GridError::SchemaViolation(format!(
                    "neighbor symmetry broken: cell {c} -> {nb} but not reverse"
                )));
            }
        }
    }
    Ok(())
}

/// Metric field names accepted by [`MpasGrid::metric_eval`]. Indexes cells
/// for cell-valued metrics and edges for edge-valued metrics.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum MpasMetricName {
    Lon,
    Lat,
    Area,
    X,
    Y,
    Z,
    NEdgesOnCell,
    LonEdge,
    LatEdge,
    DcEdge,
    DvEdge,
}

impl MpasMetricName {
    /// Parse the canonical wire-form name (e.g. `"lon_edge"`).
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "lon" => MpasMetricName::Lon,
            "lat" => MpasMetricName::Lat,
            "area" => MpasMetricName::Area,
            "x" => MpasMetricName::X,
            "y" => MpasMetricName::Y,
            "z" => MpasMetricName::Z,
            "n_edges_on_cell" => MpasMetricName::NEdgesOnCell,
            "lon_edge" => MpasMetricName::LonEdge,
            "lat_edge" => MpasMetricName::LatEdge,
            "dc_edge" => MpasMetricName::DcEdge,
            "dv_edge" => MpasMetricName::DvEdge,
            _ => return None,
        })
    }
}

/// Materialized MPAS grid.
pub struct MpasGrid {
    r: f64,
    dtype: Dtype,
    ghosts: u32,
    loader: Option<MpasLoader>,
    mesh: MpasMeshData,
}

impl std::fmt::Debug for MpasGrid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MpasGrid")
            .field("r", &self.r)
            .field("dtype", &self.dtype)
            .field("ghosts", &self.ghosts)
            .field("loader", &self.loader)
            .field("n_cells", &self.mesh.n_cells)
            .field("n_edges", &self.mesh.n_edges)
            .field("n_vertices", &self.mesh.n_vertices)
            .field("max_edges", &self.mesh.max_edges)
            .finish()
    }
}

impl MpasGrid {
    pub fn radius(&self) -> f64 {
        self.r
    }
    pub fn ghosts(&self) -> u32 {
        self.ghosts
    }
    pub fn n_cells(&self) -> usize {
        self.mesh.n_cells
    }
    pub fn n_edges(&self) -> usize {
        self.mesh.n_edges
    }
    pub fn n_vertices(&self) -> usize {
        self.mesh.n_vertices
    }
    pub fn max_edges(&self) -> usize {
        self.mesh.max_edges
    }
    pub fn loader(&self) -> Option<&MpasLoader> {
        self.loader.as_ref()
    }
    pub fn mesh(&self) -> &MpasMeshData {
        &self.mesh
    }
    pub fn topology(&self) -> &'static str {
        "unstructured"
    }

    fn check_cell(&self, c: usize) -> Result<()> {
        if c >= self.mesh.n_cells {
            return Err(GridError::InvalidOption(
                "cell",
                format!("invalid cell index {c} (expected 0..{})", self.mesh.n_cells),
            ));
        }
        Ok(())
    }

    fn check_edge(&self, e: usize) -> Result<()> {
        if e >= self.mesh.n_edges {
            return Err(GridError::InvalidOption(
                "edge",
                format!("invalid edge index {e} (expected 0..{})", self.mesh.n_edges),
            ));
        }
        Ok(())
    }

    /// Cell-center `(lon, lat)` in radians.
    pub fn cell_centers(&self, c: usize) -> Result<(f64, f64)> {
        self.check_cell(c)?;
        Ok((self.mesh.lon_cell[c], self.mesh.lat_cell[c]))
    }

    /// Cell-center cartesian `(x, y, z)` (same units as `R`).
    pub fn cell_center_cart(&self, c: usize) -> Result<[f64; 3]> {
        self.check_cell(c)?;
        Ok([
            self.mesh.x_cell[c],
            self.mesh.y_cell[c],
            self.mesh.z_cell[c],
        ])
    }

    /// Valid (non-boundary) neighbor cell indices for `c`.
    pub fn neighbors(&self, c: usize) -> Result<Vec<u32>> {
        self.check_cell(c)?;
        let k = self.mesh.n_edges_on_cell[c] as usize;
        let me = self.mesh.max_edges;
        let mut out = Vec::with_capacity(k);
        for j in 0..k {
            let nb = self.mesh.cells_on_cell[c * me + j];
            if nb >= 0 {
                out.push(nb as u32);
            }
        }
        Ok(out)
    }

    /// Raw neighbor slots for cell `c` (length `n_edges_on_cell[c]`),
    /// preserving boundary sentinels as `None`.
    pub fn neighbor_slots(&self, c: usize) -> Result<Vec<Option<u32>>> {
        self.check_cell(c)?;
        let k = self.mesh.n_edges_on_cell[c] as usize;
        let me = self.mesh.max_edges;
        let mut out = Vec::with_capacity(k);
        for j in 0..k {
            let nb = self.mesh.cells_on_cell[c * me + j];
            out.push(if nb < 0 { None } else { Some(nb as u32) });
        }
        Ok(out)
    }

    /// Cell area (meters²).
    pub fn cell_area(&self, c: usize) -> Result<f64> {
        self.check_cell(c)?;
        Ok(self.mesh.area_cell[c])
    }

    /// Edge length (Voronoi arc length `dv_edge`, meters).
    pub fn edge_length(&self, e: usize) -> Result<f64> {
        self.check_edge(e)?;
        Ok(self.mesh.dv_edge[e])
    }

    /// Cell-center-to-cell-center distance across edge `e` (`dc_edge`, meters).
    pub fn cell_distance(&self, e: usize) -> Result<f64> {
        self.check_edge(e)?;
        Ok(self.mesh.dc_edge[e])
    }

    /// Sum of `area_cell` across the mesh.
    pub fn total_area(&self) -> f64 {
        self.mesh.area_cell.iter().sum()
    }

    /// Evaluate a metric field by index. For cell-valued metrics `i` is a
    /// cell index; for edge-valued metrics (`lon_edge`, `lat_edge`,
    /// `dc_edge`, `dv_edge`) it is an edge index.
    pub fn metric_eval(&self, name: MpasMetricName, i: usize) -> Result<f64> {
        use MpasMetricName::*;
        match name {
            Lon => {
                self.check_cell(i)?;
                Ok(self.mesh.lon_cell[i])
            }
            Lat => {
                self.check_cell(i)?;
                Ok(self.mesh.lat_cell[i])
            }
            Area => {
                self.check_cell(i)?;
                Ok(self.mesh.area_cell[i])
            }
            X => {
                self.check_cell(i)?;
                Ok(self.mesh.x_cell[i])
            }
            Y => {
                self.check_cell(i)?;
                Ok(self.mesh.y_cell[i])
            }
            Z => {
                self.check_cell(i)?;
                Ok(self.mesh.z_cell[i])
            }
            NEdgesOnCell => {
                self.check_cell(i)?;
                Ok(self.mesh.n_edges_on_cell[i] as f64)
            }
            LonEdge => {
                self.check_edge(i)?;
                Ok(self.mesh.lon_edge[i])
            }
            LatEdge => {
                self.check_edge(i)?;
                Ok(self.mesh.lat_edge[i])
            }
            DcEdge => {
                self.check_edge(i)?;
                Ok(self.mesh.dc_edge[i])
            }
            DvEdge => {
                self.check_edge(i)?;
                Ok(self.mesh.dv_edge[i])
            }
        }
    }

    /// Evaluate a metric by its canonical wire-form name.
    pub fn metric_eval_by_name(&self, name: &str, i: usize) -> Result<f64> {
        let m = MpasMetricName::from_name(name).ok_or_else(|| {
            GridError::InvalidOption("metric_name", format!("unknown metric: {name:?}"))
        })?;
        self.metric_eval(m, i)
    }

    fn provenance(&self) -> Value {
        let mut obj = Map::new();
        obj.insert("binding".into(), json!("rust"));
        obj.insert("binding_version".into(), json!(env!("CARGO_PKG_VERSION")));
        obj.insert("family".into(), json!("mpas"));
        obj.insert("version".into(), json!(MPAS_FAMILY_VERSION));
        obj.insert("dtype".into(), json!(dtype_str(self.dtype)));
        match &self.loader {
            Some(l) => {
                obj.insert("loader".into(), l.to_json());
            }
            None => {
                obj.insert("loader".into(), Value::Null);
            }
        }
        Value::Object(obj)
    }
}

fn dtype_str(d: Dtype) -> &'static str {
    match d {
        Dtype::F64 => "float64",
        Dtype::F32 => "float32",
    }
}

impl Grid for MpasGrid {
    fn family(&self) -> &'static str {
        "mpas"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// Declarative §6-schema-valid lowering. No inline geometry arrays.
    fn to_esm(&self) -> Value {
        let loader = match &self.loader {
            Some(l) => l.to_json(),
            None => Value::Null,
        };
        let mut options = Map::new();
        options.insert("R".into(), json!(self.r));
        options.insert("loader".into(), loader);

        json!({
            "family": "mpas",
            "version": MPAS_FAMILY_VERSION,
            "dtype": dtype_str(self.dtype),
            "topology": "unstructured",
            "ghosts": self.ghosts,
            "n_cells": self.mesh.n_cells,
            "n_edges": self.mesh.n_edges,
            "n_vertices": self.mesh.n_vertices,
            "max_edges": self.mesh.max_edges,
            "options": Value::Object(options),
            "provenance": self.provenance(),
            "schema_version": MPAS_FAMILY_VERSION,
        })
    }
}

/// Callback that turns a loader path into an [`MpasMeshData`]. See §10.
pub type ReaderFn = Box<dyn FnOnce(&str) -> Result<MpasMeshData>>;

/// Builder for a [`MpasGrid`]. Required options are enforced at
/// [`Builder::build`] time so partial builders can be reused in tests.
#[must_use]
#[derive(Default)]
pub struct Builder {
    loader: Option<MpasLoader>,
    mesh: Option<MpasMeshData>,
    r: Option<f64>,
    dtype: Option<Dtype>,
    ghosts: Option<u32>,
    reader_fn: Option<ReaderFn>,
}

impl std::fmt::Debug for Builder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Builder")
            .field("loader", &self.loader)
            .field("mesh_set", &self.mesh.is_some())
            .field("r", &self.r)
            .field("dtype", &self.dtype)
            .field("ghosts", &self.ghosts)
            .field("reader_fn_set", &self.reader_fn.is_some())
            .finish()
    }
}

impl Builder {
    /// Loader ref per GRIDS_API §10. Required if `mesh` is not provided.
    pub fn loader(mut self, loader: MpasLoader) -> Self {
        self.loader = Some(loader);
        self
    }

    /// Directly supply an in-memory mesh. Required if `loader` is not
    /// provided, or when using a loader for provenance without going
    /// through a `reader_fn`.
    pub fn mesh(mut self, mesh: MpasMeshData) -> Self {
        self.mesh = Some(mesh);
        self
    }

    /// Sphere radius (meters). Defaults to `6.371e6`.
    pub fn r(mut self, r: f64) -> Self {
        self.r = Some(r);
        self
    }

    /// Element precision. Defaults to [`Dtype::F64`].
    pub fn dtype(mut self, dtype: Dtype) -> Self {
        self.dtype = Some(dtype);
        self
    }

    /// Halo cell width. Must be `0` for loader-backed grids.
    pub fn ghosts(mut self, ghosts: u32) -> Self {
        self.ghosts = Some(ghosts);
        self
    }

    /// Caller-supplied reader for path-based loading. NetCDF I/O is not
    /// bundled with `earthsci_grids` per GRIDS_API §10 — consumers plug in
    /// their own reader.
    pub fn reader_fn<F>(mut self, f: F) -> Self
    where
        F: FnOnce(&str) -> Result<MpasMeshData> + 'static,
    {
        self.reader_fn = Some(Box::new(f));
        self
    }

    /// Validate options and construct the grid.
    pub fn build(self) -> Result<MpasGrid> {
        let r = self.r.unwrap_or(DEFAULT_EARTH_RADIUS_M);
        if !(r.is_finite() && r > 0.0) {
            return Err(GridError::InvalidOption(
                "R",
                format!("must be a positive finite number, got {r}"),
            ));
        }
        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);
        if ghosts != 0 {
            return Err(GridError::InvalidOption(
                "ghosts",
                format!("must be 0 for loader-backed grids, got {ghosts}"),
            ));
        }

        let (mesh, loader, strict) = match (self.mesh, self.loader, self.reader_fn) {
            (Some(mesh), loader_opt, _reader) => {
                if let Some(ref l) = loader_opt {
                    l.validate()?;
                }
                // When a mesh is supplied directly, default to strict checks
                // unless the (optional) loader declares lenient mode.
                let strict = loader_opt
                    .as_ref()
                    .map(|l| l.check.as_str() != "lenient")
                    .unwrap_or(true);
                (mesh, loader_opt, strict)
            }
            (None, Some(loader), Some(reader_fn)) => {
                loader.validate()?;
                let mesh = reader_fn(&loader.path)?;
                let strict = loader.check.as_str() != "lenient";
                (mesh, Some(loader), strict)
            }
            (None, Some(_loader), None) => {
                return Err(GridError::InvalidOption(
                    "reader_fn",
                    "path-based loading requires a reader_fn(path) -> MpasMeshData; \
                     NetCDF I/O is not bundled with earthsci_grids per GRIDS_API §10"
                        .to_string(),
                ));
            }
            (None, None, _) => {
                return Err(GridError::MissingOption("loader_or_mesh"));
            }
        };

        check_mesh(&mesh, strict)?;

        Ok(MpasGrid {
            r,
            dtype,
            ghosts,
            loader,
            mesh,
        })
    }
}

/// Entry point per GRIDS_API §2.5.
pub fn builder() -> Builder {
    Builder::default()
}

// ------------------- Tests -------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Build a tiny closed 2-cell, 1-edge mesh covering the whole sphere
    /// (two hemispheres separated by one great-circle edge). Trivial
    /// geometry suffices for accessor-plumbing tests; §4 conformance is
    /// exercised cross-binding at pinned query points, not here.
    fn two_cell_mesh() -> MpasMeshData {
        let r = 1.0_f64;
        let half = 2.0 * PI * r * r; // 2π R² per hemisphere (total 4π R²)
        let max_edges = 1_usize;
        let input = MpasMeshInput {
            // Cell 0 at (0, π/4), cell 1 at (π, -π/4).
            lon_cell: vec![0.0, PI],
            lat_cell: vec![PI / 4.0, -PI / 4.0],
            area_cell: vec![half, half],
            n_edges_on_cell: vec![1, 1],
            cells_on_cell: vec![1, 0], // each neighbors the other
            edges_on_cell: vec![0, 0], // the one shared edge
            lon_edge: vec![0.0],
            lat_edge: vec![0.0],
            cells_on_edge: vec![0, 1],
            dc_edge: vec![PI / 2.0 * r],
            dv_edge: vec![PI / 2.0 * r],
            max_edges,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: Some(2),
            r: Some(r),
        };
        MpasMeshData::from_input(input).expect("valid 2-cell mesh")
    }

    /// 3-cell mesh around the equator with one boundary slot per cell,
    /// for accessor tests that exercise the sentinel path.
    fn three_cell_open_mesh() -> MpasMeshData {
        let r = 1.0_f64;
        let max_edges = 2_usize;
        // Three cells at lon = 0, 2π/3, 4π/3, lat = 0, each connected to
        // its two neighbors. We deliberately mark the second slot of cell
        // 0 as -1 (boundary) to test the sentinel pathway.
        let input = MpasMeshInput {
            lon_cell: vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0],
            lat_cell: vec![0.0, 0.0, 0.0],
            area_cell: vec![4.0 * PI * r * r / 3.0; 3],
            n_edges_on_cell: vec![2, 2, 2],
            // cell 0: [1, -1]; cell 1: [0, 2]; cell 2: [1, -1]
            // (asymmetric — tests that lenient mode accepts it).
            cells_on_cell: vec![1, NO_NEIGHBOR, 0, 2, 1, NO_NEIGHBOR],
            edges_on_cell: vec![0, NO_NEIGHBOR, 0, 1, 1, NO_NEIGHBOR],
            lon_edge: vec![PI / 3.0, PI],
            lat_edge: vec![0.0, 0.0],
            cells_on_edge: vec![0, 1, 1, 2],
            dc_edge: vec![1.0, 1.0],
            dv_edge: vec![1.0, 1.0],
            max_edges,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: Some(0),
            r: Some(r),
        };
        MpasMeshData::from_input(input).expect("valid open 3-cell mesh")
    }

    #[test]
    fn builder_requires_mesh_or_loader() {
        let err = builder().build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("loader_or_mesh")));
    }

    #[test]
    fn builder_rejects_nonpositive_r() {
        let err = builder().mesh(two_cell_mesh()).r(-1.0).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));

        let err = builder()
            .mesh(two_cell_mesh())
            .r(f64::NAN)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("R", _)));
    }

    #[test]
    fn builder_rejects_nonzero_ghosts() {
        let err = builder()
            .mesh(two_cell_mesh())
            .ghosts(1)
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("ghosts", _)));
    }

    #[test]
    fn loader_path_requires_reader_fn() {
        let err = builder()
            .loader(MpasLoader::new("/some/path.nc"))
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("reader_fn", _)));
    }

    #[test]
    fn invalid_loader_reader_rejected() {
        let bad = MpasLoader {
            path: "/x.nc".to_string(),
            reader: "bogus".to_string(),
            check: "strict".to_string(),
        };
        let err = builder()
            .loader(bad)
            .mesh(two_cell_mesh())
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("loader.reader", _)));
    }

    #[test]
    fn invalid_loader_check_rejected() {
        let bad = MpasLoader {
            path: "/x.nc".to_string(),
            reader: "auto".to_string(),
            check: "fuzzy".to_string(),
        };
        let err = builder()
            .loader(bad)
            .mesh(two_cell_mesh())
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("loader.check", _)));
    }

    #[test]
    fn empty_loader_path_rejected() {
        let bad = MpasLoader::new("");
        let err = builder()
            .loader(bad)
            .mesh(two_cell_mesh())
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("loader.path", _)));
    }

    #[test]
    fn defaults_match_contract() {
        // Builder-level R defaults to Earth radius regardless of the test
        // mesh's internal R — per GRIDS_API §2.2 the generator's own R is
        // the authoritative value used for provenance and conformance.
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        assert_eq!(g.radius(), DEFAULT_EARTH_RADIUS_M);
        assert_eq!(g.dtype(), Dtype::F64);
        assert_eq!(g.ghosts(), 0);
        assert_eq!(g.family(), "mpas");
        assert_eq!(g.topology(), "unstructured");
        assert_eq!(g.n_cells(), 2);
        assert_eq!(g.n_edges(), 1);
        assert_eq!(g.max_edges(), 1);
        assert!(g.loader().is_none());
    }

    #[test]
    fn reader_fn_invokes_with_loader_path() {
        let g = builder()
            .loader(MpasLoader::new("mem://fixture"))
            .reader_fn(|path| {
                assert_eq!(path, "mem://fixture");
                Ok(two_cell_mesh())
            })
            .build()
            .unwrap();
        assert_eq!(g.n_cells(), 2);
        assert_eq!(g.loader().unwrap().path, "mem://fixture");
    }

    #[test]
    fn reader_fn_error_propagates() {
        let err = builder()
            .loader(MpasLoader::new("bad://ref"))
            .reader_fn(|_| Err(GridError::Conformance("synthetic".to_string())))
            .build()
            .unwrap_err();
        assert!(matches!(err, GridError::Conformance(_)));
    }

    // --- accessors -----------------------------------------------------

    #[test]
    fn cell_center_returns_lon_lat() {
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        let (lon, lat) = g.cell_centers(0).unwrap();
        assert_eq!(lon, 0.0);
        assert_eq!(lat, PI / 4.0);
    }

    #[test]
    fn cell_center_cart_on_sphere_when_derived() {
        // Two-cell mesh omits x/y/z so defaults derive from (lon, lat, R=1).
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        for c in 0..g.n_cells() {
            let p = g.cell_center_cart(c).unwrap();
            let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                (n - 1.0).abs() < 1e-14,
                "cell {c} center off unit sphere: |p|={n}"
            );
        }
    }

    #[test]
    fn neighbors_drop_boundary_sentinel() {
        // Disable strict mode via a lenient loader so the asymmetric
        // fixture is accepted.
        let loader = MpasLoader {
            path: "test://mesh".to_string(),
            reader: "auto".to_string(),
            check: "lenient".to_string(),
        };
        let g = builder()
            .mesh(three_cell_open_mesh())
            .loader(loader)
            .build()
            .unwrap();
        assert_eq!(g.neighbors(0).unwrap(), vec![1]);
        assert_eq!(g.neighbors(1).unwrap(), vec![0, 2]);
        assert_eq!(g.neighbors(2).unwrap(), vec![1]);

        // Raw slots still expose the sentinel.
        let slots = g.neighbor_slots(0).unwrap();
        assert_eq!(slots.len(), 2);
        assert_eq!(slots[0], Some(1));
        assert_eq!(slots[1], None);
    }

    #[test]
    fn cell_out_of_range_rejected() {
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        assert!(g.cell_centers(2).is_err());
        assert!(g.neighbors(99).is_err());
        assert!(g.metric_eval(MpasMetricName::Area, 42).is_err());
    }

    #[test]
    fn edge_out_of_range_rejected() {
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        assert!(g.edge_length(1).is_err());
        assert!(g.cell_distance(99).is_err());
        assert!(g.metric_eval(MpasMetricName::DcEdge, 42).is_err());
    }

    #[test]
    fn total_area_sums_cells() {
        // Per-cell areas are stored on the mesh (not recomputed from R), so
        // compare total_area to the sum of the mesh's own area_cell.
        let mesh = two_cell_mesh();
        let expected: f64 = mesh.area_cell.iter().sum();
        let g = builder().mesh(mesh).build().unwrap();
        assert!((g.total_area() - expected).abs() < 1e-12);
    }

    #[test]
    fn metric_eval_covers_all_names() {
        let loader = MpasLoader {
            path: "test://mesh".to_string(),
            reader: "auto".to_string(),
            check: "lenient".to_string(),
        };
        let g = builder()
            .mesh(three_cell_open_mesh())
            .loader(loader)
            .build()
            .unwrap();
        assert_eq!(
            g.metric_eval(MpasMetricName::Lon, 1).unwrap(),
            2.0 * PI / 3.0
        );
        assert_eq!(g.metric_eval(MpasMetricName::Lat, 1).unwrap(), 0.0);
        assert!(g.metric_eval(MpasMetricName::Area, 0).unwrap() > 0.0);
        assert!(g.metric_eval(MpasMetricName::X, 2).unwrap().is_finite());
        assert!(g.metric_eval(MpasMetricName::Y, 2).unwrap().is_finite());
        assert!(g.metric_eval(MpasMetricName::Z, 2).unwrap().is_finite());
        assert_eq!(g.metric_eval(MpasMetricName::NEdgesOnCell, 0).unwrap(), 2.0);
        assert_eq!(g.metric_eval(MpasMetricName::LonEdge, 0).unwrap(), PI / 3.0);
        assert_eq!(g.metric_eval(MpasMetricName::LatEdge, 0).unwrap(), 0.0);
        assert_eq!(g.metric_eval(MpasMetricName::DcEdge, 0).unwrap(), 1.0);
        assert_eq!(g.metric_eval(MpasMetricName::DvEdge, 0).unwrap(), 1.0);
    }

    #[test]
    fn metric_eval_by_name_rejects_unknown() {
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        let err = g.metric_eval_by_name("not_a_metric", 0).unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("metric_name", _)));
    }

    // --- strict/lenient mesh check -------------------------------------

    #[test]
    fn strict_check_catches_asymmetry() {
        // three_cell_open_mesh is deliberately asymmetric (cell 0's
        // neighbor is 1, but cell 1's first slot points to 0 — which is
        // actually symmetric there; the boundary sentinel at cell 0 slot
        // 1 is what breaks strict-only requirements in general. Build a
        // separate fixture that breaks symmetry without a sentinel.)
        let mesh = MpasMeshData::from_input(MpasMeshInput {
            lon_cell: vec![0.0, 1.0],
            lat_cell: vec![0.0, 0.0],
            area_cell: vec![1.0, 1.0],
            n_edges_on_cell: vec![1, 1],
            // cell 0 claims neighbor 1; cell 1 claims neighbor 0 via an
            // index we'll corrupt to something that doesn't link back.
            cells_on_cell: vec![1, NO_NEIGHBOR], // cell 1's slot is -1
            edges_on_cell: vec![0, NO_NEIGHBOR],
            lon_edge: vec![0.5],
            lat_edge: vec![0.0],
            cells_on_edge: vec![0, 1],
            dc_edge: vec![1.0],
            dv_edge: vec![1.0],
            max_edges: 1,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: Some(0),
            r: Some(1.0),
        })
        .unwrap();

        // Strict build (no loader -> default strict) must fail.
        let err = builder().mesh(mesh.clone()).build().unwrap_err();
        assert!(matches!(err, GridError::SchemaViolation(_)));

        // Lenient loader accepts the mesh.
        let loader = MpasLoader {
            path: "test://mesh".to_string(),
            reader: "auto".to_string(),
            check: "lenient".to_string(),
        };
        let g = builder().mesh(mesh).loader(loader).build().unwrap();
        assert_eq!(g.n_cells(), 2);
    }

    #[test]
    fn check_mesh_catches_out_of_range_neighbor() {
        let mesh = MpasMeshData::from_input(MpasMeshInput {
            lon_cell: vec![0.0, 1.0],
            lat_cell: vec![0.0, 0.0],
            area_cell: vec![1.0, 1.0],
            n_edges_on_cell: vec![1, 1],
            // Cell 0 points to cell 99 — out of range.
            cells_on_cell: vec![99, 0],
            edges_on_cell: vec![0, 0],
            lon_edge: vec![0.0],
            lat_edge: vec![0.0],
            cells_on_edge: vec![0, 1],
            dc_edge: vec![1.0],
            dv_edge: vec![1.0],
            max_edges: 1,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: Some(0),
            r: Some(1.0),
        })
        .unwrap();
        let err = check_mesh(&mesh, false).unwrap_err();
        assert!(matches!(err, GridError::SchemaViolation(_)));
    }

    #[test]
    fn from_input_rejects_length_mismatch() {
        // lat_cell too short.
        let err = MpasMeshData::from_input(MpasMeshInput {
            lon_cell: vec![0.0, 1.0],
            lat_cell: vec![0.0], // wrong length
            area_cell: vec![1.0, 1.0],
            n_edges_on_cell: vec![0, 0],
            cells_on_cell: vec![NO_NEIGHBOR, NO_NEIGHBOR],
            edges_on_cell: vec![NO_NEIGHBOR, NO_NEIGHBOR],
            lon_edge: vec![],
            lat_edge: vec![],
            cells_on_edge: vec![],
            dc_edge: vec![],
            dv_edge: vec![],
            max_edges: 1,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: None,
            r: Some(1.0),
        })
        .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("lat_cell", _)));
    }

    #[test]
    fn from_input_rejects_zero_max_edges() {
        let err = MpasMeshData::from_input(MpasMeshInput {
            lon_cell: vec![],
            lat_cell: vec![],
            area_cell: vec![],
            n_edges_on_cell: vec![],
            cells_on_cell: vec![],
            edges_on_cell: vec![],
            lon_edge: vec![],
            lat_edge: vec![],
            cells_on_edge: vec![],
            dc_edge: vec![],
            dv_edge: vec![],
            max_edges: 0,
            x_cell: None,
            y_cell: None,
            z_cell: None,
            n_vertices: None,
            r: None,
        })
        .unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("max_edges", _)));
    }

    // --- to_esm --------------------------------------------------------

    #[test]
    fn to_esm_is_declarative() {
        let g = builder()
            .mesh(two_cell_mesh())
            .loader(MpasLoader::new("mem://fixture"))
            .build()
            .unwrap();
        let v = g.to_esm();
        assert_eq!(v["family"], "mpas");
        assert_eq!(v["version"], MPAS_FAMILY_VERSION);
        assert_eq!(v["dtype"], "float64");
        assert_eq!(v["topology"], "unstructured");
        assert_eq!(v["ghosts"], 0);
        assert_eq!(v["n_cells"], 2);
        assert_eq!(v["n_edges"], 1);
        assert_eq!(v["options"]["loader"]["path"], "mem://fixture");
        assert_eq!(v["options"]["loader"]["reader"], "auto");
        assert_eq!(v["options"]["loader"]["check"], "strict");
        assert_eq!(v["provenance"]["binding"], "rust");
        assert_eq!(v["provenance"]["family"], "mpas");

        // No inline geometry arrays per mayor's 2026-04-20 correction.
        let wire = serde_json::to_string(&v).unwrap();
        assert!(!wire.contains("lon_cell"));
        assert!(!wire.contains("area_cell"));
        assert!(!wire.contains("cells_on_cell"));
    }

    #[test]
    fn to_esm_loader_null_when_mesh_only() {
        let g = builder().mesh(two_cell_mesh()).build().unwrap();
        let v = g.to_esm();
        assert!(v["options"]["loader"].is_null());
        assert!(v["provenance"]["loader"].is_null());
    }

    #[test]
    fn dtype_f32_propagates_to_esm() {
        let g = builder()
            .mesh(two_cell_mesh())
            .dtype(Dtype::F32)
            .build()
            .unwrap();
        assert_eq!(g.dtype(), Dtype::F32);
        assert_eq!(g.to_esm()["dtype"], "float32");
    }
}
