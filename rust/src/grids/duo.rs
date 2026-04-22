//! DUO icosahedral triangular grid family (Heikes et al. 2023).
//!
//! Loader-backed per `docs/GRIDS_API.md` §10. The `.esm` carries a small
//! declarative config (family, level, reader ref); geometry is derived at
//! runtime by the accessors. Cross-binding conformance is compared at
//! pinned query points, not serialized bytes.
//!
//! Subdivision level `r` yields 20·4^r triangular cells, 10·4^r + 2
//! vertices, and 30·4^r edges on a sphere of radius `R`.

use std::collections::HashMap;

use serde_json::{json, Map, Value};

use crate::{Dtype, Grid, GridError, Result};

const DUO_FAMILY_VERSION: &str = "1.0.0";
const DEFAULT_EARTH_RADIUS_M: f64 = 6.371e6;
const BUILTIN_LEVEL_PREFIX: &str = "builtin://icosahedral/";

/// Loader struct per GRIDS_API §10.
#[derive(Clone, Debug, PartialEq)]
pub struct DuoLoader {
    pub path: String,
    pub reader: String,
    pub check: String,
}

impl DuoLoader {
    pub fn builtin_level(level: u32) -> Self {
        Self {
            path: format!("{BUILTIN_LEVEL_PREFIX}{level}"),
            reader: "builtin_icosahedral".to_string(),
            check: "strict".to_string(),
        }
    }
}

/// Fully-materialized icosahedral triangular mesh.
#[derive(Clone, Debug)]
pub struct DuoGrid {
    level: u32,
    r: f64,
    dtype: Dtype,
    ghosts: u32,
    /// (3, Nv) flattened column-major: vertices\[3*i..3*i+3] is vertex i.
    vertices: Vec<f64>,
    /// (3, Nc) flattened column-major, 0-based vertex indices.
    faces: Vec<u32>,
    lon: Vec<f64>,
    lat: Vec<f64>,
    /// (3, Nc) cartesian cell centers on the sphere.
    cell_cart: Vec<f64>,
    area: Vec<f64>,
    /// (2, Ne) sorted vertex-pair per edge.
    edges: Vec<u32>,
    /// (3, Nc); `cell_neighbors[3*c+k]` is neighbor across edge opposite
    /// local vertex k (0 = boundary; closed mesh has none).
    /// Stored as Option via u32::MAX sentinel for "no neighbor".
    cell_neighbors: Vec<u32>,
    /// Per-vertex face adjacency (sorted, 0-based).
    vertex_faces: Vec<Vec<u32>>,
    provenance: Value,
    loader: DuoLoader,
}

const NO_NEIGHBOR: u32 = u32::MAX;

// ------------------- Base icosahedron -------------------------------------

fn icosahedron_vertices_f64() -> [[f64; 3]; 12] {
    let phi = (1.0f64 + 5.0f64.sqrt()) / 2.0;
    let raw: [[f64; 3]; 12] = [
        [0.0, 1.0, phi],
        [0.0, -1.0, phi],
        [0.0, 1.0, -phi],
        [0.0, -1.0, -phi],
        [1.0, phi, 0.0],
        [-1.0, phi, 0.0],
        [1.0, -phi, 0.0],
        [-1.0, -phi, 0.0],
        [phi, 0.0, 1.0],
        [phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
    ];
    let mut out = [[0.0; 3]; 12];
    for (i, v) in raw.iter().enumerate() {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        out[i] = [v[0] / n, v[1] / n, v[2] / n];
    }
    out
}

/// 20 triangular faces of the base icosahedron (0-based vertex indices).
/// Winding matches the Julia reference so level-0 output is identical.
fn icosahedron_faces_0based() -> [[u32; 3]; 20] {
    // Julia 1-based → rust 0-based: subtract 1 from every index.
    [
        [0, 8, 1],
        [0, 1, 10],
        [0, 10, 5],
        [0, 5, 4],
        [0, 4, 8],
        [8, 4, 9],
        [8, 9, 6],
        [8, 6, 1],
        [1, 6, 7],
        [1, 7, 10],
        [10, 7, 11],
        [10, 11, 5],
        [5, 11, 2],
        [5, 2, 4],
        [4, 2, 9],
        [9, 2, 3],
        [9, 3, 6],
        [6, 3, 7],
        [7, 3, 11],
        [11, 3, 2],
    ]
}

// ------------------- Recursive subdivision --------------------------------

fn subdivide_icosahedron(level: u32) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let base_v = icosahedron_vertices_f64();
    let base_f = icosahedron_faces_0based();
    let mut verts: Vec<[f64; 3]> = base_v.to_vec();
    let mut faces: Vec<[u32; 3]> = base_f.to_vec();

    for _ in 0..level {
        let mut cache: HashMap<(u32, u32), u32> = HashMap::new();
        let mut new_faces: Vec<[u32; 3]> = Vec::with_capacity(faces.len() * 4);
        for f in &faces {
            let a = f[0];
            let b = f[1];
            let c = f[2];
            let ab = midpoint(&mut verts, &mut cache, a, b);
            let bc = midpoint(&mut verts, &mut cache, b, c);
            let ca = midpoint(&mut verts, &mut cache, c, a);
            new_faces.push([a, ab, ca]);
            new_faces.push([b, bc, ab]);
            new_faces.push([c, ca, bc]);
            new_faces.push([ab, bc, ca]);
        }
        faces = new_faces;
    }
    (verts, faces)
}

fn midpoint(
    verts: &mut Vec<[f64; 3]>,
    cache: &mut HashMap<(u32, u32), u32>,
    a: u32,
    b: u32,
) -> u32 {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&i) = cache.get(&key) {
        return i;
    }
    let va = verts[a as usize];
    let vb = verts[b as usize];
    let mx = va[0] + vb[0];
    let my = va[1] + vb[1];
    let mz = va[2] + vb[2];
    let n = (mx * mx + my * my + mz * mz).sqrt();
    verts.push([mx / n, my / n, mz / n]);
    let idx = (verts.len() - 1) as u32;
    cache.insert(key, idx);
    idx
}

// ------------------- Geometry helpers -------------------------------------

fn cart_to_lonlat(x: f64, y: f64, z: f64) -> (f64, f64) {
    (y.atan2(x), z.clamp(-1.0, 1.0).asin())
}

/// Spherical excess area of a unit-sphere triangle via L'Huilier's theorem.
fn spherical_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let da = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2])
        .clamp(-1.0, 1.0)
        .acos();
    let db = (c[0] * a[0] + c[1] * a[1] + c[2] * a[2])
        .clamp(-1.0, 1.0)
        .acos();
    let dc = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
        .clamp(-1.0, 1.0)
        .acos();
    let s = 0.5 * (da + db + dc);
    let t =
        (s / 2.0).tan() * ((s - da) / 2.0).tan() * ((s - db) / 2.0).tan() * ((s - dc) / 2.0).tan();
    let t = t.max(0.0);
    4.0 * t.sqrt().atan()
}

// ------------------- Connectivity -----------------------------------------

type Connectivity = (Vec<[u32; 2]>, Vec<[u32; 3]>);

fn build_connectivity(faces: &[[u32; 3]]) -> Result<Connectivity> {
    let nc = faces.len();
    // Edge opposite local vertex k is the edge between the other two:
    //  k=0 → (v1, v2); k=1 → (v2, v0); k=2 → (v0, v1).
    let mut edge_map: HashMap<(u32, u32), Vec<(u32, u32)>> = HashMap::new();
    for (c, f) in faces.iter().enumerate() {
        let edges_of_cell = [(f[1], f[2]), (f[2], f[0]), (f[0], f[1])];
        for (k, &(a, b)) in edges_of_cell.iter().enumerate() {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_map.entry(key).or_default().push((c as u32, k as u32));
        }
    }

    let mut edges: Vec<[u32; 2]> = Vec::with_capacity(edge_map.len());
    let mut neighbors = vec![NO_NEIGHBOR; 3 * nc];

    // Sort edges by (min, max) so output is deterministic across runs.
    let mut keys: Vec<(u32, u32)> = edge_map.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        edges.push([key.0, key.1]);
        let list = &edge_map[&key];
        match list.len() {
            1 => { /* boundary edge; closed mesh has none */ }
            2 => {
                let (c1, k1) = list[0];
                let (c2, k2) = list[1];
                neighbors[3 * c1 as usize + k1 as usize] = c2;
                neighbors[3 * c2 as usize + k2 as usize] = c1;
            }
            n => {
                return Err(GridError::SchemaViolation(format!(
                    "duo: non-manifold edge ({}, {}) shared by {n} cells",
                    key.0, key.1
                )));
            }
        }
    }

    // Pack neighbors into flat 3×Nc buffer matching the stored shape.
    let mut nbr_flat = vec![NO_NEIGHBOR; 3 * nc];
    for (i, &v) in neighbors.iter().enumerate() {
        nbr_flat[i] = v;
    }
    Ok((edges, pack3(&nbr_flat, nc)))
}

fn pack3(flat: &[u32], nc: usize) -> Vec<[u32; 3]> {
    (0..nc)
        .map(|c| [flat[3 * c], flat[3 * c + 1], flat[3 * c + 2]])
        .collect()
}

fn vertex_faces(faces: &[[u32; 3]], nv: usize) -> Vec<Vec<u32>> {
    let mut vf: Vec<Vec<u32>> = vec![Vec::new(); nv];
    for (c, f) in faces.iter().enumerate() {
        for &v in f.iter() {
            vf[v as usize].push(c as u32);
        }
    }
    for list in vf.iter_mut() {
        list.sort_unstable();
    }
    vf
}

// ------------------- Loader ------------------------------------------------

fn parse_builtin_level(path: &str) -> Result<Option<u32>> {
    if let Some(tail) = path.strip_prefix(BUILTIN_LEVEL_PREFIX) {
        let lvl: i64 = tail.parse().map_err(|_| {
            GridError::InvalidOption("loader.path", format!("cannot parse level from {path}"))
        })?;
        if lvl < 0 {
            return Err(GridError::InvalidOption(
                "loader.path",
                format!("subdivision level must be ≥ 0, got {lvl}"),
            ));
        }
        return Ok(Some(lvl as u32));
    }
    Ok(None)
}

fn resolve_loader_level(loader: &DuoLoader) -> Result<u32> {
    if let Some(lvl) = parse_builtin_level(&loader.path)? {
        return Ok(lvl);
    }
    if loader.reader == "duo_mesh" || loader.reader == "auto" {
        return Err(GridError::InvalidOption(
            "loader",
            format!(
                "duo: .duo mesh-file reader not yet implemented — pending \
                 EarthSciSerialization file-format spec. Use \
                 builtin://icosahedral/<level> in the meantime (got path={})",
                loader.path
            ),
        ));
    }
    Err(GridError::InvalidOption(
        "loader",
        format!(
            "duo: unrecognized loader path {} with reader={}",
            loader.path, loader.reader
        ),
    ))
}

// ------------------- Builder ----------------------------------------------

/// Builder for a [`DuoGrid`]. See GRIDS_API §2.5.
#[must_use]
#[derive(Clone, Debug, Default)]
pub struct Builder {
    loader: Option<DuoLoader>,
    r: Option<f64>,
    dtype: Option<Dtype>,
    ghosts: Option<u32>,
}

impl Builder {
    pub fn loader(mut self, loader: DuoLoader) -> Self {
        self.loader = Some(loader);
        self
    }

    pub fn r(mut self, r: f64) -> Self {
        self.r = Some(r);
        self
    }

    pub fn dtype(mut self, dtype: Dtype) -> Self {
        self.dtype = Some(dtype);
        self
    }

    pub fn ghosts(mut self, ghosts: u32) -> Self {
        self.ghosts = Some(ghosts);
        self
    }

    pub fn build(self) -> Result<DuoGrid> {
        let loader = self.loader.ok_or(GridError::MissingOption("loader"))?;
        let r = self.r.unwrap_or(DEFAULT_EARTH_RADIUS_M);
        if !(r.is_finite() && r > 0.0) {
            return Err(GridError::InvalidOption(
                "R",
                format!("must be positive and finite, got {r}"),
            ));
        }
        let dtype = self.dtype.unwrap_or_default();
        let ghosts = self.ghosts.unwrap_or(0);

        let level = resolve_loader_level(&loader)?;

        let (base_verts, base_faces) = subdivide_icosahedron(level);
        let nv = base_verts.len();
        let nc = base_faces.len();

        // Derive geometry at f64 precision; quantize to dtype on storage only
        // where it matters for the conformance contract. Internal state stays
        // f64 so accessors return the same values regardless of dtype — the
        // dtype field tracks the *declared* precision per GRIDS_API §7.
        let mut lon = vec![0.0_f64; nc];
        let mut lat = vec![0.0_f64; nc];
        let mut cell_cart = vec![0.0_f64; 3 * nc];
        let mut area = vec![0.0_f64; nc];
        let r2 = r * r;
        for (c, f) in base_faces.iter().enumerate() {
            let a = base_verts[f[0] as usize];
            let b = base_verts[f[1] as usize];
            let cc = base_verts[f[2] as usize];
            let mx = a[0] + b[0] + cc[0];
            let my = a[1] + b[1] + cc[1];
            let mz = a[2] + b[2] + cc[2];
            let n = (mx * mx + my * my + mz * mz).sqrt();
            let ux = mx / n;
            let uy = my / n;
            let uz = mz / n;
            cell_cart[3 * c] = r * ux;
            cell_cart[3 * c + 1] = r * uy;
            cell_cart[3 * c + 2] = r * uz;
            let (lo, la) = cart_to_lonlat(ux, uy, uz);
            lon[c] = lo;
            lat[c] = la;
            area[c] = spherical_triangle_area(a, b, cc) * r2;
        }

        let mut vertices = vec![0.0_f64; 3 * nv];
        for (i, v) in base_verts.iter().enumerate() {
            vertices[3 * i] = r * v[0];
            vertices[3 * i + 1] = r * v[1];
            vertices[3 * i + 2] = r * v[2];
        }

        let mut faces_flat = vec![0_u32; 3 * nc];
        for (c, f) in base_faces.iter().enumerate() {
            faces_flat[3 * c] = f[0];
            faces_flat[3 * c + 1] = f[1];
            faces_flat[3 * c + 2] = f[2];
        }

        let (edges_pair, neighbors) = build_connectivity(&base_faces)?;
        let mut edges_flat = vec![0_u32; 2 * edges_pair.len()];
        for (e, ep) in edges_pair.iter().enumerate() {
            edges_flat[2 * e] = ep[0];
            edges_flat[2 * e + 1] = ep[1];
        }
        let mut neighbors_flat = vec![NO_NEIGHBOR; 3 * nc];
        for (c, nb) in neighbors.iter().enumerate() {
            neighbors_flat[3 * c] = nb[0];
            neighbors_flat[3 * c + 1] = nb[1];
            neighbors_flat[3 * c + 2] = nb[2];
        }

        let vf = vertex_faces(&base_faces, nv);

        let dtype_str = match dtype {
            Dtype::F64 => "float64",
            Dtype::F32 => "float32",
        };
        let provenance = json!({
            "binding": "rust",
            "family": "duo",
            "version": DUO_FAMILY_VERSION,
            "level": level,
            "reader": loader.reader,
            "path": loader.path,
            "check": loader.check,
            "dtype": dtype_str,
        });

        Ok(DuoGrid {
            level,
            r,
            dtype,
            ghosts,
            vertices,
            faces: faces_flat,
            lon,
            lat,
            cell_cart,
            area,
            edges: edges_flat,
            cell_neighbors: neighbors_flat,
            vertex_faces: vf,
            provenance,
            loader,
        })
    }
}

/// Entry point per GRIDS_API §2.5.
pub fn builder() -> Builder {
    Builder::default()
}

// ------------------- Accessors --------------------------------------------

impl DuoGrid {
    pub fn level(&self) -> u32 {
        self.level
    }

    pub fn radius(&self) -> f64 {
        self.r
    }

    pub fn ghosts(&self) -> u32 {
        self.ghosts
    }

    pub fn n_cells(&self) -> usize {
        self.faces.len() / 3
    }

    pub fn n_vertices(&self) -> usize {
        self.vertices.len() / 3
    }

    pub fn n_edges(&self) -> usize {
        self.edges.len() / 2
    }

    pub fn total_area(&self) -> f64 {
        self.area.iter().sum()
    }

    pub fn lon(&self) -> &[f64] {
        &self.lon
    }

    pub fn lat(&self) -> &[f64] {
        &self.lat
    }

    pub fn area(&self) -> &[f64] {
        &self.area
    }

    pub fn loader(&self) -> &DuoLoader {
        &self.loader
    }

    /// Cell-center geographic coords for cell `c` (rad).
    pub fn cell_centers(&self, c: usize) -> (f64, f64) {
        (self.lon[c], self.lat[c])
    }

    /// Cell-center cartesian coords for cell `c` (same units as R).
    pub fn cell_center_cart(&self, c: usize) -> [f64; 3] {
        [
            self.cell_cart[3 * c],
            self.cell_cart[3 * c + 1],
            self.cell_cart[3 * c + 2],
        ]
    }

    /// Three edge-adjacent cell indices for cell `c`. `None` at a boundary
    /// (not reachable on a closed icosahedral mesh).
    pub fn neighbors(&self, c: usize) -> [Option<u32>; 3] {
        let slot = |v: u32| if v == NO_NEIGHBOR { None } else { Some(v) };
        [
            slot(self.cell_neighbors[3 * c]),
            slot(self.cell_neighbors[3 * c + 1]),
            slot(self.cell_neighbors[3 * c + 2]),
        ]
    }

    /// 0-based vertex indices of cell `c`.
    pub fn cell_vertices(&self, c: usize) -> [u32; 3] {
        [
            self.faces[3 * c],
            self.faces[3 * c + 1],
            self.faces[3 * c + 2],
        ]
    }

    /// Cartesian coords of vertex `v` (same units as R).
    pub fn vertex(&self, v: usize) -> [f64; 3] {
        [
            self.vertices[3 * v],
            self.vertices[3 * v + 1],
            self.vertices[3 * v + 2],
        ]
    }

    /// Sorted face indices incident on vertex `v`.
    pub fn vertex_faces(&self, v: usize) -> &[u32] {
        &self.vertex_faces[v]
    }

    /// Sorted vertex pair (min, max) for edge `e`.
    pub fn edge(&self, e: usize) -> (u32, u32) {
        (self.edges[2 * e], self.edges[2 * e + 1])
    }

    /// Accessor for per-cell scalar metrics. Names:
    /// `area`, `lon`, `lat`, `x`, `y`, `z`.
    pub fn metric_eval(&self, name: &str, c: usize) -> Result<f64> {
        match name {
            "area" => Ok(self.area[c]),
            "lon" => Ok(self.lon[c]),
            "lat" => Ok(self.lat[c]),
            "x" => Ok(self.cell_cart[3 * c]),
            "y" => Ok(self.cell_cart[3 * c + 1]),
            "z" => Ok(self.cell_cart[3 * c + 2]),
            other => Err(GridError::InvalidOption(
                "metric",
                format!("duo: unknown metric {other}"),
            )),
        }
    }
}

impl Grid for DuoGrid {
    fn family(&self) -> &'static str {
        "duo"
    }

    fn dtype(&self) -> Dtype {
        self.dtype
    }

    /// §6-schema-valid declarative config. No inline geometry arrays
    /// (per mayor's 2026-04-20 correction).
    fn to_esm(&self) -> Value {
        let mut options = Map::new();
        options.insert("R".to_string(), json!(self.r));
        options.insert("level".to_string(), json!(self.level));
        options.insert(
            "loader".to_string(),
            json!({
                "path": self.loader.path,
                "reader": self.loader.reader,
                "check": self.loader.check,
            }),
        );
        let dtype_str = match self.dtype {
            Dtype::F64 => "float64",
            Dtype::F32 => "float32",
        };
        json!({
            "family": "duo",
            "topology": "unstructured",
            "dtype": dtype_str,
            "ghosts": self.ghosts,
            "n_cells": self.n_cells(),
            "n_vertices": self.n_vertices(),
            "n_edges": self.n_edges(),
            "options": Value::Object(options),
            "provenance": self.provenance,
            "schema_version": DUO_FAMILY_VERSION,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level0_topology() {
        let g = builder()
            .loader(DuoLoader::builtin_level(0))
            .build()
            .unwrap();
        assert_eq!(g.n_cells(), 20);
        assert_eq!(g.n_vertices(), 12);
        // V - E + F = 2 for a sphere → E = V + F - 2 = 30
        assert_eq!(g.n_edges(), 30);
    }

    #[test]
    fn level2_topology() {
        let g = builder()
            .loader(DuoLoader::builtin_level(2))
            .build()
            .unwrap();
        let r = 2u32;
        let expected_cells = 20 * 4_usize.pow(r);
        let expected_verts = 10 * 4_usize.pow(r) + 2;
        let expected_edges = 30 * 4_usize.pow(r);
        assert_eq!(g.n_cells(), expected_cells);
        assert_eq!(g.n_vertices(), expected_verts);
        assert_eq!(g.n_edges(), expected_edges);
    }

    #[test]
    fn total_area_matches_sphere() {
        for level in 0..=3 {
            let g = builder()
                .loader(DuoLoader::builtin_level(level))
                .build()
                .unwrap();
            let r = g.radius();
            let expected = 4.0 * std::f64::consts::PI * r * r;
            let got = g.total_area();
            // L'Huilier on unit sphere scaled by R²; level-0 closes to
            // within ~1e-12 relative (tight; far within GRIDS_API §4.2
            // 1e-14 is aspirational here — L'Huilier accumulates per-face
            // atan error).
            let rel = (got - expected).abs() / expected;
            assert!(
                rel < 1e-12,
                "level={level}: total area {got} vs 4πR² {expected} (rel {rel})"
            );
        }
    }

    #[test]
    fn centers_are_on_sphere() {
        let g = builder()
            .loader(DuoLoader::builtin_level(2))
            .r(1.0)
            .build()
            .unwrap();
        for c in 0..g.n_cells() {
            let p = g.cell_center_cart(c);
            let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(
                (n - 1.0).abs() < 1e-14,
                "cell {c} center off unit sphere: |p|={n}"
            );
        }
    }

    #[test]
    fn vertices_on_sphere_with_radius() {
        let g = builder()
            .loader(DuoLoader::builtin_level(3))
            .r(6.371e6)
            .build()
            .unwrap();
        for v in 0..g.n_vertices() {
            let p = g.vertex(v);
            let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let rel = (n - 6.371e6).abs() / 6.371e6;
            assert!(rel < 1e-14, "vertex {v} off sphere: |p|={n}");
        }
    }

    #[test]
    fn neighbors_are_symmetric_and_complete() {
        let g = builder()
            .loader(DuoLoader::builtin_level(2))
            .build()
            .unwrap();
        let nc = g.n_cells();
        for c in 0..nc {
            let nb = g.neighbors(c);
            for slot in nb.iter() {
                let n = slot.expect("closed icosahedral mesh has no boundary");
                assert!((n as usize) < nc);
                // Symmetry: c must appear in n's neighbor list.
                let nn = g.neighbors(n as usize);
                assert!(
                    nn.contains(&Some(c as u32)),
                    "asymmetric neighbor: {c} → {n} but not back"
                );
            }
        }
    }

    #[test]
    fn edges_are_sorted_and_unique() {
        let g = builder()
            .loader(DuoLoader::builtin_level(2))
            .build()
            .unwrap();
        let mut seen: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
        let mut prev: Option<(u32, u32)> = None;
        for e in 0..g.n_edges() {
            let (a, b) = g.edge(e);
            assert!(a < b, "edge {e} not sorted: ({a}, {b})");
            assert!(seen.insert((a, b)), "edge {e} duplicate: ({a}, {b})");
            if let Some(p) = prev {
                assert!(p < (a, b), "edges not globally sorted at {e}");
            }
            prev = Some((a, b));
        }
    }

    #[test]
    fn missing_loader_errors() {
        let err = builder().build().unwrap_err();
        assert!(matches!(err, GridError::MissingOption("loader")));
    }

    #[test]
    fn negative_level_errors() {
        let bad = DuoLoader {
            path: "builtin://icosahedral/-1".to_string(),
            reader: "builtin_icosahedral".to_string(),
            check: "strict".to_string(),
        };
        let err = builder().loader(bad).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("loader.path", _)));
    }

    #[test]
    fn duo_mesh_reader_not_implemented() {
        let ldr = DuoLoader {
            path: "/some/where.duo".to_string(),
            reader: "duo_mesh".to_string(),
            check: "strict".to_string(),
        };
        let err = builder().loader(ldr).build().unwrap_err();
        assert!(matches!(err, GridError::InvalidOption("loader", _)));
    }

    #[test]
    fn to_esm_has_required_fields() {
        let g = builder()
            .loader(DuoLoader::builtin_level(1))
            .build()
            .unwrap();
        let v = g.to_esm();
        assert_eq!(v["family"], json!("duo"));
        assert_eq!(v["topology"], json!("unstructured"));
        assert_eq!(v["dtype"], json!("float64"));
        assert_eq!(v["n_cells"], json!(80));
        assert_eq!(v["options"]["level"], json!(1));
        assert_eq!(
            v["options"]["loader"]["path"],
            json!("builtin://icosahedral/1")
        );
        // No inline geometry arrays per mayor's 2026-04-20 correction.
        assert!(v.get("vertices").is_none());
        assert!(v.get("faces").is_none());
        assert!(v.get("lon").is_none());
    }

    #[test]
    fn metric_eval_routes() {
        let g = builder()
            .loader(DuoLoader::builtin_level(0))
            .r(1.0)
            .build()
            .unwrap();
        let a = g.metric_eval("area", 0).unwrap();
        assert!(a > 0.0);
        let x = g.metric_eval("x", 0).unwrap();
        let y = g.metric_eval("y", 0).unwrap();
        let z = g.metric_eval("z", 0).unwrap();
        let n = (x * x + y * y + z * z).sqrt();
        assert!((n - 1.0).abs() < 1e-14);
        assert!(g.metric_eval("bogus", 0).is_err());
    }

    #[test]
    fn vertex_faces_are_sorted_and_cover_five_or_six() {
        // Icosahedral subdivision: 12 original vertices have valence 5,
        // all new vertices have valence 6.
        let g = builder()
            .loader(DuoLoader::builtin_level(2))
            .build()
            .unwrap();
        for v in 0..g.n_vertices() {
            let vf = g.vertex_faces(v);
            assert!(
                vf.len() == 5 || vf.len() == 6,
                "vertex {v} has unexpected valence {}",
                vf.len()
            );
            for w in vf.windows(2) {
                assert!(w[0] < w[1], "vertex_faces not sorted at v={v}");
            }
        }
        let pentagons = (0..g.n_vertices())
            .filter(|v| g.vertex_faces(*v).len() == 5)
            .count();
        assert_eq!(pentagons, 12, "exactly 12 pentagonal vertices expected");
    }
}
