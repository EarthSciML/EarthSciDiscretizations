//! DUO icosahedral primal-cell geometry via the landed ESS AST evaluator
//! (`eval_coeff`, `crate::rule_eval`).
//!
//! Rust mirror of `src/primal_geometry_faq.jl` (D2a / esd-heg.6). Declarative
//! companion: `discretizations/grids/duo/faq/primal_geometry.esm`. That document
//! expresses the per-cell geometry of the DUO builder — the cell-center centroid
//! (normalized mean of the three corner unit vectors), the geographic `lon`/`lat`
//! (`atan2`/`asin`), and the spherical-excess cell area (L'Huilier) — as the RFC
//! `semiring-faq-unified-ir` §8.1 geometry-FAQ pattern.
//!
//! This file is the thin **evaluation bridge**: it shapes the unit vertices +
//! faces into the per-cell scalar inputs and routes every piece of grid
//! ARITHMETIC through ESS's single AST pathway (`eval_coeff`). ESD hosts no shadow
//! evaluator (AGENTS.md single-pathway invariant).
//!
//! The transcription mirrors the imperative per-cell loop (`grids::duo`)
//! bit-for-bit: squares are products (`x*x`, not `^`); the centroid sum, the
//! squared norm, and the three side-arc dot products fold in corner / space order;
//! the L'Huilier tan-product is left-folded; `area = (4·atan(√t)) · (R·R)` keeps
//! the imperative association `a_unit * r2`. The single intentional divergence is
//! dropping the `clamp` / `max` clipping guards: on a valid icosahedral mesh every
//! `acos` / `asin` argument and the L'Huilier radicand already lie in range, so the
//! guards never fire and the bytes match (the Rust `eval_coeff` and imperative
//! paths share the same `f64` libm).

use serde_json::{json, Value};

use crate::{eval_coeff, Bindings};

/// Materialized DUO primal-cell geometry. All arrays are `f64` (the grid stores
/// internal state at `f64`; the declared `dtype` tracks precision per GRIDS_API §7).
pub(crate) struct PrimalGeometry {
    /// (3, Nc) flattened column-major R-scaled cell-center cartesian coords.
    pub cell_cart: Vec<f64>,
    /// (Nc,) cell-center longitude (rad).
    pub lon: Vec<f64>,
    /// (Nc,) cell-center latitude (rad).
    pub lat: Vec<f64>,
    /// (Nc,) spherical-triangle area (R² scaled).
    pub area: Vec<f64>,
}

// AST nodes mirroring the equation bodies of
// `discretizations/grids/duo/faq/primal_geometry.esm`. Scalar corner leaves:
// a1..a3 (corner 1 = face vertex 1), b1..b3 (corner 2), c1..c3 (corner 3), R.
fn mk(op: &str, args: Vec<Value>) -> Value {
    json!({"op": op, "args": args})
}

/// All six per-cell metric ASTs, built once per call.
struct GeomAst {
    cart_x: Value,
    cart_y: Value,
    cart_z: Value,
    lon: Value,
    lat: Value,
    area: Value,
}

fn build_geom_ast() -> GeomAst {
    // centroid: m = a+b+c per component; n = sqrt(mx*mx+my*my+mz*mz); u = m/n.
    let mx = mk("+", vec![json!("a1"), json!("b1"), json!("c1")]);
    let my = mk("+", vec![json!("a2"), json!("b2"), json!("c2")]);
    let mz = mk("+", vec![json!("a3"), json!("b3"), json!("c3")]);
    let nsq = mk(
        "+",
        vec![
            mk("*", vec![mx.clone(), mx.clone()]),
            mk("*", vec![my.clone(), my.clone()]),
            mk("*", vec![mz.clone(), mz.clone()]),
        ],
    ); // squares as *, not ^
    let nrm = mk("sqrt", vec![nsq]);
    let ux = mk("/", vec![mx.clone(), nrm.clone()]);
    let uy = mk("/", vec![my.clone(), nrm.clone()]);
    let uz = mk("/", vec![mz.clone(), nrm]);

    // cell-center cartesian (R-scaled), geographic lon/lat of the unit centroid.
    let cart_x = mk("*", vec![json!("R"), ux.clone()]);
    let cart_y = mk("*", vec![json!("R"), uy.clone()]);
    let cart_z = mk("*", vec![json!("R"), uz.clone()]);
    let lon = mk("atan2", vec![uy, ux]);
    let lat = mk("asin", vec![uz]); // no clamp

    // area (L'Huilier): da=acos(b·c), db=acos(c·a), dc=acos(a·b); s=(da+db+dc)/2;
    // t=tan(s/2)·tan((s-da)/2)·tan((s-db)/2)·tan((s-dc)/2); area=(4·atan(√t))·(R·R).
    let dot_bc = mk(
        "+",
        vec![
            mk("*", vec![json!("b1"), json!("c1")]),
            mk("*", vec![json!("b2"), json!("c2")]),
            mk("*", vec![json!("b3"), json!("c3")]),
        ],
    );
    let dot_ca = mk(
        "+",
        vec![
            mk("*", vec![json!("c1"), json!("a1")]),
            mk("*", vec![json!("c2"), json!("a2")]),
            mk("*", vec![json!("c3"), json!("a3")]),
        ],
    );
    let dot_ab = mk(
        "+",
        vec![
            mk("*", vec![json!("a1"), json!("b1")]),
            mk("*", vec![json!("a2"), json!("b2")]),
            mk("*", vec![json!("a3"), json!("b3")]),
        ],
    );
    let da = mk("acos", vec![dot_bc]);
    let db = mk("acos", vec![dot_ca]);
    let dc = mk("acos", vec![dot_ab]);
    let s = mk(
        "*",
        vec![
            json!(0.5),
            mk("+", vec![da.clone(), db.clone(), dc.clone()]),
        ],
    );
    let half = |x: Value| mk("/", vec![x, json!(2.0)]);
    let t = mk(
        "*",
        vec![
            mk("tan", vec![half(s.clone())]),
            mk("tan", vec![half(mk("-", vec![s.clone(), da]))]),
            mk("tan", vec![half(mk("-", vec![s.clone(), db]))]),
            mk("tan", vec![half(mk("-", vec![s, dc]))]),
        ],
    );
    let area = mk(
        "*",
        vec![
            mk("*", vec![json!(4.0), mk("atan", vec![mk("sqrt", vec![t])])]),
            mk("*", vec![json!("R"), json!("R")]),
        ],
    ); // no max guard

    GeomAst {
        cart_x,
        cart_y,
        cart_z,
        lon,
        lat,
        area,
    }
}

/// Materialize the DUO primal-cell geometry from the declarative FAQ, routing all
/// arithmetic through the landed ESS AST evaluator (`eval_coeff`).
///
/// `verts_unit` are `(Nv)` cartesian **unit-sphere** vertex coords (NOT R-scaled),
/// as produced by `duo_subdivide_faq`; `faces` are `(Nc)` 0-based vertex-index
/// triples per cell; `r` is the sphere radius. Byte-identical (`f64`) to the
/// imperative cell-geometry loop (`grids::duo`).
pub(crate) fn primal_geometry_faq(
    verts_unit: &[[f64; 3]],
    faces: &[[u32; 3]],
    r: f64,
) -> PrimalGeometry {
    let ast = build_geom_ast();
    let nc = faces.len();
    let mut cell_cart = vec![0.0_f64; 3 * nc];
    let mut lon = vec![0.0_f64; nc];
    let mut lat = vec![0.0_f64; nc];
    let mut area = vec![0.0_f64; nc];

    for (c, f) in faces.iter().enumerate() {
        let av = verts_unit[f[0] as usize];
        let bv = verts_unit[f[1] as usize];
        let cv = verts_unit[f[2] as usize];
        // One bindings map per cell drives every co-derived metric.
        let mut b = Bindings::new();
        b.insert("a1".to_string(), av[0]);
        b.insert("a2".to_string(), av[1]);
        b.insert("a3".to_string(), av[2]);
        b.insert("b1".to_string(), bv[0]);
        b.insert("b2".to_string(), bv[1]);
        b.insert("b3".to_string(), bv[2]);
        b.insert("c1".to_string(), cv[0]);
        b.insert("c2".to_string(), cv[1]);
        b.insert("c3".to_string(), cv[2]);
        b.insert("R".to_string(), r);

        cell_cart[3 * c] = eval_coeff(&ast.cart_x, &b).expect("duo geom faq: cart_x");
        cell_cart[3 * c + 1] = eval_coeff(&ast.cart_y, &b).expect("duo geom faq: cart_y");
        cell_cart[3 * c + 2] = eval_coeff(&ast.cart_z, &b).expect("duo geom faq: cart_z");
        lon[c] = eval_coeff(&ast.lon, &b).expect("duo geom faq: lon");
        lat[c] = eval_coeff(&ast.lat, &b).expect("duo geom faq: lat");
        area[c] = eval_coeff(&ast.area, &b).expect("duo geom faq: area");
    }

    PrimalGeometry {
        cell_cart,
        lon,
        lat,
        area,
    }
}
