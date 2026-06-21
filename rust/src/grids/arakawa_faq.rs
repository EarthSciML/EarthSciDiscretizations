//! Arakawa staggered-grid construction via the landed M1 elementwise FAQ path
//! (`eval_coeff`, `crate::rule_eval`).
//!
//! Rust mirror of `src/arakawa_faq.jl` (esd-3we, S4). Declarative companion:
//! `discretizations/grids/arakawa/rules/arakawa_construction.esm`. The
//! construction is a semiring-FAQ (RFC `semiring-faq-unified-ir` §5.1/§5.2:
//! cartesian interval index-sets + elementwise `arrayop` bodies — the M1 path, NO
//! value-invention, NO recursion).
//!
//! Arakawa staggering over a cartesian base is the cartesian product of two
//! staggered 1-D affine maps per axis: a **center** map `center(i) = lo +
//! (i-0.5)*d` and a **node** map `node(i) = lo + (i-1)*d`. The four staggered
//! locations are pure structural gathers of those two maps (row-major). This
//! module routes every piece of ARITHMETIC — `dx`/`dy`, the two axis maps, the
//! cell-volume product — through ESS's AST pathway (`eval_coeff`); the structural
//! arrays (the cartesian-product gather into the four locations, the row-major
//! neighbor linearization with the `-1` off-grid sentinel, the boundary masks,
//! and the static A–E variable-location table) are pure index logic.
//! Byte-identical to the imperative accessors and to the committed
//! Julia-reference `tests/conformance/grids/arakawa/construction/golden.json`.
//!
//! Takes `(&CartesianBase, Stagger)` rather than the `ArakawaGrid` (whose base is
//! a `Box<dyn BaseGrid>` that does not expose the raw extents the affine maps
//! need); the geometry/neighbors/boundary are stagger-independent and the
//! variable-location table is the only stagger-dependent output.

use serde_json::{json, Value};

use super::arakawa::{BaseGrid, CartesianBase, Location, Stagger};
use crate::{eval_coeff, Bindings};

/// Materialized arakawa construction. Per-cell arrays are over the `(nx, ny)`
/// cell-center set (row-major); each location's coordinate pair is flattened
/// row-major over that location's shape. Mirrors the `construction/golden.json`
/// contract.
pub struct ArakawaConstruction {
    /// Base interior cell counts.
    pub nx: usize,
    /// Base interior cell counts.
    pub ny: usize,
    /// Affine x spacing.
    pub dx: f64,
    /// Affine y spacing.
    pub dy: f64,
    /// Cell measure `dx*dy` (uniform).
    pub cell_volume: f64,
    /// Cell-center `(xs, ys)` flattened row-major over `(nx, ny)`.
    pub coords_cell_center: (Vec<f64>, Vec<f64>),
    /// U-edge `(xs, ys)` over `(nx+1, ny)`.
    pub coords_u_edge: (Vec<f64>, Vec<f64>),
    /// V-edge `(xs, ys)` over `(nx, ny+1)`.
    pub coords_v_edge: (Vec<f64>, Vec<f64>),
    /// Corner `(xs, ys)` over `(nx+1, ny+1)`.
    pub coords_corner: (Vec<f64>, Vec<f64>),
    /// X −1 neighbor over the cell-center set (`-1` off-grid).
    pub neighbor_x_minus: Vec<i64>,
    /// X +1 neighbor.
    pub neighbor_x_plus: Vec<i64>,
    /// Y −1 neighbor.
    pub neighbor_y_minus: Vec<i64>,
    /// Y +1 neighbor.
    pub neighbor_y_plus: Vec<i64>,
    /// X lower boundary mask.
    pub boundary_x_lower: Vec<bool>,
    /// X upper boundary mask.
    pub boundary_x_upper: Vec<bool>,
    /// Y lower boundary mask.
    pub boundary_y_lower: Vec<bool>,
    /// Y upper boundary mask.
    pub boundary_y_upper: Vec<bool>,
    /// The stagger this construction was built for.
    pub stagger: Stagger,
    /// Whether the stagger is rotated (E).
    pub rotated: bool,
    /// `(h, u, v)` variable locations for this stagger.
    pub variable_locations: (Location, Location, Location),
}

fn mk(op: &str, args: Vec<Value>) -> Value {
    json!({"op": op, "args": args})
}

struct ArkAst {
    d: Value,      // (hi - lo) / n
    node: Value,   // lo + (i - 1)*d
    center: Value, // lo + (i - 0.5)*d
    prod: Value,   // a * b
}

fn build_ark_ast() -> ArkAst {
    ArkAst {
        d: mk(
            "/",
            vec![mk("-", vec![json!("hi"), json!("lo")]), json!("n")],
        ),
        node: mk(
            "+",
            vec![
                json!("lo"),
                mk("*", vec![mk("-", vec![json!("i"), json!(1)]), json!("d")]),
            ],
        ),
        center: mk(
            "+",
            vec![
                json!("lo"),
                mk("*", vec![mk("-", vec![json!("i"), json!(0.5)]), json!("d")]),
            ],
        ),
        prod: mk("*", vec![json!("a"), json!("b")]),
    }
}

fn spacing(ast: &ArkAst, lo: f64, hi: f64, n: usize) -> f64 {
    let mut b = Bindings::new();
    b.insert("hi".to_string(), hi);
    b.insert("lo".to_string(), lo);
    b.insert("n".to_string(), n as f64);
    eval_coeff(&ast.d, &b).expect("arakawa faq: spacing")
}

/// `n` cell-centre coordinates `center[i] = lo + (i - 0.5)*d` (1-based `i`).
fn center_map(ast: &ArkAst, n: usize, lo: f64, d: f64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let mut b = Bindings::new();
            b.insert("lo".to_string(), lo);
            b.insert("i".to_string(), (i + 1) as f64);
            b.insert("d".to_string(), d);
            eval_coeff(&ast.center, &b).expect("arakawa faq: center")
        })
        .collect()
}

/// `np1` cell-face/corner coordinates `node[i] = lo + (i - 1)*d` (1-based `i`).
fn node_map(ast: &ArkAst, np1: usize, lo: f64, d: f64) -> Vec<f64> {
    (0..np1)
        .map(|i| {
            let mut b = Bindings::new();
            b.insert("lo".to_string(), lo);
            b.insert("i".to_string(), (i + 1) as f64);
            b.insert("d".to_string(), d);
            eval_coeff(&ast.node, &b).expect("arakawa faq: node")
        })
        .collect()
}

/// Cartesian-product gather of two axis maps into a location's row-major
/// `(xs, ys)` arrays. `xmap` has length `ni` (the location's first-dim size),
/// `ymap` length `nj`; `out[j*ni + i] = (xmap[i], ymap[j])`.
fn gather(xmap: &[f64], ymap: &[f64], ni: usize, nj: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0_f64; ni * nj];
    let mut ys = vec![0.0_f64; ni * nj];
    for j in 0..nj {
        for i in 0..ni {
            xs[j * ni + i] = xmap[i];
            ys[j * ni + i] = ymap[j];
        }
    }
    (xs, ys)
}

/// Row-major neighbor over the `(nx, ny)` cell-centre set for one axis
/// (`axis_x` ⇒ x) and `offset` (∓1), `-1` off-grid.
fn axis_neighbor(nx: usize, ny: usize, axis_x: bool, offset: i64) -> Vec<i64> {
    let mut out = vec![0i64; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let ii = if axis_x { i as i64 + offset } else { i as i64 };
            let jj = if axis_x { j as i64 } else { j as i64 + offset };
            out[j * nx + i] = if ii >= 0 && (ii as usize) < nx && jj >= 0 && (jj as usize) < ny {
                (jj as usize * nx + ii as usize) as i64
            } else {
                -1
            };
        }
    }
    out
}

/// Row-major boundary mask over the `(nx, ny)` cell-centre set.
fn axis_boundary(nx: usize, ny: usize, axis_x: bool, lower: bool) -> Vec<bool> {
    let target = if lower {
        0
    } else if axis_x {
        nx - 1
    } else {
        ny - 1
    };
    let mut out = vec![false; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let v = if axis_x { i } else { j };
            out[j * nx + i] = v == target;
        }
    }
    out
}

/// Static A–E variable-location table, reconstructed from the stagger label
/// INDEPENDENT of the imperative `Stagger::variable_locations` so the conformance
/// test is a genuine cross-check. `h` is always CellCenter.
fn variable_locations(s: Stagger) -> (Location, Location, Location) {
    let (u, v) = match s {
        Stagger::A => (Location::CellCenter, Location::CellCenter),
        Stagger::B => (Location::Corner, Location::Corner),
        Stagger::C => (Location::UEdge, Location::VEdge),
        Stagger::D => (Location::VEdge, Location::UEdge),
        Stagger::E => (Location::Corner, Location::Corner),
    };
    (Location::CellCenter, u, v)
}

/// Materialize the arakawa construction (over a cartesian base) from the
/// declarative FAQ, routing all arithmetic through `eval_coeff`.
pub fn arakawa_construction_faq(base: &CartesianBase, stagger: Stagger) -> ArakawaConstruction {
    let ast = build_ark_ast();
    let nx = base.nx();
    let ny = base.ny();
    let dx = spacing(&ast, base.xlo(), base.xhi(), nx);
    let dy = spacing(&ast, base.ylo(), base.yhi(), ny);

    let center_x = center_map(&ast, nx, base.xlo(), dx);
    let center_y = center_map(&ast, ny, base.ylo(), dy);
    let node_x = node_map(&ast, nx + 1, base.xlo(), dx);
    let node_y = node_map(&ast, ny + 1, base.ylo(), dy);

    let coords_cell_center = gather(&center_x, &center_y, nx, ny);
    let coords_u_edge = gather(&node_x, &center_y, nx + 1, ny);
    let coords_v_edge = gather(&center_x, &node_y, nx, ny + 1);
    let coords_corner = gather(&node_x, &node_y, nx + 1, ny + 1);

    let cell_volume = {
        let mut b = Bindings::new();
        b.insert("a".to_string(), dx);
        b.insert("b".to_string(), dy);
        eval_coeff(&ast.prod, &b).expect("arakawa faq: cell_volume")
    };

    ArakawaConstruction {
        nx,
        ny,
        dx,
        dy,
        cell_volume,
        coords_cell_center,
        coords_u_edge,
        coords_v_edge,
        coords_corner,
        neighbor_x_minus: axis_neighbor(nx, ny, true, -1),
        neighbor_x_plus: axis_neighbor(nx, ny, true, 1),
        neighbor_y_minus: axis_neighbor(nx, ny, false, -1),
        neighbor_y_plus: axis_neighbor(nx, ny, false, 1),
        boundary_x_lower: axis_boundary(nx, ny, true, true),
        boundary_x_upper: axis_boundary(nx, ny, true, false),
        boundary_y_lower: axis_boundary(nx, ny, false, true),
        boundary_y_upper: axis_boundary(nx, ny, false, false),
        stagger,
        rotated: stagger == Stagger::E,
        variable_locations: variable_locations(stagger),
    }
}
