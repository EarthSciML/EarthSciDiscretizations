//! Cartesian structured-grid construction via the landed M1 elementwise FAQ
//! path (`eval_coeff`, `crate::rule_eval`).
//!
//! Rust mirror of `src/cartesian_faq.jl` (esd-3we.1 / S1). Declarative
//! companion: `discretizations/grids/cartesian/rules/cartesian_construction.esm`.
//! That document expresses the cartesian construction as a semiring-FAQ
//! (RFC `semiring-faq-unified-ir` §5.1/§5.2: cartesian interval index-sets +
//! elementwise `arrayop` bodies — the M1 path, NO value-invention, all CONST
//! cadence). It is the structured-grid analogue of the DUO value-invention path
//! in `duo_topology_faq.rs`.
//!
//! This module is the thin **evaluation bridge**: it shapes a `CartesianGrid`
//! into the per-axis elementwise inputs and routes every piece of grid
//! ARITHMETIC — the affine coordinate map, the cell-center/width derivations,
//! and the cell-volume product — through ESS's AST pathway (`eval_coeff`). ESD
//! hosts no shadow evaluator (AGENTS.md single-pathway invariant; GRIDS_API
//! §4.3), so the determinism contract — and therefore the cross-binding
//! byte-identity of the output — lives entirely in `earthsci_toolkit`. The only
//! ESD-side logic is the *structural* part — the cartesian-product flattening
//! (column-major), the neighbor linearization with the `-1` off-grid sentinel,
//! the boundary masks, and the identity-metric Jacobian — exactly as
//! `cartesian_faq.jl` materializes those arrays structurally while routing the
//! arithmetic through ESS.
//!
//! Byte-identity: the Rust `eval_coeff` and the imperative builder share the
//! same `f64` libm and association order (the affine map `lo + i*dx` and the
//! width difference fold identically), so the FAQ output is bit-for-bit equal to
//! the imperative `CartesianGrid` arrays and to the committed Julia-reference
//! `tests/conformance/grids/cartesian/construction/golden.json`.

use serde_json::{json, Value};

use super::cartesian::CartesianGrid;
use crate::{eval_coeff, Bindings};

/// Materialized cartesian construction. Float arrays are `f64`; the declared
/// `dtype` tracks precision per GRIDS_API §7. The fields mirror the imperative
/// `CartesianGrid` trait arrays and the `construction/golden.json` contract.
///
/// Public (rather than `pub(crate)` like the DUO front-door bridges) because the
/// structured-grid FAQ is a **parallel construction oracle**: the imperative
/// builder remains the front-door this bead (esd-dru, the py/rust W-step), and
/// the FAQ is exercised by external conformance/parity tests. The S5 reroute
/// (esd-3we.5) wires it into the constructor and tightens visibility.
pub struct CartesianConstruction {
    /// Per-axis edge coordinates; `edges[d]` has length `n[d] + 1`.
    pub edges: Vec<Vec<f64>>,
    /// Per-axis cell centers; `centers[d]` has length `n[d]`.
    pub centers: Vec<Vec<f64>>,
    /// Per-axis cell widths; `widths[d]` has length `n[d]`.
    pub widths: Vec<Vec<f64>>,
    /// Per-axis cell centers flattened column-major over the cell set (one value
    /// per cell, length `n_cells`).
    pub cell_centers: Vec<Vec<f64>>,
    /// Per-axis cell widths flattened column-major (length `n_cells`).
    pub cell_widths: Vec<Vec<f64>>,
    /// Per-cell measure = product of axis widths (length `n_cells`).
    pub cell_volume: Vec<f64>,
    /// Identity-metric Jacobian = `cell_volume`.
    pub metric_jacobian: Vec<f64>,
    /// Per-axis flattened linear id (0-based, column-major) of the offset −1
    /// neighbor, with `-1` for off-grid.
    pub neighbor_minus: Vec<Vec<i64>>,
    /// Per-axis flattened linear id of the offset +1 neighbor, `-1` off-grid.
    pub neighbor_plus: Vec<Vec<i64>>,
    /// Per-axis first-cell-of-axis boundary mask (length `n_cells`).
    pub boundary_lower: Vec<Vec<bool>>,
    /// Per-axis last-cell-of-axis boundary mask (length `n_cells`).
    pub boundary_upper: Vec<Vec<bool>>,
}

// AST nodes mirroring the equation bodies of
// `discretizations/grids/cartesian/rules/cartesian_construction.esm`. The leaf
// `i` is the 1-based interface/cell index (binding-neutral with the Julia
// const), so a uniform edge is `lo + (i-1)*dx`.
fn mk(op: &str, args: Vec<Value>) -> Value {
    json!({"op": op, "args": args})
}

/// The four elementwise construction ASTs, built once per call.
struct CartAst {
    dx: Value,     // (hi - lo) / n
    edge: Value,   // lo + (i - 1)*dx
    center: Value, // (e_lo + e_hi) / 2
    width: Value,  // e_hi - e_lo
    prod: Value,   // a * b
}

fn build_cart_ast() -> CartAst {
    CartAst {
        dx: mk(
            "/",
            vec![mk("-", vec![json!("hi"), json!("lo")]), json!("n")],
        ),
        edge: mk(
            "+",
            vec![
                json!("lo"),
                mk("*", vec![mk("-", vec![json!("i"), json!(1)]), json!("dx")]),
            ],
        ),
        center: mk(
            "/",
            vec![mk("+", vec![json!("e_lo"), json!("e_hi")]), json!(2)],
        ),
        width: mk("-", vec![json!("e_hi"), json!("e_lo")]),
        prod: mk("*", vec![json!("a"), json!("b")]),
    }
}

/// Construct one axis's `(edges, centers, widths)` via the M1 elementwise FAQ.
///
/// A uniform axis is generated by the affine map `edges[i] = lo + i*dx`
/// (0-based `i`) with `dx = (hi - lo)/n`; a non-uniform axis takes
/// `supplied_edges` as DATA. In both cases `centers[i] = (edges[i] +
/// edges[i+1])/2` and `widths[i] = edges[i+1] - edges[i]`. Every arithmetic step
/// is evaluated by ESS so the result is bit-identical to the other bindings.
fn faq_axis(
    ast: &CartAst,
    uniform: bool,
    n: usize,
    lo: f64,
    hi: f64,
    supplied_edges: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut edges = vec![0.0_f64; n + 1];
    if uniform {
        let mut b = Bindings::new();
        b.insert("hi".to_string(), hi);
        b.insert("lo".to_string(), lo);
        b.insert("n".to_string(), n as f64);
        let dx = eval_coeff(&ast.dx, &b).expect("cartesian faq: dx");
        for (i, e) in edges.iter_mut().enumerate() {
            let mut eb = Bindings::new();
            eb.insert("lo".to_string(), lo);
            eb.insert("i".to_string(), (i + 1) as f64); // 1-based binding
            eb.insert("dx".to_string(), dx);
            *e = eval_coeff(&ast.edge, &eb).expect("cartesian faq: edge");
        }
    } else {
        edges.copy_from_slice(&supplied_edges[..n + 1]);
    }
    let mut centers = vec![0.0_f64; n];
    let mut widths = vec![0.0_f64; n];
    for i in 0..n {
        let mut b = Bindings::new();
        b.insert("e_lo".to_string(), edges[i]);
        b.insert("e_hi".to_string(), edges[i + 1]);
        centers[i] = eval_coeff(&ast.center, &b).expect("cartesian faq: center");
        widths[i] = eval_coeff(&ast.width, &b).expect("cartesian faq: width");
    }
    (edges, centers, widths)
}

/// Column-major strides for shape `n` (axis 0 fastest), matching the Julia
/// `LinearIndices` / `CartesianIndices` layout.
fn col_major_strides(n: &[usize]) -> Vec<usize> {
    let mut s = vec![1usize; n.len()];
    for d in 1..n.len() {
        s[d] = s[d - 1] * n[d - 1];
    }
    s
}

/// Column-major multi-index of flat cell `k` (axis 0 fastest).
fn unflatten(k: usize, n: &[usize], strides: &[usize]) -> Vec<usize> {
    (0..n.len()).map(|d| (k / strides[d]) % n[d]).collect()
}

/// Materialize the cartesian grid construction from the declarative FAQ, routing
/// all arithmetic through the landed ESS AST evaluator (`eval_coeff`).
///
/// Byte-identical (`f64`) to the imperative `CartesianGrid` arrays and to the
/// committed Julia-reference construction golden.
pub fn cartesian_construction_faq(grid: &CartesianGrid) -> CartesianConstruction {
    let ast = build_cart_ast();
    let ndim = grid.ndim();
    let n: Vec<usize> = grid.n().to_vec();
    let nc: usize = n.iter().product();
    let strides = col_major_strides(&n);

    // --- per-axis elementwise FAQ (affine coords / centers / widths) ---------
    let mut ax_edges: Vec<Vec<f64>> = Vec::with_capacity(ndim);
    let mut ax_centers: Vec<Vec<f64>> = Vec::with_capacity(ndim);
    let mut ax_widths: Vec<Vec<f64>> = Vec::with_capacity(ndim);
    for (d, &nd) in n.iter().enumerate() {
        let (lo, hi) = grid.extent()[d];
        let (e, c, w) = faq_axis(&ast, grid.uniform()[d], nd, lo, hi, &grid.edges()[d]);
        ax_edges.push(e);
        ax_centers.push(c);
        ax_widths.push(w);
    }

    // --- per-cell flattened coordinate gathers (column-major structure) ------
    let cell_centers: Vec<Vec<f64>> = (0..ndim)
        .map(|d| {
            (0..nc)
                .map(|k| ax_centers[d][(k / strides[d]) % n[d]])
                .collect()
        })
        .collect();
    let cell_widths: Vec<Vec<f64>> = (0..ndim)
        .map(|d| {
            (0..nc)
                .map(|k| ax_widths[d][(k / strides[d]) % n[d]])
                .collect()
        })
        .collect();

    // --- cell_volume = Π_d widths_d, the product routed through ESS ----------
    let cell_volume: Vec<f64> = (0..nc)
        .map(|k| {
            let mut acc = cell_widths[0][k];
            for cw_d in cell_widths.iter().take(ndim).skip(1) {
                let mut b = Bindings::new();
                b.insert("a".to_string(), acc);
                b.insert("b".to_string(), cw_d[k]);
                acc = eval_coeff(&ast.prod, &b).expect("cartesian faq: cell_volume");
            }
            acc
        })
        .collect();
    let metric_jacobian = cell_volume.clone();

    // --- neighbor maps (0-based linear ids, -1 sentinel) --------------------
    let neighbors = |offset: i64| -> Vec<Vec<i64>> {
        (0..ndim)
            .map(|d| {
                (0..nc)
                    .map(|k| {
                        let ix = unflatten(k, &n, &strides);
                        let new_d = ix[d] as i64 + offset;
                        if new_d >= 0 && (new_d as usize) < n[d] {
                            let mut lin = 0i64;
                            for (p, &ixp) in ix.iter().enumerate().take(ndim) {
                                let coord = if p == d { new_d as usize } else { ixp };
                                lin += (coord * strides[p]) as i64;
                            }
                            lin
                        } else {
                            -1
                        }
                    })
                    .collect()
            })
            .collect()
    };
    let neighbor_minus = neighbors(-1);
    let neighbor_plus = neighbors(1);

    // --- boundary masks -----------------------------------------------------
    let boundary = |lower: bool| -> Vec<Vec<bool>> {
        (0..ndim)
            .map(|d| {
                (0..nc)
                    .map(|k| {
                        let ix = (k / strides[d]) % n[d];
                        if lower {
                            ix == 0
                        } else {
                            ix == n[d] - 1
                        }
                    })
                    .collect()
            })
            .collect()
    };
    let boundary_lower = boundary(true);
    let boundary_upper = boundary(false);

    CartesianConstruction {
        edges: ax_edges,
        centers: ax_centers,
        widths: ax_widths,
        cell_centers,
        cell_widths,
        cell_volume,
        metric_jacobian,
        neighbor_minus,
        neighbor_plus,
        boundary_lower,
        boundary_upper,
    }
}
