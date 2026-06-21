//! Vertical (1-D column) structured-grid construction via the landed M1
//! elementwise FAQ path (`eval_coeff`, `crate::rule_eval`).
//!
//! Rust mirror of `src/vertical_faq.jl` (esd-3we, S2). Declarative companion:
//! `discretizations/grids/vertical/rules/vertical_construction.esm`. The
//! construction is a semiring-FAQ (RFC `semiring-faq-unified-ir` §5.1/§5.2: two
//! interval index-sets — `interfaces` (nz+1 half-levels) and `layers` (nz cells)
//! — plus elementwise `arrayop` bodies; the M1 path, NO value-invention).
//!
//! This module is the thin **evaluation bridge**: every piece of vertical grid
//! ARITHMETIC — the per-coordinate level synthesis, the cell-center/width
//! derivations, the cell-volume, and the named layer metrics — flows through
//! ESS's single AST pathway (`eval_coeff`). The structural integer arrays
//! (±k neighbor maps with the `-1` off-column sentinel, boundary masks) are pure
//! index logic. Byte-identical to the imperative `VerticalGrid` accessors and to
//! the committed Julia-reference
//! `tests/conformance/grids/vertical/construction/golden.json` (strict 0.0
//! tolerance — no transcendentals).
//!
//! Level synthesis mirrors the imperative `vertical::build` branch-for-branch:
//! `eta` is the divide-add `ak[k]/p0 + bk[k]`; uniform `sigma` /
//! `hybrid_sigma_theta` is the affine `1 - k/nz` (used only when it reproduces
//! the stored levels bit-for-bit; an explicit-level grid carries them as DATA);
//! `z` / `theta` / `z_star` supply their levels verbatim.

use serde_json::{json, Value};

use super::vertical::{VerticalCoordinate, VerticalGrid};
use crate::{eval_coeff, Bindings};

/// Materialized vertical construction. Optional metric fields are `None` when
/// the coordinate kind does not define them (matching the imperative
/// `metric_eval` availability). Mirrors the `construction/golden.json` contract.
pub struct VerticalConstruction {
    /// Interface (half-level) coordinates, length `nz + 1`.
    pub levels: Vec<f64>,
    /// Cell centers, length `nz`.
    pub centers: Vec<f64>,
    /// Cell thicknesses, length `nz`.
    pub widths: Vec<f64>,
    /// Per-layer measure (= `widths`).
    pub cell_volume: Vec<f64>,
    /// `metric_eval("dz")` over the column (= `widths`).
    pub metric_dz: Vec<f64>,
    /// `metric_eval("z")` over the column (= `centers`).
    pub metric_z: Vec<f64>,
    /// `metric_eval("sigma")` — `Some` for sigma-like coordinates.
    pub metric_sigma: Option<Vec<f64>>,
    /// `metric_eval("pressure")` — `Some` when hybrid ak/bk are present.
    pub metric_pressure: Option<Vec<f64>>,
    /// `metric_eval("ak")` — `Some` when ak is present.
    pub metric_ak: Option<Vec<f64>>,
    /// `metric_eval("bk")` — `Some` when bk is present.
    pub metric_bk: Option<Vec<f64>>,
    /// Offset −1 layer neighbor (0-based, `-1` off-column).
    pub neighbor_minus: Vec<i64>,
    /// Offset +1 layer neighbor (0-based, `-1` off-column).
    pub neighbor_plus: Vec<i64>,
    /// Bottom-layer boundary mask.
    pub boundary_lower: Vec<bool>,
    /// Top-layer boundary mask.
    pub boundary_upper: Vec<bool>,
    /// Hybrid `ak` coefficients (empty when the kind has none).
    pub ak: Vec<f64>,
    /// Hybrid `bk` coefficients (empty when the kind has none).
    pub bk: Vec<f64>,
    /// Reference surface pressure.
    pub p0: f64,
}

fn mk(op: &str, args: Vec<Value>) -> Value {
    json!({"op": op, "args": args})
}

/// The five elementwise construction ASTs, built once per call.
struct VertAst {
    sigma: Value,    // 1 - k/nz   (k = 0-based interface index)
    eta: Value,      // ak/p0 + bk
    avg2: Value,     // (a + b) / 2
    width: Value,    // abs(b - a)
    pressure: Value, // ((ak_lo + bk_lo*p0) + (ak_hi + bk_hi*p0)) / 2
}

fn build_vert_ast() -> VertAst {
    VertAst {
        sigma: mk("-", vec![json!(1), mk("/", vec![json!("k"), json!("nz")])]),
        eta: mk(
            "+",
            vec![mk("/", vec![json!("ak"), json!("p0")]), json!("bk")],
        ),
        avg2: mk("/", vec![mk("+", vec![json!("a"), json!("b")]), json!(2)]),
        width: mk("abs", vec![mk("-", vec![json!("b"), json!("a")])]),
        pressure: mk(
            "/",
            vec![
                mk(
                    "+",
                    vec![
                        mk(
                            "+",
                            vec![json!("ak_lo"), mk("*", vec![json!("bk_lo"), json!("p0")])],
                        ),
                        mk(
                            "+",
                            vec![json!("ak_hi"), mk("*", vec![json!("bk_hi"), json!("p0")])],
                        ),
                    ],
                ),
                json!(2),
            ],
        ),
    }
}

/// Synthesize the interface (half-level) coordinates via the M1 FAQ, branching
/// on the coordinate kind exactly as the imperative builder does.
fn faq_levels(ast: &VertAst, grid: &VerticalGrid) -> Vec<f64> {
    let nz = grid.n_cells();
    match grid.coordinate() {
        VerticalCoordinate::Eta => {
            let p0 = grid.p0();
            (0..=nz)
                .map(|k| {
                    let mut b = Bindings::new();
                    b.insert("ak".to_string(), grid.ak()[k]);
                    b.insert("p0".to_string(), p0);
                    b.insert("bk".to_string(), grid.bk()[k]);
                    eval_coeff(&ast.eta, &b).expect("vertical faq: eta level")
                })
                .collect()
        }
        VerticalCoordinate::Sigma | VerticalCoordinate::HybridSigmaTheta => {
            let cand: Vec<f64> = (0..=nz)
                .map(|k| {
                    let mut b = Bindings::new();
                    b.insert("k".to_string(), k as f64);
                    b.insert("nz".to_string(), nz as f64);
                    eval_coeff(&ast.sigma, &b).expect("vertical faq: sigma level")
                })
                .collect();
            if cand.as_slice() == grid.levels() {
                cand
            } else {
                grid.levels().to_vec()
            }
        }
        _ => grid.levels().to_vec(),
    }
}

/// Elementwise midpoint / absolute-difference over the `layers` interval.
fn faq_centers_widths(ast: &VertAst, levels: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = levels.len() - 1;
    let mut centers = vec![0.0_f64; n];
    let mut widths = vec![0.0_f64; n];
    for k in 0..n {
        let mut b = Bindings::new();
        b.insert("a".to_string(), levels[k]);
        b.insert("b".to_string(), levels[k + 1]);
        centers[k] = eval_coeff(&ast.avg2, &b).expect("vertical faq: center");
        widths[k] = eval_coeff(&ast.width, &b).expect("vertical faq: width");
    }
    (centers, widths)
}

/// Per-layer two-point average of an interface-defined coefficient.
fn faq_layer_avg(ast: &VertAst, coeff: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let mut b = Bindings::new();
            b.insert("a".to_string(), coeff[k]);
            b.insert("b".to_string(), coeff[k + 1]);
            eval_coeff(&ast.avg2, &b).expect("vertical faq: layer avg")
        })
        .collect()
}

/// Materialize the vertical column construction from the declarative FAQ.
pub fn vertical_construction_faq(grid: &VerticalGrid) -> VerticalConstruction {
    let ast = build_vert_ast();
    let levels = faq_levels(&ast, grid);
    let (centers, widths) = faq_centers_widths(&ast, &levels);
    let n = widths.len();

    let cell_volume = widths.clone();
    let metric_dz = widths.clone();
    let metric_z = centers.clone();

    let sigma_like = matches!(
        grid.coordinate(),
        VerticalCoordinate::Sigma | VerticalCoordinate::HybridSigmaTheta | VerticalCoordinate::Eta
    );
    let metric_sigma = if sigma_like {
        Some(centers.clone())
    } else {
        None
    };

    let has_ak = !grid.ak().is_empty();
    let has_bk = !grid.bk().is_empty();
    let p0 = grid.p0();

    let metric_pressure = if has_ak && has_bk {
        Some(
            (0..n)
                .map(|k| {
                    let mut b = Bindings::new();
                    b.insert("ak_lo".to_string(), grid.ak()[k]);
                    b.insert("bk_lo".to_string(), grid.bk()[k]);
                    b.insert("ak_hi".to_string(), grid.ak()[k + 1]);
                    b.insert("bk_hi".to_string(), grid.bk()[k + 1]);
                    b.insert("p0".to_string(), p0);
                    eval_coeff(&ast.pressure, &b).expect("vertical faq: pressure")
                })
                .collect(),
        )
    } else {
        None
    };
    let metric_ak = if has_ak {
        Some(faq_layer_avg(&ast, grid.ak(), n))
    } else {
        None
    };
    let metric_bk = if has_bk {
        Some(faq_layer_avg(&ast, grid.bk(), n))
    } else {
        None
    };

    let neighbor_minus: Vec<i64> = (0..n)
        .map(|k| if k > 0 { (k - 1) as i64 } else { -1 })
        .collect();
    let neighbor_plus: Vec<i64> = (0..n)
        .map(|k| if k + 1 < n { (k + 1) as i64 } else { -1 })
        .collect();
    let boundary_lower: Vec<bool> = (0..n).map(|k| k == 0).collect();
    let boundary_upper: Vec<bool> = (0..n).map(|k| k == n - 1).collect();

    VerticalConstruction {
        levels,
        centers,
        widths,
        cell_volume,
        metric_dz,
        metric_z,
        metric_sigma,
        metric_pressure,
        metric_ak,
        metric_bk,
        neighbor_minus,
        neighbor_plus,
        boundary_lower,
        boundary_upper,
        ak: grid.ak().to_vec(),
        bk: grid.bk().to_vec(),
        p0,
    }
}
