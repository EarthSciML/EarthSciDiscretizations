//! Lat-lon structured-grid construction via the landed M1 elementwise FAQ path
//! (`eval_coeff`, `crate::rule_eval`).
//!
//! Rust mirror of `src/latlon_faq.jl` (esd-3we, S3). Declarative companion:
//! `discretizations/grids/latlon/rules/latlon_construction.esm`. The construction
//! is a semiring-FAQ (RFC `semiring-faq-unified-ir` §5.1/§5.2: interval
//! index-sets + elementwise `arrayop` bodies — the M1 path, NO value-invention).
//!
//! This module is the thin **evaluation bridge**: every piece of grid ARITHMETIC
//! — the per-row affine longitude map (`dlon = 2π/nlon`, `lon = lon_start +
//! (i-½)·dlon`), the spherical-rectangle area `R²·dlon·(sin φ_n − sin φ_s)`, and
//! the closed-form curvilinear metric (`g_λλ = R²cos²φ`, `g_φφ = R²`, `ginv`,
//! Jacobian `R²|cos φ|`, and the lone non-vanishing derivative `∂g_λλ/∂φ =
//! −2R²cos φ sin φ`) — flows through ESS's single AST pathway (`eval_coeff`,
//! which dispatches `sin`/`cos`/`abs`). The structural integer arrays (the
//! ragged-row-major flattening, periodic-longitude neighbor linearization, the
//! reduced-Gaussian nearest-center N/S remap, boundary masks) are pure index
//! logic. Byte/ULP-identical to the imperative `LatLonGrid` accessors and to the
//! committed Julia-reference
//! `tests/conformance/grids/latlon/construction/golden.json`.

use serde_json::{json, Value};

use super::lat_lon::LatLonGrid;
use crate::{eval_coeff, Bindings};

/// Materialized lat-lon construction. All per-cell arrays use the flat
/// ragged-row-major layout (cell `(j, i)` at flat index `row_off[j] + i`).
/// The metric is stored as the golden's named scalar components.
pub struct LatLonConstruction {
    /// Per-cell longitude (rad).
    pub cell_centers_lon: Vec<f64>,
    /// Per-cell latitude (rad) — the row-center value.
    pub cell_centers_lat: Vec<f64>,
    /// Per-cell longitude width `dlon`.
    pub cell_widths_lon: Vec<f64>,
    /// Per-cell latitude width (lat-edge difference).
    pub cell_widths_lat: Vec<f64>,
    /// Per-cell spherical-rectangle area.
    pub cell_volume: Vec<f64>,
    /// 1-D latitude edges (length `nlat + 1`).
    pub lat_edges: Vec<f64>,
    /// 1-D latitude centers (length `nlat`).
    pub lat_centers: Vec<f64>,
    /// `g_λλ = R²cos²φ` per cell.
    pub metric_g_lonlon: Vec<f64>,
    /// `g_φφ = R²` per cell.
    pub metric_g_latlat: Vec<f64>,
    /// `ginv_λλ = 1/g_λλ` per cell.
    pub metric_ginv_lonlon: Vec<f64>,
    /// Jacobian `R²|cos φ|` per cell.
    pub metric_jacobian: Vec<f64>,
    /// `∂g_λλ/∂φ = −2R²cos φ sin φ` per cell.
    pub dg_lonlon_dlat: Vec<f64>,
    /// Periodic longitude −1 neighbor (flat linear id).
    pub neighbor_lon_minus: Vec<i64>,
    /// Periodic longitude +1 neighbor (flat linear id).
    pub neighbor_lon_plus: Vec<i64>,
    /// Latitude −1 neighbor (nearest-center remap; `-1` at the pole).
    pub neighbor_lat_minus: Vec<i64>,
    /// Latitude +1 neighbor (`-1` at the pole).
    pub neighbor_lat_plus: Vec<i64>,
    /// First-latitude-row boundary mask.
    pub boundary_lat_lower: Vec<bool>,
    /// Last-latitude-row boundary mask.
    pub boundary_lat_upper: Vec<bool>,
}

fn mk(op: &str, args: Vec<Value>) -> Value {
    json!({"op": op, "args": args})
}

/// The elementwise construction ASTs, built once per call. `2π` is
/// `std::f64::consts::TAU`, the nearest f64 to 2π — bit-identical to the Julia
/// const `6.283185307179586` and the imperative `2 * pi` — so `dlon = 2π/n` is
/// bit-for-bit equal across bindings.
struct LlAst {
    dlon: Value,       // 2π / n
    lon_center: Value, // lon_start + (i - 0.5)*dlon
    lat_width: Value,  // e_hi - e_lo
    area: Value,       // ((R*R)*dlon) * (sin(lat_n) - sin(lat_s))
    cos: Value,        // cos(lat)
    gll: Value,        // ((R*R)*cos_lat) * cos_lat
    gpp: Value,        // R*R
    inv: Value,        // 1 / x
    jac: Value,        // (R*R) * abs(cos(lat))
    dgll: Value,       // ((((-2)*R)*R)*cos(phi)) * sin(phi)
}

fn build_ll_ast() -> LlAst {
    LlAst {
        dlon: mk("/", vec![json!(std::f64::consts::TAU), json!("n")]),
        lon_center: mk(
            "+",
            vec![
                json!("lon_start"),
                mk(
                    "*",
                    vec![mk("-", vec![json!("i"), json!(0.5)]), json!("dlon")],
                ),
            ],
        ),
        lat_width: mk("-", vec![json!("e_hi"), json!("e_lo")]),
        area: mk(
            "*",
            vec![
                mk(
                    "*",
                    vec![mk("*", vec![json!("R"), json!("R")]), json!("dlon")],
                ),
                mk(
                    "-",
                    vec![
                        mk("sin", vec![json!("lat_n")]),
                        mk("sin", vec![json!("lat_s")]),
                    ],
                ),
            ],
        ),
        cos: mk("cos", vec![json!("lat")]),
        gll: mk(
            "*",
            vec![
                mk(
                    "*",
                    vec![mk("*", vec![json!("R"), json!("R")]), json!("cos_lat")],
                ),
                json!("cos_lat"),
            ],
        ),
        gpp: mk("*", vec![json!("R"), json!("R")]),
        inv: mk("/", vec![json!(1.0), json!("x")]),
        jac: mk(
            "*",
            vec![
                mk("*", vec![json!("R"), json!("R")]),
                mk("abs", vec![mk("cos", vec![json!("lat")])]),
            ],
        ),
        dgll: mk(
            "*",
            vec![
                mk(
                    "*",
                    vec![
                        mk(
                            "*",
                            vec![mk("*", vec![json!(-2.0), json!("R")]), json!("R")],
                        ),
                        mk("cos", vec![json!("phi")]),
                    ],
                ),
                mk("sin", vec![json!("phi")]),
            ],
        ),
    }
}

/// Reduced-Gaussian nearest-center N/S column remap (0-based), mirroring the
/// imperative `_latlon_map_i` (the relational rank/join materialized
/// structurally).
fn latlon_map_i(i: usize, from_n: usize, to_n: usize) -> usize {
    if from_n == to_n {
        return i;
    }
    let frac = (i as f64 + 0.5) / from_n as f64;
    let k = (frac * to_n as f64).floor() as i64;
    k.clamp(0, to_n as i64 - 1) as usize
}

/// Materialize the lat-lon grid construction from the declarative FAQ, routing
/// all arithmetic through the landed ESS AST evaluator (`eval_coeff`).
pub fn lat_lon_construction_faq(grid: &LatLonGrid) -> LatLonConstruction {
    let ast = build_ll_ast();
    let nlat = grid.nlat();
    let per_row = grid.nlon_per_row();
    let nc = grid.n_cells();
    let r = grid.r();
    let lon_start = grid.lon_start();
    let lat_edges = grid.lat_edges();
    let lat_centers = grid.lat_centers();

    // Row offsets (0-based flat start of each row).
    let mut row_off = vec![0usize; nlat + 1];
    for j in 0..nlat {
        row_off[j + 1] = row_off[j] + per_row[j];
    }

    let mut cc_lon = vec![0.0_f64; nc];
    let mut cc_lat = vec![0.0_f64; nc];
    let mut cw_lon = vec![0.0_f64; nc];
    let mut cw_lat = vec![0.0_f64; nc];
    let mut cell_volume = vec![0.0_f64; nc];
    let mut g_lonlon = vec![0.0_f64; nc];
    let mut g_latlat = vec![0.0_f64; nc];
    let mut ginv_lonlon = vec![0.0_f64; nc];
    let mut jacobian = vec![0.0_f64; nc];
    let mut dg = vec![0.0_f64; nc];

    // --- per-row elementwise FAQ (affine lon, trig metric, area) ------------
    for j in 0..nlat {
        let n_i = per_row[j];
        let base = row_off[j];
        let dlon = {
            let mut b = Bindings::new();
            b.insert("n".to_string(), n_i as f64);
            eval_coeff(&ast.dlon, &b).expect("latlon faq: dlon")
        };
        let lat_c = lat_centers[j];
        let lat_s = lat_edges[j];
        let lat_n = lat_edges[j + 1];
        let lat_w = {
            let mut b = Bindings::new();
            b.insert("e_hi".to_string(), lat_n);
            b.insert("e_lo".to_string(), lat_s);
            eval_coeff(&ast.lat_width, &b).expect("latlon faq: lat width")
        };
        let mut rb = Bindings::new();
        rb.insert("R".to_string(), r);
        rb.insert("dlon".to_string(), dlon);
        rb.insert("lat_n".to_string(), lat_n);
        rb.insert("lat_s".to_string(), lat_s);
        let area = eval_coeff(&ast.area, &rb).expect("latlon faq: area");
        let cos_lat = {
            let mut b = Bindings::new();
            b.insert("lat".to_string(), lat_c);
            eval_coeff(&ast.cos, &b).expect("latlon faq: cos")
        };
        let g_ll = {
            let mut b = Bindings::new();
            b.insert("R".to_string(), r);
            b.insert("cos_lat".to_string(), cos_lat);
            eval_coeff(&ast.gll, &b).expect("latlon faq: g_ll")
        };
        let g_pp = {
            let mut b = Bindings::new();
            b.insert("R".to_string(), r);
            eval_coeff(&ast.gpp, &b).expect("latlon faq: g_pp")
        };
        let inv_ll = if g_ll > 0.0 {
            let mut b = Bindings::new();
            b.insert("x".to_string(), g_ll);
            eval_coeff(&ast.inv, &b).expect("latlon faq: inv_ll")
        } else {
            f64::INFINITY
        };
        let jac = {
            let mut b = Bindings::new();
            b.insert("R".to_string(), r);
            b.insert("lat".to_string(), lat_c);
            eval_coeff(&ast.jac, &b).expect("latlon faq: jacobian")
        };
        let dgll = {
            let mut b = Bindings::new();
            b.insert("R".to_string(), r);
            b.insert("phi".to_string(), lat_c);
            eval_coeff(&ast.dgll, &b).expect("latlon faq: dg")
        };
        for i in 0..n_i {
            let k = base + i;
            let lon = {
                let mut b = Bindings::new();
                b.insert("lon_start".to_string(), lon_start);
                b.insert("i".to_string(), (i + 1) as f64); // 1-based binding
                b.insert("dlon".to_string(), dlon);
                eval_coeff(&ast.lon_center, &b).expect("latlon faq: lon")
            };
            cc_lon[k] = lon;
            cc_lat[k] = lat_c;
            cw_lon[k] = dlon;
            cw_lat[k] = lat_w;
            cell_volume[k] = area;
            g_lonlon[k] = g_ll;
            g_latlat[k] = g_pp;
            ginv_lonlon[k] = inv_ll;
            jacobian[k] = jac;
            dg[k] = dgll;
        }
    }

    // --- neighbor maps (flat ragged-row-major linear ids) -------------------
    let mut nbr_lon_minus = vec![0i64; nc];
    let mut nbr_lon_plus = vec![0i64; nc];
    for j in 0..nlat {
        let n_i = per_row[j];
        let base = row_off[j];
        for i in 0..n_i {
            let im = (i + n_i - 1) % n_i; // periodic longitude wrap
            let ip = (i + 1) % n_i;
            nbr_lon_minus[base + i] = (base + im) as i64;
            nbr_lon_plus[base + i] = (base + ip) as i64;
        }
    }

    let nbr_lat = |offset: i64| -> Vec<i64> {
        let mut out = vec![0i64; nc];
        for j in 0..nlat {
            let n_i = per_row[j];
            let base = row_off[j];
            let jj = j as i64 + offset;
            if jj >= 0 && (jj as usize) < nlat {
                let jt = jj as usize;
                let n_t = per_row[jt];
                let base_t = row_off[jt];
                for i in 0..n_i {
                    out[base + i] = (base_t + latlon_map_i(i, n_i, n_t)) as i64;
                }
            } else {
                for i in 0..n_i {
                    out[base + i] = -1; // pole: pole_policy = none
                }
            }
        }
        out
    };
    let neighbor_lat_minus = nbr_lat(-1);
    let neighbor_lat_plus = nbr_lat(1);

    // --- boundary masks (lon wraps → all-false; lat marks first/last row) ---
    let bd_lat = |target: usize| -> Vec<bool> {
        let mut out = vec![false; nc];
        for j in 0..nlat {
            let n_i = per_row[j];
            let base = row_off[j];
            let mark = j == target;
            for i in 0..n_i {
                out[base + i] = mark;
            }
        }
        out
    };
    let boundary_lat_lower = bd_lat(0);
    let boundary_lat_upper = bd_lat(nlat - 1);

    LatLonConstruction {
        cell_centers_lon: cc_lon,
        cell_centers_lat: cc_lat,
        cell_widths_lon: cw_lon,
        cell_widths_lat: cw_lat,
        cell_volume,
        lat_edges: lat_edges.to_vec(),
        lat_centers: lat_centers.to_vec(),
        metric_g_lonlon: g_lonlon,
        metric_g_latlat: g_latlat,
        metric_ginv_lonlon: ginv_lonlon,
        metric_jacobian: jacobian,
        dg_lonlon_dlat: dg,
        neighbor_lon_minus: nbr_lon_minus,
        neighbor_lon_plus: nbr_lon_plus,
        neighbor_lat_minus,
        neighbor_lat_plus,
        boundary_lat_lower,
        boundary_lat_upper,
    }
}
