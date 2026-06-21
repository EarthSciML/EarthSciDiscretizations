//! Cross-language conformance for the cartesian **construction FAQ** (esd-dru).
//!
//! Builds each fixture grid imperatively, re-derives its construction arrays
//! through the elementwise FAQ bridge (`cartesian_construction_faq`, routed via
//! the ESS `eval_coeff` evaluator), and asserts the result is
//!   1. byte/ULP-identical to the imperative `CartesianGrid` accessors (parity),
//!   2. equal to the committed Julia-reference
//!      `tests/conformance/grids/cartesian/construction/golden.json` (the
//!      cross-binding contract — floats within the fixture tolerance, the
//!      integer neighbor maps / boundary masks exactly).
//!
//! Golden arrays are stored as compact JSON strings (the dense byte form shared
//! with `tests/conformance/grids/duo/topology/golden.json`); neighbor ids are
//! 0-based with a `-1` off-grid sentinel and boundary masks are `0`/`1`.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::cartesian::{self, CartesianGrid};
use earthsci_grids::grids::cartesian_faq::cartesian_construction_faq;
use earthsci_grids::Dtype;
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/cartesian")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

fn dtype_from(s: &str) -> Dtype {
    match s {
        "float64" => Dtype::F64,
        "float32" => Dtype::F32,
        other => panic!("unsupported dtype in fixture: {other}"),
    }
}

fn build_grid(opts: &Value) -> CartesianGrid {
    let mut b = cartesian::builder();
    if let Some(dtype) = opts.get("dtype").and_then(|v| v.as_str()) {
        b = b.dtype(dtype_from(dtype));
    }
    if let Some(g) = opts.get("ghosts").and_then(|v| v.as_u64()) {
        b = b.ghosts(g as usize);
    }
    if let Some(edges) = opts.get("edges").and_then(|v| v.as_array()) {
        let edges: Vec<Vec<f64>> = edges
            .iter()
            .map(|axis| {
                axis.as_array()
                    .expect("edges axis is an array")
                    .iter()
                    .map(|x| x.as_f64().expect("edge is a number"))
                    .collect()
            })
            .collect();
        b = b.edges(edges);
    } else {
        if let Some(nx) = opts.get("nx").and_then(|v| v.as_u64()) {
            b = b.nx(nx as usize);
        }
        if let Some(ny) = opts.get("ny").and_then(|v| v.as_u64()) {
            b = b.ny(ny as usize);
        }
        if let Some(nz) = opts.get("nz").and_then(|v| v.as_u64()) {
            b = b.nz(nz as usize);
        }
        if let Some(extent) = opts.get("extent").and_then(|v| v.as_array()) {
            let extent: Vec<(f64, f64)> = extent
                .iter()
                .map(|pair| {
                    let arr = pair.as_array().expect("extent axis is array");
                    (
                        arr[0].as_f64().expect("lo is a number"),
                        arr[1].as_f64().expect("hi is a number"),
                    )
                })
                .collect();
            b = b.extent(extent);
        }
    }
    b.build().expect("cartesian fixture builds")
}

/// Parse a golden compact-JSON-string into a nested `f64` array.
fn parse_axis_floats(s: &str) -> Vec<Vec<f64>> {
    let v: Value = serde_json::from_str(s).expect("golden float-axis string parses");
    v.as_array()
        .expect("axis array")
        .iter()
        .map(|axis| {
            axis.as_array()
                .expect("axis is array")
                .iter()
                .map(|x| x.as_f64().expect("float"))
                .collect()
        })
        .collect()
}

fn parse_floats(s: &str) -> Vec<f64> {
    let v: Value = serde_json::from_str(s).expect("golden flat-float string parses");
    v.as_array()
        .expect("array")
        .iter()
        .map(|x| x.as_f64().expect("float"))
        .collect()
}

fn parse_axis_ints(s: &str) -> Vec<Vec<i64>> {
    let v: Value = serde_json::from_str(s).expect("golden int-axis string parses");
    v.as_array()
        .expect("axis array")
        .iter()
        .map(|axis| {
            axis.as_array()
                .expect("axis is array")
                .iter()
                .map(|x| x.as_i64().expect("int"))
                .collect()
        })
        .collect()
}

/// Column-major strides for shape `n` (axis 0 fastest).
fn col_major_strides(n: &[usize]) -> Vec<usize> {
    let mut s = vec![1usize; n.len()];
    for d in 1..n.len() {
        s[d] = s[d - 1] * n[d - 1];
    }
    s
}

fn unflatten(k: usize, n: &[usize], strides: &[usize]) -> Vec<usize> {
    (0..n.len()).map(|d| (k / strides[d]) % n[d]).collect()
}

#[test]
fn cartesian_construction_faq_matches_golden_and_imperative() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let golden_doc = read_json(&hdir.join("construction/golden.json"));
    let rel_tol = golden_doc["tolerance"]["relative"].as_f64().unwrap();

    let by_name: std::collections::HashMap<String, &Value> = golden_doc["fixtures"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| (f["name"].as_str().unwrap().to_string(), f))
        .collect();

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let grid = build_grid(&fixture["opts"]);
        let ndim = grid.ndim();
        let n: Vec<usize> = grid.n().to_vec();
        let nc: usize = n.iter().product();
        let strides = col_major_strides(&n);
        let faq = cartesian_construction_faq(&grid);
        let gl = by_name
            .get(name)
            .unwrap_or_else(|| panic!("golden missing fixture {name}"));

        // --- counts ---
        assert_eq!(ndim, gl["ndim"].as_u64().unwrap() as usize, "{name}: ndim");
        assert_eq!(
            nc as u64,
            gl["n_cells"].as_u64().unwrap(),
            "{name}: n_cells"
        );

        // --- parity: FAQ == imperative accessors, bit-for-bit (same f64 libm) ---
        for d in 0..ndim {
            assert_eq!(faq.edges[d], grid.edges()[d], "{name}: edges[{d}] parity");
            assert_eq!(
                faq.centers[d],
                grid.centers()[d],
                "{name}: centers[{d}] parity"
            );
            assert_eq!(
                faq.widths[d],
                grid.widths()[d],
                "{name}: widths[{d}] parity"
            );
        }
        for k in 0..nc {
            let idx = unflatten(k, &n, &strides);
            let cc = grid.cell_center(&idx).unwrap();
            let cw = grid.cell_widths(&idx).unwrap();
            for d in 0..ndim {
                assert_eq!(faq.cell_centers[d][k], cc[d], "{name}: cell_centers parity");
                assert_eq!(faq.cell_widths[d][k], cw[d], "{name}: cell_widths parity");
            }
            assert_eq!(
                faq.cell_volume[k],
                grid.cell_volume(&idx).unwrap(),
                "{name}: cell_volume parity"
            );
        }

        // --- byte/ULP identity against the cross-binding golden ---
        let g_edges = parse_axis_floats(gl["edges"].as_str().unwrap());
        let g_centers = parse_axis_floats(gl["centers"].as_str().unwrap());
        let g_widths = parse_axis_floats(gl["widths"].as_str().unwrap());
        for d in 0..ndim {
            for (a, b) in faq.edges[d].iter().zip(g_edges[d].iter()) {
                assert!(close_rel(*a, *b, rel_tol), "{name}: edges[{d}] golden");
            }
            for (a, b) in faq.centers[d].iter().zip(g_centers[d].iter()) {
                assert!(close_rel(*a, *b, rel_tol), "{name}: centers[{d}] golden");
            }
            for (a, b) in faq.widths[d].iter().zip(g_widths[d].iter()) {
                assert!(close_rel(*a, *b, rel_tol), "{name}: widths[{d}] golden");
            }
        }

        let g_vol = parse_floats(gl["cell_volume"].as_str().unwrap());
        let g_jac = parse_floats(gl["metric_jacobian"].as_str().unwrap());
        for (a, b) in faq.cell_volume.iter().zip(g_vol.iter()) {
            assert!(close_rel(*a, *b, rel_tol), "{name}: cell_volume golden");
        }
        for (a, b) in faq.metric_jacobian.iter().zip(g_jac.iter()) {
            assert!(close_rel(*a, *b, rel_tol), "{name}: metric_jacobian golden");
        }

        let g_nm = parse_axis_ints(gl["neighbor_minus"].as_str().unwrap());
        let g_np = parse_axis_ints(gl["neighbor_plus"].as_str().unwrap());
        let g_bl = parse_axis_ints(gl["boundary_lower"].as_str().unwrap());
        let g_bu = parse_axis_ints(gl["boundary_upper"].as_str().unwrap());
        for d in 0..ndim {
            assert_eq!(
                faq.neighbor_minus[d], g_nm[d],
                "{name}: neighbor_minus[{d}]"
            );
            assert_eq!(faq.neighbor_plus[d], g_np[d], "{name}: neighbor_plus[{d}]");
            let bl: Vec<i64> = faq.boundary_lower[d].iter().map(|&x| x as i64).collect();
            let bu: Vec<i64> = faq.boundary_upper[d].iter().map(|&x| x as i64).collect();
            assert_eq!(bl, g_bl[d], "{name}: boundary_lower[{d}]");
            assert_eq!(bu, g_bu[d], "{name}: boundary_upper[{d}]");
        }
    }
}
