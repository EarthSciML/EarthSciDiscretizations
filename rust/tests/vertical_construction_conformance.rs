//! Cross-language conformance for the vertical **construction FAQ** (esd-dru).
//!
//! Builds each fixture grid imperatively, re-derives its construction arrays
//! through the elementwise FAQ bridge (`vertical_construction_faq`, routed via
//! the ESS `eval_coeff` evaluator), and asserts the result is
//!   1. byte-identical to the imperative `VerticalGrid` accessors (parity),
//!   2. equal to the committed Julia-reference
//!      `tests/conformance/grids/vertical/construction/golden.json`.
//!
//! The vertical golden carries `tolerance.relative = 0.0` (no transcendentals —
//! strict byte equality). Arrays are stored as compact JSON strings; neighbor
//! ids are 0-based with a `-1` off-column sentinel and boundary masks are 0/1.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::grids::vertical_faq::vertical_construction_faq;
use earthsci_grids::vertical::{self, VerticalCoordinate, VerticalGrid};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/vertical")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

fn f64_array(v: &Value, field: &str) -> Vec<f64> {
    v.as_array()
        .unwrap_or_else(|| panic!("{field} must be an array"))
        .iter()
        .map(|x| x.as_f64().expect("numeric"))
        .collect()
}

fn build_grid(opts: &Value) -> VerticalGrid {
    let coordinate =
        VerticalCoordinate::from_name(opts["coordinate"].as_str().expect("opts.coordinate string"))
            .expect("known coordinate");
    let mut b = vertical::builder().coordinate(coordinate);
    if let Some(nz) = opts.get("nz") {
        b = b.nz(nz.as_u64().expect("nz integer") as usize);
    }
    if let Some(levels) = opts.get("levels") {
        b = b.levels(f64_array(levels, "levels"));
    }
    if let Some(ak) = opts.get("ak") {
        b = b.ak(f64_array(ak, "ak"));
    }
    if let Some(bk) = opts.get("bk") {
        b = b.bk(f64_array(bk, "bk"));
    }
    if let Some(p0) = opts.get("p0") {
        b = b.p0(p0.as_f64().expect("p0 numeric"));
    }
    b.build().expect("vertical fixture builds")
}

/// Parse a golden compact-JSON-string of a flat float array using Rust's
/// correctly-rounded `str::parse::<f64>()`. NOT `serde_json`, whose number
/// parser mis-rounds some 17-digit decimals (e.g. `0.9683012089999999`) by 1
/// ULP — which would spuriously fail this golden's strict `0.0` tolerance even
/// though the FAQ computation is bit-identical to the Julia reference.
fn parse_floats(s: &str) -> Vec<f64> {
    s.trim()
        .trim_start_matches('[')
        .trim_end_matches(']')
        .split(',')
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<f64>().expect("golden float token parses"))
        .collect()
}

/// Parse a golden compact-JSON-string into a flat `i64` array.
fn parse_ints(s: &str) -> Vec<i64> {
    let v: Value = serde_json::from_str(s).expect("golden int string parses");
    v.as_array()
        .expect("array")
        .iter()
        .map(|x| x.as_i64().expect("int"))
        .collect()
}

fn assert_floats_eq(got: &[f64], golden_str: &str, tol: f64, what: &str) {
    let g = parse_floats(golden_str);
    assert_eq!(got.len(), g.len(), "{what}: length");
    for (a, b) in got.iter().zip(g.iter()) {
        assert!(close_rel(*a, *b, tol), "{what}: {a} vs {b}");
    }
}

#[test]
fn vertical_construction_faq_matches_golden_and_imperative() {
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
        let faq = vertical_construction_faq(&grid);
        let gl = by_name
            .get(name)
            .unwrap_or_else(|| panic!("golden missing fixture {name}"));
        let n = faq.widths.len();

        // --- counts ---
        assert_eq!(n as u64, gl["n_cells"].as_u64().unwrap(), "{name}: n_cells");
        assert_eq!(
            faq.levels.len() as u64,
            gl["n_vertices"].as_u64().unwrap(),
            "{name}: n_vertices"
        );

        // --- parity: FAQ == imperative accessors, bit-for-bit ---
        assert_eq!(
            faq.levels.as_slice(),
            grid.levels(),
            "{name}: levels parity"
        );
        assert_eq!(
            faq.centers.as_slice(),
            grid.centers(),
            "{name}: centers parity"
        );
        assert_eq!(
            faq.widths.as_slice(),
            grid.widths(),
            "{name}: widths parity"
        );
        for k in 0..n {
            assert_eq!(
                faq.metric_dz[k],
                grid.metric_eval_by_name("dz", k).unwrap(),
                "{name}: metric dz parity"
            );
            assert_eq!(
                faq.metric_z[k],
                grid.metric_eval_by_name("z", k).unwrap(),
                "{name}: metric z parity"
            );
        }

        // --- byte identity against the cross-binding golden ---
        assert_floats_eq(
            &faq.levels,
            gl["levels"].as_str().unwrap(),
            rel_tol,
            "levels",
        );
        assert_floats_eq(
            &faq.centers,
            gl["centers"].as_str().unwrap(),
            rel_tol,
            "centers",
        );
        assert_floats_eq(
            &faq.widths,
            gl["widths"].as_str().unwrap(),
            rel_tol,
            "widths",
        );
        assert_floats_eq(
            &faq.cell_volume,
            gl["cell_volume"].as_str().unwrap(),
            rel_tol,
            "cell_volume",
        );
        assert_floats_eq(
            &faq.metric_dz,
            gl["metric_dz"].as_str().unwrap(),
            rel_tol,
            "metric_dz",
        );
        assert_floats_eq(
            &faq.metric_z,
            gl["metric_z"].as_str().unwrap(),
            rel_tol,
            "metric_z",
        );

        // Optional metrics — present in the golden iff the coordinate defines them.
        if let Some(ms) = gl.get("metric_sigma") {
            assert_floats_eq(
                faq.metric_sigma.as_ref().expect("faq metric_sigma"),
                ms.as_str().unwrap(),
                rel_tol,
                "metric_sigma",
            );
        }
        if let Some(mp) = gl.get("metric_pressure") {
            assert_floats_eq(
                faq.metric_pressure.as_ref().expect("faq metric_pressure"),
                mp.as_str().unwrap(),
                rel_tol,
                "metric_pressure",
            );
        }
        if let Some(ma) = gl.get("metric_ak") {
            assert_floats_eq(
                faq.metric_ak.as_ref().expect("faq metric_ak"),
                ma.as_str().unwrap(),
                rel_tol,
                "metric_ak",
            );
        }
        if let Some(mb) = gl.get("metric_bk") {
            assert_floats_eq(
                faq.metric_bk.as_ref().expect("faq metric_bk"),
                mb.as_str().unwrap(),
                rel_tol,
                "metric_bk",
            );
        }
        if let Some(ak) = gl.get("ak") {
            assert_floats_eq(&faq.ak, ak.as_str().unwrap(), rel_tol, "ak");
        }
        if let Some(bk) = gl.get("bk") {
            assert_floats_eq(&faq.bk, bk.as_str().unwrap(), rel_tol, "bk");
        }
        assert!(
            close_rel(faq.p0, gl["p0"].as_f64().unwrap(), rel_tol),
            "{name}: p0"
        );

        // --- structural integer arrays (exact) ---
        assert_eq!(
            faq.neighbor_minus,
            parse_ints(gl["neighbor_minus"].as_str().unwrap()),
            "{name}: neighbor_minus"
        );
        assert_eq!(
            faq.neighbor_plus,
            parse_ints(gl["neighbor_plus"].as_str().unwrap()),
            "{name}: neighbor_plus"
        );
        let bl: Vec<i64> = faq.boundary_lower.iter().map(|&x| x as i64).collect();
        let bu: Vec<i64> = faq.boundary_upper.iter().map(|&x| x as i64).collect();
        assert_eq!(
            bl,
            parse_ints(gl["boundary_lower"].as_str().unwrap()),
            "{name}: boundary_lower"
        );
        assert_eq!(
            bu,
            parse_ints(gl["boundary_upper"].as_str().unwrap()),
            "{name}: boundary_upper"
        );
    }
}
