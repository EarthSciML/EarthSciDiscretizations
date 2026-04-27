//! Cross-binding rule-eval conformance test for the `ppm_reconstruction`
//! finite-volume rule.
//!
//! Drives `tests/conformance/rules/ppm_reconstruction/fixtures.json`,
//! which lists three layers of checks:
//!
//! 1. **`sub_stencil_cases`** — the named sub-stencils parse to the
//!    offsets and coefficient floats the harness pins. Coefficients are
//!    obtained by walking the rule AST with the canonical
//!    [`earthsci_grids::eval_coeff`] evaluator (no PPM-specific shadow
//!    parser).
//! 2. **`output_kind_cases`** — sub-stencil application
//!    (`Σ coeff_k · q(offset_k)`) reproduces hand-pinned values for the
//!    right-edge interpolation (CW84 eq. 1.6) and the parabola formula
//!    is hit at its endpoints.
//! 3. **`convergence_case`** — the canonical MMS sweep from
//!    `tests/fixtures/ppm_reconstruction/{input,expected}.esm`. The
//!    same inlined dispatch (sub_stencil → parabola) runs against the
//!    smooth-periodic profile the Julia reference uses, and the
//!    measured minimum convergence order must clear `expected_min_order`.
//!
//! The PPM-specific dispatch + parabola formula live here in the test —
//! per `AGENTS.md`, ESD does not carry per-rule shadow evaluators in its
//! public API.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::eval_coeff;
use serde_json::Value;

fn repo_root() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .into()
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0_f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

/// One entry in a PPM sub-stencil: cell-relative offset + coefficient
/// float resolved from the rule AST.
#[derive(Debug, Clone)]
struct StencilEntry {
    offset: i64,
    coeff: f64,
}

/// Walk `rule.discretizations.ppm_reconstruction.stencil[name]` and
/// resolve each entry's `selector.offset` (i64) and `coeff` AST (via
/// the canonical scalar evaluator) into a flat `(offset, coeff)` list.
fn parse_sub_stencil(rule_value: &Value, name: &str) -> Vec<StencilEntry> {
    let spec = rule_value
        .get("discretizations")
        .and_then(|d| d.get("ppm_reconstruction"))
        .unwrap_or_else(|| panic!("rule JSON missing discretizations.ppm_reconstruction"));
    let arr = spec
        .get("stencil")
        .and_then(|s| s.get(name))
        .and_then(Value::as_array)
        .unwrap_or_else(|| panic!("rule missing stencil.{name} array"));

    let bindings: HashMap<String, f64> = HashMap::new();
    arr.iter()
        .enumerate()
        .map(|(i, entry)| {
            let offset = entry
                .get("selector")
                .and_then(|s| s.get("offset"))
                .and_then(Value::as_i64)
                .unwrap_or_else(|| panic!("{name}[{i}]: selector.offset missing or not integer"));
            let coeff_node = entry
                .get("coeff")
                .unwrap_or_else(|| panic!("{name}[{i}]: missing coeff field"));
            let coeff = eval_coeff(coeff_node, &bindings)
                .unwrap_or_else(|e| panic!("{name}[{i}]: coeff AST eval failed: {e}"));
            StencilEntry { offset, coeff }
        })
        .collect()
}

/// Apply a sub-stencil at cell `i` of a periodic field via
/// `Σ coeff_k · q[(i + offset_k) mod n]`.
fn apply_sub_stencil<F>(entries: &[StencilEntry], q_at: F) -> f64
where
    F: Fn(i64) -> f64,
{
    entries.iter().map(|e| e.coeff * q_at(e.offset)).sum()
}

/// Parabolic reconstruction `a(ξ)` from CW84 eqs. (1.5), (1.7), (1.10).
/// Cell-average preservation: `(a_L + a_R)/2 + a_6/6 = q_i`.
fn parabola(a_l: f64, a_r: f64, q_i: f64, xi: f64) -> f64 {
    let da = a_r - a_l;
    let a6 = 6.0 * (q_i - 0.5 * (a_l + a_r));
    a_l + xi * (da + a6 * (1.0 - xi))
}

#[test]
fn ppm_rule_evals_match_conformance_fixture() {
    let root = repo_root();
    let spec = read_json(&root.join("tests/conformance/rules/ppm_reconstruction/fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let rule_path = spec["rule_path"].as_str().unwrap();
    let rule_value = read_json(&root.join(rule_path));
    let q_left_edge = parse_sub_stencil(&rule_value, "q_left_edge");
    let q_right_edge = parse_sub_stencil(&rule_value, "q_right_edge");

    // (1) sub_stencil structural cases: offsets + coefficient floats.
    for case in spec["sub_stencil_cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let sub_name = case["sub_stencil"].as_str().unwrap();
        let entries = match sub_name {
            "q_left_edge" => &q_left_edge,
            "q_right_edge" => &q_right_edge,
            other => panic!("unknown sub_stencil `{other}` in fixture"),
        };

        let expected_offsets: Vec<i64> = case["expected_offsets"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap())
            .collect();
        let expected_coeffs: Vec<f64> = case["expected_coeffs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();

        assert_eq!(
            entries.len(),
            expected_offsets.len(),
            "{name}: stencil length mismatch"
        );
        for (k, entry) in entries.iter().enumerate() {
            assert_eq!(
                entry.offset, expected_offsets[k],
                "{name}: offset[{k}] {} vs {}",
                entry.offset, expected_offsets[k]
            );
            assert!(
                close_rel(entry.coeff, expected_coeffs[k], rel_tol),
                "{name}: coeff[{k}] {} vs {}",
                entry.coeff,
                expected_coeffs[k]
            );
        }
    }

    // (2) output_kind dispatch cases.
    for case in spec["output_kind_cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let kind = case["output_kind"].as_str().unwrap();
        match kind {
            "q_right_edge" | "q_left_edge" => {
                let offsets: Vec<i64> = case["window_offsets"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_i64().unwrap())
                    .collect();
                let values: Vec<f64> = case["window_values"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap())
                    .collect();
                let q_at = |off: i64| {
                    let idx = offsets
                        .iter()
                        .position(|o| *o == off)
                        .unwrap_or_else(|| panic!("{name}: window has no offset {off}"));
                    values[idx]
                };
                let entries = if kind == "q_right_edge" {
                    &q_right_edge
                } else {
                    &q_left_edge
                };
                let got = apply_sub_stencil(entries, q_at);
                let expected = case["expected"].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: {got} vs {expected}"
                );
            }
            "q_parabola" => {
                for sample in case["samples"].as_array().unwrap() {
                    let a_l = sample["a_L"].as_f64().unwrap();
                    let a_r = sample["a_R"].as_f64().unwrap();
                    let q_i = sample["q_i"].as_f64().unwrap();
                    let xi = sample["xi"].as_f64().unwrap();
                    let expected = sample["expected"].as_f64().unwrap();
                    let got = parabola(a_l, a_r, q_i, xi);
                    assert!(
                        close_rel(got, expected, rel_tol),
                        "{name}: parabola(a_L={a_l}, a_R={a_r}, q_i={q_i}, xi={xi}) = {got} vs {expected}"
                    );
                }
            }
            other => panic!("{name}: unknown output_kind `{other}`"),
        }
    }
}

/// Manufactured solution from the canonical PPM fixture:
/// `f(x) = sin(2πx) + 0.3 cos(4πx)` and its analytical antiderivative
/// `F(x)` (used for exact per-cell averages — no quadrature error).
fn mms_f(x: f64) -> f64 {
    use std::f64::consts::PI;
    (2.0 * PI * x).sin() + 0.3 * (4.0 * PI * x).cos()
}

fn mms_big_f(x: f64) -> f64 {
    use std::f64::consts::PI;
    -(2.0 * PI * x).cos() / (2.0 * PI) + 0.3 * (4.0 * PI * x).sin() / (4.0 * PI)
}

fn cell_average(a: f64, b: f64) -> f64 {
    (mms_big_f(b) - mms_big_f(a)) / (b - a)
}

#[test]
fn ppm_convergence_sweep_matches_expected_min_order() {
    let root = repo_root();
    let spec = read_json(&root.join("tests/conformance/rules/ppm_reconstruction/fixtures.json"));
    let conv = &spec["convergence_case"];

    let input = read_json(&root.join(conv["fixture_input"].as_str().unwrap()));
    let expected = read_json(&root.join(conv["fixture_expected"].as_str().unwrap()));

    assert_eq!(input["rule"].as_str().unwrap(), "ppm_reconstruction");
    assert_eq!(expected["rule"].as_str().unwrap(), "ppm_reconstruction");

    let min_order = expected["expected_min_order"].as_f64().unwrap();
    let samples_per_cell = input["samples_per_cell"].as_u64().unwrap() as usize;
    let grids: Vec<usize> = input["grids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|g| g["n"].as_u64().unwrap() as usize)
        .collect();

    let rule_path = spec["rule_path"].as_str().unwrap();
    let rule_value = read_json(&root.join(rule_path));
    let q_left_edge = parse_sub_stencil(&rule_value, "q_left_edge");
    let q_right_edge = parse_sub_stencil(&rule_value, "q_right_edge");

    let linf_error = |n: usize| -> f64 {
        let dx = 1.0 / n as f64;
        // Cell averages on [(i)·dx, (i+1)·dx].
        let q: Vec<f64> = (0..n)
            .map(|i| cell_average(i as f64 * dx, (i + 1) as f64 * dx))
            .collect();
        let q_at_cell = |i: usize, off: i64| -> f64 {
            // Periodic wrap: i in [0, n), offset in [-2, 2].
            let idx = (i as i64 + off).rem_euclid(n as i64) as usize;
            q[idx]
        };

        let mut err: f64 = 0.0;
        for i in 0..n {
            let q_at = |off: i64| q_at_cell(i, off);
            let a_l = apply_sub_stencil(&q_left_edge, q_at);
            let a_r = apply_sub_stencil(&q_right_edge, q_at);
            let q_i = q_at(0);
            for k in 1..=samples_per_cell {
                let xi = (k as f64 - 0.5) / samples_per_cell as f64;
                let recon = parabola(a_l, a_r, q_i, xi);
                let x = i as f64 * dx + xi * dx;
                let truth = mms_f(x);
                err = err.max((recon - truth).abs());
            }
        }
        err
    };

    let errors: Vec<f64> = grids.iter().map(|&n| linf_error(n)).collect();
    for (n, e) in grids.iter().zip(errors.iter()) {
        assert!(
            e.is_finite() && *e > 0.0,
            "n={n}: error must be finite + positive, got {e}"
        );
    }

    // Errors decrease monotonically.
    for i in 0..(errors.len() - 1) {
        assert!(
            errors[i + 1] < errors[i],
            "non-monotone refinement: errors[{i}]={} >= errors[{}]={}",
            errors[i],
            i + 1,
            errors[i + 1]
        );
    }

    let orders: Vec<f64> = errors.windows(2).map(|w| (w[0] / w[1]).log2()).collect();
    let measured_min = orders.iter().copied().fold(f64::INFINITY, f64::min);

    assert!(
        measured_min >= min_order,
        "measured min convergence order {measured_min} < expected {min_order} \
         (errors: {errors:?}, orders: {orders:?})"
    );
}
