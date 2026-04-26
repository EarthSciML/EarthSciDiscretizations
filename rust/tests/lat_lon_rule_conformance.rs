//! Cross-binding rule-eval conformance test for lat-lon FD coefficients.
//!
//! Verifies that the Rust ESS-AST evaluator (`earthsci_grids::eval_coeff`)
//! reproduces the golden coefficient values published under
//! `tests/conformance/grids/latlon/golden/<fixture>.json` → `rule_evals`.
//! The reference golden is authored by `regenerate_golden.py` per
//! `docs/GRIDS_API.md` §4.3 (Python is the temporary reference binding;
//! Julia takes over once dsc-1ts lands the Julia lat-lon accessor) and
//! the cross-binding agreement tolerance lives in `fixtures.json`
//! → `tolerance.relative`.
//!
//! What this test exercises:
//!  1. The lat-lon grid accessor surfaces the per-cell bindings
//!     (`R`, `cos_lat`, `dlon`, `dlat`) needed by the
//!     `centered_2nd_uniform_latlon` rule, byte-for-byte against the
//!     reference.
//!  2. `earthsci_grids::eval_coeff` walks the ESS expression-AST and
//!     evaluates each stencil entry's `coeff` to the same float as the
//!     reference evaluator.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::{eval_coeff, lat_lon, Bindings};
use lat_lon::{LatLonVariant, PolePolicy};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/latlon")
}

fn repo_root() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .into()
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("failed to read {}: {e}", path.display());
    });
    serde_json::from_str(&text).unwrap_or_else(|e| {
        panic!("failed to parse {}: {e}", path.display());
    })
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

fn parse_variant(s: &str) -> LatLonVariant {
    match s {
        "regular" => LatLonVariant::Regular,
        "reduced_gaussian" => LatLonVariant::ReducedGaussian,
        other => panic!("unknown variant: {other}"),
    }
}

fn parse_pole_policy(s: &str) -> PolePolicy {
    match s {
        "none" => PolePolicy::None,
        "average" => PolePolicy::Average,
        "fold" => PolePolicy::Fold,
        other => panic!("unknown pole policy: {other}"),
    }
}

fn build_grid(opts: &Value) -> lat_lon::LatLonGrid {
    let variant = parse_variant(opts["variant"].as_str().unwrap());
    let r = opts["R"].as_f64().unwrap();
    let ghosts = opts["ghosts"].as_u64().unwrap() as usize;
    let pole_policy = parse_pole_policy(opts["pole_policy"].as_str().unwrap());

    let mut builder = lat_lon::builder()
        .variant(variant)
        .r(r)
        .ghosts(ghosts)
        .pole_policy(pole_policy);

    match variant {
        LatLonVariant::Regular => {
            builder = builder
                .nlon(opts["nlon"].as_u64().unwrap() as usize)
                .nlat(opts["nlat"].as_u64().unwrap() as usize);
        }
        LatLonVariant::ReducedGaussian => {
            let nlon_per_row: Vec<usize> = opts["nlon_per_row"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_u64().unwrap() as usize)
                .collect();
            builder = builder.nlon_per_row(nlon_per_row);
            if let Some(edges) = opts.get("lat_edges").and_then(|v| v.as_array()) {
                let lat_edges: Vec<f64> = edges.iter().map(|v| v.as_f64().unwrap()).collect();
                builder = builder.lat_edges(lat_edges);
            }
        }
    }

    builder.build().unwrap()
}

/// Build the `centered_2nd_uniform_latlon` bindings for cell `(j, i)`.
///
/// Mirrors `_bindings_for` in `regenerate_golden.py`: `dlon` is the
/// per-row angular cell width (`2π / nlon_per_row[j]`); `dlat` is
/// `lat_edges[j+1] - lat_edges[j]` (uniform for `Regular`, ragged for
/// `ReducedGaussian`).
fn bindings_for(grid: &lat_lon::LatLonGrid, j: usize, i: usize) -> Bindings {
    let (_lon, lat) = grid.cell_center(j, i).expect("cell center");
    let nlon = grid.nlon(j).expect("nlon for row");
    let dlon = 2.0 * std::f64::consts::PI / nlon as f64;
    let edges = grid.lat_edges();
    let dlat = edges[j + 1] - edges[j];
    let mut b: HashMap<String, f64> = HashMap::new();
    b.insert("R".into(), grid.r());
    b.insert("cos_lat".into(), lat.cos());
    b.insert("dlon".into(), dlon);
    b.insert("dlat".into(), dlat);
    b
}

fn load_rule(rel_path: &str) -> (String, Value) {
    let path = repo_root().join(rel_path);
    let value = read_json(&path);
    let discs = value["discretizations"]
        .as_object()
        .unwrap_or_else(|| panic!("{}: missing `discretizations` object", path.display()));
    if discs.len() != 1 {
        panic!(
            "{}: expected exactly one discretization, got {:?}",
            path.display(),
            discs.keys().collect::<Vec<_>>()
        );
    }
    let (name, rule) = discs.iter().next().unwrap();
    (name.clone(), rule.clone())
}

#[test]
fn lat_lon_rule_evals_match_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let grid = build_grid(&fixture["opts"]);
        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        let qps = fixture["query_points"].as_array().unwrap();

        let rule_blocks = golden["rule_evals"].as_array().unwrap_or_else(|| {
            panic!("{name}: golden missing `rule_evals` array (regenerate the golden)")
        });
        assert!(
            !rule_blocks.is_empty(),
            "{name}: golden `rule_evals` is empty"
        );

        for block in rule_blocks {
            let rule_path = block["rule_path"].as_str().unwrap();
            let (rule_name, rule) = load_rule(rule_path);
            assert_eq!(
                rule_name,
                block["rule"].as_str().unwrap(),
                "{name}: rule name mismatch for {rule_path}"
            );

            let stencil = rule["stencil"].as_array().unwrap();
            let golden_coeffs = block["stencil_coeffs"].as_array().unwrap();
            let golden_bindings = block["bindings_per_qp"].as_array().unwrap();

            assert_eq!(
                qps.len(),
                golden_coeffs.len(),
                "{name}/{rule_name}: golden coeff rows != query points"
            );

            for (k, qp) in qps.iter().enumerate() {
                let j = qp[0].as_u64().unwrap() as usize;
                let i = qp[1].as_u64().unwrap() as usize;
                let bindings = bindings_for(&grid, j, i);

                // (1) Bindings extracted from the Rust grid match those the
                //     reference binding wrote into the golden.
                let g_b = golden_bindings[k].as_object().unwrap();
                for (key, val) in g_b {
                    let expected = val.as_f64().unwrap();
                    let got = bindings.get(key).copied().unwrap_or_else(|| {
                        panic!("{name}/{rule_name} qp[{k}]=({j},{i}): missing binding `{key}`")
                    });
                    assert!(
                        close_rel(got, expected, rel_tol),
                        "{name}/{rule_name} qp[{k}]=({j},{i}): binding `{key}` {got} vs {expected}"
                    );
                }

                // (2) Each stencil entry's AST coefficient evaluates to the
                //     same float as the reference walker emitted.
                let golden_row = golden_coeffs[k].as_array().unwrap();
                assert_eq!(
                    stencil.len(),
                    golden_row.len(),
                    "{name}/{rule_name} qp[{k}]: stencil length != golden coeff length"
                );
                for (s_idx, entry) in stencil.iter().enumerate() {
                    let got = eval_coeff(&entry["coeff"], &bindings).unwrap_or_else(|e| {
                        panic!(
                            "{name}/{rule_name} qp[{k}]=({j},{i}) stencil[{s_idx}]: \
                             eval_coeff failed: {e}"
                        )
                    });
                    let expected = golden_row[s_idx].as_f64().unwrap();
                    assert!(
                        close_rel(got, expected, rel_tol),
                        "{name}/{rule_name} qp[{k}]=({j},{i}) stencil[{s_idx}]: \
                         eval {got} vs golden {expected}"
                    );
                }
            }
        }
    }
}
