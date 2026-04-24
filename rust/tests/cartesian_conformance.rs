//! Cross-language conformance test for the cartesian grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs
//! matching the committed Julia-reference golden within the
//! `docs/GRIDS_API.md` §4.2 tolerance. Corpus lives at
//! `tests/conformance/grids/cartesian/`.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::cartesian::{self, CartesianGrid, MetricName};
use earthsci_grids::Dtype;
use serde_json::Value;

fn harness_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR is .../rust; harness lives at ../tests/conformance/...
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/cartesian")
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

fn qp_to_indices(qp: &Value) -> Vec<usize> {
    qp.as_array()
        .expect("query point is an array")
        .iter()
        .map(|v| v.as_u64().expect("qp index is a non-negative integer") as usize)
        .collect()
}

#[test]
fn cartesian_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let width_names = ["dx", "dy", "dz"];
    let face_names = ["face_area_x", "face_area_y", "face_area_z"];

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let grid = build_grid(&fixture["opts"]);
        let ndim = grid.ndim();

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        assert_eq!(ndim, golden["ndim"].as_u64().unwrap() as usize);
        assert_eq!(grid.n_cells() as u64, golden["n_cells"].as_u64().unwrap());

        let qps = fixture["query_points"].as_array().unwrap();
        assert_eq!(qps.len(), golden["cell_centers"].as_array().unwrap().len());

        for (k, qp) in qps.iter().enumerate() {
            let idx = qp_to_indices(qp);

            let c = grid.cell_center(&idx).unwrap();
            let gc = golden["cell_centers"][k].as_array().unwrap();
            for d in 0..ndim {
                assert!(
                    close_rel(c[d], gc[d].as_f64().unwrap(), rel_tol),
                    "{name}: cell_centers[{d}] at qp[{k}]={idx:?}"
                );
            }

            let w = grid.cell_widths(&idx).unwrap();
            let gw = golden["cell_widths"][k].as_array().unwrap();
            for d in 0..ndim {
                assert!(
                    close_rel(w[d], gw[d].as_f64().unwrap(), rel_tol),
                    "{name}: cell_widths[{d}] at qp[{k}]={idx:?}"
                );
            }

            let v = grid.cell_volume(&idx).unwrap();
            assert!(
                close_rel(v, golden["cell_volume"][k].as_f64().unwrap(), rel_tol),
                "{name}: cell_volume at qp[{k}]={idx:?}"
            );

            // Neighbors already come back in axis-major, low-before-high
            // order from the Rust accessor; the golden uses the same layout.
            let nbrs = grid.neighbors(&idx).unwrap();
            let gnbrs = golden["neighbors"][k].as_array().unwrap();
            assert_eq!(
                nbrs.len(),
                gnbrs.len(),
                "{name}: neighbor count at qp[{k}]={idx:?}"
            );
            for (got, exp) in nbrs.iter().zip(gnbrs.iter()) {
                assert_eq!(got.axis as u64, exp["axis"].as_u64().unwrap());
                assert_eq!(got.side as i64, exp["side"].as_i64().unwrap());
                let eidx: Vec<usize> = exp["index"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_u64().unwrap() as usize)
                    .collect();
                assert_eq!(got.index, eidx);
            }

            // Scalar metrics
            let got_vol = grid.metric_eval(MetricName::Volume, &idx).unwrap();
            assert!(close_rel(
                got_vol,
                golden["metric_volume"][k].as_f64().unwrap(),
                rel_tol,
            ));
            let got_jac = grid.metric_eval(MetricName::Jacobian, &idx).unwrap();
            assert!(close_rel(
                got_jac,
                golden["metric_jacobian"][k].as_f64().unwrap(),
                rel_tol,
            ));

            for d in 0..ndim {
                let wgot = grid.metric_eval_by_name(width_names[d], &idx).unwrap();
                let wexp = golden[format!("metric_{}", width_names[d])][k]
                    .as_f64()
                    .unwrap();
                assert!(
                    close_rel(wgot, wexp, rel_tol),
                    "{name}: metric {} at qp[{k}]={idx:?}: {wgot} vs {wexp}",
                    width_names[d]
                );
                let fgot = grid.metric_eval_by_name(face_names[d], &idx).unwrap();
                let fexp = golden[format!("metric_{}", face_names[d])][k]
                    .as_f64()
                    .unwrap();
                assert!(
                    close_rel(fgot, fexp, rel_tol),
                    "{name}: metric {} at qp[{k}]={idx:?}: {fgot} vs {fexp}",
                    face_names[d]
                );
            }

            let g = grid.metric_g(&idx).unwrap();
            let gg = golden["metric_g"][k].as_array().unwrap();
            for i in 0..ndim {
                let row = gg[i].as_array().unwrap();
                for j in 0..ndim {
                    assert!(
                        close_rel(g[i][j], row[j].as_f64().unwrap(), rel_tol),
                        "{name}: metric_g[{i}][{j}] at qp[{k}]={idx:?}"
                    );
                }
            }
        }
    }
}
