//! Cross-language conformance test for the cubed_sphere grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs
//! matching the committed Julia-reference golden within the
//! `docs/GRIDS_API.md` §4.2 ULP tolerance.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::cubed_sphere::{self, Direction, MetricName};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR is .../rust; harness lives at ../tests/conformance/...
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/cubed_sphere")
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

fn metric_name(name: &str) -> MetricName {
    MetricName::from_name(name).unwrap_or_else(|| panic!("unknown metric {name}"))
}

#[test]
fn cubed_sphere_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let metric_names = [
        "J",
        "g_xixi",
        "g_etaeta",
        "g_xieta",
        "ginv_xixi",
        "ginv_etaeta",
        "ginv_xieta",
    ];
    let dir_keys = ["W", "E", "S", "N"];
    let dirs = [Direction::W, Direction::E, Direction::S, Direction::N];

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let opts = &fixture["opts"];
        let nc = opts["Nc"].as_u64().unwrap() as usize;
        let r = opts["R"].as_f64().unwrap();
        let ghosts = opts["ghosts"].as_u64().unwrap() as usize;

        let grid = cubed_sphere::builder()
            .nc(nc)
            .r(r)
            .ghosts(ghosts)
            .build()
            .unwrap();

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        assert_eq!(grid.n_cells() as u64, golden["n_cells"].as_u64().unwrap());

        let qps = fixture["query_points"].as_array().unwrap();
        for (k, qp) in qps.iter().enumerate() {
            let p = qp[0].as_u64().unwrap() as usize;
            let i = qp[1].as_u64().unwrap() as usize;
            let j = qp[2].as_u64().unwrap() as usize;

            let (lon, lat) = grid.cell_centers(p, i, j).unwrap();
            let gc = &golden["cell_centers"][k];
            assert!(
                close_rel(lon, gc["lon"].as_f64().unwrap(), rel_tol),
                "{name}: lon mismatch at qp[{k}]=({p},{i},{j}): {lon} vs {}",
                gc["lon"]
            );
            assert!(
                close_rel(lat, gc["lat"].as_f64().unwrap(), rel_tol),
                "{name}: lat mismatch at qp[{k}]=({p},{i},{j}): {lat} vs {}",
                gc["lat"]
            );

            for (dkey, dir) in dir_keys.iter().zip(dirs.iter()) {
                let (np, ni, nj) = grid.neighbor(p, i, j, *dir).unwrap();
                let golden_nb = &golden[format!("neighbors_{dkey}")][k];
                let gp = golden_nb[0].as_u64().unwrap() as usize;
                let gi = golden_nb[1].as_u64().unwrap() as usize;
                let gj = golden_nb[2].as_u64().unwrap() as usize;
                assert_eq!(
                    (np, ni, nj),
                    (gp, gi, gj),
                    "{name}: neighbor {dkey} mismatch at qp[{k}]=({p},{i},{j})"
                );
            }

            for m in metric_names {
                let got = grid.metric_eval(metric_name(m), p, i, j).unwrap();
                let expected = golden[format!("metric_{m}")][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: metric {m} mismatch at qp[{k}]=({p},{i},{j}): {got} vs {expected}"
                );
            }

            let got_area = grid.metric_eval(MetricName::Area, p, i, j).unwrap();
            let expected_area = golden["area"][k].as_f64().unwrap();
            assert!(
                close_rel(got_area, expected_area, rel_tol),
                "{name}: area mismatch at qp[{k}]=({p},{i},{j}): {got_area} vs {expected_area}"
            );
        }
    }
}
