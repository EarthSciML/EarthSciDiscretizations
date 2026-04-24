//! Cross-language conformance test for the DUO grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs
//! matching the committed Julia-reference golden within the
//! `docs/GRIDS_API.md` §4.2 ULP tolerance. Corpus lives at
//! `tests/conformance/grids/duo/`.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::duo::{self, DuoLoader};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR is .../rust; harness lives at ../tests/conformance/...
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/duo")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

#[test]
fn duo_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let metric_names = ["area", "lon", "lat", "x", "y", "z"];

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let opts = &fixture["opts"];
        let level = opts["level"].as_u64().unwrap() as u32;
        let r = opts["R"].as_f64().unwrap();
        let ghosts = opts["ghosts"].as_u64().unwrap() as u32;

        let grid = duo::builder()
            .loader(DuoLoader::builtin_level(level))
            .r(r)
            .ghosts(ghosts)
            .build()
            .unwrap();

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        assert_eq!(grid.n_cells() as u64, golden["n_cells"].as_u64().unwrap());
        assert_eq!(
            grid.n_vertices() as u64,
            golden["n_vertices"].as_u64().unwrap()
        );
        assert_eq!(grid.n_edges() as u64, golden["n_edges"].as_u64().unwrap());

        let qps = fixture["query_points"].as_array().unwrap();
        assert_eq!(qps.len(), golden["cell_centers"].as_array().unwrap().len());

        for (k, qp) in qps.iter().enumerate() {
            let c = qp.as_u64().unwrap() as usize;

            let (lon, lat) = grid.cell_centers(c);
            let gc = &golden["cell_centers"][k];
            assert!(
                close_rel(lon, gc["lon"].as_f64().unwrap(), rel_tol),
                "{name}: lon mismatch at qp[{k}]=c{c}: {lon} vs {}",
                gc["lon"]
            );
            assert!(
                close_rel(lat, gc["lat"].as_f64().unwrap(), rel_tol),
                "{name}: lat mismatch at qp[{k}]=c{c}: {lat} vs {}",
                gc["lat"]
            );

            // Rust returns [Option<u32>; 3]; golden uses -1 for "no neighbor".
            let got_nb: Vec<i64> = grid
                .neighbors(c)
                .iter()
                .map(|n| match n {
                    Some(v) => *v as i64,
                    None => -1,
                })
                .collect();
            let expected_nb: Vec<i64> = golden["neighbors"][k]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(got_nb, expected_nb, "{name}: neighbors mismatch at qp[{k}]=c{c}");

            for m in metric_names {
                let got = grid.metric_eval(m, c).unwrap();
                let expected = golden[format!("metric_{m}")][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: metric {m} mismatch at qp[{k}]=c{c}: {got} vs {expected}"
                );
            }
        }
    }
}
