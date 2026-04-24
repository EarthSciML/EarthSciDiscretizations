//! Cross-language conformance test for the mpas grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs
//! matching the committed reference golden within the `docs/GRIDS_API.md`
//! §4.2 ULP tolerance. The mpas reference binding is interim Python until
//! bead dsc-7j0 (Julia mpas) lands — see
//! `tests/conformance/grids/mpas/README.md` for rationale. The tolerance
//! and compare-at-query-points protocol are unchanged.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::mpas::{
    self, MpasLoader, MpasMeshData, MpasMeshInput, MpasMetricName,
};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/mpas")
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

fn as_f64_vec(v: &Value, label: &str) -> Vec<f64> {
    v.as_array()
        .unwrap_or_else(|| panic!("{label}: not an array"))
        .iter()
        .map(|x| x.as_f64().unwrap_or_else(|| panic!("{label}: not a number")))
        .collect()
}

fn as_i32_vec(v: &Value, label: &str) -> Vec<i32> {
    v.as_array()
        .unwrap_or_else(|| panic!("{label}: not an array"))
        .iter()
        .map(|x| {
            x.as_i64()
                .unwrap_or_else(|| panic!("{label}: not an integer"))
                as i32
        })
        .collect()
}

fn as_usize(v: &Value, label: &str) -> usize {
    v.as_u64()
        .unwrap_or_else(|| panic!("{label}: not a non-negative integer")) as usize
}

fn build_mesh(fixture: &Value) -> MpasMeshData {
    let m = &fixture["mesh"];
    let r = fixture["opts"]["R"].as_f64().unwrap();
    let input = MpasMeshInput {
        lon_cell: as_f64_vec(&m["lon_cell"], "lon_cell"),
        lat_cell: as_f64_vec(&m["lat_cell"], "lat_cell"),
        area_cell: as_f64_vec(&m["area_cell"], "area_cell"),
        n_edges_on_cell: as_i32_vec(&m["n_edges_on_cell"], "n_edges_on_cell"),
        cells_on_cell: as_i32_vec(&m["cells_on_cell"], "cells_on_cell"),
        edges_on_cell: as_i32_vec(&m["edges_on_cell"], "edges_on_cell"),
        lon_edge: as_f64_vec(&m["lon_edge"], "lon_edge"),
        lat_edge: as_f64_vec(&m["lat_edge"], "lat_edge"),
        cells_on_edge: as_i32_vec(&m["cells_on_edge"], "cells_on_edge"),
        dc_edge: as_f64_vec(&m["dc_edge"], "dc_edge"),
        dv_edge: as_f64_vec(&m["dv_edge"], "dv_edge"),
        max_edges: as_usize(&m["max_edges"], "max_edges"),
        x_cell: Some(as_f64_vec(&m["x_cell"], "x_cell")),
        y_cell: Some(as_f64_vec(&m["y_cell"], "y_cell")),
        z_cell: Some(as_f64_vec(&m["z_cell"], "z_cell")),
        n_vertices: Some(as_usize(&m["n_vertices"], "n_vertices")),
        r: Some(r),
    };
    MpasMeshData::from_input(input).expect("valid mesh")
}

const CELL_METRIC_NAMES: &[(&str, MpasMetricName)] = &[
    ("lon", MpasMetricName::Lon),
    ("lat", MpasMetricName::Lat),
    ("area", MpasMetricName::Area),
    ("x", MpasMetricName::X),
    ("y", MpasMetricName::Y),
    ("z", MpasMetricName::Z),
    ("n_edges_on_cell", MpasMetricName::NEdgesOnCell),
];

const EDGE_METRIC_NAMES: &[(&str, MpasMetricName)] = &[
    ("lon_edge", MpasMetricName::LonEdge),
    ("lat_edge", MpasMetricName::LatEdge),
    ("dc_edge", MpasMetricName::DcEdge),
    ("dv_edge", MpasMetricName::DvEdge),
];

#[test]
fn mpas_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let opts = &fixture["opts"];
        let r = opts["R"].as_f64().unwrap();
        let ghosts = opts["ghosts"].as_u64().unwrap() as u32;

        let mesh = build_mesh(fixture);
        // Use a lenient loader so the harness accepts asymmetric future
        // fixtures without reissuing golden. The tetra mesh is symmetric,
        // so strict would also pass.
        let loader = MpasLoader {
            path: format!("fixture://{name}"),
            reader: "mpas_mesh".into(),
            check: "lenient".into(),
        };
        let grid = mpas::builder()
            .mesh(mesh)
            .loader(loader)
            .r(r)
            .ghosts(ghosts)
            .build()
            .unwrap_or_else(|e| panic!("{name}: build failed: {e:?}"));

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        assert_eq!(
            grid.n_cells() as u64,
            golden["n_cells"].as_u64().unwrap(),
            "{name}: n_cells mismatch"
        );
        assert_eq!(
            grid.n_edges() as u64,
            golden["n_edges"].as_u64().unwrap(),
            "{name}: n_edges mismatch"
        );
        assert_eq!(
            grid.max_edges() as u64,
            golden["max_edges"].as_u64().unwrap(),
            "{name}: max_edges mismatch"
        );

        let query_cells = fixture["query_cells"].as_array().unwrap();
        for (k, qc) in query_cells.iter().enumerate() {
            let c = qc.as_u64().unwrap() as usize;

            let (lon, lat) = grid.cell_centers(c).unwrap();
            let gc = &golden["cell_centers"][k];
            assert!(
                close_rel(lon, gc["lon"].as_f64().unwrap(), rel_tol),
                "{name}: lon mismatch at cell {c}"
            );
            assert!(
                close_rel(lat, gc["lat"].as_f64().unwrap(), rel_tol),
                "{name}: lat mismatch at cell {c}"
            );

            let cart = grid.cell_center_cart(c).unwrap();
            let gx = &golden["cell_centers_cart"][k];
            assert!(
                close_rel(cart[0], gx["x"].as_f64().unwrap(), rel_tol),
                "{name}: cart x mismatch at cell {c}"
            );
            assert!(
                close_rel(cart[1], gx["y"].as_f64().unwrap(), rel_tol),
                "{name}: cart y mismatch at cell {c}"
            );
            assert!(
                close_rel(cart[2], gx["z"].as_f64().unwrap(), rel_tol),
                "{name}: cart z mismatch at cell {c}"
            );

            let got_nb: Vec<i64> = grid
                .neighbors(c)
                .unwrap()
                .into_iter()
                .map(|n| n as i64)
                .collect();
            let exp_nb: Vec<i64> = golden["neighbors"][k]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(got_nb, exp_nb, "{name}: neighbors mismatch at cell {c}");

            let got_area = grid.cell_area(c).unwrap();
            let exp_area = golden["cell_area"][k].as_f64().unwrap();
            assert!(
                close_rel(got_area, exp_area, rel_tol),
                "{name}: cell_area mismatch at cell {c}"
            );

            for (key, metric) in CELL_METRIC_NAMES {
                let got = grid.metric_eval(*metric, c).unwrap();
                let expected = golden["cell_metrics"][*key][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: cell metric {key} mismatch at cell {c}: {got} vs {expected}"
                );
            }
        }

        let query_edges = fixture["query_edges"].as_array().unwrap();
        for (k, qe) in query_edges.iter().enumerate() {
            let e = qe.as_u64().unwrap() as usize;
            let got_len = grid.edge_length(e).unwrap();
            let exp_len = golden["edge_length"][k].as_f64().unwrap();
            assert!(
                close_rel(got_len, exp_len, rel_tol),
                "{name}: edge_length mismatch at edge {e}"
            );
            let got_dist = grid.cell_distance(e).unwrap();
            let exp_dist = golden["cell_distance"][k].as_f64().unwrap();
            assert!(
                close_rel(got_dist, exp_dist, rel_tol),
                "{name}: cell_distance mismatch at edge {e}"
            );
            for (key, metric) in EDGE_METRIC_NAMES {
                let got = grid.metric_eval(*metric, e).unwrap();
                let expected = golden["edge_metrics"][*key][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: edge metric {key} mismatch at edge {e}: {got} vs {expected}"
                );
            }
        }
    }
}
