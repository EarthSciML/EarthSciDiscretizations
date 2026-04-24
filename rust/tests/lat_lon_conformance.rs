//! Cross-language conformance test for the lat_lon grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs
//! matching the committed reference golden within the `docs/GRIDS_API.md`
//! §4.2 ULP tolerance.
//!
//! Julia is the nominal reference binding (§4.3) but its lat_lon runtime is
//! still tracked by bead dsc-1ts at the time of writing; the corpus is
//! produced by the Python binding as documented in
//! `tests/conformance/grids/latlon/README.md`. The contract (query points,
//! tolerances, neighbor encoding) is binding-neutral — the golden is just
//! committed output, regardless of which binding emitted it.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::lat_lon::{self, Direction, LatLonVariant, MetricName, PolePolicy};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR is .../rust; harness lives at ../tests/conformance/...
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/latlon")
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

fn golden_neighbor(value: &Value) -> Option<(usize, usize)> {
    if value.is_null() {
        return None;
    }
    let arr = value.as_array().unwrap();
    let gj = arr[0].as_u64().unwrap() as usize;
    let gi = arr[1].as_u64().unwrap() as usize;
    Some((gj, gi))
}

#[test]
fn lat_lon_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let metric_names = [
        "J",
        "g_lonlon",
        "g_latlat",
        "g_lonlat",
        "ginv_lonlon",
        "ginv_latlat",
        "ginv_lonlat",
    ];
    let dir_keys = ["W", "E", "S", "N"];
    let dirs = [Direction::W, Direction::E, Direction::S, Direction::N];

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let grid = build_grid(&fixture["opts"]);

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));
        assert_eq!(grid.n_cells() as u64, golden["n_cells"].as_u64().unwrap());

        let qps = fixture["query_points"].as_array().unwrap();
        for (k, qp) in qps.iter().enumerate() {
            let j = qp[0].as_u64().unwrap() as usize;
            let i = qp[1].as_u64().unwrap() as usize;

            let (lon, lat) = grid.cell_center(j, i).unwrap();
            let gc = &golden["cell_centers"][k];
            assert!(
                close_rel(lon, gc["lon"].as_f64().unwrap(), rel_tol),
                "{name}: lon at qp[{k}]=({j},{i}): {lon} vs {}",
                gc["lon"]
            );
            assert!(
                close_rel(lat, gc["lat"].as_f64().unwrap(), rel_tol),
                "{name}: lat at qp[{k}]=({j},{i}): {lat} vs {}",
                gc["lat"]
            );

            for (dkey, dir) in dir_keys.iter().zip(dirs.iter()) {
                let got = grid.neighbor(j, i, *dir).unwrap();
                let expected = golden_neighbor(&golden[format!("neighbors_{dkey}")][k]);
                assert_eq!(
                    got, expected,
                    "{name}: neighbor {dkey} at qp[{k}]=({j},{i})"
                );
            }

            for m in metric_names {
                let got = grid.metric_eval(metric_name(m), j, i).unwrap();
                let expected = golden[format!("metric_{m}")][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: metric {m} at qp[{k}]=({j},{i}): {got} vs {expected}"
                );
            }

            let got_area = grid.metric_eval(MetricName::Area, j, i).unwrap();
            let expected_area = golden["area"][k].as_f64().unwrap();
            assert!(
                close_rel(got_area, expected_area, rel_tol),
                "{name}: area at qp[{k}]=({j},{i}): {got_area} vs {expected_area}"
            );
        }
    }
}
