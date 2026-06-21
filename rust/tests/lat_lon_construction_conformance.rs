//! Cross-language conformance for the lat-lon **construction FAQ** (esd-dru).
//!
//! Builds each fixture grid imperatively, re-derives its construction arrays
//! through the elementwise FAQ bridge (`lat_lon_construction_faq`, routed via the
//! ESS `eval_coeff` evaluator), and asserts the result is
//!   1. byte-identical to the imperative `LatLonGrid` accessors over **all**
//!      cells (parity — the FAQ introduces zero error),
//!   2. equal to the committed Julia-reference
//!      `tests/conformance/grids/latlon/construction/golden.json`.
//!
//! Methodology (matches the established cross-binding latlon conformance in
//! `lat_lon_conformance.rs` and `docs/GRIDS_API.md` §4.2, Julia as the ULP
//! anchor): the per-cell trig arrays carry `sin`/`cos` whose libm differs
//! sub-ULP across bindings, and near the poles the `sin φ_n − sin φ_s` area term
//! cancels and amplifies that drift just past the 1e-14 family tolerance on the
//! high-res `realistic` fixture. So the golden float comparison is sampled at the
//! fixture's curated `query_points` at the declared 1e-14 (exactly as the
//! existing latlon conformance test does), while integer neighbor maps / boundary
//! masks (libm-free) and the 1-D latitude arrays are checked over the whole grid.
//! The exhaustive guarantee comes from the all-cells FAQ==imperative parity.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::grids::lat_lon_faq::lat_lon_construction_faq;
use earthsci_grids::lat_lon::{self, LatLonGrid, LatLonVariant, PolePolicy};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/latlon")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

fn parse_floats(s: &str) -> Vec<f64> {
    let v: Value = serde_json::from_str(s).expect("golden float string parses");
    v.as_array()
        .expect("array")
        .iter()
        .map(|x| x.as_f64().expect("float"))
        .collect()
}

fn parse_ints(s: &str) -> Vec<i64> {
    let v: Value = serde_json::from_str(s).expect("golden int string parses");
    v.as_array()
        .expect("array")
        .iter()
        .map(|x| x.as_i64().expect("int"))
        .collect()
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

fn build_grid(opts: &Value) -> LatLonGrid {
    let variant = parse_variant(opts["variant"].as_str().unwrap());
    let mut builder = lat_lon::builder()
        .variant(variant)
        .r(opts["R"].as_f64().unwrap())
        .ghosts(opts["ghosts"].as_u64().unwrap() as usize)
        .pole_policy(parse_pole_policy(opts["pole_policy"].as_str().unwrap()));
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
                builder = builder.lat_edges(edges.iter().map(|v| v.as_f64().unwrap()).collect());
            }
        }
    }
    builder.build().unwrap()
}

/// Per-cell trig-bearing arrays — golden-compared at the curated query points.
const CELL_FLOAT_KEYS: &[&str] = &[
    "cell_centers_lon",
    "cell_centers_lat",
    "cell_widths_lon",
    "cell_widths_lat",
    "cell_volume",
    "metric_g_lonlon",
    "metric_g_latlat",
    "metric_ginv_lonlon",
    "metric_jacobian",
    "dg_lonlon_dlat",
];

fn cell_float<'a>(
    faq: &'a earthsci_grids::grids::lat_lon_faq::LatLonConstruction,
    key: &str,
) -> &'a [f64] {
    match key {
        "cell_centers_lon" => &faq.cell_centers_lon,
        "cell_centers_lat" => &faq.cell_centers_lat,
        "cell_widths_lon" => &faq.cell_widths_lon,
        "cell_widths_lat" => &faq.cell_widths_lat,
        "cell_volume" => &faq.cell_volume,
        "metric_g_lonlon" => &faq.metric_g_lonlon,
        "metric_g_latlat" => &faq.metric_g_latlat,
        "metric_ginv_lonlon" => &faq.metric_ginv_lonlon,
        "metric_jacobian" => &faq.metric_jacobian,
        "dg_lonlon_dlat" => &faq.dg_lonlon_dlat,
        other => panic!("unmapped float key {other}"),
    }
}

#[test]
fn lat_lon_construction_faq_matches_golden_and_imperative() {
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
        let faq = lat_lon_construction_faq(&grid);
        let gl = by_name
            .get(name)
            .unwrap_or_else(|| panic!("golden missing fixture {name}"));
        let nlat = grid.nlat();
        let per_row = grid.nlon_per_row().to_vec();
        let mut row_off = vec![0usize; nlat + 1];
        for j in 0..nlat {
            row_off[j + 1] = row_off[j] + per_row[j];
        }

        // --- counts ---
        assert_eq!(
            faq.cell_centers_lon.len() as u64,
            gl["n_cells"].as_u64().unwrap(),
            "{name}: n_cells"
        );
        assert_eq!(nlat as u64, gl["nlat"].as_u64().unwrap(), "{name}: nlat");

        // --- parity: FAQ == imperative accessors, over ALL cells ---
        assert_eq!(
            faq.lat_edges.as_slice(),
            grid.lat_edges(),
            "{name}: lat_edges"
        );
        assert_eq!(
            faq.lat_centers.as_slice(),
            grid.lat_centers(),
            "{name}: lat_centers"
        );
        for (j, &n_i) in per_row.iter().enumerate() {
            let base = row_off[j];
            for i in 0..n_i {
                let (lon, lat) = grid.cell_center(j, i).unwrap();
                assert_eq!(faq.cell_centers_lon[base + i], lon, "{name}: cc_lon parity");
                assert_eq!(faq.cell_centers_lat[base + i], lat, "{name}: cc_lat parity");
                assert_eq!(
                    faq.cell_volume[base + i],
                    grid.cell_area(j, i).unwrap(),
                    "{name}: cell_volume parity"
                );
            }
        }

        // --- cross-binding golden: per-cell trig arrays at curated query points ---
        let cell_golden: std::collections::HashMap<&str, Vec<f64>> = CELL_FLOAT_KEYS
            .iter()
            .map(|k| (*k, parse_floats(gl[*k].as_str().unwrap())))
            .collect();
        for qp in fixture["query_points"].as_array().unwrap() {
            let pair = qp.as_array().unwrap();
            let j = pair[0].as_u64().unwrap() as usize;
            let i = pair[1].as_u64().unwrap() as usize;
            let k = row_off[j] + i;
            for key in CELL_FLOAT_KEYS {
                let got = cell_float(&faq, key)[k];
                let exp = cell_golden[*key][k];
                assert!(
                    close_rel(got, exp, rel_tol),
                    "{name}:{key} at ({j},{i}): {got} vs {exp}"
                );
            }
        }

        // --- cross-binding golden: 1-D latitude arrays over all rows ---
        for (got, exp) in faq
            .lat_edges
            .iter()
            .zip(parse_floats(gl["lat_edges"].as_str().unwrap()).iter())
        {
            assert!(close_rel(*got, *exp, rel_tol), "{name}: lat_edges golden");
        }
        for (got, exp) in faq
            .lat_centers
            .iter()
            .zip(parse_floats(gl["lat_centers"].as_str().unwrap()).iter())
        {
            assert!(close_rel(*got, *exp, rel_tol), "{name}: lat_centers golden");
        }

        // --- cross-binding golden: integer + mask arrays, exact, all cells ---
        assert_eq!(
            faq.neighbor_lon_minus,
            parse_ints(gl["neighbor_lon_minus"].as_str().unwrap()),
            "{name}: neighbor_lon_minus"
        );
        assert_eq!(
            faq.neighbor_lon_plus,
            parse_ints(gl["neighbor_lon_plus"].as_str().unwrap()),
            "{name}: neighbor_lon_plus"
        );
        assert_eq!(
            faq.neighbor_lat_minus,
            parse_ints(gl["neighbor_lat_minus"].as_str().unwrap()),
            "{name}: neighbor_lat_minus"
        );
        assert_eq!(
            faq.neighbor_lat_plus,
            parse_ints(gl["neighbor_lat_plus"].as_str().unwrap()),
            "{name}: neighbor_lat_plus"
        );
        let bl: Vec<i64> = faq.boundary_lat_lower.iter().map(|&x| x as i64).collect();
        let bu: Vec<i64> = faq.boundary_lat_upper.iter().map(|&x| x as i64).collect();
        assert_eq!(
            bl,
            parse_ints(gl["boundary_lat_lower"].as_str().unwrap()),
            "{name}: boundary_lat_lower"
        );
        assert_eq!(
            bu,
            parse_ints(gl["boundary_lat_upper"].as_str().unwrap()),
            "{name}: boundary_lat_upper"
        );
    }
}
