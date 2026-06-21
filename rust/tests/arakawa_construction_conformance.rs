//! Cross-language conformance for the arakawa **construction FAQ** (esd-dru).
//!
//! Builds each fixture's cartesian base, re-derives the staggered construction
//! through the elementwise FAQ bridge (`arakawa_construction_faq`, routed via the
//! ESS `eval_coeff` evaluator), and asserts the result is
//!   1. byte-identical to the imperative `CartesianBase` coordinate accessors over
//!      **all** cells of every location (parity — the FAQ introduces zero error),
//!   2. equal to the committed Julia-reference
//!      `tests/conformance/grids/arakawa/construction/golden.json`.
//!
//! The golden samples the four staggered location coordinate arrays, the
//! cell-center neighbor maps, and boundary masks at curated `points`, and carries
//! the full static A–E variable-location / shape table.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::arakawa::{BaseGrid, CartesianBase, Location, Stagger};
use earthsci_grids::grids::arakawa_faq::{arakawa_construction_faq, ArakawaConstruction};
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/arakawa")
}

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()))
}

fn close_rel(a: f64, b: f64, tol: f64) -> bool {
    let scale = 1.0f64.max(a.abs()).max(b.abs());
    (a - b).abs() <= tol * scale
}

fn parse_ints(s: &str) -> Vec<i64> {
    let v: Value = serde_json::from_str(s).expect("int string parses");
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_i64().unwrap())
        .collect()
}

const LOCATIONS: &[&str] = &["cell_center", "u_edge", "v_edge", "corner"];
const STAGGERS: &[&str] = &["A", "B", "C", "D", "E"];

fn loc_shape(loc: &str, nx: usize, ny: usize) -> (usize, usize) {
    match loc {
        "cell_center" => (nx, ny),
        "u_edge" => (nx + 1, ny),
        "v_edge" => (nx, ny + 1),
        "corner" => (nx + 1, ny + 1),
        other => panic!("unknown location {other}"),
    }
}

fn loc_name(loc: Location) -> &'static str {
    match loc {
        Location::CellCenter => "cell_center",
        Location::UEdge => "u_edge",
        Location::VEdge => "v_edge",
        Location::Corner => "corner",
    }
}

fn faq_coords<'a>(faq: &'a ArakawaConstruction, loc: &str) -> &'a (Vec<f64>, Vec<f64>) {
    match loc {
        "cell_center" => &faq.coords_cell_center,
        "u_edge" => &faq.coords_u_edge,
        "v_edge" => &faq.coords_v_edge,
        "corner" => &faq.coords_corner,
        other => panic!("unknown location {other}"),
    }
}

fn base_coord(base: &CartesianBase, loc: &str, i: usize, j: usize) -> (f64, f64) {
    match loc {
        "cell_center" => base.cell_center(i, j),
        "u_edge" => base.x_edge(i, j),
        "v_edge" => base.y_edge(i, j),
        "corner" => base.corner(i, j),
        other => panic!("unknown location {other}"),
    }
}

#[test]
fn arakawa_construction_faq_matches_golden_and_imperative() {
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
        let bo = &fixture["opts"]["base"];
        assert_eq!(bo["family"].as_str().unwrap(), "cartesian");
        let base = CartesianBase::new(
            bo["xlo"].as_f64().unwrap(),
            bo["xhi"].as_f64().unwrap(),
            bo["ylo"].as_f64().unwrap(),
            bo["yhi"].as_f64().unwrap(),
            bo["nx"].as_u64().unwrap() as usize,
            bo["ny"].as_u64().unwrap() as usize,
        )
        .unwrap();
        let gl = by_name
            .get(name)
            .unwrap_or_else(|| panic!("golden missing fixture {name}"));
        let nx = base.nx();
        let ny = base.ny();
        let faq = arakawa_construction_faq(&base, Stagger::A); // geometry stagger-independent

        // --- counts + scalars ---
        assert_eq!(faq.nx as u64, gl["nx"].as_u64().unwrap(), "{name}: nx");
        assert_eq!(faq.ny as u64, gl["ny"].as_u64().unwrap(), "{name}: ny");
        assert_eq!(
            (nx * ny) as u64,
            gl["n_cells"].as_u64().unwrap(),
            "{name}: n_cells"
        );
        assert!(
            close_rel(faq.dx, gl["dx"].as_f64().unwrap(), rel_tol),
            "{name}: dx"
        );
        assert!(
            close_rel(faq.dy, gl["dy"].as_f64().unwrap(), rel_tol),
            "{name}: dy"
        );
        assert!(
            close_rel(
                faq.cell_volume,
                gl["cell_volume"].as_f64().unwrap(),
                rel_tol
            ),
            "{name}: cell_volume"
        );

        // --- parity: FAQ == imperative base accessors, all cells/locations ---
        for loc in LOCATIONS {
            let (ni, nj) = loc_shape(loc, nx, ny);
            let (xs, ys) = faq_coords(&faq, loc);
            for j in 0..nj {
                for i in 0..ni {
                    let (ex, ey) = base_coord(&base, loc, i, j);
                    assert_eq!(xs[j * ni + i], ex, "{name}:{loc} x parity ({i},{j})");
                    assert_eq!(ys[j * ni + i], ey, "{name}:{loc} y parity ({i},{j})");
                }
            }
        }

        // --- cross-binding golden: coords at curated points ---
        for loc in LOCATIONS {
            let (ni, _) = loc_shape(loc, nx, ny);
            let pts = gl["coords"]["points"][*loc].as_array().unwrap();
            let gxy: Value =
                serde_json::from_str(gl["coords"][*loc].as_str().unwrap()).expect("coords parse");
            let (xs, ys) = faq_coords(&faq, loc);
            for (k, pt) in pts.iter().enumerate() {
                let i = pt[0].as_u64().unwrap() as usize;
                let j = pt[1].as_u64().unwrap() as usize;
                let flat = j * ni + i;
                let gx = gxy[k][0].as_f64().unwrap();
                let gy = gxy[k][1].as_f64().unwrap();
                assert!(close_rel(xs[flat], gx, rel_tol), "{name}:{loc} x golden");
                assert!(close_rel(ys[flat], gy, rel_tol), "{name}:{loc} y golden");
            }
        }

        // --- cross-binding golden: boundary masks at curated points ---
        let bd = &gl["boundary"];
        let bd_pts = bd["points"].as_array().unwrap();
        for (key, arr) in [
            ("x_lower", &faq.boundary_x_lower),
            ("x_upper", &faq.boundary_x_upper),
            ("y_lower", &faq.boundary_y_lower),
            ("y_upper", &faq.boundary_y_upper),
        ] {
            let exp = parse_ints(bd[key].as_str().unwrap());
            for (k, pt) in bd_pts.iter().enumerate() {
                let i = pt[0].as_u64().unwrap() as usize;
                let j = pt[1].as_u64().unwrap() as usize;
                assert_eq!(arr[j * nx + i] as i64, exp[k], "{name}: boundary {key}");
            }
        }

        // --- cross-binding golden: neighbor maps at curated points ---
        let ng = &gl["neighbor_indices"];
        let ng_pts = ng["points"].as_array().unwrap();
        for (key, arr) in [
            ("x_minus", &faq.neighbor_x_minus),
            ("x_plus", &faq.neighbor_x_plus),
            ("y_minus", &faq.neighbor_y_minus),
            ("y_plus", &faq.neighbor_y_plus),
        ] {
            let exp = parse_ints(ng[key].as_str().unwrap());
            for (k, pt) in ng_pts.iter().enumerate() {
                let i = pt[0].as_u64().unwrap() as usize;
                let j = pt[1].as_u64().unwrap() as usize;
                assert_eq!(arr[j * nx + i], exp[k], "{name}: neighbor {key}");
            }
        }

        // --- cross-binding golden: static A–E variable-location / shape table ---
        for s in STAGGERS {
            let stagger = Stagger::from_name(s).unwrap();
            let faq_s = arakawa_construction_faq(&base, stagger);
            let gs = &gl["staggers"][*s];
            assert_eq!(
                faq_s.rotated,
                gs["rotated"].as_bool().unwrap(),
                "{name}: {s} rotated"
            );
            let (h, u, v) = faq_s.variable_locations;
            assert_eq!(loc_name(h), gs["variable_locations"]["h"].as_str().unwrap());
            assert_eq!(loc_name(u), gs["variable_locations"]["u"].as_str().unwrap());
            assert_eq!(loc_name(v), gs["variable_locations"]["v"].as_str().unwrap());
            for loc in LOCATIONS {
                let (ni, nj) = loc_shape(loc, nx, ny);
                let exp = gs["location_shapes"][*loc].as_array().unwrap();
                assert_eq!(
                    ni as u64,
                    exp[0].as_u64().unwrap(),
                    "{name}: {s} loc_shape {loc}"
                );
                assert_eq!(
                    nj as u64,
                    exp[1].as_u64().unwrap(),
                    "{name}: {s} loc_shape {loc}"
                );
            }
            for (var, vloc) in [("h", h), ("u", u), ("v", v)] {
                let (ni, nj) = loc_shape(loc_name(vloc), nx, ny);
                let exp = gs["variable_shapes"][var].as_array().unwrap();
                assert_eq!(
                    ni as u64,
                    exp[0].as_u64().unwrap(),
                    "{name}: {s} var_shape {var}"
                );
                assert_eq!(
                    nj as u64,
                    exp[1].as_u64().unwrap(),
                    "{name}: {s} var_shape {var}"
                );
            }
        }
    }
}
