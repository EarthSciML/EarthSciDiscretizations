//! Cross-language conformance test for the arakawa grid family.
//!
//! Verifies that the Rust accessor runtime produces accessor outputs matching
//! the committed Julia-reference golden within the tolerance declared in
//! `tests/conformance/grids/arakawa/fixtures.json`.

use std::fs;
use std::path::{Path, PathBuf};

use earthsci_grids::arakawa::{self, ArakawaGrid, CartesianBase, Location, Stagger, Variable};
use earthsci_grids::Dtype;
use serde_json::Value;

fn harness_dir() -> PathBuf {
    let base: PathBuf = env!("CARGO_MANIFEST_DIR").into();
    base.parent()
        .expect("rust crate has a parent directory")
        .join("tests/conformance/grids/arakawa")
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

fn location_name(loc: Location) -> &'static str {
    match loc {
        Location::CellCenter => "cell_center",
        Location::UEdge => "u_edge",
        Location::VEdge => "v_edge",
        Location::Corner => "corner",
    }
}

fn stagger_from_name(name: &str) -> Stagger {
    Stagger::from_name(name).unwrap_or_else(|| panic!("unknown stagger {name}"))
}

fn build_grid(base_opts: &Value, stagger: Stagger, dtype: Dtype, ghosts: usize) -> ArakawaGrid {
    assert_eq!(
        base_opts["family"], "cartesian",
        "only cartesian base is supported"
    );
    let base = CartesianBase::new(
        base_opts["xlo"].as_f64().unwrap(),
        base_opts["xhi"].as_f64().unwrap(),
        base_opts["ylo"].as_f64().unwrap(),
        base_opts["yhi"].as_f64().unwrap(),
        base_opts["nx"].as_u64().unwrap() as usize,
        base_opts["ny"].as_u64().unwrap() as usize,
    )
    .expect("valid cartesian base");
    arakawa::builder()
        .base(base)
        .stagger(stagger)
        .dtype(dtype)
        .ghosts(ghosts)
        .build()
        .expect("arakawa grid builds")
}

fn coord_at(grid_a: &ArakawaGrid, grid_c: &ArakawaGrid, loc: Location, i: usize, j: usize)
    -> (f64, f64)
{
    match loc {
        Location::CellCenter => grid_a.cell_centers(i, j).unwrap(),
        Location::Corner => grid_a.corners(i, j).unwrap(),
        Location::UEdge => grid_c.u_face(i, j).unwrap(),
        Location::VEdge => grid_c.v_face(i, j).unwrap(),
    }
}

fn neighbor_to_wire(n: Option<(usize, usize)>) -> Option<[usize; 2]> {
    n.map(|(i, j)| [i, j])
}

fn golden_neighbor(v: &Value) -> Option<[usize; 2]> {
    if v.is_null() {
        None
    } else {
        let arr = v.as_array().unwrap();
        Some([
            arr[0].as_u64().unwrap() as usize,
            arr[1].as_u64().unwrap() as usize,
        ])
    }
}

#[test]
fn arakawa_matches_golden() {
    let hdir = harness_dir();
    let spec = read_json(&hdir.join("fixtures.json"));
    let rel_tol = spec["tolerance"]["relative"].as_f64().unwrap();

    let dir_keys = ["W", "E", "S", "N"];
    let locations = [
        Location::CellCenter,
        Location::UEdge,
        Location::VEdge,
        Location::Corner,
    ];

    for fixture in spec["fixtures"].as_array().unwrap() {
        let name = fixture["name"].as_str().unwrap();
        let opts = &fixture["opts"];
        let base_opts = &opts["base"];
        let dtype_str = opts["dtype"].as_str().unwrap();
        let dtype = match dtype_str {
            "float64" => Dtype::F64,
            "float32" => Dtype::F32,
            other => panic!("unknown dtype {other}"),
        };
        let ghosts = opts["ghosts"].as_u64().unwrap() as usize;

        let golden = read_json(&hdir.join(format!("golden/{name}.json")));

        let g_a = build_grid(base_opts, Stagger::A, dtype, ghosts);
        let g_c = build_grid(base_opts, Stagger::C, dtype, ghosts);

        let expected_cells = golden["n_cells"].as_u64().unwrap() as usize;
        assert_eq!(g_a.n_cells(), expected_cells, "{name}: n_cells");

        // ---- coord + neighbor tables (stagger-independent) -------------
        for loc in locations {
            let lname = location_name(loc);
            let coords = &golden["coords"][lname];
            let nbrs = &golden["neighbors"][lname];
            let points = coords["points"].as_array().unwrap();

            for (k, qp) in points.iter().enumerate() {
                let i = qp[0].as_u64().unwrap() as usize;
                let j = qp[1].as_u64().unwrap() as usize;
                let (x, y) = coord_at(&g_a, &g_c, loc, i, j);
                let ex_xy = &coords["xy"][k];
                let ex_x = ex_xy[0].as_f64().unwrap();
                let ex_y = ex_xy[1].as_f64().unwrap();
                assert!(
                    close_rel(x, ex_x, rel_tol),
                    "{name}: {lname} x mismatch at qp[{k}]=({i},{j}): {x} vs {ex_x}"
                );
                assert!(
                    close_rel(y, ex_y, rel_tol),
                    "{name}: {lname} y mismatch at qp[{k}]=({i},{j}): {y} vs {ex_y}"
                );

                let [w, e, s, n] = g_a.neighbors(loc, i, j).unwrap();
                let got = [w, e, s, n].map(neighbor_to_wire);
                for (idx, dkey) in dir_keys.iter().enumerate() {
                    let expected = golden_neighbor(&nbrs[dkey][k]);
                    assert_eq!(
                        got[idx], expected,
                        "{name}: neighbor {dkey} {lname} mismatch at qp[{k}]=({i},{j})"
                    );
                }
            }
        }

        // ---- metric_eval (stagger-independent for cartesian) ----------
        let metrics = &golden["metrics"];
        let mpoints = metrics["points"].as_array().unwrap();
        for (k, qp) in mpoints.iter().enumerate() {
            let i = qp[0].as_u64().unwrap() as usize;
            let j = qp[1].as_u64().unwrap() as usize;
            for mname in ["dx", "dy", "area"] {
                let got = g_a.metric_eval_by_name(mname, i, j).unwrap();
                let expected = metrics[mname][k].as_f64().unwrap();
                assert!(
                    close_rel(got, expected, rel_tol),
                    "{name}: metric {mname} mismatch at qp[{k}]=({i},{j}): {got} vs {expected}"
                );
            }
        }

        // ---- per-stagger variable locations + shapes ------------------
        for sname_v in fixture["staggers"].as_array().unwrap() {
            let sname = sname_v.as_str().unwrap();
            let stagger = stagger_from_name(sname);
            let g_s = build_grid(base_opts, stagger, dtype, ghosts);
            let stab = &golden["staggers"][sname];

            let expected_rotated = stab["rotated"].as_bool().unwrap();
            assert_eq!(
                stagger == Stagger::E,
                expected_rotated,
                "{name}: stagger={sname} rotated flag"
            );

            for (var, vname) in [(Variable::H, "h"), (Variable::U, "u"), (Variable::V, "v")] {
                let got_loc = location_name(g_s.variable_location(var));
                let ex_loc = stab["variable_locations"][vname].as_str().unwrap();
                assert_eq!(
                    got_loc, ex_loc,
                    "{name}: stagger={sname} var={vname} location"
                );
                let (ni, nj) = g_s.variable_shape(var);
                let ex = stab["variable_shapes"][vname].as_array().unwrap();
                let ex_ni = ex[0].as_u64().unwrap() as usize;
                let ex_nj = ex[1].as_u64().unwrap() as usize;
                assert_eq!(
                    (ni, nj),
                    (ex_ni, ex_nj),
                    "{name}: stagger={sname} var={vname} shape"
                );
            }

            for loc in locations {
                let lname = location_name(loc);
                let (ni, nj) = g_s.location_shape(loc);
                let ex = stab["location_shapes"][lname].as_array().unwrap();
                let ex_ni = ex[0].as_u64().unwrap() as usize;
                let ex_nj = ex[1].as_u64().unwrap() as usize;
                assert_eq!(
                    (ni, nj),
                    (ex_ni, ex_nj),
                    "{name}: stagger={sname} location={lname} shape"
                );
            }
        }
    }
}
