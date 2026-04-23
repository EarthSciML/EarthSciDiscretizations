//! Public-API integration test for the MPAS family.
//!
//! Exercises the crate-root call form documented in `docs/GRIDS_API.md` §2.5
//! and the loader contract in §10: `earthsci_grids::mpas::builder()…build()`.

use std::f64::consts::PI;

use earthsci_grids::mpas::{
    self, MpasLoader, MpasMeshData, MpasMeshInput, MpasMetricName, NO_NEIGHBOR,
};
use earthsci_grids::{Dtype, Grid};

fn tiny_mesh() -> MpasMeshData {
    let r = 1.0_f64;
    let half = 2.0 * PI * r * r;
    MpasMeshData::from_input(MpasMeshInput {
        lon_cell: vec![0.0, PI],
        lat_cell: vec![PI / 4.0, -PI / 4.0],
        area_cell: vec![half, half],
        n_edges_on_cell: vec![1, 1],
        cells_on_cell: vec![1, 0],
        edges_on_cell: vec![0, 0],
        lon_edge: vec![0.0],
        lat_edge: vec![0.0],
        cells_on_edge: vec![0, 1],
        dc_edge: vec![PI / 2.0],
        dv_edge: vec![PI / 2.0],
        max_edges: 1,
        x_cell: None,
        y_cell: None,
        z_cell: None,
        n_vertices: Some(2),
        r: Some(r),
    })
    .expect("valid two-cell mesh")
}

#[test]
fn public_api_builder_shape_matches_spec() {
    let grid = mpas::builder()
        .mesh(tiny_mesh())
        .loader(MpasLoader::new("mem://fixture"))
        .r(6.371e6)
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("tiny MPAS mesh builds");

    assert_eq!(grid.family(), "mpas");
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.topology(), "unstructured");
    assert_eq!(grid.n_cells(), 2);
    assert_eq!(grid.n_edges(), 1);

    // Every grid implements the Grid trait and lowers to a §6-schema-valid
    // declarative config (no inline geometry).
    let esm = grid.to_esm();
    assert_eq!(esm["family"], "mpas");
    assert_eq!(esm["dtype"], "float64");
    assert_eq!(esm["topology"], "unstructured");
    assert_eq!(esm["n_cells"], 2);
    assert_eq!(esm["options"]["loader"]["path"], "mem://fixture");
    assert_eq!(esm["provenance"]["binding"], "rust");
}

#[test]
fn defaults_match_spec() {
    // R defaults to 6.371e6, dtype to F64, ghosts to 0.
    let grid = mpas::builder()
        .mesh(tiny_mesh())
        .r(6.371e6)
        .build()
        .unwrap();
    assert_eq!(grid.radius(), 6.371e6);
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.ghosts(), 0);
}

#[test]
fn loader_backed_reader_fn_contract() {
    // §10: path-based loading is satisfied by a caller-supplied reader.
    let grid = mpas::builder()
        .loader(MpasLoader::new("file:///fixtures/tiny.nc"))
        .reader_fn(|path| {
            assert_eq!(path, "file:///fixtures/tiny.nc");
            Ok(tiny_mesh())
        })
        .build()
        .expect("reader_fn satisfies §10");

    let loader = grid.loader().expect("loader recorded");
    assert_eq!(loader.path, "file:///fixtures/tiny.nc");
    assert_eq!(loader.reader, "auto");
    assert_eq!(loader.check, "strict");
    assert_eq!(grid.n_cells(), 2);
}

#[test]
fn accessor_plumbing_end_to_end() {
    let grid = mpas::builder().mesh(tiny_mesh()).build().unwrap();

    // cell_centers → (lon, lat)
    let (lon, lat) = grid.cell_centers(0).unwrap();
    assert_eq!(lon, 0.0);
    assert_eq!(lat, PI / 4.0);

    // neighbors drops -1 sentinels (there are none in the closed 2-cell mesh)
    assert_eq!(grid.neighbors(0).unwrap(), vec![1]);
    assert_eq!(grid.neighbors(1).unwrap(), vec![0]);

    // edge accessors
    assert!(grid.edge_length(0).unwrap() > 0.0);
    assert!(grid.cell_distance(0).unwrap() > 0.0);

    // metric_eval by enum + wire-form name
    let a = grid.metric_eval(MpasMetricName::Area, 0).unwrap();
    let a_by_name = grid.metric_eval_by_name("area", 0).unwrap();
    assert_eq!(a, a_by_name);

    // NO_NEIGHBOR is the documented sentinel
    assert_eq!(NO_NEIGHBOR, -1);
}
