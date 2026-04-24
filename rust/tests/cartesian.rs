//! Public-API integration test for the cartesian family.
//!
//! Exercises the crate-root call form documented in `docs/GRIDS_API.md` §2.5:
//! `earthsci_grids::cartesian::builder().<opt>(…).build()`.

use earthsci_grids::cartesian::{self, CartesianNeighbor, MetricName};
use earthsci_grids::{Dtype, Grid};

#[test]
fn public_api_builder_shape_matches_spec() {
    let grid = cartesian::builder()
        .nx(8)
        .ny(4)
        .extent(vec![(0.0, 1.0), (0.0, 0.5)])
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("uniform 2D cartesian builds");

    assert_eq!(grid.family(), "cartesian");
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.ndim(), 2);
    assert_eq!(grid.n(), &[8, 4]);
    assert_eq!(grid.n_cells(), 8 * 4);
    assert_eq!(grid.topology(), "rectilinear");
    assert_eq!(grid.uniform(), &[true, true]);

    // Every grid implements the Grid trait and lowers to a §6-schema-valid
    // declarative config.
    let esm = grid.to_esm();
    assert_eq!(esm["family"], "cartesian");
    assert_eq!(esm["dtype"], "float64");
    assert_eq!(esm["ndim"], 2);
    assert_eq!(esm["n_cells"], 8 * 4);
    assert_eq!(esm["provenance"]["binding"], "rust");
}

#[test]
fn public_api_1d_uniform() {
    let grid = cartesian::builder()
        .nx(4)
        .extent(vec![(0.0, 1.0)])
        .build()
        .expect("1D uniform builds");
    assert_eq!(grid.ndim(), 1);
    assert_eq!(grid.n_cells(), 4);
    let c = grid.cell_center(&[2]).unwrap();
    assert!((c[0] - 0.625).abs() < 1e-15);
}

#[test]
fn public_api_3d_uniform() {
    let grid = cartesian::builder()
        .nx(4)
        .ny(2)
        .nz(5)
        .extent(vec![(0.0, 1.0), (0.0, 0.5), (0.0, 0.25)])
        .build()
        .expect("3D uniform builds");
    assert_eq!(grid.ndim(), 3);
    assert_eq!(grid.n_cells(), 4 * 2 * 5);
    let v = grid.cell_volume(&[0, 0, 0]).unwrap();
    // dx=0.25, dy=0.25, dz=0.05 -> 0.003125
    assert!((v - 0.25 * 0.25 * 0.05).abs() < 1e-15);
}

#[test]
fn public_api_nonuniform_edges() {
    let grid = cartesian::builder()
        .edges(vec![vec![0.0, 0.1, 0.5, 1.0]])
        .build()
        .expect("non-uniform 1D builds from edges");
    assert_eq!(grid.ndim(), 1);
    assert_eq!(grid.n(), &[3]);
    assert_eq!(grid.uniform(), &[false]);
    let esm = grid.to_esm();
    assert!(esm.get("edges").is_some());
}

#[test]
fn public_api_neighbors_at_corner() {
    let grid = cartesian::builder()
        .nx(4)
        .ny(2)
        .extent(vec![(0.0, 1.0), (0.0, 0.5)])
        .build()
        .unwrap();
    let ns = grid.neighbors(&[0, 0]).unwrap();
    // corner of 2D grid has exactly two neighbors, both on the high side.
    assert_eq!(ns.len(), 2);
    assert!(ns.iter().all(|n: &CartesianNeighbor| n.side == 1));
}

#[test]
fn public_api_metric_volume_equals_cell_volume() {
    let grid = cartesian::builder()
        .nx(4)
        .ny(2)
        .extent(vec![(0.0, 1.0), (0.0, 0.5)])
        .build()
        .unwrap();
    let m = grid.metric_eval(MetricName::Volume, &[1, 0]).unwrap();
    let v = grid.cell_volume(&[1, 0]).unwrap();
    assert!((m - v).abs() < 1e-15);
}

#[test]
fn defaults_match_spec() {
    let grid = cartesian::builder()
        .nx(4)
        .extent(vec![(0.0, 1.0)])
        .build()
        .unwrap();
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.ghosts(), 0);
}

#[test]
fn missing_required_nx_errors() {
    let err = cartesian::builder()
        .extent(vec![(0.0, 1.0)])
        .build()
        .unwrap_err();
    assert!(matches!(
        err,
        earthsci_grids::GridError::MissingOption("nx")
    ));
}

#[test]
fn missing_required_extent_errors() {
    let err = cartesian::builder().nx(4).build().unwrap_err();
    assert!(matches!(
        err,
        earthsci_grids::GridError::MissingOption("extent")
    ));
}

#[test]
fn edges_and_uniform_mutually_exclusive() {
    let err = cartesian::builder()
        .nx(4)
        .edges(vec![vec![0.0, 1.0]])
        .build()
        .unwrap_err();
    assert!(matches!(
        err,
        earthsci_grids::GridError::InvalidOption("edges", _)
    ));
}
