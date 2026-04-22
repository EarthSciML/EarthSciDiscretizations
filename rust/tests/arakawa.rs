//! Public-API integration test for the arakawa family.
//!
//! Exercises the crate-root call form documented in `docs/GRIDS_API.md` §2.5:
//! `earthsci_grids::arakawa::builder().<opt>(…).build()`.

use earthsci_grids::arakawa::{self, CartesianBase, Location, Stagger, Variable};
use earthsci_grids::{Dtype, Grid, GridError};

#[test]
fn public_builder_call_form() {
    let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 8, 8).unwrap();
    let grid = arakawa::builder()
        .base(base)
        .stagger(Stagger::C)
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("arakawa C-grid builds");

    assert_eq!(grid.family(), "arakawa");
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.n_cells(), 64);
    assert_eq!(grid.stagger(), Stagger::C);
    assert_eq!(grid.variable_shape(Variable::U), (9, 8));
    assert_eq!(grid.variable_shape(Variable::V), (8, 9));
}

#[test]
fn to_esm_schema_fields_present() {
    let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
    let grid = arakawa::builder()
        .base(base)
        .stagger(Stagger::B)
        .build()
        .unwrap();
    let doc = grid.to_esm();
    for key in [
        "family",
        "version",
        "dtype",
        "topology",
        "ghosts",
        "n_cells",
        "stagger",
        "rotated",
        "base",
        "provenance",
    ] {
        assert!(doc.get(key).is_some(), "missing key: {key}");
    }
    assert_eq!(doc["family"], "arakawa");
    assert_eq!(doc["stagger"], "B");
    assert_eq!(doc["rotated"], false);
}

#[test]
fn missing_required_base_errors() {
    let err = arakawa::builder().stagger(Stagger::A).build().unwrap_err();
    assert!(matches!(err, GridError::MissingOption("base")));
}

#[test]
fn missing_required_stagger_errors() {
    let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
    let err = arakawa::builder().base(base).build().unwrap_err();
    assert!(matches!(err, GridError::MissingOption("stagger")));
}

#[test]
fn metric_eval_by_name_works_across_the_public_api() {
    let base = CartesianBase::new(0.0, 2.0, 0.0, 6.0, 4, 3).unwrap();
    let grid = arakawa::builder()
        .base(base)
        .stagger(Stagger::C)
        .build()
        .unwrap();
    assert!((grid.metric_eval_by_name("dx", 0, 0).unwrap() - 0.5).abs() < 1e-14);
    assert!((grid.metric_eval_by_name("dy", 0, 0).unwrap() - 2.0).abs() < 1e-14);
    assert!((grid.metric_eval_by_name("area", 3, 2).unwrap() - 1.0).abs() < 1e-14);
}

#[test]
fn uedge_neighbors_include_cross_column_stepping() {
    let base = CartesianBase::new(0.0, 1.0, 0.0, 1.0, 4, 4).unwrap();
    let grid = arakawa::builder()
        .base(base)
        .stagger(Stagger::C)
        .build()
        .unwrap();
    // UEdge shape is (5, 4). Interior UEdge cell has all four neighbors.
    let [w, e, s, n] = grid.neighbors(Location::UEdge, 2, 2).unwrap();
    assert_eq!(w, Some((1, 2)));
    assert_eq!(e, Some((3, 2)));
    assert_eq!(s, Some((2, 1)));
    assert_eq!(n, Some((2, 3)));
}
