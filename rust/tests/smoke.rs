//! Smoke test: crate compiles, public types are reachable.

use earthsci_grids::{Dtype, GridError};

#[test]
fn dtype_default_is_f64() {
    assert_eq!(Dtype::default(), Dtype::F64);
}

#[test]
fn grid_error_displays() {
    let e = GridError::MissingOption("Nc");
    assert!(format!("{e}").contains("Nc"));
}
