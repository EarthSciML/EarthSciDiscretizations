//! Integration smoke tests for the cubed_sphere grid family (public API).
//!
//! Verifies the call form documented in `docs/GRIDS_API.md` §2.5 works from
//! outside the crate and produces a grid with the §7 common-minimum fields.

use earthsci_grids::{cubed_sphere, Dtype, Grid};

#[test]
fn public_builder_call_form() {
    let g = cubed_sphere::builder()
        .nc(8)
        .r(6.371e6)
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("builder should succeed with required Nc");

    assert_eq!(g.family(), "cubed_sphere");
    assert_eq!(g.dtype(), Dtype::F64);
    assert_eq!(g.n_cells(), 6 * 8 * 8);
}

#[test]
fn to_esm_schema_fields_present() {
    let g = cubed_sphere::builder().nc(4).build().unwrap();
    let doc = g.to_esm();
    for key in [
        "family",
        "version",
        "dtype",
        "topology",
        "generator",
        "params",
        "provenance",
    ] {
        assert!(doc.get(key).is_some(), "missing key: {key}");
    }
}

#[test]
fn missing_required_nc_errors() {
    let err = cubed_sphere::builder().build().unwrap_err();
    assert!(matches!(
        err,
        earthsci_grids::GridError::MissingOption("Nc")
    ));
}
