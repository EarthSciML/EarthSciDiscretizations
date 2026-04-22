//! Public-API integration test for the lat-lon family.
//!
//! Exercises the crate-root call form documented in `docs/GRIDS_API.md` §2.5:
//! `earthsci_grids::lat_lon::builder().<opt>(…).build()`.

use earthsci_grids::lat_lon::{self, LatLonVariant, PolePolicy};
use earthsci_grids::{Dtype, Grid};

#[test]
fn public_api_builder_shape_matches_spec() {
    let grid = lat_lon::builder()
        .nlon(72)
        .nlat(36)
        .r(6.371e6)
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("regular lat-lon builds");

    assert_eq!(grid.family(), "lat_lon");
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.variant(), LatLonVariant::Regular);
    assert_eq!(grid.n_cells(), 72 * 36);
    assert_eq!(grid.nlon_uniform(), Some(72));
    assert_eq!(grid.nlat(), 36);
    assert_eq!(grid.pole_policy(), PolePolicy::None);

    // Every grid implements the Grid trait and lowers to a §6-schema-valid
    // declarative config.
    let esm = grid.to_esm();
    assert_eq!(esm["family"], "lat_lon");
    assert_eq!(esm["dtype"], "float64");
    assert_eq!(esm["variant"], "regular");
    assert_eq!(esm["params"]["nlon"], 72);
    assert_eq!(esm["params"]["nlat"], 36);
    assert_eq!(esm["provenance"]["binding"], "rust");
}

#[test]
fn reduced_gaussian_public_api() {
    let grid = lat_lon::builder()
        .variant(LatLonVariant::ReducedGaussian)
        .nlon_per_row(vec![4, 8, 12, 12, 8, 4])
        .r(6.371e6)
        .build()
        .expect("reduced-Gaussian builds with schedule only");
    assert_eq!(grid.variant(), LatLonVariant::ReducedGaussian);
    assert_eq!(grid.nlat(), 6);
    assert_eq!(grid.n_cells(), 4 + 8 + 12 + 12 + 8 + 4);
    assert!(grid.nlon_uniform().is_none());

    let esm = grid.to_esm();
    assert_eq!(esm["variant"], "reduced_gaussian");
    assert_eq!(esm["params"]["nlat"], 6);
    assert!(esm["params"]["nlon_per_row"].is_array());
}

#[test]
fn defaults_match_spec() {
    let grid = lat_lon::builder().nlon(4).nlat(4).build().unwrap();
    assert_eq!(grid.r(), 6.371e6);
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.ghosts(), 0);
    assert_eq!(grid.pole_policy(), PolePolicy::None);
}

#[test]
fn missing_required_nlon_errors() {
    let err = lat_lon::builder().nlat(4).build().unwrap_err();
    assert!(matches!(
        err,
        earthsci_grids::GridError::MissingOption("nlon")
    ));
}
