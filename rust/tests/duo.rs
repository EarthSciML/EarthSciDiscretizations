//! Public-API integration test for the DUO family.
//!
//! Exercises the crate-root call form documented in `docs/GRIDS_API.md` §2.5:
//! `earthsci_grids::duo::builder().<opt>(…).build()`.

use earthsci_grids::duo::{self, DuoLoader};
use earthsci_grids::{Dtype, Grid};

#[test]
fn public_api_builder_shape_matches_spec() {
    let grid = duo::builder()
        .loader(DuoLoader::builtin_level(1))
        .r(6.371e6)
        .dtype(Dtype::F64)
        .ghosts(0)
        .build()
        .expect("level-1 icosahedron builds");

    assert_eq!(grid.family(), "duo");
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.level(), 1);
    assert_eq!(grid.n_cells(), 80);
    assert_eq!(grid.n_vertices(), 42);
    assert_eq!(grid.n_edges(), 120);

    // Every grid implements the Grid trait and lowers to a §6-schema-valid
    // declarative config.
    let esm = grid.to_esm();
    assert_eq!(esm["family"], "duo");
    assert_eq!(esm["dtype"], "float64");
    assert_eq!(esm["options"]["level"], 1);
    assert_eq!(esm["provenance"]["binding"], "rust");
}

#[test]
fn defaults_match_spec() {
    // R defaults to 6.371e6, dtype defaults to F64, ghosts to 0.
    let grid = duo::builder()
        .loader(DuoLoader::builtin_level(0))
        .build()
        .unwrap();
    assert_eq!(grid.radius(), 6.371e6);
    assert_eq!(grid.dtype(), Dtype::F64);
    assert_eq!(grid.ghosts(), 0);
}
