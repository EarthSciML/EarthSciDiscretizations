//! Per-family grid accessor modules.
//!
//! Each family lands as its own submodule (e.g., `grids::cartesian`) per
//! `docs/GRIDS_API.md` §2.5. The crate root re-exports families via
//! `pub use grids::<family>` so the documented call form
//! `earthsci_grids::<family>::builder()` works at the crate level.
