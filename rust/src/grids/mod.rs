//! Per-family grid accessor modules.
//!
//! Each family lands as its own submodule (e.g., `grids::cartesian`) per
//! `docs/GRIDS_API.md` §2.5. The crate root re-exports families via
//! `pub use grids::<family>` so the documented call form
//! `earthsci_grids::<family>::builder()` works at the crate level.

pub mod arakawa;
pub mod cartesian;
pub mod duo;
pub mod lat_lon;
pub mod mpas;
pub mod vertical;

// DUO front-door construction bridges (esd-un6 / W2): D3 subdivision, D2a primal
// geometry, and D1a value-invention topology routed through the ESS evaluator and
// relational primitives. Crate-internal — consumed by `duo::Builder::build`.
pub(crate) mod duo_primal_geometry_faq;
pub(crate) mod duo_subdivide_faq;
pub(crate) mod duo_topology_faq;
