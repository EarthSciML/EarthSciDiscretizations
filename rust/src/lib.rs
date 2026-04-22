//! earthsci_grids — Grid accessor runtime for the EarthSci model stack.
//!
//! Conforms to the cross-binding contract in `docs/GRIDS_API.md`.
//! Per-family grid accessor modules land under [`grids`] as Phase 3
//! implementation beads complete.

pub mod grids;

pub use error::{GridError, Result};
pub use grids::duo;
pub use traits::{Dtype, Grid};

// Re-export per-family modules at the crate root so the documented
// `earthsci_grids::<family>::builder()` call form works.
pub use grids::arakawa;
pub use grids::cubed_sphere;
pub use grids::lat_lon;

mod error {
    use thiserror::Error;

    #[derive(Debug, Error)]
    pub enum GridError {
        #[error("missing required option: {0}")]
        MissingOption(&'static str),
        #[error("invalid option {0}: {1}")]
        InvalidOption(&'static str, String),
        #[error("schema violation: {0}")]
        SchemaViolation(String),
        #[error("loader I/O error: {0}")]
        Io(#[from] std::io::Error),
        #[error("conformance tolerance exceeded: {0}")]
        Conformance(String),
    }

    pub type Result<T> = std::result::Result<T, GridError>;
}

mod traits {
    #[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
    pub enum Dtype {
        #[default]
        F64,
        F32,
    }

    /// Every family's returned grid struct implements this trait.
    /// See `docs/GRIDS_API.md` §3.3.
    pub trait Grid {
        fn family(&self) -> &'static str;
        fn dtype(&self) -> Dtype;
        fn to_esm(&self) -> serde_json::Value;
    }
}
