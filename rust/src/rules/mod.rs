//! Per-rule evaluator modules.
//!
//! Each rule that requires more than the scalar AST walk in
//! [`crate::rule_eval`] (sub-stencil dispatch, output-kind dispatch,
//! reconstruction post-processing, …) lands as its own submodule under
//! `rules::<rule_name>`. The crate root re-exports the canonical entry
//! points so binding-level call sites match `docs/GRIDS_API.md` §4.3.

pub mod ppm_reconstruction;
