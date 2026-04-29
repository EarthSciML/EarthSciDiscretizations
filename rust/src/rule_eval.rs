//! Thin passthrough to the canonical EarthSciSerialization scalar evaluator.
//!
//! ESD rule files (`discretizations/<family>/*.json`) carry stencil
//! coefficients as ESS `ExpressionNode` AST trees per `esm-spec.md` §4.2.
//! Per `AGENTS.md`, ESD does not carry a shadow evaluator: rule authoring
//! and evaluation flow through the ESS tree-walk evaluator. This module
//! is the single Rust entry point that imports
//! `earthsci_toolkit::evaluate`, so the rest of ESD imports one function
//! regardless of how the ESS public surface evolves (mirrors
//! `src/rule_eval.jl`).
//!
//! See `EarthSciSerialization/esm-spec.md` for the closed op registry.

use std::collections::HashMap;

use earthsci_toolkit::{evaluate, Expr};
use serde_json::Value;
use thiserror::Error;

/// Variable bindings consumed by [`eval_coeff`]. Keys are AST variable
/// names (e.g. `"R"`, `"cos_lat"`, `"dlon"`); values are the per-cell
/// floats supplied by the binding's grid accessor.
pub type Bindings = HashMap<String, f64>;

/// Errors raised while evaluating a rule-AST coefficient node.
#[derive(Debug, Error)]
pub enum RuleEvalError {
    /// The JSON does not deserialize into an ESS `Expr` (missing field,
    /// wrong shape, unrecognized literal kind).
    #[error("schema violation: {0}")]
    Schema(String),
    /// One or more variable references in the AST have no entry in the
    /// supplied `bindings`.
    #[error("unbound variables: {0:?}")]
    UnboundVariables(Vec<String>),
    /// Evaluation failed for a reason other than unbound variables (e.g.
    /// domain error, unsupported op for the scalar evaluator).
    #[error("evaluation failed")]
    EvalFailed,
}

/// Evaluate a JSON-decoded ESS expression-AST coefficient node by
/// delegating to `earthsci_toolkit::evaluate`.
///
/// `node` is what `serde_json` produces from a rule file's `"coeff"`,
/// `"stencil"` entry, or any nested `{"op": ..., "args": ...}` block:
/// a JSON number, a JSON string (variable reference), or an object
/// matching ESS's `ExpressionNode` schema.
pub fn eval_coeff(node: &Value, bindings: &Bindings) -> Result<f64, RuleEvalError> {
    let expr: Expr =
        serde_json::from_value(node.clone()).map_err(|e| RuleEvalError::Schema(e.to_string()))?;
    evaluate(&expr, bindings).map_err(|unbound| {
        if unbound.is_empty() {
            RuleEvalError::EvalFailed
        } else {
            RuleEvalError::UnboundVariables(unbound)
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn b(pairs: &[(&str, f64)]) -> Bindings {
        pairs.iter().map(|(k, v)| ((*k).to_string(), *v)).collect()
    }

    #[test]
    fn centered_2nd_uniform_latlon_first_coeff() {
        // -1 / (2 * R * cos_lat * dlon), lifted from
        // discretizations/finite_difference/centered_2nd_uniform_latlon.json.
        let r = 1.0;
        let cos_lat = 0.5_f64.cos();
        let dlon = std::f64::consts::PI / 8.0;
        let bindings = b(&[("R", r), ("cos_lat", cos_lat), ("dlon", dlon)]);
        let ast = json!({
            "op": "/",
            "args": [-1, {"op": "*", "args": [2, "R", "cos_lat", "dlon"]}]
        });
        let got = eval_coeff(&ast, &bindings).unwrap();
        let expected = -1.0 / (2.0 * r * cos_lat * dlon);
        assert!((got - expected).abs() <= 1e-15 * expected.abs().max(1.0));
    }

    #[test]
    fn unbound_variable_surfaces_name() {
        let bindings = Bindings::new();
        let err = eval_coeff(&json!("y"), &bindings).unwrap_err();
        match err {
            RuleEvalError::UnboundVariables(names) => assert_eq!(names, vec!["y".to_string()]),
            other => panic!("expected UnboundVariables, got {other:?}"),
        }
    }
}
