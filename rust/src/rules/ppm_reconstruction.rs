//! Rust evaluator for the `ppm_reconstruction` declarative rule.
//!
//! Mirrors the Julia reference (CW84 eqs. 1.6-1.10) by consuming the
//! ESS-emitted AST coefficients directly from
//! `discretizations/finite_volume/ppm_reconstruction.json`. Two layers
//! of dispatch are exposed, matching the rule's structure:
//!
//! * **`SubStencil`** picks among the named stencil entries
//!   (`q_left_edge`, `q_right_edge`). Each is a 4-point edge
//!   interpolation — coefficients live in the rule AST and are walked
//!   by [`crate::rule_eval::eval_coeff`].
//! * **`OutputKind`** picks the produced quantity: a left-edge value, a
//!   right-edge value, or the parabolic reconstruction sampled at a
//!   sub-cell coordinate `ξ ∈ [0, 1]`. The parabola formula
//!   `a(ξ) = a_L + ξ·(δa + a₆·(1-ξ))` (CW84 eqs. 1.5, 1.7, 1.10) is
//!   pinned by the rule's `parabola` documentation block.
//!
//! Coefficient ASTs in the canonical rule are pure rationals
//! (`{"op": "/", "args": [-1, 12]}` etc.), so evaluation needs no
//! per-cell bindings. The neighbour-window indexer abstraction lets
//! callers apply the stencil over any 1D slice (periodic, ghosted, …)
//! without baking a boundary policy into the evaluator.

use std::collections::HashMap;

use serde_json::Value;
use thiserror::Error;

use crate::rule_eval::{eval_coeff, RuleEvalError};

/// Sub-stencil selector — picks one of the rule's named stencil
/// entries. The PPM rule defines `q_left_edge` (interface i-1/2) and
/// `q_right_edge` (interface i+1/2).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SubStencil {
    QLeftEdge,
    QRightEdge,
}

/// Output-kind selector — picks the quantity the evaluator produces.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum OutputKind {
    /// Left-edge value `q_{i-1/2}` — applies `q_left_edge` sub-stencil.
    QLeftEdge,
    /// Right-edge value `q_{i+1/2}` — applies `q_right_edge`.
    QRightEdge,
    /// Parabolic reconstruction `a(ξ)` sampled at `xi ∈ [0, 1]`.
    QParabola { xi: f64 },
}

/// Errors raised while parsing or evaluating the PPM rule.
#[derive(Debug, Error)]
pub enum PpmRuleError {
    #[error("schema violation: {0}")]
    Schema(String),
    #[error("missing sub-stencil `{0}`")]
    MissingSubStencil(&'static str),
    #[error("coefficient AST evaluation failed: {0}")]
    Coeff(#[from] RuleEvalError),
    #[error("xi out of [0, 1]: {0}")]
    ParabolaArg(f64),
}

/// One entry in a PPM sub-stencil: a relative offset along the
/// reconstruction axis and the AST coefficient applied to the
/// neighbour cell average at that offset.
#[derive(Debug, Clone)]
pub struct StencilEntry {
    pub offset: i64,
    pub coeff: f64,
}

/// Parsed `ppm_reconstruction` rule.
///
/// The two named sub-stencils are resolved into pre-evaluated
/// `(offset, coeff)` pairs at construction time so the hot path is
/// just a multiply-add over the neighbour window.
#[derive(Debug, Clone)]
pub struct PpmRule {
    q_left_edge: Vec<StencilEntry>,
    q_right_edge: Vec<StencilEntry>,
}

impl PpmRule {
    /// Parse the rule from a JSON-decoded
    /// `discretizations/finite_volume/ppm_reconstruction.json`. Accepts
    /// either the full file shape (`{"discretizations":
    /// {"ppm_reconstruction": {...}}}`) or the inner spec object.
    pub fn from_json(value: &Value) -> Result<Self, PpmRuleError> {
        let spec = if let Some(discs) = value.get("discretizations") {
            discs.get("ppm_reconstruction").ok_or_else(|| {
                PpmRuleError::Schema(
                    "discretizations object missing `ppm_reconstruction` entry".into(),
                )
            })?
        } else {
            value
        };
        let stencil = spec
            .get("stencil")
            .and_then(Value::as_object)
            .ok_or_else(|| PpmRuleError::Schema("rule missing `stencil` object".into()))?;

        let q_left_edge = parse_sub_stencil(stencil, "q_left_edge")?;
        let q_right_edge = parse_sub_stencil(stencil, "q_right_edge")?;
        Ok(Self {
            q_left_edge,
            q_right_edge,
        })
    }

    /// Return the parsed entries for a sub-stencil.
    pub fn entries(&self, sub: SubStencil) -> &[StencilEntry] {
        match sub {
            SubStencil::QLeftEdge => &self.q_left_edge,
            SubStencil::QRightEdge => &self.q_right_edge,
        }
    }

    /// Apply a sub-stencil to a neighbour-window indexer.
    ///
    /// `q_at` receives the relative offset (e.g. `-2..=1`) and returns
    /// the cell-averaged value at that neighbour. The caller owns the
    /// boundary policy (periodic wrap, ghost padding, …).
    pub fn apply_sub_stencil<F>(&self, sub: SubStencil, q_at: F) -> f64
    where
        F: Fn(i64) -> f64,
    {
        self.entries(sub)
            .iter()
            .map(|e| e.coeff * q_at(e.offset))
            .sum()
    }

    /// Top-level dispatch on `output_kind`.
    ///
    /// Edge outputs apply the matching sub-stencil; the parabola output
    /// applies both edge sub-stencils to recover `a_L` and `a_R`,
    /// reads `q_i` from `q_at(0)`, and samples the cell parabola at
    /// `xi`.
    pub fn evaluate<F>(&self, output: OutputKind, q_at: F) -> Result<f64, PpmRuleError>
    where
        F: Fn(i64) -> f64,
    {
        match output {
            OutputKind::QLeftEdge => Ok(self.apply_sub_stencil(SubStencil::QLeftEdge, &q_at)),
            OutputKind::QRightEdge => Ok(self.apply_sub_stencil(SubStencil::QRightEdge, &q_at)),
            OutputKind::QParabola { xi } => {
                if !(0.0..=1.0).contains(&xi) {
                    return Err(PpmRuleError::ParabolaArg(xi));
                }
                let a_l = self.apply_sub_stencil(SubStencil::QLeftEdge, &q_at);
                let a_r = self.apply_sub_stencil(SubStencil::QRightEdge, &q_at);
                let q_i = q_at(0);
                Ok(parabola(a_l, a_r, q_i, xi))
            }
        }
    }
}

/// Parabolic reconstruction `a(ξ)` from CW84 eqs. (1.5), (1.7), (1.10).
///
/// `xi` is the normalised cell coordinate in `[0, 1]` (0 = left edge,
/// 1 = right edge); `a_L`, `a_R` are the limited or unlimited edge
/// values; `q_i` is the cell average. Cell-average preservation:
/// `(a_L + a_R)/2 + a_6/6 = q_i`.
pub fn parabola(a_l: f64, a_r: f64, q_i: f64, xi: f64) -> f64 {
    let da = a_r - a_l;
    let a6 = 6.0 * (q_i - 0.5 * (a_l + a_r));
    a_l + xi * (da + a6 * (1.0 - xi))
}

fn parse_sub_stencil(
    stencil: &serde_json::Map<String, Value>,
    key: &'static str,
) -> Result<Vec<StencilEntry>, PpmRuleError> {
    let arr = stencil
        .get(key)
        .and_then(Value::as_array)
        .ok_or(PpmRuleError::MissingSubStencil(key))?;
    let bindings = HashMap::new();
    let mut out = Vec::with_capacity(arr.len());
    for (idx, entry) in arr.iter().enumerate() {
        let offset = entry
            .get("selector")
            .and_then(|s| s.get("offset"))
            .and_then(Value::as_i64)
            .ok_or_else(|| {
                PpmRuleError::Schema(format!(
                    "{key}[{idx}]: selector.offset missing or not an integer"
                ))
            })?;
        let coeff_node = entry
            .get("coeff")
            .ok_or_else(|| PpmRuleError::Schema(format!("{key}[{idx}]: missing `coeff` field")))?;
        let coeff = eval_coeff(coeff_node, &bindings)?;
        out.push(StencilEntry { offset, coeff });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn rule_json() -> Value {
        json!({
            "discretizations": {
                "ppm_reconstruction": {
                    "stencil": {
                        "q_left_edge": [
                            { "selector": { "offset": -2 }, "coeff": { "op": "/", "args": [-1, 12] } },
                            { "selector": { "offset": -1 }, "coeff": { "op": "/", "args": [ 7, 12] } },
                            { "selector": { "offset":  0 }, "coeff": { "op": "/", "args": [ 7, 12] } },
                            { "selector": { "offset":  1 }, "coeff": { "op": "/", "args": [-1, 12] } }
                        ],
                        "q_right_edge": [
                            { "selector": { "offset": -1 }, "coeff": { "op": "/", "args": [-1, 12] } },
                            { "selector": { "offset":  0 }, "coeff": { "op": "/", "args": [ 7, 12] } },
                            { "selector": { "offset":  1 }, "coeff": { "op": "/", "args": [ 7, 12] } },
                            { "selector": { "offset":  2 }, "coeff": { "op": "/", "args": [-1, 12] } }
                        ]
                    }
                }
            }
        })
    }

    #[test]
    fn parses_named_sub_stencils() {
        let rule = PpmRule::from_json(&rule_json()).unwrap();
        let left: Vec<i64> = rule
            .entries(SubStencil::QLeftEdge)
            .iter()
            .map(|e| e.offset)
            .collect();
        let right: Vec<i64> = rule
            .entries(SubStencil::QRightEdge)
            .iter()
            .map(|e| e.offset)
            .collect();
        assert_eq!(left, vec![-2, -1, 0, 1]);
        assert_eq!(right, vec![-1, 0, 1, 2]);
        // Coefficients are simple rationals — check 7/12 lands without
        // round-off.
        assert!(rule
            .entries(SubStencil::QLeftEdge)
            .iter()
            .any(|e| (e.coeff - 7.0 / 12.0).abs() < 1e-15));
    }

    #[test]
    fn edge_value_matches_cw84_eq_1_6() {
        // CW84 eq. (1.6): q_{i+1/2} = (-q_{i-1} + 7 q_i + 7 q_{i+1} - q_{i+2}) / 12.
        // Picking arbitrary integers makes the expected value exact.
        let rule = PpmRule::from_json(&rule_json()).unwrap();
        let window = [3.0_f64, 5.0, 11.0, -2.0];
        // window[k] corresponds to offset k - 1 in q_right_edge ordering.
        let q_at = |off: i64| match off {
            -1 => window[0],
            0 => window[1],
            1 => window[2],
            2 => window[3],
            other => panic!("unexpected offset {other}"),
        };
        let got = rule
            .evaluate(OutputKind::QRightEdge, q_at)
            .expect("eval succeeds");
        let expected = (-window[0] + 7.0 * window[1] + 7.0 * window[2] - window[3]) / 12.0;
        assert!((got - expected).abs() < 1e-13 * expected.abs().max(1.0));
    }

    #[test]
    fn parabola_preserves_cell_average() {
        // A truly cell-average-preserving (a_L, a_R, q_i) triple has
        // (a_L + a_R)/2 + a_6/6 = q_i. Pick a_L, a_R freely — the
        // formula's a_6 is constructed precisely to enforce this.
        let (a_l, a_r, q_i) = (0.4_f64, 1.7, 1.0);
        // ∫_0^1 a(ξ) dξ = a_L + (da/2) + a_6/6.
        let da = a_r - a_l;
        let a6 = 6.0 * (q_i - 0.5 * (a_l + a_r));
        let avg = a_l + da / 2.0 + a6 / 6.0;
        assert!((avg - q_i).abs() < 1e-15);

        // Direct call into evaluator's parabola path.
        assert_eq!(parabola(a_l, a_r, q_i, 0.0), a_l);
        assert!((parabola(a_l, a_r, q_i, 1.0) - a_r).abs() < 1e-15);
    }

    #[test]
    fn evaluate_dispatch_routes_through_named_sub_stencil() {
        // q_left_edge applied at cell i is the right-edge interpolation
        // shifted one cell left: it samples offsets -2..1.
        let rule = PpmRule::from_json(&rule_json()).unwrap();
        let q = |off: i64| (off + 5) as f64; // arbitrary linear field
        let got = rule.evaluate(OutputKind::QLeftEdge, q).unwrap();
        let expected = (-q(-2) + 7.0 * q(-1) + 7.0 * q(0) - q(1)) / 12.0;
        assert!((got - expected).abs() < 1e-14);
    }

    #[test]
    fn parabola_xi_out_of_range_errors() {
        let rule = PpmRule::from_json(&rule_json()).unwrap();
        let q = |_off: i64| 1.0_f64;
        let err = rule
            .evaluate(OutputKind::QParabola { xi: 1.5 }, q)
            .unwrap_err();
        assert!(matches!(err, PpmRuleError::ParabolaArg(_)));
    }

    #[test]
    fn missing_sub_stencil_errors() {
        let bad = json!({"stencil": {"q_left_edge": []}});
        let err = PpmRule::from_json(&bad).unwrap_err();
        assert!(matches!(
            err,
            PpmRuleError::MissingSubStencil("q_right_edge")
        ));
    }
}
