//! Scalar evaluator for ESS expression-AST coefficient nodes.
//!
//! ESD rule files (`discretizations/<family>/*.json`) carry stencil
//! coefficients as ESS `ExpressionNode` AST trees per `esm-spec.md` §4.2.
//! This module evaluates a JSON-decoded AST against a bindings table and
//! returns a `f64`, mirroring the reference Julia evaluator
//! (`EarthSciSerialization.evaluate(::OpExpr, ::Dict{String,Float64})`).
//! The supported scalar op set is the closed registry implemented by ESS:
//! arithmetic (`+ - * / ^`), elementary functions (`sin cos tan exp log
//! sqrt abs`), spec-pinned constants (`π`/`pi`, `e`), and the `const`
//! literal carrier. Array, calculus, conditional, comparison, logical,
//! `fn`, and `enum` ops are out of scope: rule files that require them
//! must extend the evaluator alongside the matching ESS contract.

use std::collections::HashMap;

use serde_json::{Map, Value};
use thiserror::Error;

/// Variable bindings consumed by [`eval_coeff`]. Keys are AST variable
/// names (e.g. `"R"`, `"cos_lat"`, `"dlon"`); values are the per-cell
/// floats supplied by the binding's grid accessor.
pub type Bindings = HashMap<String, f64>;

/// Errors raised while walking a rule-AST coefficient node.
#[derive(Debug, Error)]
pub enum RuleEvalError {
    /// A string variable reference appeared in the AST but has no entry
    /// in `bindings`.
    #[error("unbound variable: {0}")]
    UnboundVariable(String),
    /// The op string is not part of the closed scalar registry this
    /// evaluator implements.
    #[error("unsupported op: {0}")]
    UnsupportedOp(String),
    /// The AST violates the ESS schema (missing field, wrong shape, …).
    #[error("schema violation: {0}")]
    Schema(String),
    /// An op was called with the wrong number of arguments.
    #[error("invalid arity for op {op}: expected {expected}, got {got}")]
    Arity {
        op: String,
        expected: String,
        got: usize,
    },
    /// A numerical operation has no real-valued result for the given
    /// inputs (e.g. `sqrt` of a negative, `log` of a non-positive).
    #[error("domain error: {0}")]
    Domain(String),
}

/// Evaluate a JSON-decoded ESS expression-AST coefficient node.
///
/// `node` is a JSON number, a JSON string (variable reference), or an
/// object carrying `"op"` and `"args"` (and optional `"value"` for
/// `const`). `bindings` maps variable names to floats.
pub fn eval_coeff(node: &Value, bindings: &Bindings) -> Result<f64, RuleEvalError> {
    match node {
        Value::Number(n) => n
            .as_f64()
            .ok_or_else(|| RuleEvalError::Schema(format!("non-finite numeric literal: {n}"))),
        Value::String(s) => bindings
            .get(s)
            .copied()
            .ok_or_else(|| RuleEvalError::UnboundVariable(s.clone())),
        Value::Object(map) => eval_op(map, bindings),
        other => Err(RuleEvalError::Schema(format!(
            "AST node must be number, string, or object; got {other}"
        ))),
    }
}

fn eval_op(map: &Map<String, Value>, bindings: &Bindings) -> Result<f64, RuleEvalError> {
    let op = map
        .get("op")
        .and_then(Value::as_str)
        .ok_or_else(|| RuleEvalError::Schema("AST object missing string `op` field".into()))?;

    if op == "const" {
        let v = map
            .get("value")
            .ok_or_else(|| RuleEvalError::Schema("`const` op missing `value` field".into()))?;
        return match v {
            Value::Number(n) => n
                .as_f64()
                .ok_or_else(|| RuleEvalError::Schema(format!("non-finite `const` value: {n}"))),
            other => Err(RuleEvalError::Schema(format!(
                "scalar `const` expected, got {other}"
            ))),
        };
    }
    if op == "enum" {
        return Err(RuleEvalError::Schema(
            "`enum` op encountered during evaluation; expected lowering to `const`".into(),
        ));
    }
    if op == "fn" {
        let name = map
            .get("name")
            .and_then(Value::as_str)
            .unwrap_or("<unnamed>");
        return Err(RuleEvalError::UnsupportedOp(format!("fn:{name}")));
    }

    let args = map
        .get("args")
        .and_then(Value::as_array)
        .ok_or_else(|| RuleEvalError::Schema(format!("op `{op}` missing `args` array")))?;
    let mut values = Vec::with_capacity(args.len());
    for a in args {
        values.push(eval_coeff(a, bindings)?);
    }
    apply_op(op, &values)
}

fn apply_op(op: &str, args: &[f64]) -> Result<f64, RuleEvalError> {
    let arity = |expected: &str| RuleEvalError::Arity {
        op: op.to_string(),
        expected: expected.to_string(),
        got: args.len(),
    };
    match op {
        // Arithmetic ----------------------------------------------------
        "+" => match args.len() {
            0 => Ok(0.0),
            1 => Ok(args[0]),
            // Iterator::sum on f64 folds left-to-right starting from 0.0,
            // matching the Julia reference (`sum(args)`).
            _ => Ok(args.iter().sum()),
        },
        "-" => match args.len() {
            1 => Ok(-args[0]),
            2 => Ok(args[0] - args[1]),
            _ => Err(arity("1 or 2")),
        },
        "*" => match args.len() {
            0 => Ok(1.0),
            1 => Ok(args[0]),
            // Iterator::product folds left-to-right starting from 1.0,
            // matching the Julia reference (`prod(args)`).
            _ => Ok(args.iter().product()),
        },
        "/" => {
            if args.len() != 2 {
                return Err(arity("2"));
            }
            Ok(args[0] / args[1])
        }
        "^" => {
            if args.len() != 2 {
                return Err(arity("2"));
            }
            Ok(args[0].powf(args[1]))
        }
        // Elementary functions -----------------------------------------
        "sin" | "cos" | "tan" | "exp" | "log" | "sqrt" | "abs" => {
            if args.len() != 1 {
                return Err(arity("1"));
            }
            let x = args[0];
            Ok(match op {
                "sin" => x.sin(),
                "cos" => x.cos(),
                "tan" => x.tan(),
                "exp" => x.exp(),
                "log" => {
                    if x <= 0.0 {
                        return Err(RuleEvalError::Domain(format!("log({x})")));
                    }
                    x.ln()
                }
                "sqrt" => {
                    if x < 0.0 {
                        return Err(RuleEvalError::Domain(format!("sqrt({x})")));
                    }
                    x.sqrt()
                }
                "abs" => x.abs(),
                _ => unreachable!(),
            })
        }
        // Constants -----------------------------------------------------
        "π" | "pi" => {
            if !args.is_empty() {
                return Err(arity("0"));
            }
            Ok(std::f64::consts::PI)
        }
        "e" => {
            if !args.is_empty() {
                return Err(arity("0"));
            }
            Ok(std::f64::consts::E)
        }
        other => Err(RuleEvalError::UnsupportedOp(other.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn b(pairs: &[(&str, f64)]) -> Bindings {
        pairs.iter().map(|(k, v)| ((*k).to_string(), *v)).collect()
    }

    #[test]
    fn number_and_variable() {
        let bindings = b(&[("x", 3.5)]);
        assert_eq!(eval_coeff(&json!(2.5), &bindings).unwrap(), 2.5);
        assert_eq!(eval_coeff(&json!(-7), &bindings).unwrap(), -7.0);
        assert_eq!(eval_coeff(&json!("x"), &bindings).unwrap(), 3.5);
    }

    #[test]
    fn unbound_variable_errors() {
        let bindings = Bindings::new();
        let err = eval_coeff(&json!("y"), &bindings).unwrap_err();
        matches!(err, RuleEvalError::UnboundVariable(name) if name == "y");
    }

    #[test]
    fn nary_arithmetic() {
        let bindings = b(&[("a", 2.0), ("b", 3.0), ("c", 5.0)]);
        // a + b + c
        let ast = json!({"op": "+", "args": ["a", "b", "c"]});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), 10.0);
        // a * b * c
        let ast = json!({"op": "*", "args": ["a", "b", "c"]});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), 30.0);
        // -1 / (2 * a * b)
        let ast = json!({
            "op": "/",
            "args": [-1, {"op": "*", "args": [2, "a", "b"]}]
        });
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), -1.0 / 12.0);
    }

    #[test]
    fn unary_minus_and_const() {
        let bindings = b(&[("x", 4.0)]);
        let ast = json!({"op": "-", "args": ["x"]});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), -4.0);
        let ast = json!({"op": "const", "value": 7.5, "args": []});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), 7.5);
    }

    #[test]
    fn elementary_functions() {
        let bindings = b(&[("x", 0.5)]);
        let ast = json!({"op": "sin", "args": ["x"]});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), 0.5f64.sin());
        let ast = json!({"op": "sqrt", "args": [{"op": "*", "args": ["x", 8]}]});
        assert_eq!(eval_coeff(&ast, &bindings).unwrap(), 4.0f64.sqrt());
    }

    #[test]
    fn unsupported_op_errors() {
        let bindings = Bindings::new();
        let ast = json!({"op": "ifelse", "args": [true, 1, 0]});
        let err = eval_coeff(&ast, &bindings).unwrap_err();
        matches!(err, RuleEvalError::UnsupportedOp(name) if name == "ifelse");
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
}
