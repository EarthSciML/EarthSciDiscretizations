/**
 * Tree-walk evaluator for ESS-emitted ExpressionNode coefficients.
 *
 * Mirrors `EarthSciSerialization.evaluate` (Julia) and the Julia entry point
 * `EarthSciDiscretizations.eval_coeff` so the TypeScript binding produces
 * numerically identical results when fed the same AST and bindings.
 *
 * Supported ops: arithmetic (`+ - * / ^`), unary `-`, `abs`, `sign`,
 * trig (`sin cos tan asin acos atan atan2`), `exp log log10 sqrt`,
 * `min max floor ceil`, comparisons (`> < >= <= == !=`), boolean
 * (`and or not`), `ifelse`, and scalar `const`.
 */

import type { ExpressionNode, Bindings, OpNode } from "./types.js";

export class UnboundVariableError extends Error {
  readonly variableName: string;
  constructor(variableName: string) {
    super(`Variable '${variableName}' not found in bindings`);
    this.variableName = variableName;
    this.name = "UnboundVariableError";
  }
}

function isOpNode(node: ExpressionNode): node is OpNode {
  return typeof node === "object" && node !== null && typeof (node as OpNode).op === "string";
}

export function evaluate(node: ExpressionNode, bindings: Bindings): number {
  if (typeof node === "number") return node;
  if (typeof node === "string") {
    const v = bindings.get(node);
    if (v === undefined) throw new UnboundVariableError(node);
    return v;
  }
  if (!isOpNode(node)) {
    throw new Error(`unrecognized AST node type: ${typeof node}`);
  }

  const { op } = node;

  if (op === "const") {
    const v = node.value;
    if (typeof v === "number") return v;
    throw new Error("`const` op with non-scalar value cannot be evaluated as a number");
  }
  if (op === "enum") {
    throw new Error("`enum` op encountered during evaluate(); enum nodes must be lowered to `const` integers first");
  }

  const args = node.args.map((a) => evaluate(a, bindings));

  switch (op) {
    case "+":
      if (args.length === 1) return args[0];
      return args.reduce((s, v) => s + v, 0);
    case "-":
      if (args.length === 1) return -args[0];
      if (args.length === 2) return args[0] - args[1];
      throw new Error(`Subtraction requires 1 or 2 arguments, got ${args.length}`);
    case "*":
      if (args.length === 1) return args[0];
      return args.reduce((p, v) => p * v, 1);
    case "/":
      if (args.length !== 2) throw new Error("Division requires exactly 2 arguments");
      if (args[1] === 0) throw new Error("Division by zero");
      return args[0] / args[1];
    case "^":
      if (args.length !== 2) throw new Error("Exponentiation requires exactly 2 arguments");
      return Math.pow(args[0], args[1]);
    case "abs":
      return Math.abs(args[0]);
    case "sign":
      return Math.sign(args[0]);
    case "sqrt":
      if (args[0] < 0) throw new Error("sqrt argument must be non-negative");
      return Math.sqrt(args[0]);
    case "exp":
      return Math.exp(args[0]);
    case "log":
      if (args[0] <= 0) throw new Error("log argument must be positive");
      return Math.log(args[0]);
    case "log10":
      if (args[0] <= 0) throw new Error("log10 argument must be positive");
      return Math.log10(args[0]);
    case "sin":
      return Math.sin(args[0]);
    case "cos":
      return Math.cos(args[0]);
    case "tan":
      return Math.tan(args[0]);
    case "asin":
      if (args[0] < -1 || args[0] > 1) throw new Error("asin argument must be in [-1, 1]");
      return Math.asin(args[0]);
    case "acos":
      if (args[0] < -1 || args[0] > 1) throw new Error("acos argument must be in [-1, 1]");
      return Math.acos(args[0]);
    case "atan":
      return Math.atan(args[0]);
    case "atan2":
      if (args.length !== 2) throw new Error("atan2 requires exactly 2 arguments");
      return Math.atan2(args[0], args[1]);
    case "floor":
      return Math.floor(args[0]);
    case "ceil":
      return Math.ceil(args[0]);
    case "min":
      if (args.length === 0) throw new Error("min requires at least 1 argument");
      return Math.min(...args);
    case "max":
      if (args.length === 0) throw new Error("max requires at least 1 argument");
      return Math.max(...args);
    case ">":
      return args[0] > args[1] ? 1 : 0;
    case "<":
      return args[0] < args[1] ? 1 : 0;
    case ">=":
      return args[0] >= args[1] ? 1 : 0;
    case "<=":
      return args[0] <= args[1] ? 1 : 0;
    case "==":
      return args[0] === args[1] ? 1 : 0;
    case "!=":
      return args[0] !== args[1] ? 1 : 0;
    case "and":
      return args.every((x) => x !== 0) ? 1 : 0;
    case "or":
      return args.some((x) => x !== 0) ? 1 : 0;
    case "not":
      return args[0] === 0 ? 1 : 0;
    case "ifelse":
      if (args.length !== 3) throw new Error("ifelse requires exactly 3 arguments");
      return args[0] !== 0 ? args[1] : args[2];
    default:
      throw new Error(`Unsupported operator: ${op}`);
  }
}

/**
 * Convenience overload that accepts a plain object of bindings instead of a
 * `Map`. Mirrors the Julia `eval_coeff(node, Dict{String,Float64})` surface.
 */
export function evalCoeff(
  node: ExpressionNode,
  bindings: Bindings | Record<string, number>,
): number {
  const map = bindings instanceof Map ? bindings : new Map(Object.entries(bindings));
  return evaluate(node, map);
}
