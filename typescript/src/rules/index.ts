/**
 * Rule evaluation runtime.
 *
 * Coefficient evaluation delegates to `earthsci-toolkit` (ESS) so the TS
 * binding shares one implementation across all languages — see Julia
 * `src/rule_eval.jl`, which similarly passes through to
 * `EarthSciSerialization.evaluate`.
 *
 * Catalog loading (`loadCatalog` etc.) is ESD-specific: ESS validates the
 * `discretizations` field of an EsmFile but has no helper for the
 * standalone `{discretizations: {...}}` per-rule JSON shape that lives
 * under `<repo>/discretizations/<family>/*.json`.
 */

import { evaluate as essEvaluate } from "earthsci-toolkit";
import type { Bindings, ExpressionNode } from "./types.js";

type EssExpr = Parameters<typeof essEvaluate>[0];

export function evaluate(node: ExpressionNode, bindings: Bindings): number {
  return essEvaluate(node as EssExpr, bindings);
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
  return essEvaluate(node as EssExpr, map);
}

export { loadCatalog, loadRuleFile, parseRuleFile, findRule } from "./load.js";
export type {
  Bindings,
  ExpressionNode,
  OpNode,
  Rule,
  RuleSelector,
  RuleStencil,
  RuleStencilEntry,
} from "./types.js";
