/**
 * Rule evaluation runtime.
 *
 * Mirrors the Julia surface in `src/rule_eval.jl`: a thin `evalCoeff` (and
 * underlying `evaluate`) that consumes ESS-emitted ExpressionNode coefficients
 * against numeric bindings, plus loaders for the on-disk discretization
 * catalog.
 */

export { evaluate, evalCoeff, UnboundVariableError } from "./expression.js";
export { loadCatalog, loadRuleFile, parseRuleFile, findRule } from "./load.js";
export {
  applyEdgeStencil,
  assertPpmRule,
  evaluateParabola,
  getSubStencil,
  reconstructCell,
  reconstructOutput,
  resolveAxis,
} from "./ppm_reconstruction.js";
export type {
  CellAverages,
  CellReconstruction,
  PpmOutputKind,
  PpmSubStencil,
} from "./ppm_reconstruction.js";
export type {
  Bindings,
  ExpressionNode,
  OpNode,
  Rule,
  RuleSelector,
  RuleStencil,
  RuleStencilEntry,
} from "./types.js";
