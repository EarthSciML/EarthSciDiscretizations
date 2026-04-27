/**
 * Test-private PPM (Colella & Woodward 1984) reconstruction helpers.
 *
 * Inlines the multi-stencil dispatch + parabola formula required by the
 * `ppm_reconstruction` conformance / convergence / unit tests. Per
 * `AGENTS.md`, ESD does not expose per-rule shadow evaluators in its
 * public API; the canonical Julia test for this rule
 * (`test/test_ppm_reconstruction_rule.jl`) inlines the same math, and
 * this module mirrors that pattern for the TypeScript binding.
 *
 * The edge stencils (`q_left_edge`, `q_right_edge`) carry their
 * coefficients as ESS expression-AST nodes, so coefficient evaluation
 * goes through the canonical generic evaluator
 * (`rules.evaluate`); only PPM-specific pieces (sub-stencil dispatch,
 * the closed-form parabola) live here.
 */

import { rules } from "../src/index.js";
import type {
  Rule,
  RuleSelector,
  RuleStencilEntry,
} from "../src/rules/index.js";

export type PpmSubStencil = "q_left_edge" | "q_right_edge";
export type PpmOutputKind = PpmSubStencil | "q_parabola";

export type CellAverages = ReadonlyArray<number>;

const DEFAULT_AXIS = "x";

function stripPlaceholder(name: string): string {
  return name.startsWith("$") ? name.slice(1) : name;
}

export function assertPpmRule(rule: Rule): void {
  if (rule.name !== "ppm_reconstruction") {
    throw new Error(
      `assertPpmRule: rule name must be 'ppm_reconstruction', got '${rule.name}'`,
    );
  }
  if (rule.grid_family !== "cartesian") {
    throw new Error(
      `assertPpmRule: ppm_reconstruction grid_family must be 'cartesian', got '${rule.grid_family}'`,
    );
  }
  if (Array.isArray(rule.stencil)) {
    throw new Error(
      "assertPpmRule: ppm_reconstruction stencil must be a multi-stencil object (q_left_edge / q_right_edge), not an array",
    );
  }
  for (const sub of ["q_left_edge", "q_right_edge"] as const) {
    const entries = rule.stencil[sub];
    if (!Array.isArray(entries) || entries.length === 0) {
      throw new Error(
        `assertPpmRule: ppm_reconstruction stencil.${sub} must be a non-empty array of entries`,
      );
    }
  }
}

export function resolveAxis(rule: Rule): string {
  const applies = rule.applies_to;
  if (
    typeof applies === "object" &&
    applies !== null &&
    "dim" in applies &&
    typeof (applies as { dim?: unknown }).dim === "string"
  ) {
    return stripPlaceholder((applies as { dim: string }).dim);
  }
  return DEFAULT_AXIS;
}

export function getSubStencil(rule: Rule, name: string): RuleStencilEntry[] {
  if (Array.isArray(rule.stencil)) {
    throw new Error(
      `getSubStencil: rule '${rule.name}' has flat stencil array; cannot dispatch to sub-stencil '${name}'`,
    );
  }
  const entries = rule.stencil[name];
  if (!Array.isArray(entries)) {
    throw new Error(
      `getSubStencil: rule '${rule.name}' has no sub-stencil named '${name}'`,
    );
  }
  return entries;
}

function selectorAxis(sel: RuleSelector | RuleSelector[]): string {
  if (Array.isArray(sel)) {
    if (sel.length === 0) {
      throw new Error("PPM stencil entry selector array is empty");
    }
    return selectorAxis(sel[0]);
  }
  if (typeof sel.axis !== "string") {
    throw new Error(
      `PPM stencil entry selector missing 'axis' string: ${JSON.stringify(sel)}`,
    );
  }
  return sel.axis;
}

function selectorOffset(sel: RuleSelector | RuleSelector[]): number {
  if (Array.isArray(sel)) {
    if (sel.length === 0) {
      throw new Error("PPM stencil entry selector array is empty");
    }
    return selectorOffset(sel[0]);
  }
  if (typeof sel.offset !== "number") {
    throw new Error(
      `PPM stencil entry selector missing 'offset' number: ${JSON.stringify(sel)}`,
    );
  }
  return sel.offset;
}

function periodicIndex(i: number, n: number): number {
  return ((i % n) + n) % n;
}

export function applyEdgeStencil(
  rule: Rule,
  subStencil: PpmSubStencil,
  q: CellAverages,
  i: number,
): number {
  if (q.length < 1) {
    throw new RangeError("applyEdgeStencil: cell-average array is empty");
  }
  if (!Number.isInteger(i)) {
    throw new RangeError(`applyEdgeStencil: i=${i} must be an integer`);
  }
  const axis = resolveAxis(rule);
  const entries = getSubStencil(rule, subStencil);
  const bindings = new Map<string, number>();
  let acc = 0;
  for (const entry of entries) {
    const entryAxis = stripPlaceholder(selectorAxis(entry.selector));
    if (entryAxis !== axis) {
      throw new Error(
        `applyEdgeStencil: sub-stencil ${subStencil} entry targets axis '${entryAxis}', expected '${axis}'`,
      );
    }
    const off = selectorOffset(entry.selector);
    const idx = periodicIndex(i + off, q.length);
    const c = rules.evaluate(entry.coeff, bindings);
    acc += c * q[idx];
  }
  return acc;
}

/**
 * Evaluate the PPM parabola at `xi ∈ [0, 1]` (CW84 eqs. 1.5, 1.7, 1.10):
 *   a(ξ) = a_L + ξ · (Δa + a_6 · (1 − ξ))
 *   Δa  = a_R − a_L
 *   a_6 = 6 · (q̄ − ½(a_L + a_R))
 */
export function evaluateParabola(
  aL: number,
  aR: number,
  qbar: number,
  xi: number,
): number {
  if (!(xi >= 0 && xi <= 1)) {
    throw new RangeError(
      `evaluateParabola: xi=${xi} must lie in [0, 1] (cell-normalised coordinate)`,
    );
  }
  const da = aR - aL;
  const a6 = 6 * (qbar - 0.5 * (aL + aR));
  return aL + xi * (da + a6 * (1 - xi));
}

export interface CellReconstruction {
  q_left_edge: number;
  q_right_edge: number;
  parabola: (xi: number) => number;
}

export function reconstructCell(
  rule: Rule,
  q: CellAverages,
  i: number,
): CellReconstruction {
  assertPpmRule(rule);
  const aL = applyEdgeStencil(rule, "q_left_edge", q, i);
  const aR = applyEdgeStencil(rule, "q_right_edge", q, i);
  const qi = q[periodicIndex(i, q.length)];
  return {
    q_left_edge: aL,
    q_right_edge: aR,
    parabola: (xi: number) => evaluateParabola(aL, aR, qi, xi),
  };
}

export function reconstructOutput(
  rule: Rule,
  output: PpmOutputKind,
  q: CellAverages,
  i: number,
  options: { xi?: number } = {},
): number {
  if (output === "q_left_edge" || output === "q_right_edge") {
    return applyEdgeStencil(rule, output, q, i);
  }
  if (output === "q_parabola") {
    if (typeof options.xi !== "number") {
      throw new Error(
        "reconstructOutput: output 'q_parabola' requires options.xi",
      );
    }
    const cell = reconstructCell(rule, q, i);
    return cell.parabola(options.xi);
  }
  throw new Error(`reconstructOutput: unknown output kind '${String(output)}'`);
}
