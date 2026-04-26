/**
 * PPM (Colella & Woodward 1984) reconstruction rule evaluator.
 *
 * Mirrors the Julia reference test in `test/test_ppm_reconstruction_rule.jl`
 * and lifts the same semantics into a reusable runtime that consumes the
 * on-disk `discretizations/finite_volume/ppm_reconstruction.json` rule.
 *
 * Why this lives next to the generic AST evaluator: the PPM rule's `stencil`
 * field is a multi-stencil mapping (`q_left_edge`, `q_right_edge`) — the
 * generic flat-array stencil applicator does not apply. Output dispatch
 * across the rule's three declared outputs (`q_left_edge`, `q_right_edge`,
 * `q_parabola`) lives here.
 *
 * The edge stencils are 4th-order interpolations of cell-averaged data
 * (CW84 eq. 1.6); the parabola sub-cell reconstruction (CW84 eqs. 1.5,
 * 1.7, 1.10) is a closed-form combination of the two edge values and
 * the cell average and is computed without a separate AST.
 */

import { evaluate } from "./expression.js";
import type {
  Bindings,
  Rule,
  RuleSelector,
  RuleStencilEntry,
} from "./types.js";

export type PpmSubStencil = "q_left_edge" | "q_right_edge";
export type PpmOutputKind = PpmSubStencil | "q_parabola";

/** A periodic 1D array of cell averages along the reconstructed axis. */
export type CellAverages = ReadonlyArray<number>;

const DEFAULT_AXIS = "x";

/** Verify the rule's shape matches what this evaluator expects. */
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

/**
 * Resolve the placeholder axis name from `applies_to.dim`. The rule uses the
 * ESS placeholder convention (`"$x"`); the bare name (`"x"`) is also accepted.
 * The leading `$`, if present, is stripped — `selectorMatchesAxis` does the
 * same on the selector side so the two compare equal regardless of which form
 * the catalog file happens to use.
 */
export function resolveAxis(rule: Rule): string {
  const applies = rule.applies_to;
  if (
    typeof applies === "object" &&
    applies !== null &&
    "dim" in applies &&
    typeof (applies as { dim?: unknown }).dim === "string"
  ) {
    const dim = (applies as { dim: string }).dim;
    return stripPlaceholder(dim);
  }
  return DEFAULT_AXIS;
}

function stripPlaceholder(name: string): string {
  return name.startsWith("$") ? name.slice(1) : name;
}

/**
 * Look up the named sub-stencil array on a multi-output rule. Throws if the
 * rule's `stencil` is in flat-array form or the sub-stencil is missing.
 */
export function getSubStencil(
  rule: Rule,
  name: string,
): RuleStencilEntry[] {
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
  // JS `%` keeps the sign of the dividend; wrap to [0, n).
  const m = ((i % n) + n) % n;
  return m;
}

/**
 * Apply a named edge sub-stencil at cell `i` of a periodic 1D field.
 *
 * Returns `Σ coeff_k * q[(i + offset_k) mod n]` where `(coeff_k, offset_k)`
 * are taken from `rule.stencil[subStencil]`. The dollar-prefixed axis name
 * declared by the rule (`$x` → `x`) is used to assert that all selectors
 * target the reconstructed axis; non-matching selectors raise.
 */
export function applyEdgeStencil(
  rule: Rule,
  subStencil: PpmSubStencil,
  q: CellAverages,
  i: number,
  bindings: Bindings = new Map(),
): number {
  if (q.length < 1) {
    throw new RangeError("applyEdgeStencil: cell-average array is empty");
  }
  if (!Number.isInteger(i)) {
    throw new RangeError(`applyEdgeStencil: i=${i} must be an integer`);
  }
  const axis = resolveAxis(rule);
  const entries = getSubStencil(rule, subStencil);
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
    const c = evaluate(entry.coeff, bindings);
    acc += c * q[idx];
  }
  return acc;
}

/**
 * Evaluate the cell-internal PPM parabola at the normalised coordinate
 * `xi ∈ [0, 1]` (CW84 eqs. 1.5, 1.7, 1.10):
 *
 *     a(ξ) = a_L + ξ · (Δa + a_6 · (1 − ξ))
 *     Δa  = a_R − a_L
 *     a_6 = 6 · (q̄ − ½(a_L + a_R))
 *
 * `aL` is the left edge of the cell (= `q_left_edge` for that cell, or
 * equivalently `q_right_edge` of the neighbour to the left), `aR` is the
 * right edge of the cell, `qbar` is the cell average.
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

/**
 * Reconstructed values for a single cell: both edge values plus a closure
 * that samples the parabola at any `xi ∈ [0, 1]`.
 */
export interface CellReconstruction {
  q_left_edge: number;
  q_right_edge: number;
  parabola: (xi: number) => number;
}

/**
 * Reconstruct cell `i` of a periodic 1D field: compute `q_left_edge[i]`,
 * `q_right_edge[i]`, and capture the parabola closure for sub-cell sampling.
 *
 * Note: per CW84, the left edge of cell `i` is the right edge of cell `i-1`,
 * which is also what `stencil.q_left_edge` evaluated at cell `i` returns.
 * Both forms agree to ULP because the rule's `q_left_edge` coefficients
 * are exactly the `q_right_edge` coefficients shifted by one.
 */
export function reconstructCell(
  rule: Rule,
  q: CellAverages,
  i: number,
  bindings: Bindings = new Map(),
): CellReconstruction {
  assertPpmRule(rule);
  const aL = applyEdgeStencil(rule, "q_left_edge", q, i, bindings);
  const aR = applyEdgeStencil(rule, "q_right_edge", q, i, bindings);
  const qi = q[periodicIndex(i, q.length)];
  return {
    q_left_edge: aL,
    q_right_edge: aR,
    parabola: (xi: number) => evaluateParabola(aL, aR, qi, xi),
  };
}

/**
 * One-shot dispatch on the rule's declared `outputs` for a single cell.
 *
 *  - `"q_left_edge"`: returns `applyEdgeStencil(rule, "q_left_edge", …)`
 *  - `"q_right_edge"`: returns `applyEdgeStencil(rule, "q_right_edge", …)`
 *  - `"q_parabola"`: requires `xi`; returns the parabola sampled at `xi`
 */
export function reconstructOutput(
  rule: Rule,
  output: PpmOutputKind,
  q: CellAverages,
  i: number,
  options: { xi?: number; bindings?: Bindings } = {},
): number {
  const bindings = options.bindings ?? new Map<string, number>();
  if (output === "q_left_edge" || output === "q_right_edge") {
    return applyEdgeStencil(rule, output, q, i, bindings);
  }
  if (output === "q_parabola") {
    if (typeof options.xi !== "number") {
      throw new Error(
        "reconstructOutput: output 'q_parabola' requires options.xi",
      );
    }
    const cell = reconstructCell(rule, q, i, bindings);
    return cell.parabola(options.xi);
  }
  throw new Error(`reconstructOutput: unknown output kind '${String(output)}'`);
}
