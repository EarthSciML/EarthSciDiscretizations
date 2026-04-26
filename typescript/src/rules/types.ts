/**
 * Types for ESS-emitted discretization rule JSON (esm-spec §7).
 *
 * A rule's `coeff` field is a JSON-decoded `ExpressionNode`: a number, a
 * variable name (string), or an `{op, args}` object whose args are themselves
 * `ExpressionNode`s. Mirrors the `ExpressionNode` carried by the Julia,
 * Python, and Rust bindings.
 */

export type ExpressionNode = number | string | OpNode;

export interface OpNode {
  op: string;
  args: ExpressionNode[];
  name?: string;
  value?: number | number[];
}

export interface RuleSelector {
  kind: string;
  axis?: string;
  offset?: number;
  [key: string]: unknown;
}

export interface RuleStencilEntry {
  selector: RuleSelector | RuleSelector[];
  coeff: ExpressionNode;
}

/**
 * A rule's `stencil` is either:
 *  - a flat array of entries (single linear combination, the common case), or
 *  - an object mapping sub-stencil names → entry arrays (multi-output rules
 *    such as `ppm_reconstruction`, where the same support set produces several
 *    distinct linear combinations like `q_left_edge` / `q_right_edge`).
 *
 * The on-disk JSON shape is preserved verbatim; downstream evaluators dispatch
 * on type via `Array.isArray`.
 */
export type RuleStencil =
  | RuleStencilEntry[]
  | Record<string, RuleStencilEntry[]>;

export interface Rule {
  name: string;
  family: string;
  applies_to: ExpressionNode;
  grid_family: string;
  combine?: string;
  accuracy?: string;
  stencil: RuleStencil;
  /** Declared output kinds for multi-output rules (e.g. PPM). */
  outputs?: string[];
  /** Reconstruction form (e.g. `"piecewise_parabolic"`). */
  form?: string;
  /** Limiter selection (e.g. `"none"`, `"monotonic"`). */
  limiter?: string;
  /** Reference citation (e.g. CW84). */
  reference?: string;
}

export type Bindings = Map<string, number>;
