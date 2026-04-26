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

export interface Rule {
  name: string;
  family: string;
  applies_to: ExpressionNode;
  grid_family: string;
  combine?: string;
  accuracy?: string;
  stencil: RuleStencilEntry[];
}

export type Bindings = Map<string, number>;
