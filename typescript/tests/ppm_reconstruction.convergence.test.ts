/**
 * MMS convergence sweep for the TypeScript ppm_reconstruction evaluator.
 *
 * Mirrors the Julia sibling test in `test/test_ppm_reconstruction_rule.jl`.
 * Drives the on-disk multi-stencil rule through `applyEdgeStencil` +
 * `evaluateParabola` on a uniform 1D periodic Cartesian grid, sweeping
 * `n ∈ {16, 32, 64, 128}`, and asserts that the L_infinity error against
 * the smooth periodic manufactured solution decays at the rule-declared
 * order.
 *
 * The convergence input/expected fixture lives at
 * `tests/fixtures/ppm_reconstruction/{input,expected}.esm`. Per CW84, the
 * unlimited PPM reconstruction is theoretically O(dx^3); the fixture's
 * `expected_min_order = 2.7` allows headroom for the pre-asymptotic regime
 * on the coarsest grid.
 */

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { rules } from "../src/index.js";
import {
  reconstructCell,
  reconstructOutput,
} from "./ppm_reconstruction_helpers.js";

const REPO_ROOT = resolve(fileURLToPath(new URL("../", import.meta.url)), "..");
const RULE_PATH = join(
  REPO_ROOT,
  "discretizations",
  "finite_volume",
  "ppm_reconstruction.json",
);
const FIXTURE_DIR = join(REPO_ROOT, "tests", "fixtures", "ppm_reconstruction");

interface ConvergenceInput {
  rule: string;
  domain: { extent: Array<[number, number]>; boundary: "periodic" };
  manufactured_solution: { expression: string };
  grids: Array<{ n: number; uniform: boolean; dx: number }>;
  norm: string;
  samples_per_cell: number;
}

interface ConvergenceExpected {
  rule: string;
  expected_min_order: number;
  theoretical_order: number;
  norm: string;
}

// f(x) = sin(2*pi*x) + 0.3*cos(4*pi*x), antiderivative used for exact cell averages.
const TWO_PI = 2 * Math.PI;
const FOUR_PI = 4 * Math.PI;
const f = (x: number): number =>
  Math.sin(TWO_PI * x) + 0.3 * Math.cos(FOUR_PI * x);
const F = (x: number): number =>
  -Math.cos(TWO_PI * x) / TWO_PI + (0.3 * Math.sin(FOUR_PI * x)) / FOUR_PI;

describe("ppm_reconstruction MMS convergence (TypeScript)", () => {
  const input: ConvergenceInput = JSON.parse(
    readFileSync(join(FIXTURE_DIR, "input.esm"), "utf8"),
  );
  const expected: ConvergenceExpected = JSON.parse(
    readFileSync(join(FIXTURE_DIR, "expected.esm"), "utf8"),
  );

  const ruleSet = rules.parseRuleFile(
    RULE_PATH,
    readFileSync(RULE_PATH, "utf8"),
  );
  const rule = ruleSet.find((r) => r.name === input.rule);
  if (!rule) throw new Error(`rule ${input.rule} missing from ${RULE_PATH}`);

  it("rule + fixtures agree on identity and norm", () => {
    expect(input.rule).toBe("ppm_reconstruction");
    expect(expected.rule).toBe("ppm_reconstruction");
    expect(expected.norm).toBe("Linf");
    expect(input.norm).toBe("Linf");
  });

  it("sweep over n ∈ {16,32,64,128} hits the declared min order", () => {
    const samplesPerCell = input.samples_per_cell;
    const grids = input.grids.map((g) => g.n);

    const linfError = (n: number): number => {
      const dx = 1.0 / n;
      const q = new Array<number>(n);
      for (let i = 0; i < n; i++) {
        const a = i * dx;
        const b = (i + 1) * dx;
        q[i] = (F(b) - F(a)) / (b - a);
      }
      let err = 0.0;
      for (let i = 0; i < n; i++) {
        const cell = reconstructCell(rule, q, i);
        for (let k = 0; k < samplesPerCell; k++) {
          const xi = (k + 0.5) / samplesPerCell;
          const x = i * dx + xi * dx;
          const diff = Math.abs(cell.parabola(xi) - f(x));
          if (diff > err) err = diff;
        }
      }
      return err;
    };

    const errors = grids.map((n) => linfError(n));
    for (const e of errors) expect(Number.isFinite(e) && e > 0).toBe(true);

    // Errors must decrease monotonically as the grid refines.
    for (let i = 0; i + 1 < errors.length; i++) {
      expect(
        errors[i + 1] < errors[i],
        `errors[${i + 1}]=${errors[i + 1]} not < errors[${i}]=${errors[i]}`,
      ).toBe(true);
    }

    const orders: number[] = [];
    for (let i = 0; i + 1 < errors.length; i++) {
      orders.push(Math.log2(errors[i] / errors[i + 1]));
    }
    const minOrder = Math.min(...orders);
    expect(
      minOrder >= expected.expected_min_order,
      `min observed order ${minOrder} < expected ${expected.expected_min_order}; orders=${orders.join(", ")}`,
    ).toBe(true);
    // Theoretical CW84 order is 3; we should be in the right ballpark.
    expect(minOrder).toBeGreaterThan(expected.theoretical_order - 0.5);
  });

  it("reconstructOutput dispatch matches reconstructCell for each output kind", () => {
    const n = 32;
    const dx = 1.0 / n;
    const q = new Array<number>(n);
    for (let i = 0; i < n; i++) {
      const a = i * dx;
      const b = (i + 1) * dx;
      q[i] = (F(b) - F(a)) / (b - a);
    }
    for (const i of [0, 5, 13, n - 1]) {
      const cell = reconstructCell(rule, q, i);
      expect(reconstructOutput(rule, "q_left_edge", q, i)).toBe(
        cell.q_left_edge,
      );
      expect(reconstructOutput(rule, "q_right_edge", q, i)).toBe(
        cell.q_right_edge,
      );
      for (const xi of [0.0, 0.25, 0.5, 0.75, 1.0]) {
        expect(reconstructOutput(rule, "q_parabola", q, i, { xi })).toBe(
          cell.parabola(xi),
        );
      }
    }
  });
});
