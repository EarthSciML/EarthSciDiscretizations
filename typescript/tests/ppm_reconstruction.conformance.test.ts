/**
 * Cross-binding conformance for the ppm_reconstruction rule evaluator.
 *
 * Loads the on-disk rule from `discretizations/finite_volume/` and the
 * shared fixture + golden under
 * `tests/conformance/grids/cartesian/rules/ppm_reconstruction/`, then
 * verifies that the TypeScript PPM evaluator reproduces:
 *
 *   1. The cell averages of the manufactured solution (sanity check on the
 *      analytical antiderivative used by the regenerator).
 *   2. `q_left_edge[i]` and `q_right_edge[i]` at every query cell, as
 *      computed by `applyEdgeStencil` against the rule's multi-stencil.
 *   3. `q_parabola(i, xi)` for every query cell × xi sample, dispatched
 *      through `reconstructOutput` (output_kind = "q_parabola").
 *
 * The Python and Rust sibling evaluators consume the same fixture + golden,
 * so passing tests across all three bindings imply cross-binding agreement
 * to the declared tolerance.
 */

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { rules } from "../src/index.js";

const REPO_ROOT = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
);
const RULE_PATH = join(
  REPO_ROOT,
  "discretizations",
  "finite_volume",
  "ppm_reconstruction.json",
);
const FIXTURE_DIR = join(
  REPO_ROOT,
  "tests",
  "conformance",
  "grids",
  "cartesian",
  "rules",
  "ppm_reconstruction",
);

interface GridSpec {
  family: "cartesian";
  ndim: 1;
  n: number;
  extent: Array<[number, number]>;
  boundary: "periodic";
  dtype: "float64";
}

interface FixtureSpec {
  rule: string;
  tolerance: { relative: number };
  manufactured_solution: {
    name: string;
    expression: string;
    antiderivative: string;
  };
  grid: GridSpec;
  query_cells: number[];
  parabola_xi: number[];
}

interface Golden {
  rule: string;
  grid: GridSpec;
  cell_averages: number[];
  query_cells: number[];
  q_left_edge: number[];
  q_right_edge: number[];
  parabola_xi: number[];
  q_parabola: number[][];
}

function closeRel(a: number, b: number, tol: number): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function buildCellAverages(grid: GridSpec): number[] {
  const n = grid.n;
  const [lo, hi] = grid.extent[0];
  const dx = (hi - lo) / n;
  const TWO_PI = 2 * Math.PI;
  const FOUR_PI = 4 * Math.PI;
  const F = (x: number): number =>
    -Math.cos(TWO_PI * x) / TWO_PI + (0.3 * Math.sin(FOUR_PI * x)) / FOUR_PI;
  const out = new Array<number>(n);
  for (let i = 0; i < n; i++) {
    const a = lo + i * dx;
    const b = lo + (i + 1) * dx;
    out[i] = (F(b) - F(a)) / (b - a);
  }
  return out;
}

describe("ppm_reconstruction rule eval cross-binding conformance", () => {
  const fixture: FixtureSpec = JSON.parse(
    readFileSync(join(FIXTURE_DIR, "fixtures.json"), "utf8"),
  );
  const golden: Golden = JSON.parse(
    readFileSync(join(FIXTURE_DIR, "golden.json"), "utf8"),
  );

  const ruleSet = rules.parseRuleFile(
    RULE_PATH,
    readFileSync(RULE_PATH, "utf8"),
  );
  const rule = ruleSet.find((r) => r.name === fixture.rule);
  if (!rule) throw new Error(`rule ${fixture.rule} missing from ${RULE_PATH}`);

  const tol = fixture.tolerance.relative;
  const cellAverages = buildCellAverages(fixture.grid);

  it("rule shape passes the PPM contract assertions", () => {
    expect(() => rules.assertPpmRule(rule)).not.toThrow();
    expect(rule.family).toBe("finite_volume");
    expect(rule.grid_family).toBe("cartesian");
    expect(rules.resolveAxis(rule)).toBe("x");
    const left = rules.getSubStencil(rule, "q_left_edge");
    const right = rules.getSubStencil(rule, "q_right_edge");
    expect(left.length).toBe(4);
    expect(right.length).toBe(4);
  });

  it("computes cell averages that match the golden", () => {
    expect(cellAverages.length).toBe(golden.cell_averages.length);
    for (let i = 0; i < cellAverages.length; i++) {
      expect(closeRel(cellAverages[i], golden.cell_averages[i], tol)).toBe(true);
    }
  });

  it("matches the golden q_left_edge / q_right_edge at every query cell", () => {
    expect(fixture.query_cells.length).toBe(golden.query_cells.length);
    for (let k = 0; k < fixture.query_cells.length; k++) {
      const i = fixture.query_cells[k];
      expect(i).toBe(golden.query_cells[k]);
      const aL = rules.applyEdgeStencil(rule, "q_left_edge", cellAverages, i);
      const aR = rules.applyEdgeStencil(rule, "q_right_edge", cellAverages, i);
      expect(
        closeRel(aL, golden.q_left_edge[k], tol),
        `q_left_edge[i=${i}]: got ${aL}, golden ${golden.q_left_edge[k]}`,
      ).toBe(true);
      expect(
        closeRel(aR, golden.q_right_edge[k], tol),
        `q_right_edge[i=${i}]: got ${aR}, golden ${golden.q_right_edge[k]}`,
      ).toBe(true);
    }
  });

  it("matches the golden q_parabola samples at every (cell, xi)", () => {
    for (let k = 0; k < fixture.query_cells.length; k++) {
      const i = fixture.query_cells[k];
      for (let s = 0; s < fixture.parabola_xi.length; s++) {
        const xi = fixture.parabola_xi[s];
        const got = rules.reconstructOutput(
          rule,
          "q_parabola",
          cellAverages,
          i,
          { xi },
        );
        const want = golden.q_parabola[k][s];
        expect(
          closeRel(got, want, tol),
          `q_parabola(i=${i}, xi=${xi}): got ${got}, golden ${want}`,
        ).toBe(true);
      }
    }
  });

  it("reconstructCell agrees with applyEdgeStencil + evaluateParabola", () => {
    for (let k = 0; k < fixture.query_cells.length; k++) {
      const i = fixture.query_cells[k];
      const cell = rules.reconstructCell(rule, cellAverages, i);
      const aL = rules.applyEdgeStencil(rule, "q_left_edge", cellAverages, i);
      const aR = rules.applyEdgeStencil(rule, "q_right_edge", cellAverages, i);
      expect(cell.q_left_edge).toBe(aL);
      expect(cell.q_right_edge).toBe(aR);
      for (const xi of fixture.parabola_xi) {
        const direct = rules.evaluateParabola(aL, aR, cellAverages[i], xi);
        expect(cell.parabola(xi)).toBe(direct);
      }
    }
  });
});
