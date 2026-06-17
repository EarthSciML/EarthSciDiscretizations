/**
 * Cross-binding conformance for the centered_2nd_uniform_latlon rule
 * evaluator.
 *
 * Loads the on-disk rule from `discretizations/finite_difference/` and the
 * shared fixture + golden under `tests/conformance/grids/latlon/rules/`, then
 * verifies the TypeScript AST evaluator reproduces the binding-independent
 * golden coefficients to within the fixture-declared tolerance.
 *
 * The Python (dsc-ve1) and Rust (dsc-k86) sibling evaluators consume the
 * same fixture + golden, so passing tests across all three bindings imply
 * cross-binding agreement.
 */

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids, rules } from "../src/index.js";
import type { ExpressionNode, OpNode } from "../src/rules/types.js";
import type { LatLonGrid, LatLonOpts } from "../src/grids/lat_lon.js";

const REPO_ROOT = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
);
const RULE_PATH = join(
  REPO_ROOT,
  "discretizations",
  "finite_difference",
  "centered_2nd_uniform_latlon.json",
);
const FIXTURE_DIR = join(
  REPO_ROOT,
  "tests",
  "conformance",
  "grids",
  "latlon",
  "rules",
  "centered_2nd_uniform_latlon",
);

interface GridSpec {
  variant: "regular";
  nlon: number;
  nlat: number;
  R: number;
  ghosts: number;
  dtype: "float64" | "float32";
  pole_policy: "none";
}

interface StencilKey {
  axis: "lon" | "lat";
  offset: number;
}

interface FixtureSpec {
  rule: string;
  tolerance: { relative: number };
  grid: GridSpec;
  query_points: Array<[number, number]>;
  stencil_order: StencilKey[];
}

interface Golden {
  rule: string;
  grid: GridSpec;
  query_points: Array<[number, number]>;
  cell_centers: Array<{ lon: number; lat: number }>;
  cos_lat: number[];
  coeffs: Record<string, number[]>;
}

function closeRel(a: number, b: number, tol: number): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function buildGrid(spec: GridSpec): LatLonGrid {
  const opts: LatLonOpts = {
    variant: spec.variant,
    nlon: spec.nlon,
    nlat: spec.nlat,
    R: spec.R,
    ghosts: spec.ghosts,
    dtype: spec.dtype,
    pole_policy: spec.pole_policy,
  };
  return grids.lat_lon(opts);
}

function tagFor(key: StencilKey): string {
  return `${key.axis}_${key.offset >= 0 ? "p" : "m"}${Math.abs(key.offset)}`;
}

// esd-t4h: the rule no longer carries a `stencil` list. centered_2nd_uniform_latlon
// lowers grad(u) to a closed `replacement` op-tree:
//   { op: "+", args: [term, term, term, term] }
//   term = { op: "*", args: [<coeff-AST>, { op: "index",
//            args: ["$u", <lat-pos>, <lon-pos>] }] }
// Each <pos> is either a bare axis name ("lat"/"lon") or an offset node
// { op: "+", args: [<axis>, <int>] }. We reconstruct the directional
// { axis, offset, coeff } triples from that AST — the real current contract.

interface StencilTerm {
  axis: "lon" | "lat";
  offset: number;
  coeff: ExpressionNode;
}

function isOpNode(node: ExpressionNode): node is OpNode {
  return typeof node === "object" && node !== null && "op" in node;
}

/** Parse an index position into `{axis, offset}`, or `null` for a bare axis. */
function axisOffset(pos: ExpressionNode): StencilTerm | null {
  if (typeof pos === "string") return null;
  if (!isOpNode(pos) || pos.op !== "+") {
    throw new Error(`unrecognized index position node: ${JSON.stringify(pos)}`);
  }
  const axis = pos.args.find(
    (a): a is "lon" | "lat" => a === "lon" || a === "lat",
  );
  const offset = pos.args.find((a): a is number => typeof a === "number");
  if (axis === undefined || offset === undefined || offset === 0) {
    throw new Error(`malformed offset node: ${JSON.stringify(pos)}`);
  }
  return { axis, offset, coeff: 0 };
}

/** Extract the directional `{axis, offset, coeff}` terms from `replacement`. */
function extractTerms(rule: rules.Rule): StencilTerm[] {
  const repl = rule.replacement;
  if (repl === undefined || !isOpNode(repl) || repl.op !== "+") {
    throw new Error(
      `centered_2nd_uniform_latlon expected a sum replacement op-tree, got ${JSON.stringify(repl)}`,
    );
  }
  return repl.args.map((term) => {
    if (!isOpNode(term) || term.op !== "*") {
      throw new Error(`expected a scaled index term, got ${JSON.stringify(term)}`);
    }
    const [coeff, idx] = term.args;
    if (!isOpNode(idx) || idx.op !== "index") {
      throw new Error(`expected an index op as the second factor, got ${JSON.stringify(idx)}`);
    }
    // index args: [field, lat-pos, lon-pos]
    const differenced = idx.args
      .slice(1)
      .map(axisOffset)
      .filter((d): d is StencilTerm => d !== null);
    if (differenced.length !== 1) {
      throw new Error(
        `each term must difference exactly one axis, got ${differenced.length}`,
      );
    }
    return { axis: differenced[0].axis, offset: differenced[0].offset, coeff };
  });
}

function findTerm(
  terms: StencilTerm[],
  axis: string,
  offset: number,
): StencilTerm {
  const found = terms.find((t) => t.axis === axis && t.offset === offset);
  if (!found) {
    throw new Error(`replacement term for axis=${axis} offset=${offset} not found`);
  }
  return found;
}

describe("centered_2nd_uniform_latlon rule eval cross-binding conformance", () => {
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
  const grid = buildGrid(fixture.grid);

  it("agrees with golden cell_centers and cos_lat at every query point", () => {
    expect(fixture.query_points.length).toBe(golden.query_points.length);
    for (let k = 0; k < fixture.query_points.length; k++) {
      const [j, i] = fixture.query_points[k];
      const { lon, lat } = grid.cell_center(j, i);
      expect(closeRel(lon, golden.cell_centers[k].lon, tol)).toBe(true);
      expect(closeRel(lat, golden.cell_centers[k].lat, tol)).toBe(true);
      expect(closeRel(Math.cos(lat), golden.cos_lat[k], tol)).toBe(true);
    }
  });

  it("evaluates each directional coefficient to match the golden", () => {
    const dlon = (2 * Math.PI) / fixture.grid.nlon;
    const dlat = Math.PI / fixture.grid.nlat;

    const terms = extractTerms(rule);
    for (const key of fixture.stencil_order) {
      const entry = findTerm(terms, key.axis, key.offset);
      const tag = tagFor(key);
      const expected = golden.coeffs[tag];
      expect(
        expected,
        `golden coeffs missing tag ${tag}`,
      ).toBeDefined();

      for (let k = 0; k < fixture.query_points.length; k++) {
        const [j, i] = fixture.query_points[k];
        const { lat } = grid.cell_center(j, i);
        const bindings = new Map<string, number>([
          ["R", fixture.grid.R],
          ["dlon", dlon],
          ["dlat", dlat],
          ["cos_lat", Math.cos(lat)],
        ]);
        const got = rules.evaluate(entry.coeff, bindings);
        expect(
          closeRel(got, expected[k], tol),
          `coeff ${tag} at qp[${k}]=(${j},${i}): got ${got}, expected ${expected[k]}`,
        ).toBe(true);
      }
    }
  });

  it("rule structural shape matches the migrated (esd-t4h) contract", () => {
    expect(rule.family).toBe("finite_difference");
    expect(rule.grid_family).toBe("latlon");
    expect(rule.combine).toBe("+");
    // Migrated off `stencil` to a closed `replacement` op-tree that sums four
    // directional terms (lon±1, lat±1).
    expect(rule.stencil).toBeUndefined();
    expect(rule.replacement).toBeDefined();
    const terms = extractTerms(rule);
    expect(terms.length).toBe(4);
    const axes = new Set(terms.map((t) => t.axis));
    expect(axes.has("lon")).toBe(true);
    expect(axes.has("lat")).toBe(true);
    const offsets = new Set(terms.map((t) => t.offset));
    expect(offsets.has(-1)).toBe(true);
    expect(offsets.has(1)).toBe(true);
  });
});
