import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type { MetricName } from "../src/grids/cubed_sphere.js";

// Harness lives at <repo>/tests/conformance/grids/cubed_sphere — the TS
// package root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "cubed_sphere",
);

interface Fixture {
  name: string;
  opts: { Nc: number; R: number; ghosts: number; dtype: "float64" | "float32" };
  query_points: Array<[number, number, number]>;
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

interface Golden {
  n_cells: number;
  cell_centers: Array<{ lon: number; lat: number }>;
  neighbors_W: Array<[number, number, number]>;
  neighbors_E: Array<[number, number, number]>;
  neighbors_S: Array<[number, number, number]>;
  neighbors_N: Array<[number, number, number]>;
  metric_J: number[];
  metric_g_xixi: number[];
  metric_g_etaeta: number[];
  metric_g_xieta: number[];
  metric_ginv_xixi: number[];
  metric_ginv_etaeta: number[];
  metric_ginv_xieta: number[];
  area: number[];
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);
const REL_TOL = spec.tolerance.relative;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

const METRIC_NAMES: MetricName[] = [
  "J",
  "g_xixi",
  "g_etaeta",
  "g_xieta",
  "ginv_xixi",
  "ginv_etaeta",
  "ginv_xieta",
];

type Dir = "W" | "E" | "S" | "N";
const DIR_EDGE: Record<Dir, "west" | "east" | "south" | "north"> = {
  W: "west",
  E: "east",
  S: "south",
  N: "north",
};
const DIR_ORDER: Dir[] = ["W", "E", "S", "N"];

describe("cubed_sphere cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the Julia-reference golden (fixture=${fixture.name})`, () => {
      const g = grids.cubed_sphere({
        Nc: fixture.opts.Nc,
        R: fixture.opts.R,
        ghosts: fixture.opts.ghosts,
        dtype: fixture.opts.dtype,
      });
      const golden: Golden = JSON.parse(
        readFileSync(join(HARNESS_DIR, "golden", `${fixture.name}.json`), "utf8"),
      );
      expect(g.n_cells).toBe(golden.n_cells);
      expect(fixture.query_points.length).toBe(golden.cell_centers.length);

      for (let k = 0; k < fixture.query_points.length; k++) {
        const [p, i, j] = fixture.query_points[k];

        const { lon, lat } = g.cell_centers(p, i, j);
        expect(closeRel(lon, golden.cell_centers[k].lon)).toBe(true);
        expect(closeRel(lat, golden.cell_centers[k].lat)).toBe(true);

        const nbrs = g.neighbors(p, i, j);
        for (const dir of DIR_ORDER) {
          const edge = DIR_EDGE[dir];
          const ref = nbrs.find((r) => r.edge === edge);
          expect(ref, `neighbor ${dir} missing at qp[${k}]`).toBeDefined();
          const expected = golden[`neighbors_${dir}` as const][k];
          expect([ref!.panel, ref!.i, ref!.j]).toEqual(expected);
        }

        for (const name of METRIC_NAMES) {
          const got = g.metric_eval(name, p, i, j);
          const expected = golden[`metric_${name}` as const][k];
          expect(
            closeRel(got, expected),
            `metric ${name} at qp[${k}]=(${p},${i},${j}): got ${got}, expected ${expected}`,
          ).toBe(true);
        }

        const gotArea = g.cell_area(p, i, j);
        expect(
          closeRel(gotArea, golden.area[k]),
          `area at qp[${k}]=(${p},${i},${j}): got ${gotArea}, expected ${golden.area[k]}`,
        ).toBe(true);
      }
    });
  }
});
