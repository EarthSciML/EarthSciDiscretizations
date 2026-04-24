import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type { DuoMetricName } from "../src/grids/duo.js";

// Harness lives at <repo>/tests/conformance/grids/duo — the TS package
// root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "duo",
);

interface Fixture {
  name: string;
  opts: { level: number; R: number; ghosts: number; dtype: "float64" | "float32" };
  query_points: number[];
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

interface Golden {
  n_cells: number;
  n_vertices: number;
  n_edges: number;
  cell_centers: Array<{ lon: number; lat: number }>;
  neighbors: Array<[number, number, number]>;
  metric_area: number[];
  metric_lon: number[];
  metric_lat: number[];
  metric_x: number[];
  metric_y: number[];
  metric_z: number[];
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);
const REL_TOL = spec.tolerance.relative;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

const METRIC_NAMES: DuoMetricName[] = ["area", "lon", "lat", "x", "y", "z"];

describe("duo cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the Julia-reference golden (fixture=${fixture.name})`, () => {
      const g = grids.duo({
        loader: { path: `builtin://icosahedral/${fixture.opts.level}` },
        R: fixture.opts.R,
        ghosts: fixture.opts.ghosts,
        dtype: fixture.opts.dtype,
      });
      const golden: Golden = JSON.parse(
        readFileSync(join(HARNESS_DIR, "golden", `${fixture.name}.json`), "utf8"),
      );

      expect(g.n_cells).toBe(golden.n_cells);
      expect(g.n_vertices).toBe(golden.n_vertices);
      expect(g.n_edges).toBe(golden.n_edges);
      expect(fixture.query_points.length).toBe(golden.cell_centers.length);

      for (let k = 0; k < fixture.query_points.length; k++) {
        const c = fixture.query_points[k];

        const { lon, lat } = g.cell_centers(c);
        expect(
          closeRel(lon, golden.cell_centers[k].lon),
          `lon at qp[${k}]=c${c}: got ${lon}, expected ${golden.cell_centers[k].lon}`,
        ).toBe(true);
        expect(
          closeRel(lat, golden.cell_centers[k].lat),
          `lat at qp[${k}]=c${c}: got ${lat}, expected ${golden.cell_centers[k].lat}`,
        ).toBe(true);

        const nbrs = g.neighbors(c);
        expect([nbrs[0], nbrs[1], nbrs[2]]).toEqual(golden.neighbors[k]);

        for (const name of METRIC_NAMES) {
          const got = g.metric_eval(name, c);
          const expected = golden[`metric_${name}` as const][k];
          expect(
            closeRel(got, expected),
            `metric ${name} at qp[${k}]=c${c}: got ${got}, expected ${expected}`,
          ).toBe(true);
        }
      }
    });
  }
});
