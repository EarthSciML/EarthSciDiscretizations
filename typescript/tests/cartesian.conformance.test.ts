import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type {
  AxisExtent,
  CartesianMetricName,
  CartesianOpts,
} from "../src/grids/cartesian.js";

// Harness lives at <repo>/tests/conformance/grids/cartesian — the TS
// package root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "cartesian",
);

interface FixtureOpts {
  nx?: number;
  ny?: number;
  nz?: number;
  extent?: Array<[number, number]>;
  edges?: number[][];
  ghosts?: number;
  dtype?: "float64" | "float32";
}

interface Fixture {
  name: string;
  opts: FixtureOpts;
  query_points: number[][];
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

interface NeighborEntry {
  axis: number;
  side: number;
  index: number[];
}

interface Golden {
  ndim: number;
  n_cells: number;
  cell_centers: number[][];
  cell_widths: number[][];
  cell_volume: number[];
  neighbors: NeighborEntry[][];
  metric_g: number[][][];
  metric_volume: number[];
  metric_jacobian: number[];
  [key: string]: unknown;
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);
const REL_TOL = spec.tolerance.relative;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function buildOpts(opts: FixtureOpts): CartesianOpts {
  const out: CartesianOpts = {};
  if (opts.dtype !== undefined) out.dtype = opts.dtype;
  if (opts.ghosts !== undefined) out.ghosts = opts.ghosts;
  if (opts.edges !== undefined) {
    out.edges = opts.edges.map((e) => e.slice());
  } else {
    if (opts.nx !== undefined) out.nx = opts.nx;
    if (opts.ny !== undefined) out.ny = opts.ny;
    if (opts.nz !== undefined) out.nz = opts.nz;
    if (opts.extent !== undefined) {
      out.extent = opts.extent.map(
        (e) => [e[0], e[1]] as unknown as AxisExtent,
      );
    }
  }
  return out;
}

const WIDTH_NAMES = ["dx", "dy", "dz"] as const satisfies readonly CartesianMetricName[];
const FACE_NAMES = [
  "face_area_x",
  "face_area_y",
  "face_area_z",
] as const satisfies readonly CartesianMetricName[];

describe("cartesian cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the Julia-reference golden (fixture=${fixture.name})`, () => {
      const grid = grids.cartesian(buildOpts(fixture.opts));
      const golden: Golden = JSON.parse(
        readFileSync(
          join(HARNESS_DIR, "golden", `${fixture.name}.json`),
          "utf8",
        ),
      );
      const ndim = grid.ndim;
      expect(ndim).toBe(golden.ndim);
      expect(grid.n_cells).toBe(golden.n_cells);
      expect(fixture.query_points.length).toBe(golden.cell_centers.length);

      for (let k = 0; k < fixture.query_points.length; k++) {
        const qp = fixture.query_points[k];

        const c = grid.cell_centers(...qp);
        const gc = golden.cell_centers[k];
        for (let d = 0; d < ndim; d++) {
          expect(
            closeRel(c[d], gc[d]),
            `cell_centers[${d}] at qp[${k}]=${JSON.stringify(qp)}`,
          ).toBe(true);
        }

        const w = grid.cell_widths(...qp);
        const gw = golden.cell_widths[k];
        for (let d = 0; d < ndim; d++) {
          expect(
            closeRel(w[d], gw[d]),
            `cell_widths[${d}] at qp[${k}]=${JSON.stringify(qp)}`,
          ).toBe(true);
        }

        const v = grid.cell_volume(...qp);
        expect(
          closeRel(v, golden.cell_volume[k]),
          `cell_volume at qp[${k}]=${JSON.stringify(qp)}`,
        ).toBe(true);

        // Neighbors: TS returns a list in axis-major / low-before-high
        // order with 0-based axis indices — identical to the golden shape.
        const nbrs = grid.neighbors(...qp);
        const gnbrs = golden.neighbors[k];
        expect(nbrs.length).toBe(gnbrs.length);
        for (let n = 0; n < nbrs.length; n++) {
          expect(nbrs[n].axis).toBe(gnbrs[n].axis);
          expect(nbrs[n].side).toBe(gnbrs[n].side);
          expect(nbrs[n].index).toEqual(gnbrs[n].index);
        }

        // Scalar metrics
        expect(
          closeRel(
            grid.metric_eval("volume", ...qp) as number,
            golden.metric_volume[k],
          ),
        ).toBe(true);
        expect(
          closeRel(
            grid.metric_eval("jacobian", ...qp) as number,
            golden.metric_jacobian[k],
          ),
        ).toBe(true);
        for (let d = 0; d < ndim; d++) {
          const wGot = grid.metric_eval(WIDTH_NAMES[d], ...qp) as number;
          const wExp = golden[`metric_${WIDTH_NAMES[d]}`] as number[];
          expect(
            closeRel(wGot, wExp[k]),
            `metric ${WIDTH_NAMES[d]} at qp[${k}]=${JSON.stringify(qp)}`,
          ).toBe(true);
          const fGot = grid.metric_eval(FACE_NAMES[d], ...qp) as number;
          const fExp = golden[`metric_${FACE_NAMES[d]}`] as number[];
          expect(
            closeRel(fGot, fExp[k]),
            `metric ${FACE_NAMES[d]} at qp[${k}]=${JSON.stringify(qp)}`,
          ).toBe(true);
        }

        const g = grid.metric_eval("g", ...qp) as number[][];
        const gg = golden.metric_g[k];
        for (let i = 0; i < ndim; i++) {
          for (let j = 0; j < ndim; j++) {
            expect(
              closeRel(g[i][j], gg[i][j]),
              `metric_g[${i}][${j}] at qp[${k}]=${JSON.stringify(qp)}`,
            ).toBe(true);
          }
        }
      }
    });
  }
});
