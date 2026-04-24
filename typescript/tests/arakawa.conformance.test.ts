import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type {
  ArakawaGrid,
  ArakawaLocation,
  ArakawaMetricName,
  ArakawaStagger,
  ArakawaVariable,
} from "../src/grids/arakawa.js";

// Harness lives at <repo>/tests/conformance/grids/arakawa — the TS package
// root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "arakawa",
);

interface BaseOpts {
  family: "cartesian";
  xlo: number;
  xhi: number;
  ylo: number;
  yhi: number;
  nx: number;
  ny: number;
}

interface FixtureOpts {
  base: BaseOpts;
  dtype: "float64" | "float32";
  ghosts: number;
}

interface Fixture {
  name: string;
  opts: FixtureOpts;
  staggers: ArakawaStagger[];
  query_points: Record<ArakawaLocation, Array<[number, number]>>;
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

type NeighborWire = [number, number] | null;

interface PerLocCoords {
  points: Array<[number, number]>;
  xy: Array<[number, number]>;
}

interface PerLocNeighbors {
  points: Array<[number, number]>;
  W: NeighborWire[];
  E: NeighborWire[];
  S: NeighborWire[];
  N: NeighborWire[];
}

interface StaggerTable {
  rotated: boolean;
  variable_locations: Record<ArakawaVariable, ArakawaLocation>;
  variable_shapes: Record<ArakawaVariable, [number, number]>;
  location_shapes: Record<ArakawaLocation, [number, number]>;
}

interface Golden {
  n_cells: number;
  coords: Record<ArakawaLocation, PerLocCoords>;
  neighbors: Record<ArakawaLocation, PerLocNeighbors>;
  metrics: {
    points: Array<[number, number]>;
    dx: number[];
    dy: number[];
    area: number[];
  };
  staggers: Record<ArakawaStagger, StaggerTable>;
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);
const REL_TOL = spec.tolerance.relative;

const LOCATIONS: ArakawaLocation[] = [
  "cell_center",
  "u_edge",
  "v_edge",
  "corner",
];
const VARIABLES: ArakawaVariable[] = ["h", "u", "v"];
const METRIC_NAMES: ArakawaMetricName[] = ["dx", "dy", "area"];
type Dir = "W" | "E" | "S" | "N";
const DIR_ORDER: Dir[] = ["W", "E", "S", "N"];
const DIR_EDGE: Record<Dir, "west" | "east" | "south" | "north"> = {
  W: "west",
  E: "east",
  S: "south",
  N: "north",
};

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function buildGrid(
  base: BaseOpts,
  stagger: ArakawaStagger,
  dtype: "float64" | "float32",
  ghosts: number,
): ArakawaGrid {
  if (base.family !== "cartesian") {
    throw new Error(`only cartesian base supported; got ${base.family}`);
  }
  const b = grids.cartesianBase({
    xlo: base.xlo,
    xhi: base.xhi,
    ylo: base.ylo,
    yhi: base.yhi,
    nx: base.nx,
    ny: base.ny,
  });
  return grids.arakawa({ base: b, stagger, dtype, ghosts });
}

function coordAt(
  gA: ArakawaGrid,
  gC: ArakawaGrid,
  loc: ArakawaLocation,
  i: number,
  j: number,
): [number, number] {
  switch (loc) {
    case "cell_center":
      return gA.cell_centers(i, j);
    case "corner":
      return gA.corners(i, j);
    case "u_edge":
      return gC.u_face(i, j);
    case "v_edge":
      return gC.v_face(i, j);
  }
}

describe("arakawa cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the Julia-reference golden (fixture=${fixture.name})`, () => {
      const opts = fixture.opts;
      const golden: Golden = JSON.parse(
        readFileSync(
          join(HARNESS_DIR, "golden", `${fixture.name}.json`),
          "utf8",
        ),
      );

      const gA = buildGrid(opts.base, "A", opts.dtype, opts.ghosts);
      const gC = buildGrid(opts.base, "C", opts.dtype, opts.ghosts);
      expect(gA.n_cells).toBe(golden.n_cells);

      // --- coord + neighbor tables (stagger-independent) ----------------
      for (const loc of LOCATIONS) {
        const coords = golden.coords[loc];
        const nbrs = golden.neighbors[loc];
        for (let k = 0; k < coords.points.length; k++) {
          const [i, j] = coords.points[k];
          const [x, y] = coordAt(gA, gC, loc, i, j);
          const [exX, exY] = coords.xy[k];
          expect(
            closeRel(x, exX),
            `${fixture.name}: ${loc} x mismatch at qp[${k}]=(${i},${j}): ${x} vs ${exX}`,
          ).toBe(true);
          expect(
            closeRel(y, exY),
            `${fixture.name}: ${loc} y mismatch at qp[${k}]=(${i},${j}): ${y} vs ${exY}`,
          ).toBe(true);

          const got = gA.neighbors(loc, i, j);
          for (const dkey of DIR_ORDER) {
            const expected = nbrs[dkey][k];
            const gotEdge = got[DIR_EDGE[dkey]];
            if (expected === null) {
              expect(
                gotEdge,
                `${fixture.name}: neighbor ${dkey} ${loc} at qp[${k}]=(${i},${j}) should be null`,
              ).toBeNull();
            } else {
              expect(
                gotEdge,
                `${fixture.name}: neighbor ${dkey} ${loc} at qp[${k}]=(${i},${j})`,
              ).toEqual([expected[0], expected[1]]);
            }
          }
        }
      }

      // --- metric_eval (stagger-independent for cartesian) --------------
      const metrics = golden.metrics;
      for (let k = 0; k < metrics.points.length; k++) {
        const [i, j] = metrics.points[k];
        for (const mname of METRIC_NAMES) {
          const got = gA.metric_eval(mname, i, j);
          const expected = metrics[mname][k];
          expect(
            closeRel(got, expected),
            `${fixture.name}: metric ${mname} mismatch at qp[${k}]=(${i},${j}): ${got} vs ${expected}`,
          ).toBe(true);
        }
      }

      // --- per-stagger variable locations + shapes ----------------------
      for (const sname of fixture.staggers) {
        const gS = buildGrid(opts.base, sname, opts.dtype, opts.ghosts);
        const stab = golden.staggers[sname];

        expect(gS.rotated, `stagger=${sname} rotated`).toBe(stab.rotated);
        for (const v of VARIABLES) {
          expect(
            gS.variable_location(v),
            `stagger=${sname} var=${v} location`,
          ).toBe(stab.variable_locations[v]);
          const got = gS.variable_shape(v);
          expect(
            [got[0], got[1]],
            `stagger=${sname} var=${v} shape`,
          ).toEqual([stab.variable_shapes[v][0], stab.variable_shapes[v][1]]);
        }
        for (const loc of LOCATIONS) {
          const got = gS.location_shape(loc);
          expect(
            [got[0], got[1]],
            `stagger=${sname} loc=${loc} shape`,
          ).toEqual([stab.location_shapes[loc][0], stab.location_shapes[loc][1]]);
        }
      }
    });
  }
});
