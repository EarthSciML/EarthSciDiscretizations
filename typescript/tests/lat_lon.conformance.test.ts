import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type {
  LatLonGrid,
  LatLonMetricName,
  LatLonVariant,
  PolePolicy,
} from "../src/grids/lat_lon.js";

// Harness lives at <repo>/tests/conformance/grids/latlon — the TS package
// root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "latlon",
);

interface FixtureOpts {
  variant: LatLonVariant;
  nlon?: number;
  nlat?: number;
  nlon_per_row?: number[];
  lat_edges?: number[];
  R: number;
  ghosts: number;
  dtype: "float64" | "float32";
  pole_policy: PolePolicy;
}

interface Fixture {
  name: string;
  opts: FixtureOpts;
  query_points: Array<[number, number]>;
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

type GoldenNeighbor = [number, number] | null;

interface Golden {
  n_cells: number;
  cell_centers: Array<{ lon: number; lat: number }>;
  neighbors_W: GoldenNeighbor[];
  neighbors_E: GoldenNeighbor[];
  neighbors_S: GoldenNeighbor[];
  neighbors_N: GoldenNeighbor[];
  metric_J: number[];
  metric_g_lonlon: number[];
  metric_g_latlat: number[];
  metric_g_lonlat: number[];
  metric_ginv_lonlon: number[];
  metric_ginv_latlat: number[];
  metric_ginv_lonlat: number[];
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

const METRIC_NAMES: LatLonMetricName[] = [
  "J",
  "g_lonlon",
  "g_latlat",
  "g_lonlat",
  "ginv_lonlon",
  "ginv_latlat",
  "ginv_lonlat",
];

type Dir = "W" | "E" | "S" | "N";
const DIR_ORDER: Dir[] = ["W", "E", "S", "N"];

function buildGrid(opts: FixtureOpts): LatLonGrid {
  if (opts.variant === "regular") {
    return grids.lat_lon({
      variant: "regular",
      nlon: opts.nlon!,
      nlat: opts.nlat!,
      R: opts.R,
      ghosts: opts.ghosts,
      dtype: opts.dtype,
      pole_policy: opts.pole_policy,
    });
  }
  return grids.lat_lon({
    variant: "reduced_gaussian",
    nlon_per_row: opts.nlon_per_row!,
    lat_edges: opts.lat_edges,
    R: opts.R,
    ghosts: opts.ghosts,
    dtype: opts.dtype,
    pole_policy: opts.pole_policy,
  });
}

describe("lat_lon cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the reference golden (fixture=${fixture.name})`, () => {
      const g = buildGrid(fixture.opts);
      const golden: Golden = JSON.parse(
        readFileSync(join(HARNESS_DIR, "golden", `${fixture.name}.json`), "utf8"),
      );
      expect(g.n_cells).toBe(golden.n_cells);
      expect(fixture.query_points.length).toBe(golden.cell_centers.length);

      for (let k = 0; k < fixture.query_points.length; k++) {
        const [j, i] = fixture.query_points[k];

        const { lon, lat } = g.cell_center(j, i);
        expect(
          closeRel(lon, golden.cell_centers[k].lon),
          `lon at qp[${k}]=(${j},${i}): got ${lon}, expected ${golden.cell_centers[k].lon}`,
        ).toBe(true);
        expect(
          closeRel(lat, golden.cell_centers[k].lat),
          `lat at qp[${k}]=(${j},${i}): got ${lat}, expected ${golden.cell_centers[k].lat}`,
        ).toBe(true);

        const nbrs = g.neighbors(j, i);
        for (const dir of DIR_ORDER) {
          const got = nbrs[dir];
          const expected = golden[`neighbors_${dir}` as const][k];
          if (expected === null) {
            expect(
              got,
              `neighbor ${dir} at qp[${k}]=(${j},${i}): expected pole (null), got ${JSON.stringify(got)}`,
            ).toBeNull();
          } else {
            expect(
              got,
              `neighbor ${dir} at qp[${k}]=(${j},${i}): expected ${JSON.stringify(expected)}, got null`,
            ).not.toBeNull();
            expect([got!.j, got!.i]).toEqual(expected);
          }
        }

        for (const name of METRIC_NAMES) {
          const got = g.metric_eval(name, j, i);
          const expected = golden[`metric_${name}` as const][k];
          expect(
            closeRel(got, expected),
            `metric ${name} at qp[${k}]=(${j},${i}): got ${got}, expected ${expected}`,
          ).toBe(true);
        }

        const gotArea = g.cell_area(j, i);
        expect(
          closeRel(gotArea, golden.area[k]),
          `area at qp[${k}]=(${j},${i}): got ${gotArea}, expected ${golden.area[k]}`,
        ).toBe(true);
      }
    });
  }
});
