import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import {
  mpasMeshData,
  type MPASMeshData,
  type MpasMetricName,
} from "../src/grids/index.js";

// Harness lives at <repo>/tests/conformance/grids/mpas — the TS package
// root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "mpas",
);

interface MeshArrays {
  max_edges: number;
  n_vertices: number;
  lon_cell: number[];
  lat_cell: number[];
  x_cell: number[];
  y_cell: number[];
  z_cell: number[];
  area_cell: number[];
  n_edges_on_cell: number[];
  cells_on_cell: number[];
  edges_on_cell: number[];
  lon_edge: number[];
  lat_edge: number[];
  cells_on_edge: number[];
  dc_edge: number[];
  dv_edge: number[];
}

interface Fixture {
  name: string;
  opts: { R: number; ghosts: number; dtype: "float64" | "float32" };
  mesh: MeshArrays;
  query_cells: number[];
  query_edges: number[];
}

interface FixtureSpec {
  tolerance: { relative: number };
  fixtures: Fixture[];
}

interface Golden {
  n_cells: number;
  n_edges: number;
  max_edges: number;
  cell_centers: Array<{ lon: number; lat: number }>;
  cell_centers_cart: Array<{ x: number; y: number; z: number }>;
  neighbors: number[][];
  cell_area: number[];
  cell_metrics: Record<string, number[]>;
  edge_length: number[];
  cell_distance: number[];
  edge_metrics: Record<string, number[]>;
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);
const REL_TOL = spec.tolerance.relative;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function buildMesh(fx: Fixture): MPASMeshData {
  return mpasMeshData({
    lon_cell: fx.mesh.lon_cell,
    lat_cell: fx.mesh.lat_cell,
    area_cell: fx.mesh.area_cell,
    n_edges_on_cell: fx.mesh.n_edges_on_cell,
    cells_on_cell: fx.mesh.cells_on_cell,
    edges_on_cell: fx.mesh.edges_on_cell,
    lon_edge: fx.mesh.lon_edge,
    lat_edge: fx.mesh.lat_edge,
    cells_on_edge: fx.mesh.cells_on_edge,
    dc_edge: fx.mesh.dc_edge,
    dv_edge: fx.mesh.dv_edge,
    max_edges: fx.mesh.max_edges,
    x_cell: fx.mesh.x_cell,
    y_cell: fx.mesh.y_cell,
    z_cell: fx.mesh.z_cell,
    n_vertices: fx.mesh.n_vertices,
    R: fx.opts.R,
  });
}

const CELL_METRICS: MpasMetricName[] = [
  "lon",
  "lat",
  "area",
  "x",
  "y",
  "z",
  "n_edges_on_cell",
];
const EDGE_METRICS: MpasMetricName[] = [
  "lon_edge",
  "lat_edge",
  "dc_edge",
  "dv_edge",
];

describe("mpas cross-language conformance", () => {
  for (const fixture of spec.fixtures) {
    it(`matches the reference golden (fixture=${fixture.name})`, () => {
      const grid = grids.mpas({
        mesh: buildMesh(fixture),
        R: fixture.opts.R,
        ghosts: fixture.opts.ghosts,
        dtype: fixture.opts.dtype,
      });

      const golden: Golden = JSON.parse(
        readFileSync(join(HARNESS_DIR, "golden", `${fixture.name}.json`), "utf8"),
      );
      expect(grid.n_cells).toBe(golden.n_cells);
      expect(grid.n_edges).toBe(golden.n_edges);
      expect(grid.max_edges).toBe(golden.max_edges);

      fixture.query_cells.forEach((c, k) => {
        const { lon, lat } = grid.cell_centers(c);
        expect(closeRel(lon, golden.cell_centers[k].lon)).toBe(true);
        expect(closeRel(lat, golden.cell_centers[k].lat)).toBe(true);

        const { x, y, z } = grid.cell_center_cart(c);
        expect(closeRel(x, golden.cell_centers_cart[k].x)).toBe(true);
        expect(closeRel(y, golden.cell_centers_cart[k].y)).toBe(true);
        expect(closeRel(z, golden.cell_centers_cart[k].z)).toBe(true);

        expect(grid.neighbors(c)).toEqual(golden.neighbors[k]);
        expect(closeRel(grid.cell_area(c), golden.cell_area[k])).toBe(true);

        for (const name of CELL_METRICS) {
          const got = grid.metric_eval(name, c);
          const expected = golden.cell_metrics[name][k];
          expect(closeRel(got, expected)).toBe(true);
        }
      });

      fixture.query_edges.forEach((e, k) => {
        expect(closeRel(grid.edge_length(e), golden.edge_length[k])).toBe(true);
        expect(closeRel(grid.cell_distance(e), golden.cell_distance[k])).toBe(true);
        for (const name of EDGE_METRICS) {
          const got = grid.metric_eval(name, e);
          const expected = golden.edge_metrics[name][k];
          expect(closeRel(got, expected)).toBe(true);
        }
      });
    });
  }
});
