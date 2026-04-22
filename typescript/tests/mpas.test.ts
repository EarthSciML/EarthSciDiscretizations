import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import {
  mpasMeshData,
  checkMesh,
  type MPASMeshData,
  type MPASGrid,
} from "../src/grids/index.js";

/**
 * Tiny synthetic global MPAS-like mesh: 4 triangular Voronoi cells at
 * tetrahedral vertex directions. Each cell has 3 neighbors (the other
 * three cells); total area = 4πR²; per-cell area = πR².
 *
 * Mirrors `test/test_mpas_grid.jl::_tetra_mesh` so the pinned query
 * points match the Julia reference binding's expected values.
 */
function tetraMesh(R = 1.0): MPASMeshData {
  const verts: Array<[number, number, number]> = [
    [1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
  ];
  const n_cells = 4;
  const max_edges = 3;
  const lon_cell = new Array<number>(n_cells);
  const lat_cell = new Array<number>(n_cells);
  const x_cell = new Array<number>(n_cells);
  const y_cell = new Array<number>(n_cells);
  const z_cell = new Array<number>(n_cells);
  for (let c = 0; c < n_cells; c++) {
    const [x, y, z] = verts[c];
    const n = Math.sqrt(x * x + y * y + z * z);
    const xh = x / n;
    const yh = y / n;
    const zh = z / n;
    lon_cell[c] = Math.atan2(yh, xh);
    lat_cell[c] = Math.asin(zh);
    x_cell[c] = R * xh;
    y_cell[c] = R * yh;
    z_cell[c] = R * zh;
  }
  const area_cell = Array<number>(n_cells).fill(Math.PI * R * R);
  const n_edges_on_cell = Array<number>(n_cells).fill(3);
  // cells_on_cell row-major (n_cells × max_edges): all other cells.
  const cells_on_cell = new Array<number>(n_cells * max_edges);
  for (let c = 0; c < n_cells; c++) {
    let j = 0;
    for (let k = 0; k < n_cells; k++) {
      if (k === c) continue;
      cells_on_cell[c * max_edges + j] = k;
      j++;
    }
  }
  // Edges: one per unordered cell pair = 6 edges.
  const n_edges = 6;
  const cells_on_edge = new Array<number>(n_edges * 2);
  const edgeLookup = new Map<string, number>();
  let e = 0;
  for (let i = 0; i < n_cells; i++) {
    for (let j = i + 1; j < n_cells; j++) {
      cells_on_edge[e * 2 + 0] = i;
      cells_on_edge[e * 2 + 1] = j;
      edgeLookup.set(`${i}_${j}`, e);
      edgeLookup.set(`${j}_${i}`, e);
      e++;
    }
  }
  const edges_on_cell = new Array<number>(n_cells * max_edges);
  for (let c = 0; c < n_cells; c++) {
    for (let j = 0; j < max_edges; j++) {
      const nb = cells_on_cell[c * max_edges + j];
      const key = c < nb ? `${c}_${nb}` : `${nb}_${c}`;
      edges_on_cell[c * max_edges + j] = edgeLookup.get(key) as number;
    }
  }
  const lon_edge = new Array<number>(n_edges);
  const lat_edge = new Array<number>(n_edges);
  const dc_edge = new Array<number>(n_edges);
  for (let ei = 0; ei < n_edges; ei++) {
    const c1 = cells_on_edge[ei * 2 + 0];
    const c2 = cells_on_edge[ei * 2 + 1];
    const v1: [number, number, number] = [
      x_cell[c1] / R,
      y_cell[c1] / R,
      z_cell[c1] / R,
    ];
    const v2: [number, number, number] = [
      x_cell[c2] / R,
      y_cell[c2] / R,
      z_cell[c2] / R,
    ];
    let d = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    if (d > 1) d = 1;
    if (d < -1) d = -1;
    dc_edge[ei] = R * Math.acos(d);
    const mx = (v1[0] + v2[0]) / 2;
    const my = (v1[1] + v2[1]) / 2;
    const mz = (v1[2] + v2[2]) / 2;
    const n = Math.sqrt(mx * mx + my * my + mz * mz);
    lon_edge[ei] = Math.atan2(my / n, mx / n);
    lat_edge[ei] = Math.asin(mz / n);
  }
  const dv_edge = dc_edge.slice();
  return mpasMeshData({
    lon_cell,
    lat_cell,
    area_cell,
    n_edges_on_cell,
    cells_on_cell,
    edges_on_cell,
    lon_edge,
    lat_edge,
    cells_on_edge,
    dc_edge,
    dv_edge,
    max_edges,
    x_cell,
    y_cell,
    z_cell,
    n_vertices: n_cells,
    R,
  });
}

describe("earthsci.grids.mpas", () => {
  it("is exposed via both namespaces", () => {
    expect(typeof grids.mpas).toBe("function");
    expect(earthsci.grids.mpas).toBe(grids.mpas);
  });

  it("constructs from an in-memory mesh", () => {
    const g: MPASGrid = grids.mpas({ mesh: tetraMesh(), R: 1 });
    expect(g.family).toBe("mpas");
    expect(g.topology).toBe("unstructured");
    expect(g.dtype).toBe("float64");
    expect(g.n_cells).toBe(4);
    expect(g.n_edges).toBe(6);
    expect(g.max_edges).toBe(3);
  });

  it("topology consistency: neighbor symmetry + edge refs (tetra mesh)", () => {
    const g = grids.mpas({ mesh: tetraMesh() });
    for (let c = 0; c < g.n_cells; c++) {
      const nbs = g.neighbors(c);
      expect(nbs.length).toBe(3);
      for (const nb of nbs) {
        expect(nb).toBeGreaterThanOrEqual(0);
        expect(nb).toBeLessThan(g.n_cells);
        expect(g.neighbors(nb).includes(c)).toBe(true);
      }
    }
    for (let e = 0; e < g.n_edges; e++) {
      const c1 = g.mesh.cells_on_edge[e * 2 + 0];
      const c2 = g.mesh.cells_on_edge[e * 2 + 1];
      expect(c1).toBeGreaterThanOrEqual(0);
      expect(c1).toBeLessThan(g.n_cells);
      expect(c2).toBeGreaterThanOrEqual(0);
      expect(c2).toBeLessThan(g.n_cells);
      expect(c1).not.toBe(c2);
    }
  });

  it("metric accuracy: 4πR² total area + edge arc = R·arccos(v1·v2)", () => {
    const R = 2.5;
    const g = grids.mpas({ mesh: tetraMesh(R), R });
    const expectedTotal = 4 * Math.PI * R * R;
    expect(Math.abs(g.total_area() - expectedTotal) / expectedTotal).toBeLessThan(1e-12);
    let acc = 0;
    for (let c = 0; c < g.n_cells; c++) acc += g.cell_area(c);
    expect(Math.abs(acc - expectedTotal) / expectedTotal).toBeLessThan(1e-12);

    for (let e = 0; e < g.n_edges; e++) {
      const c1 = g.mesh.cells_on_edge[e * 2 + 0];
      const c2 = g.mesh.cells_on_edge[e * 2 + 1];
      const p1 = g.cell_center_cart(c1);
      const p2 = g.cell_center_cart(c2);
      let d = (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z) / (R * R);
      if (d > 1) d = 1;
      if (d < -1) d = -1;
      const expected = R * Math.acos(d);
      expect(Math.abs(g.edge_length(e) - expected) / R).toBeLessThan(1e-12);
      expect(Math.abs(g.cell_distance(e) - expected) / R).toBeLessThan(1e-12);
    }
  });

  it("pinned accessor outputs (cross-binding reference)", () => {
    // Vertex 0 = (1,1,1)/√3; matches Julia test's pinned query point.
    const g = grids.mpas({ mesh: tetraMesh(1.0), R: 1.0 });
    const { lon, lat } = g.cell_centers(0);
    expect(Math.abs(lon - Math.PI / 4)).toBeLessThan(1e-14);
    expect(Math.abs(lat - Math.asin(1 / Math.sqrt(3)))).toBeLessThan(1e-14);
    const { x, y, z } = g.cell_center_cart(0);
    expect(Math.abs(Math.sqrt(x * x + y * y + z * z) - 1)).toBeLessThan(1e-14);
    expect(g.metric_eval("lon", 0)).toBe(lon);
    expect(g.metric_eval("area", 0)).toBe(g.cell_area(0));
    expect(g.metric_eval("dc_edge", 0)).toBe(g.cell_distance(0));
    expect(() => g.metric_eval("nonsense" as "area", 0)).toThrow(RangeError);
  });

  it("cell centers lie on the sphere of radius R (within 1e-12)", () => {
    const R = 6.371e6;
    const g = grids.mpas({ mesh: tetraMesh(R), R });
    for (let c = 0; c < g.n_cells; c++) {
      const { x, y, z } = g.cell_center_cart(c);
      const r = Math.sqrt(x * x + y * y + z * z);
      expect(Math.abs(r - R) / R).toBeLessThan(1e-12);
    }
  });

  it("rejects invalid cell / edge indices", () => {
    const g = grids.mpas({ mesh: tetraMesh() });
    expect(() => g.cell_centers(-1)).toThrow(RangeError);
    expect(() => g.cell_centers(g.n_cells)).toThrow(RangeError);
    expect(() => g.neighbors(3.5)).toThrow(RangeError);
    expect(() => g.cell_area(-1)).toThrow(RangeError);
    expect(() => g.edge_length(-1)).toThrow(RangeError);
    expect(() => g.edge_length(g.n_edges)).toThrow(RangeError);
    expect(() => g.metric_eval("lon_edge", g.n_edges)).toThrow(RangeError);
  });

  it("rejects invalid generator options", () => {
    expect(() => grids.mpas(undefined as unknown as { mesh: MPASMeshData })).toThrow(TypeError);
    expect(() => grids.mpas({})).toThrow(TypeError);
    const mesh = tetraMesh();
    expect(() => grids.mpas({ mesh, R: -1 })).toThrow(RangeError);
    expect(() => grids.mpas({ mesh, R: NaN })).toThrow(RangeError);
    expect(() => grids.mpas({ mesh, ghosts: 2 })).toThrow(RangeError);
    expect(() => grids.mpas({ mesh, ghosts: -1 })).toThrow(RangeError);
    expect(() =>
      grids.mpas({ mesh, dtype: "bogus" as "float64" }),
    ).toThrow(RangeError);
    // path-based loader without readerFn
    expect(() =>
      grids.mpas({ loader: { path: "/tmp/fake.nc", reader: "nc4" } }),
    ).toThrow(TypeError);
    // invalid loader fields
    expect(() =>
      grids.mpas({
        mesh,
        loader: { path: "x", reader: "bogus" as "nc4" },
      }),
    ).toThrow(RangeError);
    expect(() =>
      grids.mpas({
        mesh,
        loader: { path: "x", check: "bogus" as "strict" },
      }),
    ).toThrow(RangeError);
    expect(() =>
      grids.mpas({
        mesh,
        loader: { path: "" },
      }),
    ).toThrow(TypeError);
  });

  it("path-based loading uses readerFn and records loader in provenance", () => {
    const mesh = tetraMesh();
    let capturedPath = "";
    const g = grids.mpas({
      loader: { path: "/tmp/fake.nc", reader: "mpas_mesh", check: "lenient" },
      readerFn: (p) => {
        capturedPath = p;
        return mesh;
      },
    });
    expect(capturedPath).toBe("/tmp/fake.nc");
    expect(g.loader).toEqual({
      path: "/tmp/fake.nc",
      reader: "mpas_mesh",
      check: "lenient",
    });
    const esm = g.toESM() as Record<string, unknown>;
    const options = esm.options as Record<string, unknown>;
    const loader = options.loader as Record<string, unknown>;
    expect(loader.path).toBe("/tmp/fake.nc");
    expect(loader.reader).toBe("mpas_mesh");
    expect(loader.check).toBe("lenient");
  });

  it("strict check catches broken neighbor symmetry", () => {
    const mesh = tetraMesh();
    // Corrupt: cell 0 forgets its first neighbor (cell 1). Cell 1 still
    // claims 0 as a neighbor → reciprocity broken.
    mesh.cells_on_cell[0] = -1;
    expect(() => checkMesh(mesh, true)).toThrow(/neighbor symmetry/);
    // Lenient mode passes (bounds are still valid).
    expect(() => checkMesh(mesh, false)).not.toThrow();
  });

  describe("toESM declarative config", () => {
    it("emits a small declarative config (no inline geometry)", () => {
      const g = grids.mpas({
        mesh: tetraMesh(1.0),
        loader: { path: "/tmp/fake.nc", reader: "mpas_mesh", check: "lenient" },
        R: 1.0,
      });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.family).toBe("mpas");
      expect(esm.version).toBe("1.0.0");
      expect(esm.dtype).toBe("float64");
      expect(esm.topology).toBe("unstructured");
      expect(esm.ghosts).toBe(0);
      expect(esm.n_cells).toBe(4);
      expect(esm.n_edges).toBe(6);
      expect(esm.max_edges).toBe(3);
      const options = esm.options as Record<string, unknown>;
      expect(options.R).toBe(1.0);
      const loader = options.loader as Record<string, unknown>;
      expect(loader.path).toBe("/tmp/fake.nc");
      expect(loader.reader).toBe("mpas_mesh");
      expect(loader.check).toBe("lenient");

      const prov = esm.provenance as Record<string, unknown>;
      expect(prov.binding).toBe("typescript");
      expect(prov.source_sha).toBe("dsc-3nw");
      expect(prov.reader_version).toBeDefined();

      // Mayor 2026-04-20 correction: no inline geometry arrays.
      for (const forbidden of [
        "cells",
        "edges",
        "vertices",
        "lon_cell",
        "lat_cell",
        "area_cell",
        "cells_on_cell",
        "cells_on_edge",
      ]) {
        expect(esm).not.toHaveProperty(forbidden);
      }

      // JSON round-trips cleanly.
      const s = JSON.stringify(esm);
      expect(JSON.parse(s)).toEqual(esm);
    });

    it("omits loader from provenance when mesh-only", () => {
      const g = grids.mpas({ mesh: tetraMesh(), R: 1.0 });
      expect(g.loader).toBeNull();
      const esm = g.toESM() as Record<string, unknown>;
      const options = esm.options as Record<string, unknown>;
      expect(options.loader).toBeNull();
    });

    it("coerces loader defaults (reader='auto', check='strict')", () => {
      const g = grids.mpas({
        mesh: tetraMesh(),
        loader: { path: "/tmp/x.nc" },
      });
      expect(g.loader).toEqual({
        path: "/tmp/x.nc",
        reader: "auto",
        check: "strict",
      });
    });
  });

  it("dtype propagation (float32)", () => {
    const g = grids.mpas({ mesh: tetraMesh(), dtype: "float32" });
    expect(g.dtype).toBe("float32");
    const esm = g.toESM() as Record<string, unknown>;
    expect(esm.dtype).toBe("float32");
  });

  it("determinism: two calls with the same mesh produce identical toESM", () => {
    const m1 = tetraMesh(1.0);
    const m2 = tetraMesh(1.0);
    const g1 = grids.mpas({ mesh: m1, R: 1.0 });
    const g2 = grids.mpas({ mesh: m2, R: 1.0 });
    expect(JSON.stringify(g1.toESM())).toBe(JSON.stringify(g2.toESM()));
    for (let c = 0; c < g1.n_cells; c++) {
      expect(g2.cell_centers(c)).toEqual(g1.cell_centers(c));
      expect(g2.neighbors(c)).toEqual(g1.neighbors(c));
      expect(g2.cell_area(c)).toBe(g1.cell_area(c));
    }
  });
});
