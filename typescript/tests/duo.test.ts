import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { DuoGrid } from "../src/grids/index.js";

const builtin = (level: number) => ({
  path: `builtin://icosahedral/${level}`,
  reader: "builtin_icosahedral",
  check: "strict",
});

describe("earthsci.grids.duo", () => {
  it("is exposed via both namespaces", () => {
    expect(typeof grids.duo).toBe("function");
    expect(earthsci.grids.duo).toBe(grids.duo);
  });

  it("level 0 reproduces the bare icosahedron (20 / 12 / 30)", () => {
    const g: DuoGrid = grids.duo({ loader: builtin(0), R: 1 });
    expect(g.family).toBe("duo");
    expect(g.topology).toBe("unstructured");
    expect(g.level).toBe(0);
    expect(g.n_cells).toBe(20);
    expect(g.n_vertices).toBe(12);
    expect(g.n_edges).toBe(30);
    // Euler characteristic on the sphere: V - E + F = 2.
    expect(g.n_vertices - g.n_edges + g.n_cells).toBe(2);
  });

  it("subdivision counts match 20·4^r / 10·4^r+2 / 30·4^r", () => {
    for (let r = 0; r <= 3; r++) {
      const g = grids.duo({ loader: builtin(r), R: 1 });
      expect(g.n_cells).toBe(20 * 4 ** r);
      expect(g.n_vertices).toBe(10 * 4 ** r + 2);
      expect(g.n_edges).toBe(30 * 4 ** r);
      expect(g.n_vertices - g.n_edges + g.n_cells).toBe(2);
    }
  });

  it("summed cell areas equal 4π R² (closed unit sphere)", () => {
    for (let r = 0; r <= 3; r++) {
      const g = grids.duo({ loader: builtin(r), R: 1 });
      let total = 0;
      for (let c = 0; c < g.n_cells; c++) total += g.cell_area(c);
      expect(Math.abs(total - 4 * Math.PI) / (4 * Math.PI)).toBeLessThan(1e-10);
    }
  });

  it("summed cell areas scale as R²", () => {
    const R = 6.371e6;
    const g = grids.duo({ loader: builtin(2), R });
    let total = 0;
    for (let c = 0; c < g.n_cells; c++) total += g.cell_area(c);
    expect(Math.abs(total - 4 * Math.PI * R * R) / (4 * Math.PI * R * R)).toBeLessThan(1e-10);
  });

  it("every cell has strictly positive area", () => {
    const g = grids.duo({ loader: builtin(2), R: 1 });
    for (let c = 0; c < g.n_cells; c++) {
      expect(g.cell_area(c)).toBeGreaterThan(0);
    }
  });

  it("neighbors are symmetric and exhaustive (closed mesh → 3 per cell)", () => {
    const g = grids.duo({ loader: builtin(2), R: 1 });
    for (let c = 0; c < g.n_cells; c++) {
      const nbs = g.neighbors(c);
      expect(nbs.length).toBe(3);
      for (const nb of nbs) {
        // Closed icosahedral mesh has no boundaries.
        expect(nb).toBeGreaterThanOrEqual(0);
        expect(nb).toBeLessThan(g.n_cells);
        // Reciprocity.
        const back = g.neighbors(nb);
        expect(back.includes(c)).toBe(true);
      }
    }
  });

  it("cell_centers lie on the sphere of radius R (within 1e-12)", () => {
    const R = 6.371e6;
    const g = grids.duo({ loader: builtin(2), R });
    for (let c = 0; c < g.n_cells; c++) {
      const x = g.metric_eval("x", c);
      const y = g.metric_eval("y", c);
      const z = g.metric_eval("z", c);
      const r = Math.sqrt(x * x + y * y + z * z);
      expect(Math.abs(r - R) / R).toBeLessThan(1e-12);
    }
  });

  it("cell_centers returns finite lon/lat in [-π,π]×[-π/2,π/2]", () => {
    const g = grids.duo({ loader: builtin(1), R: 1 });
    for (let c = 0; c < g.n_cells; c++) {
      const { lon, lat } = g.cell_centers(c);
      expect(Number.isFinite(lon)).toBe(true);
      expect(Number.isFinite(lat)).toBe(true);
      expect(lon).toBeGreaterThanOrEqual(-Math.PI);
      expect(lon).toBeLessThanOrEqual(Math.PI);
      expect(lat).toBeGreaterThanOrEqual(-Math.PI / 2);
      expect(lat).toBeLessThanOrEqual(Math.PI / 2);
    }
  });

  it("metric_eval exposes area/lon/lat/x/y/z and rejects unknowns", () => {
    const g = grids.duo({ loader: builtin(1), R: 6.371e6 });
    for (const c of [0, 5, g.n_cells - 1]) {
      expect(g.metric_eval("area", c)).toBe(g.cell_area(c));
      const { lon, lat } = g.cell_centers(c);
      expect(g.metric_eval("lon", c)).toBe(lon);
      expect(g.metric_eval("lat", c)).toBe(lat);
    }
    expect(() => g.metric_eval("bogus" as "area", 0)).toThrow(RangeError);
  });

  it("rejects invalid cell indices", () => {
    const g = grids.duo({ loader: builtin(1), R: 1 });
    expect(() => g.cell_centers(-1)).toThrow(RangeError);
    expect(() => g.cell_centers(g.n_cells)).toThrow(RangeError);
    expect(() => g.neighbors(3.5)).toThrow(RangeError);
    expect(() => g.metric_eval("area", -1)).toThrow(RangeError);
    expect(() => g.cell_area(g.n_cells + 1)).toThrow(RangeError);
  });

  it("rejects missing / invalid options per §9", () => {
    expect(() => grids.duo(undefined as unknown as { loader: { path: string } })).toThrow(TypeError);
    expect(() =>
      grids.duo({} as unknown as { loader: { path: string } }),
    ).toThrow(TypeError);
    // loader.path required
    expect(() =>
      grids.duo({ loader: {} as unknown as { path: string } }),
    ).toThrow(TypeError);
    // duo_mesh reader not yet implemented
    expect(() =>
      grids.duo({ loader: { path: "/tmp/nope.duo", reader: "duo_mesh" } }),
    ).toThrow(RangeError);
    // unknown reader + non-builtin path
    expect(() =>
      grids.duo({ loader: { path: "bogus://x", reader: "unknown_reader" } }),
    ).toThrow(RangeError);
    // negative subdivision level
    expect(() => grids.duo({ loader: builtin(-1) })).toThrow(RangeError);
    // bad R
    expect(() => grids.duo({ loader: builtin(0), R: -1 })).toThrow(RangeError);
    // bad ghosts
    expect(() => grids.duo({ loader: builtin(0), ghosts: -1 })).toThrow(RangeError);
    // bad dtype
    expect(() =>
      grids.duo({ loader: builtin(0), dtype: "bogus" as "float64" }),
    ).toThrow(RangeError);
  });

  it("determinism: two calls with the same opts produce identical grids", () => {
    const g1 = grids.duo({ loader: builtin(2), R: 1 });
    const g2 = grids.duo({ loader: builtin(2), R: 1 });
    for (let c = 0; c < g1.n_cells; c++) {
      expect(g2.cell_centers(c)).toEqual(g1.cell_centers(c));
      expect(g2.neighbors(c)).toEqual(g1.neighbors(c));
      expect(g2.cell_area(c)).toBe(g1.cell_area(c));
    }
    expect(JSON.stringify(g1.toESM())).toBe(JSON.stringify(g2.toESM()));
  });

  describe("toESM declarative config", () => {
    it("emits a small declarative config (no inline geometry)", () => {
      const g = grids.duo({ loader: builtin(1), R: 6.371e6 });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.family).toBe("duo");
      expect(esm.version).toBe("1.0.0");
      expect(esm.dtype).toBe("float64");
      expect(esm.topology).toBe("unstructured");
      expect(esm.n_cells).toBe(g.n_cells);
      expect(esm.n_vertices).toBe(g.n_vertices);
      expect(esm.n_edges).toBe(g.n_edges);

      const options = esm.options as Record<string, unknown>;
      expect(options.R).toBe(6.371e6);
      expect(options.level).toBe(1);
      const loader = options.loader as Record<string, unknown>;
      expect(loader.path).toBe("builtin://icosahedral/1");
      expect(loader.reader).toBe("builtin_icosahedral");
      expect(loader.check).toBe("strict");

      expect(esm.provenance).toBeDefined();
      // No inline geometry per 2026-04-20 mayor correction.
      for (const forbidden of [
        "cells",
        "vertices",
        "vertices_xyz",
        "edges",
        "edges_list",
        "faces",
        "lon",
        "lat",
      ]) {
        expect(esm).not.toHaveProperty(forbidden);
      }

      // JSON round-trips cleanly.
      const s = JSON.stringify(esm);
      expect(JSON.parse(s)).toEqual(esm);
    });

    it("defaults loader.reader='auto' and loader.check='strict' are preserved", () => {
      const g = grids.duo({ loader: { path: "builtin://icosahedral/0" }, R: 1 });
      // auto reader on a builtin path resolves to level 0 without touching the
      // unimplemented .duo file reader.
      expect(g.level).toBe(0);
      const esm = g.toESM() as Record<string, unknown>;
      const options = esm.options as Record<string, unknown>;
      const loader = options.loader as Record<string, unknown>;
      expect(loader.reader).toBe("auto");
      expect(loader.check).toBe("strict");
    });
  });

  it("accepts dtype 'float32' and records it in toESM / grid", () => {
    const g = grids.duo({ loader: builtin(1), R: 1, dtype: "float32" });
    expect(g.dtype).toBe("float32");
    const esm = g.toESM() as Record<string, unknown>;
    expect(esm.dtype).toBe("float32");
  });
});
