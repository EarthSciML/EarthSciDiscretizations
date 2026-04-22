import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { CartesianGrid } from "../src/grids/index.js";

const REL_TOL = 1e-14;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

describe("earthsci.grids.cartesian", () => {
  it("is exposed via both grids and earthsci.grids namespaces", () => {
    expect(typeof grids.cartesian).toBe("function");
    expect(earthsci.grids.cartesian).toBe(grids.cartesian);
  });

  describe("construction — uniform", () => {
    it("builds a 1D uniform grid from nx + extent", () => {
      const g: CartesianGrid = grids.cartesian({
        nx: 10,
        extent: [[0, 1]],
      });
      expect(g.family).toBe("cartesian");
      expect(g.topology).toBe("rectilinear");
      expect(g.dtype).toBe("float64");
      expect(g.ndim).toBe(1);
      expect(g.n).toEqual([10]);
      expect(g.ghosts).toBe(0);
      expect(g.n_cells).toBe(10);
      expect(g.uniform).toEqual([true]);
      expect(g.edges[0].length).toBe(11);
      expect(g.centers[0].length).toBe(10);
      expect(g.edges[0][0]).toBeCloseTo(0, 14);
      expect(g.edges[0][10]).toBeCloseTo(1, 14);
      for (const w of g.widths[0]) expect(closeRel(w, 0.1, 1e-12)).toBe(true);
    });

    it("builds a 2D uniform grid", () => {
      const g = grids.cartesian({
        nx: 4,
        ny: 8,
        extent: [
          [0, 2],
          [0, 4],
        ],
      });
      expect(g.ndim).toBe(2);
      expect(g.n).toEqual([4, 8]);
      expect(g.extent).toEqual([
        [0, 2],
        [0, 4],
      ]);
      expect(g.n_cells).toBe(32);
      for (const w of g.widths[0]) expect(closeRel(w, 0.5, 1e-12)).toBe(true);
      for (const w of g.widths[1]) expect(closeRel(w, 0.5, 1e-12)).toBe(true);
    });

    it("builds a 3D uniform grid", () => {
      const g = grids.cartesian({
        nx: 2,
        ny: 3,
        nz: 4,
        extent: [
          [0, 1],
          [0, 1.5],
          [0, 2],
        ],
      });
      expect(g.ndim).toBe(3);
      expect(g.n).toEqual([2, 3, 4]);
      expect(g.uniform).toEqual([true, true, true]);
      expect(g.n_cells).toBe(24);
    });
  });

  describe("construction — non-uniform", () => {
    it("builds a 2D grid from explicit edges", () => {
      const xe = [0, 0.1, 0.3, 0.7, 1];
      const ye = [-1, 0, 0.5, 1];
      const g = grids.cartesian({ edges: [xe, ye] });
      expect(g.ndim).toBe(2);
      expect(g.n).toEqual([4, 3]);
      expect(g.uniform).toEqual([false, false]);
      expect(g.extent).toEqual([
        [0, 1],
        [-1, 1],
      ]);
      const expectedX = [0.1, 0.2, 0.4, 0.3];
      for (let i = 0; i < expectedX.length; i++) {
        expect(closeRel(g.widths[0][i], expectedX[i], 1e-12)).toBe(true);
      }
      const expectedY = [1, 0.5, 0.5];
      for (let i = 0; i < expectedY.length; i++) {
        expect(closeRel(g.widths[1][i], expectedY[i], 1e-12)).toBe(true);
      }
    });

    it("flags per-axis uniformity correctly when edges happen to be uniform", () => {
      const g = grids.cartesian({ edges: [[0, 0.25, 0.5, 0.75, 1]] });
      expect(g.uniform).toEqual([true]);
    });
  });

  describe("options — dtype, ghosts", () => {
    it("dtype=float32 is honored", () => {
      const g = grids.cartesian({
        nx: 4,
        extent: [[0, 1]],
        dtype: "float32",
      });
      expect(g.dtype).toBe("float32");
    });

    it("ghosts is honored", () => {
      const g = grids.cartesian({ nx: 4, extent: [[0, 1]], ghosts: 3 });
      expect(g.ghosts).toBe(3);
    });
  });

  describe("accessor contract", () => {
    it("cell_centers derives from extent (uniform 1D)", () => {
      const g = grids.cartesian({ nx: 4, extent: [[0, 1]] });
      expect(g.cell_centers(0)).toEqual([0.125]);
      expect(g.cell_centers(1)).toEqual([0.375]);
      expect(g.cell_centers(2)).toEqual([0.625]);
      expect(g.cell_centers(3)).toEqual([0.875]);
    });

    it("cell_centers in 2D", () => {
      const g = grids.cartesian({
        nx: 2,
        ny: 2,
        extent: [
          [0, 2],
          [10, 14],
        ],
      });
      expect(g.cell_centers(0, 0)).toEqual([0.5, 11]);
      expect(g.cell_centers(1, 1)).toEqual([1.5, 13]);
    });

    it("cell_volume = product of widths", () => {
      const g3 = grids.cartesian({
        nx: 2,
        ny: 2,
        nz: 2,
        extent: [
          [0, 2],
          [0, 4],
          [0, 8],
        ],
      });
      expect(closeRel(g3.cell_volume(0, 0, 0), 1 * 2 * 4, 1e-12)).toBe(true);
      const gNu = grids.cartesian({ edges: [[0, 0.5, 1.5]] });
      expect(closeRel(gNu.cell_volume(0), 0.5, 1e-12)).toBe(true);
      expect(closeRel(gNu.cell_volume(1), 1.0, 1e-12)).toBe(true);
    });

    it("cell_widths returns per-axis widths", () => {
      const g = grids.cartesian({ edges: [[0, 0.1, 0.4]] });
      expect(g.cell_widths(0)).toEqual([0.1]);
      expect(closeRel(g.cell_widths(1)[0], 0.3, 1e-12)).toBe(true);
    });

    it("neighbors omit out-of-bounds sides", () => {
      const g = grids.cartesian({
        nx: 3,
        ny: 3,
        extent: [
          [0, 1],
          [0, 1],
        ],
      });
      const nb = g.neighbors(1, 1);
      expect(nb.length).toBe(4);
      const byKey: Record<string, number[]> = {};
      for (const r of nb) byKey[`${r.axis}:${r.side}`] = r.index;
      expect(byKey["0:-1"]).toEqual([0, 1]);
      expect(byKey["0:1"]).toEqual([2, 1]);
      expect(byKey["1:-1"]).toEqual([1, 0]);
      expect(byKey["1:1"]).toEqual([1, 2]);

      const corner = g.neighbors(0, 0);
      expect(corner.length).toBe(2);
      const hasSides = new Set(corner.map((r) => `${r.axis}:${r.side}`));
      expect(hasSides.has("0:1")).toBe(true);
      expect(hasSides.has("1:1")).toBe(true);
      expect(hasSides.has("0:-1")).toBe(false);
      expect(hasSides.has("1:-1")).toBe(false);
    });

    it("metric_eval: volume / dx / face_area / jacobian / g", () => {
      const g = grids.cartesian({
        nx: 2,
        ny: 4,
        extent: [
          [0, 2],
          [0, 4],
        ],
      });
      expect(closeRel(g.metric_eval("volume", 0, 0) as number, 1, 1e-12)).toBe(true);
      expect(closeRel(g.metric_eval("dx", 0, 0) as number, 1, 1e-12)).toBe(true);
      expect(closeRel(g.metric_eval("dy", 0, 0) as number, 1, 1e-12)).toBe(true);
      expect(closeRel(g.metric_eval("face_area_x", 0, 0) as number, 1, 1e-12)).toBe(true);
      expect(closeRel(g.metric_eval("face_area_y", 0, 0) as number, 1, 1e-12)).toBe(true);
      expect(g.metric_eval("jacobian", 0, 0)).toBe(1);
      expect(g.metric_eval("g", 0, 0)).toEqual([
        [1, 0],
        [0, 1],
      ]);
    });

    it("metric_eval rejects unknown name / axis not present", () => {
      const g1 = grids.cartesian({ nx: 2, extent: [[0, 1]] });
      expect(() => g1.metric_eval("dy", 0)).toThrow(RangeError);
      expect(() => g1.metric_eval("dz", 0)).toThrow(RangeError);
      expect(() => g1.metric_eval("face_area_y", 0)).toThrow(RangeError);
      expect(() =>
        g1.metric_eval("bogus" as "volume", 0),
      ).toThrow(RangeError);
    });

    it("rejects out-of-range or wrong-arity indices", () => {
      const g = grids.cartesian({
        nx: 2,
        ny: 2,
        extent: [
          [0, 1],
          [0, 1],
        ],
      });
      expect(() => g.cell_centers(-1, 0)).toThrow(RangeError);
      expect(() => g.cell_centers(0, 2)).toThrow(RangeError);
      expect(() => g.cell_centers(0.5, 0)).toThrow(RangeError);
      // wrong arity
      expect(() => g.cell_centers(0)).toThrow(RangeError);
      expect(() => g.cell_centers(0, 0, 0)).toThrow(RangeError);
    });
  });

  describe("topology consistency", () => {
    it("edges bound centers, widths sum to extent span (uniform)", () => {
      for (const nx of [1, 5, 16]) {
        for (const ext of [
          [0, 1] as [number, number],
          [-3, 7.5] as [number, number],
        ]) {
          const g = grids.cartesian({ nx, extent: [ext] });
          for (let i = 0; i < nx; i++) {
            expect(g.edges[0][i]).toBeLessThan(g.centers[0][i]);
            expect(g.centers[0][i]).toBeLessThan(g.edges[0][i + 1]);
          }
          const sum = g.widths[0].reduce((p, v) => p + v, 0);
          expect(closeRel(sum, ext[1] - ext[0], 1e-12)).toBe(true);
          expect(closeRel(g.edges[0][0], ext[0], 1e-12)).toBe(true);
          expect(closeRel(g.edges[0][nx], ext[1], 1e-12)).toBe(true);
        }
      }
    });

    it("partition of unity on non-uniform 2D grid", () => {
      const g = grids.cartesian({
        edges: [
          [0, 0.3, 1],
          [0, 0.5, 2],
        ],
      });
      let total = 0;
      for (let i = 0; i < g.n[0]; i++) {
        for (let j = 0; j < g.n[1]; j++) {
          total += g.cell_volume(i, j);
        }
      }
      expect(closeRel(total, 1 * 2, 1e-12)).toBe(true);
    });
  });

  describe("metric accuracy", () => {
    it("uniform dx is exact", () => {
      const nx = 7;
      const extent: [number, number] = [-1, 1];
      const g = grids.cartesian({ nx, extent: [extent] });
      const expected = (extent[1] - extent[0]) / nx;
      for (let i = 0; i < nx; i++) {
        expect(
          closeRel(g.metric_eval("dx", i) as number, expected, 1e-14),
        ).toBe(true);
      }
    });

    it("non-uniform widths match edges exactly", () => {
      const xe = [0, 1, 1.5, 4, 4.25];
      const g = grids.cartesian({ edges: [xe] });
      for (let i = 0; i + 1 < xe.length; i++) {
        expect(
          closeRel(g.metric_eval("dx", i) as number, xe[i + 1] - xe[i], 1e-14),
        ).toBe(true);
      }
    });
  });

  describe("error contract — §9", () => {
    it("missing required options raise TypeError", () => {
      expect(() => grids.cartesian({} as Parameters<typeof grids.cartesian>[0])).toThrow(TypeError);
      expect(() => grids.cartesian({ nx: 4 })).toThrow(TypeError); // missing extent
      expect(() =>
        grids.cartesian({ nx: 4, ny: 2 } as Parameters<typeof grids.cartesian>[0]),
      ).toThrow(TypeError); // 2D missing extent
      expect(() =>
        grids.cartesian({
          nx: 4,
          nz: 2,
          extent: [
            [0, 1],
            [0, 1],
          ],
        } as Parameters<typeof grids.cartesian>[0]),
      ).toThrow(TypeError); // ny missing for 3D
    });

    it("invalid extent / negative ghosts raise RangeError", () => {
      expect(() => grids.cartesian({ nx: 4, extent: [[1, 0]] })).toThrow(RangeError);
      expect(() => grids.cartesian({ nx: 0, extent: [[0, 1]] })).toThrow(RangeError);
      expect(() => grids.cartesian({ nx: 3.5, extent: [[0, 1]] })).toThrow(RangeError);
      expect(() =>
        grids.cartesian({ nx: 4, extent: [[0, 1]], ghosts: -1 }),
      ).toThrow(RangeError);
    });

    it("non-uniform edges must be strictly increasing", () => {
      expect(() =>
        grids.cartesian({ edges: [[0, 0.5, 0.5, 1]] }),
      ).toThrow(RangeError); // duplicate
      expect(() =>
        grids.cartesian({ edges: [[0, 0.5, 0.25, 1]] }),
      ).toThrow(RangeError); // decreasing
      expect(() => grids.cartesian({ edges: [[0]] })).toThrow(RangeError); // too short
    });

    it("cannot mix `edges` with `nx`/`extent`", () => {
      expect(() =>
        grids.cartesian({ nx: 4, edges: [[0, 1, 2]] }),
      ).toThrow(TypeError);
      expect(() =>
        grids.cartesian({ extent: [[0, 1]], edges: [[0, 1, 2]] }),
      ).toThrow(TypeError);
    });

    it("non-uniform ndim cap (≤ 3)", () => {
      expect(() =>
        grids.cartesian({
          edges: [
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
          ],
        }),
      ).toThrow(RangeError);
    });

    it("rejects invalid dtype", () => {
      expect(() =>
        grids.cartesian({
          nx: 4,
          extent: [[0, 1]],
          dtype: "bogus" as "float64",
        }),
      ).toThrow(RangeError);
    });
  });

  describe("toESM declarative config", () => {
    it("emits §6 fields for uniform (no edges payload)", () => {
      const g = grids.cartesian({
        nx: 4,
        ny: 2,
        extent: [
          [0, 1],
          [0, 0.5],
        ],
        ghosts: 2,
      });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.family).toBe("cartesian");
      expect(esm.topology).toBe("rectilinear");
      expect(esm.version).toBe("1.0.0");
      expect(esm.dtype).toBe("float64");
      expect(esm.ndim).toBe(2);
      expect(esm.ghosts).toBe(2);
      expect(esm.n_cells).toBe(8);
      expect(esm.n).toEqual([4, 2]);
      expect(esm.uniform).toEqual([true, true]);
      expect(esm.extent).toEqual([
        [0, 1],
        [0, 0.5],
      ]);
      expect(esm).not.toHaveProperty("edges"); // uniform → no edges payload
      expect(esm).not.toHaveProperty("cells");
      expect(esm).not.toHaveProperty("vertices");
      const prov = esm.provenance as Record<string, unknown>;
      expect(prov.binding).toBe("typescript");
      expect(prov.binding_version).toBeDefined();
      expect(prov.math_lib).toBeDefined();

      // JSON round-trips cleanly.
      const s = JSON.stringify(esm);
      expect(JSON.parse(s)).toEqual(esm);
    });

    it("includes edges for non-uniform axes", () => {
      const xe = [0, 0.1, 0.4, 1];
      const g = grids.cartesian({ edges: [xe] });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.uniform).toEqual([false]);
      const e = esm.edges as number[][];
      expect(e).toBeDefined();
      expect(e[0]).toEqual(xe);
    });

    it("reflects dtype in lowering", () => {
      const g = grids.cartesian({
        nx: 2,
        extent: [[0, 1]],
        dtype: "float32",
      });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.dtype).toBe("float32");
    });

    it("is deterministic (same opts → same config bytes)", () => {
      const opts = {
        nx: 8,
        ny: 8,
        extent: [
          [0, 1],
          [0, 1],
        ] as AxisExtentLike,
        ghosts: 1,
      };
      const a = grids.cartesian(opts).toESM();
      const b = grids.cartesian(opts).toESM();
      expect(JSON.stringify(a)).toBe(JSON.stringify(b));
    });
  });
});

type AxisExtentLike = [number, number][];
