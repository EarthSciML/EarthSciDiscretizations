import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { CubedSphereGrid } from "../src/grids/index.js";

const REL_TOL = 1e-14;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

describe("earthsci.grids.cubed_sphere", () => {
  it("is exposed via both grids and earthsci.grids namespaces", () => {
    expect(typeof grids.cubed_sphere).toBe("function");
    expect(earthsci.grids.cubed_sphere).toBe(grids.cubed_sphere);
  });

  it("returns a Grid with the contracted shape", () => {
    const g: CubedSphereGrid = grids.cubed_sphere({ Nc: 4 });
    expect(g.family).toBe("cubed_sphere");
    expect(g.topology).toBe("block_structured");
    expect(g.generator).toBe("gnomonic_c6");
    expect(g.dtype).toBe("float64");
    expect(g.Nc).toBe(4);
    expect(g.R).toBeCloseTo(6.371e6);
    expect(g.ghosts).toBe(0);
    expect(g.n_panels).toBe(6);
    expect(g.n_cells).toBe(6 * 4 * 4);
  });

  it("honors option overrides (Nc, R, dtype, ghosts)", () => {
    const g = grids.cubed_sphere({ Nc: 8, R: 1, dtype: "float32", ghosts: 2 });
    expect(g.Nc).toBe(8);
    expect(g.R).toBe(1);
    expect(g.dtype).toBe("float32");
    expect(g.ghosts).toBe(2);
  });

  it("rejects missing / invalid options per §9 error contract", () => {
    expect(() => grids.cubed_sphere({} as unknown as { Nc: number })).toThrow(TypeError);
    expect(() => grids.cubed_sphere({ Nc: 0 })).toThrow(RangeError);
    expect(() => grids.cubed_sphere({ Nc: 3.5 })).toThrow(RangeError);
    expect(() => grids.cubed_sphere({ Nc: 4, R: -1 })).toThrow(RangeError);
    expect(() => grids.cubed_sphere({ Nc: 4, ghosts: -1 })).toThrow(RangeError);
    expect(() =>
      grids.cubed_sphere({ Nc: 4, dtype: "bogus" as "float64" }),
    ).toThrow(RangeError);
  });

  describe("topology consistency", () => {
    const g = grids.cubed_sphere({ Nc: 6, R: 1 });

    it("panel centers of panels 2 and 5 are the poles", () => {
      // Nc even → no cell exactly at ξ=0; use edge-symmetric reflection.
      // Instead check the mean of the four center-ish cells: lat → ±π/2 as Nc→∞.
      // Directly verify via (ξ=0, η=0) evaluation is covered by the math; use
      // a direct call to confirm panel 2 = north, panel 5 = south at origin.
      const fine = grids.cubed_sphere({ Nc: 1, R: 1 });
      expect(fine.cell_centers(2, 0, 0).lat).toBeCloseTo(Math.PI / 2, 12);
      expect(fine.cell_centers(5, 0, 0).lat).toBeCloseTo(-Math.PI / 2, 12);
      expect(fine.cell_centers(0, 0, 0).lat).toBeCloseTo(0, 12);
      expect(fine.cell_centers(0, 0, 0).lon).toBeCloseTo(0, 12);
    });

    it("each interior cell has 4 neighbors (W/E/S/N)", () => {
      for (let p = 0; p < 6; p++) {
        for (let i = 1; i < g.Nc - 1; i++) {
          for (let j = 1; j < g.Nc - 1; j++) {
            const nbs = g.neighbors(p, i, j);
            expect(nbs.length).toBe(4);
            const edges = nbs.map((n) => n.edge).sort();
            expect(edges).toEqual(["east", "north", "south", "west"]);
          }
        }
      }
    });

    it("boundary cells also produce 4 neighbors (cross-panel)", () => {
      for (let p = 0; p < 6; p++) {
        for (let i of [0, g.Nc - 1]) {
          for (let j = 0; j < g.Nc; j++) {
            expect(g.neighbors(p, i, j).length).toBe(4);
          }
        }
        for (let j of [0, g.Nc - 1]) {
          for (let i = 0; i < g.Nc; i++) {
            expect(g.neighbors(p, i, j).length).toBe(4);
          }
        }
      }
    });

    it("cross-panel neighbors cross back to the source panel (involutive)", () => {
      // For every boundary cell, the neighbor across the edge must have the
      // source cell in its own neighbor set (symmetry of adjacency).
      for (let p = 0; p < 6; p++) {
        const border: Array<[number, number]> = [];
        for (let k = 0; k < g.Nc; k++) {
          border.push([0, k], [g.Nc - 1, k], [k, 0], [k, g.Nc - 1]);
        }
        for (const [i, j] of border) {
          const nbs = g.neighbors(p, i, j);
          for (const nb of nbs) {
            if (nb.panel === p) continue;
            const back = g.neighbors(nb.panel, nb.i, nb.j);
            const found = back.some((b) => b.panel === p && b.i === i && b.j === j);
            expect(
              found,
              `panel ${p} (${i},${j}) → panel ${nb.panel} (${nb.i},${nb.j}) did not symmetrize`,
            ).toBe(true);
          }
        }
      }
    });
  });

  describe("metric accuracy", () => {
    it("J at the panel center equals R² for any panel (analytic identity)", () => {
      const R = 7;
      const g = grids.cubed_sphere({ Nc: 1, R });
      for (let p = 0; p < 6; p++) {
        const J = g.metric_eval("J", p, 0, 0);
        expect(closeRel(J, R * R)).toBe(true);
      }
    });

    it("metric inverse satisfies g·g⁻¹ = I at cell centers", () => {
      const g = grids.cubed_sphere({ Nc: 5, R: 6.371e6 });
      for (const [i, j] of [
        [0, 0],
        [2, 2],
        [4, 4],
        [0, 4],
      ] as const) {
        const gxx = g.metric_eval("g_xixi", 0, i, j);
        const gee = g.metric_eval("g_etaeta", 0, i, j);
        const gxe = g.metric_eval("g_xieta", 0, i, j);
        const ixx = g.metric_eval("ginv_xixi", 0, i, j);
        const iee = g.metric_eval("ginv_etaeta", 0, i, j);
        const ixe = g.metric_eval("ginv_xieta", 0, i, j);
        // g_ab · g^{bc} = δ_a^c; off-diagonal (xx·xe-inv_xieta + xe·ee-inv) = 0 etc.
        const p11 = gxx * ixx + gxe * ixe;
        const p22 = gxe * ixe + gee * iee;
        const p12 = gxx * ixe + gxe * iee;
        expect(closeRel(p11, 1, 1e-12)).toBe(true);
        expect(closeRel(p22, 1, 1e-12)).toBe(true);
        expect(Math.abs(p12)).toBeLessThan(1e-6 * Math.max(gxx, gee));
      }
    });

    it("summed cell areas equal 4π R² (surface area of sphere)", () => {
      const R = 1;
      const Nc = 8;
      const g = grids.cubed_sphere({ Nc, R });
      let total = 0;
      for (let p = 0; p < 6; p++) {
        for (let i = 0; i < Nc; i++) {
          for (let j = 0; j < Nc; j++) {
            total += g.cell_area(p, i, j);
          }
        }
      }
      // Spherical-excess formula is exact up to floating-point round-off.
      expect(Math.abs(total - 4 * Math.PI * R * R) / (4 * Math.PI)).toBeLessThan(1e-12);
    });

    it("cell areas are strictly positive", () => {
      const g = grids.cubed_sphere({ Nc: 4, R: 2 });
      for (let p = 0; p < 6; p++) {
        for (let i = 0; i < 4; i++) {
          for (let j = 0; j < 4; j++) {
            expect(g.cell_area(p, i, j)).toBeGreaterThan(0);
          }
        }
      }
    });

    it("rejects out-of-range panel / cell indices", () => {
      const g = grids.cubed_sphere({ Nc: 3 });
      expect(() => g.cell_centers(-1, 0, 0)).toThrow(RangeError);
      expect(() => g.cell_centers(6, 0, 0)).toThrow(RangeError);
      expect(() => g.cell_centers(0, 3, 0)).toThrow(RangeError);
      expect(() => g.cell_centers(0, 0, -1)).toThrow(RangeError);
      expect(() => g.metric_eval("bogus" as "J", 0, 0, 0)).toThrow(RangeError);
    });
  });

  describe("toESM declarative config", () => {
    it("returns a small declarative config (no per-cell geometry)", () => {
      const g = grids.cubed_sphere({ Nc: 48, R: 6.371e6, ghosts: 3 });
      const esm = g.toESM() as Record<string, unknown>;
      expect(esm.family).toBe("cubed_sphere");
      expect(esm.version).toBe("1.0.0");
      expect(esm.dtype).toBe("float64");
      expect(esm.topology).toBe("block_structured");
      expect(esm.generator).toBe("gnomonic_c6");
      expect(esm.n_panels).toBe(6);
      expect(esm.n_cells).toBe(6 * 48 * 48);
      const options = esm.options as Record<string, unknown>;
      expect(options.Nc).toBe(48);
      expect(options.R).toBe(6.371e6);
      expect(options.ghosts).toBe(3);

      // Must not inline cells/vertices/edges arrays per the 2026-04-20 mayor
      // correction (declarative config, not serialized geometry).
      expect(esm).not.toHaveProperty("cells");
      expect(esm).not.toHaveProperty("vertices");
      expect(esm).not.toHaveProperty("edges");
      expect(esm).not.toHaveProperty("lon");
      expect(esm).not.toHaveProperty("lat");

      // JSON round-trips cleanly.
      const s = JSON.stringify(esm);
      expect(JSON.parse(s)).toEqual(esm);
    });

    it("is deterministic (same opts → same config)", () => {
      const a = grids.cubed_sphere({ Nc: 16, R: 6.371e6 }).toESM();
      const b = grids.cubed_sphere({ Nc: 16, R: 6.371e6 }).toESM();
      expect(JSON.stringify(a)).toBe(JSON.stringify(b));
    });
  });
});
