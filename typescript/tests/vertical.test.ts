import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { VerticalGrid } from "../src/grids/index.js";

const REL_TOL = 1e-14;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function allClose(a: readonly number[], b: readonly number[], tol = 1e-12) {
  expect(a.length).toBe(b.length);
  for (let i = 0; i < a.length; i++) expect(closeRel(a[i], b[i], tol)).toBe(true);
}

describe("earthsci.grids.vertical", () => {
  // --------------------------------------------------------------- API
  describe("namespace wiring", () => {
    it("exposes vertical via both grids and earthsci.grids", () => {
      expect(typeof grids.vertical).toBe("function");
      expect(earthsci.grids.vertical).toBe(grids.vertical);
    });

    it("rejects missing options object", () => {
      // @ts-expect-error - exercising runtime guard
      expect(() => grids.vertical()).toThrow(TypeError);
    });

    it("rejects missing coordinate", () => {
      // @ts-expect-error - missing required field
      expect(() => grids.vertical({})).toThrow(TypeError);
    });

    it("rejects unknown coordinate", () => {
      expect(() =>
        // @ts-expect-error - outside the literal union
        grids.vertical({ coordinate: "foo", nz: 4 }),
      ).toThrow(RangeError);
    });

    it("rejects invalid dtype", () => {
      expect(() =>
        // @ts-expect-error - outside the literal union
        grids.vertical({ coordinate: "sigma", nz: 4, dtype: "float16" }),
      ).toThrow(RangeError);
    });

    it("rejects negative ghosts", () => {
      expect(() =>
        grids.vertical({ coordinate: "sigma", nz: 4, ghosts: -1 }),
      ).toThrow(RangeError);
    });

    it("rejects non-integer ghosts", () => {
      expect(() =>
        grids.vertical({ coordinate: "sigma", nz: 4, ghosts: 1.5 }),
      ).toThrow(RangeError);
    });

    it("applies default dtype / ghosts / p0", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(g.dtype).toBe("float64");
      expect(g.ghosts).toBe(0);
      expect(g.p0).toBe(1.0e5);
      expect(g.family).toBe("vertical");
      expect(g.topology).toBe("column");
      expect(g.ndim).toBe(1);
    });
  });

  // ------------------------------------------------------------ sigma
  describe("sigma coordinate", () => {
    it("builds a uniform column from nz", () => {
      const g: VerticalGrid = grids.vertical({ coordinate: "sigma", nz: 10 });
      expect(g.coordinate).toBe("sigma");
      expect(g.nz).toBe(10);
      expect(g.levels.length).toBe(11);
      expect(g.centers.length).toBe(10);
      expect(g.widths.length).toBe(10);
      expect(closeRel(g.levels[0], 1.0)).toBe(true);
      expect(closeRel(g.levels[10], 0.0)).toBe(true);
      for (const w of g.widths) expect(closeRel(w, 0.1, 1e-12)).toBe(true);
    });

    it("accepts explicit decreasing levels in [0, 1]", () => {
      const lv = [1.0, 0.9, 0.7, 0.3, 0.0];
      const g = grids.vertical({ coordinate: "sigma", levels: lv });
      expect(g.nz).toBe(4);
      allClose(g.levels, lv);
      allClose(g.centers, [0.95, 0.8, 0.5, 0.15]);
      allClose(g.widths, [0.1, 0.2, 0.4, 0.3]);
    });

    it("rejects non-decreasing levels", () => {
      expect(() =>
        grids.vertical({ coordinate: "sigma", levels: [1.0, 0.5, 0.7, 0.0] }),
      ).toThrow(RangeError);
    });

    it("rejects levels outside [0, 1]", () => {
      expect(() =>
        grids.vertical({ coordinate: "sigma", levels: [1.2, 0.5, 0.0] }),
      ).toThrow(RangeError);
      expect(() =>
        grids.vertical({ coordinate: "sigma", levels: [1.0, 0.5, -0.1] }),
      ).toThrow(RangeError);
    });

    it("requires nz or levels", () => {
      expect(() => grids.vertical({ coordinate: "sigma" })).toThrow(TypeError);
    });

    it("checks nz/levels consistency", () => {
      expect(() =>
        grids.vertical({
          coordinate: "sigma",
          nz: 10,
          levels: [1.0, 0.5, 0.0],
        }),
      ).toThrow(RangeError);
    });
  });

  // ------------------------------------------------------- z / theta / z_star
  describe("z / theta / z_star coordinates", () => {
    it("builds a z column from strictly-increasing levels", () => {
      const z = [0.0, 100.0, 300.0, 700.0, 1500.0];
      const g = grids.vertical({ coordinate: "z", levels: z });
      expect(g.coordinate).toBe("z");
      expect(g.nz).toBe(4);
      allClose(g.levels, z);
      allClose(g.centers, [50.0, 200.0, 500.0, 1100.0]);
      allClose(g.widths, [100.0, 200.0, 400.0, 800.0]);
    });

    it("requires levels for z", () => {
      expect(() => grids.vertical({ coordinate: "z", nz: 10 })).toThrow(
        TypeError,
      );
    });

    it("rejects non-monotone z levels", () => {
      expect(() =>
        grids.vertical({ coordinate: "z", levels: [0.0, 100.0, 50.0] }),
      ).toThrow(RangeError);
    });

    it("builds a theta column", () => {
      const g = grids.vertical({
        coordinate: "theta",
        levels: [280.0, 300.0, 340.0, 400.0],
      });
      expect(g.nz).toBe(3);
      allClose(g.centers, [290.0, 320.0, 370.0]);
      allClose(g.widths, [20.0, 40.0, 60.0]);
    });

    it("builds a z_star column", () => {
      const g = grids.vertical({
        coordinate: "z_star",
        levels: [0.0, 0.25, 0.75, 1.0],
      });
      expect(g.nz).toBe(3);
      allClose(g.widths, [0.25, 0.5, 0.25]);
    });
  });

  // -------------------------------------------------------------- eta
  describe("eta coordinate", () => {
    it("synthesizes sigma from ak/p0 + bk", () => {
      const p0 = 1.0e5;
      const ak = [0.0, 1000.0, 5000.0, 100.0];
      const bk = [1.0, 0.7, 0.3, 0.0];
      const g = grids.vertical({ coordinate: "eta", ak, bk, p0 });
      expect(g.coordinate).toBe("eta");
      expect(g.nz).toBe(3);
      expect(g.ak.length).toBe(4);
      expect(g.bk.length).toBe(4);
      const expected = ak.map((a, i) => a / p0 + bk[i]);
      allClose(g.levels, expected);
    });

    it("requires ak and bk", () => {
      expect(() =>
        grids.vertical({ coordinate: "eta", bk: [1.0, 0.5, 0.0] }),
      ).toThrow(TypeError);
      expect(() =>
        grids.vertical({ coordinate: "eta", ak: [0.0, 100.0, 200.0] }),
      ).toThrow(TypeError);
    });

    it("rejects ak/bk length mismatch", () => {
      expect(() =>
        grids.vertical({
          coordinate: "eta",
          ak: [0.0, 1000.0, 100.0],
          bk: [1.0, 0.5, 0.25, 0.0],
        }),
      ).toThrow(RangeError);
    });

    it("rejects non-decreasing synthesized sigma", () => {
      expect(() =>
        grids.vertical({
          coordinate: "eta",
          ak: [0.0, 1000.0, 500.0, 100.0],
          bk: [1.0, 0.5, 0.8, 0.0],
        }),
      ).toThrow(RangeError);
    });
  });

  // ---------------------------------------------- hybrid_sigma_theta
  describe("hybrid_sigma_theta coordinate", () => {
    it("defaults to a uniform sigma-like column from nz", () => {
      const g = grids.vertical({ coordinate: "hybrid_sigma_theta", nz: 5 });
      expect(g.nz).toBe(5);
      expect(closeRel(g.levels[0], 1.0)).toBe(true);
      expect(closeRel(g.levels[5], 0.0)).toBe(true);
    });

    it("accepts a transition in (0, 1)", () => {
      const g = grids.vertical({
        coordinate: "hybrid_sigma_theta",
        nz: 4,
        transition: 0.4,
      });
      expect(g.transition).toBe(0.4);
    });

    it("rejects transition out of (0, 1)", () => {
      expect(() =>
        grids.vertical({
          coordinate: "hybrid_sigma_theta",
          nz: 4,
          transition: 1.5,
        }),
      ).toThrow(RangeError);
      expect(() =>
        grids.vertical({
          coordinate: "hybrid_sigma_theta",
          nz: 4,
          transition: 0.0,
        }),
      ).toThrow(RangeError);
    });

    it("accepts optional ak / bk", () => {
      const levels = [1.0, 0.75, 0.5, 0.25, 0.0];
      const ak = [0.0, 500.0, 2000.0, 800.0, 100.0];
      const bk = [1.0, 0.7, 0.4, 0.2, 0.0];
      const g = grids.vertical({
        coordinate: "hybrid_sigma_theta",
        levels,
        ak,
        bk,
        p0: 1.0e5,
      });
      expect(g.ak.length).toBe(5);
      expect(g.bk.length).toBe(5);
    });
  });

  // --------------------------------------------------------- accessors
  describe("accessors", () => {
    it("cell_centers: bulk and scalar", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0, 700.0],
      });
      const bulk = g.cell_centers() as number[];
      expect(bulk.length).toBe(3);
      const expected = [50.0, 200.0, 500.0];
      for (let k = 0; k < 3; k++) {
        expect(closeRel(g.cell_centers(k) as number, expected[k])).toBe(true);
      }
    });

    it("cell_widths: bulk and scalar", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0, 700.0],
      });
      const bulk = g.cell_widths() as number[];
      expect(bulk.length).toBe(3);
      const expected = [100.0, 200.0, 400.0];
      for (let k = 0; k < 3; k++) {
        expect(closeRel(g.cell_widths(k) as number, expected[k])).toBe(true);
      }
    });

    it("cell accessors reject out-of-range layer", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(() => g.cell_centers(4)).toThrow(RangeError);
      expect(() => g.cell_centers(-1)).toThrow(RangeError);
      expect(() => g.cell_widths(10)).toThrow(RangeError);
    });

    it("neighbors at boundaries and interior", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 5 });
      expect(g.neighbors(0)).toEqual({ up: 1 });
      expect(g.neighbors(4)).toEqual({ down: 3 });
      expect(g.neighbors(2)).toEqual({ down: 1, up: 3 });
    });

    it("neighbors rejects out-of-range", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(() => g.neighbors(4)).toThrow(RangeError);
      expect(() => g.neighbors(-1)).toThrow(RangeError);
    });

    it("exposes n_cells / n_vertices / n_edges", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 7 });
      expect(g.n_cells).toBe(7);
      expect(g.n_vertices).toBe(8);
      expect(g.n_edges).toBe(7);
    });
  });

  // -------------------------------------------------------- metric_eval
  describe("metric_eval", () => {
    it("dz / z match cell accessors", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0, 700.0],
      });
      for (let k = 0; k < g.nz; k++) {
        expect(
          closeRel(g.metric_eval("dz", k), g.cell_widths(k) as number),
        ).toBe(true);
        expect(
          closeRel(g.metric_eval("z", k), g.cell_centers(k) as number),
        ).toBe(true);
      }
    });

    it("sigma valid for sigma-like coordinates", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      for (let k = 0; k < g.nz; k++) {
        expect(
          closeRel(g.metric_eval("sigma", k), g.cell_centers(k) as number),
        ).toBe(true);
      }
    });

    it("sigma rejected for z", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0],
      });
      expect(() => g.metric_eval("sigma", 0)).toThrow(RangeError);
    });

    it("pressure averages interface values for eta", () => {
      const p0 = 1.0e5;
      const ak = [0.0, 1000.0, 5000.0, 100.0];
      const bk = [1.0, 0.7, 0.3, 0.0];
      const g = grids.vertical({ coordinate: "eta", ak, bk, p0 });
      for (let k = 0; k < g.nz; k++) {
        const pLo = ak[k] + bk[k] * p0;
        const pHi = ak[k + 1] + bk[k + 1] * p0;
        expect(
          closeRel(g.metric_eval("pressure", k), 0.5 * (pLo + pHi), 1e-12),
        ).toBe(true);
      }
    });

    it("pressure requires hybrid coefficients", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(() => g.metric_eval("pressure", 0)).toThrow(RangeError);
    });

    it("ak / bk averages across layer interfaces", () => {
      const ak = [0.0, 1000.0, 5000.0, 100.0];
      const bk = [1.0, 0.7, 0.3, 0.0];
      const g = grids.vertical({ coordinate: "eta", ak, bk });
      for (let k = 0; k < g.nz; k++) {
        expect(
          closeRel(g.metric_eval("ak", k), 0.5 * (ak[k] + ak[k + 1])),
        ).toBe(true);
        expect(
          closeRel(g.metric_eval("bk", k), 0.5 * (bk[k] + bk[k + 1])),
        ).toBe(true);
      }
    });

    it("rejects unknown metric name", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(() =>
        // @ts-expect-error - outside the literal union
        g.metric_eval("not_a_metric", 0),
      ).toThrow(RangeError);
    });

    it("rejects out-of-range layer", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      expect(() => g.metric_eval("dz", 4)).toThrow(RangeError);
    });
  });

  // --------------------------------------------------------------- toESM
  describe("toESM (declarative lowering)", () => {
    it("emits schema-shaped declarative config for sigma", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4, ghosts: 2 });
      const doc = g.toESM() as Record<string, unknown>;
      expect(doc.family).toBe("vertical");
      expect(doc.topology).toBe("column");
      expect(doc.ndim).toBe(1);
      expect(doc.dtype).toBe("float64");
      expect(doc.ghosts).toBe(2);
      expect(doc.n_cells).toBe(4);
      expect(doc.n_vertices).toBe(5);
      const options = doc.options as Record<string, unknown>;
      expect(options.coordinate).toBe("sigma");
      expect(options.nz).toBe(4);
      expect((options.levels as number[]).length).toBe(5);
      expect("ak" in options).toBe(false);
      expect("p0" in options).toBe(false);
      expect(doc.provenance).toBeTruthy();
    });

    it("includes hybrid coeffs for eta", () => {
      const ak = [0.0, 1000.0, 5000.0, 100.0];
      const bk = [1.0, 0.7, 0.3, 0.0];
      const g = grids.vertical({ coordinate: "eta", ak, bk, p0: 1.0e5 });
      const doc = g.toESM() as Record<string, unknown>;
      const options = doc.options as Record<string, unknown>;
      expect(options.coordinate).toBe("eta");
      expect(options.ak).toEqual(ak);
      expect(options.bk).toEqual(bk);
      expect(options.p0).toBe(1.0e5);
    });

    it("omits derived arrays from wire form", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0, 700.0],
      });
      const wire = JSON.stringify(g.toESM());
      expect(wire.includes("centers")).toBe(false);
      expect(wire.includes("widths")).toBe(false);
    });

    it("JSON round-trips", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 5 });
      const parsed = JSON.parse(JSON.stringify(g.toESM())) as Record<
        string,
        unknown
      >;
      expect((parsed.options as Record<string, unknown>).nz).toBe(5);
      expect(parsed.family).toBe("vertical");
    });

    it("provenance carries typescript binding identity", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 4 });
      const doc = g.toESM() as Record<string, unknown>;
      const prov = doc.provenance as Record<string, unknown>;
      expect(prov.binding).toBe("typescript");
      expect(prov.family).toBe("vertical");
      expect(prov.coordinate).toBe("sigma");
    });

    it("matches the python-emitted sigma_uniform_n16 fixture shape", () => {
      const g = grids.vertical({ coordinate: "sigma", nz: 16 });
      const doc = g.toESM() as Record<string, unknown>;
      expect(doc.n_cells).toBe(16);
      expect(doc.n_vertices).toBe(17);
      expect(doc.schema_version).toBe("1.0.0");
      const options = doc.options as Record<string, unknown>;
      const levels = options.levels as number[];
      expect(levels.length).toBe(17);
      expect(closeRel(levels[0], 1.0)).toBe(true);
      expect(closeRel(levels[16], 0.0)).toBe(true);
      const expected: number[] = [];
      for (let k = 0; k <= 16; k++) expected.push(1 - k / 16);
      allClose(levels, expected);
    });
  });

  // --------------------------------------------------------------- dtype
  describe("dtype propagation", () => {
    it("float32 flows to grid + provenance", () => {
      const g = grids.vertical({
        coordinate: "z",
        levels: [0.0, 100.0, 300.0],
        dtype: "float32",
      });
      expect(g.dtype).toBe("float32");
      const doc = g.toESM() as Record<string, unknown>;
      expect(doc.dtype).toBe("float32");
      expect((doc.provenance as Record<string, unknown>).dtype).toBe("float32");
    });
  });
});
