import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type {
  ArakawaGrid,
  ArakawaStagger,
} from "../src/grids/index.js";

const REL_TOL = 1e-14;

function closeRel(a: number, b: number, tol = REL_TOL): boolean {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= tol * scale;
}

function closePair(
  p: readonly [number, number],
  q: readonly [number, number],
  tol = REL_TOL,
): boolean {
  return closeRel(p[0], q[0], tol) && closeRel(p[1], q[1], tol);
}

function baseCartesian(nx = 4, ny = 4) {
  return grids.cartesianBase({
    xlo: 0,
    xhi: 1,
    ylo: 0,
    yhi: 1,
    nx,
    ny,
  });
}

function makeGrid(stagger: ArakawaStagger, nx = 4, ny = 4): ArakawaGrid {
  return grids.arakawa({ base: baseCartesian(nx, ny), stagger });
}

describe("earthsci.grids.arakawa", () => {
  it("is exposed via both grids and earthsci.grids namespaces", () => {
    expect(typeof grids.arakawa).toBe("function");
    expect(earthsci.grids.arakawa).toBe(grids.arakawa);
    expect(typeof grids.cartesianBase).toBe("function");
  });

  describe("cartesianBase validation", () => {
    it("rejects non-positive nx", () => {
      expect(() =>
        grids.cartesianBase({
          xlo: 0,
          xhi: 1,
          ylo: 0,
          yhi: 1,
          nx: 0,
          ny: 4,
        }),
      ).toThrow(RangeError);
    });

    it("rejects reversed extent", () => {
      expect(() =>
        grids.cartesianBase({
          xlo: 1,
          xhi: 0,
          ylo: 0,
          yhi: 1,
          nx: 4,
          ny: 4,
        }),
      ).toThrow(RangeError);
    });

    it("rejects non-finite extent", () => {
      expect(() =>
        grids.cartesianBase({
          xlo: Number.NaN,
          xhi: 1,
          ylo: 0,
          yhi: 1,
          nx: 4,
          ny: 4,
        }),
      ).toThrow(RangeError);
    });

    it("exposes declared properties", () => {
      const b = grids.cartesianBase({
        xlo: -1,
        xhi: 3,
        ylo: 0,
        yhi: 2,
        nx: 4,
        ny: 2,
      });
      expect(b.family).toBe("cartesian");
      expect(b.nx).toBe(4);
      expect(b.ny).toBe(2);
      expect(b.xlo).toBe(-1);
      expect(b.xhi).toBe(3);
      expect(b.dx()).toBe(1);
      expect(b.dy()).toBe(1);
    });
  });

  describe("option contract", () => {
    it("rejects missing base", () => {
      // @ts-expect-error — validating runtime error path
      expect(() => grids.arakawa({ stagger: "C" })).toThrow(TypeError);
    });

    it("rejects missing stagger", () => {
      // @ts-expect-error — validating runtime error path
      expect(() => grids.arakawa({ base: baseCartesian() })).toThrow(
        TypeError,
      );
    });

    it("rejects invalid stagger label", () => {
      expect(() =>
        grids.arakawa({
          base: baseCartesian(),
          // @ts-expect-error — validating runtime error path
          stagger: "Q",
        }),
      ).toThrow(RangeError);
    });

    it("rejects invalid dtype", () => {
      expect(() =>
        grids.arakawa({
          base: baseCartesian(),
          stagger: "C",
          // @ts-expect-error — validating runtime error path
          dtype: "float16",
        }),
      ).toThrow(RangeError);
    });

    it("rejects negative ghosts", () => {
      expect(() =>
        grids.arakawa({
          base: baseCartesian(),
          stagger: "C",
          ghosts: -1,
        }),
      ).toThrow(RangeError);
    });

    it("rejects non-integer ghosts", () => {
      expect(() =>
        grids.arakawa({
          base: baseCartesian(),
          stagger: "C",
          ghosts: 1.5,
        }),
      ).toThrow(RangeError);
    });
  });

  describe("stagger → location table", () => {
    it("matches the API spec", () => {
      expect(grids.variableLocations("A")).toEqual([
        "cell_center",
        "cell_center",
        "cell_center",
      ]);
      expect(grids.variableLocations("B")).toEqual([
        "cell_center",
        "corner",
        "corner",
      ]);
      expect(grids.variableLocations("C")).toEqual([
        "cell_center",
        "u_edge",
        "v_edge",
      ]);
      expect(grids.variableLocations("D")).toEqual([
        "cell_center",
        "v_edge",
        "u_edge",
      ]);
      // E: topologically like B
      expect(grids.variableLocations("E")).toEqual([
        "cell_center",
        "corner",
        "corner",
      ]);
    });
  });

  describe("locationShape", () => {
    it("produces the standard shapes", () => {
      expect(grids.locationShape("cell_center", 10, 20)).toEqual([10, 20]);
      expect(grids.locationShape("u_edge", 10, 20)).toEqual([11, 20]);
      expect(grids.locationShape("v_edge", 10, 20)).toEqual([10, 21]);
      expect(grids.locationShape("corner", 10, 20)).toEqual([11, 21]);
    });
  });

  describe("defaults", () => {
    it("match the contract", () => {
      const gr = makeGrid("C");
      expect(gr.family).toBe("arakawa");
      expect(gr.topology).toBe("block_structured");
      expect(gr.dtype).toBe("float64");
      expect(gr.ghosts).toBe(0);
      expect(gr.nx).toBe(4);
      expect(gr.ny).toBe(4);
      expect(gr.n_cells).toBe(16);
      expect(gr.stagger).toBe("C");
      expect(gr.rotated).toBe(false);
    });
  });

  describe("C-grid accessors recover Cartesian coordinates", () => {
    it("exactly reproduces centers / faces / corners", () => {
      const gr = grids.arakawa({
        base: grids.cartesianBase({
          xlo: 0,
          xhi: 1,
          ylo: 0,
          yhi: 1,
          nx: 10,
          ny: 10,
        }),
        stagger: "C",
      });

      expect(closePair(gr.cell_centers(0, 0), [0.05, 0.05])).toBe(true);
      expect(closePair(gr.cell_centers(9, 9), [0.95, 0.95])).toBe(true);

      // u at u-faces: (0, 0) is the western boundary at mid-y of row 0.
      expect(closePair(gr.u_face(0, 0), [0.0, 0.05])).toBe(true);
      // u(10, 0) is the eastern boundary of the last column.
      expect(closeRel(gr.u_face(10, 0)[0], 1.0)).toBe(true);

      // v at v-faces.
      expect(closePair(gr.v_face(0, 0), [0.05, 0.0])).toBe(true);
      expect(closeRel(gr.v_face(0, 10)[1], 1.0)).toBe(true);

      // Corners.
      expect(gr.corners(0, 0)).toEqual([0, 0]);
      expect(gr.corners(10, 10)).toEqual([1, 1]);
    });
  });

  describe("A-grid: colocation", () => {
    it("puts h, u, v all at cell centers", () => {
      const gr = makeGrid("A");
      expect(gr.variable_shape("h")).toEqual([4, 4]);
      expect(gr.variable_shape("u")).toEqual([4, 4]);
      expect(gr.variable_shape("v")).toEqual([4, 4]);
      expect(gr.u_face(1, 2)).toEqual(gr.cell_centers(1, 2));
      expect(gr.v_face(1, 2)).toEqual(gr.cell_centers(1, 2));
    });
  });

  describe("B-grid: u, v at corners", () => {
    it("matches the corner shape", () => {
      const gr = makeGrid("B");
      expect(gr.variable_shape("h")).toEqual([4, 4]);
      expect(gr.variable_shape("u")).toEqual([5, 5]);
      expect(gr.variable_shape("v")).toEqual([5, 5]);
      expect(gr.u_face(0, 0)).toEqual(gr.corners(0, 0));
      expect(gr.v_face(2, 1)).toEqual(gr.corners(2, 1));
    });
  });

  describe("D-grid: swaps C's u/v faces", () => {
    it("swaps shapes and locations vs. C", () => {
      const c = grids.arakawa({
        base: grids.cartesianBase({
          xlo: 0,
          xhi: 2,
          ylo: 0,
          yhi: 1,
          nx: 4,
          ny: 4,
        }),
        stagger: "C",
      });
      const d = grids.arakawa({
        base: grids.cartesianBase({
          xlo: 0,
          xhi: 2,
          ylo: 0,
          yhi: 1,
          nx: 4,
          ny: 4,
        }),
        stagger: "D",
      });
      expect(d.variable_shape("u")).toEqual(c.variable_shape("v"));
      expect(d.variable_shape("v")).toEqual(c.variable_shape("u"));
      expect(d.u_face(0, 0)).toEqual(c.v_face(0, 0));
      expect(d.v_face(0, 0)).toEqual(c.u_face(0, 0));
    });
  });

  describe("E-grid: rotated B", () => {
    it("shares topology with B and flags rotated", () => {
      const e = makeGrid("E");
      const b = makeGrid("B");
      expect(e.rotated).toBe(true);
      expect(b.rotated).toBe(false);
      expect(e.variable_shape("h")).toEqual(b.variable_shape("h"));
      expect(e.variable_shape("u")).toEqual(b.variable_shape("u"));
      expect(e.variable_shape("v")).toEqual(b.variable_shape("v"));
    });
  });

  describe("neighbors", () => {
    it("returns all four for an interior cell", () => {
      const gr = makeGrid("C");
      const nb = gr.neighbors("cell_center", 1, 1);
      expect(nb.west).toEqual([0, 1]);
      expect(nb.east).toEqual([2, 1]);
      expect(nb.south).toEqual([1, 0]);
      expect(nb.north).toEqual([1, 2]);
    });

    it("returns null on domain boundaries", () => {
      const gr = makeGrid("C");
      const nb = gr.neighbors("cell_center", 0, 0);
      expect(nb.west).toBeNull();
      expect(nb.south).toBeNull();
      expect(nb.east).toEqual([1, 0]);
      expect(nb.north).toEqual([0, 1]);
    });

    it("respects UEdge shape on the last column", () => {
      // UEdge shape is (5, 4); (4, 3) is the top-right corner of UEdge.
      const gr = makeGrid("C");
      const nb = gr.neighbors("u_edge", 4, 3);
      expect(nb.east).toBeNull();
      expect(nb.north).toBeNull();
      expect(nb.west).toEqual([3, 3]);
      expect(nb.south).toEqual([4, 2]);
    });

    it("rejects out-of-range indices", () => {
      const gr = makeGrid("C");
      expect(() => gr.neighbors("cell_center", 4, 0)).toThrow(RangeError);
      expect(() => gr.neighbors("cell_center", 0, 4)).toThrow(RangeError);
    });
  });

  describe("metric_eval", () => {
    it("returns uniform Cartesian metrics", () => {
      const gr = grids.arakawa({
        base: grids.cartesianBase({
          xlo: 0,
          xhi: 2,
          ylo: 0,
          yhi: 6,
          nx: 4,
          ny: 3,
        }),
        stagger: "C",
      });
      expect(closeRel(gr.metric_eval("dx", 0, 0), 0.5)).toBe(true);
      expect(closeRel(gr.metric_eval("dy", 0, 0), 2.0)).toBe(true);
      expect(closeRel(gr.metric_eval("area", 3, 2), 1.0)).toBe(true);
    });

    it("rejects unknown metric name", () => {
      const gr = makeGrid("C");
      expect(() =>
        // @ts-expect-error — validating runtime error path
        gr.metric_eval("not_a_metric", 0, 0),
      ).toThrow(RangeError);
    });

    it("rejects out-of-range cell", () => {
      const gr = makeGrid("C");
      expect(() => gr.metric_eval("dx", 4, 0)).toThrow(RangeError);
    });
  });

  describe("bounds checks on accessors", () => {
    it("rejects out-of-range centers/corners", () => {
      const gr = makeGrid("C");
      expect(() => gr.cell_centers(4, 0)).toThrow(RangeError);
      expect(() => gr.cell_centers(-1, 0)).toThrow(RangeError);
      expect(() => gr.corners(5, 0)).toThrow(RangeError);
    });

    it("rejects non-integer indices", () => {
      const gr = makeGrid("C");
      expect(() => gr.cell_centers(0.5, 0)).toThrow(RangeError);
    });
  });

  describe("toESM: declarative config", () => {
    it("emits a small declarative .esm without inline geometry", () => {
      const gr = grids.arakawa({
        base: grids.cartesianBase({
          xlo: 0,
          xhi: 1,
          ylo: 0,
          yhi: 1,
          nx: 5,
          ny: 5,
        }),
        stagger: "C",
        ghosts: 2,
      });
      const doc = gr.toESM() as Record<string, unknown>;
      expect(doc.family).toBe("arakawa");
      expect(doc.version).toBe("1.0.0");
      expect(doc.dtype).toBe("float64");
      expect(doc.topology).toBe("block_structured");
      expect(doc.stagger).toBe("C");
      expect(doc.ghosts).toBe(2);
      expect(doc.n_cells).toBe(25);
      expect(doc.rotated).toBe(false);

      const baseDoc = doc.base as Record<string, unknown>;
      expect(baseDoc.family).toBe("cartesian");
      expect(baseDoc.nx).toBe(5);
      expect(baseDoc.ny).toBe(5);

      // Per §0 correction: no inline geometry arrays.
      const wire = JSON.stringify(doc);
      expect(wire).not.toContain('"cells":[');
      expect(wire).not.toContain('"edges":[');
      expect(wire).not.toContain('"vertices":[');
      expect(wire.length).toBeLessThan(2_000);
    });

    it("sets rotated=true for E and false for B", () => {
      const e = makeGrid("E");
      const b = makeGrid("B");
      expect((e.toESM() as Record<string, unknown>).rotated).toBe(true);
      expect((b.toESM() as Record<string, unknown>).rotated).toBe(false);
    });

    it("propagates dtype=float32", () => {
      const gr = grids.arakawa({
        base: baseCartesian(),
        stagger: "C",
        dtype: "float32",
      });
      expect(gr.dtype).toBe("float32");
      expect((gr.toESM() as Record<string, unknown>).dtype).toBe("float32");
    });

    it("identifies the TypeScript binding in provenance", () => {
      const gr = makeGrid("C");
      const doc = gr.toESM() as Record<string, unknown>;
      const prov = doc.provenance as Record<string, unknown>;
      expect(prov.binding).toBe("typescript");
      expect(prov.api_version).toBe("1.0.0");
      expect(prov.stagger).toBe("C");
    });
  });

  describe("custom base grids plug in via the interface", () => {
    it("accepts any object that satisfies ArakawaBaseGrid", () => {
      // A degenerate 2-cell × 2-cell base that returns constants; enough to
      // exercise that arakawa doesn't depend on CartesianBase specifically.
      const fake = {
        nx: 2,
        ny: 2,
        dx: () => 7,
        dy: () => 11,
        cellCenter: (): [number, number] => [0, 0],
        xEdge: (): [number, number] => [1, 0],
        yEdge: (): [number, number] => [0, 1],
        corner: (): [number, number] => [2, 2],
        toESM: () => ({ family: "fake", nx: 2, ny: 2 }),
      };
      const gr = grids.arakawa({ base: fake, stagger: "C" });
      expect(gr.n_cells).toBe(4);
      expect(gr.metric_eval("dx", 0, 0)).toBe(7);
      expect(gr.metric_eval("dy", 0, 0)).toBe(11);
      expect(gr.metric_eval("area", 0, 0)).toBe(77);
      expect(gr.corners(2, 2)).toEqual([2, 2]);
    });
  });
});
