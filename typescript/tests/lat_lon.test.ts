import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { LatLonGrid } from "../src/grids/index.js";

const PI = Math.PI;
const HALF_PI = PI / 2;

function regular(nlon: number, nlat: number, R = 6.371e6): LatLonGrid {
  return grids.lat_lon({ nlon, nlat, R });
}

describe("earthsci.grids.lat_lon — surface", () => {
  it("is exposed via both grids and earthsci.grids namespaces", () => {
    expect(typeof grids.lat_lon).toBe("function");
    expect(earthsci.grids.lat_lon).toBe(grids.lat_lon);
  });

  it("returns a Grid with the contracted shape (regular defaults)", () => {
    const g = regular(8, 4);
    expect(g.family).toBe("lat_lon");
    expect(g.topology).toBe("rectilinear");
    expect(g.variant).toBe("regular");
    expect(g.dtype).toBe("float64");
    expect(g.nlat).toBe(4);
    expect(g.nlon_uniform).toBe(8);
    expect(g.nlon_per_row).toEqual([8, 8, 8, 8]);
    expect(g.R).toBe(6.371e6);
    expect(g.ghosts).toBe(0);
    expect(g.pole_policy).toBe("none");
    expect(g.lon_start).toBeCloseTo(-PI, 15);
    expect(g.n_cells).toBe(8 * 4);
  });

  it("honors option overrides (R, dtype, ghosts, lon_start)", () => {
    const g = grids.lat_lon({
      nlon: 6,
      nlat: 3,
      R: 1,
      dtype: "float32",
      ghosts: 2,
      lon_start: 0,
    });
    expect(g.R).toBe(1);
    expect(g.dtype).toBe("float32");
    expect(g.ghosts).toBe(2);
    expect(g.lon_start).toBe(0);
  });
});

describe("earthsci.grids.lat_lon — error contract (§9)", () => {
  it("rejects missing required options for regular variant", () => {
    expect(() => grids.lat_lon({} as unknown as { nlon: 1; nlat: 1 })).toThrow(
      TypeError,
    );
    expect(() =>
      grids.lat_lon({ nlat: 4 } as unknown as { nlon: 1; nlat: 1 }),
    ).toThrow(TypeError);
    expect(() =>
      grids.lat_lon({ nlon: 4 } as unknown as { nlon: 1; nlat: 1 }),
    ).toThrow(TypeError);
  });

  it("rejects non-positive-integer nlon / nlat", () => {
    expect(() => grids.lat_lon({ nlon: 0, nlat: 4 })).toThrow(RangeError);
    expect(() => grids.lat_lon({ nlon: 4, nlat: 0 })).toThrow(RangeError);
    expect(() => grids.lat_lon({ nlon: 3.5, nlat: 4 })).toThrow(RangeError);
  });

  it("rejects invalid R / lon_start / ghosts / dtype", () => {
    expect(() => grids.lat_lon({ nlon: 4, nlat: 4, R: -1 })).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({ nlon: 4, nlat: 4, R: Number.NaN }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({ nlon: 4, nlat: 4, lon_start: Number.POSITIVE_INFINITY }),
    ).toThrow(RangeError);
    expect(() => grids.lat_lon({ nlon: 4, nlat: 4, ghosts: -1 })).toThrow(
      RangeError,
    );
    expect(() =>
      grids.lat_lon({
        nlon: 4,
        nlat: 4,
        dtype: "bogus" as "float64",
      }),
    ).toThrow(RangeError);
  });

  it("rejects unimplemented pole_policy values at build time", () => {
    expect(() =>
      grids.lat_lon({ nlon: 4, nlat: 4, pole_policy: "average" }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({ nlon: 4, nlat: 4, pole_policy: "fold" }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({
        nlon: 4,
        nlat: 4,
        pole_policy: "wrap" as "none",
      }),
    ).toThrow(RangeError);
  });

  it("rejects nlon_per_row in regular variant and nlon in reduced_gaussian", () => {
    expect(() =>
      grids.lat_lon({ nlon: 4, nlat: 2, nlon_per_row: [4, 4] }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({
        variant: "reduced_gaussian",
        nlon: 4,
        nlon_per_row: [4, 8],
      }),
    ).toThrow(RangeError);
  });

  it("requires nlon_per_row for reduced_gaussian and validates length vs nlat", () => {
    expect(() =>
      grids.lat_lon({ variant: "reduced_gaussian", nlat: 4 }),
    ).toThrow(TypeError);
    expect(() =>
      grids.lat_lon({
        variant: "reduced_gaussian",
        nlat: 3,
        nlon_per_row: [4, 8],
      }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({
        variant: "reduced_gaussian",
        nlon_per_row: [4, 0, 8],
      }),
    ).toThrow(RangeError);
  });

  it("rejects malformed lat_edges / lat_centers", () => {
    expect(() =>
      grids.lat_lon({ nlon: 2, nlat: 2, lat_edges: [0, 0, 0.5] }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({ nlon: 2, nlat: 1, lat_edges: [-2, 2] }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({ nlon: 2, nlat: 2, lat_edges: [-HALF_PI, 0] }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({
        nlon: 2,
        nlat: 2,
        lat_edges: [-HALF_PI, 0, HALF_PI],
        lat_centers: [-HALF_PI - 0.1, 1],
      }),
    ).toThrow(RangeError);
    expect(() =>
      grids.lat_lon({
        nlon: 2,
        nlat: 2,
        lat_edges: [-HALF_PI, 0, HALF_PI],
        lat_centers: [-1],
      }),
    ).toThrow(RangeError);
  });
});

describe("earthsci.grids.lat_lon — topology (regular)", () => {
  it("interior neighbors are local", () => {
    const g = regular(8, 6);
    const n = g.neighbors(2, 3);
    expect(n.W).toEqual({ j: 2, i: 2 });
    expect(n.E).toEqual({ j: 2, i: 4 });
    expect(n.S).toEqual({ j: 1, i: 3 });
    expect(n.N).toEqual({ j: 3, i: 3 });
  });

  it("longitude wraps periodically (W of i=0 and E of i=N-1)", () => {
    const g = regular(8, 4);
    expect(g.neighbors(2, 0).W).toEqual({ j: 2, i: 7 });
    expect(g.neighbors(2, 7).E).toEqual({ j: 2, i: 0 });
  });

  it("south pole row has no S neighbor; north pole row has no N neighbor", () => {
    const g = regular(8, 4);
    const south = g.neighbors(0, 3);
    expect(south.S).toBeNull();
    expect(south.N).toEqual({ j: 1, i: 3 });
    const north = g.neighbors(3, 3);
    expect(north.N).toBeNull();
    expect(north.S).toEqual({ j: 2, i: 3 });
  });

  it("every neighbor index is in range across the grid", () => {
    const g = regular(5, 4);
    for (let j = 0; j < g.nlat; j++) {
      for (let i = 0; i < g.nlon(j); i++) {
        const ns = g.neighbors(j, i);
        for (const c of [ns.W, ns.E, ns.S, ns.N]) {
          if (c === null) continue;
          expect(c.j).toBeGreaterThanOrEqual(0);
          expect(c.j).toBeLessThan(g.nlat);
          expect(c.i).toBeGreaterThanOrEqual(0);
          expect(c.i).toBeLessThan(g.nlon(c.j));
        }
      }
    }
  });

  it("longitudinal periodicity: the W of E neighbor of (j,i) is (j,i)", () => {
    const g = regular(6, 5);
    for (let j = 0; j < g.nlat; j++) {
      for (let i = 0; i < g.nlon(j); i++) {
        const e = g.neighbors(j, i).E;
        expect(g.neighbors(e.j, e.i).W).toEqual({ j, i });
      }
    }
  });

  it("cell index validation rejects out-of-range j and i", () => {
    const g = regular(4, 4);
    expect(() => g.cell_center(4, 0)).toThrow(RangeError);
    expect(() => g.cell_center(0, 4)).toThrow(RangeError);
    expect(() => g.neighbors(-1, 0)).toThrow(RangeError);
    expect(() => g.neighbor(0, 0, "X" as "W")).toThrow(RangeError);
  });
});

describe("earthsci.grids.lat_lon — topology (reduced_gaussian)", () => {
  it("rejects N/S neighbor at correct rounded column when widths differ", () => {
    const g = grids.lat_lon({
      variant: "reduced_gaussian",
      nlon_per_row: [4, 8, 4],
    });
    // Row 1 cell (1,5): center frac (5+0.5)/8 = 0.6875 -> floor(0.6875*4)=2 in row 0.
    expect(g.neighbors(1, 5).S).toEqual({ j: 0, i: 2 });
    // Row 1 cell (1,0): center frac 0.0625 -> floor(0.0625*4)=0 in row 0.
    expect(g.neighbors(1, 0).S).toEqual({ j: 0, i: 0 });
  });

  it("row offsets and total cell counts are consistent", () => {
    const g = grids.lat_lon({
      variant: "reduced_gaussian",
      nlon_per_row: [4, 8, 4],
    });
    expect(g.n_cells).toBe(16);
    expect(g.nlat).toBe(3);
    expect(g.nlon_uniform).toBeNull();
    expect(g.row_offset(0)).toBe(0);
    expect(g.row_offset(1)).toBe(4);
    expect(g.row_offset(2)).toBe(12);
    expect(g.row_offset(3)).toBe(16);
  });

  it("nlat defaults to nlon_per_row.length when omitted", () => {
    const g = grids.lat_lon({
      variant: "reduced_gaussian",
      nlon_per_row: [4, 8, 4, 4],
    });
    expect(g.nlat).toBe(4);
  });
});

describe("earthsci.grids.lat_lon — centers and edges (metric accuracy)", () => {
  it("default lat edges span -π/2 to π/2 with nlat+1 entries", () => {
    const g = regular(4, 6);
    expect(g.lat_edges.length).toBe(7);
    expect(g.lat_edges[0]).toBeCloseTo(-HALF_PI, 15);
    expect(g.lat_edges[6]).toBeCloseTo(HALF_PI, 15);
  });

  it("lon_edges and lon_centers tile the row uniformly", () => {
    const g = grids.lat_lon({ nlon: 4, nlat: 2, lon_start: 0 });
    expect(g.lon_edges(0)).toEqual([0, PI / 2, PI, (3 * PI) / 2, 2 * PI]);
    expect(g.lon_centers(0)).toEqual([
      PI / 4,
      (3 * PI) / 4,
      (5 * PI) / 4,
      (7 * PI) / 4,
    ]);
  });

  it("cell_center matches the cell_centers_bulk flatten layout", () => {
    const g = regular(6, 4);
    const { lon, lat } = g.cell_centers_bulk();
    expect(lon.length).toBe(g.n_cells);
    expect(lat.length).toBe(g.n_cells);
    for (const j of [0, 1, 3]) {
      for (const i of [0, 2, 5]) {
        const c = g.cell_center(j, i);
        const idx = g.row_offset(j) + i;
        expect(c.lon).toBeCloseTo(lon[idx], 14);
        expect(c.lat).toBeCloseTo(lat[idx], 14);
      }
    }
  });

  it("with even nlat the equatorial-edge sits on lat=0 (rectangle symmetry)", () => {
    const g = regular(4, 2);
    // edges = [-π/2, 0, π/2]; row 0 center = -π/4, row 1 center = π/4
    expect(g.cell_center(0, 0).lat).toBeCloseTo(-PI / 4, 15);
    expect(g.cell_center(1, 0).lat).toBeCloseTo(PI / 4, 15);
  });
});

describe("earthsci.grids.lat_lon — area accuracy", () => {
  it("sum of cell areas = 4πR² for a regular grid", () => {
    const R = 6.371e6;
    const g = grids.lat_lon({ nlon: 72, nlat: 36, R });
    const total = g.area_bulk().reduce((a, b) => a + b, 0);
    const expected = 4 * PI * R * R;
    expect(Math.abs(total - expected) / expected).toBeLessThan(1e-12);
  });

  it("sum of cell areas = 4π on a unit sphere", () => {
    const g = grids.lat_lon({ nlon: 10, nlat: 8, R: 1 });
    const total = g.area_bulk().reduce((a, b) => a + b, 0);
    expect(Math.abs(total - 4 * PI)).toBeLessThan(1e-12);
  });

  it("reduced-Gaussian total area also = 4πR²", () => {
    const g = grids.lat_lon({
      variant: "reduced_gaussian",
      nlon_per_row: [4, 8, 12, 12, 8, 4],
      R: 1,
    });
    const total = g.area_bulk().reduce((a, b) => a + b, 0);
    expect(Math.abs(total - 4 * PI)).toBeLessThan(1e-12);
  });

  it("every cell has positive area", () => {
    const g = regular(6, 4);
    for (const a of g.area_bulk()) {
      expect(a).toBeGreaterThan(0);
    }
  });

  it("scalar cell_area matches area_bulk per cell", () => {
    const g = regular(6, 4);
    const bulk = g.area_bulk();
    for (const j of [0, 1, 3]) {
      for (const i of [0, 2, 5]) {
        const idx = g.row_offset(j) + i;
        expect(g.cell_area(j, i)).toBeCloseTo(bulk[idx], 12);
      }
    }
  });
});

describe("earthsci.grids.lat_lon — metric tensor", () => {
  it("g_lonlat and ginv_lonlat are zero everywhere", () => {
    const g = regular(4, 4);
    for (let j = 0; j < g.nlat; j++) {
      expect(g.metric_eval("g_lonlat", j, 0)).toBe(0);
      expect(g.metric_eval("ginv_lonlat", j, 0)).toBe(0);
    }
  });

  it("g_latlat = R² (independent of lat)", () => {
    const R = 3;
    const g = grids.lat_lon({ nlon: 4, nlat: 5, R });
    for (let j = 0; j < g.nlat; j++) {
      expect(g.metric_eval("g_latlat", j, 0)).toBeCloseTo(R * R, 14);
    }
  });

  it("g_lonlon ≈ R² near the equator (large nlat)", () => {
    const R = 3;
    const g = grids.lat_lon({ nlon: 4, nlat: 400, R });
    const jEq = 200; // lat center close to 0
    expect(
      Math.abs(g.metric_eval("g_lonlon", jEq, 0) - R * R) / (R * R),
    ).toBeLessThan(1e-4);
  });

  it("inverse metric components are exact inverses", () => {
    const g = grids.lat_lon({ nlon: 6, nlat: 6, R: 1 });
    for (const j of [0, 2, 5]) {
      const gLL = g.metric_eval("g_lonlon", j, 0);
      const gPP = g.metric_eval("g_latlat", j, 0);
      expect(g.metric_eval("ginv_lonlon", j, 0) * gLL).toBeCloseTo(1, 12);
      expect(g.metric_eval("ginv_latlat", j, 0) * gPP).toBeCloseTo(1, 12);
    }
  });

  it("J = R²|cos(lat)| > 0 in the interior, area metric matches cell_area", () => {
    const g = regular(6, 6);
    for (let j = 1; j < g.nlat - 1; j++) {
      for (let i = 0; i < g.nlon(j); i++) {
        expect(g.metric_eval("J", j, i)).toBeGreaterThan(0);
        expect(g.metric_eval("area", j, i)).toBeCloseTo(g.cell_area(j, i), 12);
      }
    }
  });

  it("rejects unknown metric names and out-of-range cells", () => {
    const g = regular(4, 4);
    expect(() =>
      g.metric_eval("not_a_metric" as "J", 0, 0),
    ).toThrow(RangeError);
    expect(() => g.metric_eval("J", 4, 0)).toThrow(RangeError);
    expect(() => g.metric_eval("J", 0, 4)).toThrow(RangeError);
  });
});

describe("earthsci.grids.lat_lon — toESM / declarative lowering", () => {
  it("regular variant emits scalar dimensions only (no geometry blob)", () => {
    const g = grids.lat_lon({ nlon: 72, nlat: 36, R: 6.371e6, ghosts: 2 });
    const doc = g.toESM() as Record<string, unknown>;
    expect(doc.family).toBe("lat_lon");
    expect(doc.version).toBe("1.0.0");
    expect(doc.dtype).toBe("float64");
    expect(doc.topology).toBe("rectilinear");
    expect(doc.variant).toBe("regular");
    expect(doc.generator).toBe("lat_lon_regular");
    const params = doc.params as Record<string, unknown>;
    expect(params.nlon).toBe(72);
    expect(params.nlat).toBe(36);
    expect(params.R).toBe(6.371e6);
    expect(params.ghosts).toBe(2);
    expect(params.pole_policy).toBe("none");
    expect(doc.provenance).toBeDefined();

    const wire = JSON.stringify(doc);
    expect(wire).not.toContain("cells");
    expect(wire).not.toContain("lon_array");
    expect(wire.length).toBeLessThan(2000);
  });

  it("reduced_gaussian variant carries the row-width schedule and lat edges", () => {
    const g = grids.lat_lon({
      variant: "reduced_gaussian",
      nlon_per_row: [4, 8, 4],
      R: 1,
    });
    const doc = g.toESM() as Record<string, unknown>;
    expect(doc.variant).toBe("reduced_gaussian");
    expect(doc.generator).toBe("lat_lon_reduced_gaussian");
    const params = doc.params as Record<string, unknown>;
    expect(params.nlat).toBe(3);
    expect(params.nlon_per_row).toEqual([4, 8, 4]);
    expect(Array.isArray(params.lat_edges)).toBe(true);
    expect((params.lat_edges as number[]).length).toBe(4);
  });

  it("toESM round-trips through JSON.stringify / JSON.parse", () => {
    const g = grids.lat_lon({ nlon: 4, nlat: 3, R: 1 });
    const text = JSON.stringify(g.toESM());
    const reparsed = JSON.parse(text) as Record<string, unknown>;
    const params = reparsed.params as Record<string, unknown>;
    expect(params.nlon).toBe(4);
    expect(params.nlat).toBe(3);
    expect(params.R).toBe(1);
  });

  it("provenance identifies the typescript binding", () => {
    const g = regular(4, 4);
    const doc = g.toESM() as Record<string, unknown>;
    const prov = doc.provenance as Record<string, unknown>;
    expect(prov.binding).toBe("typescript");
    expect(prov.api_version).toBe("1.0.0");
    expect(prov.generator).toBe("lat_lon_regular");
  });

  it("dtype=float32 propagates to toESM", () => {
    const g = grids.lat_lon({ nlon: 4, nlat: 4, dtype: "float32" });
    expect(g.dtype).toBe("float32");
    expect((g.toESM() as Record<string, unknown>).dtype).toBe("float32");
  });
});
