/**
 * Cubed-sphere grid (gnomonic, 6 panels × Nc × Nc) — accessor runtime.
 *
 * Per `docs/GRIDS_API.md` §2.6 and the 2026-04-20 mayor correction in bead
 * dsc-suu: `toESM()` emits a small declarative config, not a serialized
 * geometry blob. Geometry is derived on demand via the accessors
 * (`cell_centers`, `neighbors`, `metric_eval`, `cell_area`).
 *
 * Panels are 0-indexed internally (0..5) mapping to the standard FV3
 * numbering: 0=front, 1=right, 2=top (north pole), 3=back, 4=left,
 * 5=bottom (south pole). Cell indices i,j are 0-indexed in [0, Nc).
 */

import type { Dtype, Grid } from "./types.js";

export interface CubedSphereOpts {
  Nc: number;
  R?: number;
  dtype?: Dtype;
  ghosts?: number;
}

export type Edge = "west" | "east" | "south" | "north";

export interface NeighborRef {
  panel: number;
  i: number;
  j: number;
  edge: Edge;
}

export type MetricName =
  | "J"
  | "g_xixi"
  | "g_etaeta"
  | "g_xieta"
  | "ginv_xixi"
  | "ginv_etaeta"
  | "ginv_xieta";

export interface CubedSphereGrid extends Grid {
  readonly family: "cubed_sphere";
  readonly topology: "block_structured";
  readonly generator: "gnomonic_c6";
  readonly Nc: number;
  readonly R: number;
  readonly ghosts: number;
  readonly n_panels: 6;
  readonly n_cells: number;

  cell_centers(panel: number, i: number, j: number): { lon: number; lat: number };
  neighbors(panel: number, i: number, j: number): NeighborRef[];
  metric_eval(name: MetricName, panel: number, i: number, j: number): number;
  cell_area(panel: number, i: number, j: number): number;
}

interface PanelNeighbor {
  panel: number;
  edge: Edge;
  reverseIndex: boolean;
}

type ConnectivityRow = Record<Edge, PanelNeighbor>;

// Panel connectivity — mirrors src/grids/panel_connectivity.jl (0-indexed).
// Julia panels 1..6 → TS panels 0..5.
const PANEL_CONNECTIVITY: ReadonlyArray<ConnectivityRow> = [
  // 0 (front)
  {
    west: { panel: 4, edge: "east", reverseIndex: false },
    east: { panel: 1, edge: "west", reverseIndex: false },
    south: { panel: 5, edge: "north", reverseIndex: false },
    north: { panel: 2, edge: "south", reverseIndex: false },
  },
  // 1 (right)
  {
    west: { panel: 0, edge: "east", reverseIndex: false },
    east: { panel: 3, edge: "west", reverseIndex: false },
    south: { panel: 5, edge: "east", reverseIndex: true },
    north: { panel: 2, edge: "east", reverseIndex: false },
  },
  // 2 (top / north pole)
  {
    west: { panel: 4, edge: "north", reverseIndex: true },
    east: { panel: 1, edge: "north", reverseIndex: false },
    south: { panel: 0, edge: "north", reverseIndex: false },
    north: { panel: 3, edge: "north", reverseIndex: true },
  },
  // 3 (back)
  {
    west: { panel: 1, edge: "east", reverseIndex: false },
    east: { panel: 4, edge: "west", reverseIndex: false },
    south: { panel: 5, edge: "south", reverseIndex: true },
    north: { panel: 2, edge: "north", reverseIndex: true },
  },
  // 4 (left)
  {
    west: { panel: 3, edge: "east", reverseIndex: false },
    east: { panel: 0, edge: "west", reverseIndex: false },
    south: { panel: 5, edge: "west", reverseIndex: false },
    north: { panel: 2, edge: "west", reverseIndex: true },
  },
  // 5 (bottom / south pole)
  {
    west: { panel: 4, edge: "south", reverseIndex: false },
    east: { panel: 1, edge: "south", reverseIndex: true },
    south: { panel: 3, edge: "south", reverseIndex: true },
    north: { panel: 0, edge: "south", reverseIndex: false },
  },
];

function gnomonicToCart(
  xi: number,
  eta: number,
  panel: number,
): [number, number, number] {
  const X = Math.tan(xi);
  const Y = Math.tan(eta);
  const D = Math.sqrt(1 + X * X + Y * Y);
  switch (panel) {
    case 0:
      return [1 / D, X / D, Y / D];
    case 1:
      return [-X / D, 1 / D, Y / D];
    case 2:
      return [-Y / D, X / D, 1 / D];
    case 3:
      return [-1 / D, -X / D, Y / D];
    case 4:
      return [X / D, -1 / D, Y / D];
    case 5:
      return [Y / D, X / D, -1 / D];
    default:
      throw new RangeError(`cubed_sphere: invalid panel ${panel} (expected 0..5)`);
  }
}

function gnomonicToLonLat(
  xi: number,
  eta: number,
  panel: number,
): { lon: number; lat: number } {
  const [x, y, z] = gnomonicToCart(xi, eta, panel);
  const clamped = z < -1 ? -1 : z > 1 ? 1 : z;
  return { lon: Math.atan2(y, x), lat: Math.asin(clamped) };
}

function gnomonicMetric(
  xi: number,
  eta: number,
  R: number,
): { J: number; g_xixi: number; g_etaeta: number; g_xieta: number } {
  const X = Math.tan(xi);
  const Y = Math.tan(eta);
  const secXi2 = 1 + X * X;
  const secEta2 = 1 + Y * Y;
  const D2 = 1 + X * X + Y * Y;
  const D2_32 = Math.pow(D2, 1.5);
  const R2 = R * R;
  return {
    J: (R2 * secXi2 * secEta2) / D2_32,
    g_xixi: (R2 * secXi2 * secXi2 * secEta2) / (D2 * D2),
    g_etaeta: (R2 * secEta2 * secEta2 * secXi2) / (D2 * D2),
    g_xieta: (-R2 * X * Y * secXi2 * secEta2) / (D2 * D2),
  };
}

function cellAreaFromCorners(
  corners: ReadonlyArray<[number, number, number]>,
  R: number,
): number {
  // Spherical-excess via inner-angle sum at 4 corners on the unit sphere.
  const n = corners.length;
  let angleSum = 0;
  for (let k = 0; k < n; k++) {
    const vPrev = corners[(k - 1 + n) % n];
    const vCurr = corners[k];
    const vNext = corners[(k + 1) % n];
    const t1 = sideTangent(vCurr, vPrev);
    const t2 = sideTangent(vCurr, vNext);
    const n1 = norm(t1);
    const n2 = norm(t2);
    if (n1 < 1e-15 || n2 < 1e-15) {
      angleSum += Math.PI / 2;
    } else {
      const c = dot(t1, t2) / (n1 * n2);
      angleSum += Math.acos(c < -1 ? -1 : c > 1 ? 1 : c);
    }
  }
  return R * R * (angleSum - (n - 2) * Math.PI);
}

type Vec3 = [number, number, number];

function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function sub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function norm(a: Vec3): number {
  return Math.sqrt(dot(a, a));
}

// t = vCurr × ((vOther - vCurr) × vCurr): tangent at vCurr along the
// great-circle arc toward vOther (matches Julia compute_cell_area).
function sideTangent(vCurr: Vec3, vOther: Vec3): Vec3 {
  return cross(vCurr, cross(sub(vOther, vCurr), vCurr));
}

function transformGhostIndex(
  nb: PanelNeighbor,
  i: number,
  j: number,
  Nc: number,
  localEdge: Edge,
): { i: number; j: number } {
  // i,j are 0-indexed in the *local* panel; returns 0-indexed indices in
  // the neighbor panel for the ghost cell just across `localEdge`.
  //
  // "along" = coordinate that varies along the shared edge (0..Nc-1);
  // "perp" = depth into the neighbor (0 = first row off the shared edge).
  let along: number;
  const perp = 0;
  if (localEdge === "west" || localEdge === "east") {
    along = j;
  } else {
    along = i;
  }
  if (nb.reverseIndex) {
    along = Nc - 1 - along;
  }
  switch (nb.edge) {
    case "west":
      return { i: perp, j: along };
    case "east":
      return { i: Nc - 1 - perp, j: along };
    case "south":
      return { i: along, j: perp };
    case "north":
      return { i: along, j: Nc - 1 - perp };
  }
}

function checkPanel(panel: number): void {
  if (!Number.isInteger(panel) || panel < 0 || panel > 5) {
    throw new RangeError(`cubed_sphere: invalid panel ${panel} (expected integer 0..5)`);
  }
}

function checkCell(i: number, j: number, Nc: number): void {
  if (!Number.isInteger(i) || i < 0 || i >= Nc) {
    throw new RangeError(`cubed_sphere: invalid i ${i} (expected integer 0..${Nc - 1})`);
  }
  if (!Number.isInteger(j) || j < 0 || j >= Nc) {
    throw new RangeError(`cubed_sphere: invalid j ${j} (expected integer 0..${Nc - 1})`);
  }
}

export function cubed_sphere(opts: CubedSphereOpts): CubedSphereGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("cubed_sphere: options object is required");
  }
  const { Nc } = opts;
  if (Nc === undefined) {
    throw new TypeError("cubed_sphere: missing required option 'Nc'");
  }
  if (!Number.isInteger(Nc) || Nc <= 0) {
    throw new RangeError(`cubed_sphere: 'Nc' must be a positive integer (got ${Nc})`);
  }
  const R = opts.R ?? 6.371e6;
  if (!Number.isFinite(R) || R <= 0) {
    throw new RangeError(`cubed_sphere: 'R' must be a positive finite number (got ${R})`);
  }
  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(`cubed_sphere: 'dtype' must be 'float64' or 'float32' (got ${String(dtype)})`);
  }
  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts < 0) {
    throw new RangeError(`cubed_sphere: 'ghosts' must be a non-negative integer (got ${ghosts})`);
  }

  const dxi = Math.PI / 2 / Nc;
  const xiCenter = (i: number): number => -Math.PI / 4 + (i + 0.5) * dxi;
  const xiEdge = (i: number): number => -Math.PI / 4 + i * dxi;

  const grid: CubedSphereGrid = {
    family: "cubed_sphere",
    topology: "block_structured",
    generator: "gnomonic_c6",
    Nc,
    R,
    dtype,
    ghosts,
    n_panels: 6,
    n_cells: 6 * Nc * Nc,

    cell_centers(panel: number, i: number, j: number) {
      checkPanel(panel);
      checkCell(i, j, Nc);
      return gnomonicToLonLat(xiCenter(i), xiCenter(j), panel);
    },

    neighbors(panel: number, i: number, j: number): NeighborRef[] {
      checkPanel(panel);
      checkCell(i, j, Nc);
      const refs: NeighborRef[] = [];
      const conn = PANEL_CONNECTIVITY[panel];
      // west
      if (i > 0) {
        refs.push({ panel, i: i - 1, j, edge: "west" });
      } else {
        const t = transformGhostIndex(conn.west, i, j, Nc, "west");
        refs.push({ panel: conn.west.panel, i: t.i, j: t.j, edge: "west" });
      }
      // east
      if (i < Nc - 1) {
        refs.push({ panel, i: i + 1, j, edge: "east" });
      } else {
        const t = transformGhostIndex(conn.east, i, j, Nc, "east");
        refs.push({ panel: conn.east.panel, i: t.i, j: t.j, edge: "east" });
      }
      // south
      if (j > 0) {
        refs.push({ panel, i, j: j - 1, edge: "south" });
      } else {
        const t = transformGhostIndex(conn.south, i, j, Nc, "south");
        refs.push({ panel: conn.south.panel, i: t.i, j: t.j, edge: "south" });
      }
      // north
      if (j < Nc - 1) {
        refs.push({ panel, i, j: j + 1, edge: "north" });
      } else {
        const t = transformGhostIndex(conn.north, i, j, Nc, "north");
        refs.push({ panel: conn.north.panel, i: t.i, j: t.j, edge: "north" });
      }
      return refs;
    },

    metric_eval(name: MetricName, panel: number, i: number, j: number): number {
      checkPanel(panel);
      checkCell(i, j, Nc);
      const m = gnomonicMetric(xiCenter(i), xiCenter(j), R);
      switch (name) {
        case "J":
          return m.J;
        case "g_xixi":
          return m.g_xixi;
        case "g_etaeta":
          return m.g_etaeta;
        case "g_xieta":
          return m.g_xieta;
        case "ginv_xixi":
        case "ginv_etaeta":
        case "ginv_xieta": {
          const det = m.g_xixi * m.g_etaeta - m.g_xieta * m.g_xieta;
          if (name === "ginv_xixi") return m.g_etaeta / det;
          if (name === "ginv_etaeta") return m.g_xixi / det;
          return -m.g_xieta / det;
        }
        default:
          throw new RangeError(`cubed_sphere: unknown metric '${String(name)}'`);
      }
    },

    cell_area(panel: number, i: number, j: number): number {
      checkPanel(panel);
      checkCell(i, j, Nc);
      const xw = xiEdge(i);
      const xe = xiEdge(i + 1);
      const ys = xiEdge(j);
      const yn = xiEdge(j + 1);
      const corners: Vec3[] = [
        gnomonicToCart(xw, ys, panel),
        gnomonicToCart(xe, ys, panel),
        gnomonicToCart(xe, yn, panel),
        gnomonicToCart(xw, yn, panel),
      ];
      return cellAreaFromCorners(corners, R);
    },

    toESM(): object {
      return {
        family: "cubed_sphere",
        version: "1.0.0",
        dtype,
        topology: "block_structured",
        generator: "gnomonic_c6",
        options: { Nc, R, ghosts },
        n_panels: 6,
        n_cells: 6 * Nc * Nc,
      };
    },
  };

  return grid;
}
