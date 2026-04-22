/**
 * DUO icosahedral triangular grid family (Heikes et al. 2023) — accessor runtime.
 *
 * Loader-backed per `docs/GRIDS_API.md` §10. `toESM()` emits a small
 * declarative config (family, level, loader ref); geometry is derived at
 * runtime by the accessors. Cross-binding conformance is compared at
 * pinned query points, not serialized bytes (bead dsc-no3, mayor
 * correction 2026-04-20).
 *
 * Subdivision level `r` yields 20·4^r triangular cells, 10·4^r + 2
 * vertices, 30·4^r edges. All vertices lie on a sphere of radius R.
 *
 * Cell indices are 0-based (TS convention). Neighbor slot k (0..2) is
 * across the edge opposite local vertex k, matching the Julia reference
 * binding's 1-based slot semantics shifted down by one.
 */

import type { Dtype, Grid } from "./types.js";

const DUO_FAMILY_VERSION = "1.0.0";
const BUILTIN_PREFIX = "builtin://icosahedral/";

export interface DuoLoader {
  path: string;
  reader?: string;
  check?: string;
}

export interface DuoOpts {
  loader: DuoLoader;
  R?: number;
  dtype?: Dtype;
  ghosts?: number;
}

export type DuoMetricName = "area" | "lon" | "lat" | "x" | "y" | "z";

export interface DuoGrid extends Grid {
  readonly family: "duo";
  readonly topology: "unstructured";
  readonly level: number;
  readonly R: number;
  readonly ghosts: number;
  readonly n_cells: number;
  readonly n_vertices: number;
  readonly n_edges: number;
  readonly loader: Required<DuoLoader>;

  cell_centers(c: number): { lon: number; lat: number };
  neighbors(c: number): [number, number, number];
  metric_eval(name: DuoMetricName, c: number): number;
  cell_area(c: number): number;
}

type Vec3 = [number, number, number];

function parseBuiltinLevel(path: string): number | null {
  if (!path.startsWith(BUILTIN_PREFIX)) return null;
  const tail = path.slice(BUILTIN_PREFIX.length);
  if (tail.length === 0 || !/^-?\d+$/.test(tail)) {
    throw new RangeError(`duo: cannot parse level from loader path '${path}'`);
  }
  const lvl = parseInt(tail, 10);
  if (lvl < 0) {
    throw new RangeError(`duo: subdivision level must be >= 0 (got ${lvl})`);
  }
  return lvl;
}

function resolveLoaderLevel(loader: Required<DuoLoader>): number {
  const lvl = parseBuiltinLevel(loader.path);
  if (lvl !== null) return lvl;
  if (loader.reader === "duo_mesh" || loader.reader === "auto") {
    throw new RangeError(
      "duo: .duo mesh-file reader not yet implemented — pending " +
        "EarthSciSerialization file-format spec. Use " +
        "builtin://icosahedral/<level> in the meantime.",
    );
  }
  throw new RangeError(
    `duo: unrecognized loader path '${loader.path}' with reader='${loader.reader}'`,
  );
}

function coerceLoader(loader: unknown): Required<DuoLoader> {
  if (loader === null || typeof loader !== "object") {
    throw new TypeError("duo: 'loader' option is required and must be an object");
  }
  const l = loader as Record<string, unknown>;
  const path = l.path;
  if (typeof path !== "string" || path.length === 0) {
    throw new TypeError("duo: loader.path is required and must be a non-empty string");
  }
  const reader = l.reader === undefined ? "auto" : l.reader;
  if (typeof reader !== "string") {
    throw new TypeError("duo: loader.reader must be a string");
  }
  const check = l.check === undefined ? "strict" : l.check;
  if (typeof check !== "string") {
    throw new TypeError("duo: loader.check must be a string");
  }
  return { path, reader, check };
}

// --- Base icosahedron -----------------------------------------------------

function icosahedronVertices(): Vec3[] {
  const phi = (1 + Math.sqrt(5)) / 2;
  const raw: Vec3[] = [
    [0, 1, phi],
    [0, -1, phi],
    [0, 1, -phi],
    [0, -1, -phi],
    [1, phi, 0],
    [-1, phi, 0],
    [1, -phi, 0],
    [-1, -phi, 0],
    [phi, 0, 1],
    [phi, 0, -1],
    [-phi, 0, 1],
    [-phi, 0, -1],
  ];
  return raw.map(([x, y, z]): Vec3 => {
    const n = Math.sqrt(x * x + y * y + z * z);
    return [x / n, y / n, z / n];
  });
}

// Mirrors src/grids/duo.jl _icosahedron_faces (1-based → 0-based here).
function icosahedronFaces(): Array<[number, number, number]> {
  const f: Array<[number, number, number]> = [
    [1, 9, 2],
    [1, 2, 11],
    [1, 11, 6],
    [1, 6, 5],
    [1, 5, 9],
    [9, 5, 10],
    [9, 10, 7],
    [9, 7, 2],
    [2, 7, 8],
    [2, 8, 11],
    [11, 8, 12],
    [11, 12, 6],
    [6, 12, 3],
    [6, 3, 5],
    [5, 3, 10],
    [10, 3, 4],
    [10, 4, 7],
    [7, 4, 8],
    [8, 4, 12],
    [12, 4, 3],
  ];
  return f.map(([a, b, c]): [number, number, number] => [a - 1, b - 1, c - 1]);
}

// --- Recursive subdivision ------------------------------------------------

function subdivideIcosahedron(level: number): {
  vertices: Vec3[];
  faces: Array<[number, number, number]>;
} {
  const verts: Vec3[] = icosahedronVertices();
  let faces = icosahedronFaces();

  for (let step = 0; step < level; step++) {
    const cache = new Map<string, number>();
    const midpoint = (a: number, b: number): number => {
      const key = a < b ? `${a}_${b}` : `${b}_${a}`;
      const cached = cache.get(key);
      if (cached !== undefined) return cached;
      const va = verts[a];
      const vb = verts[b];
      const mx = va[0] + vb[0];
      const my = va[1] + vb[1];
      const mz = va[2] + vb[2];
      const n = Math.sqrt(mx * mx + my * my + mz * mz);
      verts.push([mx / n, my / n, mz / n]);
      const idx = verts.length - 1;
      cache.set(key, idx);
      return idx;
    };

    const next: Array<[number, number, number]> = [];
    for (const [a, b, c] of faces) {
      const ab = midpoint(a, b);
      const bc = midpoint(b, c);
      const ca = midpoint(c, a);
      next.push([a, ab, ca]);
      next.push([b, bc, ab]);
      next.push([c, ca, bc]);
      next.push([ab, bc, ca]);
    }
    faces = next;
  }
  return { vertices: verts, faces };
}

// --- Geometry helpers -----------------------------------------------------

function cartToLonLat(x: number, y: number, z: number): { lon: number; lat: number } {
  const clamped = z < -1 ? -1 : z > 1 ? 1 : z;
  return { lon: Math.atan2(y, x), lat: Math.asin(clamped) };
}

// L'Huilier's theorem on the unit sphere; stable on near-degenerate triangles.
function sphericalTriangleArea(a: Vec3, b: Vec3, c: Vec3): number {
  const dot = (u: Vec3, v: Vec3): number => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const acosClamp = (t: number): number => Math.acos(t < -1 ? -1 : t > 1 ? 1 : t);
  const da = acosClamp(dot(b, c));
  const db = acosClamp(dot(c, a));
  const dc = acosClamp(dot(a, b));
  const s = 0.5 * (da + db + dc);
  let t =
    Math.tan(s / 2) *
    Math.tan((s - da) / 2) *
    Math.tan((s - db) / 2) *
    Math.tan((s - dc) / 2);
  if (t < 0) t = 0;
  return 4 * Math.atan(Math.sqrt(t));
}

// --- Connectivity ---------------------------------------------------------

function buildConnectivity(faces: ReadonlyArray<[number, number, number]>): {
  edges: Array<[number, number]>;
  cellNeighbors: Int32Array; // shape (Nc * 3), slot k of cell c at index c*3 + k
} {
  const Nc = faces.length;
  const edgeMap = new Map<string, Array<[number, number]>>();
  for (let c = 0; c < Nc; c++) {
    const [v1, v2, v3] = faces[c];
    // Edge opposite vertex k=0 is (v2,v3); k=1 → (v3,v1); k=2 → (v1,v2).
    const sides: Array<[number, number]> = [
      [v2, v3],
      [v3, v1],
      [v1, v2],
    ];
    for (let k = 0; k < 3; k++) {
      const [a, b] = sides[k];
      const key = a < b ? `${a}_${b}` : `${b}_${a}`;
      const list = edgeMap.get(key);
      if (list === undefined) {
        edgeMap.set(key, [[c, k]]);
      } else {
        list.push([c, k]);
      }
    }
  }

  const edges: Array<[number, number]> = [];
  const cellNeighbors = new Int32Array(Nc * 3).fill(-1);
  for (const [key, list] of edgeMap) {
    const [sa, sb] = key.split("_");
    edges.push([parseInt(sa, 10), parseInt(sb, 10)]);
    if (list.length === 2) {
      const [[c1, k1], [c2, k2]] = list;
      cellNeighbors[c1 * 3 + k1] = c2;
      cellNeighbors[c2 * 3 + k2] = c1;
    } else if (list.length > 2) {
      throw new Error(`duo: non-manifold edge (${key}) shared by ${list.length} cells`);
    }
  }
  return { edges, cellNeighbors };
}

// --- Generator ------------------------------------------------------------

export function duo(opts: DuoOpts): DuoGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("duo: options object is required");
  }
  if (opts.loader === undefined) {
    throw new TypeError("duo: missing required option 'loader'");
  }
  const loader = coerceLoader(opts.loader);
  const level = resolveLoaderLevel(loader);

  const R = opts.R ?? 6.371e6;
  if (!Number.isFinite(R) || R <= 0) {
    throw new RangeError(`duo: 'R' must be a positive finite number (got ${R})`);
  }
  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(`duo: 'dtype' must be 'float64' or 'float32' (got ${String(dtype)})`);
  }
  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts < 0) {
    throw new RangeError(`duo: 'ghosts' must be a non-negative integer (got ${ghosts})`);
  }

  const { vertices, faces } = subdivideIcosahedron(level);
  const Nv = vertices.length;
  const Nc = faces.length;

  // Per-cell centroid (normalized) + area.
  const cellCart = new Float64Array(Nc * 3);
  const lon = new Float64Array(Nc);
  const lat = new Float64Array(Nc);
  const area = new Float64Array(Nc);
  const R2 = R * R;
  for (let c = 0; c < Nc; c++) {
    const [ai, bi, ci] = faces[c];
    const va = vertices[ai];
    const vb = vertices[bi];
    const vc = vertices[ci];
    const mx = va[0] + vb[0] + vc[0];
    const my = va[1] + vb[1] + vc[1];
    const mz = va[2] + vb[2] + vc[2];
    const n = Math.sqrt(mx * mx + my * my + mz * mz);
    const ux = mx / n;
    const uy = my / n;
    const uz = mz / n;
    cellCart[c * 3 + 0] = R * ux;
    cellCart[c * 3 + 1] = R * uy;
    cellCart[c * 3 + 2] = R * uz;
    const ll = cartToLonLat(ux, uy, uz);
    lon[c] = ll.lon;
    lat[c] = ll.lat;
    area[c] = sphericalTriangleArea(va, vb, vc) * R2;
  }

  const { edges, cellNeighbors } = buildConnectivity(faces);

  const checkCell = (c: number): void => {
    if (!Number.isInteger(c) || c < 0 || c >= Nc) {
      throw new RangeError(`duo: invalid cell index ${c} (expected integer 0..${Nc - 1})`);
    }
  };

  const grid: DuoGrid = {
    family: "duo",
    topology: "unstructured",
    dtype,
    level,
    R,
    ghosts,
    n_cells: Nc,
    n_vertices: Nv,
    n_edges: edges.length,
    loader,

    cell_centers(c: number) {
      checkCell(c);
      return { lon: lon[c], lat: lat[c] };
    },

    neighbors(c: number): [number, number, number] {
      checkCell(c);
      return [
        cellNeighbors[c * 3 + 0],
        cellNeighbors[c * 3 + 1],
        cellNeighbors[c * 3 + 2],
      ];
    },

    metric_eval(name: DuoMetricName, c: number): number {
      checkCell(c);
      switch (name) {
        case "area":
          return area[c];
        case "lon":
          return lon[c];
        case "lat":
          return lat[c];
        case "x":
          return cellCart[c * 3 + 0];
        case "y":
          return cellCart[c * 3 + 1];
        case "z":
          return cellCart[c * 3 + 2];
        default:
          throw new RangeError(`duo: metric_eval: unknown metric '${String(name)}'`);
      }
    },

    cell_area(c: number): number {
      checkCell(c);
      return area[c];
    },

    toESM(): object {
      return {
        family: "duo",
        version: DUO_FAMILY_VERSION,
        topology: "unstructured",
        dtype,
        ghosts,
        n_cells: Nc,
        n_vertices: Nv,
        n_edges: edges.length,
        options: {
          R,
          level,
          loader: {
            path: loader.path,
            reader: loader.reader,
            check: loader.check,
          },
        },
        provenance: {
          binding: "typescript",
          family: "duo",
          version: DUO_FAMILY_VERSION,
          level,
          reader: loader.reader,
          path: loader.path,
          check: loader.check,
          dtype,
        },
        schema_version: DUO_FAMILY_VERSION,
      };
    },
  };

  return grid;
}
