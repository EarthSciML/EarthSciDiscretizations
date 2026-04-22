/**
 * MPAS unstructured Voronoi grid (loader-backed, accessor runtime).
 *
 * Implements the `mpas` family per `docs/GRIDS_API.md` §2.6 / §10 and the
 * mayor's 2026-04-20 correction on bead dsc-3nw: `toESM()` emits a small
 * declarative config (family, dimensions, loader ref), not a serialized
 * geometry blob. Geometry is derived on demand by the accessors
 * (`cell_centers`, `neighbors`, `cell_area`, `edge_length`, `metric_eval`).
 * Cross-binding conformance is compared at pinned query points.
 *
 * No new runtime dependencies are introduced. Path-based loading requires
 * the caller to supply a `readerFn` that translates a filesystem path into
 * an `MPASMeshData`; this keeps NetCDF I/O out of the default dependency
 * graph while satisfying the §10 loader contract. In-memory construction
 * via `mpasMeshData({...})` is the primary path for tests and host-built
 * meshes.
 *
 * Index convention: cells, edges, and adjacency slots are 0-based (TS
 * convention). The sentinel for "no neighbor" in adjacency arrays is
 * `-1`. This differs from the Julia reference (1-based with `0`-sentinel
 * per NetCDF convention) but does not affect cross-binding conformance,
 * which compares accessor outputs at pinned query points rather than
 * serialized adjacency arrays.
 */

import type { Dtype, Grid } from "./types.js";

const MPAS_FAMILY_VERSION = "1.0.0";
const MPAS_SOURCE_SHA = "dsc-3nw";
const MPAS_READER_VERSION = "0.1.0";

export type MpasReader = "auto" | "nc4" | "mpas_mesh";
export type MpasCheck = "strict" | "lenient";

export interface MpasLoader {
  path: string;
  reader?: MpasReader;
  check?: MpasCheck;
}

export interface MpasMeshInput {
  lon_cell: ArrayLike<number>;
  lat_cell: ArrayLike<number>;
  area_cell: ArrayLike<number>;
  n_edges_on_cell: ArrayLike<number>;
  /** Row-major (n_cells × max_edges). Entry -1 denotes "no neighbor". */
  cells_on_cell: ArrayLike<number>;
  /** Row-major (n_cells × max_edges). Entry -1 denotes "no edge". */
  edges_on_cell: ArrayLike<number>;
  lon_edge: ArrayLike<number>;
  lat_edge: ArrayLike<number>;
  /** Row-major (n_edges × 2). Entry -1 denotes "external boundary". */
  cells_on_edge: ArrayLike<number>;
  dc_edge: ArrayLike<number>;
  dv_edge: ArrayLike<number>;
  max_edges: number;
  x_cell?: ArrayLike<number>;
  y_cell?: ArrayLike<number>;
  z_cell?: ArrayLike<number>;
  n_vertices?: number;
  R?: number;
}

export interface MPASMeshData {
  readonly n_cells: number;
  readonly n_edges: number;
  readonly n_vertices: number;
  readonly max_edges: number;
  readonly lon_cell: Float64Array;
  readonly lat_cell: Float64Array;
  readonly x_cell: Float64Array;
  readonly y_cell: Float64Array;
  readonly z_cell: Float64Array;
  readonly area_cell: Float64Array;
  readonly n_edges_on_cell: Int32Array;
  readonly cells_on_cell: Int32Array;
  readonly edges_on_cell: Int32Array;
  readonly lon_edge: Float64Array;
  readonly lat_edge: Float64Array;
  readonly cells_on_edge: Int32Array;
  readonly dc_edge: Float64Array;
  readonly dv_edge: Float64Array;
}

export type MpasMetricName =
  | "lon"
  | "lat"
  | "area"
  | "x"
  | "y"
  | "z"
  | "n_edges_on_cell"
  | "lon_edge"
  | "lat_edge"
  | "dc_edge"
  | "dv_edge";

export interface MpasOpts {
  loader?: MpasLoader;
  mesh?: MPASMeshData;
  R?: number;
  dtype?: Dtype;
  ghosts?: number;
  /** Path-based loading requires a caller-supplied reader — §10 contract. */
  readerFn?: (path: string) => MPASMeshData;
}

export interface MPASGrid extends Grid {
  readonly family: "mpas";
  readonly topology: "unstructured";
  readonly R: number;
  readonly ghosts: number;
  readonly n_cells: number;
  readonly n_edges: number;
  readonly n_vertices: number;
  readonly max_edges: number;
  readonly loader: Required<MpasLoader> | null;
  readonly mesh: MPASMeshData;

  cell_centers(c: number): { lon: number; lat: number };
  cell_center_cart(c: number): { x: number; y: number; z: number };
  neighbors(c: number): number[];
  cell_area(c: number): number;
  edge_length(e: number): number;
  cell_distance(e: number): number;
  metric_eval(name: MpasMetricName, i: number): number;
  total_area(): number;
}

function toFloat64(a: ArrayLike<number>, label: string, n: number): Float64Array {
  if (a.length !== n) {
    throw new RangeError(`mpas: ${label} length ${a.length} != ${n}`);
  }
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) out[i] = Number(a[i]);
  return out;
}

function toInt32(a: ArrayLike<number>, label: string, n: number): Int32Array {
  if (a.length !== n) {
    throw new RangeError(`mpas: ${label} length ${a.length} != ${n}`);
  }
  const out = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const v = a[i];
    if (!Number.isInteger(v)) {
      throw new TypeError(`mpas: ${label}[${i}] must be an integer (got ${v})`);
    }
    out[i] = v;
  }
  return out;
}

/**
 * Validated constructor for the in-memory MPAS mesh. Derives cartesian
 * cell coordinates from (lon, lat, R) when not supplied.
 */
export function mpasMeshData(input: MpasMeshInput): MPASMeshData {
  const n_cells = input.lon_cell.length;
  const n_edges = input.lon_edge.length;
  const max_edges = input.max_edges;
  if (!Number.isInteger(max_edges) || max_edges <= 0) {
    throw new RangeError(`mpas: max_edges must be a positive integer (got ${max_edges})`);
  }
  const R = input.R ?? 6.371e6;

  const lon_cell = toFloat64(input.lon_cell, "lon_cell", n_cells);
  const lat_cell = toFloat64(input.lat_cell, "lat_cell", n_cells);
  const area_cell = toFloat64(input.area_cell, "area_cell", n_cells);
  const n_edges_on_cell = toInt32(input.n_edges_on_cell, "n_edges_on_cell", n_cells);
  const cells_on_cell = toInt32(
    input.cells_on_cell,
    "cells_on_cell",
    n_cells * max_edges,
  );
  const edges_on_cell = toInt32(
    input.edges_on_cell,
    "edges_on_cell",
    n_cells * max_edges,
  );
  const lon_edge = toFloat64(input.lon_edge, "lon_edge", n_edges);
  const lat_edge = toFloat64(input.lat_edge, "lat_edge", n_edges);
  const cells_on_edge = toInt32(input.cells_on_edge, "cells_on_edge", n_edges * 2);
  const dc_edge = toFloat64(input.dc_edge, "dc_edge", n_edges);
  const dv_edge = toFloat64(input.dv_edge, "dv_edge", n_edges);

  const x_cell =
    input.x_cell !== undefined
      ? toFloat64(input.x_cell, "x_cell", n_cells)
      : new Float64Array(n_cells);
  const y_cell =
    input.y_cell !== undefined
      ? toFloat64(input.y_cell, "y_cell", n_cells)
      : new Float64Array(n_cells);
  const z_cell =
    input.z_cell !== undefined
      ? toFloat64(input.z_cell, "z_cell", n_cells)
      : new Float64Array(n_cells);

  if (input.x_cell === undefined) {
    for (let c = 0; c < n_cells; c++) {
      const clat = Math.cos(lat_cell[c]);
      x_cell[c] = R * clat * Math.cos(lon_cell[c]);
    }
  }
  if (input.y_cell === undefined) {
    for (let c = 0; c < n_cells; c++) {
      const clat = Math.cos(lat_cell[c]);
      y_cell[c] = R * clat * Math.sin(lon_cell[c]);
    }
  }
  if (input.z_cell === undefined) {
    for (let c = 0; c < n_cells; c++) {
      z_cell[c] = R * Math.sin(lat_cell[c]);
    }
  }

  return {
    n_cells,
    n_edges,
    n_vertices: input.n_vertices ?? 0,
    max_edges,
    lon_cell,
    lat_cell,
    x_cell,
    y_cell,
    z_cell,
    area_cell,
    n_edges_on_cell,
    cells_on_cell,
    edges_on_cell,
    lon_edge,
    lat_edge,
    cells_on_edge,
    dc_edge,
    dv_edge,
  };
}

function coerceLoader(loader: MpasLoader): Required<MpasLoader> {
  if (typeof loader.path !== "string" || loader.path.length === 0) {
    throw new TypeError("mpas: loader.path is required and must be a non-empty string");
  }
  const reader = loader.reader ?? "auto";
  if (reader !== "auto" && reader !== "nc4" && reader !== "mpas_mesh") {
    throw new RangeError(
      `mpas: loader.reader must be 'auto' | 'nc4' | 'mpas_mesh' (got '${String(reader)}')`,
    );
  }
  const check = loader.check ?? "strict";
  if (check !== "strict" && check !== "lenient") {
    throw new RangeError(
      `mpas: loader.check must be 'strict' | 'lenient' (got '${String(check)}')`,
    );
  }
  return { path: loader.path, reader, check };
}

/**
 * Strict-mode validation: bounds-check adjacency arrays and enforce
 * neighbor-link reciprocity. Lenient mode skips the reciprocity check.
 */
export function checkMesh(m: MPASMeshData, strict: boolean): void {
  const { n_cells, n_edges, max_edges } = m;
  for (let c = 0; c < n_cells; c++) {
    const k = m.n_edges_on_cell[c];
    if (!Number.isInteger(k) || k < 0 || k > max_edges) {
      throw new Error(
        `mpas: n_edges_on_cell[${c}]=${k} out of [0, ${max_edges}]`,
      );
    }
    for (let j = 0; j < k; j++) {
      const nb = m.cells_on_cell[c * max_edges + j];
      if (nb < -1 || nb >= n_cells) {
        throw new Error(
          `mpas: cells_on_cell[${c},${j}]=${nb} out of [-1, ${n_cells - 1}]`,
        );
      }
      const e = m.edges_on_cell[c * max_edges + j];
      if (e < -1 || e >= n_edges) {
        throw new Error(
          `mpas: edges_on_cell[${c},${j}]=${e} out of [-1, ${n_edges - 1}]`,
        );
      }
    }
  }
  for (let e = 0; e < n_edges; e++) {
    for (let s = 0; s < 2; s++) {
      const c = m.cells_on_edge[e * 2 + s];
      if (c < -1 || c >= n_cells) {
        throw new Error(
          `mpas: cells_on_edge[${e},${s}]=${c} out of [-1, ${n_cells - 1}]`,
        );
      }
    }
  }
  if (!strict) return;
  for (let c = 0; c < n_cells; c++) {
    const k = m.n_edges_on_cell[c];
    for (let j = 0; j < k; j++) {
      const nb = m.cells_on_cell[c * max_edges + j];
      if (nb < 0) continue;
      const kb = m.n_edges_on_cell[nb];
      let found = false;
      for (let jj = 0; jj < kb; jj++) {
        if (m.cells_on_cell[nb * max_edges + jj] === c) {
          found = true;
          break;
        }
      }
      if (!found) {
        throw new Error(
          `mpas: neighbor symmetry broken: cell ${c} -> ${nb} but not reverse`,
        );
      }
    }
  }
}

export function mpas(opts: MpasOpts): MPASGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("mpas: options object is required");
  }
  const R = opts.R ?? 6.371e6;
  if (!Number.isFinite(R) || R <= 0) {
    throw new RangeError(`mpas: 'R' must be a positive finite number (got ${R})`);
  }
  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(
      `mpas: 'dtype' must be 'float64' or 'float32' (got '${String(dtype)}')`,
    );
  }
  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts !== 0) {
    throw new RangeError(
      `mpas: 'ghosts' must be 0 for loader-backed grids (got ${ghosts})`,
    );
  }

  let mesh: MPASMeshData;
  let loaderRecord: Required<MpasLoader> | null = null;

  if (opts.mesh !== undefined) {
    mesh = opts.mesh;
    if (opts.loader !== undefined) {
      loaderRecord = coerceLoader(opts.loader);
    }
  } else if (opts.loader !== undefined) {
    loaderRecord = coerceLoader(opts.loader);
    if (opts.readerFn === undefined) {
      throw new TypeError(
        "mpas: path-based loading requires a 'readerFn(path) => MPASMeshData'. " +
          "NetCDF I/O is not bundled with @earthsci/grids per GRIDS_API.md §10; " +
          "pass a readerFn from your consumer package.",
      );
    }
    const produced = opts.readerFn(loaderRecord.path);
    if (produced === null || typeof produced !== "object") {
      throw new TypeError("mpas: readerFn must return an MPASMeshData object");
    }
    mesh = produced;
  } else {
    throw new TypeError(
      "mpas: provide either 'mesh' (MPASMeshData) or 'loader' (MpasLoader)",
    );
  }

  const strict = loaderRecord === null ? false : loaderRecord.check === "strict";
  checkMesh(mesh, strict);

  const checkCell = (c: number): void => {
    if (!Number.isInteger(c) || c < 0 || c >= mesh.n_cells) {
      throw new RangeError(
        `mpas: invalid cell index ${c} (expected integer 0..${mesh.n_cells - 1})`,
      );
    }
  };
  const checkEdge = (e: number): void => {
    if (!Number.isInteger(e) || e < 0 || e >= mesh.n_edges) {
      throw new RangeError(
        `mpas: invalid edge index ${e} (expected integer 0..${mesh.n_edges - 1})`,
      );
    }
  };

  const grid: MPASGrid = {
    family: "mpas",
    topology: "unstructured",
    dtype,
    R,
    ghosts,
    n_cells: mesh.n_cells,
    n_edges: mesh.n_edges,
    n_vertices: mesh.n_vertices,
    max_edges: mesh.max_edges,
    loader: loaderRecord,
    mesh,

    cell_centers(c: number) {
      checkCell(c);
      return { lon: mesh.lon_cell[c], lat: mesh.lat_cell[c] };
    },

    cell_center_cart(c: number) {
      checkCell(c);
      return { x: mesh.x_cell[c], y: mesh.y_cell[c], z: mesh.z_cell[c] };
    },

    neighbors(c: number): number[] {
      checkCell(c);
      const k = mesh.n_edges_on_cell[c];
      const out: number[] = [];
      for (let j = 0; j < k; j++) {
        const nb = mesh.cells_on_cell[c * mesh.max_edges + j];
        if (nb >= 0) out.push(nb);
      }
      return out;
    },

    cell_area(c: number): number {
      checkCell(c);
      return mesh.area_cell[c];
    },

    edge_length(e: number): number {
      checkEdge(e);
      return mesh.dv_edge[e];
    },

    cell_distance(e: number): number {
      checkEdge(e);
      return mesh.dc_edge[e];
    },

    metric_eval(name: MpasMetricName, i: number): number {
      switch (name) {
        case "lon":
          checkCell(i);
          return mesh.lon_cell[i];
        case "lat":
          checkCell(i);
          return mesh.lat_cell[i];
        case "area":
          checkCell(i);
          return mesh.area_cell[i];
        case "x":
          checkCell(i);
          return mesh.x_cell[i];
        case "y":
          checkCell(i);
          return mesh.y_cell[i];
        case "z":
          checkCell(i);
          return mesh.z_cell[i];
        case "n_edges_on_cell":
          checkCell(i);
          return mesh.n_edges_on_cell[i];
        case "lon_edge":
          checkEdge(i);
          return mesh.lon_edge[i];
        case "lat_edge":
          checkEdge(i);
          return mesh.lat_edge[i];
        case "dc_edge":
          checkEdge(i);
          return mesh.dc_edge[i];
        case "dv_edge":
          checkEdge(i);
          return mesh.dv_edge[i];
        default:
          throw new RangeError(
            `mpas: metric_eval: unknown metric '${String(name)}'`,
          );
      }
    },

    total_area(): number {
      let s = 0;
      for (let c = 0; c < mesh.n_cells; c++) s += mesh.area_cell[c];
      return s;
    },

    toESM(): object {
      const loader = loaderRecord === null
        ? null
        : {
            path: loaderRecord.path,
            reader: loaderRecord.reader,
            check: loaderRecord.check,
          };
      return {
        family: "mpas",
        version: MPAS_FAMILY_VERSION,
        dtype,
        topology: "unstructured",
        ghosts,
        n_cells: mesh.n_cells,
        n_edges: mesh.n_edges,
        n_vertices: mesh.n_vertices,
        max_edges: mesh.max_edges,
        options: {
          R,
          loader,
        },
        provenance: {
          binding: "typescript",
          binding_version: "@earthsci/grids 0.1.0",
          family: "mpas",
          version: MPAS_FAMILY_VERSION,
          source_sha: MPAS_SOURCE_SHA,
          reader_version: MPAS_READER_VERSION,
          loader,
          dtype,
        },
        schema_version: MPAS_FAMILY_VERSION,
      };
    },
  };

  return grid;
}
