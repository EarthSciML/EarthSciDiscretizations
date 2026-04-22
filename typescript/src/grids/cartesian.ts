/**
 * Cartesian grid family (1D / 2D / 3D, uniform + non-uniform) — accessor runtime.
 *
 * Per `docs/GRIDS_API.md` §2.6 and the 2026-04-20 mayor correction in bead
 * dsc-suu: `toESM()` emits a small declarative config (family + dims +
 * extents + optional non-uniform edges), NOT a serialized geometry blob.
 * Geometry is derived on demand via the accessors (`cell_centers`,
 * `cell_widths`, `cell_volume`, `neighbors`, `metric_eval`).
 *
 * Mirrors the Julia reference binding in `src/grids/cartesian.jl`. Indices
 * are 0-based in TS (Julia is 1-based); semantics and accessor outputs are
 * otherwise identical.
 */

import type { Dtype, Grid } from "./types.js";

const CARTESIAN_API_VERSION = "1.0.0";
const CARTESIAN_BINDING_VERSION = "1.0.0";

export type AxisExtent = readonly [number, number];

export interface CartesianOpts {
  nx?: number;
  ny?: number;
  nz?: number;
  extent?: ReadonlyArray<AxisExtent>;
  edges?: ReadonlyArray<ReadonlyArray<number>>;
  dtype?: Dtype;
  ghosts?: number;
}

export type CartesianMetricName =
  | "volume"
  | "jacobian"
  | "g"
  | "dx"
  | "dy"
  | "dz"
  | "face_area_x"
  | "face_area_y"
  | "face_area_z";

export interface CartesianNeighborRef {
  axis: number; // 0-based axis index
  side: -1 | 1;
  index: number[];
}

export interface CartesianGrid extends Grid {
  readonly family: "cartesian";
  readonly topology: "rectilinear";
  readonly ndim: 1 | 2 | 3;
  readonly n: number[]; // per-axis cell counts
  readonly extent: AxisExtent[]; // per-axis (lo, hi)
  readonly edges: number[][]; // per-axis edge arrays (length n[d]+1)
  readonly centers: number[][]; // per-axis cell centers
  readonly widths: number[][]; // per-axis cell widths
  readonly uniform: boolean[];
  readonly ghosts: number;
  readonly n_cells: number;

  cell_centers(...idx: number[]): number[];
  cell_widths(...idx: number[]): number[];
  cell_volume(...idx: number[]): number;
  neighbors(...idx: number[]): CartesianNeighborRef[];
  metric_eval(
    name: CartesianMetricName,
    ...idx: number[]
  ): number | number[][];
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function isPositiveInteger(v: unknown): v is number {
  return typeof v === "number" && Number.isInteger(v) && v >= 1;
}

function validateExtentAxis(
  lo: unknown,
  hi: unknown,
  axis: number,
): [number, number] {
  if (typeof lo !== "number" || !Number.isFinite(lo)) {
    throw new RangeError(
      `cartesian: extent[${axis}] lo must be a finite number (got ${String(lo)})`,
    );
  }
  if (typeof hi !== "number" || !Number.isFinite(hi)) {
    throw new RangeError(
      `cartesian: extent[${axis}] hi must be a finite number (got ${String(hi)})`,
    );
  }
  if (!(hi > lo)) {
    throw new RangeError(
      `cartesian: extent[${axis}] must satisfy hi > lo (got lo=${lo}, hi=${hi})`,
    );
  }
  return [lo, hi];
}

function normalizeExtent(
  extent: ReadonlyArray<AxisExtent>,
  ndim: number,
): [number, number][] {
  if (!Array.isArray(extent)) {
    throw new TypeError(
      `cartesian: 'extent' must be an array of (lo, hi) pairs (got ${typeof extent})`,
    );
  }
  if (extent.length !== ndim) {
    throw new RangeError(
      `cartesian: 'extent' has ${extent.length} axes, expected ${ndim}`,
    );
  }
  const out: [number, number][] = [];
  for (let d = 0; d < ndim; d++) {
    const e = extent[d];
    if (!Array.isArray(e) || e.length !== 2) {
      throw new TypeError(
        `cartesian: extent[${d}] must be a 2-element (lo, hi); got ${JSON.stringify(e)}`,
      );
    }
    out.push(validateExtentAxis(e[0], e[1], d));
  }
  return out;
}

function uniformAxis(
  n: number,
  lo: number,
  hi: number,
): { edges: number[]; centers: number[]; widths: number[] } {
  const dx = (hi - lo) / n;
  const edges: number[] = new Array(n + 1);
  const centers: number[] = new Array(n);
  const widths: number[] = new Array(n);
  for (let i = 0; i <= n; i++) {
    edges[i] = lo + i * dx;
  }
  for (let i = 0; i < n; i++) {
    centers[i] = (edges[i] + edges[i + 1]) / 2;
    widths[i] = edges[i + 1] - edges[i];
  }
  return { edges, centers, widths };
}

function nonUniformAxis(
  edgesIn: ReadonlyArray<number>,
  axis: number,
): { n: number; edges: number[]; centers: number[]; widths: number[] } {
  if (!Array.isArray(edgesIn)) {
    throw new TypeError(
      `cartesian: edges[${axis}] must be an array of numbers`,
    );
  }
  if (edgesIn.length < 2) {
    throw new RangeError(
      `cartesian: edges[${axis}] must have ≥ 2 entries; got ${edgesIn.length}`,
    );
  }
  const edges: number[] = new Array(edgesIn.length);
  for (let i = 0; i < edgesIn.length; i++) {
    const e = edgesIn[i];
    if (typeof e !== "number" || !Number.isFinite(e)) {
      throw new RangeError(
        `cartesian: edges[${axis}][${i}] must be a finite number (got ${String(e)})`,
      );
    }
    edges[i] = e;
  }
  for (let i = 0; i + 1 < edges.length; i++) {
    if (!(edges[i + 1] > edges[i])) {
      throw new RangeError(
        `cartesian: edges[${axis}] must be strictly increasing (violation at index ${i})`,
      );
    }
  }
  const n = edges.length - 1;
  const centers: number[] = new Array(n);
  const widths: number[] = new Array(n);
  for (let i = 0; i < n; i++) {
    centers[i] = (edges[i] + edges[i + 1]) / 2;
    widths[i] = edges[i + 1] - edges[i];
  }
  return { n, edges, centers, widths };
}

function isUniform(widths: ReadonlyArray<number>): boolean {
  if (widths.length <= 1) return true;
  const w0 = widths[0];
  const tol = Math.max(Number.EPSILON, Number.EPSILON * Math.abs(w0)) * 8;
  for (let i = 1; i < widths.length; i++) {
    if (Math.abs(widths[i] - w0) > tol) return false;
  }
  return true;
}

function checkIdxArity(idx: readonly number[], ndim: number): void {
  if (idx.length !== ndim) {
    throw new RangeError(
      `cartesian: expected ${ndim} cell index argument(s), got ${idx.length}`,
    );
  }
}

function checkIdxRange(idx: readonly number[], n: readonly number[]): void {
  for (let d = 0; d < idx.length; d++) {
    const v = idx[d];
    if (!Number.isInteger(v) || v < 0 || v >= n[d]) {
      throw new RangeError(
        `cartesian: index[${d}]=${v} out of range [0, ${n[d] - 1}]`,
      );
    }
  }
}

// ---------------------------------------------------------------------------
// Generator
// ---------------------------------------------------------------------------

export function cartesian(opts: CartesianOpts): CartesianGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("cartesian: options object is required");
  }

  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(
      `cartesian: 'dtype' must be 'float64' or 'float32' (got ${String(dtype)})`,
    );
  }

  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts < 0) {
    throw new RangeError(
      `cartesian: 'ghosts' must be a non-negative integer (got ${ghosts})`,
    );
  }

  const hasN = opts.nx !== undefined || opts.ny !== undefined || opts.nz !== undefined;
  const hasExtent = opts.extent !== undefined;

  let ndim: 1 | 2 | 3;
  let n: number[];
  let extentN: [number, number][];
  let axisEdges: number[][];
  let axisCenters: number[][];
  let axisWidths: number[][];
  let uniform: boolean[];

  if (opts.edges !== undefined) {
    if (hasN || hasExtent) {
      throw new TypeError(
        "cartesian: pass either 'edges' (non-uniform) OR 'nx/ny/nz'+'extent' (uniform), not both",
      );
    }
    if (!Array.isArray(opts.edges) || opts.edges.length === 0) {
      throw new TypeError(
        "cartesian: 'edges' must be a non-empty array of edge arrays",
      );
    }
    const k = opts.edges.length;
    if (k < 1 || k > 3) {
      throw new RangeError(
        `cartesian: ndim must be in 1..3 (got ${k} edge arrays)`,
      );
    }
    ndim = k as 1 | 2 | 3;
    n = new Array(ndim);
    extentN = new Array(ndim);
    axisEdges = new Array(ndim);
    axisCenters = new Array(ndim);
    axisWidths = new Array(ndim);
    uniform = new Array(ndim);
    for (let d = 0; d < ndim; d++) {
      const axis = nonUniformAxis(opts.edges[d], d);
      n[d] = axis.n;
      extentN[d] = [axis.edges[0], axis.edges[axis.edges.length - 1]];
      axisEdges[d] = axis.edges;
      axisCenters[d] = axis.centers;
      axisWidths[d] = axis.widths;
      uniform[d] = isUniform(axis.widths);
    }
  } else {
    if (opts.nx === undefined) {
      throw new TypeError("cartesian: missing required option 'nx'");
    }
    if (!isPositiveInteger(opts.nx)) {
      throw new RangeError(
        `cartesian: 'nx' must be a positive integer (got ${opts.nx})`,
      );
    }
    // Determine ndim from which of nx/ny/nz were provided.
    if (opts.nz !== undefined && opts.ny === undefined) {
      throw new TypeError(
        "cartesian: missing required option 'ny' for 3D grid",
      );
    }
    if (opts.ny === undefined) {
      ndim = 1;
      n = [opts.nx];
    } else if (opts.nz === undefined) {
      if (!isPositiveInteger(opts.ny)) {
        throw new RangeError(
          `cartesian: 'ny' must be a positive integer (got ${opts.ny})`,
        );
      }
      ndim = 2;
      n = [opts.nx, opts.ny];
    } else {
      if (!isPositiveInteger(opts.ny)) {
        throw new RangeError(
          `cartesian: 'ny' must be a positive integer (got ${opts.ny})`,
        );
      }
      if (!isPositiveInteger(opts.nz)) {
        throw new RangeError(
          `cartesian: 'nz' must be a positive integer (got ${opts.nz})`,
        );
      }
      ndim = 3;
      n = [opts.nx, opts.ny, opts.nz];
    }
    if (!hasExtent) {
      throw new TypeError("cartesian: missing required option 'extent'");
    }
    extentN = normalizeExtent(opts.extent as ReadonlyArray<AxisExtent>, ndim);

    axisEdges = new Array(ndim);
    axisCenters = new Array(ndim);
    axisWidths = new Array(ndim);
    uniform = new Array(ndim);
    for (let d = 0; d < ndim; d++) {
      const a = uniformAxis(n[d], extentN[d][0], extentN[d][1]);
      axisEdges[d] = a.edges;
      axisCenters[d] = a.centers;
      axisWidths[d] = a.widths;
      uniform[d] = true;
    }
  }

  const nCells = n.reduce((p, v) => p * v, 1);

  const grid: CartesianGrid = {
    family: "cartesian",
    topology: "rectilinear",
    dtype,
    ndim,
    n,
    extent: extentN,
    edges: axisEdges,
    centers: axisCenters,
    widths: axisWidths,
    uniform,
    ghosts,
    n_cells: nCells,

    cell_centers(...idx: number[]): number[] {
      checkIdxArity(idx, ndim);
      checkIdxRange(idx, n);
      const out: number[] = new Array(ndim);
      for (let d = 0; d < ndim; d++) out[d] = axisCenters[d][idx[d]];
      return out;
    },

    cell_widths(...idx: number[]): number[] {
      checkIdxArity(idx, ndim);
      checkIdxRange(idx, n);
      const out: number[] = new Array(ndim);
      for (let d = 0; d < ndim; d++) out[d] = axisWidths[d][idx[d]];
      return out;
    },

    cell_volume(...idx: number[]): number {
      checkIdxArity(idx, ndim);
      checkIdxRange(idx, n);
      let v = 1;
      for (let d = 0; d < ndim; d++) v *= axisWidths[d][idx[d]];
      return v;
    },

    neighbors(...idx: number[]): CartesianNeighborRef[] {
      checkIdxArity(idx, ndim);
      checkIdxRange(idx, n);
      const out: CartesianNeighborRef[] = [];
      for (let d = 0; d < ndim; d++) {
        if (idx[d] > 0) {
          const nbr = idx.slice();
          nbr[d] = idx[d] - 1;
          out.push({ axis: d, side: -1, index: nbr });
        }
        if (idx[d] < n[d] - 1) {
          const nbr = idx.slice();
          nbr[d] = idx[d] + 1;
          out.push({ axis: d, side: 1, index: nbr });
        }
      }
      return out;
    },

    metric_eval(
      name: CartesianMetricName,
      ...idx: number[]
    ): number | number[][] {
      checkIdxArity(idx, ndim);
      checkIdxRange(idx, n);
      switch (name) {
        case "volume": {
          let v = 1;
          for (let d = 0; d < ndim; d++) v *= axisWidths[d][idx[d]];
          return v;
        }
        case "jacobian":
          return 1;
        case "g": {
          const m: number[][] = new Array(ndim);
          for (let i = 0; i < ndim; i++) {
            const row: number[] = new Array(ndim);
            for (let j = 0; j < ndim; j++) row[j] = i === j ? 1 : 0;
            m[i] = row;
          }
          return m;
        }
        case "dx":
          return axisWidths[0][idx[0]];
        case "dy":
          if (ndim < 2) {
            throw new RangeError(
              `cartesian: metric 'dy' not defined for ${ndim}D grid`,
            );
          }
          return axisWidths[1][idx[1]];
        case "dz":
          if (ndim < 3) {
            throw new RangeError(
              `cartesian: metric 'dz' not defined for ${ndim}D grid`,
            );
          }
          return axisWidths[2][idx[2]];
        case "face_area_x":
          return faceArea(0, idx, ndim, axisWidths);
        case "face_area_y":
          if (ndim < 2) {
            throw new RangeError(
              `cartesian: metric 'face_area_y' not defined for ${ndim}D grid`,
            );
          }
          return faceArea(1, idx, ndim, axisWidths);
        case "face_area_z":
          if (ndim < 3) {
            throw new RangeError(
              `cartesian: metric 'face_area_z' not defined for ${ndim}D grid`,
            );
          }
          return faceArea(2, idx, ndim, axisWidths);
        default:
          throw new RangeError(
            `cartesian: unknown metric '${String(name)}'`,
          );
      }
    },

    toESM(): object {
      const base: Record<string, unknown> = {
        family: "cartesian",
        version: CARTESIAN_API_VERSION,
        topology: "rectilinear",
        dtype,
        ndim,
        ghosts,
        n_cells: nCells,
        n: n.slice(),
        extent: extentN.map((e) => [e[0], e[1]]),
        uniform: uniform.slice(),
        provenance: {
          binding: "typescript",
          binding_version: CARTESIAN_BINDING_VERSION,
          api_version: CARTESIAN_API_VERSION,
          platform: "js",
          runtime: "node",
          math_lib: "libm-js",
          source_sha: "",
        },
      };
      if (uniform.some((u) => !u)) {
        base.edges = axisEdges.map((eArr, d) =>
          uniform[d] ? [] : eArr.slice(),
        );
      }
      return base;
    },
  };

  return grid;
}

function faceArea(
  axis: number,
  idx: readonly number[],
  ndim: number,
  widths: readonly number[][],
): number {
  // Product of widths on all axes EXCEPT `axis`. In 1D this is empty → 1.
  let a = 1;
  for (let d = 0; d < ndim; d++) {
    if (d === axis) continue;
    a *= widths[d][idx[d]];
  }
  return a;
}
