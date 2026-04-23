/**
 * Lat-lon grid family (regular + reduced-Gaussian) — accessor runtime.
 *
 * Per `docs/GRIDS_API.md` §2.6, §3.4, §7 and the 2026-04-20 mayor correction
 * in bead dsc-suu: `toESM()` emits a small declarative config (family +
 * variant + dimensions + per-row schedule for reduced-Gaussian + optional
 * latitude edges), NOT a serialized geometry blob. Geometry is derived on
 * demand via the accessors (`cell_center`, `cell_centers_bulk`, `neighbors`,
 * `cell_area`, `metric_eval`).
 *
 * Mirrors the Rust binding in `rust/src/grids/lat_lon.rs`. Indices are
 * 0-based; semantics, defaults, and accessor outputs are identical.
 *
 * Polar-singularity handling: only `pole_policy="none"` is implemented in
 * this phase (boundary cells at the poles have no N/S neighbor). The
 * "average" and "fold" policies are declared in the schema but rejected at
 * build time.
 */

import type { Dtype, Grid } from "./types.js";

const LAT_LON_API_VERSION = "1.0.0";
const LAT_LON_BINDING_VERSION = "1.0.0";

const PI = Math.PI;
const HALF_PI = PI / 2;
const DEFAULT_R = 6.371e6;
const DEFAULT_LON_START = -PI;

export type LatLonVariant = "regular" | "reduced_gaussian";
export type PolePolicy = "none" | "average" | "fold";

/** Cardinal direction for face neighbors. */
export type Direction = "W" | "E" | "S" | "N";
export const DIRECTIONS: readonly Direction[] = ["W", "E", "S", "N"];

export interface LatLonOpts {
  variant?: LatLonVariant;
  nlon?: number;
  nlat?: number;
  nlon_per_row?: ReadonlyArray<number>;
  lat_edges?: ReadonlyArray<number>;
  lat_centers?: ReadonlyArray<number>;
  lon_start?: number;
  R?: number;
  ghosts?: number;
  dtype?: Dtype;
  pole_policy?: PolePolicy;
}

export type LatLonMetricName =
  | "J"
  | "g_lonlon"
  | "g_latlat"
  | "g_lonlat"
  | "ginv_lonlon"
  | "ginv_latlat"
  | "ginv_lonlat"
  | "area";

/**
 * Neighbor cell. `null` is returned at the poles under `pole_policy="none"`
 * to signal a genuine boundary rather than a wrap.
 */
export interface LatLonNeighborCell {
  j: number;
  i: number;
}

export type LatLonNeighborSet = {
  W: LatLonNeighborCell;
  E: LatLonNeighborCell;
  S: LatLonNeighborCell | null;
  N: LatLonNeighborCell | null;
};

export interface LatLonGrid extends Grid {
  readonly family: "lat_lon";
  readonly topology: "rectilinear";
  readonly variant: LatLonVariant;
  readonly nlat: number;
  readonly nlon_per_row: number[];
  readonly nlon_uniform: number | null;
  readonly R: number;
  readonly ghosts: number;
  readonly pole_policy: PolePolicy;
  readonly lon_start: number;
  readonly lat_edges: number[]; // length nlat+1
  readonly lat_centers: number[]; // length nlat
  readonly n_cells: number;

  nlon(j: number): number;
  lon_edges(j: number): number[];
  lon_centers(j: number): number[];
  cell_center(j: number, i: number): { lon: number; lat: number };
  cell_centers_bulk(): { lon: number[]; lat: number[] };
  row_offset(j: number): number;
  neighbors(j: number, i: number): LatLonNeighborSet;
  neighbor(j: number, i: number, dir: Direction): LatLonNeighborCell | null;
  cell_area(j: number, i: number): number;
  area_bulk(): number[];
  metric_eval(name: LatLonMetricName, j: number, i: number): number;
}

// ---------------------------------------------------------------------------
// Validation helpers
// ---------------------------------------------------------------------------

function isPositiveInteger(v: unknown): v is number {
  return typeof v === "number" && Number.isInteger(v) && v >= 1;
}

function isFiniteNumber(v: unknown): v is number {
  return typeof v === "number" && Number.isFinite(v);
}

function copyFiniteIntArray(
  src: ReadonlyArray<number>,
  optName: string,
): number[] {
  const out: number[] = new Array(src.length);
  for (let k = 0; k < src.length; k++) {
    const v = src[k];
    if (!Number.isInteger(v) || v < 1) {
      throw new RangeError(
        `lat_lon: ${optName}[${k}]=${String(v)} must be a positive integer`,
      );
    }
    out[k] = v;
  }
  return out;
}

function copyFiniteFloatArray(
  src: ReadonlyArray<number>,
  optName: string,
): number[] {
  const out: number[] = new Array(src.length);
  for (let k = 0; k < src.length; k++) {
    const v = src[k];
    if (!isFiniteNumber(v)) {
      throw new RangeError(
        `lat_lon: ${optName}[${k}]=${String(v)} must be a finite number`,
      );
    }
    out[k] = v;
  }
  return out;
}

function validateLatEdges(edges: ReadonlyArray<number>, nlat: number): void {
  if (edges.length !== nlat + 1) {
    throw new RangeError(
      `lat_lon: 'lat_edges' length ${edges.length} does not match nlat+1=${nlat + 1}`,
    );
  }
  for (let k = 0; k < nlat; k++) {
    if (!(edges[k + 1] > edges[k])) {
      throw new RangeError(
        `lat_lon: 'lat_edges' must be strictly increasing (violation at index ${k})`,
      );
    }
  }
  // Allow a tiny epsilon overshoot to tolerate floating-point representation
  // of ±π/2 from Gaussian-quadrature schemes.
  if (edges[0] < -HALF_PI - 1e-12 || edges[nlat] > HALF_PI + 1e-12) {
    throw new RangeError(
      `lat_lon: 'lat_edges' must lie in [-pi/2, pi/2]`,
    );
  }
}

function validateLatCenters(
  centers: ReadonlyArray<number>,
  edges: ReadonlyArray<number>,
  nlat: number,
): void {
  if (centers.length !== nlat) {
    throw new RangeError(
      `lat_lon: 'lat_centers' length ${centers.length} does not match nlat=${nlat}`,
    );
  }
  for (let k = 0; k < nlat; k++) {
    const c = centers[k];
    if (!(c >= edges[k] && c <= edges[k + 1])) {
      throw new RangeError(
        `lat_lon: 'lat_centers'[${k}]=${c} outside enclosing edges [${edges[k]}, ${edges[k + 1]}]`,
      );
    }
  }
}

/**
 * Map a column index `i` in a row of width `from` to the nearest-center
 * column in a row of width `to`. Mirrors the Rust `map_i` helper for
 * cross-binding parity in reduced-Gaussian neighbor lookup.
 */
function mapI(i: number, from: number, to: number): number {
  if (from === to) return i;
  const frac = (i + 0.5) / from;
  const k = Math.floor(frac * to);
  return Math.min(k, to - 1);
}

// ---------------------------------------------------------------------------
// Generator
// ---------------------------------------------------------------------------

export function lat_lon(opts: LatLonOpts): LatLonGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("lat_lon: options object is required");
  }

  const variant: LatLonVariant = opts.variant ?? "regular";
  if (variant !== "regular" && variant !== "reduced_gaussian") {
    throw new RangeError(
      `lat_lon: 'variant' must be 'regular' or 'reduced_gaussian' (got ${String(variant)})`,
    );
  }

  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(
      `lat_lon: 'dtype' must be 'float64' or 'float32' (got ${String(dtype)})`,
    );
  }

  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts < 0) {
    throw new RangeError(
      `lat_lon: 'ghosts' must be a non-negative integer (got ${ghosts})`,
    );
  }

  const polePolicy: PolePolicy = opts.pole_policy ?? "none";
  if (polePolicy !== "none" && polePolicy !== "average" && polePolicy !== "fold") {
    throw new RangeError(
      `lat_lon: 'pole_policy' must be 'none', 'average', or 'fold' (got ${String(polePolicy)})`,
    );
  }
  if (polePolicy !== "none") {
    throw new RangeError(
      `lat_lon: pole_policy='${polePolicy}' is declared but not implemented in this phase`,
    );
  }

  const R = opts.R ?? DEFAULT_R;
  if (!isFiniteNumber(R) || !(R > 0)) {
    throw new RangeError(
      `lat_lon: 'R' must be a positive finite number (got ${String(R)})`,
    );
  }

  const lonStart = opts.lon_start ?? DEFAULT_LON_START;
  if (!isFiniteNumber(lonStart)) {
    throw new RangeError(
      `lat_lon: 'lon_start' must be a finite number (got ${String(lonStart)})`,
    );
  }

  let nlat: number;
  let nlonPerRow: number[];

  if (variant === "regular") {
    if (opts.nlon_per_row !== undefined) {
      throw new RangeError(
        "lat_lon: 'nlon_per_row' is not allowed for variant='regular'",
      );
    }
    if (opts.nlon === undefined) {
      throw new TypeError("lat_lon: missing required option 'nlon'");
    }
    if (opts.nlat === undefined) {
      throw new TypeError("lat_lon: missing required option 'nlat'");
    }
    if (!isPositiveInteger(opts.nlon)) {
      throw new RangeError(
        `lat_lon: 'nlon' must be a positive integer (got ${opts.nlon})`,
      );
    }
    if (!isPositiveInteger(opts.nlat)) {
      throw new RangeError(
        `lat_lon: 'nlat' must be a positive integer (got ${opts.nlat})`,
      );
    }
    nlat = opts.nlat;
    nlonPerRow = new Array(nlat).fill(opts.nlon);
  } else {
    if (opts.nlon !== undefined) {
      throw new RangeError(
        "lat_lon: 'nlon' is not allowed for variant='reduced_gaussian'; use 'nlon_per_row'",
      );
    }
    if (opts.nlon_per_row === undefined) {
      throw new TypeError(
        "lat_lon: missing required option 'nlon_per_row' for variant='reduced_gaussian'",
      );
    }
    if (!Array.isArray(opts.nlon_per_row) || opts.nlon_per_row.length === 0) {
      throw new RangeError(
        "lat_lon: 'nlon_per_row' must be a non-empty array of positive integers",
      );
    }
    nlonPerRow = copyFiniteIntArray(opts.nlon_per_row, "nlon_per_row");
    nlat = opts.nlat ?? nlonPerRow.length;
    if (!isPositiveInteger(nlat)) {
      throw new RangeError(
        `lat_lon: 'nlat' must be a positive integer (got ${nlat})`,
      );
    }
    if (nlonPerRow.length !== nlat) {
      throw new RangeError(
        `lat_lon: 'nlon_per_row' length ${nlonPerRow.length} does not match nlat=${nlat}`,
      );
    }
  }

  let latEdges: number[];
  if (opts.lat_edges !== undefined) {
    if (!Array.isArray(opts.lat_edges)) {
      throw new TypeError("lat_lon: 'lat_edges' must be an array of numbers");
    }
    latEdges = copyFiniteFloatArray(opts.lat_edges, "lat_edges");
    validateLatEdges(latEdges, nlat);
  } else {
    const dlat = PI / nlat;
    latEdges = new Array(nlat + 1);
    for (let k = 0; k <= nlat; k++) latEdges[k] = -HALF_PI + k * dlat;
  }

  let latCenters: number[];
  if (opts.lat_centers !== undefined) {
    if (!Array.isArray(opts.lat_centers)) {
      throw new TypeError("lat_lon: 'lat_centers' must be an array of numbers");
    }
    latCenters = copyFiniteFloatArray(opts.lat_centers, "lat_centers");
    validateLatCenters(latCenters, latEdges, nlat);
  } else {
    latCenters = new Array(nlat);
    for (let k = 0; k < nlat; k++) {
      latCenters[k] = 0.5 * (latEdges[k] + latEdges[k + 1]);
    }
  }

  // Pre-compute row offsets for the ragged bulk layout.
  const rowOffsets: number[] = new Array(nlat + 1);
  rowOffsets[0] = 0;
  for (let j = 0; j < nlat; j++) {
    rowOffsets[j + 1] = rowOffsets[j] + nlonPerRow[j];
  }
  const nCells = rowOffsets[nlat];

  const nlonUniform = variant === "regular" ? nlonPerRow[0] : null;

  function checkCell(j: number, i: number): void {
    if (!Number.isInteger(j) || j < 0 || j >= nlat) {
      throw new RangeError(
        `lat_lon: j=${j} out of range [0, ${nlat - 1}]`,
      );
    }
    const n = nlonPerRow[j];
    if (!Number.isInteger(i) || i < 0 || i >= n) {
      throw new RangeError(
        `lat_lon: i=${i} out of range [0, ${n - 1}] for row ${j}`,
      );
    }
  }

  function checkRow(j: number): number {
    if (!Number.isInteger(j) || j < 0 || j >= nlat) {
      throw new RangeError(
        `lat_lon: j=${j} out of range [0, ${nlat - 1}]`,
      );
    }
    return nlonPerRow[j];
  }

  const grid: LatLonGrid = {
    family: "lat_lon",
    topology: "rectilinear",
    variant,
    dtype,
    nlat,
    nlon_per_row: nlonPerRow,
    nlon_uniform: nlonUniform,
    R,
    ghosts,
    pole_policy: polePolicy,
    lon_start: lonStart,
    lat_edges: latEdges,
    lat_centers: latCenters,
    n_cells: nCells,

    nlon(j: number): number {
      return checkRow(j);
    },

    lon_edges(j: number): number[] {
      const n = checkRow(j);
      const dlon = (2 * PI) / n;
      const out: number[] = new Array(n + 1);
      for (let k = 0; k <= n; k++) out[k] = lonStart + k * dlon;
      return out;
    },

    lon_centers(j: number): number[] {
      const n = checkRow(j);
      const dlon = (2 * PI) / n;
      const out: number[] = new Array(n);
      for (let k = 0; k < n; k++) out[k] = lonStart + (k + 0.5) * dlon;
      return out;
    },

    cell_center(j: number, i: number): { lon: number; lat: number } {
      checkCell(j, i);
      const n = nlonPerRow[j];
      const dlon = (2 * PI) / n;
      return {
        lon: lonStart + (i + 0.5) * dlon,
        lat: latCenters[j],
      };
    },

    cell_centers_bulk(): { lon: number[]; lat: number[] } {
      const lon: number[] = new Array(nCells);
      const lat: number[] = new Array(nCells);
      let p = 0;
      for (let j = 0; j < nlat; j++) {
        const n = nlonPerRow[j];
        const dlon = (2 * PI) / n;
        const latC = latCenters[j];
        for (let i = 0; i < n; i++) {
          lon[p] = lonStart + (i + 0.5) * dlon;
          lat[p] = latC;
          p++;
        }
      }
      return { lon, lat };
    },

    row_offset(j: number): number {
      if (!Number.isInteger(j) || j < 0 || j > nlat) {
        throw new RangeError(
          `lat_lon: row offset query j=${j} out of range [0, ${nlat}]`,
        );
      }
      return rowOffsets[j];
    },

    neighbors(j: number, i: number): LatLonNeighborSet {
      checkCell(j, i);
      const n_i = nlonPerRow[j];
      const w: LatLonNeighborCell = { j, i: i === 0 ? n_i - 1 : i - 1 };
      const e: LatLonNeighborCell = { j, i: i + 1 === n_i ? 0 : i + 1 };
      const s: LatLonNeighborCell | null =
        j === 0
          ? null
          : { j: j - 1, i: mapI(i, n_i, nlonPerRow[j - 1]) };
      const n: LatLonNeighborCell | null =
        j + 1 === nlat
          ? null
          : { j: j + 1, i: mapI(i, n_i, nlonPerRow[j + 1]) };
      return { W: w, E: e, S: s, N: n };
    },

    neighbor(
      j: number,
      i: number,
      dir: Direction,
    ): LatLonNeighborCell | null {
      const ns = this.neighbors(j, i);
      switch (dir) {
        case "W":
          return ns.W;
        case "E":
          return ns.E;
        case "S":
          return ns.S;
        case "N":
          return ns.N;
        default:
          throw new RangeError(
            `lat_lon: unknown direction '${String(dir)}'`,
          );
      }
    },

    cell_area(j: number, i: number): number {
      checkCell(j, i);
      const n = nlonPerRow[j];
      const dlon = (2 * PI) / n;
      const latS = latEdges[j];
      const latN = latEdges[j + 1];
      return R * R * dlon * (Math.sin(latN) - Math.sin(latS));
    },

    area_bulk(): number[] {
      const out: number[] = new Array(nCells);
      let p = 0;
      for (let j = 0; j < nlat; j++) {
        const n = nlonPerRow[j];
        const dlon = (2 * PI) / n;
        const latS = latEdges[j];
        const latN = latEdges[j + 1];
        const a = R * R * dlon * (Math.sin(latN) - Math.sin(latS));
        for (let i = 0; i < n; i++) {
          out[p++] = a;
        }
      }
      return out;
    },

    metric_eval(name: LatLonMetricName, j: number, i: number): number {
      checkCell(j, i);
      if (name === "area") return this.cell_area(j, i);
      const lat = latCenters[j];
      const cosLat = Math.cos(lat);
      const r2 = R * R;
      const gLL = r2 * cosLat * cosLat;
      const gPP = r2;
      switch (name) {
        case "J":
          return r2 * Math.abs(cosLat);
        case "g_lonlon":
          return gLL;
        case "g_latlat":
          return gPP;
        case "g_lonlat":
          return 0;
        case "ginv_lonlon":
          return gLL > 0 ? 1 / gLL : Infinity;
        case "ginv_latlat":
          return 1 / gPP;
        case "ginv_lonlat":
          return 0;
        default:
          throw new RangeError(
            `lat_lon: unknown metric '${String(name)}'`,
          );
      }
    },

    toESM(): object {
      const params: Record<string, unknown> =
        variant === "regular"
          ? {
              nlon: nlonPerRow[0],
              nlat,
              R,
              ghosts,
              pole_policy: polePolicy,
              lon_start: lonStart,
            }
          : {
              nlat,
              nlon_per_row: nlonPerRow.slice(),
              lat_edges: latEdges.slice(),
              R,
              ghosts,
              pole_policy: polePolicy,
              lon_start: lonStart,
            };
      return {
        family: "lat_lon",
        version: LAT_LON_API_VERSION,
        dtype,
        topology: "rectilinear",
        variant,
        generator:
          variant === "regular"
            ? "lat_lon_regular"
            : "lat_lon_reduced_gaussian",
        params,
        provenance: {
          binding: "typescript",
          binding_version: LAT_LON_BINDING_VERSION,
          api_version: LAT_LON_API_VERSION,
          platform: "js",
          runtime: "node",
          math_lib: "libm-js",
          generator:
            variant === "regular"
              ? "lat_lon_regular"
              : "lat_lon_reduced_gaussian",
          source_sha: "",
        },
      };
    },
  };

  return grid;
}
