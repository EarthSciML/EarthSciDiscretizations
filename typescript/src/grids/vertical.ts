/**
 * Vertical grid family (1D column) — accessor runtime.
 *
 * Per `docs/GRIDS_API.md` §2.6 and the 2026-04-20 mayor correction in bead
 * dsc-0k6: `toESM()` emits a small declarative config (family + coordinate
 * kind + interface levels + optional hybrid coefficients), NOT a serialized
 * geometry blob. Centers and widths are derived on demand from the
 * interface `levels` via pure arithmetic.
 *
 * Rhymes with the Python reference in
 * `python/src/earthsci_toolkit/grids/vertical.py`. Indices are 0-based.
 *
 * Supported coordinate kinds:
 *
 * | Kind                   | Value domain                  | Required options        |
 * |------------------------|-------------------------------|-------------------------|
 * | `"sigma"`              | [0, 1]; 1 = surface, 0 = top  | `nz` or `levels`        |
 * | `"eta"`                | hybrid sigma-pressure (NCAR)  | `ak`, `bk`              |
 * | `"z"`                  | geometric altitude (m)        | `levels` (increasing)   |
 * | `"theta"`              | potential temperature (K)     | `levels` (increasing)   |
 * | `"hybrid_sigma_theta"` | blended sigma→theta           | `nz` or `levels`        |
 * | `"z_star"`             | generalized height            | `levels` (increasing)   |
 */

import type { Dtype, Grid } from "./types.js";

const VERTICAL_API_VERSION = "1.0.0";
const VERTICAL_BINDING_VERSION = "0.1.0";

export type VerticalCoordinate =
  | "sigma"
  | "eta"
  | "z"
  | "theta"
  | "hybrid_sigma_theta"
  | "z_star";

const COORDINATES: ReadonlyArray<VerticalCoordinate> = [
  "sigma",
  "eta",
  "z",
  "theta",
  "hybrid_sigma_theta",
  "z_star",
];

export type VerticalMetricName =
  | "dz"
  | "z"
  | "sigma"
  | "pressure"
  | "ak"
  | "bk";

const METRIC_NAMES: ReadonlyArray<VerticalMetricName> = [
  "dz",
  "z",
  "sigma",
  "pressure",
  "ak",
  "bk",
];

export interface VerticalOpts {
  coordinate: VerticalCoordinate;
  nz?: number;
  levels?: ReadonlyArray<number>;
  ak?: ReadonlyArray<number>;
  bk?: ReadonlyArray<number>;
  p0?: number;
  transition?: number;
  dtype?: Dtype;
  ghosts?: number;
}

export interface VerticalNeighbors {
  down?: number;
  up?: number;
}

export interface VerticalGrid extends Grid {
  readonly family: "vertical";
  readonly topology: "column";
  readonly coordinate: VerticalCoordinate;
  readonly nz: number;
  readonly levels: number[];
  readonly centers: number[];
  readonly widths: number[];
  readonly ak: number[];
  readonly bk: number[];
  readonly p0: number;
  readonly transition: number | null;
  readonly ghosts: number;
  readonly ndim: 1;
  readonly n_cells: number;
  readonly n_vertices: number;
  readonly n_edges: number;

  cell_centers(k?: number): number | number[];
  cell_widths(k?: number): number | number[];
  neighbors(k: number): VerticalNeighbors;
  metric_eval(name: VerticalMetricName, k: number): number;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function isNonNegativeInteger(v: unknown): v is number {
  return typeof v === "number" && Number.isInteger(v) && v >= 0;
}

function isPositiveFiniteNumber(v: unknown): v is number {
  return typeof v === "number" && Number.isFinite(v) && v > 0;
}

function coerceFiniteArray(
  values: ReadonlyArray<number>,
  label: string,
): number[] {
  const out: number[] = new Array(values.length);
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (typeof v !== "number" || !Number.isFinite(v)) {
      throw new RangeError(
        `vertical: ${label}[${i}] must be a finite number (got ${String(v)})`,
      );
    }
    out[i] = v;
  }
  return out;
}

function coerceLevels(
  levels: ReadonlyArray<number>,
  label: string,
  mustDecrease: boolean,
  domain: [number, number] | null,
): number[] {
  if (!Array.isArray(levels)) {
    throw new TypeError(
      `vertical: \`${label}\` must be an array of numbers (got ${typeof levels})`,
    );
  }
  if (levels.length < 2) {
    throw new RangeError(
      `vertical: \`${label}\` must have >= 2 entries (got ${levels.length})`,
    );
  }
  const arr = coerceFiniteArray(levels, label);
  if (domain !== null) {
    const [lo, hi] = domain;
    for (let i = 0; i < arr.length; i++) {
      if (arr[i] < lo || arr[i] > hi) {
        throw new RangeError(
          `vertical: \`${label}\` entries must lie in [${lo}, ${hi}] (index ${i} = ${arr[i]})`,
        );
      }
    }
  }
  for (let i = 0; i + 1 < arr.length; i++) {
    const d = arr[i + 1] - arr[i];
    if (mustDecrease ? !(d < 0) : !(d > 0)) {
      throw new RangeError(
        `vertical: \`${label}\` must be strictly ${
          mustDecrease ? "decreasing" : "increasing"
        } (violation at index ${i})`,
      );
    }
  }
  return arr;
}

function coerceHybrid(
  coeffs: ReadonlyArray<number> | undefined,
  expectedLen: number,
  label: string,
): number[] {
  if (coeffs === undefined || coeffs === null) {
    throw new TypeError(`vertical: \`${label}\` is required`);
  }
  if (!Array.isArray(coeffs)) {
    throw new TypeError(
      `vertical: \`${label}\` must be an array of numbers (got ${typeof coeffs})`,
    );
  }
  if (coeffs.length !== expectedLen) {
    throw new RangeError(
      `vertical: \`${label}\` must have length nz+1 = ${expectedLen}; got ${coeffs.length}`,
    );
  }
  return coerceFiniteArray(coeffs, label);
}

function uniformSigmaLevels(nz: number): number[] {
  const out: number[] = new Array(nz + 1);
  for (let k = 0; k <= nz; k++) out[k] = 1 - k / nz;
  return out;
}

function centersAndWidths(levels: ReadonlyArray<number>): {
  centers: number[];
  widths: number[];
} {
  const n = levels.length - 1;
  const centers: number[] = new Array(n);
  const widths: number[] = new Array(n);
  for (let k = 0; k < n; k++) {
    centers[k] = 0.5 * (levels[k] + levels[k + 1]);
    widths[k] = Math.abs(levels[k + 1] - levels[k]);
  }
  return { centers, widths };
}

function checkLayer(k: unknown, nz: number): void {
  if (typeof k !== "number" || !Number.isInteger(k)) {
    throw new TypeError(
      `vertical: k must be an integer (got ${String(k)})`,
    );
  }
  if (k < 0 || k >= nz) {
    throw new RangeError(
      `vertical: k=${k} out of range [0, ${nz})`,
    );
  }
}

// ---------------------------------------------------------------------------
// Generator
// ---------------------------------------------------------------------------

export function vertical(opts: VerticalOpts): VerticalGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("vertical: options object is required");
  }

  const coordinate = opts.coordinate;
  if (coordinate === undefined || coordinate === null) {
    throw new TypeError("vertical: required option `coordinate` is missing");
  }
  if (typeof coordinate !== "string") {
    throw new TypeError(
      `vertical: \`coordinate\` must be a string (got ${typeof coordinate})`,
    );
  }
  if (!COORDINATES.includes(coordinate as VerticalCoordinate)) {
    throw new RangeError(
      `vertical: unknown coordinate '${coordinate}'; expected one of ${COORDINATES.join(
        ", ",
      )}`,
    );
  }

  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(
      `vertical: \`dtype\` must be 'float64' or 'float32' (got ${String(dtype)})`,
    );
  }

  const ghosts = opts.ghosts ?? 0;
  if (!isNonNegativeInteger(ghosts)) {
    throw new RangeError(
      `vertical: \`ghosts\` must be a non-negative integer (got ${String(ghosts)})`,
    );
  }

  const p0 = opts.p0 ?? 1.0e5;
  if (!isPositiveFiniteNumber(p0)) {
    throw new RangeError(
      `vertical: \`p0\` must be a positive finite number (got ${String(p0)})`,
    );
  }

  let nz = opts.nz;
  if (nz !== undefined) {
    if (!Number.isInteger(nz) || (nz as number) < 1) {
      throw new RangeError(
        `vertical: \`nz\` must be an integer >= 1 (got ${String(nz)})`,
      );
    }
  }

  let levels: number[];
  let ak: number[] = [];
  let bk: number[] = [];
  let transition: number | null = null;

  if (coordinate === "sigma") {
    if (opts.levels !== undefined) {
      levels = coerceLevels(opts.levels, "levels", true, [0, 1]);
      if (nz !== undefined && levels.length !== nz + 1) {
        throw new RangeError(
          `vertical: sigma: nz=${nz} inconsistent with levels length ${levels.length}`,
        );
      }
    } else {
      if (nz === undefined) {
        throw new TypeError("vertical: 'sigma' requires `nz` or `levels`");
      }
      levels = uniformSigmaLevels(nz);
    }
  } else if (
    coordinate === "z" ||
    coordinate === "theta" ||
    coordinate === "z_star"
  ) {
    if (opts.levels === undefined) {
      throw new TypeError(
        `vertical: '${coordinate}' requires explicit \`levels\``,
      );
    }
    levels = coerceLevels(opts.levels, "levels", false, null);
    if (nz !== undefined && levels.length !== nz + 1) {
      throw new RangeError(
        `vertical: ${coordinate}: nz=${nz} inconsistent with levels length ${levels.length}`,
      );
    }
  } else if (coordinate === "eta") {
    if (opts.ak === undefined) {
      throw new TypeError("vertical: 'eta' requires `ak` (length nz+1)");
    }
    if (opts.bk === undefined) {
      throw new TypeError("vertical: 'eta' requires `bk` (length nz+1)");
    }
    if (!Array.isArray(opts.ak) || !Array.isArray(opts.bk)) {
      throw new TypeError(
        "vertical: 'eta' `ak` and `bk` must be arrays of numbers",
      );
    }
    if (opts.ak.length !== opts.bk.length) {
      throw new RangeError(
        `vertical: 'eta' ak/bk must have equal length; got ${opts.ak.length} vs ${opts.bk.length}`,
      );
    }
    const nzEff = nz !== undefined ? nz : opts.ak.length - 1;
    if (nzEff < 1) {
      throw new RangeError(`vertical: nz must be >= 1 (got ${nzEff})`);
    }
    ak = coerceHybrid(opts.ak, nzEff + 1, "ak");
    bk = coerceHybrid(opts.bk, nzEff + 1, "bk");
    const sigma: number[] = new Array(ak.length);
    for (let i = 0; i < ak.length; i++) sigma[i] = ak[i] / p0 + bk[i];
    for (let i = 0; i + 1 < sigma.length; i++) {
      if (!(sigma[i + 1] - sigma[i] < 0)) {
        throw new RangeError(
          "vertical: 'eta' synthesized sigma (ak/p0 + bk) must be strictly decreasing",
        );
      }
    }
    levels = sigma;
    nz = nzEff;
  } else {
    // hybrid_sigma_theta
    if (opts.levels === undefined && nz === undefined) {
      throw new TypeError(
        "vertical: 'hybrid_sigma_theta' requires `nz` or `levels`",
      );
    }
    if (opts.levels !== undefined) {
      levels = coerceLevels(opts.levels, "levels", true, [0, 1]);
      if (nz !== undefined && levels.length !== nz + 1) {
        throw new RangeError(
          `vertical: hybrid_sigma_theta: nz=${nz} inconsistent with levels length ${levels.length}`,
        );
      }
    } else {
      levels = uniformSigmaLevels(nz as number);
    }
    if (opts.transition !== undefined) {
      if (
        typeof opts.transition !== "number" ||
        !Number.isFinite(opts.transition) ||
        !(opts.transition > 0 && opts.transition < 1)
      ) {
        throw new RangeError(
          `vertical: hybrid_sigma_theta \`transition\` must be in (0, 1); got ${String(opts.transition)}`,
        );
      }
      transition = opts.transition;
    }
    ak =
      opts.ak !== undefined ? coerceHybrid(opts.ak, levels.length, "ak") : [];
    bk =
      opts.bk !== undefined ? coerceHybrid(opts.bk, levels.length, "bk") : [];
  }

  const nzResolved = levels.length - 1;
  const { centers, widths } = centersAndWidths(levels);

  const grid: VerticalGrid = {
    family: "vertical",
    topology: "column",
    dtype,
    coordinate,
    nz: nzResolved,
    levels,
    centers,
    widths,
    ak,
    bk,
    p0,
    transition,
    ghosts,
    ndim: 1,
    n_cells: nzResolved,
    n_vertices: nzResolved + 1,
    n_edges: nzResolved,

    cell_centers(k?: number): number | number[] {
      if (k === undefined) return centers.slice();
      checkLayer(k, nzResolved);
      return centers[k];
    },

    cell_widths(k?: number): number | number[] {
      if (k === undefined) return widths.slice();
      checkLayer(k, nzResolved);
      return widths[k];
    },

    neighbors(k: number): VerticalNeighbors {
      checkLayer(k, nzResolved);
      const out: VerticalNeighbors = {};
      if (k > 0) out.down = k - 1;
      if (k < nzResolved - 1) out.up = k + 1;
      return out;
    },

    metric_eval(name: VerticalMetricName, k: number): number {
      if (!METRIC_NAMES.includes(name)) {
        throw new RangeError(`vertical: unknown metric '${String(name)}'`);
      }
      checkLayer(k, nzResolved);
      switch (name) {
        case "dz":
          return widths[k];
        case "z":
          return centers[k];
        case "sigma":
          if (
            coordinate === "sigma" ||
            coordinate === "hybrid_sigma_theta" ||
            coordinate === "eta"
          ) {
            return centers[k];
          }
          throw new RangeError(
            `vertical: 'sigma' undefined for coordinate '${coordinate}'`,
          );
        case "pressure":
          if (ak.length === 0 || bk.length === 0) {
            throw new RangeError(
              `vertical: 'pressure' requires hybrid ak/bk (coordinate '${coordinate}' has none)`,
            );
          }
          return 0.5 * (ak[k] + bk[k] * p0 + ak[k + 1] + bk[k + 1] * p0);
        case "ak":
          if (ak.length === 0) {
            throw new RangeError(
              "vertical: 'ak' unavailable (no hybrid coefficients)",
            );
          }
          return 0.5 * (ak[k] + ak[k + 1]);
        case "bk":
          if (bk.length === 0) {
            throw new RangeError(
              "vertical: 'bk' unavailable (no hybrid coefficients)",
            );
          }
          return 0.5 * (bk[k] + bk[k + 1]);
      }
    },

    toESM(): object {
      const options: Record<string, unknown> = {
        coordinate,
        nz: nzResolved,
        levels: levels.slice(),
      };
      if (ak.length > 0) options.ak = ak.slice();
      if (bk.length > 0) options.bk = bk.slice();
      if (
        coordinate === "eta" ||
        coordinate === "hybrid_sigma_theta" ||
        ak.length > 0 ||
        bk.length > 0
      ) {
        options.p0 = p0;
      }
      if (transition !== null) {
        options.transition = transition;
      }
      return {
        family: "vertical",
        topology: "column",
        dtype,
        ndim: 1,
        ghosts,
        n_cells: nzResolved,
        n_vertices: nzResolved + 1,
        n_edges: nzResolved,
        options,
        provenance: {
          binding: "typescript",
          binding_version: VERTICAL_BINDING_VERSION,
          family: "vertical",
          version: VERTICAL_API_VERSION,
          coordinate,
          dtype,
        },
        schema_version: VERTICAL_API_VERSION,
      };
    },
  };

  return grid;
}
