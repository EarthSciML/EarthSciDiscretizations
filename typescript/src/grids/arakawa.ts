/**
 * Arakawa staggering accessor runtime.
 *
 * Conforms to the cross-binding contract in `docs/GRIDS_API.md` §2.6, §3.4, §7.
 *
 * Arakawa staggering is a *transform* over an underlying base grid. Given a
 * base grid (Cartesian or lat-lon) plus a stagger label (A/B/C/D/E), the
 * runtime provides on-demand accessors for cell centers, u-face, v-face, and
 * corner locations per the staggering convention.
 *
 * Per the 2026-04-20 scope correction, `toESM()` emits a small declarative
 * config (family + base-grid ref + stagger + dimensions + extents), NOT a
 * serialized geometry blob. Geometry is derived from that config by the
 * accessors on demand.
 *
 * Stagger conventions (2-D, horizontal):
 *
 *   - A: h, u, v colocated at cell centers.
 *   - B: h at cell centers; u, v colocated at corners.
 *   - C: h at cell centers; u at u-faces (east-west); v at v-faces (north-south).
 *   - D: h at cell centers; u at v-faces; v at u-faces (swapped from C).
 *   - E: rotated B-grid; topologically equivalent to B with a 45° rotation flag
 *     in the lowered config.
 *
 * Mirrors the Julia reference binding (`src/grids/arakawa.jl`) and the Rust
 * binding (`rust/src/grids/arakawa.rs`). Indices are 0-based in TS, matching
 * `cartesian.ts` and `cubed_sphere.ts`.
 */

import type { Dtype, Grid } from "./types.js";

const ARAKAWA_API_VERSION = "1.0.0";
const ARAKAWA_BINDING_VERSION = "1.0.0";

// ---------------------------------------------------------------------------
// Stagger + Location + Variable vocabularies (snake_case on the wire).
// ---------------------------------------------------------------------------

export type ArakawaStagger = "A" | "B" | "C" | "D" | "E";

export type ArakawaLocation = "cell_center" | "u_edge" | "v_edge" | "corner";

export type ArakawaVariable = "h" | "u" | "v";

export type ArakawaMetricName = "dx" | "dy" | "area";

const STAGGERS: readonly ArakawaStagger[] = ["A", "B", "C", "D", "E"];
const LOCATIONS: readonly ArakawaLocation[] = [
  "cell_center",
  "u_edge",
  "v_edge",
  "corner",
];
const VARIABLES: readonly ArakawaVariable[] = ["h", "u", "v"];
const METRICS: readonly ArakawaMetricName[] = ["dx", "dy", "area"];

/** Per-stagger `(h_loc, u_loc, v_loc)` variable-location triple. */
export function variableLocations(
  s: ArakawaStagger,
): [ArakawaLocation, ArakawaLocation, ArakawaLocation] {
  switch (s) {
    case "A":
      return ["cell_center", "cell_center", "cell_center"];
    case "B":
      return ["cell_center", "corner", "corner"];
    case "C":
      return ["cell_center", "u_edge", "v_edge"];
    case "D":
      return ["cell_center", "v_edge", "u_edge"];
    case "E":
      // E: topologically like B; 45° rotation carried in `rotated` flag.
      return ["cell_center", "corner", "corner"];
  }
}

/** Shape `(ni, nj)` of `loc` on a base grid of `(nx, ny)` interior cells. */
export function locationShape(
  loc: ArakawaLocation,
  nx: number,
  ny: number,
): [number, number] {
  switch (loc) {
    case "cell_center":
      return [nx, ny];
    case "u_edge":
      return [nx + 1, ny];
    case "v_edge":
      return [nx, ny + 1];
    case "corner":
      return [nx + 1, ny + 1];
  }
}

// ---------------------------------------------------------------------------
// Base grid abstraction.
//
// Phase 1 (cartesian) and Phase 2/3 (lat-lon) grid types exist but have their
// own options shapes. Arakawa is parameterised over an `ArakawaBaseGrid`; a
// minimal `CartesianBase` is provided here so the accessors are exercisable
// today. Future base-grid families implement the same interface and plug in
// without touching this file.
// ---------------------------------------------------------------------------

/** Primitives every arakawa base grid must expose. */
export interface ArakawaBaseGrid {
  /** Interior cells along x. */
  readonly nx: number;
  /** Interior cells along y. */
  readonly ny: number;
  /** Uniform cell spacing in x. */
  dx(): number;
  /** Uniform cell spacing in y. */
  dy(): number;
  /** `(x, y)` of cell center `(i, j)` (0-based). */
  cellCenter(i: number, j: number): [number, number];
  /** `(x, y)` of the u-face location at `(i, j)` (0-based, `0..nx × 0..ny-1`). */
  xEdge(i: number, j: number): [number, number];
  /** `(x, y)` of the v-face location at `(i, j)` (0-based, `0..nx-1 × 0..ny`). */
  yEdge(i: number, j: number): [number, number];
  /** `(x, y)` of the cell corner at `(i, j)` (0-based, `0..nx × 0..ny`). */
  corner(i: number, j: number): [number, number];
  /** Declarative `.esm`-shaped summary of the base grid. */
  toESM(): object;
}

export interface CartesianBaseOpts {
  xlo: number;
  xhi: number;
  ylo: number;
  yhi: number;
  nx: number;
  ny: number;
}

export interface CartesianBase extends ArakawaBaseGrid {
  readonly family: "cartesian";
  readonly xlo: number;
  readonly xhi: number;
  readonly ylo: number;
  readonly yhi: number;
}

function isPositiveInteger(v: unknown): v is number {
  return typeof v === "number" && Number.isInteger(v) && v >= 1;
}

function isFiniteNumber(v: unknown): v is number {
  return typeof v === "number" && Number.isFinite(v);
}

/** Build a validated uniform Cartesian base grid. */
export function cartesianBase(opts: CartesianBaseOpts): CartesianBase {
  if (opts === undefined || opts === null) {
    throw new TypeError("cartesianBase: options object is required");
  }
  const { xlo, xhi, ylo, yhi, nx, ny } = opts;
  if (!isPositiveInteger(nx)) {
    throw new RangeError(
      `cartesianBase: 'nx' must be a positive integer (got ${String(nx)})`,
    );
  }
  if (!isPositiveInteger(ny)) {
    throw new RangeError(
      `cartesianBase: 'ny' must be a positive integer (got ${String(ny)})`,
    );
  }
  if (!isFiniteNumber(xlo) || !isFiniteNumber(xhi)) {
    throw new RangeError(
      `cartesianBase: 'xlo' and 'xhi' must be finite numbers (got xlo=${String(xlo)}, xhi=${String(xhi)})`,
    );
  }
  if (!isFiniteNumber(ylo) || !isFiniteNumber(yhi)) {
    throw new RangeError(
      `cartesianBase: 'ylo' and 'yhi' must be finite numbers (got ylo=${String(ylo)}, yhi=${String(yhi)})`,
    );
  }
  if (!(xhi > xlo)) {
    throw new RangeError(
      `cartesianBase: 'xhi' must be > 'xlo' (got xlo=${xlo}, xhi=${xhi})`,
    );
  }
  if (!(yhi > ylo)) {
    throw new RangeError(
      `cartesianBase: 'yhi' must be > 'ylo' (got ylo=${ylo}, yhi=${yhi})`,
    );
  }

  const dxVal = (xhi - xlo) / nx;
  const dyVal = (yhi - ylo) / ny;

  const base: CartesianBase = {
    family: "cartesian",
    xlo,
    xhi,
    ylo,
    yhi,
    nx,
    ny,
    dx(): number {
      return dxVal;
    },
    dy(): number {
      return dyVal;
    },
    cellCenter(i: number, j: number): [number, number] {
      return [xlo + (i + 0.5) * dxVal, ylo + (j + 0.5) * dyVal];
    },
    xEdge(i: number, j: number): [number, number] {
      return [xlo + i * dxVal, ylo + (j + 0.5) * dyVal];
    },
    yEdge(i: number, j: number): [number, number] {
      return [xlo + (i + 0.5) * dxVal, ylo + j * dyVal];
    },
    corner(i: number, j: number): [number, number] {
      return [xlo + i * dxVal, ylo + j * dyVal];
    },
    toESM(): object {
      return {
        family: "cartesian",
        nx,
        ny,
        extent: [
          [xlo, ylo],
          [xhi, yhi],
        ],
      };
    },
  };
  return base;
}

// ---------------------------------------------------------------------------
// ArakawaGrid interface + options.
// ---------------------------------------------------------------------------

export interface ArakawaOpts {
  base: ArakawaBaseGrid;
  stagger: ArakawaStagger;
  dtype?: Dtype;
  ghosts?: number;
}

export interface ArakawaNeighbors {
  readonly west: [number, number] | null;
  readonly east: [number, number] | null;
  readonly south: [number, number] | null;
  readonly north: [number, number] | null;
}

export interface ArakawaGrid extends Grid {
  readonly family: "arakawa";
  readonly topology: "block_structured";
  readonly stagger: ArakawaStagger;
  readonly rotated: boolean;
  readonly ghosts: number;
  readonly base: ArakawaBaseGrid;
  readonly nx: number;
  readonly ny: number;
  readonly n_cells: number;

  /** Shape `(ni, nj)` of the named location on this grid. */
  location_shape(loc: ArakawaLocation): [number, number];
  /** Location where `var` lives under this grid's stagger. */
  variable_location(variable: ArakawaVariable): ArakawaLocation;
  /** Shape of the location where `var` lives under this stagger. */
  variable_shape(variable: ArakawaVariable): [number, number];

  /** `(x, y)` of cell center `(i, j)` (0-based). */
  cell_centers(i: number, j: number): [number, number];
  /** `(x, y)` of the u-variable location at `(i, j)` under this stagger. */
  u_face(i: number, j: number): [number, number];
  /** `(x, y)` of the v-variable location at `(i, j)` under this stagger. */
  v_face(i: number, j: number): [number, number];
  /** `(x, y)` of cell corner `(i, j)` (0-based, `0..nx × 0..ny`). */
  corners(i: number, j: number): [number, number];
  /** Generic `(x, y)` of `(loc, i, j)`. */
  coord(loc: ArakawaLocation, i: number, j: number): [number, number];

  /**
   * Four axial neighbours of `(loc, i, j)` as `{west, east, south, north}`,
   * each `[i', j']` or `null` at a domain boundary. Stays on the same
   * location (no cross-location stepping — that's the caller's job per
   * stagger semantics).
   */
  neighbors(loc: ArakawaLocation, i: number, j: number): ArakawaNeighbors;

  /** Evaluate a named metric at cell `(i, j)`. */
  metric_eval(name: ArakawaMetricName, i: number, j: number): number;
}

// ---------------------------------------------------------------------------
// Generator.
// ---------------------------------------------------------------------------

function checkBounds(
  loc: ArakawaLocation,
  i: number,
  j: number,
  nx: number,
  ny: number,
): void {
  const [ni, nj] = locationShape(loc, nx, ny);
  if (!Number.isInteger(i) || i < 0 || i >= ni) {
    throw new RangeError(
      `arakawa: i=${String(i)} out of range [0, ${ni}) for location '${loc}'`,
    );
  }
  if (!Number.isInteger(j) || j < 0 || j >= nj) {
    throw new RangeError(
      `arakawa: j=${String(j)} out of range [0, ${nj}) for location '${loc}'`,
    );
  }
}

function coordAt(
  base: ArakawaBaseGrid,
  loc: ArakawaLocation,
  i: number,
  j: number,
): [number, number] {
  switch (loc) {
    case "cell_center":
      return base.cellCenter(i, j);
    case "u_edge":
      return base.xEdge(i, j);
    case "v_edge":
      return base.yEdge(i, j);
    case "corner":
      return base.corner(i, j);
  }
}

function isBaseGridLike(v: unknown): v is ArakawaBaseGrid {
  if (v === null || typeof v !== "object") return false;
  const b = v as Record<string, unknown>;
  return (
    typeof b.nx === "number" &&
    typeof b.ny === "number" &&
    typeof b.dx === "function" &&
    typeof b.dy === "function" &&
    typeof b.cellCenter === "function" &&
    typeof b.xEdge === "function" &&
    typeof b.yEdge === "function" &&
    typeof b.corner === "function" &&
    typeof b.toESM === "function"
  );
}

export function arakawa(opts: ArakawaOpts): ArakawaGrid {
  if (opts === undefined || opts === null) {
    throw new TypeError("arakawa: options object is required");
  }

  if (!isBaseGridLike(opts.base)) {
    throw new TypeError(
      "arakawa: missing or invalid required option 'base' (expected an ArakawaBaseGrid)",
    );
  }
  if (opts.stagger === undefined) {
    throw new TypeError("arakawa: missing required option 'stagger'");
  }
  if (!STAGGERS.includes(opts.stagger)) {
    throw new RangeError(
      `arakawa: 'stagger' must be one of ${STAGGERS.join("/")} (got ${String(opts.stagger)})`,
    );
  }

  const dtype: Dtype = opts.dtype ?? "float64";
  if (dtype !== "float64" && dtype !== "float32") {
    throw new RangeError(
      `arakawa: 'dtype' must be 'float64' or 'float32' (got ${String(dtype)})`,
    );
  }

  const ghosts = opts.ghosts ?? 0;
  if (!Number.isInteger(ghosts) || ghosts < 0) {
    throw new RangeError(
      `arakawa: 'ghosts' must be a non-negative integer (got ${String(ghosts)})`,
    );
  }

  const base = opts.base;
  if (!isPositiveInteger(base.nx)) {
    throw new RangeError(
      `arakawa: base.nx must be a positive integer (got ${String(base.nx)})`,
    );
  }
  if (!isPositiveInteger(base.ny)) {
    throw new RangeError(
      `arakawa: base.ny must be a positive integer (got ${String(base.ny)})`,
    );
  }

  const stagger: ArakawaStagger = opts.stagger;
  const rotated = stagger === "E";
  const [hLoc, uLoc, vLoc] = variableLocations(stagger);

  const grid: ArakawaGrid = {
    family: "arakawa",
    topology: "block_structured",
    dtype,
    stagger,
    rotated,
    ghosts,
    base,
    get nx(): number {
      return base.nx;
    },
    get ny(): number {
      return base.ny;
    },
    get n_cells(): number {
      return base.nx * base.ny;
    },

    location_shape(loc: ArakawaLocation): [number, number] {
      if (!LOCATIONS.includes(loc)) {
        throw new RangeError(
          `arakawa: unknown location '${String(loc)}'; expected one of ${LOCATIONS.join("/")}`,
        );
      }
      return locationShape(loc, base.nx, base.ny);
    },

    variable_location(variable: ArakawaVariable): ArakawaLocation {
      switch (variable) {
        case "h":
          return hLoc;
        case "u":
          return uLoc;
        case "v":
          return vLoc;
        default:
          throw new RangeError(
            `arakawa: unknown variable '${String(variable)}'; expected one of ${VARIABLES.join("/")}`,
          );
      }
    },

    variable_shape(variable: ArakawaVariable): [number, number] {
      return locationShape(
        this.variable_location(variable),
        base.nx,
        base.ny,
      );
    },

    cell_centers(i: number, j: number): [number, number] {
      checkBounds("cell_center", i, j, base.nx, base.ny);
      return base.cellCenter(i, j);
    },

    u_face(i: number, j: number): [number, number] {
      checkBounds(uLoc, i, j, base.nx, base.ny);
      return coordAt(base, uLoc, i, j);
    },

    v_face(i: number, j: number): [number, number] {
      checkBounds(vLoc, i, j, base.nx, base.ny);
      return coordAt(base, vLoc, i, j);
    },

    corners(i: number, j: number): [number, number] {
      checkBounds("corner", i, j, base.nx, base.ny);
      return base.corner(i, j);
    },

    coord(loc: ArakawaLocation, i: number, j: number): [number, number] {
      if (!LOCATIONS.includes(loc)) {
        throw new RangeError(
          `arakawa: unknown location '${String(loc)}'; expected one of ${LOCATIONS.join("/")}`,
        );
      }
      checkBounds(loc, i, j, base.nx, base.ny);
      return coordAt(base, loc, i, j);
    },

    neighbors(
      loc: ArakawaLocation,
      i: number,
      j: number,
    ): ArakawaNeighbors {
      if (!LOCATIONS.includes(loc)) {
        throw new RangeError(
          `arakawa: unknown location '${String(loc)}'; expected one of ${LOCATIONS.join("/")}`,
        );
      }
      checkBounds(loc, i, j, base.nx, base.ny);
      const [ni, nj] = locationShape(loc, base.nx, base.ny);
      return {
        west: i > 0 ? [i - 1, j] : null,
        east: i + 1 < ni ? [i + 1, j] : null,
        south: j > 0 ? [i, j - 1] : null,
        north: j + 1 < nj ? [i, j + 1] : null,
      };
    },

    metric_eval(name: ArakawaMetricName, i: number, j: number): number {
      if (!METRICS.includes(name)) {
        throw new RangeError(
          `arakawa: unknown metric '${String(name)}'; expected one of ${METRICS.join("/")}`,
        );
      }
      checkBounds("cell_center", i, j, base.nx, base.ny);
      const dxV = base.dx();
      const dyV = base.dy();
      switch (name) {
        case "dx":
          return dxV;
        case "dy":
          return dyV;
        case "area":
          return dxV * dyV;
      }
    },

    toESM(): object {
      return {
        family: "arakawa",
        version: ARAKAWA_API_VERSION,
        dtype,
        topology: "block_structured",
        ghosts,
        n_cells: base.nx * base.ny,
        stagger,
        rotated,
        base: base.toESM(),
        provenance: {
          binding: "typescript",
          binding_version: ARAKAWA_BINDING_VERSION,
          api_version: ARAKAWA_API_VERSION,
          platform: "js",
          runtime: "node",
          math_lib: "libm-js",
          source_sha: "",
          stagger,
        },
      };
    },
  };

  return grid;
}
