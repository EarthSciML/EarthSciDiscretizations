/**
 * Shared types for the grids namespace.
 * See `docs/GRIDS_API.md` §3.4, §7.
 */

export type Dtype = "float64" | "float32";

export interface Grid {
  readonly family: string;
  readonly dtype: Dtype;
  toESM(): object;
}
