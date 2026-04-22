/**
 * Per-family grid accessor modules.
 *
 * Each family lands as its own module (e.g., `./cubed_sphere.ts`) and is
 * re-exported here per `docs/GRIDS_API.md` §2.6.
 */

export type { Dtype, Grid } from "./types.js";
export { cartesian } from "./cartesian.js";
export type {
  CartesianOpts,
  CartesianGrid,
  CartesianMetricName,
  CartesianNeighborRef,
  AxisExtent,
} from "./cartesian.js";
export { cubed_sphere } from "./cubed_sphere.js";
export type {
  CubedSphereOpts,
  CubedSphereGrid,
  Edge,
  NeighborRef,
  MetricName,
} from "./cubed_sphere.js";
export { duo } from "./duo.js";
export type {
  DuoOpts,
  DuoGrid,
  DuoLoader,
  DuoMetricName,
} from "./duo.js";
