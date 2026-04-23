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
export { lat_lon } from "./lat_lon.js";
export type {
  LatLonOpts,
  LatLonGrid,
  LatLonVariant,
  LatLonMetricName,
  LatLonNeighborCell,
  LatLonNeighborSet,
  PolePolicy,
  Direction,
} from "./lat_lon.js";
export { duo } from "./duo.js";
export type {
  DuoOpts,
  DuoGrid,
  DuoLoader,
  DuoMetricName,
} from "./duo.js";
export { mpas, mpasMeshData, checkMesh } from "./mpas.js";
export type {
  MpasOpts,
  MPASGrid,
  MpasLoader,
  MpasReader,
  MpasCheck,
  MpasMeshInput,
  MPASMeshData,
  MpasMetricName,
} from "./mpas.js";
export {
  arakawa,
  cartesianBase,
  variableLocations,
  locationShape,
} from "./arakawa.js";
export type {
  ArakawaOpts,
  ArakawaGrid,
  ArakawaStagger,
  ArakawaLocation,
  ArakawaVariable,
  ArakawaMetricName,
  ArakawaNeighbors,
  ArakawaBaseGrid,
  CartesianBase,
  CartesianBaseOpts,
} from "./arakawa.js";
