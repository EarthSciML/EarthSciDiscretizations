/**
 * @earthsci/grids — Grid accessor runtime for the EarthSci model stack.
 *
 * Conforms to the cross-binding contract in `docs/GRIDS_API.md`.
 */

import * as grids from "./grids/index.js";

export { grids };
export type { Dtype, Grid } from "./grids/types.js";

export const earthsci = { grids };
