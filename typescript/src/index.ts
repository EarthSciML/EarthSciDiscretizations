/**
 * @earthsci/grids — Grid accessor runtime for the EarthSci model stack.
 *
 * Conforms to the cross-binding contract in `docs/GRIDS_API.md`.
 */

import * as grids from "./grids/index.js";
import * as rules from "./rules/index.js";

export { grids, rules };
export type { Dtype, Grid } from "./grids/types.js";

export const earthsci = { grids, rules };
