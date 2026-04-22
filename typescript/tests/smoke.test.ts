import { describe, it, expect } from "vitest";
import { earthsci, grids } from "../src/index.js";
import type { Dtype, Grid } from "../src/index.js";

describe("@earthsci/grids smoke", () => {
  it("exposes the grids namespace via both named and earthsci export", () => {
    expect(earthsci.grids).toBe(grids);
  });

  it("exports Dtype and Grid types", () => {
    const dtype: Dtype = "float64";
    expect(dtype).toBe("float64");

    const stub: Grid = {
      family: "stub",
      dtype: "float64",
      toESM: () => ({}),
    };
    expect(stub.family).toBe("stub");
  });
});
