/**
 * Unit tests for the ppm_reconstruction rule evaluator (sub_stencil and
 * output_kind dispatch, parabola formula, periodic indexing).
 */

import { describe, it, expect } from "vitest";

import { parseRuleFile } from "../src/rules/index.js";
import type { Rule } from "../src/rules/index.js";
import {
  applyEdgeStencil,
  assertPpmRule,
  evaluateParabola,
  getSubStencil,
  reconstructCell,
  reconstructOutput,
  resolveAxis,
} from "./ppm_reconstruction_helpers.js";

const PPM_JSON = JSON.stringify({
  discretizations: {
    ppm_reconstruction: {
      applies_to: { op: "reconstruct", args: ["$q"], dim: "$x" },
      grid_family: "cartesian",
      form: "piecewise_parabolic",
      limiter: "none",
      accuracy: "O(dx^3)",
      outputs: ["q_left_edge", "q_right_edge", "q_parabola"],
      stencil: {
        q_left_edge: [
          {
            selector: { kind: "cartesian", axis: "$x", offset: -2 },
            coeff: { op: "/", args: [-1, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: -1 },
            coeff: { op: "/", args: [7, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: 0 },
            coeff: { op: "/", args: [7, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: 1 },
            coeff: { op: "/", args: [-1, 12] },
          },
        ],
        q_right_edge: [
          {
            selector: { kind: "cartesian", axis: "$x", offset: -1 },
            coeff: { op: "/", args: [-1, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: 0 },
            coeff: { op: "/", args: [7, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: 1 },
            coeff: { op: "/", args: [7, 12] },
          },
          {
            selector: { kind: "cartesian", axis: "$x", offset: 2 },
            coeff: { op: "/", args: [-1, 12] },
          },
        ],
      },
    },
  },
});

function loadRule(): Rule {
  const set = parseRuleFile(
    "/tmp/finite_volume/ppm_reconstruction.json",
    PPM_JSON,
  );
  if (set.length !== 1 || set[0].name !== "ppm_reconstruction") {
    throw new Error("test rule fixture failed to parse");
  }
  return set[0];
}

describe("ppm_reconstruction rule contract", () => {
  it("assertPpmRule accepts the canonical multi-stencil shape", () => {
    expect(() => assertPpmRule(loadRule())).not.toThrow();
  });

  it("assertPpmRule rejects a flat-stencil rule", () => {
    const bad: Rule = {
      ...loadRule(),
      stencil: [
        {
          selector: { kind: "cartesian", axis: "$x", offset: 0 },
          coeff: 1,
        },
      ],
    };
    expect(() => assertPpmRule(bad)).toThrow(/multi-stencil/);
  });

  it("resolveAxis strips the '$' from applies_to.dim", () => {
    expect(resolveAxis(loadRule())).toBe("x");
  });

  it("getSubStencil returns named entries and throws on missing name", () => {
    const rule = loadRule();
    expect(getSubStencil(rule, "q_left_edge").length).toBe(4);
    expect(getSubStencil(rule, "q_right_edge").length).toBe(4);
    expect(() => getSubStencil(rule, "q_unknown")).toThrow(/sub-stencil/);
  });
});

describe("evaluateParabola (CW84 eqs. 1.5,1.7,1.10)", () => {
  it("reproduces a_L at xi=0 and a_R at xi=1", () => {
    expect(evaluateParabola(2, 5, 3, 0)).toBeCloseTo(2, 15);
    expect(evaluateParabola(2, 5, 3, 1)).toBeCloseTo(5, 15);
  });

  it("integrates back to the cell average over xi in [0,1] (Simpson, 3 nodes)", () => {
    // Simpson's rule with 3 nodes is exact for quadratics, which the parabola is.
    const aL = -0.4;
    const aR = 1.7;
    const qbar = 0.65;
    const v0 = evaluateParabola(aL, aR, qbar, 0);
    const v1 = evaluateParabola(aL, aR, qbar, 0.5);
    const v2 = evaluateParabola(aL, aR, qbar, 1);
    const integral = (v0 + 4 * v1 + v2) / 6;
    expect(integral).toBeCloseTo(qbar, 14);
  });

  it("rejects xi outside [0,1]", () => {
    expect(() => evaluateParabola(0, 1, 0.5, -0.01)).toThrow(/xi=/);
    expect(() => evaluateParabola(0, 1, 0.5, 1.01)).toThrow(/xi=/);
  });
});

describe("applyEdgeStencil and reconstructOutput dispatch", () => {
  // 8-cell periodic field with predictable values.
  const q = [1, 2, 3, 4, 5, 6, 7, 8];

  it("computes q_right_edge[i] = (-q[i-1] + 7*q[i] + 7*q[i+1] - q[i+2])/12", () => {
    const rule = loadRule();
    for (let i = 0; i < q.length; i++) {
      const im1 = q[(i - 1 + q.length) % q.length];
      const i0 = q[i];
      const ip1 = q[(i + 1) % q.length];
      const ip2 = q[(i + 2) % q.length];
      const expected = (-im1 + 7 * i0 + 7 * ip1 - ip2) / 12;
      expect(applyEdgeStencil(rule, "q_right_edge", q, i)).toBeCloseTo(
        expected,
        14,
      );
    }
  });

  it("computes q_left_edge as q_right_edge of the previous cell", () => {
    const rule = loadRule();
    for (let i = 0; i < q.length; i++) {
      const left = applyEdgeStencil(rule, "q_left_edge", q, i);
      const prevRight = applyEdgeStencil(
        rule,
        "q_right_edge",
        q,
        (i - 1 + q.length) % q.length,
      );
      expect(left).toBeCloseTo(prevRight, 14);
    }
  });

  it("reconstructCell.parabola samples the right edge values", () => {
    const rule = loadRule();
    const cell = reconstructCell(rule, q, 3);
    expect(cell.parabola(0)).toBeCloseTo(cell.q_left_edge, 14);
    expect(cell.parabola(1)).toBeCloseTo(cell.q_right_edge, 14);
  });

  it("reconstructOutput dispatches q_parabola only when xi is provided", () => {
    const rule = loadRule();
    expect(() => reconstructOutput(rule, "q_parabola", q, 0)).toThrow(
      /options\.xi/,
    );
    const v = reconstructOutput(rule, "q_parabola", q, 0, { xi: 0.5 });
    const cell = reconstructCell(rule, q, 0);
    expect(v).toBe(cell.parabola(0.5));
  });

  it("rejects non-integer cell indices", () => {
    const rule = loadRule();
    expect(() => applyEdgeStencil(rule, "q_left_edge", q, 1.5)).toThrow(
      /must be an integer/,
    );
  });

  it("rejects evaluation against an empty field", () => {
    const rule = loadRule();
    expect(() => applyEdgeStencil(rule, "q_left_edge", [], 0)).toThrow(/empty/);
  });
});
