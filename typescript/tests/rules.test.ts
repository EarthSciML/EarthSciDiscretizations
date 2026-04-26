import { describe, it, expect } from "vitest";

import {
  evaluate,
  evalCoeff,
  UnboundVariableError,
  parseRuleFile,
  findRule,
} from "../src/rules/index.js";
import type { ExpressionNode } from "../src/rules/index.js";

describe("rule expression evaluator", () => {
  it("evaluates numeric literals", () => {
    expect(evaluate(3.14, new Map())).toBe(3.14);
    expect(evaluate(-2, new Map())).toBe(-2);
  });

  it("evaluates variables via bindings", () => {
    const b = new Map([["x", 5]]);
    expect(evaluate("x", b)).toBe(5);
  });

  it("throws UnboundVariableError for missing variables", () => {
    expect(() => evaluate("missing", new Map())).toThrow(UnboundVariableError);
  });

  it("evaluates arithmetic ops including unary minus and n-ary +/*", () => {
    expect(evaluate({ op: "+", args: [1, 2, 3, 4] }, new Map())).toBe(10);
    expect(evaluate({ op: "*", args: [2, 3, 4] }, new Map())).toBe(24);
    expect(evaluate({ op: "-", args: [10, 3] }, new Map())).toBe(7);
    expect(evaluate({ op: "-", args: [5] }, new Map())).toBe(-5);
    expect(evaluate({ op: "/", args: [1, 4] }, new Map())).toBe(0.25);
    expect(evaluate({ op: "^", args: [2, 10] }, new Map())).toBe(1024);
  });

  it("evaluates the centered_2nd_uniform_latlon lon coefficient", () => {
    // -1 / (2 * R * cos_lat * dlon)
    const node: ExpressionNode = {
      op: "/",
      args: [
        -1,
        { op: "*", args: [2, "R", "cos_lat", "dlon"] },
      ],
    };
    const R = 1.0;
    const cos_lat = Math.cos(Math.PI / 8);
    const dlon = Math.PI / 4;
    const expected = -1 / (2 * R * cos_lat * dlon);
    const got = evalCoeff(node, { R, cos_lat, dlon });
    expect(got).toBe(expected);
  });

  it("rejects division by zero", () => {
    expect(() => evaluate({ op: "/", args: [1, 0] }, new Map())).toThrow(
      "Division by zero",
    );
  });

  it("supports trig and ifelse for completeness", () => {
    expect(evaluate({ op: "sin", args: [0] }, new Map())).toBe(0);
    expect(evaluate({ op: "cos", args: [0] }, new Map())).toBe(1);
    expect(evaluate({ op: "ifelse", args: [1, 42, 99] }, new Map())).toBe(42);
    expect(evaluate({ op: "ifelse", args: [0, 42, 99] }, new Map())).toBe(99);
  });
});

describe("rule loader", () => {
  it("parses the centered_2nd_uniform_latlon rule from the on-disk catalog", () => {
    // Hand-craft the same JSON content the loader would read so the test
    // doesn't depend on filesystem layout in unit-test scope.
    const content = JSON.stringify({
      discretizations: {
        centered_2nd_uniform_latlon: {
          applies_to: { op: "grad", args: ["$u"], dim: "$k" },
          grid_family: "latlon",
          combine: "+",
          accuracy: "O(h^2)",
          stencil: [
            {
              selector: { kind: "latlon", axis: "lon", offset: -1 },
              coeff: {
                op: "/",
                args: [-1, { op: "*", args: [2, "R", "cos_lat", "dlon"] }],
              },
            },
          ],
        },
      },
    });
    const rules = parseRuleFile(
      "/tmp/finite_difference/centered_2nd_uniform_latlon.json",
      content,
    );
    expect(rules.length).toBe(1);
    const rule = findRule(rules, "centered_2nd_uniform_latlon");
    expect(rule).toBeDefined();
    expect(rule!.family).toBe("finite_difference");
    expect(rule!.grid_family).toBe("latlon");
    expect(rule!.stencil.length).toBe(1);
  });
});
