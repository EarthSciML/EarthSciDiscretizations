import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { grids } from "../src/index.js";
import type { VerticalOpts } from "../src/grids/vertical.js";

// Harness lives at <repo>/tests/conformance/grids/vertical — the TS
// package root is <repo>/typescript, so walk up one level.
const HARNESS_DIR = resolve(
  fileURLToPath(new URL("../", import.meta.url)),
  "..",
  "tests",
  "conformance",
  "grids",
  "vertical",
);

interface Fixture {
  name: string;
  opts: VerticalOpts;
}

interface FixtureSpec {
  tolerance: { relative: number };
  golden_dir: string;
  provenance_strip_fields: string[];
  fixtures: Fixture[];
}

const spec: FixtureSpec = JSON.parse(
  readFileSync(join(HARNESS_DIR, "fixtures.json"), "utf8"),
);

const GOLDEN_DIR = resolve(HARNESS_DIR, spec.golden_dir);
const STRIP_FIELDS = new Set(spec.provenance_strip_fields);

// Produce a JSON string with keys sorted recursively and 2-space indent.
// Matches Python's `json.dumps(..., sort_keys=True, indent=2)` so the same
// golden file canonicalizes identically from both bindings. The only
// additional discipline required is for numbers: Python emits `100000.0`
// while the default `JSON.stringify` emits `100000`. See `formatNumber`.
function sortValue(v: unknown): unknown {
  if (Array.isArray(v)) return v.map(sortValue);
  if (v !== null && typeof v === "object") {
    const entries = Object.entries(v as Record<string, unknown>).sort(
      ([a], [b]) => (a < b ? -1 : a > b ? 1 : 0),
    );
    const out: Record<string, unknown> = {};
    for (const [k, val] of entries) out[k] = sortValue(val);
    return out;
  }
  return v;
}

// Render a finite number using Python's `repr(float)` rules so a ts-side
// JSON dump matches the python-generated golden byte-for-byte. Python
// always renders floats with a decimal point (`100000.0`, not `100000`)
// and uses `repr` (shortest round-trip) for non-integer values.
// JS's `Number.prototype.toString()` is spec'd to emit the same shortest
// round-trip representation as Python's `repr` for non-integer finite
// doubles, so the only gap is the trailing ".0" on integer-valued floats.
function formatNumber(n: number): string {
  if (!Number.isFinite(n)) {
    throw new Error(
      `vertical conformance canonicalization: non-finite number ${n}`,
    );
  }
  const s = n.toString();
  // Integer-valued floats: append ".0" to match python's repr.
  if (/^-?\d+$/.test(s)) return `${s}.0`;
  return s;
}

function serialize(v: unknown, indent: number): string {
  if (v === null) return "null";
  if (typeof v === "boolean") return v ? "true" : "false";
  if (typeof v === "number") return formatNumber(v);
  if (typeof v === "string") return JSON.stringify(v);
  if (Array.isArray(v)) {
    if (v.length === 0) return "[]";
    const inner = " ".repeat(indent + 2);
    const parts = v.map((x) => inner + serialize(x, indent + 2));
    return `[\n${parts.join(",\n")}\n${" ".repeat(indent)}]`;
  }
  if (typeof v === "object") {
    const entries = Object.entries(v as Record<string, unknown>);
    if (entries.length === 0) return "{}";
    const inner = " ".repeat(indent + 2);
    const parts = entries.map(
      ([k, val]) =>
        `${inner}${JSON.stringify(k)}: ${serialize(val, indent + 2)}`,
    );
    return `{\n${parts.join(",\n")}\n${" ".repeat(indent)}}`;
  }
  throw new Error(
    `vertical conformance canonicalization: unsupported value ${typeof v}`,
  );
}

function canonicalize(doc: Record<string, unknown>): string {
  const stripped: Record<string, unknown> = { ...doc };
  const prov = stripped.provenance;
  if (prov !== undefined && prov !== null && typeof prov === "object") {
    const filtered: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(prov as Record<string, unknown>)) {
      if (!STRIP_FIELDS.has(k)) filtered[k] = v;
    }
    stripped.provenance = filtered;
  }
  return serialize(sortValue(stripped), 0) + "\n";
}

describe("vertical cross-language conformance", () => {
  it("fixtures.json declares byte-equality (tolerance.relative === 0)", () => {
    expect(spec.tolerance.relative).toBe(0);
  });

  for (const fixture of spec.fixtures) {
    it(`matches the committed golden (fixture=${fixture.name})`, () => {
      const g = grids.vertical(fixture.opts);
      const emitted = canonicalize(g.toESM() as Record<string, unknown>);

      const goldenText = readFileSync(
        join(GOLDEN_DIR, `${fixture.name}.esm`),
        "utf8",
      );
      const golden = JSON.parse(goldenText) as Record<string, unknown>;
      const expected = canonicalize(golden);

      expect(emitted).toBe(expected);
    });
  }
});
