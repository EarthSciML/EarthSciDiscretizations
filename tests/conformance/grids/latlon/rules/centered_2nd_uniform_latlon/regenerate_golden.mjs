#!/usr/bin/env node
/**
 * Regenerate the golden file for the centered_2nd_uniform_latlon rule
 * conformance fixture.
 *
 * The script computes coefficients from the closed-form formulas the rule's
 * AST encodes — NOT via any binding's AST evaluator — so the golden is a
 * binding-independent reference. Each binding's evaluator must reproduce
 * these values to within the fixture-declared tolerance.
 *
 * Coefficients (per discretizations/finite_difference/centered_2nd_uniform_latlon.json):
 *   lon, offset ±1: ±1 / (2 * R * cos_lat * dlon)
 *   lat, offset ±1: ±1 / (2 * R * dlat)
 *
 * Run from the repo root:
 *
 *   node tests/conformance/grids/latlon/rules/centered_2nd_uniform_latlon/regenerate_golden.mjs
 */

import { readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve, join } from "node:path";

const HERE = dirname(fileURLToPath(import.meta.url));
const FIXTURE_PATH = join(HERE, "fixtures.json");
const GOLDEN_PATH = join(HERE, "golden.json");

const fixture = JSON.parse(readFileSync(FIXTURE_PATH, "utf8"));
const grid = fixture.grid;

const dlon = (2 * Math.PI) / grid.nlon;
const dlat = Math.PI / grid.nlat;
const lonStart = -Math.PI;

function cellCenter(j, i) {
  const lat = -Math.PI / 2 + (j + 0.5) * dlat;
  const lon = lonStart + (i + 0.5) * dlon;
  return { lon, lat };
}

const cell_centers = fixture.query_points.map(([j, i]) => cellCenter(j, i));
const cos_lat = cell_centers.map((c) => Math.cos(c.lat));

const coeffs = {};
for (const key of fixture.stencil_order) {
  const tag = `${key.axis}_${key.offset >= 0 ? "p" : "m"}${Math.abs(key.offset)}`;
  coeffs[tag] = fixture.query_points.map(([j, _i], k) => {
    void _i;
    const cl = cos_lat[k];
    const sign = key.offset > 0 ? 1 : -1;
    if (key.axis === "lon") return sign / (2 * grid.R * cl * dlon);
    if (key.axis === "lat") return sign / (2 * grid.R * dlat);
    throw new Error(`unsupported axis ${key.axis}`);
  });
}

const out = {
  rule: "centered_2nd_uniform_latlon",
  grid,
  query_points: fixture.query_points,
  cell_centers,
  cos_lat,
  coeffs,
};

writeFileSync(GOLDEN_PATH, JSON.stringify(out, null, 2) + "\n");
console.log(`wrote ${resolve(GOLDEN_PATH)}`);
