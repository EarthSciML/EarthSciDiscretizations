#!/usr/bin/env node
/**
 * Regenerate the cross-binding golden for the ppm_reconstruction rule
 * conformance fixture.
 *
 * The script computes coefficients from the closed-form formulas the rule's
 * AST encodes — NOT via any binding's AST evaluator — so the golden is a
 * binding-independent reference. Each binding's PPM evaluator must reproduce
 * these values to within the fixture-declared tolerance.
 *
 * Closed-form (per discretizations/finite_volume/ppm_reconstruction.json,
 * CW84 eqs. 1.5-1.10):
 *   q_left_edge[i]  = (-q[i-2] + 7*q[i-1] + 7*q[i] -    q[i+1]) / 12
 *   q_right_edge[i] = (-q[i-1] + 7*q[i]   + 7*q[i+1] - q[i+2]) / 12
 *   parabola(i, xi) = a_L + xi * (da + a_6 * (1 - xi))
 *     da  = a_R - a_L
 *     a_6 = 6 * (q[i] - 0.5*(a_L + a_R))
 *
 * Cell averages come from the analytical antiderivative of the manufactured
 * solution (see fixtures.json), so quadrature noise is removed and the only
 * libm-induced variability is sin/cos at cell edges.
 *
 * Run from the repo root:
 *   node tests/conformance/grids/cartesian/rules/ppm_reconstruction/regenerate_golden.mjs
 */

import { readFileSync, writeFileSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const FIXTURE_PATH = join(HERE, "fixtures.json");
const GOLDEN_PATH = join(HERE, "golden.json");

const fixture = JSON.parse(readFileSync(FIXTURE_PATH, "utf8"));

const n = fixture.grid.n;
const [lo, hi] = fixture.grid.extent[0];
const dx = (hi - lo) / n;

// f(x) = sin(2*pi*x) + 0.3*cos(4*pi*x)
// F(x) = -cos(2*pi*x)/(2*pi) + 0.3*sin(4*pi*x)/(4*pi)
const TWO_PI = 2 * Math.PI;
const FOUR_PI = 4 * Math.PI;
const F = (x) =>
  -Math.cos(TWO_PI * x) / TWO_PI + (0.3 * Math.sin(FOUR_PI * x)) / FOUR_PI;
const cellAverage = (a, b) => (F(b) - F(a)) / (b - a);

const cellAverages = new Array(n);
for (let i = 0; i < n; i++) {
  cellAverages[i] = cellAverage(lo + i * dx, lo + (i + 1) * dx);
}

const wrap = (i) => ((i % n) + n) % n;

const qLeftEdgeAll = new Array(n);
const qRightEdgeAll = new Array(n);
for (let i = 0; i < n; i++) {
  const qm2 = cellAverages[wrap(i - 2)];
  const qm1 = cellAverages[wrap(i - 1)];
  const q0 = cellAverages[wrap(i)];
  const qp1 = cellAverages[wrap(i + 1)];
  const qp2 = cellAverages[wrap(i + 2)];
  qLeftEdgeAll[i] = (-qm2 + 7 * qm1 + 7 * q0 - qp1) / 12;
  qRightEdgeAll[i] = (-qm1 + 7 * q0 + 7 * qp1 - qp2) / 12;
}

const qLeftEdge = [];
const qRightEdge = [];
const qParabola = [];
for (const i of fixture.query_cells) {
  const aL = qLeftEdgeAll[i];
  const aR = qRightEdgeAll[i];
  qLeftEdge.push(aL);
  qRightEdge.push(aR);
  const qi = cellAverages[i];
  const da = aR - aL;
  const a6 = 6 * (qi - 0.5 * (aL + aR));
  const row = fixture.parabola_xi.map((xi) => aL + xi * (da + a6 * (1 - xi)));
  qParabola.push(row);
}

const out = {
  rule: fixture.rule,
  grid: fixture.grid,
  cell_averages: cellAverages,
  query_cells: fixture.query_cells,
  q_left_edge: qLeftEdge,
  q_right_edge: qRightEdge,
  parabola_xi: fixture.parabola_xi,
  q_parabola: qParabola,
};

writeFileSync(GOLDEN_PATH, JSON.stringify(out, null, 2) + "\n");
console.log(`wrote ${resolve(GOLDEN_PATH)}`);
