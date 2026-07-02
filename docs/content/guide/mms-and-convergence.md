---
title: "MMS and convergence"
description: "How inline §6.6.5 tests, convergence manifests, and committed goldens pin correctness and observed order of accuracy."
weight: 3
---

Two questions get tested for every rule: *is the discretization right* (MMS,
at one resolution) and *does it converge at the claimed order* (the sweep).
Both start from a manufactured solution and neither ever computes a number in
this repository — all numerics run through official ESS binding runners.

## Inline MMS tests (esm-spec §6.6.5)

Each driver under `problems/` is a complete model with `tests` inline in the
`.esm` file itself. A test declares a `time_span` and `assertions`; each
assertion compares a state variable against a **reference expression** at a
given time under a `reduce` norm:

```json
{
  "variable": "u", "time": 0.1,
  "expected": 0.0, "tolerance": { "abs": 1e-4 },
  "reduce": "L2_error",
  "reference": { "op": "*", "args": [
    { "op": "exp", "args": [ { "op": "*", "args": [-0.01, 9.869604401089358, 0.1] } ] },
    { "op": "cos", "args": [ { "op": "*", "args": [3.141592653589793,
      { "op": "apply_expression_template", "args": [], "name": "x_coord", "bindings": {} } ] } ] } ] }
}
```

The method of manufactured solutions (Roache 2002) picks a solution first and
derives the problem from it. The solution must be **compatible with the rule's
boundary condition**: `problems/heat_1d_zero_grad.esm` uses
`u(x,t) = exp(−απ²t)·cos(πx)` because `cos(πx)` has zero gradient at both ends
— it satisfies the continuous PDE exactly, so no forcing term is needed and
every recorded error is pure spatial-discretization error. The derivation and
the BC-compatibility argument live in the test's `description` (that is the
required home for them, per AGENTS.md §4).

References are written explicitly over the grid's geometry templates
(`x_coord`), so no implicit index-to-coordinate mapping is assumed anywhere.

A `tests/conformance/simulation/<case>/manifest.json` routes each binding's
official simulation pathway at the problem's default resolution; the
assertions and tolerances live in the problem file itself.

## Convergence manifests

`tests/conformance/convergence/<case>/manifest.json` declares the sweep:

```json
{
  "problem": "../../../../problems/heat_1d_zero_grad.esm",
  "rule": "../../../../grids/cartesian_uniform_1d/rules/central_D2_zero_grad_bc.esm",
  "resolutions": [ { "n": 16, "bindings": { "N": 16 } },
                   { "n": 32, "bindings": { "N": 32 } }, "…" ],
  "norms": ["L2_error", "Linf_error"],
  "expected_order": 2.0, "order_tolerance": 0.2,
  "reference_binding": "julia",
  "golden": "golden/errors.json",
  "status": "pending-golden"
}
```

The problem file is loaded once per resolution with the **loader-API
metaparameter binding** (esm-spec §9.7.6, binding site 4). No per-resolution
fixture files exist or are needed — if you feel the need for a fixture
generator, the metaparameter mechanism is being misused.

## Goldens and the observed order

`scripts/regenerate-goldens.sh` (Julia-reference-only) runs the sweep through
the canonical pipeline and commits the error norms to `golden/errors.json`:

```json
{ "case": "heat_1d_zero_grad", "binding": "julia",
  "errors": [ { "n": 16, "L2_error": 1.2e-4, "Linf_error": 1.9e-4 }, "…" ] }
```

From then on everything *reads* that file:

- `scripts/check_convergence_order.py` (stdlib-only, runs in the fast CI
  validate job) computes the pairwise observed orders
  `p_k = log(e_{k−1}/e_k) / log(n_k/n_{k−1})` and asserts
  `|median(p) − expected_order| ≤ order_tolerance`.
- Other bindings are compared against the golden values under
  `tolerances.error_vs_golden_*` — they must reproduce the *same errors*, not
  merely converge.
- The docs build renders the observed-order table and the log-log convergence
  plot on the rule's grid page from the same file. Nothing is ever recomputed
  for display.

Until the golden lands, the manifest carries `"status": "pending-golden"`:
the order check reports SKIP (a missing golden without that marker is a
failure), and the grid page shows a pending callout in place of the table.
