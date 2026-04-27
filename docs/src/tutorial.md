# Tutorial: Authoring a rule

This tutorial walks through writing a new discretization rule end-to-end
in the closed-AST lowering pattern. The running example is
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}) — the
two-point centered finite difference for ∂u/∂x on a uniform Cartesian
axis. It is the smallest rule in the catalog that exercises every step,
and it is currently the canonical linear exemplar (commit `7b26ffd`).

By the end you will have:

1. Pattern-matched a §4.2 PDE operator with metavariables.
2. Lowered it to a closed `arrayop` expression in the §4.2 op vocabulary.
3. Delegated boundary handling to the domain's `boundary_conditions`
   block.
4. Validated the rule with a Layer-A canonical-form fixture.
5. Set up the Layer-B convergence sweep with a documented escape hatch
   for the in-flight ESS prerequisite.

## Step 1 — Match a §4.2 PDE operator

A rule's `applies_to` clause is a pattern that the ESS rewriter unifies
against the equation tree. The `op` field MUST be one of the §4.2 PDE
operators: `D`, `grad`, `div`, `laplacian` (plus the pointwise-math
vocabulary for non-spatial transformations).

Do **not** invent off-spec match keys. Names like `advect`,
`reconstruct`, `flux`, and `limit` are forbidden as `applies_to.op`
values — those are *implementation details of a scheme*, not PDE
operators. Every advection scheme matches a `D` (or, in flux form, a
`div`); every reconstruction is part of the `arrayop` body of the rule
that emits the reconstructed array.

For our example, the operator is the spatial gradient:

```json
"applies_to": {
  "op":   "grad",
  "args": ["$u"],
  "dim":  "$x"
}
```

`$u` and `$x` are pattern metavariables. They bind to the actual field
name and axis at rewrite time, so the same rule fires for
`grad(temperature, dim=x)` and `grad(salinity, dim=x)` alike. If a rule
applies to multiple PDE ops, write one rule file per pattern — the
catalog favours small, single-purpose rules.

## Step 2 — Lower to a closed `arrayop`

The `replacement` is a single AST node in the §4.2 op vocabulary.
`arrayop` is the usual top-level shape: it carries an `output_idx`
(symbolic indices spanning the result), an `expr` (the scalar body
evaluated at each index point), and `args` (the input operands).

For the centered difference, the body is `(u[$x+1] - u[$x-1]) / (2·dx)`:

```json
"replacement": {
  "op": "arrayop",
  "output_idx": ["$x"],
  "expr": {
    "op": "/",
    "args": [
      { "op": "-", "args": [
        { "op": "index", "args": ["$u", { "op": "+", "args": ["$x", 1] }] },
        { "op": "index", "args": ["$u", { "op": "-", "args": ["$x", 1] }] }
      ]},
      { "op": "*", "args": [2, "dx"] }
    ]
  },
  "args": ["$u"]
}
```

A few authoring rules that follow from "stay in the §4.2 vocabulary":

- **No scheme-specific kernels.** Every `op` in the body is one of the
  §4.2 ops listed in [Operators](@ref). If you need a limiter, write
  `min` / `max` / `ifelse`. If you need a polynomial reconstruction,
  write `+` / `*` / `index`. Reviewers reject `fn` nodes that
  re-implement a clamp under a custom name.
- **No `bc:*` ops in the body.** The lowering is the *interior* closed
  form; boundaries come in at Step 3.
- **No host-language code anywhere.** Rules ship JSON. The reference
  evaluator implements §4.2 once; rules compose against it.

Larger rules follow the same shape. A flux-form lowering uses an
`arrayop` whose `output_idx` ranges over edges; a limiter is a
`broadcast` of `min`/`max`/`ifelse` over operand arrays; a 5-point
Laplacian is an `arrayop` summing five `index` nodes with the
appropriate coefficients in the body. The full op alphabet is in
[Operators](@ref).

## Step 3 — Delegate boundary conditions to the domain

Boundary handling is declared once per field on the domain's
`boundary_conditions` block (`esm-spec.md` §11.5). The rule itself stays
BC-agnostic — its lowering covers only the interior. A separate set of
**downstream BC rewrite rules** consumes the domain's BC list and
rewrites the index expressions at the boundary cells:

| Domain BC | Index transformation applied to `$u[$x ± 1]` |
|---|---|
| `periodic` | wrap-around: `mod($x ± 1 + N, N)` (see [`periodic_bc`]({{< ref "/rules/periodic_bc" >}})) |
| `dirichlet` / `constant` | boundary cell reads the prescribed value |
| `neumann` / `zero_gradient` | mirror the in-range neighbor (clamp the index) |
| `robin` | mixed coefficient row at the boundary |

Concretely, this means **do not** write `ifelse` branches in your
`arrayop` body to special-case the boundary cells. Let the lowering be
the one-line interior formula; the BC rewriter will rewrite the index
expressions at the boundary at lowering time.

## Step 4 — Layer-A canonical fixture (validate the JSON)

Layer A is the canonical-form round-trip: load the rule JSON, walk the
AST, re-serialize, and compare. It catches schema violations,
metavariable typos, and accidental host-language leakage. Every rule
ships at least one Layer-A fixture covering the lowering.

The Layer-A fixture for our running example lives next to the rule:

```
discretizations/finite_difference/centered_2nd_uniform/fixtures/canonical/
```

Run the catalog tests locally to exercise it:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The catalog tests are the *defensive* layer — they confirm the rule
parses, the metavariables resolve, and the lowering walks cleanly. They
should pass before you proceed to Step 5.

## Step 5 — Layer-B convergence sweep (and the ESS prerequisite)

Layer B is an MMS (manufactured-solution) convergence sweep that
verifies the rule's empirical order of accuracy on a refinement
sequence. Each rule ships an `input.esm` that names the manufactured
solution and a sweep of grid sizes; the harness measures L∞ / L₂ error
and fits a slope.

There is an **in-flight prerequisite** worth flagging: the ESS
`mms_evaluator` currently dispatches via the legacy `spec['stencil']`
kernels rather than walking the closed `arrayop` lowering. Until ESS
gains an AST-walker dispatch path (tracked as ESS `esm-4gw`), Layer B
cannot exercise a closed-AST rule against the symbolic harness.

The right escape hatch for this **specific** prerequisite is
`applicable: false` with a `skip_reason` that names the blocker:

```json
{
  "rule": "centered_2nd_uniform",
  "manufactured_solution": "sin(2*pi*x) on [0,1] periodic; derivative 2*pi*cos(2*pi*x)",
  "sampling": "cell_center",
  "grids": [
    { "n": 16 }, { "n": 32 }, { "n": 64 }, { "n": 128 }
  ],
  "applicable": false,
  "skip_reason": "rule rewritten as closed arrayop replacement (dsc-rar); ESS mms_evaluator currently dispatches via spec['stencil'] kernels — Layer-B re-enables once ESS gains an AST-walker dispatch (see follow-up bead)."
}
```

`applicable: false` is **only** for blockers of this kind: an upstream
prerequisite that, once landed, lets the *unmodified* fixture re-enable.
It is **not** a general escape hatch for "the convergence test is hard
to write", "we haven't verified the order yet", or "the MMS choice is
wrong". A `skip_reason` that doesn't name a tracked blocker should not
land.

When ESS lands `esm-4gw`, the fixture flips to `applicable: true`
without other edits and the sweep exercises the rule's `arrayop`
lowering through the same evaluator the canonical tests use.

## Step 6 — Document and link

Each rule has a doc page under `docs/content/rules/<rule>.md`. Follow
the structure of
[`centered_2nd_uniform`]({{< ref "/rules/centered_2nd_uniform" >}}) —
overview, `applies_to` and `replacement` AST, BC handoff table,
truncation derivation, convergence figure / status. Cross-link to
related rules and to `esm-spec.md` §4.2 / §11.5 for the definitive
operator and BC vocabulary.

The catalog landing page at
[`docs/content/rules/_index.md`]({{< ref "/rules" >}}) advertises the
closed-AST lowering pattern as the default. New rules should match its
framing; rules predating the migration carry a "legacy form" note on
their page until they are rewritten.

## Adapt to your scheme

To author a different rule, swap each piece:

- Step 1 — choose the §4.2 PDE op your scheme discretizes. (`D` for
  most time-dependent schemes; `grad` / `div` / `laplacian` for spatial
  operators.)
- Step 2 — write the closed `arrayop` body. Use `index`, arithmetic,
  `min` / `max` / `ifelse`, and (for nonlinear or weighted schemes)
  `broadcast`. No new ops.
- Step 3 — leave BC handling to the domain. Do not embed BC switches in
  the body.
- Step 4 — ship a canonical fixture. Make it pass.
- Step 5 — ship the convergence fixture. Use `applicable: false` with a
  named-blocker `skip_reason` only if a tracked prerequisite genuinely
  blocks the sweep; otherwise the fixture must be `applicable: true`
  and the slope must match the declared accuracy.
- Step 6 — write the doc page in the same shape as
  `centered_2nd_uniform.md` and cross-link.

For a contributor-oriented walk through the surrounding repository
infrastructure (paths, CI layers, registration), see the
[Add a new discretization rule]({{< ref "/tutorials/add-a-rule" >}})
companion tutorial — note that until that tutorial migrates to the new
pedagogy it still describes the legacy stencil/coefficient form on some
steps.
