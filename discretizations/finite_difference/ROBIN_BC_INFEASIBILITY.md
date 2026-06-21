# robin_bc (§9.2/§9.4) — DECLARATIVE-OR-FAIL verdict (G9 / Robin BC numeric)

**Bead:** esd-6g4.8 (G9) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Shares root cause with:** esd-6g4.9 (G10, nonzero-Neumann) — gaps #2 and #3
below apply verbatim to it; gap #1 is Robin-specific (3 coefficients vs 1 value).
**ESS source pin:** EarthSciSerialization GitHub `main` @
`6ebe038f13d7d3b0e1481d07afbfc835dd3b8238` — the rev ESD `Project.toml`
`[sources]` resolves (`{rev = "main", subdir = "packages/EarthSciSerialization.jl"}`).
All ESS line numbers below are at that commit. The BC-relevant `discretize.jl`,
`parse.jl`, and `types.jl` are byte-identical in the depot tree that
`Pkg.instantiate()` resolved for the empirical probes (`EarthSciSerialization`
v0.6.0).

---

## 0. Verdict (TL;DR)

> **DECLARATIVE — INFEASIBLE over the existing ESS engine.** The Robin boundary
> condition `a·u + b·∂u/∂n = g` cannot be authored as a declarative rule that the
> **existing** ESS engine *evaluates* in the production discretization pipeline
> (`build_ode_problem`). Per the bead's `DECLARATIVE-OR-FAIL` contract, the
> correct outcome is this precise gap report + escalation — **no firing rule is
> committed and no imperative ghost operator is hand-coded.** The existing
> `robin_bc.json` stays as a non-firing schema placeholder; the
> `@test_throws RuleEngineError` (`E_BC_UNSUPPORTED`) tests in
> `test/test_bc_ic_goldens.jl` correctly pin current behavior.

The Robin rule is **spec'd and schema-valid** (the ESM round-trips losslessly and
the rule JSON parses). What is missing is **engine plumbing in ESS**, all three
pieces of which are *imperative engine changes* of the same class as the G8/G11
ESS beads (`ess-tox`, `ess-gj4` — "GENERIC engine support, the one code-not-fixture
exception"), **not** anything expressible as declarative JSON over the existing
vocabulary.

The bead title reads *"MAYOR-HELD on G8"* and the parent epic assumed *"Once the
matcher (G8) fires the rule, add … fixtures."* That assumption is **incomplete**:
the G8 matcher (`ess-tox`) is necessary but **not sufficient** for Robin. G8
delivered kind/side *discrimination* only; it did not deliver coefficient
transport, a general-ghost lowering, or `$h` threading — and its companion
"retire the imperative ghost shadow" sub-goal was **reverted** (see §2).

---

## 1. What the bead asked

> "Robin BC is schema + `@test_throws E_BC_UNSUPPORTED` only. Once the matcher
> (G8) fires the rule, add numeric UNIT (ghost coeffs) + INTEGRATION (MMS solve)
> fixtures. … DECLARATIVE-OR-FAIL: the rule MUST be declarative JSON
> (arrayop/stencil composed over the EXISTING ESS engine vocabulary). If it
> CANNOT be expressed declaratively, do NOT implement it — STOP and report the
> precise gap for escalation (NO imperative/hand-coded operator)."

Acceptance criteria (verbatim): *"Robin BC declarative rule fires + numeric UNIT
+ INTEGRATION fixtures green; OR declarative-infeasibility reported."* This
document is the second branch.

---

## 2. What G8 (`ess-tox`) actually delivered — and what it didn't

G8 = ESS bead `ess-tox` ("BC engine kind/side matcher"), **closed**. It delivered
**kind/side discrimination** for synthetic `bc` nodes:

- `parse.jl:157-168` exposes the authored `kind`/`side` keys on the parsed
  `OpExpr` as the generic match fields `fn`/`dim`.
- `canonicalize.jl:284-291` emits `fn` in the canonical form so
  `bc(u,dirichlet,xmin)` and `bc(u,robin,xmin)` no longer collapse to the same
  canonical JSON (the matcher is now kind-aware).

That is the *whole* of what landed for matching. Crucially, `ess-tox`'s
acceptance also promised *"imperative ghost shadow retired/superseded"* — and
**that half was reverted**:

| ESS bead | Status | Effect on the production BC path |
|---|---|---|
| `ess-e7u` | closed (`609bceae`) | **Retired** `_apply_nonperiodic_bcs!`, `_check_bc_supported`, `_instantiate_bc_cell`, `_substitute_ghost`, … on the plan that makearray-region lowering (`ess-j8t`) + a `bind-guards-interface` would replace them. *"That `bind-guards-interface` was never built."* → ESD trunk went **50-red**. |
| `ess-1nm` | closed | **Restored** `_apply_periodic_folding!` / `_apply_nonperiodic_bcs!` / `_instantiate_bc_cell` / `_check_bc_supported` to fix the regression. |

Net state on ESS `main` today: the **imperative** `_apply_nonperiodic_bcs!` is
the production ghost-substitution path, and the makearray-region path
(`_substitute_ghost` / `_materialize_regions!`) is **absent** (`git grep` on
`origin/main` → no matches). The declarative BC rules do **not** drive production
lowering.

---

## 3. The three precise gaps (why declarative-infeasible)

### Gap 1 — No coefficient transport (Robin-specific)

A Robin BC carries **three** independent scalars (`α`, `β`, `γ`). ESS has nowhere
to put them:

- `types.jl:97-165` — `OpExpr` is a struct with a **fixed field set**
  (`op, args, wrt, dim, …, fn, name, value, table, …`). There is **no** field for
  `robin_alpha/beta/gamma`, and no generic "extra attributes" map.
- `parse.jl` reads only those known keys for a `bc` node (`kind`→`fn`,
  `side`→`dim`, `value`→`value`). `robin_alpha/beta/gamma` are **silently
  dropped**. Empirical proof (§4): parsing
  `{op:bc, kind:robin, side:xmin, args:[u], robin_alpha:1, robin_beta:1, robin_gamma:1.5}`
  canonicalizes to `{"args":["u"],"dim":"xmin","fn":"robin","op":"bc"}` — the
  coefficients are gone.
- The production BC-rewrite entry point `_discretize_bc!` (`discretize.jl:881-942`)
  builds the synthetic wrapper from only `variable/kind/side/value`
  (`discretize.jl:898-913`): `{op:bc, fn:kind, dim:side, args:[variable, value?]}`.
  Robin has no `value`, so the wrapper is `{op:bc, fn:robin, dim:xmin, args:[u]}`
  — it **never forwards** the coefficients.

So a declarative pattern (`robin_alpha:"$a"`, …) has **nothing to bind** — the
fields don't survive parsing and aren't on the wrapper. `0` occurrences of
`robin` exist in the entire ESS `src/` tree.

### Gap 2 — No general-ghost lowering in production (shared with G10)

Even if a rule could fire, the production ghost substitution is imperative and
kind-restricted:

- `_apply_nonperiodic_bcs!` (`discretize.jl:767`, called from `discretize.jl:355`)
  and its helper `_instantiate_bc_cell` (`discretize.jl:706-709`) hard-code exactly
  two ghost policies: **dirichlet** → ghost = the BC `value` (`_deep_native(value)`)
  and **zero-flux neumann** → mirror the in-range cell (`e = 1 - e` / `2N + 1 - e`).
- `_check_bc_supported` (`discretize.jl:746-758`) **throws `E_BC_UNSUPPORTED`** for
  every other kind, including Robin and *nonzero* Neumann.

There is **no** mechanism to apply a rule-produced or coefficient-parameterized
ghost expression. The Robin ghost is
`u_ghost = (2·h·γ + (2·β − α·h)·u₀) / (α·h + 2·β)` — a function of `u₀`, `h`, `α`,
`β`, `γ`. The existing lift can substitute a constant (dirichlet) or a mirror
(zero-neumann), not this.

### Gap 3 — `$h` not threaded for the BC-rewrite path (shared with G10)

The Robin/Neumann ghost needs the grid spacing `$h`. The guard machinery exists —
`bind_side_spacing` (`rule_engine.jl:629/779`, from ESS bead `ess-bps`) binds
`$h = 1/N` from `ctx.grids[grid].{spatial_dims, dim_sizes}` — but:

- **No ESD BC rule declares a `guards` block at all** (`grep -L '"guards"'` over
  `discretizations/**/*.json`). `robin_bc.json` / `neumann_bc.json` mention
  `bind_side_spacing` only in *prose*; `$h` is referenced in the replacement but
  never bound.
- The production `_discretize_bc!` wrapper-rewrite path does not thread grid
  context into the rewrite of the synthetic `bc` node, so even a correctly-authored
  guard could not resolve `grid`/`$h` there.

---

## 4. Empirical evidence (reproducible)

Run from the repo root after `Pkg.instantiate()`:

```julia
import EarthSciDiscretizations
using EarthSciDiscretizations: build_ode_problem
using EarthSciSerialization: RuleEngineError, parse_expression, canonical_json
repo = dirname(dirname(pathof(EarthSciDiscretizations)))

# INTEGRATION path — build_ode_problem on the Robin (and nonzero-Neumann) fixture:
for name in ("bc_robin_1d", "bc_neumann_1d_nonzero")
    esm = joinpath(repo, "test", "fixtures", "bc_ic", "$name.esm")
    gdd = joinpath(repo, "discretizations", "gdd", "cartesian_1d_nonperiodic_n8.gdd.json")
    try build_ode_problem(esm; grid_ref=gdd) catch e
        println(name, " → ", e.code, ": ", e.message) end
end

# Coefficient transport — does parse_expression keep robin_alpha/beta/gamma?
node = Dict("op"=>"bc","kind"=>"robin","side"=>"xmin","args"=>["u"],
            "robin_alpha"=>1.0,"robin_beta"=>1.0,"robin_gamma"=>1.5)
println(canonical_json(parse_expression(node)))
```

Output (observed):

```
bc_robin_1d → E_BC_UNSUPPORTED: boundary condition kind 'robin' on 'u' along 'x'
    is not supported by the arrayop lift (supported: dirichlet, zero-flux neumann)
bc_neumann_1d_nonzero → E_BC_UNSUPPORTED: nonzero-Neumann boundary condition on
    'u' along 'x' is not yet supported by the arrayop lift (requires
    grid-spacing-aware ghost extrapolation)
{"args":["u"],"dim":"xmin","fn":"robin","op":"bc"}   ← α/β/γ dropped
```

---

## 5. Why no alternative declarative encoding works

- **Coefficients-in-`args`** (`args:[$u,$a,$b,$g]` instead of named fields): a rule
  could match such a node *in isolation*, but `_discretize_bc!` never populates
  those args from the ESM BC dict (Gap 1), so it cannot fire in the production
  pipeline. (It also cannot be exercised by the walker's Layer-A rewrite harness:
  `test/walk_esd_tests.jl::_expr_from_json` lifts only `op/args/wrt/dim`, dropping
  `kind`/`fn`, so the kind-discriminating pattern would not match there either.)
  Authoring such a rule would validate a fiction disconnected from production.
- **Region/algebraic-equation encoding** (`_discretize_equation_bc!`, the `region`
  path): a boundary region equation defines a *separate* equation; it does not
  substitute the ghost into the interior Laplacian's boundary stencil. The lowering
  that *could* have unified them — makearray-region (`_substitute_ghost`, `ess-j8t`)
  — is **retired** (§2). So the interior stencil at `i=1` still resolves its ghost
  read through the imperative lift (Gap 2), which throws.

There is no spatial-AST encoding of "substitute this `u₀`/`h`/coefficient-dependent
ghost into the boundary stencil read" that the *existing* engine evaluates.

---

## 6. Target shape once ESS gains support (for the future implementer)

The declarative ESD rule that *should* fire — pending the ESS work in §7 — is a
kind/side-discriminated `bc` rewrite with a spacing guard:

```jsonc
// robin_bc.json (target; requires ESS Gaps 1–3 closed)
{
  "rules": { "robin_bc": {
    "pattern":  { "op": "bc", "kind": "robin", "side": "$side",
                  "args": ["$u", "$a", "$b", "$g"] },     // α,β,γ transported as args
    "guards":   [ { "name": "bind_side_spacing", "pvar": "$h",
                    "side": "$side", "grid": "$grid" },
                  { "name": "var_has_grid", "var": "$u", "pvar": "$grid" } ],
    "replacement": { /* u_ghost = (2hγ + (2β − αh)u₀)/(αh + 2β), as today */ }
  } }
}
```

The exact `replacement`/`produces` AST already in the committed `robin_bc.json` is
analytically correct (Dirichlet limit `b=0,a=1` → `2u_bc−u₀`; Neumann limit
`a=0,b=1` → `u₀+h·g`); only the *binding surface* (pattern args + guards) and the
ESS plumbing are missing.

---

## 7. Required ESS engine work (front-door change) — escalation

Filed as cross-rig ESS bead **`ess-lhi`** ("ESS: Robin + nonzero-Neumann BC
coefficient transport + general-ghost lowering (unblocks G9/G10)"). Three
changes, all "GENERIC engine support":

1. **Transport Robin/Neumann coefficients** from the ESM BC dict through
   `_discretize_bc!` onto the synthetic `bc` wrapper (e.g. append `α,β,γ` to
   `args`, or add typed `OpExpr` fields), so a declarative pattern can bind them.
2. **General-ghost lowering**: teach the non-periodic BC lift
   (`_apply_nonperiodic_bcs!` / `_instantiate_bc_cell`) to apply a
   rule-produced / coefficient-parameterized ghost expression (function of `u₀`,
   `h`, coefficients) instead of throwing `E_BC_UNSUPPORTED` — OR revive the
   makearray-region `_substitute_ghost` path that `ess-e7u` planned and `ess-1nm`'s
   restore note flags as the intended successor.
3. **Thread `$h` / grid context** into the `bc`-wrapper rewrite so a
   `bind_side_spacing` guard resolves.

This is **distinct** from the in-progress `ess-gj4` (generic *const_array*
boundary policy for the G11 covariant-Laplacian metric gathers): that concerns
periodic-wrap/edge-extend for CONST gathers, not coefficient-bearing *state-var*
ghost extrapolation.

---

## 8. Disposition

- **Verdict:** DECLARATIVE — INFEASIBLE over the existing ESS engine (§0).
- **Action taken:** report only. **No firing rule, no UNIT/INTEGRATION fixtures,
  no imperative executor** committed — per `DECLARATIVE-OR-FAIL`. The committed
  `robin_bc.json` remains a documented non-firing placeholder; its rewrite fixture
  keeps `applicable:false` (skip_reason refreshed to the accurate gap); the
  `@test_throws E_BC_UNSUPPORTED` goldens stay (they pin current behavior).
- **Escalation:** cross-rig ESS engine bead **`ess-lhi`** (§7) blocking esd-6g4.8
  (G9) **and** esd-6g4.9 (G10); strategic disposition routed to the mayor via the
  epic's coverage-drain bead `esd-6g4.15` (G16).
