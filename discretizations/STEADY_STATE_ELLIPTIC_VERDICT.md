# Steady-state / elliptic operator — DECLARATIVE-OR-FAIL verdict

**Bead:** esd-6g4.5 (G5) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Oracle source pin:** EarthSciSerialization GitHub `main` @
`fd2e588318b0a8343540d4ed31a6532d99301ab0` — the rev ESD `Project.toml`
`[sources]` resolves (verified equal to the local refinery ESS checkout HEAD).
All ESS line numbers below are at that commit. ESD line numbers are at this
branch's base (`origin/main`).

This file is a verdict doc only. It lives at the `discretizations/` root, so it
is invisible to `load_rules` (which globs `*.json` inside family subdirs) — it
adds no rule, no fixture, and no walker/conformance surface.

---

## 0. Verdict (TL;DR)

> **DECLARATIVE — FEASIBLE, and already expressible with ZERO new rules.** The
> discrete elliptic **operator** and **residual** ARE expressible over the
> EXISTING ESS engine. They require **no new rule, no new fixture, no engine
> code, and no imperative solver.** The discrete elliptic operator is the
> existing `laplacian_*` / `nonlinear_laplacian_*` family — those rules are
> **solve-mode-agnostic**, and the *same* rule produces a steady elliptic
> residual when the operator is placed on the **RHS of an algebraic
> (non-time-derivative) equation**. ESS then classifies the system as a DAE
> (`system_class: "dae"`) and hands the *solve* to an external DAE/nonlinear
> assembler.
>
> The **only** out-of-scope piece is the elliptic **solve** (inverting the
> discrete operator for `u`). That is explicitly an ESS RFC §2 non-goal
> ("implicit solvers … out of scope") and is barred by the bead ("Do NOT add an
> imperative elliptic solver") and the RHS-only library principle. No "solver
> rule type" exists; the roadmap already tracks it as `not_in_repo`.

**Net:** there is nothing to author. "Elliptic" is not a missing *rule* — it is
a *solve mode* that the existing operator catalog already feeds, demonstrated
empirically in §3 and pinned by a new regression test
(`test/test_elliptic_residual_discretization.jl`).

---

## 1. The bead's two questions, disentangled

The bead asks for "the discrete elliptic **operator/residual** declaratively,
NOT a solver," and flags a tension with the RHS-only library principle. Three
distinct objects hide inside "elliptic," and they live at different layers:

| Object | Math | Layer | Status |
|---|---|---|---|
| Elliptic **operator** | `L[u] = ∇²u` (or `∇·(a∇u)`) | discretization **rule** | **Already cataloged** (`laplacian_2nd_uniform_cartesian.json`, `nonlinear_laplacian_*`, `spherical_laplacian_*`, `covariant_laplacian_cubed_sphere`) |
| Elliptic **residual** | `R[u] = ∇²u − f = 0` | **equation** (model) | **Expressible for free** via existing rule + algebraic equation (§3) |
| Elliptic **solve** | find `u` s.t. `R[u]=0` | external integrator/solver | **Out of scope** (RFC §2; §6) |

The discretization-rule layer owns only the first. The second is an
*equation-structure* choice the model author makes; the third is downstream of
ESS entirely.

---

## 2. The operator is already cataloged — and solve-mode-agnostic

A discretization rule rewrites *one operator subtree* into a stencil; it never
assembles `L u = f`. The existing 2nd-order Cartesian Laplacian
(`discretizations/finite_difference/laplacian_2nd_uniform_cartesian.json`) is the
canonical discrete elliptic operator. Critically, **the rule does not know or
care whether it is used in a parabolic or an elliptic problem** — the solve mode
is decided entirely by the *equation* the operator sits in:

- `∂u/∂t = ∇²u` (LHS is `D(u, wrt=t)`) → **parabolic**, `system_class: "ode"`.
- `f = ∇²u` (LHS is a plain variable) → **elliptic/steady**, `system_class: "dae"`.

Both rewrite through the *identical* Laplacian rule. This is the precise content
of the "RHS-only tension": the operator must appear on the **RHS** of the
equation to be discretized (§4), which is exactly the natural residual form
`f = ∇²u` ≡ `∇²u − f = 0`.

---

## 3. Empirical demonstration (attempt-in-earnest)

Built from the Laplacian rule's own canonical fixture
(`…/laplacian_2nd_uniform_cartesian/fixtures/canonical/input.esm`, a 4×4 periodic
Cartesian grid), swapping the ODE equation for a steady Poisson equation and
adding a forcing field `f`. Driven through `EarthSciSerialization.discretize`
(the canonical pipeline; ESD's `eval_coeff` is a passthrough). Reproduced by the
committed regression test.

| Case | Equation given to `discretize` | `system_class` | `algebraic_eqn_count` | Laplacian rewritten → stencil? |
|---|---|---|---|---|
| Baseline (ODE) | `∂u/∂t = ∇²u` (`lhs: D(u,wrt=t)`) | `ode` | 0 | **yes** |
| **Form A** (operator on RHS) | `f = ∇²u` (`lhs: "f"`) | **`dae`** | **1** | **yes** ✅ |
| Form B (operator on LHS) | `∇²u = f` (`lhs: laplacian(u)`) | `dae` | 1 | **no** ❌ |

Plus the DAE contract: **Form A under `dae_support=false` throws
`E_NO_DAE_SUPPORT`** — correct, a binding that hands off only to an ODE
integrator cannot accept the elliptic (algebraic) system.

**Form A is the feasible path.** The discrete output equation is, verbatim:

```jsonc
{ "lhs": "f",
  "rhs": { "op": "+", "args": [
      // (1/dx²)·u[i-1,j] + (-2/dx²)·u[i,j] + (1/dx²)·u[i+1,j]
      // + (1/dy²)·u[i,j-1] + (-2/dy²)·u[i,j] + (1/dy²)·u[i,j+1]
      …  // the 5-point discrete Laplacian, periodic-folded
  ] } }
```

i.e. the discrete elliptic **residual** `f = ∇²ₕu` (≡ `∇²ₕu − f = 0`), produced
by the existing rule, classified algebraic, **with the solve left external.**
Form B confirms the RHS-only constraint: an operator on the equation LHS is *not*
rewritten (see §4).

---

## 4. Why it works — engine mechanism (file:line @ pin)

`EarthSciSerialization.jl/src/discretize.jl`:

- **`_discretize_equation!` (L864–877):** *"Rewrite RHS. The LHS is of the form
  `D(x, wrt=t)` or `x`; we canonicalize it"* — the rule engine rewrites **only
  the RHS** (`_rewrite_or_passthrough!`); the LHS is merely canonicalized
  (`_canonicalize_value`). This is exactly why Form A (operator on RHS) is
  discretized and Form B (operator on LHS) is not. The RHS-only rewrite is by
  design, and it *is* the RHS-only principle in code.
- **`_apply_dae_contract!` (L1159; docstring from L1141):** an equation is
  *differential* **iff** its LHS is `D(x, wrt=<independent var>)`; **every other
  equation is algebraic**. It stamps
  `metadata.system_class = total_algebraic > 0 ? "dae" : "ode"` (L1191) and
  `metadata.dae_info`, and throws `E_NO_DAE_SUPPORT` (L1199) when an algebraic
  equation is present but `dae_support=false`. Its docstring pins the contract to
  RFC §12: *"hand to a DAE assembler, or abort."* (The `_is_algebraic_equation`
  predicate is at L1228.)
- `discretize(…; dae_support=…)` signature + DAE doc: L33, L67–77.

The DAE/algebraic path is already covered by ESD conformance
(`tests/conformance/discretization/dae_missing/{pure_ode_baseline,mixed_dae_observed}.json`
+ `test/test_dae_discretization_conformance.jl`). This verdict adds the missing
combination: a **spatial operator rule firing inside an algebraic equation**.

---

## 5. The RHS-only tension, resolved

The library exposes the discrete elliptic operator/residual **declaratively and
for free**: the operator rule is solve-mode-agnostic, and the residual emerges
when the model author writes the operator on the RHS of an algebraic equation.
The "RHS-only" requirement is not a limitation to work around — it is the very
mechanism that yields the residual form `∇²ₕu − f = 0`. ESS classifies the result
as a DAE and defers the solve. **No rule encodes "elliptic"; no rule encodes a
solve mode.** That separation is correct and is what makes the operator catalog
reusable across parabolic and elliptic problems.

---

## 6. What is out of scope — the solve

Authoring an elliptic *solver* (or a rule that bundles one) is barred and
unnecessary:

- **ESS RFC §2 (Non-goals, discretization.md:36–41):** "Time integration. …
  Time-stepping schemes, **implicit solvers**, operator-splitting … are out of
  scope." The solve is handed to MTK / SUNDIALS / nonlinear solvers downstream.
- **Bead instruction:** "Do NOT add an imperative elliptic solver."
- **No solver rule type exists.** `docs/rule-catalog-roadmap.md:302–303` lists
  `multigrid_geometric` ("Standard elliptic solver") and `gmres_krylov` as
  `family: solver`, `not_in_repo`, gap **"needs solver rule type"** — i.e. an ESS
  RFC schema extension, not a rule-authoring task.
- **The residual side-channel (`produces: {kind: algebraic}`, RFC §9.4,
  discretization.md:2128–2135) is *not needed* here and is not yet implemented**
  anyway: `discretize.jl`'s DAE docstring says such equations are handled "once
  that [rule] lands," and `discretizations/finite_difference/dirichlet_bc.json:4`
  notes `produces` is "not yet processed by the current rule engine; documented
  for forward compatibility." The elliptic residual does **not** rely on it —
  authoring the equation directly as algebraic (Form A) already works today.

A standalone "elliptic operator" rule would either (a) duplicate the existing
Laplacian, or (b) misuse the unimplemented `produces: algebraic` to bake a
solve-mode decision into the operator catalog — a category error. Neither is
warranted.

---

## 7. Scope, handoff, and reproduction

- **No rule or fixture authored** (correct: the operator already exists; the
  residual is an equation-layer composition). No imperative code added.
- **Regression test added:** `test/test_elliptic_residual_discretization.jl`
  pins §3 — the existing Laplacian rule, applied to a steady Poisson model,
  yields `system_class: "dae"` with the operator rewritten, errors
  `E_NO_DAE_SUPPORT` without DAE support, and stays `ode` in the time-dependent
  form (solve-mode-agnosticism). Runs under `Pkg.test` (TestItemRunner).
- **Open decision for the mayor / ESS spec owners (not a polecat task):** whether
  ESD/ESS should ever grow a first-class *solver* or *elliptic-solve* rule type.
  That is a spec change (RFC §7 schema + engine), explicitly out of scope for
  this bead and the RHS-only library. The mayor audit already flags this as a
  spec/MOL-gap (`esd-coverage-gaps-2026-06-21.md`: "steady-state/elliptic ✗ —
  absent from spec+ESD (also a MOL-gap)").
- Collected by the G16 drain bead (esd-6g4.15) as **REPORTED — feasible, no new
  rule required**.

---

## 8. Summary mapping — every elliptic ingredient → existing mechanism

| Elliptic ingredient | Declarative expression | Existing mechanism | New code? |
|---|---|---|---|
| Discrete operator `∇²ₕu` | `laplacian($u)` rewrite → 5-pt stencil | `laplacian_2nd_uniform_cartesian.json` + rule engine | **No** |
| Variable-coeff `∇·(a∇u)` | `nonlinear_laplacian_*` | existing rules | **No** |
| Residual `∇²ₕu − f = 0` | algebraic equation `f = ∇²u` (operator on RHS) | `_discretize_equation!` RHS rewrite (discretize.jl:864–877) | **No** |
| Steady/elliptic classification | LHS not `D(·,wrt=t)` ⇒ algebraic | `_apply_dae_contract!` → `system_class:"dae"` (discretize.jl:1159, stamp L1191) | **No** |
| Forcing field `f` | a `parameter` array variable in the model | base ESM variable | **No** |
| The solve (invert `L_h`) | — | external DAE/nonlinear solver | **Out of scope** (RFC §2) |

**Verdict: DECLARATIVE-FEASIBLE — already expressible; no rule, no fixture, no
solver to add.**
