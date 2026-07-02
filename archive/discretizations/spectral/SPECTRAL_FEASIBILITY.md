# Spectral family — declarative feasibility verdict

| | |
|---|---|
| **Bead** | esd-6g4.4 — G4: spectral family rules (Fourier/Chebyshev/spherical-harmonic) — DECLARATIVE-OR-FAIL |
| **Epic** | esd-6g4 — ESD operator×grid coverage (DECLARATIVE-OR-FAIL) |
| **Label** | `operator-grid-coverage` |
| **Outcome** | Verdict (no imperative solver added; no non-general rule shipped) |

**Oracle source pin.** The findings below were derived against the exact
EarthSciSerialization (ESS) source that this repo's `Project.toml` `[sources]`
(`rev = "main"`) resolves to via the instantiated `Manifest.toml`:

- ESS commit **`60c1ca89423a95418d73517d01e353d5b1551715`** ("fix(schema): sync
  language package copies after ess-z31 OQ1 edit", 2026-06-14)
- subdir-tree (`packages/EarthSciSerialization.jl`) **`8542c29a1534c7a76751e0dd35baa668eb7eb949`**
  (the `git-tree-sha1` pinned in `Manifest.toml`)
- ESS `origin/main` HEAD at verdict time: `f137bf749f926e03d506cb7e65843c4b5f3e7e32`
  (has advanced past the pinned tree; re-verify if the Manifest is re-resolved —
  see [[ess-feature-check-github-main]]).

`src/<file>.jl` citations are ESS files at the pinned tree unless prefixed
`ESD:`.

---

## 0. Verdict (TL;DR)

Spectral operators are **global** (every output point depends on every input
point). The only declarative shape that can express that over the existing ESS
engine is a **dense matrix contraction** — an `arrayop` with `reduce:"+"` whose
body sums a differentiation-matrix entry `D[i,j]` against `u[j]` over a
reduction index `j`. This was the hypothesis the bead named ("a dense transform
/ DFT-matrix contraction arrayop"). It was tested empirically. Results:

1. **The contraction mechanism works.** A dense `reduce` `arrayop` over a state
   vector executes end-to-end through the canonical pipeline
   (`discretize → ArrayOp → eval`) to machine precision. ✅
2. **Fourier and Chebyshev differentiation matrices are fully AST-expressible**
   — entries built inline from `cos`/`sin`/`ifelse`/`==` and index arithmetic,
   **no precomputed matrix, no imperative code**. Both reproduce the analytic
   derivative to **machine precision** (Fourier L∞ ≈ 3e-16; Chebyshev ≈ 4e-15).
   This **overturns the "likely infeasible" prior** for differentiation. ✅
3. **BUT no _general_ spectral rule is shippable today**, for two reasons:
   - **GAP A (differentiation blocker).** A spectral reduction must span the
     full grid (`j ∈ 1..N`). Reduction-range bounds must be **concrete
     compile-time integers**; they cannot reference the dynamic grid dimension
     size. `discretize` injects grid sizes into *output* ranges only — a rule's
     *reduction* range passes through unchanged — so a dense-reduction rule can
     only be written for **one hardcoded N**, not as a grid-agnostic rule.
   - **GAP B (convergence-harness mismatch).** Spectral accuracy is
     *exponential*, not an algebraic order; the Layer-B MMS runner measures
     algebraic order, so even a fixed-N spectral rule cannot be validated by the
     standard convergence harness.
4. **Spherical harmonics are fundamentally infeasible — GAP C.** The transform
   needs associated Legendre functions `Pₗᵐ(cosθ)`, which are **recurrence-
   defined** with no elementary closed form. The AST has **no recurrence/loop/
   fold**, so the SHT matrices cannot be assembled in-AST; supplying them
   precomputed is math-out-of-the-AST, which the declarative mandate forbids
   (ESD `AGENTS.md`).

**Net:** differentiation (Fourier/Chebyshev) is **declaratively expressible**
and proven correct, but blocked from being a *general, testable* rule by a
single, well-scoped ESS limitation (GAP A) plus a test-harness gap (GAP B);
spherical-harmonic operators are **declaratively infeasible** (GAP C). The
correct DECLARATIVE-OR-FAIL outcome is this report + the precise ESS engine
support that would unblock differentiation (§5). No imperative spectral solver
was added.

Reproducible evidence: [`feasibility_probe/`](feasibility_probe/) (runs against
an ESS-only environment; see its README).

---

## 1. What works (positive findings)

### 1.1 The engine vocabulary is pointwise + an arrayop-level reduce

The scalar op table is pointwise only — `+ - * / ^`, comparisons, `and/or/not`,
`ifelse`, `sin cos tan asin acos atan exp log sqrt abs sign floor ceil min max`,
`pi`/`e` (`src/tree_walk.jl:625-746`). There is **no scalar reduction op**.
Reduction lives at the `arrayop` level via two fields parsed onto every
`arrayop`: `reduce` and `ranges` (`src/parse.jl:61,128,180`; round-tripped in
`test/array_ops_test.jl:344-360`). Semantics: an index that appears in the body
but **not** in `output_idx` is reduced over with the `reduce` op — i.e. a
contraction. ESS's own test states the matrix-multiply shape
`C[i,j]=Σ_k A[i,k]·B[k,j]` (`test/array_ops_test.jl:346-351`).

### 1.2 Precedent: a reduce-arrayop rule already ships in ESD

`ESD:discretizations/finite_difference/nn_diffusion_duo.json` is a `reduce:"+"`
arrayop rule whose body double-indexes a 2-D array by `[output_idx, reduce_idx]`
(`edges_on_face[i,k]`). A spectral operator has the **same shape**, except the
reduction spans the whole grid and `D[i,j]` is a differentiation matrix instead
of mesh geometry.

### 1.3 The dense contraction executes (proven)

Through the official tree-walk evaluator (`build_evaluator`, the runner the
ESD Layer-B per-topology runners use — `ESD:src/ode_problem.jl:84,100`,
`ESD:test/walk_esd_tests.jl:1228-1230`):

- `du[i] = Σ_j u[j]` (dense sum over a **state** vector) → correct.
- `du[i] = Σ_j j·u[j]` → correct (the reduction index is usable as a **scalar
  value** in body arithmetic, not only as a subscript).

### 1.4 Fourier & Chebyshev matrices, inline in the AST (proven, machine precision)

Periodic Fourier first-derivative matrix (N even, `xⱼ=2π(j-1)/N`):

```
D[i,j] = ifelse(i==j, 0, 0.5 · (-1)^(i-j) · cot((i-j)π/N))
```

expressed with `(-1)^m = cos(πm)` and `cot θ = cos θ / sin θ` — all existing
ops. Contracted `out[i]=Σ_j D[i,j]·u[j]`, it reproduces `d/dx sin(x)=cos(x)`
to **L∞ ≈ 3.3e-16** (N=8) and `d/dx sin(2x)` to **≈ 4.4e-15** (N=16). The
diagonal singularity is harmless: `ifelse` selects the `0` branch and the
unused `cot(0)` value is discarded (`src/tree_walk.jl` ifelse arm ~`694-700`).

Chebyshev–Gauss–Lobatto first-derivative matrix (`xₖ=cos((k-1)π/(N-1))`),
including endpoint `cᵢ` factors, the `±(2(N-1)²+1)/6` corners and the
`-xᵢ/(2(1-xᵢ²))` interior diagonal — all via `cos`/`ifelse`/`==`/rational ops —
differentiates `x³` exactly: **L∞ ≈ 3.6e-15** (N=6). No precomputed matrix.

The full Fourier rule was also driven through the **real rule→discretize
pipeline** (`discretize(esm; lift_1d_arrayop=true)` then `build_evaluator`):
the `grad` pattern is rewritten to the reduce-arrayop and evaluates to
**L∞ ≈ 3.3e-16** for a concrete grid size.

---

## 2. GAP A — reduction-range bounds cannot track the grid size (differentiation blocker)

A general rule must reduce over `j ∈ 1..N` where **N = grid dimension size**,
unknown at authoring time. Empirically:

- `discretize` resolves **output** ranges to concrete integers from the grid
  (e.g. the shipped `laplacian_2nd_uniform_cartesian` canonical golden has
  `ranges:{i:[1,4],j:[1,4]}` from `size:4`; output index←dimension mapping at
  `src/rule_engine.jl:1246-1248`).
- A rule's **reduction** range is passed through **unchanged**: declaring
  `ranges:{j:[1,"Nx"]}` (parameter) or `ranges:{j:[1,"x"]}` (dimension name)
  survives `discretize` verbatim — neither is resolved to the grid size.
- At evaluation the unroller requires concrete integers:
  `_expand_int_range(r::Vector{Int})` (`src/tree_walk.jl:907`); a symbolic
  bound reaches `_eval_const_int(::VarExpr)` and throws
  **`E_TREEWALK_UNBOUND_LOOP_VAR`** (`src/tree_walk.jl:924`). Observed:
  `E_TREEWALK_UNBOUND_LOOP_VAR: Nx` and `: x`.
- A parameter cannot help: the loop is unrolled at build time, before
  parameters are bound, so a runtime `p.Nx` is not available to size the loop.

A *dimension-keyed* range **is** resolved to a size — but only inside the
`scheme` / `multi_output_stencil` expander (`src/scheme_expansion.jl:750-763`,
`:1015-1028`, reading `dim_sizes`), which handles only the `stencil` /
`multi_output_stencil` kinds (see [[esd-declarative-or-fail-scheme-kinds]]).
A dense contraction is neither kind, so it never reaches that resolver. There is
**no token** that sizes a *reduction* range in an `arrayop` replacement.

**Consequence:** a declarative spectral differentiation rule can be written and
proven correct **only for a single hardcoded N**. That is not a usable grid-
agnostic discretization rule, so none is shipped (per the mandate: do not ship
something that cannot be expressed *as a rule*).

---

## 3. GAP B — convergence harness assumes algebraic order

Layer-B (`run_mms_convergence`, `ESD:test/walk_esd_tests.jl`) computes
`log2(errₙ/err₂ₙ)` and asserts a minimum **algebraic** order. Spectral methods
converge **exponentially** for smooth data — error collapses to machine
precision at modest N, so the algebraic-order statistic is meaningless (it is
either undefined or astronomically large, and round-off dominates the ratio).
Validating a spectral rule needs a **property-based** Layer-B variant
(e.g. "L∞ error < tol at moderate N", mirroring the limiter property runner
`run_layer_limiter`, `ESD:test/walk_esd_tests.jl:183`), which does not exist.
This is a *testing* gap, secondary to GAP A, but it independently blocks a
convergence-tested spectral rule.

---

## 4. GAP C — spherical harmonics are fundamentally infeasible

A spherical-harmonic operator (e.g. the Laplace–Beltrami operator, diagonal
`-l(l+1)/R²` in spectral space) is applied as
`inverse-SHT ∘ diag(eigenvalue) ∘ forward-SHT`. The transforms require the
**associated Legendre functions `Pₗᵐ(cosθ)`** evaluated at the grid latitudes.

- `Pₗᵐ` has **no elementary closed form** in `(l, m, θ)`; it is defined by a
  three-term recurrence in `l` (equivalently an `l`-th derivative via Rodrigues).
- The ESS AST has **no recurrence, no loop, no fold over a variable count, and
  no special-function op** for `Pₗᵐ` (op table, `src/tree_walk.jl:625-746`).
  A fixed `arrayop`/`reduce` cannot synthesize an `l`-step recurrence.
- The composite grid-point operator `L[i,j] = Σ_{l,m} Pₗᵐ(θᵢ)Pₗᵐ(θⱼ)·(…)` has
  no closed form either — it is exactly the sum over `Pₗᵐ` that cannot be
  written in-AST.
- Supplying `Pₗᵐ` (or `L`) as a precomputed `const` matrix moves the
  discretization math out of the AST, which the declarative mandate forbids
  ("Do NOT write … helper functions that implement math a rule should own",
  ESD `AGENTS.md`). (Such an inline `const` matrix is not even indexable in the
  evaluator backends — only state arrays and runner-supplied named arrays are.)

Even the Gauss–Legendre grid **nodes** are roots of Legendre polynomials (no
closed form); those could be admitted as grid coordinate primitives, but the
`Pₗᵐ` **basis values inside the transform** cannot. Spherical harmonics are
therefore declaratively infeasible — not for lack of a contraction mechanism,
but because the basis itself is not AST-expressible.

---

## 5. Recommended generic engine support (to unblock differentiation)

GAP A is small and well-scoped, and it lives in **ESS** (the engine), not ESD
(which carries no evaluator — ESD `AGENTS.md` single-pathway rule). The ESS-vs-
ESD ownership of any sweep/range machinery is a **mayor** decision
([[esd-declarative-or-fail-scheme-kinds]]). Recommended ESS enhancement:

> Allow an `arrayop` **reduction**-range bound to reference a grid dimension
> (by dimension name, or a `{dim_size: <name>}` token), resolved to the concrete
> size at `discretize` time — the same resolution already applied to *output*
> ranges and to `multi_output_stencil` ranges (`src/scheme_expansion.jl:756`).

With that one change, the proven Fourier/Chebyshev rules in the appendix become
**general** declarative rules. GAP B (a property-based spectral Layer-B runner)
would then make them convergence-testable. Both are follow-ups for escalation
(this verdict feeds the G-series drain bead esd-6g4.15); no imperative operator
should be added to force a pass in the meantime.

---

## Appendix — proven rule ASTs (proof-of-concept, blocked on GAP A)

These are correct and execute to machine precision at a fixed grid size `N`.
They are **not** committed as catalog rules because the reduction bound `N`
cannot yet track the grid (GAP A). `Nx` denotes the grid-size parameter; the
literal `N` in a bound marks where GAP A forces a hardcode.

**Fourier first derivative** (`grad`, cartesian, periodic):

```jsonc
{ "applies_to": { "op": "grad", "args": ["$u"], "dim": "x" },
  "grid_family": "cartesian",
  "replacement": {
    "op": "arrayop", "reduce": "+",
    "ranges": { "j": [1, N] },            // GAP A: must be the grid size
    "expr": { "op": "*", "args": [
      { "op": "ifelse", "args": [
        { "op": "==", "args": ["i", "j"] }, 0.0,
        { "op": "*", "args": [ 0.5,
          { "op": "cos", "args": [ { "op": "*", "args": [
            { "op": "pi", "args": [] }, { "op": "-", "args": ["i","j"] } ] } ] },
          { "op": "/", "args": [
            { "op": "cos", "args": [ { "op": "*", "args": [
              { "op": "-", "args": ["i","j"] },
              { "op": "/", "args": [ { "op": "pi", "args": [] }, "Nx" ] } ] } ] },
            { "op": "sin", "args": [ { "op": "*", "args": [
              { "op": "-", "args": ["i","j"] },
              { "op": "/", "args": [ { "op": "pi", "args": [] }, "Nx" ] } ] } ] } ] } ] } ] },
      { "op": "index", "args": ["$u", "j"] } ] } } }
```

**Chebyshev first derivative** (Chebyshev–Gauss–Lobatto, `xₖ=cos((k-1)π/(N-1))`):
entries `D[i,j] = (cᵢ/cⱼ)(-1)^(i+j)/(xᵢ-xⱼ)` off-diagonal,
`-xᵢ/(2(1-xᵢ²))` interior diagonal, `±(2(N-1)²+1)/6` corners, with
`cᵢ = ifelse(i==1 ∨ i==N, 2, 1)`, `(-1)^(i+j)=cos(π(i+j))`,
`xᵢ=cos((i-1)π/(N-1))` — all existing ops; same `arrayop`/`reduce` envelope as
above. Full builder: [`feasibility_probe/probe3.jl`](feasibility_probe/probe3.jl).
