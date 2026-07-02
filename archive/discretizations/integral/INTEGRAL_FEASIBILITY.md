# Integral / PIDE family — declarative feasibility verdict

| | |
|---|---|
| **Bead** | esd-6g4.14 — G15: integral/PIDE as a declarative catalog rule — DECLARATIVE-OR-FAIL |
| **Epic** | esd-6g4 — ESD operator×grid coverage (DECLARATIVE-OR-FAIL) |
| **Label** | `operator-grid-coverage` |
| **Outcome** | **FEASIBLE** — rule shipped (`midpoint_1d.json`) + Layer-A canonical fixture (byte-exact). No imperative operator added. |

**Oracle source pin.** Derived against the exact EarthSciSerialization (ESS)
source this repo's `Project.toml` `[sources]` (`rev = "main"`) resolves to via a
fresh `Pkg.instantiate()`:

- ESS commit **`fd2e588318b0a8343540d4ed31a6532d99301ab0`** ("refactor(schema):
  retire cross_metric scheme kind (RFC §7.6) …", 2026-06-21)
- subdir-tree (`packages/EarthSciSerialization.jl`) **`f2b411a6617fdb85539f2abddec7bb4523149d19`**
  (the `git-tree-sha1` pinned in the resolved `Manifest.toml`)

This pin is **newer** than the spectral verdict's
([[esd-dense-contraction-reduce-arrayop-gap]]; `SPECTRAL_FEASIBILITY.md` tree
`8542c29a…`), so the feasibility below was established **empirically against the
current engine** (`feasibility_probe/`). `src/<file>.jl` citations are ESS files
at this pinned tree unless prefixed `ESD:`.

---

## 0. Verdict (TL;DR)

**FEASIBLE.** A domain integral `∫u dV` is a **weighted sum-reduction over grid
cells** — `Σ_j w_j u_j` — which is exactly a `reduce:"+"` `arrayop` (the
"quadrature = weighted sum_product arrayop" the bead named). It ships as the
declarative catalog rule [`midpoint_1d.json`](midpoint_1d.json), proven through
the official `discretize → ArrayOp → eval` pipeline (`feasibility_probe/`, all
`err = 0.0`):

1. The sum_product mechanism executes and reproduces the midpoint integral to
   machine precision (Check 2). ✅
2. **The naive grid-spanning reduction bound fails** — tying `j ∈ 1..N` to the
   grid via a parameter `"Nx"` or dimension name `"x"` throws
   `E_TREEWALK_UNBOUND_LOOP_VAR` (Check 3). This is **GAP A** — the limitation
   the spectral verdict declared a hard blocker.
3. **★ GAP A is bypassed by supplying the size as DATA. ★** A reduction bound
   `index(size_x, 1)`, with `size_x = [N]` a host-supplied **const_array**,
   resolves at build-time unrolling and makes the reduction **track the grid
   size grid-agnostically** — `err = 0.0` at N = 8/16/32 (Check 4). This is the
   **same host-supplied-grid-data pattern `nn_diffusion` already uses** for mesh
   connectivity (its per-cell bound `index(n_edges_on_cell, i)`). It needs **no
   ESS change**: const_arrays are available during the unroll, whereas
   parameters bind later (§2).
4. The full **rule path** — the shipped `midpoint_1d` pattern/replacement
   rewriting the `integral` op → the sum_product arrayop via `discretize`,
   then evaluating — reproduces the integral exactly (Check 5). The Layer-A
   canonical fixture byte-matches the rewrite (`run_layer_a` → PASS, 1043 bytes).
   ✅

**No new rule or op was invented** — the shape is the existing `reduce` arrayop
(C. Tessum's bead note: "same sum_product rule, no special-casing … NO new
rule/op vs the structured integral"). 

**Cross-cutting correction.** The bypass is a property of the **arrayop
reduction-bound resolver** (`_resolve_scalar_arrayop` /
`_expand_int_range_dyn`, `src/tree_walk.jl:2078`), which is **operator-agnostic**
— it never inspects the arrayop body. So the **same `index(size, 1)`
const_array bound resolves identically for the spectral family's dense
contraction** `out[i] = Σ_j D[i,j] u_j`. GAP A is therefore **not** the hard
blocker the spectral verdict ([[esd-dense-contraction-reduce-arrayop-gap]])
concluded — the Fourier/Chebyshev rules it proved correct-but-blocked are
shippable by the same mechanism. Recommended re-confirmation feeds esd-6g4.15
(§5).

---

## 1. What works (proven, `feasibility_probe/probe.jl`)

### 1.1 The `integral` op is already grid-agnostic host quadrature

`build_ode_problem` carries **no** integral handling
(`ESD:src/ode_problem.jl` — it calls `discretize`/`build_evaluator`). The
quadrature lives in the ESS evaluator: the `integral` op
(`src/tree_walk.jl:2215-2243`) reads the integrand state variable's tracked
extent `array_var_info[vname] → (lo1, hi1)` and expands to
`d{x}·(u[lo1]+…+u[hi1])`, sized from the integrand's own grid bounds — correct
for **any N** (Check 1; the PIDE unit test `ESD:test/test_build_ode_problem.jl:502`).
This is engine-builtin host quadrature, not a catalog rule.

### 1.2 The scalar `reduce` arrayop reproduces the integral (the mechanism)

A scalar (`output_idx`-empty) `reduce:"+"` arrayop in equation-RHS position is
expanded inline by `_resolve_scalar_arrayop` (`src/tree_walk.jl:2062-2110`): the
contraction index is unrolled and the terms combined. With a **concrete** bound
it evaluates `du ≡ -Σ_j h·u_j = -1.0` to `err = 0.0` (Check 2).

### 1.3 The const_array-sized bound makes it grid-agnostic (the key)

Bound `index(size_x, 1)` with `size_x = [N]` a const_array → `err = 0.0` at
N = 8/16/32, **one grid-agnostic rule** (Check 4). Precedent for host-supplied
grid arrays: `ESD:discretizations/finite_difference/nn_diffusion_duo.json`.

### 1.4 The shipped rule rewrites and evaluates correctly

`midpoint_1d` (pattern `{op:integral, args:[$u], var:x, …}` → replacement
`dx · Σ_{j=1..index(size_x,1)} u[j]`) is applied by `discretize` and evaluates
to `err = 0.0` (Check 5). Its Layer-A canonical fixture byte-matches (§4).

---

## 2. Why the const_array bound bypasses GAP A (and the symbolic bound does not)

A quadrature reduction spans the whole grid (`j ∈ 1..N`, N = dimension size).
The arrayop contraction bound is expanded by `_expand_int_range` (const range)
or `_expand_int_range_dyn(…, idx_env, const_arrays)` (`src/tree_walk.jl:2078-2081`).
The bound must reduce to a concrete `Int`, **but it is evaluated against
`const_arrays` at unroll time**:

- A bare grid-size **symbol** (`"Nx"` param, `"x"` dim) reaches
  `_eval_const_int(::VarExpr)` with no binding and throws
  `E_TREEWALK_UNBOUND_LOOP_VAR` (`src/tree_walk.jl:1424-1427`) — the GAP A
  failure (Check 3). A **parameter** also fails: the loop unrolls at build time,
  before parameters bind.
- `index(size_x, 1)` is **not** a bare symbol — it is an `index` into a
  const_array, which `_eval_const_int` resolves to `const_arrays["size_x"][1] =
  N` **at build time** (const_arrays are threaded into the unroller). So the
  reduction is sized from grid DATA, not from a symbol the unroller can't see.
  This is exactly how `nn_diffusion`'s per-cell bound `index(n_edges_on_cell, i)`
  works (`src/tree_walk.jl:2070-2076`) — only here the count is a global grid
  size indexed by the literal `1`.

The resolver **never inspects the arrayop body**, so this is operator-agnostic:
the identical bound resolves a spectral differentiation contraction
`out[i] = Σ_j D[i,j] u_j` the same way. GAP A blocked only the *bound*, and the
bound is suppliable as data.

---

## 3. The shipped rule (`midpoint_1d.json`)

Scope: the **full-domain midpoint/Euler quadrature of a 1-D integral on a
uniform cartesian grid** — `∫u dx ≈ dx · Σ_{j=1}^{N} u_j`. The replacement is
`{op:*, args:["dx", {arrayop reduce:+ over j∈[1, index(size_x,1)] of index($u,j)}]}`.

**Matching caveat (documented limitation).** The ESS rule matcher
(`src/rule_engine.jl:286-318`, `_match`) keys on `op` + `args` + `wrt/dim/fn`
and **ignores** the integral op's `lower`/`upper`/`var` fields. So this rule
fires on **any** `integral($u)` whose first arg unifies, regardless of the
integration bounds — it treats every match as a full-domain integral. That is
correct for the domain integral / PIDE reduction term this bead targets; a
**partial**-integral rule (summing a sub-range) would need ESS to expose
`lower`/`upper` in matching, or a `var`-keyed selector. Noted for esd-6g4.15.

The weight is the cell width `dx` (a parameter, as in the PIDE GDD); the grid
size enters as the `size_x` const_array (§5).

---

## 4. Testing — Layer-A proves the rewrite; Layer-B MMS is N/A for a quadrature

- **Layer A (canonical):** `fixtures/canonical/` drives `discretize(input)` and
  byte-compares the canonical document to `expected.esm`. `run_layer_a` →
  **PASS** (1043-byte match): the `integral` op rewrites to the weighted
  sum_product arrayop. This is the declarative-rewrite acceptance signature, and
  it needs no eval / no const_arrays, so it **cannot regress the existing PIDE
  path** (a catalog rule applies only when a GDD references it via
  `ESD:_inject_rules!`; the PIDE test's GDD does not). Regenerate `expected.esm`
  the way the walker compares it:
  `canonical_doc_json(EarthSciSerialization.discretize(input))`
  (`ESD:test/walk_esd_tests.jl` `apply_rule_and_diff` / `canonical_doc_json`).
- **Layer B (MMS convergence): SKIP — semantically inapplicable.** The Layer-B
  runner measures the **algebraic order of a differential operator** on a
  manufactured solution. A quadrature is not a differential operator: midpoint
  is `O(h²)` in *quadrature* error, the property a dedicated property-test would
  assert (cf. the limiter runner `run_layer_limiter`), not the MMS-order runner.
  This mirrors the spectral verdict's GAP B (harness-shape mismatch). Eval
  correctness is instead proven by `feasibility_probe/` Check 5 (`err = 0.0`).

`load_rules` is unaffected structurally: the new `integral/` family contributes
exactly one rule (`midpoint_1d`); the probe/`*.md` files carry no `*.json` so
they are invisible to the catalog glob.

---

## 5. Production wiring & generalisations (→ esd-6g4.15 / G16)

- **Use the rule in real `build_ode_problem` PIDE solves.** The host must supply
  `size_<dim> = [N_dim]` const_arrays for cartesian grids (a few lines beside
  the existing `coord_<dim>` supply, analogous to `_unstructured_grid_const_arrays`).
  This is **GENERIC, harmless engine support**: `size_x` is referenced only when
  a GDD opts into `midpoint_1d`, so the ESS-builtin `integral` op path (§1.1) is
  untouched. Layer-A already proves the rewrite without it.
- **Nonuniform grids:** replace the scalar `dx` with a per-cell weight array
  `index(w, j)` (host-supplied cell widths) — same envelope.
- **Integral OPERATOR** `(Ku)_i = Σ_j w_j K(x_i,x_j) u_j` (field output): the
  **global/dense** kernel is the same `index(size_x,1)`-bounded contraction with
  an AST-expressible elementary kernel `K` (Gaussian, polynomial — built inline
  exactly as the spectral verdict built Fourier/Chebyshev matrices); the
  **compact-support** kernel is the `nn_diffusion` per-cell-bound shape
  (`feasibility_probe/` Check 6, proven), whose only remaining need is the
  spatial-join / value-invention broad-phase (`src/tree_walk.jl:739`,
  `materialize_value_invention`) to build the support connectivity from a radius.
- **Spectral re-confirmation (high value).** §0/§2: the same const_array bound
  bypasses the spectral GAP A. Re-running `discretizations/spectral/feasibility_probe/`
  with an `index(size,1)` bound should turn the Fourier/Chebyshev "proved but
  blocked" rules into shippable ones. Recommended for esd-6g4.15.

This bead's outcome for the drain (esd-6g4.15): the integral/PIDE gap is
**CLOSED** (rule + Layer-A fixture, walker-green), with the production-wiring and
operator generalisations scoped above. No imperative operator was added.

---

## Appendix — the shipped rule, and the operator generalisations

**Domain integral (shipped, `midpoint_1d.json`)** — `du/dt = -∫u dx`:

```jsonc
{ "applies_to": { "op": "integral", "args": ["$u"], "var": "x",
                  "lower": {"op":"const","value":0.0}, "upper": {"op":"const","value":1.0} },
  "grid_family": "cartesian",
  "replacement": { "op": "*", "args": [
    "dx",
    { "op": "arrayop", "reduce": "+",
      "ranges": { "j": [1, { "op": "index", "args": ["size_x", 1] }] },   // GAP A bypass: size as DATA
      "expr": { "op": "index", "args": ["$u", "j"] },
      "args": ["size_x"] } ] } }
```

**Compact-support integral operator** (proven shape, Check 6) — the
`nn_diffusion` envelope, `(Ku)_i = Σ_{k≤n_supp[i]} w[i,k]·u[supp[i,k]]`:

```jsonc
{ "op": "arrayop", "reduce": "+", "output_idx": ["i"],
  "ranges": { "k": [1, { "op": "index", "args": ["n_supp", "i"] }] },     // per-cell bound (variable valence)
  "expr": { "op": "*", "args": [
    { "op": "index", "args": ["w", "i", "k"] },
    { "op": "index", "args": ["u", { "op": "index", "args": ["supp", "i", "k"] }] } ] },
  "args": ["w", "supp", "n_supp"] }   // supp/w from the spatial-join front-door for a radius-r kernel
```
