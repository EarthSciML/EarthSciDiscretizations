# Integral / PIDE feasibility probe (esd-6g4.14)

Reproducible evidence for the verdict in
[`../INTEGRAL_FEASIBILITY.md`](../INTEGRAL_FEASIBILITY.md). This script drives
the **official ESS Julia evaluator** (`build_evaluator` / `discretize`) — the
same canonical pipeline the ESD walker uses — so the results reflect the real
engine, not a shadow evaluator.

This directory contains **no `*.json` rule files**, so it is invisible to
`load_rules` (`ESD:src/rules.jl` globs `*.json` in family subdirs only) and has
zero walker / rule-catalog impact.

## Running

`Project.toml` here pulls **only** EarthSciSerialization (the tree-walk path
needs no ModelingToolkit), so instantiation is seconds, not minutes:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. probe.jl
```

Pin note: the verdict was taken against ESS commit
`fd2e588318b0a8343540d4ed31a6532d99301ab0` (subdir-tree `f2b411a6…`). The
`[sources]` rev is `main`; re-confirm against `origin/main` if it has advanced.

## What the probe shows (verdict: FEASIBLE)

| check | shows |
|---|---|
| 1 | the ESS **`integral` op** is grid-agnostic host quadrature — `err = 0.0` at N = 8/16/32 (the existing capability) |
| 2 | a scalar `reduce` arrayop with a **concrete** bound reproduces the midpoint integral, `err = 0.0` (the sum_product **mechanism**) |
| 3 | the **naive grid-tied symbolic bound** (param `"Nx"` / dim `"x"`) throws `E_TREEWALK_UNBOUND_LOOP_VAR` — the apparent **GAP A** blocker |
| 4 | **★ the const_array bound `index(size_x,1)` BYPASSES GAP A ★** — the reduction tracks the grid size as DATA, grid-agnostic, `err = 0.0` at N = 8/16/32 |
| 5 | the full **RULE path**: the shipped `midpoint_1d` pattern/replacement rewrites the `integral` op → arrayop via `discretize` and evaluates, `err = 0.0` |
| 6 | the **compact-support** operator shape (per-cell bound `index(n_supp,i)`, the `nn_diffusion` envelope) executes, `err = 0.0` — the integral-OPERATOR generalisation |

Expected output:

```
--- Check 1: ESS `integral` op grid-agnostic (existing host quadrature) ---
[1: integral OP      N=8] PASS  du≡-1.0  (Linf err=0.0)
[1: integral OP      N=16] PASS  du≡-1.0  (Linf err=0.0)
[1: integral OP      N=32] PASS  du≡-1.0  (Linf err=0.0)
--- Check 2: scalar reduce-arrayop integral, CONCRETE bound (mechanism) ---
[2: arrayop concrete  N=8] PASS  du≡-1.0  (Linf err=0.0)
[2: arrayop concrete  N=16] PASS  du≡-1.0  (Linf err=0.0)
--- Check 3: NAIVE grid-tied symbolic bound fails (the apparent GAP A) ---
[3: arrayop param-Nx  N=8 ] eval THREW: E_TREEWALK_UNBOUND_LOOP_VAR: Nx
[3: arrayop dim-x     N=16] eval THREW: E_TREEWALK_UNBOUND_LOOP_VAR: x
--- Check 4: *** const_array bound index(size_x,1) BYPASSES GAP A *** ---
[4: const-array bound N=8] PASS  du≡-1.0  (Linf err=0.0)
[4: const-array bound N=16] PASS  du≡-1.0  (Linf err=0.0)
[4: const-array bound N=32] PASS  du≡-1.0  (Linf err=0.0)
--- Check 5: full RULE path (shipped midpoint_1d rewrites the integral op) ---
[5: rule path        N=8] PASS  du≡-1.0  (Linf err=0.0)
[5: rule path        N=16] PASS  du≡-1.0  (Linf err=0.0)
[5: rule path        N=32] PASS  du≡-1.0  (Linf err=0.0)
--- Check 6: compact-support operator shape (integral OPERATOR generalisation) ---
[6: compact OP       N=8] PASS  du≡-0.375  (Linf err=0.0)
[6: compact OP       N=16] PASS  du≡-0.1875  (Linf err=0.0)
```
