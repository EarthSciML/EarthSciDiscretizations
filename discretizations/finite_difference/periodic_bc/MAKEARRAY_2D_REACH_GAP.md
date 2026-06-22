# periodic_bc 2-D makearray reach gap (esd-7mj) — RESOLVED by ess-wg0

**Status: RESOLVED.** The `periodic_bc` rule is COMPLETE and CORRECT. The 1-D
integration wrap is byte-exact (`test/test_periodic_bc_rule.jl`, `maxerr == 0.0`).
The 2-D corner integration was blocked by a **pre-existing ESS bug** in
`_scan_stencil_reach!`, unrelated to the periodic rule; that bug was fixed in ESS
**ess-wg0** (the two-arm reach scan below), which is now on ESS `main`. The 2-D
corner assertions in `test/test_periodic_bc_rule.jl` are now plain `@test` (they
wrap byte-exact, corners included). This document is retained as the historical
finding for the 2-D portion (DECLARATIVE-OR-FAIL Gate D).

## What works

`periodic_bc.json` authors the periodic wrap as a makearray-region BC ghost rule
in the landed neumann/robin fn/dim `bc`-wrapper encoding:

```
bc(periodic, $side, [$u])  ->  index($u, $N - 1)      ($N = bind_side_dim_size)
```

The ess-hjg makearray lowering (`_reindex_ghost`) re-indexes the 0-based offset
`L = N-1` to the wrapped OPPOSITE end per side: `xmin -> 1+(N-1) = N -> u[N]`,
`xmax -> N-(N-1) = 1 -> u[1]`. No symbolic `Nx`, no `mod`, no new engine
primitive — only the existing `bind_side_dim_size` binding. The 1-D centered-FD
Laplacian on a NON-periodic grid with model-level `kind:"periodic"` BCs wraps
exactly (`du[1]` reads `u[N]`, `du[N]` reads `u[1]`).

## What is blocked (and why it is NOT the periodic rule)

A 2-D model whose interior stencil is the standard `laplacian_2nd_uniform_cartesian`
rule does NOT splice ANY boundary ghost — every out-of-range read falls back to
the zero-ghost convention. The rule fires correctly (all four BCs receive
`index(u, N-1)`); the failure is upstream, in the makearray box-shrink:

1. The Laplacian rule authors its neighbour reads as additive offsets,
   `index(u, i + (-1), j)` and `index(u, i + 1, j)`.
2. ESS `canonicalize` reorders the commutative `+` to put the literal first, so
   the lowered body carries `index(u, {op:+, args:[-1, "i"]}, "j")` — a
   **constant-first** offset.
3. `EarthSciSerialization._scan_stencil_reach!` (discretize.jl) only detects
   **variable-first** offsets:
   ```julia
   if length(aa) == 2 && aa[1] isa AbstractString && aa[2] isa Number
   ```
   The constant-first form `[-1, "i"]` is missed, so the per-axis reach is
   computed as 0.
4. With reach 0, `_apply_makearray_bcs!` finds no bounded side
   (`any(bounded) || return`) and returns WITHOUT wrapping the body in a
   makearray. No boundary regions are emitted, so no ghost (periodic OR
   neumann/robin/dirichlet-nonzero) is spliced in 2-D.

The 1-D `d2` rule (`centered_2nd_deriv_uniform`) escapes the bug only because it
authors one neighbour as a SUBTRACTION, `index(u, x - 1)` — non-commutative, so
canonicalize preserves the variable-first arg order and reach is detected (the
`max` over reads makes the single detected read sufficient).

The pre-existing 2-D BC goldens do not catch this because they never need a
non-trivial 2-D ghost: `bc_dirichlet_2d_nonperiodic` uses Dirichlet-0, whose
ghost EQUALS the zero-ghost fallback, and `bc_dirichlet_2d_periodic_x` wraps x
through the grid-level periodic-folding path (`_apply_periodic_folding!`), not
the makearray ghost path.

## Fix (ESS bead ess-wg0 — LANDED)

`_scan_stencil_reach!` was made canonicalization-robust by detecting the
constant-first offset symmetrically:

```julia
elseif length(aa) == 2 && aa[1] isa Number && aa[2] isa AbstractString
    k = abs(Int(aa[1])); v = String(aa[2])
    haskey(reach, v) && (reach[v] = max(reach[v], k))
```

With this two-arm reach scan the 2-D periodic Laplacian wraps byte-exactly against
the grid-periodic ground truth (corners included), through the periodic rule
unchanged. The fix is general — it also unblocks nonzero-Dirichlet / Neumann /
Robin ghosts on any 2-D stencil whose rule authors additive offsets.

Landed in ESS bead **ess-wg0** (now on ESS `main`). The 2-D corner assertions in
`test/test_periodic_bc_rule.jl` are now plain `@test` (esd-2bo flipped them from
`@test_broken` once the ESS reach fix merged).
