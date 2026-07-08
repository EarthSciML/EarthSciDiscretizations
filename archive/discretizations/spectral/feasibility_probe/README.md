# Spectral feasibility probe (esd-6g4.4)

Reproducible evidence for the verdict in
[`../SPECTRAL_FEASIBILITY.md`](../SPECTRAL_FEASIBILITY.md). These scripts drive
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
julia --project=. probe.jl       # mechanism + Fourier 1st-deriv (machine precision)
julia --project=. probe2.jl      # rule -> discretize -> eval pipeline (machine precision)
julia --project=. probe3.jl      # Chebyshev 1st-deriv + Fourier generality (N=16)
julia --project=. gapA_demo.jl   # GAP A: symbolic reduction-range bound -> E_TREEWALK_UNBOUND_LOOP_VAR
```

Pin note: the verdict was taken against ESS commit
`60c1ca89423a95418d73517d01e353d5b1551715` (subdir-tree `8542c29a…`). The
`[sources]` rev is `main`; re-confirm against `origin/main` if it has advanced.

## What each probe shows

| script | shows |
|---|---|
| `probe.jl` | dense `reduce` arrayop executes over a state vector; reduction index usable as a scalar; **Fourier 1st-derivative matrix inline in the AST → L∞≈3e-16** |
| `probe2.jl` | a Fourier `grad` **rule** is applied by `discretize` and evaluates to machine precision for a concrete grid size (full rule→ArrayOp→eval path) |
| `probe3.jl` | **Chebyshev** 1st-derivative inline in the AST → exact (L∞≈4e-15); Fourier generality at N=16, wavenumber 2 |
| `gapA_demo.jl` | **GAP A**: a reduction range tied to the grid size (`"Nx"` / `"x"`) throws `E_TREEWALK_UNBOUND_LOOP_VAR` — reduction bounds must be concrete integers, so no _general_ rule is expressible |
