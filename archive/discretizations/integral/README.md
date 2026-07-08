# Integral / PIDE Rules

Rule files for quadrature operators: domain integrals (`∫u dV`, the reduction
term of a partial integro-differential equation) and integral operators
(`(Ku)(x) = ∫K(x,y)u(y)dy`).

## Rules

- `midpoint_1d.json` — full-domain midpoint/Euler quadrature of a 1-D integral
  on a uniform cartesian grid: `∫u dx ≈ dx · Σ_j u_j`, expressed as a
  `reduce:"+"` sum_product `arrayop`. The reduction tracks the grid size via a
  host-supplied `size_x` const_array (`index(size_x, 1)`) — the same
  host-supplied-grid-data pattern `nn_diffusion` uses for mesh connectivity.

## Feasibility verdict

The declarative feasibility of this family was established under bead esd-6g4.14
(DECLARATIVE-OR-FAIL). Verdict: **FEASIBLE** — see
[`INTEGRAL_FEASIBILITY.md`](INTEGRAL_FEASIBILITY.md).

A quadrature is a weighted sum-reduction over grid cells, expressible as a
`reduce`-`arrayop` sum_product. The naive grid-spanning reduction bound fails
(the "GAP A" that blocked the spectral family), **but supplying the
integration-dimension size as a const_array bypasses it** — the reduction tracks
the grid size as data, no ESS change. Proven to machine precision through the
real `discretize → ArrayOp → eval` pipeline (`feasibility_probe/`). The shipped
rule's Layer-A canonical rewrite is byte-exact; Layer-B MMS convergence is N/A
for a quadrature (it has no differential-operator order). The bypass is
operator-agnostic, so it also unblocks the spectral family's GAP A (see the
verdict §5). Read the verdict before adding integral rules.
