# MPAS `flux` operator — DECLARATIVE-OR-FAIL verdict (esd-6g4.2 / G2)

**Verdict: a standalone MPAS `flux` rule is NOT expressible in the closed
declarative vocabulary — by design. The numerical flux is REALIZED as the
`arrayop` body of the operators that produce or consume it, both shipped in
G2/G1.** This is the "OR per-op declarative-infeasibility reported" branch of
the G2 acceptance criterion.

G2 ("MPAS gradient + flux + advection") ships two full declarative rules —
`gradient_mpas` and `advection_mpas` (this directory) — plus this verdict for
the third. The companion `divergence_mpas` (G1, esd-6g4.1) is the flux-consuming
operator. Together they cover the MPAS C-grid mimetic-operator family.

## Why `flux` is not a top-level rule

The closed authoring vocabulary (`docs/content/guide/operators.md`, ESS
`esm-spec.md` §4.2) admits exactly four PDE operators in `applies_to.op`:

> `D`, `grad`, `div`, `laplacian`.

and explicitly lists `flux` as **forbidden** as an op name, with the prescribed
encoding:

> `flux` → flux assembly is an `arrayop` whose `output_idx` ranges over edges
> and whose body is the numerical-flux formula in pointwise math + `ifelse`.

So a numerical flux is the **body of an `arrayop`**, lowered from a `grad`/`div`
match — never a separately-`applies_to`-matchable rule. There is no PDE op a
standalone flux rule could pattern-match against; authoring one would require an
off-spec `flux` op, which the reviewer guidance forbids ("if the only
justification for a new op name is *this scheme calls it that*, it is the wrong
move").

## The two concrete MPAS edge fluxes ARE shipped — as arrayop bodies

| Flux | Closed form | Where it is realized |
|------|-------------|----------------------|
| Diffusive / gradient flux | `F_e = −K·(∇φ)_e = −K·(φ[c2]−φ[c1])/dc_edge[e]` | The **edge-output body of `gradient_mpas`** (this directory). For constant `K` the flux is the gradient times a const-folded scalar; `gradient_mpas` ships the `∇φ` half with its own Layer-A byte contract and an O(h²) `unstructured_gradient` Layer-B convergence proof. |
| Advective flux | `F_e = u_e·q̂_e`, `q̂_e = (q[c1]+q[c2])/2` | The **reduction body of `advection_mpas`** (this directory): the divergence sums the per-edge advective flux `u_e·q̂_e`. The flux sub-expression is the inner `index($u,…)·(index($q,c1)+index($q,c2))/2` of the `arrayop`. |
| Edge-normal flux (consumed) | `(div F)_i = (1/A_i)Σ_k s_{i,k}·dv_e·F_e` | `divergence_mpas` (G1) — the flux-CONSUMING operator. |

## Why a standalone advective-flux rule additionally fails

Even setting aside the "no `flux` op" prohibition, a rule that *emitted* the
advective flux `F_e = u_e·q̂_e` as an edge state would have to match a
velocity×scalar product `$u·$q` at edges. That product is **commutative**, and
ESS canonicalizes commutative `*` operands by **variable name**
(`canonicalize.jl::_sort_args!`), while the matcher
(`rule_engine.jl::_match`) is strictly positional with no commutative
backtracking. The edge-velocity (direct edge gather) and cell-scalar (centered
cell→edge interpolation) operands therefore land in a name-determined order and
get mis-gathered for ~half of all variable namings (verified: a positional
`div($u·$q)` rule indexes the cell scalar at edge positions).

`advection_mpas` resolves this by shipping **two `var_location_is`-guarded
variants** (operand-order mirrors) so exactly one fires for any naming — but
those guards hang off the `div` op it matches. A standalone flux has no
top-level op to hang them on, so the disambiguation has nowhere to live. The
flux is correctly expressible only *inside* a `div`/`grad` lowering, which is
exactly what the closed-vocabulary design intends.

## Bottom line

- `gradient_mpas` — **FEASIBLE, shipped** (full Layer-A + Layer-B O(h²)).
- `advection_mpas` — **FEASIBLE, shipped** (full Layer-A byte contract; Layer-B
  convergence deferred — the unlimited centered C-grid advection is L∞-low-order
  on the non-centroidal Voronoi dual, see the convergence fixture `skip_reason`).
- `flux` — **not a standalone rule**: realized as the `arrayop` body of the
  above, per the closed-vocabulary `operators.md` design. No imperative operator,
  no off-spec op.
