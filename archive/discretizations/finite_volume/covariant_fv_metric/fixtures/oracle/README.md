# Covariant FV operator — static oracle golden fixtures

These JSON files are **frozen numeric outputs of the ESS `grid_assembly.jl`
oracle**, captured in ESD (esd-zk9.1 / R1) so they survive the planned deletion
of `EarthSciSerialization` `grid_assembly.jl` / `grid_assembly_symbolic.jl`
(ess-4g1). They are the **conformance target** for the R2 declarative
arrayop-einsum rule (esd-zk9.2): the rule's output must be tolerance-identical
to these arrays.

## Provenance

| Field | Value |
|-------|-------|
| Oracle | `EarthSciSerialization.precompute_laplacian_stencil` / `precompute_gradient_stencil` (`src/grid_assembly.jl`) |
| ESS commit | `0e935a3384a6f13d396598a7dfebe2197e20b22c` (GitHub `main`, the rev ESD's `Project.toml` `[sources]` pins) |
| Generator | `generate_oracle_goldens.jl` (frozen copy in this dir; `generate_oracle_goldens.Project.toml` is the minimal ESS-only env) |
| Determinism | Byte-identical on re-run (Julia Float64 shortest-round-trippable printing; pure deterministic arithmetic) |

The generator defines a minimal `AbstractCurvilinearGrid` mock (mirroring ESS
`test/grid_assembly_symbolic_test.jl::SymTestCartesianGrid`) so the metric —
including a **non-orthogonal** `g^{ξη} ≠ 0` that a real orthogonal lat-lon grid
would zero out — is controlled exactly. The script will NOT run after
`grid_assembly.jl` is deleted; it is retained as a provenance record only.

## Cases

### `latlon_small.json` — orthogonal spherical (the named example)
- 8 (lon) × 6 (lat) = 48 cells, R = 1. ξ = lon ∈ [0,2π) **periodic**;
  η = lat ∈ (−π/2, π/2) **non-periodic** (poles).
- Metric: `g^{lonlon}=1/(R²cos²φ)`, `g^{latlat}=1/R²`, `g^{lonlat}=0`,
  `J=R²cosφ` (matches ESD `LatLonGrid` `metric_ginv`/`metric_jacobian`).
- **All 4 corner weights are exactly `0.0`** — lat-lon is orthogonal, so the
  cross-derivative `2g^{ξη}∂²/∂ξ∂η` vanishes and no corner gather is needed.
  The **η metric-derivative correction** (`∂(Jg^{ηη})/∂η = −sinφ ≠ 0`) IS
  exercised.
- Boundary: lon wraps (cell 1's W neighbor = cell 8); lat poles use the oracle's
  **sentinel→self** policy (south row's S-neighbor index = self; north row's
  N-neighbor index = self). See `grid.neighbors` columns.

### `nonorthogonal_small.json` — non-orthogonal curvilinear (corner-exercising)
- 8 × 8 = 64 cells, **periodic in both axes**. ξ, η ∈ [0,2π).
- Synthetic smooth SPD metric: `g^{ξξ}=1+0.2cos ξ`, `g^{ηη}=1+0.2sin η`,
  `g^{ξη}=0.25 sin(ξ+η)`, `J=1/√det(g⁻¹)` (physically consistent).
- **Exercises every weight term**, including the 4 NE/NW/SE/SW corner
  cross-derivative weights (`max|corner weight| ≈ 0.2026`) and the cross/
  orthogonal metric-derivative corrections. Verified: corner weights equal
  `±cross_d2 = ±2g^{ξη}/(4dξdη)` with the `+,−,−,+` (NE,NW,SE,SW) sign pattern
  from `grid_assembly.jl:186-189`.

## JSON schema (both files)

```
meta              : {name, ess_commit, oracle, grid_kind, note, [R]}
grid              : {Nx, Ny, periodic_xi, periodic_eta, dxi, deta,
                     flat_index, xi_centers[N], eta_centers[N]}
metric            : {ginv_xixi[N], ginv_xieta[N], ginv_etaeta[N], jacobian_J[N]}
gradient_target   : symbol passed to coord_jacobian
coord_jacobian    : {dxi_dt1[N], deta_dt1[N], dxi_dt2[N], deta_dt2[N]}
field_phi         : [N]        analytic test field φ the stencils were applied to
laplacian         : {stencil_columns[9]=C,E,W,N,S,NE,NW,SE,SW,
                     weights[N][9], neighbors[N][9] (1-based flat),
                     applied_result[N] = Σ_k weights[c,k]·φ[neighbors[c,k]]}
gradient          : {stencil_columns[5]=C,E,W,N,S,
                     weights_t1[N][5], weights_t2[N][5], neighbors[N][5],
                     applied_result_t1[N], applied_result_t2[N]}
checksums         : {sum_du_lap, sumabs_du_lap, maxabs_du_lap,
                     sumabs_grad_t1, sumabs_grad_t2}
```

Flat indexing: `c = i + (j-1)·Nx`, `i ∈ 1:Nx` (ξ), `j ∈ 1:Ny` (η); neighbor
indices are **1-based**. The `applied_result` arrays let R2 verify the full
operator end-to-end (not just the weight tables) against a known field.

See `../../ORACLE_CHARACTERIZATION.md` for the operator math, line-level
oracle provenance, and the declarative-feasibility verdict.
