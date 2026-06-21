# MPAS-SCVT variable-resolution reference fixture

Bead `esd-e5m.5` (D5, conformance criterion **(c)**) of the declarative MPAS-SCVT
mesh generator (epic `esd-e5m`).

Variable resolution is the headline capability of the SCVT generator over the
fixed icosahedral DUO dual: a density function `ρ(x, y, z)` pulls more (smaller)
cells into the region where `ρ` is larger. `golden.json` pins the converged
variable-resolution mesh the canonical `build_scvt_mesh` pipeline produces, so the
var-res result stays reproducible and its density-shift behaviour is locked
against regression.

## Reference instance

| field | value |
|-------|-------|
| seed | 42-generator level-1 icosahedral (`duo_subdivide_faq(Float64, 1)`) |
| density | `ρ(x, y, z) = 2 + z` (`z` = unit-sphere z; `ρ ∈ [1, 3]`) |
| background | level 3 (1280 quadrature points — finer than the 42 generators) |
| convergence | `tol = 1e-12` (discrete centroidal fixed point) |

`ρ` is largest toward the north pole (`z → 1`), so the converged mesh has
**smaller cells there** — the variable-resolution signal the fixture pins (the
test asserts `area_cell` is negatively correlated with `z_cell`). The
uniform-density (`ρ ≡ 1` → CVT) mesh on the same seed/background is pinned
alongside so the test asserts the density genuinely reshaped the mesh.

## What is pinned, and to what contract

The integer **structure** (`n_cells` / `n_edges` / `n_vertices` / `max_edges` /
`n_edges_on_cell`) is exact. The mesh **geometry** (`area_cell`, `lat/lon_cell`,
`x/y/z_cell`) rides the **§5.8 tolerance contract** — it derives from the
libm-dependent spherical-excess areas and the centroid sums — so the test
compares to tolerance, not byte identity. The Julia reference is itself bitwise
reproducible run to run (the discrete Lloyd fixed point is deterministic), so the
golden is stable.

The **deterministic integer topology** of the var-res mesh is governed by the
same determinism contract as every other SCVT mesh and is already pinned
cross-binding by the spherical-topology leaf golden
([`../topology_leaf/`](../topology_leaf/)); it is not re-pinned here.

## VERDICT — why this is an INTERNAL reference (no external jigsaw fixture)

Criterion (c) is phrased "(optional) offline jigsaw / MPAS-Tools
variable-resolution reference fixture checked in." The standard external
references — **jigsaw-geo-python** and the **MPAS-Tools** `MpasMeshConverter` /
`mesh_definition_tools` SCVT generators — are offline command-line mesh
generators that are **not available in this environment** (no network, not
vendored, heavy non-Julia toolchains). Per the epic's DECLARATIVE-OR-FAIL culture
(an unavailable external dependency is reported, not forced), this fixture is an
**internal reference**: the canonical `build_scvt_mesh` pipeline's own
reproducible variable-resolution output.

This is sufficient for the conformance obligation because the var-res mesh's
correctness is independently anchored without the external tool:

1. **The density integral is analytically pinned** at the D1 layer
   (`test/test_mpas_scvt_background_faq.jl`): the discrete `ρ = 2 + z` quadrature
   reproduces the continuous density-weighted centroid `∫(2+z)z dA / ∫(2+z) dA =
   1/6` over the unit sphere — the variable-resolution weighting is correct at
   the source.
2. **The CVT fixed-point property** holds for the var-res mesh too
   (`test/test_mpas_scvt_conformance.jl`): each generator is the
   density-weighted centroid of its Voronoi cell at convergence.
3. **The density-shift direction** is asserted against this fixture (smaller
   cells where `ρ` is larger), and the mesh tiles the whole sphere
   (`Σ area = 4πR²`) and passes `check_mesh(strict)`.

Should jigsaw / MPAS-Tools become available, a follow-up can add an external
cross-reference fixture alongside this internal one (the tolerance contract
already covers the cross-generator drift).

## Regenerate

After an intentional, reviewed change to the SCVT pipeline:

```
julia --project=. tests/conformance/grids/mpas/scvt/variable_resolution/regenerate_golden.jl
```
