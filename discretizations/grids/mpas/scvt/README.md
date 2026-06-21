# MPAS-SCVT declarative mesh generator

The Spherical Centroidal Voronoi Tessellation (SCVT) generator builds an MPAS
Voronoi mesh — including **variable resolution** via a density function `ρ` — as
a declarative per-iteration FAQ step driven by an external fixed-point loop, plus
a one-time spherical-Delaunay topology leaf (epic `esd-e5m`). The pieces:

| Bead | Piece | Artifact |
|------|-------|----------|
| `esd-e5m.1` (**this dir**) | **D1** background quadrature mesh + density sampling | `background_quadrature.esm` |
| `esd-e5m.2` | D2 one Lloyd iteration (assign → centroid → move) | *(step `.esm`, MAYOR-held on ESS E1+E2)* |
| `esd-e5m.3` | D3 spherical-topology leaf (wrap S2B) | *(leaf op, MAYOR-held on S2B)* |
| `esd-e5m.4` | D4 `build_scvt_mesh` host driver (external fixed-point loop) | *(host driver)* |
| `esd-e5m.5` | D5 conformance drain | *(conformance)* |

---

# Background quadrature mesh + density sampling as a CONST integrand FAQ (esd-e5m.1 / D1)

`background_quadrature.esm` declares the **fixed point set the SCVT step
integrates over**. One Lloyd iteration computes, per generator `g`, the
density-weighted centroid

```
c_g = ( Σ_{p→g} ρ_p dA_p x_p ) / ( Σ_{p→g} ρ_p dA_p )
```

over the background quadrature points `p` assigned to `g`, then moves `g` to
`c_g` (re-normalized to the sphere). This document declares the **static,
generator-independent half** of that integral — the quadrature point set, the
sampled density, and the per-point density-weighted integrands the step sums. The
grouped reduction over generators (and the `argmin` nearest-generator assignment)
depends on the *moving* generators and lives in the per-iteration step (D2).

## The point set is DUO-subdivision-based (reuses esd-heg.3)

Each background quadrature point **is** one DUO icosahedral primal cell of the
level-N mesh — the same mesh the
[`subdivide_fold.esm`](../../duo/faq/subdivide_fold.esm) (esd-heg.4) fold builds
by folding the one-level value-invention refine pass
([`subdivide_refine.esm`](../../duo/faq/subdivide_refine.esm), **esd-heg.3**) over
the base icosahedron seed N times. For each cell:

- **`bg_coord`** — the quadrature point — is the unit-sphere primal-cell
  **centroid** (`primal_geometry.esm` `cell_cart / R`, esd-heg.6);
- **`area`** — the quadrature weight — is the spherical-triangle **area** `dA`
  (`primal_geometry.esm` `area`, the L'Huilier spherical excess, esd-heg.6).

Both are therefore **CONST parameters** fed by the upstream subdivision + geometry
passes (value-free shape; concrete data supplied by the conformance harness /
the D4 driver), exactly as `subdivide_fold.esm` declares its `vert_raw` / `face_vert`
seed as parameters fed by the harness.

## Density sampling

**`rho`** is the density evaluated at each quadrature point — a CONST parameter:
`ρ_p = f(bg_coord_p)` for a user density function `f`, sampled through the
single-pathway evaluator (`eval_coeff`) exactly as the expression-IC rule samples
its RHS at grid cells ([`discretizations/ic/expression_ic.json`](../../../ic/expression_ic.json)).

- `ρ ≡ 1` is the **uniform-density regression**: SCVT reduces to a CVT, whose
  fixed point is the quasi-uniform icosahedral-dual MPAS mesh (the `x1.*` family).
- a non-constant `ρ` drives **variable-resolution** refinement (denser where `ρ`
  is larger). The D4 driver supplies the sampled field.

## The FAQ (the new declarative content)

Two elementwise `sum_product` aggregates (no contraction) over the cell set give
the per-point integrands the step sums:

1. **`bg_mass = ρ · area`** — the density-weighted measure `ρ_p dA_p` (the
   centroid **denominator** integrand).
2. **`bg_moment = bg_mass · bg_coord = ρ · area · bg_coord`** — the
   density-weighted position `ρ_p dA_p x_p` (the centroid **numerator**
   integrand). Reusing `bg_mass` keeps the `ρ·area` product byte-identical between
   numerator and denominator.

Products use `*` (not `^`), so the float ops are bit-stable and `ρ ≡ 1` reproduces
`area` exactly. The document is a **value-free structural spec** and validates
against the ESS `esm-schema.json`.

## What "byte-identical; matches the quadrature" means here

`fixtures/canonical/background_quadrature_level{0,1}.json` pin, for both the
uniform density (`ρ ≡ 1`) and a sampled latitude-graded density (`ρ = 2 + z`,
`z` the unit-centroid z-component, `ρ ∈ [1,3]`), the exact Float64 quadrature
points, weights, sampled densities, and integrands every binding must reproduce.
The contract is bitwise Float64 identity. Regenerate with
[`regenerate_fixtures.jl`](regenerate_fixtures.jl), which drives the same
single-pathway evaluator (per `AGENTS.md`: golden regeneration drives the
canonical pipeline).

## Proof

[`test/test_mpas_scvt_background_faq.jl`](../../../../test/test_mpas_scvt_background_faq.jl)
drives the **landed ESS engine** (`eval_coeff` — the single-pathway passthrough,
no shadow evaluator) and proves, across levels 0–2:

- `ρ ≡ 1` reduces the integrand to the bare quadrature weight, **bit-for-bit**
  (`bg_mass == area`, `bg_moment == area · centroid`);
- density sampling is exactly `ρ = 2 + z`, always positive;
- the quadrature weights **tile the sphere**: `Σ bg_mass = 4πR²` (the central
  correctness property — the background mesh integrates the constant `1` to the
  full sphere area);
- the density-weighted centroid `Σ bg_moment / Σ bg_mass` sits at the **sphere
  centre** under uniform density and **shifts toward the dense hemisphere** under
  `ρ = 2 + z` — to the analytic value `1/6` = `∫(2+z)z dA / ∫(2+z) dA` over the
  unit sphere — confirming the discrete quadrature reproduces the continuous
  density-weighted integral;
- the schema-valid declarative document and the cross-binding canonical-byte
  contract (reproduced bit-for-bit by both the FAQ and the imperative grid area).
