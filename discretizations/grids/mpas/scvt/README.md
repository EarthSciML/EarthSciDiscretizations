# MPAS-SCVT declarative mesh generator

The Spherical Centroidal Voronoi Tessellation (SCVT) generator builds an MPAS
Voronoi mesh — including **variable resolution** via a density function `ρ` — as
a declarative per-iteration FAQ step driven by an external fixed-point loop, plus
a one-time spherical-Delaunay topology leaf (epic `esd-e5m`). The pieces:

| Bead | Piece | Artifact |
|------|-------|----------|
| `esd-e5m.1` (**this dir**) | **D1** background quadrature mesh + density sampling | `background_quadrature.esm` |
| `esd-e5m.2` (**this dir**) | **D2** one Lloyd iteration (assign → centroid → move) | `lloyd_step.esm`, `LLOYD_STEP_CONTRACT.md` |
| `esd-e5m.3` (**this dir**) | **D3** spherical-topology leaf (wrap S2B) + determinism/tolerance contract | `topology_leaf.esm`, `TOPOLOGY_LEAF_CONTRACT.md` |
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

---

# One-iteration Lloyd step (assign → centroid → move) (esd-e5m.2 / D2)

`lloyd_step.esm` declares **one** Lloyd / SCVT iteration over the fixed D1
background quadrature point set, evaluated end-to-end by the EarthSciSerialization
**value-invention front-door** (`materialize_value_invention`) — composing the two
MAYOR-held ESS keystones, **E1** (`argmin` arg-witness, `ess-os1`) and **E2**
(grouped `sum_product` over a derived key, `ess-2u5`), both now landed.

## The three declarative step pieces

Given the background points `bg_coord`, their density-weighted measure `bg_mass`
(both from D1), and the current generators `gen`, one iteration is:

1. **Assign (E1).** `assign[c] = argmin_g dist²(bg_coord[c], gen[g])` — per
   background point the nearest-generator **index**, a squared-Euclidean distance
   (so the metric is `*`/`−`/`+` only, no `sqrt`) with the §5.7 rule-6
   **smallest-generator-id tie-break**, making the integer buffer byte-identical
   across bindings.
2. **Centroid (E2).** `den[g] = Σ_{c:assign[c]=g} bg_mass[c]` and
   `num_{x,y,z}[g] = Σ_{c:assign[c]=g} bg_mass[c]·bg_coord[c,·]` — grouped
   `sum_product` reductions keyed on the data-dependent `assign` buffer, run
   through `Relational.group_aggregate`. One scalar reduction per **space**
   component, because a grouped reduction emits a single output index.
3. **Move.** `centroid_{x,y,z}[g] = num_{x,y,z}[g] / den[g]` — the derived
   elementwise division: the generator's next (pre-projection) position. An
   unattended generator (`den = 0`) folds to `NaN`; the host keeps its old seed.

All eight buffers are value-invention — materialised **once** at build time, off
the per-step hot path, and dropped from the ODE.

## DECLARATIVE-OR-FAIL: where the boundary falls

Per the epic rule the per-iteration **step** is the declarative FAQ above; the
**loop is host (RHS-only)** and the topology is the D3 leaf. The final **sphere
re-projection** that turns the centroid into the next generator —
`gen_next = R · centroid / |centroid|` — needs a Euclidean norm (`sqrt`), which is
**outside the value-invention build-time op set** (`index`/`skolem`/`floor`/`ceil`/
`+`/`−`/`*`/`/`); the relational front-door carries no transcendental geometry, the
same boundary at which D3 draws its convex-hull leaf. So the projection is the
**host loop's RHS-only step** (it reads the centroid, writes the next generators,
and re-invokes the build), documented in
[`LLOYD_STEP_CONTRACT.md`](LLOYD_STEP_CONTRACT.md) — not emitted as a buffer. The
centroid `num/den` **is** the planar Lloyd move (exactly as the landed E2 fixture
frames it); the projection is the only sphere-specific addition.

## What "byte/tolerance-identical across bindings" means here

- **Determinism (integer).** `assign` is byte-identical across bindings (§5.7
  rule 6 tie-break) and density-independent (a pure distance argmin).
- **Tolerance (float).** `den` / `num_*` / `centroid_*` ride the §5.8 tolerance
  contract across bindings (the grouped float sums derive from the libm-dependent
  spherical-excess `bg_mass`); the Julia reference reproduces them bit-for-bit.

The full contract is **[`LLOYD_STEP_CONTRACT.md`](LLOYD_STEP_CONTRACT.md)**.

## Proof

[`test/test_mpas_scvt_lloyd_step_faq.jl`](../../../../test/test_mpas_scvt_lloyd_step_faq.jl)
drives `materialize_value_invention` end-to-end and proves: a **planar
regression** reduces the 3-D step bit-for-bit to the landed ESS E2 result
(`centroid_x == [0, 1.125, 2.0]`); the **level-0 composition** (D1's 20 points × 12
icosahedral generators) has `assign` equal to an independent nearest-generator
recomputation, conserves mass (`Σ den = Σ bg_mass`), and folds unattended
generators to `NaN`; a **CVT regression** shows refinement attends every generator
and shrinks the one-step move (the icosahedral generators are a centroidal fixed
point); and the schema-valid declaration + the level-0 canonical byte golden
([`fixtures/canonical/lloyd_step_level0.json`](fixtures/canonical/lloyd_step_level0.json),
regenerate with [`regenerate_lloyd_step_fixtures.jl`](regenerate_lloyd_step_fixtures.jl))
are reproduced by the front-door for both `ρ ≡ 1` and `ρ = 2 + z`.

---

# Spherical-topology leaf (esd-e5m.3 / D3)

The one-time, post-convergence topology operation. Given the converged SCVT
generators it emits the spherical **Delaunay** triangulation and its dual
**Voronoi** connectivity (`cells_on_cell` / `edges_on_cell` / the circumcentre-
ring `vertices_on_cell`). It runs ONCE at convergence — never in the Lloyd
recurrence — and the D4 driver assembles its output into `MpasMeshData`.

## The leaf boundary

```
generators ─[LEAF: spherical Delaunay = 3-D convex hull of unit directions]─▶
    triangles + circumcentres ─[FAQ: voronoi_dual_topology_faq, esd-heg.2]─▶
        cells_on_cell / edges_on_cell / vertices_on_cell
```

Per the epic's DECLARATIVE-OR-FAIL rule the SCVT *step* (D2) is a declarative
FAQ, but the topology is **"an allowed irreducible leaf (like
`intersect_polygon`)."** Only the `generators → triangulation` step is
irreducible (an iterative convex-hull construction whose robustness lives in how
the orientation arithmetic is evaluated); everything below it is the **landed**
declarative dual-topology FAQ (`voronoi_dual_topology_faq`,
[`../../duo/rules/voronoi_dual_topology.esm`](../../duo/rules/voronoi_dual_topology.esm)),
which previously hard-wired the DUO primal triangulation and now consumes the
leaf's arbitrary-generator triangulation.

## Canonical executor + cross-binding consistency

The triangulation is emitted in production by the **s2bindings.rs S2B FFI**
(`SphericalDelaunay`, bead `s2b-s7b`), which uses an exact orientation predicate
so the integer connectivity is byte-identical across bindings even at cospherical
degeneracies. The Julia reference (`scvt_spherical_delaunay` /
`scvt_voronoi_connectivity`, [`src/grids/mpas_scvt_topology.jl`](../../../../src/grids/mpas_scvt_topology.jl))
computes the same hull with a Float64 predicate — exact for the non-degenerate
inputs the contract requires — and emits the identical canonical ordering, so the
two agree byte-for-byte on a well-posed mesh (the `intersect_polygon` pattern:
native per-binding kernels under one contract).

## Determinism + tolerance contract

The full normative contract is **[`TOPOLOGY_LEAF_CONTRACT.md`](TOPOLOGY_LEAF_CONTRACT.md)**:

- **Determinism (integer connectivity).** Triangles CCW-from-outside, rotated
  smallest-index-first, lexicographically sorted; dual rings CCW rotated to the
  smallest neighbour, `cells_on_cell` / `vertices_on_cell` index-aligned;
  undirected edges the canonical skolem `(min, max)`; Euler `n_tri = 2n − 4` hard-
  failed otherwise. Byte-identical across bindings (`CONFORMANCE_SPEC.md` §5.5 /
  §5.7).
- **Tolerance (geometry).** The circumcentre coordinates (and downstream areas /
  edge lengths) ride the §5.8 tolerance contract (`atol ≈ 1e-15 R²`), not byte
  identity.

The declarative companion is [`topology_leaf.esm`](topology_leaf.esm) (schema-
valid; declares the leaf's emitted connectivity + the circumcentre-equidistance
tolerance witness).

## Proof

[`test/test_mpas_scvt_topology_leaf.jl`](../../../../test/test_mpas_scvt_topology_leaf.jl)
proves, on the octahedron (6 → 8) and icosahedral level-1 (42 → 80) seeds: exact
connectivity + the byte golden
([`tests/conformance/grids/mpas/scvt/topology_leaf/`](../../../../tests/conformance/grids/mpas/scvt/topology_leaf/));
structural invariants (Euler, two-triangles-per-edge, canonical ordering,
adjacency symmetry, S2B vertex alignment); **determinism** (permuting the
generators leaves the triangulation invariant); the **ρ ≡ 1 regression**
(`cells_on_cell` byte-identical to the imperative `_duo_voronoi_dual` — uniform-
density SCVT recovers the quasi-uniform icosahedral-dual mesh); the tolerance-
geometry invariant; and the schema-valid declaration.
