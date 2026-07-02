# MPAS-SCVT mesh generator host driver — RHS-only fixed-point contract

**Bead** `esd-e5m.4` (D4) · **epic** `esd-e5m` · **reference** `build_scvt_mesh` /
`scvt_lloyd_solve` / `scvt_background_quadrature` (`src/grids/mpas_scvt.jl`) ·
**step** [`lloyd_step.esm`](lloyd_step.esm) (D2) · **leaf** `scvt_voronoi_connectivity`
([`TOPOLOGY_LEAF_CONTRACT.md`](TOPOLOGY_LEAF_CONTRACT.md), D3) · **background**
[`background_quadrature.esm`](background_quadrature.esm) (D1) · **proof**
[`test/test_mpas_scvt_build.jl`](../../../../test/test_mpas_scvt_build.jl)

The MPAS-SCVT mesh generator is a per-iteration declarative FAQ **step** (D2)
driven by an external Lloyd fixed-point **loop**, capped — once, at convergence —
by the spherical-topology **leaf** (D3). This document is the normative contract
for the **loop**: the host driver that ties D1 + D2 + D3 into an `MpasMeshData`.
It is the grid-construction analogue of `build_ode_problem` (`src/ode_problem.jl`):
the engine materialises the RHS (here, ONE Lloyd iteration), and the host owns the
solve (here, the fixed-point iteration).

## 1. The driver boundary — what is host, what is declarative, what is a leaf

```
 generators₀, ρ, level (D1)
   │   bg_coord, bg_mass  =  scvt_background_quadrature(level; density=ρ)   (D1, eval_coeff)
   ▼
 ┌─[LOOP, host RHS-only — NO loop in the IR]──────────────────────────────────────┐
 │   vi      = materialize_value_invention(lloyd_step.esm, {bg_coord,bg_mass,gen}) │  ◀ D2 STEP (declarative)
 │   centroid = vi.centroid_{x,y,z}             (density-weighted Lloyd centroid)   │
 │   gen_next = R · centroid / |centroid|        (sphere re-projection, HOST √)     │
 │   until   max_g ‖gen_next[g] − gen[g]‖ < tol  (HOST convergence test)            │
 └─────────────────────────────────────────────────────────────────────────────────┘
   │   gen (converged generators on the sphere of radius R)
   ▼
   conn = scvt_voronoi_connectivity(gen; R)       (D3 LEAF, ONCE — never in the loop)
   │   + dual GEOMETRY (duo_dual_geometry_faq, the shared FAQ stack)
   ▼
   MpasMeshData
```

Per the epic's **DECLARATIVE-OR-FAIL** rule the per-iteration **step** is a
declarative semiring-FAQ, the **loop is host (RHS-only)**, and the **topology is an
allowed irreducible leaf**. The boundary is drawn at exactly two host operations:

- **Host (the loop).** The **convergence test** `max_g ‖gen_next−gen‖ < tol` and the
  **sphere re-projection** `gen_next = R·centroid/|centroid|`. The projection needs
  a Euclidean norm `|centroid| = √(Σ_d centroid_d²)`; `sqrt` is **outside the
  build-time value-invention op set** (`LLOYD_STEP_CONTRACT.md` §1) — the relational
  front-door deliberately carries no transcendental geometry — so the re-projection
  lives with the loop that consumes it, exactly the boundary at which D3 draws its
  convex-hull leaf. **No recurrence is ever lowered into the IR**: the engine
  evaluates ONE Lloyd iteration (the `lloyd_step.esm` step buffers leave the ODE at
  setup), and the host iterates it. `build_scvt_mesh` calls
  `materialize_value_invention` once per host iteration with `gen` as the only
  varying parameter.
- **Declarative (the step).** `assign → den/num → centroid` is materialised by the
  value-invention front-door (D2). Each host iteration re-evaluates the SAME
  declarative step against the updated generators.
- **Leaf (the topology).** `generators → Voronoi connectivity` runs **once** at
  convergence (D3), never in the recurrence.

## 2. The fixed-point contract

- **Unit-space iteration.** The step's argmin metric and centroid live on the unit
  sphere (`bg_coord` is the unit-sphere quadrature set, `LLOYD_STEP_CONTRACT.md`),
  so the loop iterates unit directions and restores the radius `R` only at output
  and for the downstream geometry. The argmin and the centroid direction are
  scale-invariant, so the converged directions are independent of whether the loop
  carries `R`.
- **Unattended generators are held.** An empty group (`den[g] = 0`, NaN centroid —
  the coarse-background degenerate case, `LLOYD_STEP_CONTRACT.md` §3/§4) keeps its
  previous position rather than moving to a NaN; the driver warns when the
  background is not finer than the generators. With a background finer than the
  generator set every generator is attended and no generator is held.
- **Convergence is a genuine fixed point.** At convergence the generators are a
  centroidal Voronoi tessellation: applying ONE more declarative step moves them by
  less than `tol` (the proof's fixed-point check). Under uniform density `ρ ≡ 1`
  the SCVT reduces to a CVT whose fixed point is the quasi-uniform icosahedral-dual
  mesh; a density function `ρ(x,y,z)` is the variable-resolution mechanism.
- **Non-convergence is reported, not hidden.** If the loop exhausts `max_iters`
  without reaching `tol`, the driver warns and emits the last iterate (a valid mesh
  for the current generators) rather than failing — the caller chooses `tol` /
  `max_iters`.

## 3. The emitted `MpasMeshData`

`build_scvt_mesh` emits the **same** `MpasMeshData` (`src/grids/mpas.jl`) the
imperative DUO Voronoi dual (`_duo_voronoi_dual`) produces, reusing the identical
dual TOPOLOGY (`voronoi_dual_topology_faq`, wrapped by the D3 leaf) and dual
GEOMETRY (`duo_dual_geometry_faq`) FAQ stack. The ONLY differences from the DUO
dual are the **source of the generators** (Lloyd convergence under `ρ`, not the
fixed icosahedral subdivision) and the **source of the triangulation** (the
spherical-Delaunay leaf, not the DUO primal faces):

- **cells** = the converged generators; `x/y/z_cell` are the generators scaled to
  `R`, `lon/lat_cell` their geographic coordinates.
- **dual vertices** = the spherical-Delaunay circumcentres (`conn.circumcenters`);
  `n_vertices = n_triangles = 2·Nc − 4` (Euler, closed mesh).
- **edges** = the Delaunay edges; `cells_on_edge` is each edge's generator pair,
  `n_edges = 3·Nc − 6`. The canonical edge numbering is
  `primal_topology_faq(conn.faces, Nc).edges` — the SAME numbering the D3 leaf's
  `edges_on_cell` references, so every per-edge array is index-aligned with the
  leaf connectivity.
- **areas** are the spherical-excess Voronoi-cell areas over the circumcentre ring;
  they conserve mass: `Σ_c area_cell[c] = 4πR²` to rounding.

The integer connectivity rides the D3 determinism contract (byte-identical across
bindings); the floating-point geometry rides the §5.8 tolerance contract.

## 4. Proof

[`test/test_mpas_scvt_build.jl`](../../../../test/test_mpas_scvt_build.jl) drives
`build_scvt_mesh` / `scvt_lloyd_solve` end-to-end and proves:

- **Uniform CVT** — a level-0 icosahedral seed + a level-2 background converges and
  emits the dodecahedron (12 pentagons, `n_edges = 30`, `n_vertices = 20`); mass is
  conserved (`Σ area = 4πR²`) and the adjacency passes `check_mesh(strict)`.
- **Icosahedral dual at level 1** — a 42-generator seed emits exactly 12 pentagons
  + 30 hexagons (`n_edges = 120`, `n_vertices = 80`).
- **Variable resolution** — a density `ρ = 2 + z` reweights the centroids, so the
  converged mesh (and its cell areas) differ from the uniform CVT while still
  tiling the whole sphere.
- **Host loop over the declarative step** — a deterministic off-CVT seed takes
  more than one iteration and converges to a genuine fixed point (one more
  declarative step moves the generators by `< 1e-9`), the behavioural proof that
  the loop is host RHS-only over the declarative step — no recurrence in the IR.
