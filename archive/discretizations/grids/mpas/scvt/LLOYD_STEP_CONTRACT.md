# MPAS-SCVT one-iteration Lloyd step — declarative boundary, determinism & tolerance contract

**Bead** `esd-e5m.2` (D2) · **epic** `esd-e5m` · **engine** EarthSciSerialization
`materialize_value_invention` (the value-invention front-door) · **keystones**
ESS `argmin` (E1, `ess-os1`) + grouped `sum_product` (E2, `ess-2u5`) ·
**declaration** [`lloyd_step.esm`](lloyd_step.esm) · **golden**
[`fixtures/canonical/lloyd_step_level0.json`](fixtures/canonical/lloyd_step_level0.json)

The MPAS-SCVT mesh generator is a per-iteration declarative FAQ **step** (this
bead) driven by an external Lloyd fixed-point **loop** (D4), capped — once, at
convergence — by the spherical-topology **leaf** (D3,
[`TOPOLOGY_LEAF_CONTRACT.md`](TOPOLOGY_LEAF_CONTRACT.md)). This document is the
normative contract for one iteration of the step: what is declarative, what is
host, and what cross-binding identity each output obeys.

## 1. The step boundary — what is declarative, what is host

One Lloyd iteration over the **fixed** background quadrature point set declared by
[`background_quadrature.esm`](background_quadrature.esm) (D1) computes, per
generator `g`, the density-weighted centroid and moves `g` to it:

```
 bg_coord, bg_mass (D1)  +  gen (current generators)
   │
   ├─[STEP, declarative FAQ — value-invention front-door]────────────────────────┐
   │   (1) ASSIGN   assign[c] = argmin_g dist2(bg_coord[c], gen[g])      (E1)     │
   │   (2) CENTROID den[g]      = Σ_{c:assign[c]=g} bg_mass[c]            (E2)     │
   │                num_{x,y,z}[g] = Σ_{c:assign[c]=g} bg_mass[c]·bg_coord[c,·]    │
   │   (3) MOVE     centroid_{x,y,z}[g] = num_{x,y,z}[g] / den[g]                  │
   └──────────────────────────────────────────────────────────────────────────────┘
       │ centroid (the next pre-projection generator position, in R³)
   ┌─[LOOP, host RHS-only]──────────────────────────────────────────────────────┐
   │   PROJECT  gen_next[g,:] = R · centroid[g,:] / |centroid[g,:]|              │
   │   repeat the build with gen := gen_next until ‖gen_next − gen‖ < tol         │
   └──────────────────────────────────────────────────────────────────────────────┘
```

Per the epic's **DECLARATIVE-OR-FAIL** rule the per-iteration **step** must be a
declarative semiring-FAQ; the **loop is host (RHS-only)** and the **topology is an
allowed irreducible leaf**. The boundary is drawn at exactly two operations:

- **Declarative (the step).** `assign → den/num → centroid` is materialised
  end-to-end by the value-invention front-door (`value_invention.{jl,rs,py}`),
  **once** at build time, off the per-step hot path (RFC `semiring-faq-unified-ir`
  §6.1), and dropped from the ODE — exactly as the E1/E2 fixtures
  ([`nearest_generator_argmin.esm`](https://github.com/EarthSciML/EarthSciSerialization)
  / `nearest_generator_centroid.esm`) are. `assign` is the E1 arg-witness
  (§5.7 rule 6); `den`/`num_*` are E2 grouped `sum_product` reductions keyed on
  the `assign` buffer (§5.5 rule 5, via `Relational.group_aggregate`);
  `centroid_*` is the E2 derived elementwise `/`. The only irreducible leaf
  *inside* the step is the E1 `argmin` itself (like `intersect_polygon`, §8.1).
- **Host (the loop).** The **sphere re-projection** `centroid → gen_next` and the
  **fixed-point iteration** are host RHS-only code (the D4 driver). The projection
  needs a Euclidean norm `|centroid| = √(Σ_d centroid_d²)`; `sqrt` is **outside
  the value-invention build-time op set** (`index`, `skolem`, `floor`, `ceil`,
  `+`, `−`, `*`, `/`, comparisons — the relational front-door deliberately carries
  no transcendental geometry, the same boundary at which D3 draws its convex-hull
  leaf and the RFC draws `acos`/`sqrt`). It is a trivial elementwise map the host
  applies between builds; it is **not** emitted as a value-invention buffer, and a
  `sqrt` derived buffer reading the centroid is correctly rejected by the engine
  (`E_TREEWALK_VI_OP`). The centroid IS the planar Lloyd move; the projection is
  the only sphere-specific addition, and it lives with the loop that consumes it.

### Why the centroid is three scalar reductions

A grouped semiring reduction emits a **single output index** (the group key `g`),
so the 3-D moment `Σ_{c→g} bg_mass[c]·bg_coord[c,:]` cannot be one
`[generators, space]` aggregate — it is the three scalar reductions
`num_x`/`num_y`/`num_z`, each a group-by over `g`, and likewise three derived
`centroid_*`. The squared distance in the `argmin` is written as the three
explicit axis terms (literal component indices `1`/`2`/`3`) for the same reason:
the arg-witness body is a flat scalar FAQ over `*`/`−`/`+`.

## 2. Determinism contract — the INTEGER assignment (byte-identical)

`assign` is a **pure function of `bg_coord` and `gen`**, emitted byte-identically
by every binding (ESS `CONFORMANCE_SPEC.md` §5.5.1 rule 6 — the same regime as
`distinct`/`rank`/`skolem`):

1. **Squared metric.** `dist2 = Σ_d (bg_coord[c,d] − gen[g,d])²` uses `*`/`−`/`+`
   only — `argmin` is invariant to the monotone square, so no `sqrt` is needed and
   the metric is a deterministic FAQ.
2. **Smallest-id tie-break.** Equidistant generators resolve to the **smallest
   generator id** (§5.7 rule 6): equal `dist2` values resolve to the lower index,
   making the integer buffer independent of candidate enumeration order. (The
   reference test reproduces `assign` with an independent ascending-`g` argmin
   using the identical left-associated arithmetic — they agree byte-for-byte.)
3. **Density-independent.** `assign` depends only on geometry, never on `bg_mass`,
   so it is identical under any density `ρ` (the golden stores it once).
4. **Empty candidate set is an error.** A point with no candidate generator has no
   witnessing index (`E_TREEWALK_VI_ARGEMPTY`) — never an ad-hoc default.
5. **Broad-phase prune (optional, at scale).** Replacing the unpruned candidate
   set with a bin-Skolem `join` (the landed `NearestGeneratorBinned` idiom) prunes
   candidates to a point's lat-lon bin; it changes performance, not the result
   (the same nearest generator, the same tie-break).

## 3. Tolerance contract — the FLOATING-POINT grouped / derived buffers

`den`, `num_{x,y,z}`, and `centroid_{x,y,z}` are floating point and ride the
**tolerance** contract (ESS `CONFORMANCE_SPEC.md` §5.8), **not** byte-identity
across bindings:

- The grouped `⊕` sums **sequentially in canonical key order** (§5.5 rule 5), so a
  single binding is bit-reproducible; across bindings the sums ride §5.8 because
  `bg_mass` itself derives from the spherical-excess area (a libm-dependent
  quantity, per D1's tolerance contract).
- **Mass is conserved**: `Σ_g den[g] = Σ_c bg_mass[c]` to rounding — the grouped
  reduction partitions all of the background mass among the generators.
- **Empty groups fold to 0̄.** A generator with no assigned point has `den = 0`
  and `num = 0` (the empty-`⊕` identity); its `centroid = 0/0 = NaN` — the host
  loop keeps that generator's previous seed. The golden serialises an unattended
  centroid as JSON `null`.
- The density-weighted centroid is a convex combination of unit-sphere points, so
  `|centroid[g]| ≤ R` (it lies inside the ball; the host projects it back out).

## 4. CVT well-posedness

Under uniform density (`ρ ≡ 1`) the SCVT reduces to a CVT, whose fixed point is
the quasi-uniform icosahedral-dual mesh. The icosahedral generators are therefore
a **centroidal Voronoi fixed point**: as the background quadrature refines, one
Lloyd step's projected displacement shrinks toward zero. A background **finer than
the generator set** attends every generator (no empty groups); the **coarse
level-0 quadrature** (20 points, 12 generators) is the degenerate edge case — each
face centroid is equidistant from its three corner generators, so the smallest-id
tie-break starves 4 of the 12 generators (the empty-group path above). Both
regimes are deterministic and pinned.

## 5. Proof

[`test/test_mpas_scvt_lloyd_step_faq.jl`](../../../../test/test_mpas_scvt_lloyd_step_faq.jl)
drives `materialize_value_invention` end-to-end and proves:

- **Planar regression** — embedding the landed E2 fixture's geometry in 3-D
  reproduces its result bit-for-bit (`assign = [1,2,2,3]`, `centroid_x =
  [0, 1.125, 2.0]`, `centroid_y = centroid_z = 0`); all eight step buffers leave
  the ODE.
- **Level-0 composition** — D1's 20 background points × the 12 icosahedral
  generators: `assign` equals an independent nearest-generator argmin (the
  smallest-id tie-break), mass is conserved, `centroid = num/den` for every
  attended generator and `NaN` for the 4 unattended ones, and the centroids lie in
  the ball; a sampled `ρ = 2 + z` reweights the per-generator mass.
- **CVT regression** — refining the background (levels 1→2) attends every
  generator and **shrinks** the one-step projected displacement (the icosahedral
  generators are a centroidal fixed point).
- **Canonical byte contract** — the schema-valid [`lloyd_step.esm`](lloyd_step.esm)
  and the level-0 golden, reproduced by the front-door: `assign` byte-identical,
  the grouped/derived floats bit-for-bit by the Julia reference (tolerance across
  bindings), for both `ρ ≡ 1` and `ρ = 2 + z`.
