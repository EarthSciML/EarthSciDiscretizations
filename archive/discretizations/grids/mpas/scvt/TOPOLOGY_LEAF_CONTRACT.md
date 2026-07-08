# MPAS-SCVT spherical-topology leaf — determinism & tolerance contract

**Bead** `esd-e5m.3` (D3) · **epic** `esd-e5m` · **canonical executor** s2bindings.rs
`s2b-s7b` · **reference** `scvt_voronoi_connectivity` (`src/grids/mpas_scvt_topology.jl`)
· **declaration** [`topology_leaf.esm`](topology_leaf.esm) · **golden**
[`tests/conformance/grids/mpas/scvt/topology_leaf/`](../../../../tests/conformance/grids/mpas/scvt/topology_leaf/)

The MPAS-SCVT mesh generator is a per-iteration declarative FAQ step (D2) driven
by an external Lloyd fixed-point loop (D4). At convergence — **once**, never in
the recurrence — the generators are turned into an MPAS Voronoi mesh by the
**spherical-topology leaf** defined here. This document is the normative contract
the leaf's output obeys; it is what makes the topology cross-binding consistent.

## 1. The leaf boundary — what is irreducible, what is declarative

```
 generators ─[LEAF: spherical Delaunay = 3-D convex hull of the unit directions]─▶
     triangles + circumcentres ─[FAQ: voronoi_dual_topology_faq, esd-heg.2]─▶
         cells_on_cell / edges_on_cell / vertices_on_cell
```

Per the epic's **DECLARATIVE-OR-FAIL** rule the SCVT *step* (D2) must be a
declarative semiring-FAQ, but the topology is **"an allowed irreducible leaf
(like `intersect_polygon`)."** The boundary is drawn at exactly one operation:

- **Irreducible (the leaf).** `generators → Delaunay triangulation`. On a sphere
  the Delaunay triangulation of a point set is exactly the set of faces of the
  points' 3-D convex hull (taken as unit vectors). This is an iterative,
  control-flow-heavy construction whose robustness lives in *how* the orientation
  arithmetic is evaluated — it **cannot** be expressed as a semiring aggregate,
  the same status the RFC gives `acos`/`sqrt` and the ESS `intersect_polygon`
  geometry leaf (RFC `semiring-faq-unified-ir` §8.1).
- **Declarative (already landed).** `triangulation → MPAS connectivity`. The
  angular ordering + bridge-vertex set-intersection that turn a triangulation
  into `cells_on_cell`/`edges_on_cell` is the value-invention FAQ
  `voronoi_dual_topology_faq` (bead `esd-heg.2`,
  [`discretizations/grids/duo/rules/voronoi_dual_topology.esm`](../../duo/rules/voronoi_dual_topology.esm)).
  It previously consumed the hard-wired DUO primal triangulation; the leaf now
  supplies the **arbitrary-generator** triangulation it consumes.

## 2. Canonical executor: the S2B FFI

In production the triangulation is emitted by the **s2bindings.rs S2B FFI**
(`SphericalDelaunay::from_lon_lat` → `voronoi_connectivity`, bead `s2b-s7b`).
S2B resolves every orientation decision with an **exact predicate** — double
precision under a Shewchuk static error filter, falling back to s2geometry's
arbitrary-precision `ExactFloat` — so the integer connectivity has no
floating-point tie ambiguity and is independent of insertion order, **even at
the cospherical degeneracies** (cocircular quads, the 12 icosahedral pentagon
seeds) where a naive predicate would diverge between bindings.

This Julia reference binding computes the same convex hull with a Float64
orientation predicate — exact for the **non-degenerate** inputs §4 requires — and
emits the **identical canonical ordering** (§3), so the Julia reference and the
S2B-backed bindings agree byte-for-byte on any well-posed mesh. This is the same
pattern `intersect_polygon` uses (GeometryOps in Julia, `s2bindings` in Rust)
under one shared contract: **the contract, not a shared implementation, is the
cross-binding glue.**

## 3. Determinism contract — the INTEGER connectivity (byte-identical)

`triangle_cell` and the dual `cells_on_cell` / `edges_on_cell` /
`vertices_on_cell` are a **pure function of the generator coordinates**, emitted
byte-identically by every binding (ESS `CONFORMANCE_SPEC.md` §5.5 / §5.7 — the
same regime as `distinct`/`rank`/`skolem`). The normative rules:

1. **Triangles (`triangle_cell`).** Each Delaunay triangle is emitted CCW as seen
   from outside the sphere (outward normal `n` satisfies `n · a > 0` for a hull
   enclosing the centre), then **rotated so its smallest cell index is first**,
   and the whole triangle list is **sorted lexicographically**. Voronoi vertex
   `t` (the circumcentre) is dual to triangle `t`.
2. **Neighbour rings (`cells_on_cell`).** Per cell, the neighbours are in CCW
   angular order, **rotated to start at the smallest neighbour index**.
3. **Vertex rings (`vertices_on_cell`).** Index-aligned with `cells_on_cell` per
   the S2B convention: `vertices_on_cell[i]` is the circumcentre of the triangle
   between neighbours `cells_on_cell[i]` and `cells_on_cell[(i + 1) mod k]`.
4. **Edges are canonical skolem tuples.** An undirected Delaunay/Voronoi edge is
   `(min(u, v), max(u, v))`; `edges_on_cell` holds the dense id of that canonical
   edge (the `voronoi_dual_topology_faq` / `primal_topology_faq` numbering), never
   a hash- or insertion-order id (§5.7 rule 4).
5. **No floats in keys.** The angular sort key is a numeric `atan2` bearing; it is
   a non-key float operation materialised in the bridge (a stable sort), never a
   relational key (§5.7 rule 1). Its cross-binding identity rides the geometry
   pipeline's fixed evaluation order, pinned by the golden.
6. **Closed-mesh / Euler check.** A closed spherical mesh of `n` generators has
   exactly `n_tri = 2n − 4` triangles; the leaf **hard-fails** any other count
   (degenerate or non-enclosing generators) rather than guessing.
7. **Index base** is pinned per binding (Julia 1-based, Rust/Python 0-based); the
   golden is serialised 0-based binding-neutral and converted at the harness
   boundary (`rank_base_pin` in `golden.json`).

## 4. Well-posedness (when byte-identity holds)

Byte-identity across bindings is guaranteed for a **well-posed** SCVT mesh: the
generators are distinct, finite, and their convex hull **encloses the sphere
centre** (as any SCVT mesh's do), with **no four exactly cospherically coplanar**
generators. Coordinates produced by the Lloyd step from lon/lat trigonometry are
never *exactly* coplanar, so a converged SCVT mesh is always well-posed. An
exactly-degenerate configuration is **reported as an error** (S2B) or trips the
Euler check (the Julia reference) — never resolved by an ad-hoc rule. At a
genuine cocircular degeneracy the byte-identity guarantee requires S2B's exact
predicate; the Float64 reference is exact only for the non-degenerate case and
defers the degenerate tail to S2B.

## 5. Tolerance contract — the FLOATING-POINT geometry

`circumcenter` (the dual Voronoi vertices) and any downstream Voronoi cell
**areas** / edge **lengths** are floating point and ride the **tolerance**
contract (ESS `CONFORMANCE_SPEC.md` §5.8), **not** byte-identity:

- Two bindings are compared only under the **same manifold** (spherical).
- A geometric quantity matches within `|x − x_ref| ≤ atol + rtol · |x_ref|`, with
  the sliver floor `atol ≈ 1e-15 · R²` for areas (`≈ 1e-15 · R` for lengths);
  sub-`atol` quantities are treated as equal-to-zero.
- The defining invariant — each circumcentre is **equidistant from its triangle's
  three generators** and lies on the sphere of radius `R` — is asserted to
  tolerance (`topology_leaf.esm` `circ_dist2`; `test_mpas_scvt_topology_leaf.jl`).

## 6. Proof

[`test/test_mpas_scvt_topology_leaf.jl`](../../../../test/test_mpas_scvt_topology_leaf.jl)
proves, on two non-degenerate seeds (the octahedron, 6 → 8; icosahedral level 1,
42 → 80):

- exact connectivity + byte-identical canonical serialisation vs the golden;
- structural invariants — Euler `2n − 4`, every edge in exactly two triangles,
  CCW-from-outside + smallest-index-first + lexicographic order, adjacency
  symmetry, and the §3.3 `vertices_on_cell` ↔ `cells_on_cell` S2B alignment;
- **determinism** — permuting the generator order leaves the triangulation (as a
  set of generator-identity triples) unchanged;
- the **ρ ≡ 1 regression** — on the icosahedral seed the leaf's `cells_on_cell`
  is byte-identical to the landed imperative `_duo_voronoi_dual`, i.e. uniform-
  density SCVT recovers exactly the quasi-uniform icosahedral-dual MPAS mesh;
- the tolerance geometry invariant (circumcentres on the sphere of radius `R`);
- the [`topology_leaf.esm`](topology_leaf.esm) declaration is schema-valid.
