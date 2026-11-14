---
title: "mpas"
description: "MPAS unstructured spherical centroidal-Voronoi (icosahedral) C-grid on the unit sphere: NCELLS Voronoi cells, NEDGES edges, NVERTICES triangle vertices, with cell valence bounded by MAXEDGES (on a closed sphere the Euler relations give NEDGES = 3*(NCELLS - 2) and NVERTICES = 2*(NCELLS - 2), but the three counts stay independent metaparameters so planar or open meshes can reuse the file)."
esd_kind: "grid"
families: ["finite_volume", "mpas"]
ops: ["div", "grad", "laplacian"]
bcs: ["none"]
tags: ["esd:grid", "family:mpas", "ndim:2", "topology:unstructured_voronoi", "staggering:c-grid", "sphere:unit", "esd:rule", "family:finite_volume", "grid:mpas", "op:div", "axes:cells", "order:1", "bc:none", "spec:9.6.8", "op:grad", "axes:edges", "op:laplacian"]
source: "grids/mpas/grid.esm"
---

MPAS unstructured spherical centroidal-Voronoi (icosahedral) C-grid on the unit sphere: NCELLS Voronoi cells, NEDGES edges, NVERTICES triangle vertices, with cell valence bounded by MAXEDGES (on a closed sphere the Euler relations give NEDGES = 3*(NCELLS - 2) and NVERTICES = 2*(NCELLS - 2), but the three counts stay independent metaparameters so planar or open meshes can reuse the file). Unlike the rectangular grids there is no coordinate axis: the iteration domains are the three entity interval sets (cells, edges, vertices) plus the ragged connectivity set `edges_of_cell` (esm-schema IndexSet kind "ragged"; RFC semiring-faq-unified-ir 5.2), which binds by NAME to two keyed factors the consuming document must expose as cell-shaped arrays: `nEdgesOnCell` (offsets: the valence of each cell) and `edgesOnCell` (values: the edge ids around each cell, counterclockwise, zero-padded to maxEdges). Geometry is likewise DATA, not formula: mesh files under grids/mpas/mesh/ (section-4.7 subsystem files with literal `const` arrays) supply the factors nEdgesOnCell, edgesOnCell, edgeSignOnCell, cellsOnEdge, verticesOnEdge, dvEdge, dcEdge, areaCell, latCell, lonCell -- and, since the 2026-07 regridding wave WIDENED the mesh contract, the VERTEX-GEOMETRY factors latVertex, lonVertex (vertex = Delaunay-triangle circumcenter, radians, lon in (-pi, pi], same convention as latCell/lonCell) and verticesOnCell (vertex ids around each cell, counterclockwise viewed from outside, zero-padded to maxEdges exactly like edgesOnCell; verticesOnCell[i][k] is the vertex SHARED by consecutive CCW edges edgesOnCell[i][k] and edgesOnCell[i][k+1 cyclic]; on a closed mesh every cell has exactly as many vertices as edges, so nEdgesOnCell doubles as the per-cell vertex count and there is no separate offsets factor) -- under exactly those names; a consuming model mounts one as a subsystem and re-exposes the factors a rule needs as bare-name observed aliases (e.g. areaCell := mesh.areaCell), since keyed factors resolve by bare name in the model scope (see problems/divergence_mpas.esm for the pattern). The widened contract makes CELL POLYGON RINGS derivable in-document: the match-less template `mpas_cell_rings_deg` gathers each cell's vertex ring into the [cells, maxEdges, ring_coord] lon-lat-DEGREE stack the conservative regridder consumes as src_poly/tgt_poly (regridding/conservative_overlap.esm; ring vertex = [lon, lat] in degrees, CCW), clamping the gather with min(k, nEdgesOnCell(i)) so PADDING SLOTS REPEAT THE CELL'S LAST REAL VERTEX -- rectangular (padded) ring storage rather than a ragged ring, chosen because the overlap-area templates (esm:std/geometry.esm overlap_area_planar / overlap_area_spherical, which take the ring slot count explicitly) accept trailing-duplicate padded rings value-identically (a repeated vertex contributes exactly zero; verified empirically for both templates through the Julia, Python and Rust conformance runners: padded and unpadded rings agree at 0 ulp as subject and as clip operand; on a mesh whose MAXEDGES equals every cell's valence, e.g. the level-0 pentagons, no padding ever occurs). The other closed-form geometry is the resolution scale: `cell_area_mean` = 4*pi/NCELLS, the mean cell area on the unit sphere (the h^2 measure for convergence sweeps, h ~ sqrt(cell_area_mean)). Orientation conventions are fixed by the mesh data (positive edge normal from cellsOnEdge[e][1] to cellsOnEdge[e][2]; positive tangent v1->v2 with t = k x n; edgeSignOnCell = +1 where the positive normal leaves the cell), per Ringler et al. (2010).

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NCELLS` | integer | `12` | Number of Voronoi cells (generators). Default: the icosahedral level-0 mesh (grids/mpas/mesh/level0.esm). |
| `NEDGES` | integer | `30` | Number of edges (cell-pair interfaces; 3*(NCELLS-2) on a closed sphere). |
| `NVERTICES` | integer | `20` | Number of mesh vertices (Delaunay-triangle circumcenters; 2*(NCELLS-2) on a closed sphere). |
| `MAXEDGES` | integer | `5` | Maximum cell valence: the padded row length of the per-cell connectivity arrays (5 for level0, 6 for level1 and all larger icosahedral meshes). |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>areaCell</code> <code>cellsOnEdge</code> <code>dcEdge</code> <code>dvEdge</code> <code>edgeSignOnCell</code> <code>edgesOnCell</code> <code>nEdgesOnCell</code></p>

The mesh geometry arrives as **data**, not formula: these are the keyed-factor arrays a consuming model must expose under exactly these bare names, normally as observed aliases of a mounted `grids/mpas/mesh/` subsystem's `const` data (see `problems/divergence_mpas.esm`). The 2026-07 regridding wave widened the contract with the vertex-geometry factors so cell polygon rings derive in-document via `mpas_cell_rings_deg`.

## Index sets

| Name | Kind | Size |
|---|---|---|
| `cells` | interval | `NCELLS` |
| `edges` | interval | `NEDGES` |
| `vertices` | interval | `NVERTICES` |
| `maxEdges` | interval | `MAXEDGES` |
| `edge_ends` | interval | `2` |
| `edges_of_cell` | ragged | `per cells (offsets: nEdgesOnCell, values: edgesOnCell)` |
| `ring_coord` | interval | `2` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `cell_area_mean` | — | `12.566370614359172 / NCELLS` | Mean Voronoi-cell area on the unit sphere, 4*pi/NCELLS: the grid's resolution measure (nominal mesh spacing h ~ sqrt(cell_area_mean)), used by convergence sweeps exactly as dx = 1/N is on cartesian_uniform_1d. |
| `mpas_cell_rings_deg` | `latVertex`, `lonVertex`, `verticesOnCell`, `nEdgesOnCell` | `sum[i, k, c] (ifelse(c == 1, lonVertex[verticesOnCell[i, min(k, nEdgesOnCell[i])]] * 57.29577951308232, latVertex[verticesOnCell[i, min(k, nEdgesOnCell[i])]] * 57.29577951308232)) where {c in ring_coord, i in cells, k in maxEdges}` | Cell polygon rings from the widened mesh contract, regridder-ready: ring[i,k,c] = [lonVertex, latVertex\](verticesOnCell[i, min(k, nEdgesOnCell(i))]) * 180/pi -- the [cells, maxEdges, ring_coord] closed CCW lon-lat-DEGREE vertex-ring stack that regridding/conservative_overlap.esm consumes as src_poly/tgt_poly (ring vertex = [lon, lat] degrees; CCW is inherited from verticesOnCell's counterclockwise order). The gather is CLAMPED, min(k, nEdgesOnCell(i)), so a cell with fewer than MAXEDGES vertices pads its ring by REPEATING ITS LAST REAL VERTEX rather than dereferencing a zero pad slot -- trailing-duplicate padded rings are accepted value-identically by the overlap-area templates (esm:std/geometry.esm overlap_area_planar / overlap_area_spherical, which take the ring slot count explicitly) (a repeated vertex contributes exactly zero; empirically verified, 0 ulp, planar and spherical), which is why the ring ships rectangular instead of ragged. The four mesh factors are TEMPLATE PARAMETERS (this template is match-less and explicitly applied, so params are the natural mechanism -- unlike the div rule, whose match-fired body must use the bare-name keyed-factor contract): bind them to the mounted mesh subsystem's factors DIRECTLY under their subsystem-scoped names (e.g. latVertex: 'mesh.latVertex'), which keeps the ring expression's references const-backed -- the build-time geometry-setup evaluator materializes aggregates over const data, and routing through bare-name VarExpr aliases would put a non-const reference inside the setup closure (the alias re-exposure pass and the geometry-setup pass do not compose in that order upstream). Radians-to-degrees (factor 57.29577951308232 = 180/pi) is part of the template because the mesh files publish radians (MPAS convention) while the geometry kernels take degrees. Build-time CONST: the rings fold into the regridder's overlap-matrix setup. |

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `fv_divergence_cell`

<p class="rule-meta"><code>op:div</code> <code>order:1</code> <code>bc:none</code> <code>axes:cells</code></p>

Source: [`grids/mpas/rules/fv_divergence_cell.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mpas/rules/fv_divergence_cell.esm)

**Rewrites** `div(F)` — match pattern:

```json
{"op": "div", "args": ["F"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>F</code> is a bare field shaped <code>[edges]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sum[i] (edgeSignOnCell[i, k] * F[edgesOnCell[i, k]] * dvEdge[edgesOnCell[i, k]] / areaCell[i]) where {i in cells, k in edges_of_cell(i)}</div>

Finite-volume cell divergence of an edge-normal flux on the MPAS Voronoi grid (the TRiSK discrete divergence, Ringler et al. 2010 eq. 21): div_i = (1/A_i) * sum_{k=1..nEdgesOnCell(i)} edgeSignOnCell(i,k) * F(e) * dvEdge(e) with e = edgesOnCell(i,k) — the exact Gauss/divergence theorem on cell i with the flux integral along each edge approximated by the midpoint rule, F being the flux component along the edge's positive normal (which points from cellsOnEdge(e,1) to cellsOnEdge(e,2)) and edgeSignOnCell flipping it outward from cell i. Matches the open-namespace rewrite-target {op: div, args: [F]} (esm-spec 9.6.8: div is author sugar with no evaluable core, and a rule MAY match it directly and emit a one-shot scheme), match-scoped with `where: {F: {shape: [edges]}}` (esm-spec 9.6.1) so it fires ONLY on the divergence of a bare field declared over this grid's `edges` index set — the axis-less div pattern would otherwise fire on every div in a consuming component. The shape constraint is the mechanism that lets two MPAS meshes coexist in one document: import this rule twice under distinct index-set renames (esm-spec 9.7.7), and each instance's `where` shape follows its renamed edge set (e.g. [meshA.edges] vs [meshB.edges]), so each div lowers by its own mesh's rule with no first-declared-wins tie-break. Because the constraint requires a bare variable of exactly that shape, a consumer that divergences a COMPOUND inline flux (e.g. div(u*h)) must first bind the flux to a declared edge-shaped observed (the divergence_mpas problem's `flux` observed is exactly that pattern). The body is a single aggregate over the cells axis whose contracted index k runs over the RAGGED index set edges_of_cell (esm-spec 4.3.1 { from, of }; RFC semiring-faq-unified-ir 5.2): per-cell dynamic bound k <= nEdgesOnCell(i) plus a gather through the edgesOnCell values factor — this is the connectivity-based stencil form that replaces the makearray region form of the rectangular grids. The division by areaCell(i) sits inside the summand, which equals (1/A_i)*sum since A_i is constant over k. There are no boundary regions and no makearray because a closed sphere has no boundary: every edge separates exactly two cells (bc:none; a rule for a limited-area MPAS mesh would instead overwrite the boundary-cell entries, and would be a different rule). The free names nEdgesOnCell, edgesOnCell, edgeSignOnCell, dvEdge, areaCell are the grid's keyed-factor contract (see grids/mpas/grid.esm): the consuming model must expose them as bare-name arrays, normally as observed aliases of a grids/mpas/mesh/ subsystem's const data. Accuracy: exactly conservative for any F (sum_i A_i*div_i = 0, each edge cancels between its two cells) and exactly annihilates discrete curls (div of F_e = (eta(v2)-eta(v1))/dvEdge(e) is zero to rounding — the TRiSK div-curl identity; problems/divergence_mpas.esm asserts it); order tag is the guaranteed first-order bound on arbitrarily distorted Voronoi meshes, though the midpoint flux is second-order on centroidal/symmetric meshes (Ringler et al. 2010, section 3.2).

**References**

- Ringler, T., Thuburn, J., Klemp, J. B., & Skamarock, W. C. (2010). A unified approach to energy conservation and potential vorticity dynamics for arbitrarily-structured C-grids. Journal of Computational Physics, 229(9), 3065-3090. Eq. (21) and the div-curl mimetic identities of section 3.
- Thuburn, J., Ringler, T. D., Skamarock, W. C., & Klemp, J. B. (2009). Numerical representation of geostrophic modes on arbitrarily structured C-grids. Journal of Computational Physics, 228(22), 8321-8335.

**Conformance fixtures**

- ast: [`fv_divergence_cell`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/fv_divergence_cell) — status `active`
- convergence: [`divergence_mpas`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/divergence_mpas) — status `active`
- simulation: [`divergence_mpas`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/divergence_mpas) — status `active`

#### Convergence — case [`divergence_mpas`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/divergence_mpas)

Error norms read from the committed golden [`tests/conformance/convergence/divergence_mpas/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/divergence_mpas/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | Linf_error | observed order |
|---|---|---|
| 977 | 7.759e-01 | — |
| 1828 | 4.140e-01 | 1.00 |
| 3590 | 2.335e-01 | 0.85 |

Expected order: **1** (± 0.25).

![Convergence of divergence_mpas (log-log)](/plots/convergence/divergence_mpas.png)

### `fv_gradient_edge`

<p class="rule-meta"><code>op:grad</code> <code>order:1</code> <code>bc:none</code> <code>axes:edges</code></p>

Source: [`grids/mpas/rules/fv_gradient_edge.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mpas/rules/fv_gradient_edge.esm)

**Rewrites** `grad(phi)` — match pattern:

```json
{"op": "grad", "args": ["phi"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>phi</code> is a bare field shaped <code>[cells]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sum[e] ((phi[cellsOnEdge[e, 2]] - phi[cellsOnEdge[e, 1]]) / dcEdge[e]) where {e in edges}</div>

Finite-volume edge-normal gradient of a cell scalar on the MPAS Voronoi grid (the TRiSK discrete gradient, Ringler et al. 2010 eq. 22): grad_e = (phi(cellsOnEdge(e,2)) - phi(cellsOnEdge(e,1))) / dcEdge(e) — the component, along each edge's positive normal, of the gradient of the cell field phi, formed as the two-point difference of the flanking cell values divided by the center-to-center arc dcEdge(e). The positive normal points from cellsOnEdge(e,1) to cellsOnEdge(e,2) (the MPAS orientation shared with edgeSignOnCell and with the divergence rule's flux convention), so the difference is taken head-minus-tail in that direction and dcEdge is the great-circle distance the difference spans. This is the exact MIMETIC DUAL of grids/mpas/rules/fv_divergence_cell: that rule maps an edge-normal flux to a cell divergence over cells; this rule maps a cell scalar to an edge-normal gradient over edges, and the two are negative adjoints of one another (see Accuracy). Matches the open-namespace rewrite-target {op: grad, args: [phi]} (esm-spec 9.6.8, 4.2: grad is author sugar with no evaluable core, and a rule MAY match it directly on the bare op and emit a one-shot scheme rather than routing through the D-sugar identity grad(f,dim:x)=D(f,wrt:x)), match-scoped with `where: {phi: {shape: [cells]}}` (esm-spec 9.6.1) so it fires ONLY on the gradient of a bare field declared over this grid's `cells` index set — the axis-less grad pattern would otherwise fire on every grad in a consuming component. The shape constraint is the mechanism that lets two MPAS meshes coexist in one document: import this rule twice under distinct index-set renames (esm-spec 9.7.7), and each instance's `where` shape follows its renamed cell set (e.g. [meshA.cells] vs [meshB.cells]), so each grad lowers by its own mesh's rule with no first-declared-wins tie-break. Because the constraint requires a bare variable of exactly that shape, a consumer that takes the gradient of a COMPOUND inline cell expression (e.g. grad(a*b)) must first bind it to a declared cell-shaped observed. The body is a single aggregate over the `edges` axis with output index e and NO contracted index — a pure connectivity map, one two-point stencil per edge — because on the MPAS C-grid the gradient is intrinsically edge-local: every edge separates exactly two cells and the edge-normal derivative reads only those two, unlike the divergence whose body contracts the ragged edges_of_cell set. There is no makearray and no boundary region because a closed sphere has no boundary: every edge has both a cellsOnEdge(e,1) and a cellsOnEdge(e,2) (bc:none; a rule for a limited-area MPAS mesh would instead special-case the boundary edges that have only one interior cell, and would be a different rule). The free names cellsOnEdge, dcEdge are the grid's keyed-factor contract (see grids/mpas/grid.esm): the consuming model must expose them as bare-name arrays, normally as observed aliases of a grids/mpas/mesh/ subsystem's const data (problems/gradient_mpas.esm is the pattern). Accuracy: exactly annihilates spatially-constant fields (phi(c2)-phi(c1) = 0 for phi = const, hence grad_e = 0 at every edge to rounding — problems/gradient_mpas.esm asserts this), which is the gradient's null-space mimetic property and the discrete statement that grad(const) = 0. The rule is the negative discrete adjoint of fv_divergence_cell under the TRiSK inner products <a,b>_cells = sum_i areaCell(i)*a_i*b_i and <u,v>_edges = sum_e dvEdge(e)*dcEdge(e)*u_e*v_e: summation by parts on the closed sphere carries no boundary term, so sum_i areaCell(i)*phi_i*div(F)_i = -sum_e dvEdge(e)*dcEdge(e)*grad(phi)_e*F_e for any cell field phi and edge field F, the discrete Gauss/Green identity the two operators share (Ringler et al. 2010, section 3). The order tag is the guaranteed first-order bound on arbitrarily distorted Voronoi meshes, though the two-point difference is second-order on centroidal/symmetric meshes where the edge midpoint bisects the cell-center arc (Ringler et al. 2010, section 3.2). A convergence-order sweep on refined MPAS meshes is a possible follow-up, out of scope for the exactness gate that ships with this rule.

**References**

- Ringler, T., Thuburn, J., Klemp, J. B., & Skamarock, W. C. (2010). A unified approach to energy conservation and potential vorticity dynamics for arbitrarily-structured C-grids. Journal of Computational Physics, 229(9), 3065-3090. Eq. (22) and the discrete Green/adjoint (grad = -div^T) identities of section 3.
- Thuburn, J., Ringler, T. D., Skamarock, W. C., & Klemp, J. B. (2009). Numerical representation of geostrophic modes on arbitrarily structured C-grids. Journal of Computational Physics, 228(22), 8321-8335.

**Conformance fixtures**

- ast: [`fv_gradient_edge`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/fv_gradient_edge) — status `active`
- simulation: [`gradient_mpas`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/gradient_mpas) — status `active`

### `fv_laplacian_cell`

<p class="rule-meta"><code>op:laplacian</code> <code>order:1</code> <code>bc:none</code> <code>axes:cells</code></p>

Source: [`grids/mpas/rules/fv_laplacian_cell.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mpas/rules/fv_laplacian_cell.esm)

**Rewrites** `laplacian(phi)` — match pattern:

```json
{"op": "laplacian", "args": ["phi"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>phi</code> is a bare field shaped <code>[cells]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sum[i] (edgeSignOnCell[i, k] * dvEdge[edgesOnCell[i, k]] * ((phi[cellsOnEdge[edgesOnCell[i, k], 2]] - phi[cellsOnEdge[edgesOnCell[i, k], 1]]) / dcEdge[edgesOnCell[i, k]]) / areaCell[i]) where {i in cells, k in edges_of_cell(i)}</div>

Finite-volume cell Laplacian of a cell scalar on the MPAS Voronoi grid (the TRiSK discrete Laplacian, Ringler et al. 2010, the div-of-grad composition of eqs. 21-22): lap_i = (1/A_i) * sum_{k=1..nEdgesOnCell(i)} edgeSignOnCell(i,k) * dvEdge(e) * (phi(cellsOnEdge(e,2)) - phi(cellsOnEdge(e,1))) / dcEdge(e) with e = edgesOnCell(i,k) -- the exact Gauss/divergence theorem on cell i applied to the edge-normal gradient of phi, i.e. div(grad phi) with the flux integral along each edge approximated by the midpoint rule. This rule is the mimetic COMPOSITION of grids/mpas/rules/fv_divergence_cell with grids/mpas/rules/fv_gradient_edge: the divergence's per-edge flux F(e) is replaced INLINE by the two-point edge-normal gradient grad_e = (phi(c2) - phi(c1)) / dcEdge(e) that fv_gradient_edge would emit (c1 = cellsOnEdge(e,1), c2 = cellsOnEdge(e,2); the positive edge normal points c1 -> c2, the MPAS orientation shared by edgeSignOnCell and by both parent rules), so the summand is the outward-signed edge-normal gradient flux edgeSignOnCell(i,k)*dvEdge(e)*grad_e and the sum is the net flux through cell i's boundary, per unit area A_i. Emitting the composed one-shot stencil (rather than lowering div and grad on separate fixpoint passes) is the open-namespace option esm-spec 9.6.8 explicitly grants: 'A rule author MAY instead match laplacian directly and emit a one-shot stencil.' Matches the open-namespace rewrite-target {op: laplacian, args: [phi]} (esm-spec 9.6.8, 4.2: laplacian is author sugar with no evaluable core -- laplacian(u) = sum_i D(D(u,x_i),x_i) -- and a rule MAY match it directly on the bare op and emit a scheme rather than routing through the nested-D identity), match-scoped with `where: {phi: {shape: [cells]}}` (esm-spec 9.6.1) so it fires ONLY on the Laplacian of a bare field declared over this grid's `cells` index set -- the axis-less laplacian pattern would otherwise fire on every laplacian in a consuming component. The shape constraint is the mechanism that lets two MPAS meshes coexist in one document: import this rule twice under distinct index-set renames (esm-spec 9.7.7), and each instance's `where` shape follows its renamed cell set (e.g. [meshA.cells] vs [meshB.cells]), so each laplacian lowers by its own mesh's rule with no first-declared-wins tie-break. Because the constraint requires a bare variable of exactly that shape, a consumer that takes the Laplacian of a COMPOUND inline cell expression (e.g. laplacian(a*b)) must first bind it to a declared cell-shaped observed. The body is a single aggregate over the `cells` axis whose contracted index k runs over the RAGGED index set edges_of_cell (esm-spec 4.3.1 { from, of }; RFC semiring-faq-unified-ir 5.2): per-cell dynamic bound k <= nEdgesOnCell(i) plus a gather through the edgesOnCell values factor -- the same connectivity-based stencil form as fv_divergence_cell, since the Laplacian IS a divergence (its output is a cell aggregate contracting the cell's ragged edge ring), only with the edge flux built in-place from phi and dcEdge instead of read from a declared edge field. The division by areaCell(i) sits inside the summand, which equals (1/A_i)*sum since A_i is constant over k, exactly as in fv_divergence_cell. There are no boundary regions and no makearray because a closed sphere has no boundary: every edge separates exactly two cells, so every cellsOnEdge(e,1)/cellsOnEdge(e,2) gather and every edgeSignOnCell entry the ragged bound reaches is well defined (bc:none; a rule for a limited-area MPAS mesh would instead overwrite the boundary-cell entries -- one-sided/Neumann/Dirichlet closure of the edges with only one interior cell -- and would be a different rule). The free names edgeSignOnCell, edgesOnCell, dvEdge, cellsOnEdge, dcEdge, areaCell are the grid's keyed-factor contract (see grids/mpas/grid.esm): the consuming model must expose them as bare-name arrays, normally as observed aliases of a grids/mpas/mesh/ subsystem's const data (problems/laplacian_mpas.esm is the pattern) -- the union of fv_divergence_cell's factors (edgeSignOnCell, edgesOnCell, dvEdge, areaCell) and fv_gradient_edge's factors (cellsOnEdge, dcEdge), since the composed stencil needs both the divergence's cell-boundary connectivity and the gradient's flanking-cell gather. Accuracy: exactly ANNIHILATES spatially-constant fields (phi(c2)-phi(c1) = 0 for phi = const at every edge, hence every summand edgeSignOnCell*dvEdge*0/dcEdge = 0 and lap_i = 0 at every cell to rounding -- the Laplacian's null-space mimetic property, the discrete statement that laplacian(const) = 0; problems/laplacian_mpas.esm asserts it). Exactly CONSERVATIVE for any cell field phi (sum_i A_i*lap_i = 0: each interior edge e appears in exactly its two cells cellsOnEdge(e,1) and cellsOnEdge(e,2) with opposite edgeSignOnCell and the same dvEdge*grad_e, so the two cells' contributions cancel and the global integral of the discrete Laplacian is zero -- the discrete divergence-theorem closure inherited from fv_divergence_cell; problems/laplacian_mpas.esm asserts it for a non-constant field). And SYMMETRIC NEGATIVE-SEMIDEFINITE under the TRiSK inner products <a,b>_cells = sum_i areaCell(i)*a_i*b_i and <u,v>_edges = sum_e dvEdge(e)*dcEdge(e)*u_e*v_e: because lap = div(grad) and fv_gradient_edge is the negative discrete adjoint of fv_divergence_cell (grad = -div^T under those inner products, the closed-sphere summation-by-parts identity of section 3), the Laplacian satisfies <phi, lap psi>_cells = -sum_e dvEdge(e)*dcEdge(e)*grad(phi)_e*grad(psi)_e, which is symmetric in phi<->psi and equals -||grad phi||^2 <= 0 at psi = phi -- the discrete Dirichlet energy, the property that makes the TRiSK Laplacian a stable diffusion operator. The order tag is the guaranteed first-order bound on arbitrarily distorted Voronoi meshes, though the composed midpoint flux is second-order on centroidal/symmetric meshes where the edge midpoint bisects the cell-center arc (Ringler et al. 2010, section 3.2). A convergence-order sweep on refined MPAS meshes is a possible follow-up, out of scope for the exactness gate that ships with this rule.

**References**

- Ringler, T., Thuburn, J., Klemp, J. B., & Skamarock, W. C. (2010). A unified approach to energy conservation and potential vorticity dynamics for arbitrarily-structured C-grids. Journal of Computational Physics, 229(9), 3065-3090. Eqs. (21)-(22) and the discrete div-grad adjoint (grad = -div^T) / Dirichlet-energy identities of section 3.
- Thuburn, J., Ringler, T. D., Skamarock, W. C., & Klemp, J. B. (2009). Numerical representation of geostrophic modes on arbitrarily structured C-grids. Journal of Computational Physics, 228(22), 8321-8335.

**Conformance fixtures**

- ast: [`fv_laplacian_cell`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/fv_laplacian_cell) — status `active`
- simulation: [`laplacian_mpas`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/laplacian_mpas) — status `active`

