---
title: "Regridding"
description: "Conservative and interpolating remap expressions (coupling transforms). Each entry is a pure ESM template library — declarative expressions, no kernels."
---

Conservative and interpolating remap expressions (coupling transforms). Each entry is a pure ESM template library — declarative expressions, no kernels.

The conservative overlap regridder is **end-to-end declarative**: a consumer constructs source and target cell rings and runs the candidate-gated overlap → normalize → apply, all from `.esm` templates with no host glue.

- **Cell-ring constructors** turn a grid spec into the `[cells, verts, coord]` vertex rings the regridder consumes — `cartesian_cell_rings.esm` from a uniform cartesian spec, and `grids/mpas/grid.esm`'s `mpas_cell_rings_deg` from MPAS mesh vertex data (the widened contract).
- **Broad phase in-library**: a `join.overlap` gate on the ring envelopes (bounding boxes the engine derives from the same rings the narrow phase clips), so the candidate set is byte-identical across bindings with no binning parameter.
- **Gated == dense** value-identically (the gated denominator broad-phases like the numerator, so partition-of-unity stays exact).
- **One file, both manifolds**: the consumer picks the narrow phase — `conservative_overlap_A_ij_planar` or `conservative_overlap_A_ij_spherical` (over unit-vector rings from `conservative_overlap_unit_rings`) — and everything downstream of A_ij is shared. Both run on Julia + Python + Rust; the two rewrite-only ports (Go, TypeScript) are scope-excluded — recorded in each case manifest, never shimmed.

End-to-end exemplars: `tests/conformance/regridding/cartesian_rings_regrid_gated_3x3_to_2x2` (grid spec → rings → broad phase → gated overlap → regrid) and `mpas_l0_to_octants_sphere` (MPAS mesh-derived rings on the sphere).

### `cartesian_cell_rings`

<p class="rule-meta"><code>method:cartesian_cell_rings</code></p>

Source: [`regridding/cartesian_cell_rings.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/regridding/cartesian_cell_rings.esm)

Declarative cell GEOMETRY for a uniform cartesian grid, in regridder-ready form: given the grid spec -- per-axis origin x0/y0 (the west-most / south-most cell-EDGE coordinate), spacing dx/dy, and sizes GX/GY -- the templates emit the per-cell representative centers (cartesian_cell_rings_lon / _lat, shape [ring_cells]) and the per-cell four-corner closed CCW rings (cartesian_cell_rings_poly, shape [ring_cells, ring_verts, ring_coord]) that regridding/conservative_overlap.esm consumes as its src_poly/tgt_poly (narrow phase) and its binning coordinates (broad phase). This is the 0.8.0 template-library rewrite of the archived model-form file (archive/regridding/cartesian_target_cells.esm), with two deliberate changes: (1) the cell axis is FLATTENED to the single ring_cells = GX*GY enumeration the regridder's src_cells/tgt_cells axes expect -- row-major from the southwest, cell c = (y-1)*GX + x with x fastest (identical to the conformance exemplars' 'i = (row-1)*3 + col'), decoded in-expression as x-1 = (c-1) - GX*floor((c-1)/GX) and y-1 = floor((c-1)/GX); (2) the Camp Fire 19x21 LCC defaults are GONE -- a concrete grid is a consumer concern, so the INTEGER sizes are load-time metaparameters (GX/GY, bound at the import edge per esm-spec 9.7.6) and the REAL-VALUED origin/spacing are template parameters (x0/y0/dx/dy, bound per apply-site per 9.6.2), the same integers-are-metaparameters / reals-are-consumer-supplied split the grids/ libraries use. All values are affine in the indices and CONST at build time: the geometry folds into the regridder's overlap-matrix setup, off any per-step hot path. Vertex order is the CCW ring (W,S)->(E,S)->(E,N)->(W,N) = ring_verts 1..4, ring_coord 1=x 2=y; vertex x = x0 + (x-1)*dx + dx*[v in {2,3}], vertex y = y0 + (y-1)*dy + dy*[v in {3,4}]; centers are x0 + (x-1)*dx + dx/2 (= xmin + (x-1)*dx). Coordinates are in the grid's OWN frame -- lon-lat degrees for a geographic grid (feed conservative_overlap_A_ij_planar for the coordinate-plane model or compose with a reprojection for 'spherical'), projected metres for an LCC fire grid (compose with reprojection/lambert_conformal.esm to reach lon-lat). TWO cartesian ring grids in one document (cartesian-to-cartesian regridding): import this file TWICE under distinct esm-spec 9.7.7 namespace prefixes -- e.g. prefix `src` (GX/GY of the source) and prefix `tgt` (GX/GY of the target). Renaming is transitive, so the ring_* index sets arrive as src.ring_cells / tgt.ring_cells (no diamond-rule collision -- differently-prefixed edges register distinctly per esm-spec 9.7.7), and the match-less templates namespace to src.cartesian_cell_rings_poly / tgt.cartesian_cell_rings_poly, applied by their DOTTED names with the {from: ring_cells} references inside their aggregate bodies renamed to each edge's own axis. Both constructed ring stacks then feed one conservative_overlap import as src_poly/tgt_poly, so grid-spec -> rings -> regrid runs with no out-of-document precomputation on either side (tests/conformance/regridding/cartesian_rings_regrid_cart2cart_3x3_to_2x2/ is the end-to-end exemplar; the single-import sibling cartesian_rings_regrid_3x3_to_2x2/ predates renaming and still supplies its source rings as const data). A bare (un-prefixed) second import WOULD collide on the merged ring_* index sets (esm-spec 9.7.5/9.7.6 diamond rule) -- the prefix is what makes the two instances coexist. Verification: the affine formulas are pinned exactly (0.0 error) by the conformance case tests/conformance/regridding/cartesian_rings_regrid_3x3_to_2x2/, whose consumer compares the emitted rings/centers against inline expected constants and then feeds the SAME emitted rings through the conservative regridder end-to-end (grid-spec -> rings -> broad phase -> narrow phase -> apply), gated on the exact conservation / partition-of-unity invariants.

<details><summary>Templates (3)</summary>

**`cartesian_cell_rings_lon`**

Per-cell representative x-coordinate (cell centre) over the flattened row-major enumeration: lon[c] = x0 + ((c-1) mod GX)*dx + dx/2, with the mod spelled (c-1) - GX*floor((c-1)/GX) from the evaluable core (no bespoke mod op). This is the regridder's broad-phase binning coordinate (conservative_overlap_{src,tgt}_bin) and any point-based consumer's cell centre. CONST at build time.

Parameters: `x0`, `dx`

<div class="math-block">sum[c] (x0 + (c - 1 - GX * floor((c - 1) / GX)) * dx + dx / 2) where {c in ring_cells}</div>

**`cartesian_cell_rings_lat`**

Per-cell representative y-coordinate (cell centre) over the flattened row-major enumeration: lat[c] = y0 + floor((c-1)/GX)*dy + dy/2. CONST at build time.

Parameters: `y0`, `dy`

<div class="math-block">sum[c] (y0 + floor((c - 1) / GX) * dy + dy / 2) where {c in ring_cells}</div>

**`cartesian_cell_rings_poly`**

Per-cell closed corner ring over the flattened row-major enumeration, CCW (W,S)->(E,S)->(E,N)->(W,N) = ring_verts 1..4, ring_coord 1=x 2=y: vertex x = x0 + ((c-1) mod GX)*dx + dx*[v in {2,3}], vertex y = y0 + floor((c-1)/GX)*dy + dy*[v in {3,4}] (the boolean bracket is the sum of == comparisons, 0 or 1 in the evaluable core). Shape [ring_cells, ring_verts, ring_coord] -- exactly the [cells, overlap_cell_verts, overlap_coord] ring stack conservative_overlap_A_ij_planar / conservative_overlap_A_ij_gated_planar take as src_poly/tgt_poly (bind the consumer's NSRC/NTGT to GX*GY at the import edges, and NVERT and the narrow phase's src_nv/tgt_nv to one consumer metaparameter equal to 4). CONST at build time; folds into the overlap-matrix setup.

Parameters: `x0`, `y0`, `dx`, `dy`

<div class="math-block">sum[c, v, k] (ifelse(k == 1, x0 + (c - 1 - GX * floor((c - 1) / GX)) * dx + dx * ((v == 2) + (v == 3)), y0 + floor((c - 1) / GX) * dy + dy * ((v == 3) + (v == 4)))) where {c in ring_cells, k in ring_coord, v in ring_verts}</div>

</details>

**References**

- EarthSciAST esm-spec.md 8.6/8.6.1 (regridding is a coupling expression; cell rings as ordinary [verts, coord] coordinate arrays) and 9.7 (template libraries and metaparameters).

### `cartesian_nearest_regrid`

<p class="rule-meta"><code>method:nearest</code></p>

Source: [`regridding/cartesian_nearest_regrid.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/regridding/cartesian_nearest_regrid.esm)

Nearest-neighbour (containing-cell) regridding from a uniform cartesian SOURCE grid onto an arbitrary target grid, factored as a template library (esm-spec 9.7) so a coupling can import it by reference instead of restating the gather inline — the companion to regridding/conservative_overlap.esm. Conservative overlap area-AVERAGES and is correct for continuous fields (elevation, temperature, wind); this GATHERS the source value unchanged and is the correct choice for CATEGORICAL fields (integer class codes — fuel model, land cover, soil type), where averaging codes 5 and 91 into 48 is physically meaningless. It is also a pure index map (no polygon clipping), so it is far cheaper than the conservative narrow phase. One template: cartesian_nearest_regrid_apply. See its description for the exact floor-quantised source-cell lookup and the determinism/tie-breaking contract.

<details><summary>Templates (1)</summary>

**`cartesian_nearest_regrid_apply`**

NEAREST-NEIGHBOUR (containing-cell) regrid onto a target grid: F_tgt[j] = F_src[c(j)], where c(j) is the flattened index of the SOURCE cell whose extent contains target cell j's centre. Unlike conservative_overlap (area-weighted AVERAGE — correct for continuous fields), this is a verbatim GATHER — the source value is copied unchanged, so CATEGORICAL fields (integer class codes: land cover, fuel model, soil type) keep exact valid codes with no meaningless between-class averages. It is also far cheaper: a pure index map, no polygon intersection. The source cell is located by floor-quantising the target centre into the source's uniform cartesian lattice: sx = clamp(floor((tgt_lon[j]-src_x0)/src_dx)+1, 1, GX_SRC), sy likewise with src_y0/src_dy/GY_SRC, and c = (sy-1)*GX_SRC + sx over the x-fastest row-major enumeration cartesian_cell_rings uses. floor makes the pick DETERMINISTIC (a centre on a cell boundary resolves to the east/north cell — no argmin tie), so the result is byte-identical across bindings. tgt_lon/tgt_lat are the per-target-cell centres (bind cartesian_cell_rings_lon/lat); src_x0/src_dx/src_y0/src_dy + the GX_SRC/GY_SRC metaparameters describe the SOURCE lattice. F_src is shaped [src_cells]; the clamp keeps a target cell outside the source extent on the nearest edge cell rather than out of bounds. For a non-uniform or unstructured source, a centroid-argmin nearest (min-distance value-at-argmin, the lev_min idiom) is the general form; this template is the uniform-cartesian case.

Parameters: `F_src`, `tgt_lon`, `tgt_lat`, `src_x0`, `src_dx`, `src_y0`, `src_dy`

<div class="math-block">sum[j] (F_src[(max(1, min(GY_SRC, (floor((tgt_lat[j] - src_y0) / src_dy) + 1))) - 1) * GX_SRC + max(1, min(GX_SRC, (floor((tgt_lon[j] - src_x0) / src_dx) + 1)))]) where {j in tgt_cells}</div>

</details>

**References**

- EarthSciAST esm-spec.md 8.6 (regridding is a coupling expression) / 9.7 (expression-template libraries); the categorical-regrid gotcha (gt-categorical) noted in components/earthsci_data/landfire.esm — the upstream EarthSciData.jl LANDFIRE loader used nearest-neighbour (BSpline(Constant())) for exactly this reason.

### `conservative_overlap`

<p class="rule-meta"><code>method:conservative_overlap</code></p>

Source: [`regridding/conservative_overlap.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/regridding/conservative_overlap.esm)

Conservative overlap-area regridding, factored as a template library (esm-spec 9.7) so a coupling can invoke it instead of restating the aggregates inline. This is the library form of exactly the computation in the ESS worked fixtures tests/valid/geometry/conservative_regrid_overlap_join.esm and tests/coupling/cross_domain_coupling.esm (esm-spec 8.6 / 8.6.1): F_tgt[j] = (1/A_j) * sum_i A_ij * F_src[i], with A_ij = the overlap area of src_i and tgt_j (the esm:std/geometry.esm templates overlap_area_planar / overlap_area_spherical, esm 2.0.0) and A_j = sum_i A_ij. The match-less templates factor the stages: (1) conservative_overlap_A_ij_planar / conservative_overlap_A_ij_spherical -- the NARROW PHASE as a dense evaluable faq whose element is the scalar overlap-area template of esm:std/geometry.esm (esm-spec 8.6.1: the template returns the scalar area with its convex clip written out in its body, so the narrow phase is a dense tensor contraction with no ragged clip-ring intermediate); (2) conservative_overlap_A_j -- the ROW-NORMALIZATION denominator, the group-by-j row-sum of the SAME computed areas with the sliver filter A_ij > atol (CONFORMANCE_SPEC 5.8.2), which is what makes partition-of-unity sum_i W_ij = 1 EXACT BY CONSTRUCTION (5.8.3: the denominator is the row-sum of the computed areas, not the true target-cell area) -- note the denominator applies the SAME filter as the numerator, so the surviving sets coincide and the invariant is exact; (3) conservative_overlap_apply -- the weighted-sum APPLY + NORMALIZE, F_tgt[j] = sum_i A_ij*F_src[i] / A_j[j] over the same filtered set (the sparse mat-vec of ConservativeRegridding.jl, W_ij = A_ij/A_j folded in). THE BROAD PHASE IS THE SPATIAL-OVERLAP JOIN GATE (projection-pushdown Phase 2a): (4) the _gated narrow-phase/assembly variants (conservative_overlap_A_ij_gated_planar / _A_j_gated / _W_ij_gated / _apply_gated) -- the SAME aggregates with the candidate join folded in as a join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]} gate. The engine builds the candidate (i,j) set ONCE at build time from the cell-ring ENVELOPES -- each ring's axis-aligned bounding box, derived by the engine from the SAME polygon stacks the narrow phase clips (no hand-authored binning coordinate, no bin key, no bin-width parameter), inflated outward by a small eps -- byte-identically in every binding. A product term is contributed only where the source and target envelopes intersect (src_poly is the QUERY side, tgt_poly the INDEXED cell side), so the contraction runs over O(#candidates) instead of O(NSRC*NTGT). DENSE OR GATED IS THE CONSUMER'S CHOICE, value-identically: with the sliver filter present the two forms agree EXACTLY -- a non-candidate pair has disjoint bounding boxes hence zero true overlap area, which the dense path computes and the filter drops while the gated path never visits it; identical surviving sets, identical sums, and adding exact-zero terms to a nonnegative-area sum is FP-exact. The broad phase is ADMISSIBLE BY CONSTRUCTION: envelope candidacy never misses a genuine positive-area overlap for any eps >= 0 (two cells whose interiors overlap have overlapping bounding boxes), so no consumer bin-width admissibility obligation remains -- disjoint cell clusters are still pruned automatically (their bounding boxes do not intersect). The gated A_j applies the SAME join gate AND the SAME filter as the gated numerator, so the surviving sets coincide and partition-of-unity stays exact by construction -- the denominator broad-phases like the numerator. Cell polygons ARRIVE AS TEMPLATE PARAMETERS (closed vertex rings shaped [cells, overlap_cell_verts, overlap_coord], CCW, ring vertex v = [lon, lat]), so the templates are grid-agnostic: any pair of grids that can present its cells as fixed-vertex-count rings flows through unchanged (regridding/cartesian_cell_rings.esm derives such rings from a uniform cartesian grid spec; grids/mpas/grid.esm's mpas_cell_rings_deg derives them from MPAS mesh vertex data). Grid sizes are load-time metaparameters (NSRC, NTGT source/target cell counts; NVERT ring length, default quads), bound at the import edge; the defaults match the 3x3 -> 2x2 conformance exemplar. The MANIFOLD is chosen by WHICH narrow-phase template a consumer applies (esm 2.0.0 removed the polygon_intersection_area leaf and its manifold field; esm-spec 4.2 'Removed ops'): conservative_overlap_A_ij_planar clips in the source/target coordinate plane (areas in the squared coordinate unit; exact for the lon-lat-rectangle coarsening exemplar and for any projected patch), and conservative_overlap_A_ij_spherical clips with great-circle edges on the unit sphere (areas in steradians; the correct model for lat-lon-on-sphere regridding, with the CONFORMANCE_SPEC 5.8.4 great-circle-edge caveat for parallel edges). All three narrow-phase templates are the same pairwise skeleton (conservative_overlap_pairs, or conservative_overlap_pairs_gated with the broad-phase join) around a per-manifold element (conservative_overlap_area_planar_ij / conservative_overlap_area_spherical_ij). The spherical template reads [cells, verts, overlap_xyz] UNIT-VECTOR stacks, which the consumer builds from its lon-lat rings with conservative_overlap_unit_rings (a faq of esm:std/geometry.esm lonlat_to_unit over the ring's own cell and vertex axes); the planar templates read the lon-lat (or projected x-y) rings directly. Both take the source and target ring slot counts (src_nv, tgt_nv) explicitly, since a ring stack's vertex axis need not be overlap_cell_verts (an MPAS source ring is padded to maxEdges); a consumer binds them to the same metaparameter that sizes that axis, so the count is stated once. The broad phase is planar-envelope only, so the gated variants are planar. Everything downstream of A_ij (A_j, W_ij, apply and their gated forms) is manifold-agnostic. Conformance contract for consumers: candidate set byte-identical (when the broad phase is used), conservation and partition-of-unity as the primary exact gates, per-pair areas/weights toleranced with the sliver floor (CONFORMANCE_SPEC 5.8.1-5.8.3), and dense/gated value-identity gated exactly where both paths are computed. Exemplars: tests/conformance/regridding/conservative_overlap_3x3_to_2x2/ (dense, planar), conservative_overlap_sphere_fan_4_to_2/ (dense, spherical), conservative_overlap_gated_3x3_to_2x2/ (gated vs dense value-identity on the exemplar geometry), conservative_overlap_gated_clusters_4_to_2/ (gated with genuine candidate pruning across disjoint clusters), cartesian_rings_regrid_3x3_to_2x2/ (grid-spec -> rings -> broad phase -> narrow phase -> apply, end-to-end declarative), and mpas_l0_to_octants_sphere/ (MPAS mesh-derived rings on the sphere).

<details><summary>Templates (14)</summary>

**`conservative_overlap_pairs`**

DENSE PAIRWISE SKELETON (narrow-phase scaffold shared by the planar and spherical A_ij templates): the [src_cells, tgt_cells] matrix whose (i, j) entry is the expression bound to `area`, evaluated with loop symbols i (source cell) and j (target cell) in scope -- bind `area` to conservative_overlap_area_planar_ij or conservative_overlap_area_spherical_ij, which read exactly those symbols. src_ring and tgt_ring are the ring stacks the element gathers from (the faq's args). Every source/target pair is visited; see conservative_overlap_pairs_gated for the broad-phase-pruned form.

Parameters: `src_ring`, `tgt_ring`, `area`

<div class="math-block">sum[i, j] (area) where {i in src_cells, j in tgt_cells}</div>

**`conservative_overlap_pairs_gated`**

GATED PAIRWISE SKELETON: conservative_overlap_pairs with the spatial-overlap broad phase folded in as a join.overlap{src_env:[src_ring], tgt_env:[tgt_ring]} gate, so `area` is evaluated only for the candidate (i, j) pairs whose ring envelopes (axis-aligned bounding boxes inflated by eps) intersect; every other entry is the additive identity 0. The envelopes are planar, so src_ring / tgt_ring are the lon-lat (or projected x-y) rings, never unit-vector stacks.

Parameters: `src_ring`, `tgt_ring`, `area`

<div class="math-block">sum[i, j] (area) where {i in src_cells, j in tgt_cells} join(overlap(src=[src_poly], tgt=[tgt_poly], eps=1.0e-9))</div>

**`conservative_overlap_area_planar_ij`**

NARROW-PHASE ELEMENT, PLANAR: the esm:std/geometry.esm overlap_area_planar of source cell i of src_poly and target cell j of tgt_poly, i.e. one entry of A_ij. It reads the loop symbols i and j of the enclosing conservative_overlap_pairs / conservative_overlap_pairs_gated faq, so it is only meaningful bound to that skeleton's `area`. src_poly and tgt_poly are [cells, verts, overlap_coord] lon-lat (or x-y) stacks; src_nv and tgt_nv are the vertex slot counts of their rings.

Parameters: `src_poly`, `tgt_poly`, `src_nv`, `tgt_nv`

<div class="math-block">overlap_area_planar&lt;cell_a=i, cell_b=j, nv_a=src_nv, nv_b=tgt_nv, ring_a=src_poly, ring_b=tgt_poly&gt;</div>

**`conservative_overlap_area_spherical_ij`**

NARROW-PHASE ELEMENT, SPHERICAL: the esm:std/geometry.esm overlap_area_spherical of source cell i of src_unit and target cell j of tgt_unit, i.e. one entry of A_ij. It reads the loop symbols i and j of the enclosing conservative_overlap_pairs / conservative_overlap_pairs_gated faq, so it is only meaningful bound to that skeleton's `area`. src_unit and tgt_unit are [cells, verts, overlap_xyz] unit-vector stacks; src_nv and tgt_nv are the vertex slot counts of their rings.

Parameters: `src_unit`, `tgt_unit`, `src_nv`, `tgt_nv`

<div class="math-block">overlap_area_spherical&lt;cell_a=i, cell_b=j, nv_a=src_nv, nv_b=tgt_nv, ring_a=src_unit, ring_b=tgt_unit&gt;</div>

**`conservative_overlap_unit_rings`**

UNIT-VECTOR RING STACK (input adapter for conservative_overlap_A_ij_spherical): the [ring_cells, ring_verts, overlap_xyz] stack of unit vectors of the lon-lat-degree ring stack lonlat_ring, each entry the esm:std/geometry.esm lonlat_to_unit of one vertex component, so the trig is evaluated once per vertex rather than once per cell pair. ring_cells and ring_verts NAME the index sets lonlat_ring is shaped over (e.g. src_cells and overlap_cell_verts, or an MPAS mesh's cells and maxEdges), so one template serves source and target stacks of any grid.

Parameters: `lonlat_ring`, `ring_cells`, `ring_verts`

<div class="math-block">sum[c, v, d] (lonlat_to_unit&lt;cell=c, comp=d, ring=lonlat_ring, vert=v&gt;) where {c in ring_cells, d in overlap_xyz, v in ring_verts}</div>

**`conservative_overlap_A_ij_planar`**

NARROW PHASE, PLANAR (overlap-weight template): the raw overlap-area matrix A_ij[i,j] = overlap_area_planar(src_poly[i], tgt_poly[j]) (esm:std/geometry.esm) as a dense evaluable faq (esm-spec 8.6.1): conservative_overlap_pairs over the conservative_overlap_area_planar_ij element, clipping with straight edges in the coordinate plane; areas in the squared coordinate unit. src_poly is shaped [src_cells, verts, overlap_coord], tgt_poly [tgt_cells, verts, overlap_coord]; each cell is a closed convex vertex ring with vertex v = [lon, lat] (or [x, y] on a projected patch), either orientation. src_nv is the slot count of each source ring (its vertex axis length), tgt_nv of each target ring; a ring shorter than its slot count repeats its last vertex, which contributes nothing. Bind src_nv / tgt_nv to the SAME metaparameter that sizes the ring stack's vertex axis, never to a separate literal: a count below the axis length silently clips each cell as a polygon with fewer vertices, and because the A_j denominator row-sums the same (wrong) areas, partition-of-unity and conservation still hold -- no invariant gate would catch it. Entries for disjoint pairs are exactly 0; near-tangent slivers are dropped downstream by the atol filter of the A_j / apply templates, never here -- the raw areas are what the CONFORMANCE_SPEC 5.8.2 tolerance gates compare.

Parameters: `src_poly`, `tgt_poly`, `src_nv`, `tgt_nv`

<div class="math-block">conservative_overlap_pairs&lt;area=conservative_overlap_area_planar_ij&lt;src_nv=src_nv, src_poly=src_poly, tgt_nv=tgt_nv, tgt_poly=tgt_poly&gt;, src_ring=src_poly, tgt_ring=tgt_poly&gt;</div>

**`conservative_overlap_A_ij_spherical`**

NARROW PHASE, SPHERICAL (overlap-weight template): the raw overlap-area matrix A_ij[i,j] = overlap_area_spherical(src_unit[i], tgt_unit[j]) (esm:std/geometry.esm) as a dense evaluable faq (esm-spec 8.6.1): conservative_overlap_pairs over the conservative_overlap_area_spherical_ij element, clipping with great-circle edges on the unit sphere; areas in steradians (multiply by R^2 for a physical area). src_unit is shaped [src_cells, verts, overlap_xyz], tgt_unit [tgt_cells, verts, overlap_xyz]: UNIT-VECTOR ring stacks, built by the consumer from its lon-lat-degree rings with conservative_overlap_unit_rings, so the conversion is computed once per vertex rather than once per pair. src_nv is the slot count of each source ring (its vertex axis length), tgt_nv of each target ring; a ring shorter than its slot count repeats its last vertex, which contributes nothing. Bind src_nv / tgt_nv to the SAME metaparameter that sizes the ring stack's vertex axis, never to a separate literal: a count below the axis length silently clips each cell as a polygon with fewer vertices, and because the A_j denominator row-sums the same (wrong) areas, partition-of-unity and conservation still hold -- no invariant gate would catch it. Agrees with the removed polygon_intersection_area(manifold spherical) within the CONFORMANCE_SPEC 5.8.2 tolerance, not bit for bit. Entries for disjoint pairs are exactly 0; slivers are dropped downstream by the atol filter.

Parameters: `src_unit`, `tgt_unit`, `src_nv`, `tgt_nv`

<div class="math-block">conservative_overlap_pairs&lt;area=conservative_overlap_area_spherical_ij&lt;src_nv=src_nv, src_unit=src_unit, tgt_nv=tgt_nv, tgt_unit=tgt_unit&gt;, src_ring=src_unit, tgt_ring=tgt_unit&gt;</div>

**`conservative_overlap_A_j`**

ROW NORMALIZATION (partition-of-unity template): the per-target-cell denominator A_j[j] = sum_i A_ij[i,j] over the surviving-overlap set (filter A_ij > atol -- the CONFORMANCE_SPEC 5.8.2 sliver floor; sub-atol areas are equal-to-zero and contribute the additive identity). Because the apply template divides by THIS row-sum of the SAME computed areas under the SAME filter, the weights W_ij = A_ij/A_j satisfy sum_i W_ij = 1 exactly by construction (5.8.3), regardless of edge-model or clipping error -- the exact invariant the conformance gate asserts. This is ConservativeRegridding.jl's dst_areas (row-sums of the intersections matrix), equation (3) of the ESS overlap-join fixture.

Parameters: `A_ij`, `atol`

<div class="math-block">sum[j] (A_ij[i, j]) where {i in src_cells, j in tgt_cells} if A_ij[i, j] &gt; atol</div>

**`conservative_overlap_W_ij`**

NORMALIZED WEIGHTS (observability template): the explicit weight matrix W_ij[i,j] = A_ij[i,j] / A_j[j] over the surviving-overlap set, under the SAME filter A_ij > atol as the A_j and apply templates (a filtered-out sliver's weight is exactly 0, the aggregate's additive identity). This is ConservativeRegridding.jl's W matrix that conservative_overlap_apply folds into its sparse mat-vec -- factored out as its own template so a conformance runner can gate the per-pair weights through the official build-time pathway (CONFORMANCE_SPEC 5.8.1-5.8.2; the EarthSciAST.jl BuildInspection setup-array surface). Because the denominator is the row-sum of the SAME computed areas under the SAME filter, sum_i W_ij = 1 exactly for every target cell with a surviving overlap (5.8.3). Bind A_j to the materialized result of conservative_overlap_A_j. W_ij is build-once diagnostic data: the apply template never reads it, so declaring it adds no work to the ODE hot path.

Parameters: `A_ij`, `A_j`, `atol`

<div class="math-block">sum[i, j] (A_ij[i, j] / A_j[j]) where {i in src_cells, j in tgt_cells} if A_ij[i, j] &gt; atol</div>

**`conservative_overlap_apply`**

WEIGHTED-SUM APPLY (+ fold-in normalize): F_tgt[j] = sum_i A_ij[i,j]*F_src[i] / A_j[j] over the same filtered surviving set as the A_j template -- i.e. F_tgt = sum_i W_ij*F_src[i] with W_ij = A_ij/A_j, the sparse mat-vec apply of ConservativeRegridding.jl (equations (4)+(5) of the ESS overlap-join fixture, the heat_flux observed of cross_domain_coupling.esm). Bind A_j to the materialized result of conservative_overlap_A_j so the denominator is the row-sum of the same computed areas: that is what makes partition-of-unity exact and conservation sum_j A_j*F_tgt[j] = sum_i A_i*F_src[i] hold wherever the target tiles the source domain (CONFORMANCE_SPEC 5.8.3). For a time-varying F_src this template is the apply-many half of the build-once/apply-many split; A_ij and A_j are build-once.

Parameters: `A_ij`, `A_j`, `F_src`, `atol`

<div class="math-block">sum[j] (A_ij[i, j] * F_src[i] / A_j[j]) where {i in src_cells, j in tgt_cells} if A_ij[i, j] &gt; atol</div>

**`conservative_overlap_A_ij_gated_planar`**

CANDIDATE-GATED NARROW PHASE: the overlap-area matrix of conservative_overlap_A_ij_planar with the SPATIAL-OVERLAP broad-phase candidate join folded in -- join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]} (projection-pushdown Phase 2a) admits a product term only where the source cell's and target cell's envelopes intersect (each ring factor's axis-aligned bounding box, computed once at build time from the SAME polygon rings the narrow phase clips; inflated outward by eps), so one polygon clip runs per CANDIDATE pair instead of per (i,j) pair: O(#candidates) build cost. src_poly is the QUERY side (src_env), tgt_poly the INDEXED cell side (tgt_env); the candidate (i,j) set is built once and consulted per contracted tuple. Non-candidate entries are the aggregate's additive identity, exactly 0. VALUE-IDENTICAL to the dense template: a pruned pair has zero true overlap area (disjoint bounding boxes => disjoint cells), which the dense path computes as 0 and the downstream sliver filter drops anyway -- identical surviving sets, identical sums. eps is a small conservative slack (the broad phase never misses a genuine positive-area overlap for eps >= 0). PLANAR only: the envelope broad phase is a planar bounding-box test, and the narrow phase is overlap_area_planar exactly as in conservative_overlap_A_ij_planar (same src_nv / tgt_nv contract) -- the two differ only in the skeleton, conservative_overlap_pairs_gated instead of conservative_overlap_pairs.

Parameters: `src_poly`, `tgt_poly`, `src_nv`, `tgt_nv`

<div class="math-block">conservative_overlap_pairs_gated&lt;area=conservative_overlap_area_planar_ij&lt;src_nv=src_nv, src_poly=src_poly, tgt_nv=tgt_nv, tgt_poly=tgt_poly&gt;, src_ring=src_poly, tgt_ring=tgt_poly&gt;</div>

**`conservative_overlap_A_j_gated`**

CANDIDATE-GATED ROW NORMALIZATION: A_j[j] = sum_i A_ij[i,j] over the CANDIDATE pairs (join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]}) that survive the sliver filter A_ij > atol. THE DENOMINATOR BROAD-PHASES LIKE THE NUMERATOR: this template applies the SAME spatial-overlap join gate and the SAME filter as conservative_overlap_A_ij_gated_planar / the apply, so the surviving-overlap sets coincide and partition-of-unity sum_i W_ij = 1 stays EXACT BY CONSTRUCTION (CONFORMANCE_SPEC 5.8.3) -- gating only the numerator, or only the denominator, would break the invariant wherever a spurious sliver or a pruned pair appears on one side only. Value-identical to conservative_overlap_A_j on a gated A_ij: the pruned entries have zero overlap area and are sub-atol, so both forms sum the same surviving terms. Bind src_poly/tgt_poly to the SAME ring stacks passed to conservative_overlap_A_ij_gated_planar (the envelope broad phase reads their bounding boxes at build time).

Parameters: `A_ij`, `src_poly`, `tgt_poly`, `atol`

<div class="math-block">sum[j] (A_ij[i, j]) where {i in src_cells, j in tgt_cells} join(overlap(src=[src_poly], tgt=[tgt_poly], eps=1.0e-9)) if A_ij[i, j] &gt; atol</div>

**`conservative_overlap_W_ij_gated`**

CANDIDATE-GATED NORMALIZED WEIGHTS: W_ij[i,j] = A_ij[i,j] / A_j[j] over the candidate pairs surviving the SAME spatial-overlap join gate (join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]}) and SAME filter as the gated A_j (a pruned or filtered pair's weight is exactly 0, the additive identity). Bind A_j to the materialized result of conservative_overlap_A_j_gated, and src_poly/tgt_poly to the same ring stacks as the numerator. Build-once diagnostic data exactly like the dense W_ij template: the conformance runner reads it from the official build-time setup-array surface; the apply never does.

Parameters: `A_ij`, `A_j`, `src_poly`, `tgt_poly`, `atol`

<div class="math-block">sum[i, j] (A_ij[i, j] / A_j[j]) where {i in src_cells, j in tgt_cells} join(overlap(src=[src_poly], tgt=[tgt_poly], eps=1.0e-9)) if A_ij[i, j] &gt; atol</div>

**`conservative_overlap_apply_gated`**

CANDIDATE-GATED WEIGHTED-SUM APPLY (+ fold-in normalize): F_tgt[j] = sum_i A_ij[i,j]*F_src[i] / A_j[j] over the candidate pairs surviving the SAME spatial-overlap join gate (join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]}) and SAME filter as the gated A_j -- the O(#overlaps) production spelling of the ESS worked fixture's equations (4)+(5), the sparse mat-vec restricted to the candidate set. Bind A_j to the materialized result of conservative_overlap_A_j_gated so the denominator is the row-sum of the same computed areas over the same surviving set: partition-of-unity exact by construction, conservation wherever the target tiles the source domain (CONFORMANCE_SPEC 5.8.3). Value-identical to conservative_overlap_apply on the same (gated) A_ij/A_j. A_ij and A_j are build-once; when F_src is build-time-known the whole apply materializes at setup (the candidate set + narrow filter select the same surviving terms as the dense apply). Bind src_poly/tgt_poly to the same ring stacks as the numerator.

Parameters: `A_ij`, `A_j`, `F_src`, `src_poly`, `tgt_poly`, `atol`

<div class="math-block">sum[j] (A_ij[i, j] * F_src[i] / A_j[j]) where {i in src_cells, j in tgt_cells} join(overlap(src=[src_poly], tgt=[tgt_poly], eps=1.0e-9)) if A_ij[i, j] &gt; atol</div>

</details>

**References**

- Jones, P. W. (1999). First- and second-order conservative remapping schemes for grids in spherical coordinates. Monthly Weather Review 127(9), 2204-2210. [`doi:10.1175/1520-0493(1999)127<2204:FASOCR>2.0.CO;2`](https://doi.org/10.1175/1520-0493%281999%29127%3C2204%3AFASOCR%3E2.0.CO%3B2)
- EarthSciAST esm-spec.md 8.6/8.6.1 (regridding is a coupling expression; the overlap-area templates of esm:std/geometry.esm) and CONFORMANCE_SPEC.md 5.8 (geometry tolerance contract).

