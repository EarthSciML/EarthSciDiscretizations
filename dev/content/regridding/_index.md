---
title: "Regridding"
description: "Conservative and interpolating remap expressions (coupling transforms). Each entry is a pure ESM template library — declarative expressions, no kernels."
---

Conservative and interpolating remap expressions (coupling transforms). Each entry is a pure ESM template library — declarative expressions, no kernels.

The conservative overlap regridder is **end-to-end declarative**: a consumer constructs source and target cell rings, derives the broad-phase bin keys, and runs the candidate-gated overlap → normalize → apply, all from `.esm` templates with no host glue.

- **Cell-ring constructors** turn a grid spec into the `[cells, verts, coord]` vertex rings the regridder consumes — `cartesian_cell_rings.esm` from a uniform cartesian spec, and `grids/mpas/grid.esm`'s `mpas_cell_rings_deg` from MPAS mesh vertex data (the widened contract).
- **Broad phase in-library**: representative binning coordinates derived from the rings themselves plus integer `skolem` bin keys, so the candidate set is byte-identical across bindings.
- **Gated == dense** value-identically under an admissible binning (the gated denominator broad-phases like the numerator, so partition-of-unity stays exact).
- **One file, every manifold**: `planar` or `spherical` is a per-invocation template parameter. Spherical clip runs on Julia + Rust today; the Python spherical case is gated on the optional `spherely` dependency (planar cases run on Julia + Python + Rust), and the two rewrite-only ports (Go, TypeScript) are scope-excluded — recorded in each case manifest, never shimmed.

End-to-end exemplars: `tests/conformance/regridding/cartesian_rings_regrid_gated_3x3_to_2x2` (grid spec → rings → broad phase → gated overlap → regrid) and `mpas_l0_to_octants_sphere` (MPAS mesh-derived rings on the sphere).

### `cartesian_cell_rings`

<p class="rule-meta"><code>method:cartesian_cell_rings</code></p>

Source: [`regridding/cartesian_cell_rings.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/regridding/cartesian_cell_rings.esm)

Declarative cell GEOMETRY for a uniform cartesian grid, in regridder-ready form: given the grid spec -- per-axis origin x0/y0 (the west-most / south-most cell-EDGE coordinate), spacing dx/dy, and sizes GX/GY -- the templates emit the per-cell representative centers (cartesian_cell_rings_lon / _lat, shape [ring_cells]) and the per-cell four-corner closed CCW rings (cartesian_cell_rings_poly, shape [ring_cells, ring_verts, ring_coord]) that regridding/conservative_overlap.esm consumes as its src_poly/tgt_poly (narrow phase) and its binning coordinates (broad phase). This is the 0.8.0 template-library rewrite of the archived model-form file (archive/regridding/cartesian_target_cells.esm), with two deliberate changes: (1) the cell axis is FLATTENED to the single ring_cells = GX*GY enumeration the regridder's src_cells/tgt_cells axes expect -- row-major from the southwest, cell c = (y-1)*GX + x with x fastest (identical to the conformance exemplars' 'i = (row-1)*3 + col'), decoded in-expression as x-1 = (c-1) - GX*floor((c-1)/GX) and y-1 = floor((c-1)/GX); (2) the Camp Fire 19x21 LCC defaults are GONE -- a concrete grid is a consumer concern, so the INTEGER sizes are load-time metaparameters (GX/GY, bound at the import edge per esm-spec 9.7.6) and the REAL-VALUED origin/spacing are template parameters (x0/y0/dx/dy, bound per apply-site per 9.6.2), the same integers-are-metaparameters / reals-are-consumer-supplied split the grids/ libraries use. All values are affine in the indices and CONST at build time: the geometry folds into the regridder's overlap-matrix setup, off any per-step hot path. Vertex order is the CCW ring (W,S)->(E,S)->(E,N)->(W,N) = ring_verts 1..4, ring_coord 1=x 2=y; vertex x = x0 + (x-1)*dx + dx*[v in {2,3}], vertex y = y0 + (y-1)*dy + dy*[v in {3,4}]; centers are x0 + (x-1)*dx + dx/2 (= xmin + (x-1)*dx). Coordinates are in the grid's OWN frame -- lon-lat degrees for a geographic grid (feed conservative_overlap at manifold 'planar' for the coordinate-plane model or compose with a reprojection for 'spherical'), projected metres for an LCC fire grid (compose with reprojection/lambert_conformal.esm to reach lon-lat). TWO cartesian ring grids in one document (cartesian-to-cartesian regridding): import this file TWICE under distinct esm-spec 9.7.7 namespace prefixes -- e.g. prefix `src` (GX/GY of the source) and prefix `tgt` (GX/GY of the target). Renaming is transitive, so the ring_* index sets arrive as src.ring_cells / tgt.ring_cells (no diamond-rule collision -- differently-prefixed edges register distinctly per esm-spec 9.7.7), and the match-less templates namespace to src.cartesian_cell_rings_poly / tgt.cartesian_cell_rings_poly, applied by their DOTTED names with the {from: ring_cells} references inside their aggregate bodies renamed to each edge's own axis. Both constructed ring stacks then feed one conservative_overlap import as src_poly/tgt_poly, so grid-spec -> rings -> regrid runs with no out-of-document precomputation on either side (tests/conformance/regridding/cartesian_rings_regrid_cart2cart_3x3_to_2x2/ is the end-to-end exemplar; the single-import sibling cartesian_rings_regrid_3x3_to_2x2/ predates renaming and still supplies its source rings as const data). A bare (un-prefixed) second import WOULD collide on the merged ring_* index sets (esm-spec 9.7.5/9.7.6 diamond rule) -- the prefix is what makes the two instances coexist. Verification: the affine formulas are pinned exactly (0.0 error) by the conformance case tests/conformance/regridding/cartesian_rings_regrid_3x3_to_2x2/, whose consumer compares the emitted rings/centers against inline expected constants and then feeds the SAME emitted rings through the conservative regridder end-to-end (grid-spec -> rings -> broad phase -> narrow phase -> apply), gated on the exact conservation / partition-of-unity invariants.

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

Per-cell closed corner ring over the flattened row-major enumeration, CCW (W,S)->(E,S)->(E,N)->(W,N) = ring_verts 1..4, ring_coord 1=x 2=y: vertex x = x0 + ((c-1) mod GX)*dx + dx*[v in {2,3}], vertex y = y0 + floor((c-1)/GX)*dy + dy*[v in {3,4}] (the boolean bracket is the sum of == comparisons, 0 or 1 in the evaluable core). Shape [ring_cells, ring_verts, ring_coord] -- exactly the [cells, overlap_cell_verts, overlap_coord] ring stack conservative_overlap_A_ij{,_gated} takes as src_poly/tgt_poly (bind the consumer's NSRC/NTGT to GX*GY and NVERT to 4 at the import edges). CONST at build time; folds into the overlap-matrix setup.

Parameters: `x0`, `y0`, `dx`, `dy`

<div class="math-block">sum[c, v, k] (ifelse(k == 1, x0 + (c - 1 - GX * floor((c - 1) / GX)) * dx + dx * ((v == 2) + (v == 3)), y0 + floor((c - 1) / GX) * dy + dy * ((v == 3) + (v == 4)))) where {c in ring_cells, k in ring_coord, v in ring_verts}</div>

</details>

**References**

- EarthSciAST esm-spec.md 8.6/8.6.1 (regridding is a coupling expression; cell rings as ordinary [verts, coord] coordinate arrays) and 9.7 (template libraries and metaparameters).

### `cartesian_nearest_regrid`

<p class="rule-meta"><code>method:nearest</code></p>

Source: [`regridding/cartesian_nearest_regrid.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/regridding/cartesian_nearest_regrid.esm)

Nearest-neighbour (containing-cell) regridding from a uniform cartesian SOURCE grid onto an arbitrary target grid, factored as a template library (esm-spec 9.7) so a coupling can import it by reference instead of restating the gather inline — the companion to regridding/conservative_overlap.esm. Conservative overlap area-AVERAGES and is correct for continuous fields (elevation, temperature, wind); this GATHERS the source value unchanged and is the correct choice for CATEGORICAL fields (integer class codes — fuel model, land cover, soil type), where averaging codes 5 and 91 into 48 is physically meaningless. It is also a pure index map (no polygon_intersection_area), so it is far cheaper than the conservative narrow phase. One template: cartesian_nearest_regrid_apply. See its description for the exact floor-quantised source-cell lookup and the determinism/tie-breaking contract.

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

Conservative overlap-area regridding, factored as a template library (esm-spec 9.7) so a coupling can invoke it instead of restating the aggregates inline. This is the library form of exactly the computation in the ESS worked fixtures tests/valid/geometry/conservative_regrid_overlap_join.esm and tests/coupling/cross_domain_coupling.esm (esm-spec 8.6 / 8.6.1): F_tgt[j] = (1/A_j) * sum_i A_ij * F_src[i], with A_ij = polygon_intersection_area(src_i, tgt_j) and A_j = sum_i A_ij. The match-less templates factor the stages: (1) conservative_overlap_A_ij -- the NARROW PHASE as a dense evaluable aggregate over the fused polygon_intersection_area leaf (esm-spec 8.6.1: the fused leaf is what makes the narrow phase a dense tensor contraction with no ragged clip-ring intermediate); (2) conservative_overlap_A_j -- the ROW-NORMALIZATION denominator, the group-by-j row-sum of the SAME computed areas with the sliver filter A_ij > atol (CONFORMANCE_SPEC 5.8.2), which is what makes partition-of-unity sum_i W_ij = 1 EXACT BY CONSTRUCTION (5.8.3: the denominator is the row-sum of the computed areas, not the true target-cell area) -- note the denominator applies the SAME filter as the numerator, so the surviving sets coincide and the invariant is exact; (3) conservative_overlap_apply -- the weighted-sum APPLY + NORMALIZE, F_tgt[j] = sum_i A_ij*F_src[i] / A_j[j] over the same filtered set (the sparse mat-vec of ConservativeRegridding.jl, W_ij = A_ij/A_j folded in). THE BROAD PHASE IS THE SPATIAL-OVERLAP JOIN GATE (projection-pushdown Phase 2a): (4) the _gated narrow-phase/assembly variants (conservative_overlap_A_ij_gated / _A_j_gated / _W_ij_gated / _apply_gated) -- the SAME aggregates with the candidate join folded in as a join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]} gate. The engine builds the candidate (i,j) set ONCE at build time from the cell-ring ENVELOPES -- each ring's axis-aligned bounding box, derived by the engine from the SAME polygon stacks the narrow phase clips (no hand-authored binning coordinate, no bin key, no bin-width parameter), inflated outward by a small eps -- via an STRtree (when the GeometryOps extension is loaded) or a dependency-free brute-force sweep, byte-identically. A product term is contributed only where the source and target envelopes intersect (src_poly is the QUERY side, tgt_poly the INDEXED cell side), so the contraction runs over O(#candidates) instead of O(NSRC*NTGT). DENSE OR GATED IS THE CONSUMER'S CHOICE, value-identically: with the sliver filter present the two forms agree EXACTLY -- a non-candidate pair has disjoint bounding boxes hence zero true overlap area, which the dense path computes and the filter drops while the gated path never visits it; identical surviving sets, identical sums, and adding exact-zero terms to a nonnegative-area sum is FP-exact. The broad phase is ADMISSIBLE BY CONSTRUCTION: envelope candidacy never misses a genuine positive-area overlap for any eps >= 0 (two cells whose interiors overlap have overlapping bounding boxes), so no consumer bin-width admissibility obligation remains -- disjoint cell clusters are still pruned automatically (their bounding boxes do not intersect). The gated A_j applies the SAME join gate AND the SAME filter as the gated numerator, so the surviving sets coincide and partition-of-unity stays exact by construction -- the denominator broad-phases like the numerator. Cell polygons ARRIVE AS TEMPLATE PARAMETERS (closed vertex rings shaped [cells, overlap_cell_verts, overlap_coord], CCW, ring vertex v = [lon, lat]), so the templates are grid-agnostic: any pair of grids that can present its cells as fixed-vertex-count rings flows through unchanged (regridding/cartesian_cell_rings.esm derives such rings from a uniform cartesian grid spec; grids/mpas/grid.esm's mpas_cell_rings_deg derives them from MPAS mesh vertex data). Grid sizes are load-time metaparameters (NSRC, NTGT source/target cell counts; NVERT ring length, default quads), bound at the import edge; the defaults match the 3x3 -> 2x2 conformance exemplar. The MANIFOLD stays a TEMPLATE PARAMETER of the A_ij templates: per the esm-spec 9.6.1 scalar-field substitution rule (a parameter name appearing as the string value of a scalar Expression-node field in a body is a substitution site, the mirror of the match-side scalar-field binding rule), the geometry kernel's manifold flag is bound at each invocation -- 'planar' (clipping in the source/target coordinate plane; exact for the lon-lat-rectangle coarsening exemplar and for any projected patch) or 'spherical' (great-circle-edge clipping on the unit sphere, areas in steradians; the correct model for lat-lon-on-sphere regridding, with the CONFORMANCE_SPEC 5.8.4 great-circle-edge caveat for parallel edges). The bound value MUST be a member of the closed set {planar, spherical, geodesic}: substitution is pure-syntactic and the enum is enforced on the expanded form (9.6.4, diagnostic geometry_manifold_invalid). One library file therefore serves every manifold; there is no separate spherical near-duplicate. Conformance contract for consumers: candidate set byte-identical (when the broad phase is used), conservation and partition-of-unity as the primary exact gates, per-pair areas/weights toleranced with the sliver floor (CONFORMANCE_SPEC 5.8.1-5.8.3), and dense/gated value-identity gated exactly where both paths are computed. Exemplars: tests/conformance/regridding/conservative_overlap_3x3_to_2x2/ (dense, manifold 'planar'), conservative_overlap_sphere_fan_4_to_2/ (dense, 'spherical'), conservative_overlap_gated_3x3_to_2x2/ (gated vs dense value-identity on the exemplar geometry), conservative_overlap_gated_clusters_4_to_2/ (gated with genuine candidate pruning across disjoint clusters), cartesian_rings_regrid_3x3_to_2x2/ (grid-spec -> rings -> broad phase -> narrow phase -> apply, end-to-end declarative), and mpas_l0_to_octants_sphere/ (MPAS mesh-derived rings on the sphere).

<details><summary>Templates (8)</summary>

**`conservative_overlap_A_ij`**

NARROW PHASE (overlap-weight template): the raw overlap-area matrix A_ij[i,j] = polygon_intersection_area(src_poly[i], tgt_poly[j]) as a dense evaluable aggregate (esm-spec 8.6.1). src_poly is shaped [src_cells, overlap_cell_verts, overlap_coord], tgt_poly [tgt_cells, overlap_cell_verts, overlap_coord]; each cell is a closed CCW vertex ring with ring vertex v = [lon, lat] (or [x, y] on a projected patch). The fused leaf hides the clip ring, so no ragged intermediate appears and the whole matrix is a dense tensor contraction. The MANIFOLD is the third template parameter, occupying the leaf's scalar manifold field as a substitution site (esm-spec 9.6.1): bind 'planar' for coordinate-plane clipping (areas in the squared coordinate unit) or 'spherical' for great-circle-edge clipping on the unit sphere (areas in steradians); the bound value must be a member of the closed set {planar, spherical, geodesic}, enforced post-expansion (9.6.4, geometry_manifold_invalid). Entries for disjoint pairs are exactly 0; near-tangent slivers are dropped downstream by the atol filter of the A_j / apply templates, never here -- the raw areas are what the CONFORMANCE_SPEC 5.8.2 tolerance gates compare.

Parameters: `src_poly`, `tgt_poly`, `manifold`

<div class="math-block">sum[i, j] (polygon_intersection_area(src_poly[i], tgt_poly[j], manifold=manifold)) where {i in src_cells, j in tgt_cells}</div>

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

**`conservative_overlap_A_ij_gated`**

CANDIDATE-GATED NARROW PHASE: the overlap-area matrix of conservative_overlap_A_ij with the SPATIAL-OVERLAP broad-phase candidate join folded in -- join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]} (projection-pushdown Phase 2a) admits a product term only where the source cell's and target cell's envelopes intersect (each ring factor's axis-aligned bounding box, computed once at build time from the SAME polygon rings the narrow phase clips; inflated outward by eps), so one polygon clip runs per CANDIDATE pair instead of per (i,j) pair: O(#candidates) build cost. src_poly is the QUERY side (src_env), tgt_poly the INDEXED cell side (tgt_env); the candidate (i,j) set is built once (STRtree when the GeometryOps extension is loaded, else brute-force) and consulted per contracted tuple. Non-candidate entries are the aggregate's additive identity, exactly 0. VALUE-IDENTICAL to the dense template: a pruned pair has zero true overlap area (disjoint bounding boxes => disjoint cells), which the dense path computes as 0 and the downstream sliver filter drops anyway -- identical surviving sets, identical sums. eps is a small conservative slack (the broad phase never misses a genuine positive-area overlap for eps >= 0). Manifold semantics identical to conservative_overlap_A_ij.

Parameters: `src_poly`, `tgt_poly`, `manifold`

<div class="math-block">sum[i, j] (polygon_intersection_area(src_poly[i], tgt_poly[j], manifold=manifold)) where {i in src_cells, j in tgt_cells} join(overlap(src=[src_poly], tgt=[tgt_poly], eps=1.0e-9))</div>

**`conservative_overlap_A_j_gated`**

CANDIDATE-GATED ROW NORMALIZATION: A_j[j] = sum_i A_ij[i,j] over the CANDIDATE pairs (join.overlap{src_env:[src_poly], tgt_env:[tgt_poly]}) that survive the sliver filter A_ij > atol. THE DENOMINATOR BROAD-PHASES LIKE THE NUMERATOR: this template applies the SAME spatial-overlap join gate and the SAME filter as conservative_overlap_A_ij_gated / the apply, so the surviving-overlap sets coincide and partition-of-unity sum_i W_ij = 1 stays EXACT BY CONSTRUCTION (CONFORMANCE_SPEC 5.8.3) -- gating only the numerator, or only the denominator, would break the invariant wherever a spurious sliver or a pruned pair appears on one side only. Value-identical to conservative_overlap_A_j on a gated A_ij: the pruned entries have zero overlap area and are sub-atol, so both forms sum the same surviving terms. Bind src_poly/tgt_poly to the SAME ring stacks passed to conservative_overlap_A_ij_gated (the envelope broad phase reads their bounding boxes at build time).

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
- EarthSciAST esm-spec.md 8.6/8.6.1 (regridding is a coupling expression; the polygon_intersection_area fused leaf) and CONFORMANCE_SPEC.md 5.8 (geometry tolerance contract).

