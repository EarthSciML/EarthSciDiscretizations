---
title: "latlon"
description: "Uniform longitude-latitude grid with consumer-supplied extent (global or regional)."
esd_kind: "grid"
families: ["finite_difference", "finite_volume", "latlon"]
ops: ["D", "D2", "laplacian"]
bcs: ["periodic", "zero_gradient"]
tags: ["esd:grid", "family:latlon", "ndim:2", "spacing:uniform", "axis:lon", "axis:lat", "esd:stencil", "family:finite_difference", "grid:latlon", "op:D2", "order:2", "stencil_width:3", "op:D", "family:finite_volume", "op:laplacian", "stencil_width:5", "order:1", "stencil_width:2", "esd:rule", "axes:lon,lat", "bc:zero_gradient", "spec:9.6.8", "bc:periodic"]
source: "grids/latlon/grid.esm"
---

Uniform longitude-latitude grid with consumer-supplied extent (global or regional). The `lon` and `lat` index sets are the iteration axes; NLON and NLAT (the point counts) are the esm-spec 9.7 metaparameters, bound at an import edge or the loader API. The real-valued geometry is NOT baked into this file: `lon0_deg`, `dlon_deg`, `lat0_deg`, `dlat_deg` (all in degrees) are free names in the coordinate/metric template bodies and in this grid's stencils/rules, resolving in the consuming component's scope at evaluation — the same keyed-factor contract grids/mpas uses for areaCell/dvEdge and grids/cartesian_uniform_1d uses for x0/dx. COORDINATE CONVENTIONS (pinned; every rule in this directory assumes them): longitude is CELL-BASED — NLON cells with centers lon_i = lon0_deg + (i - 1/2) dlon_deg, i = 1..NLON, so the axis covers [lon0_deg, lon0_deg + NLON dlon_deg) and closes into a periodic zonal circle exactly when NLON dlon_deg = 360; latitude is POINT-BASED and endpoint-inclusive — NLAT points lat_j = lat0_deg + (j - 1) dlat_deg, j = 1..NLAT, covering [lat0_deg, lat0_deg + (NLAT - 1) dlat_deg]. GLOBAL RECIPE (the default 144 x 91 sizes match the esm-spec 9.6.8 worked example): the consuming model defines lon0_deg = 0 (parameter), dlon_deg = 360/NLON, lat0_deg = -90 (parameter), dlat_deg = 180/(NLAT - 1) — the lat points then run pole to pole with lat_1 = -90 and lat_NLAT = +90. Spell the spacings as observed variables whose expressions divide by the METAPARAMETER NAME (e.g. dlon_deg = {op: /, args: [360, NLON]}), which esm-spec 9.7.6 substitutes as an integer literal at load, so a loader-API rebinding of NLON/NLAT — a convergence sweep — keeps the geometry consistent automatically. A model that binds NLON/NLAT at the import edge closes them before they reach the model's own scope and must spell the matching literals instead. REGIONAL RECIPE: any other origin/spacing (e.g. lon0_deg = -130, dlon_deg = 0.25 over CONUS); the zonal circle then does NOT close, so use the zero-gradient (or a future one-sided) lon rule, never the periodic one. PHYSICAL-UNITS RECIPE: every rule here differentiates PER DEGREE of arc; on a sphere of radius R_sphere (a further free name, meters, consumer-supplied) the eastward and northward physical derivatives are df/dx_east = (180/pi) D(f, lon) / (R_sphere * coslat) and df/dy_north = (180/pi) D(f, lat) / R_sphere — compose these metric factors in the consuming equation (the `coslat` template below supplies cos(lat) over the lat axis; 180/pi = 57.29577951308232). Keeping the metric composition in the consumer keeps the rule count low and the rules exact on the index space. 

FACE (NODE) AXES, FOR THE CONSERVATIVE FLUX-FORM OPERATORS. Alongside the two CELL/POINT axes this grid declares two N+1 EDGE axes — `lon_nodes` (size NLON+1) and `lat_nodes` (size NLAT+1), esm-spec 9.7.6 metaparameter-expression sizes folded at load. They carry the face-staggered horizontal winds and the face-length metrics of the flux-form transport operators, exactly as `lev_nodes` does in the vertical on grids/latlon3d. FACE INDEXING (pinned, and identical to the vertical convention): FACE j IS THE LOWER (southern) EDGE OF CELL j, so cell j is bounded by face j to the south and face j+1 to the north; faces 1 and NLAT+1 are the two domain-end walls. lon_nodes follows the same rule (face i is the western edge of cell i) and is declared here so the contract is frozen before the zonal flux-form work lands; nothing in this grid reads it yet. Declaring `lat_nodes` as a distinct index-set NAME is also what makes the flux-form rules DISCRIMINABLE: an esm-spec 9.6.1 `where` shape constraint matches on the index-set NAME list, so a face wind of shape [lon, lat_nodes, lev] can never be mistaken for a cell tracer of shape [lon, lat, lev]. 

HALF-POLAR FINITE VOLUMES FROM POINT LATITUDES (the `lat_edge` derivation). Latitude is POINT-based and carries no edges of its own, so a finite-volume operator has to manufacture them. The rule pinned here: put each interior edge HALFWAY BETWEEN the two neighbouring points and CLAMP at the domain ends, phie[j] = max(lat0_deg, min(latN_deg, lat0_deg + (j - 1.5) dlat_deg)) for j = 1..NLAT+1, with latN_deg = lat0_deg + (NLAT - 1) dlat_deg. The min/max make it branch-free. On the GLOBAL recipe (lat0_deg = -90, dlat_deg = 180/(NLAT-1)) the j = 1 edge clamps to -90 and the j = NLAT+1 edge clamps to +90, so the two end cells span only HALF a nominal grid interval: this reproduces GEOS-Chem's HALF-POLAR CAPS as a consequence of the clamp rather than as a special case, and it does so on the regional recipe too (there the clamp simply halves the two end rows). Interior cells span a full dlat_deg. 

DERIVED LATITUDE GEOMETRY (the keyed-factor contract, mirroring grids/latlon3d's pe -> lev_thickness/lev_center). The CONSUMER declares one observed `phie` of shape [lat_nodes] by applying `lat_edge` (which takes no arguments — it folds from the free names lat0_deg/dlat_deg and the NLAT metaparameter), then feeds that SAME `phie` to the three derived templates: `coslat_edge` (shape [lat_nodes]) = cos(rad(phie[j])), the FACE-LENGTH weight; `latcell_area` (shape [lat]) = sin(rad(phie[j+1])) - sin(rad(phie[j])), the EXACT spherical cell-area weight (an exact integral of cos(phi) dphi, not a midpoint approximation — this is what makes the flux-form divergence conservative to machine precision); and `latcell_dphi` (shape [lat]) = rad(phie[j+1]) - rad(phie[j]), the cell angular width, a HALF-width at the two caps. rad(x) = x pi/180. Angles go in as DEGREES and come out as RADIANS (widths) or dimensionless (areas, cosines); on a sphere of radius R the physical cell area is R^2 dlon_rad latcell_area[j] and the physical face length is R dlon_rad coslat_edge[j]. 

AGGREGATE-INDEX HYGIENE, AND WHY THE FOUR TEMPLATES USE FOUR DIFFERENT INDEX NAMES (gje / gjc / gja / gjd). Template application is inlined NON-HYGIENICALLY: an aggregate index in an applied template body is NOT renamed apart from an identically-named index in the surrounding (or the supplied) expression. Two consequences, both load-bearing. (1) None of these names may be `j` (nor `i`/`k`), or a derived observed inlined into a consuming stencil would capture that stencil's own aggregate index — the footgun documented in grids/cartesian_nonuniform_1d and in grids/latlon3d's `gk`. (2) More subtly, the PRODUCER (`lat_edge`, index `gje`) and each CONSUMER (`coslat_edge` `gjc`, `latcell_area` `gja`, `latcell_dphi` `gjd`) must not share an index name EITHER: when the consumer's body is inlined it gathers the producer's inlined body at an OFFSET (phie[gja + 1]), and if the producer aggregated over that same name the offset is silently swallowed — the gather degenerates to phie[gja] and latcell_area returns an array of EXACT ZEROS (which then divides). That failure is silent (no error, no NaN — just zeros, and Inf downstream), so the four distinct names are a correctness requirement, not a style choice. A consuming model must likewise not name its own aggregate index gje/gjc/gja/gjd. 

THE POLE FACES ARE NEVER READ. On the global recipe coslat_edge[1] and coslat_edge[NLAT+1] are cos(-/+pi/2), which in IEEE double is 6.1e-17, NOT zero — so a flux-form operator that MULTIPLIED the pole-wall flux by the pole face length would leak a tiny but nonzero amount of mass across the pole. The flux-form lat rules therefore OMIT the two pole-wall flux terms from the expression entirely (the same literal-absence discipline the vertical no-flux stencils use at the ground and the model lid), which makes the polar no-flux wall exact in floating point. A consequence: whatever a consumer stores in the pole faces of a [lon, lat_nodes, ...] meridional wind is IGNORED, so a nonzero value there is harmless rather than a silent mass leak.

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NLON` | integer | `144` | Number of cells along longitude (cell-based axis; the global recipe closes the zonal circle with dlon_deg = 360/NLON). |
| `NLAT` | integer | `91` | Number of points along latitude (point-based, endpoint-inclusive axis; the global recipe is pole-to-pole with dlat_deg = 180/(NLAT - 1)). |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>R_sphere</code> <code>dlat_deg</code> <code>dlon_deg</code> <code>lat0_deg</code> <code>lon0_deg</code></p>

Region-generic: the same files serve a **global** or a **regional** grid depending on the free-name values the consumer supplies. GLOBAL recipe — `lon0_deg = 0`, `dlon_deg = 360/NLON`, `lat0_deg = −90`, `dlat_deg = 180/(NLAT − 1)` (pole-to-pole; the zonal circle closes, so the periodic-lon rule applies). REGIONAL recipe — any other origin/spacing (e.g. `lon0_deg = −130`, `dlon_deg = 0.25` over CONUS); the circle does not close, so use the zero-gradient lon rule. Spell the spacings as observeds dividing by the metaparameter names (`NLON`/`NLAT`) so a convergence sweep stays consistent. A further free name `R_sphere` (sphere radius, m) enters only the consumer's physical-units metric composition (per-degree derivatives → per-metre; the `coslat` template supplies cos(lat)).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `lon` | interval | `NLON` |
| `lat` | interval | `NLAT` |
| `lon_nodes` | interval | `{'op': '+', 'args': ['NLON', 1]}` |
| `lat_nodes` | interval | `{'op': '+', 'args': ['NLAT', 1]}` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `lon_coord` | — | `sum[i] (lon0_deg + (i - 0.5) * dlon_deg) where {i in lon}` | Cell-center longitudes lon_i = lon0_deg + (i - 1/2) dlon_deg (degrees), shaped over the lon index set. lon0_deg and dlon_deg are free names supplied by the consuming model (see metadata.description for the contract and the global recipe). |
| `lat_coord` | — | `sum[j] (lat0_deg + (j - 1) * dlat_deg) where {j in lat}` | Point latitudes lat_j = lat0_deg + (j - 1) dlat_deg (degrees), shaped over the lat index set (endpoint-inclusive; the global recipe runs pole to pole). lat0_deg and dlat_deg are free names supplied by the consuming model (see metadata.description). |
| `coslat` | — | `cos(3.141592653589793 / 180 * lat_coord<>)` | Metric factor cos(lat_j) over the lat index set: the cosine of lat_coord converted to radians inside (pi/180 with pi = 3.141592653589793). Used by consumers composing physical-units zonal derivatives, df/dx_east = (180/pi) D(f, lon) / (R_sphere * coslat) — see metadata.description. Inherits lat_coord's free names lat0_deg and dlat_deg. |
| `deg2rad` | `angle_deg` | `angle_deg * (3.141592653589793 / 180)` | Degree-to-radian helper: angle_deg * pi/180 (pi = 3.141592653589793). Invoke with any scalar or array expression binding; pure sugar for metric-factor compositions. |
| `lat_edge` | — | `sum[gje] (max(lat0_deg, min((lat0_deg + (NLAT - 1) * dlat_deg), (lat0_deg + (gje - 1.5) * dlat_deg)))) where {gje in lat_nodes}` | Latitude CELL EDGES in DEGREES, over the lat_nodes (NLAT+1 face) axis: phie[gje] = max(lat0_deg, min(latN_deg, lat0_deg + (gje - 1.5) dlat_deg)), latN_deg = lat0_deg + (NLAT - 1) dlat_deg. Each interior edge is the midpoint of the two neighbouring lat POINTS; the two end edges CLAMP to the domain ends, which halves the two end cells and so reproduces GEOS-Chem's half-polar caps on the global recipe (see metadata.description). Branch-free (min/max, no ifelse). Face gje is the LOWER (southern) edge of cell gje. Assumes dlat_deg > 0. Takes no arguments: it folds from the free names lat0_deg and dlat_deg plus the NLAT metaparameter. THE CONSUMER declares an observed `phie` of shape [lat_nodes] by applying this template, then feeds that same `phie` to coslat_edge / latcell_area / latcell_dphi — the keyed-factor contract grids/latlon3d uses for pe -> lev_thickness. This is the PRODUCER template, and its aggregate index `gje` is deliberately DISTINCT from the consumers' (gjc/gja/gjd): non-hygienic inlining would otherwise let a consumer's offset gather phie[gj+1] be captured by this body's own index and silently lose the +1 (see the AGGREGATE-INDEX HYGIENE paragraph in metadata.description). |
| `coslat_edge` | `phie` | `sum[gjc] (cos(3.141592653589793 / 180 * phie[gjc])) where {gjc in lat_nodes}` | FACE-LENGTH metric over the lat_nodes (NLAT+1 face) axis: coslat_edge[gjc] = cos(pi/180 * phie[gjc]), where phie is the consumer's degree-valued edge array (shape [lat_nodes], from lat_edge). This is the zonal length of the latitude face, in units of R dlon_rad. On the global recipe the two POLE faces have coslat_edge = cos(-/+pi/2) = 6.1e-17 in IEEE double, NOT 0 — which is precisely why the flux-form lat rules OMIT the pole-wall flux terms rather than multiplying them by this near-zero factor (see metadata.description). Aggregate index `gjc`, distinct from the producer's `gje` and from its sibling consumers' `gja`/`gjd` (hygiene; see metadata.description). |
| `latcell_area` | `phie` | `sum[gja] (sin(3.141592653589793 / 180 * phie[gja + 1]) - sin(3.141592653589793 / 180 * phie[gja])) where {gja in lat}` | EXACT spherical cell-area weight over the lat (cell) axis: latcell_area[gja] = sin(pi/180 * phie[gja+1]) - sin(pi/180 * phie[gja]), gathering the consumer's NLAT+1-node degree-valued edge array phie. This is the exact integral of cos(phi) dphi across the cell, not a midpoint approximation, so the flux-form meridional divergence (F[j+1] coslat_edge[j+1] - F[j] coslat_edge[j]) / latcell_area[j] telescopes to an exact global mass budget. Positive when dlat_deg > 0. The physical cell area on a sphere of radius R is R^2 dlon_rad latcell_area[gja]. At the two half-polar caps latcell_area ~ dphi^2/8 — SECOND order small in the grid spacing while the adjacent face is only FIRST order small (coslat_edge ~ dphi/2); that O(1/dphi) amplification is why near-pole convergence requires the product V dq/dphi to vanish like cos(phi) at the pole (see grids/latlon3d/rules/upwind1_flux_D_lat_polar). Aggregate index `gja`, distinct from the producer's `gje` — with a shared name the offset gather phie[gja+1] is captured and silently degenerates to phie[gja], returning EXACT ZEROS (see the AGGREGATE-INDEX HYGIENE paragraph in metadata.description). |
| `latcell_dphi` | `phie` | `sum[gjd] (3.141592653589793 / 180 * (phie[gjd + 1] - phie[gjd])) where {gjd in lat}` | Cell angular WIDTH in RADIANS over the lat (cell) axis: latcell_dphi[gjd] = pi/180 * (phie[gjd+1] - phie[gjd]), gathering the consumer's NLAT+1-node degree-valued edge array phie. Interior cells are a full dlat_deg wide; the two end (polar-cap) cells are HALF that, a direct consequence of the lat_edge clamp. The physical north-south cell width on a sphere of radius R is R latcell_dphi[gjd]. Needed by the ZONAL flux-form operator (which weights the zonal divergence by the cell width) and by any consumer wanting the meridional spacing; the MERIDIONAL flux-form operator does NOT need it, since latcell_area already carries the exact area weight. Aggregate index `gjd`, distinct from the producer's `gje` for the offset-gather hygiene reason in metadata.description. |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `central_D2_lat_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/central_D2_lat_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/central_D2_lat_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i, j + 1] + -2 * f[i, j] + f[i, j - 1]) / (dlat_deg * dlat_deg)) where {i in lon, j in 2:NLAT - 1}</div>

Interior-only 3-point centered second derivative along latitude on the latlon grid: (f[i,j+1] - 2 f[i,j] + f[i,j-1]) / dlat_deg^2 over all i, j in [2, NLAT-1] — a PER-DEGREE second derivative (the diffusion/Laplacian companion to central_D_lat_interior; see grids/latlon/grid.esm for the physical-units metric recipe, of which the full metric spherical Laplacian with cos(lat) weighting is a documented FOLLOW-UP, not this operator). Match-less named template consumed by boundary-condition layers (rules/). `dlat_deg` is a free name resolving to the consuming model's meridional spacing variable (the grid's consumer-supplied geometry contract). The 3-point centered weights [1, -2, 1] are the Fornberg formula for the 2nd derivative on a uniform stencil; O(dlat_deg^2) on this uniform axis.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2_lon_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/central_D2_lon_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/central_D2_lon_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i + 1, j] + -2 * f[i, j] + f[i - 1, j]) / (dlon_deg * dlon_deg)) where {i in 2:NLON - 1, j in lat}</div>

Interior-only 3-point centered second derivative along longitude on the latlon grid: (f[i+1,j] - 2 f[i,j] + f[i-1,j]) / dlon_deg^2 over i in [2, NLON-1], all j — a PER-DEGREE second derivative (the diffusion/Laplacian companion to central_D_lon_interior; see grids/latlon/grid.esm for the physical-units metric recipe, of which the full metric spherical Laplacian with cos(lat) weighting is a documented FOLLOW-UP, not this operator). Match-less named template consumed by boundary-condition layers (rules/). `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (the grid's consumer-supplied geometry contract). The 3-point centered weights [1, -2, 1] are the Fornberg formula for the 2nd derivative on a uniform stencil; O(dlon_deg^2) on this uniform axis.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D_lat_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/central_D_lat_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/central_D_lat_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i, j + 1] - f[i, j - 1]) / (2 * dlat_deg)) where {i in lon, j in 2:NLAT - 1}</div>

Interior-only centered first derivative along latitude on the latlon grid: (f[i,j+1] - f[i,j-1]) / (2 dlat_deg) over all i, j in [2, NLAT-1] — a PER-DEGREE derivative (see grids/latlon/grid.esm for the physical-units metric recipe). Match-less named template consumed by boundary-condition layers (rules/). `dlat_deg` is a free name resolving to the consuming model's meridional spacing variable (the grid's consumer-supplied geometry contract). O(dlat_deg^2) on this uniform axis.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D_lon_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/central_D_lon_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/central_D_lon_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i + 1, j] - f[i - 1, j]) / (2 * dlon_deg)) where {i in 2:NLON - 1, j in lat}</div>

Interior-only centered first derivative along longitude on the latlon grid: (f[i+1,j] - f[i-1,j]) / (2 dlon_deg) over i in [2, NLON-1], all j — a PER-DEGREE derivative (see grids/latlon/grid.esm for the physical-units metric recipe). Match-less named template consumed by boundary-condition layers (rules/). `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (the grid's consumer-supplied geometry contract). O(dlon_deg^2) on this uniform axis.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `laplace_beltrami_interior`

<p class="rule-meta"><code>op:laplacian</code> <code>order:2</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/laplace_beltrami_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/laplace_beltrami_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] ((u[i + 1, j] + -2 * u[i, j] + u[i - 1, j]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + (cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + 0.5) * dlat_deg)) * (u[i, j + 1] - u[i, j]) - cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + -0.5) * dlat_deg)) * (u[i, j] - u[i, j - 1])) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 2:NLON - 1, j in 2:NLAT - 1}</div>

Interior-only metric spherical Laplacian (Laplace-Beltrami) on the latlon grid: the doubly-interior 5-point cross value lap u = (1/(R_sphere^2 cos^2 phi)) d2u/dlambda^2 + (1/(R_sphere^2 cos phi)) d/dphi(cos phi du/dphi) over i in [2,NLON-1], j in [2,NLAT-1] (lambda=lon, phi=lat in RADIANS; dlambda=(pi/180)dlon_deg, dphi=(pi/180)dlat_deg, pi=3.141592653589793). ZONAL term: the centered per-radian second difference (u[i+1,j]-2u[i,j]+u[i-1,j])/(R_sphere^2 cos^2 phi_j dlambda^2) -- the per-degree D2 of central_D2_lon_interior scaled by the (180/pi)^2/(R_sphere^2 cos^2 phi_j) metric weight (the cos(lat) weighting that distinguishes this from the bare per-degree operator). MERIDIONAL term: the conservative (mimetic) flux-form divergence of the cos-weighted gradient, (cos(phi_{j+1/2})(u[i,j+1]-u[i,j]) - cos(phi_{j-1/2})(u[i,j]-u[i,j-1])) /(R_sphere^2 cos phi_j dphi^2), with the cell-face latitudes phi_{j+-1/2} at the half-index nodes -- the standard 2nd-order finite-volume Laplace-Beltrami stencil that annihilates constants and is discretely conservative in the meridional flux. cos(phi) at cell centers and faces is computed IN THIS BODY from phi_j=(pi/180)(lat0_deg+(j-1)dlat_deg) and phi_{j+-1/2}=(pi/180)(lat0_deg+(j-1+-1/2)dlat_deg): free-name geometry (lat0_deg, dlat_deg, dlon_deg, R_sphere) resolves in the consuming model's scope at evaluation (the grid's consumer-supplied geometry contract, grids/latlon/grid.esm). VALID ONLY on a REGIONAL spherical patch with lat strictly inside (-90,90) so cos phi>0 everywhere (a spherical cap/band); the global pole-cap treatment of the 1/cos^2 zonal singularity at phi=+-90 is a documented FOLLOW-UP. Match-less named template consumed by the boundary-condition rule layer (rules/laplace_beltrami_lon_periodic_lat_zero_grad). O(dlon_deg^2)+O(dlat_deg^2) on this uniform index space.

**References**

- Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992). A standard test set for numerical approximations to the shallow water equations in spherical geometry. Journal of Computational Physics, 102(1), 211-224. (Laplace-Beltrami / covariant Laplacian on the sphere.)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Chapter 2 (finite-volume flux-form divergence of a variable-coefficient gradient) and the [1,-2,1] second difference.

### `upwind1_D_lon_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/latlon/stencils/upwind1_D_lon_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/stencils/upwind1_D_lon_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i, j] - f[i - 1, j]) / dlon_deg) where {i in 2:NLON, j in lat}</div>

Interior-only first-order upwind (backward) first derivative along longitude on the latlon grid: (f[i,j] - f[i-1,j]) / dlon_deg over i in [2, NLON], all j — a PER-DEGREE derivative (see grids/latlon/grid.esm for the physical-units metric recipe), for zonal advection with non-negative angular velocity (flow toward increasing longitude). Match-less named template consumed by boundary-condition layers (rules/) — e.g. rules/upwind1_D_lon_periodic.esm adds the i=1 periodic wrap. This is the one-sided upwind companion to the centered stencils/central_D_lon_interior.esm; it keeps the full interior span [2, NLON] because only the i=1 column lacks an upstream neighbor (contrast the centered stencil, which drops both end columns). `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (the grid's consumer-supplied geometry contract). O(dlon_deg) on this uniform axis. Sign convention: upwind for flow toward +lon; a westward variant belongs in a separate stencil.

**References**

- LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge University Press. Section 4.8 (upwind methods).

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `central_D2_lat_zero_grad_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:lon,lat</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/latlon/rules/central_D2_lat_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/central_D2_lat_zero_grad_bc.esm)

**Rewrites** `D(D(f)/Dlat)/Dlat` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "lat"}], "wrt": "lat"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLON, 2:NLAT - 1] = central_D2_lat_interior&lt;f=f&gt;, [1:NLON, 1:1] = sum[i, j] ((f[i, 2] - f[i, 1]) / (dlat_deg * dlat_deg)) where {i in lon, j in 1:1}, [1:NLON, NLAT:NLAT] = sum[i, j] ((f[i, NLAT - 1] - f[i, NLAT]) / (dlat_deg * dlat_deg)) where {i in lon, j in NLAT:NLAT})</div>

Centered second derivative (per-degree diffusion / 1-D Laplacian) along latitude with zero-gradient (homogeneous Neumann) pole rows on the latlon grid. Compound match on D(D(f, wrt: lat), wrt: lat) at priority 10, so it fires on the whole second derivative before any plain-D rule could lower the inner meridional derivative (esm-spec 9.6.3 compound precedence); a PER-DEGREE operator (the physical northward diffusion composes the (180/pi)^2 / R_sphere^2 metric factor in the consuming equation — see grids/latlon/grid.esm; the full metric spherical Laplacian with cos(lat) weighting is a documented FOLLOW-UP, not this rule). Interior rows [2, NLAT-1]: the imported 3-point centered stencil (f[i,j+1] - 2 f[i,j] + f[i,j-1]) / dlat_deg^2. Pole rows j=1 and j=NLAT: ghost elimination with f[i,0]=f[i,1] and f[i,NLAT+1]=f[i,NLAT] (the same zero-gradient ghost closure as cartesian_uniform_1d/rules/central_D2_zero_grad_bc, lifted to the lat axis), giving (f[i,2]-f[i,1])/dlat_deg^2 and (f[i,NLAT-1]-f[i,NLAT])/dlat_deg^2; the three regions tile the axes, so the discretized second derivative is fully defined with its BC (esm-spec 9.6.8). The face values are FULL-RANK aggregates (singleton output range on the lat axis) so every binding's simulator fills the (NLON x 1) regions without rank promotion. POLE-TREATMENT ACCURACY (important, and different from the cell-centered cartesian case): on the GLOBAL recipe the lat axis is POINT-BASED with lat_1 = -90 and lat_NLAT = +90 exactly AT the poles (not half a cell inside, as cartesian_uniform_1d's cell centers are). At a node-coincident boundary the ghost-elimination face (f[i,2]-f[i,1])/dlat_deg^2 is the finite-volume Neumann Laplacian evaluated half a cell off the node, so its leading truncation is O(1) unless the field's first THREE lat-derivatives vanish at the pole (a plain zero-gradient / cos(k*lat) eigenmode, whose second derivative is nonzero at the pole, degrades the observed order to 1; see the de-risk note in problems/latlon_diffusion_meridional_zero_grad.esm). For a pole-supercompatible field (all of d/dlat, d2/dlat2, d3/dlat3 zero at the poles — e.g. cos^4((pi/180) lat)) the face residual is O(dlat_deg^2) and the global scheme is cleanly 2nd order; such fields are exactly those a regular scalar on the sphere produces near the polar cap. This is a deliberately CRUDE pole treatment (as for the first-derivative central_D_lat_zero_grad_bc): it neither couples the pole rows across longitude nor conserves anything through the pole; proper pole handling (averaged polar caps, or an unstructured cap) is a different rule/grid. On a REGIONAL grid the same rule is an ordinary zero-gradient (homogeneous Neumann) frame boundary. `dlat_deg` is a free name resolving to the consuming model's meridional spacing variable (consumer-supplied geometry contract). Match scoping: the pattern is the compound `wrt: lat` twice plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1), so it fires only on the second meridional derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(D(., wrt: lat), wrt: lat) reusing the axis name, and, under import-edge renaming (esm-spec 9.7.7), the wrt literals and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (Neumann boundary conditions by ghost-point elimination) and Chapter 2 (the [1, -2, 1] second-difference operator).
- Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992). A standard test set for numerical approximations to the shallow water equations in spherical geometry. Journal of Computational Physics, 102(1), 211-224. (Context for pole treatments on latitude-longitude grids.)

**Conformance fixtures**

- ast: [`central_D2_lat_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2_lat_zero_grad_bc) — status `active`
- convergence: [`latlon_diffusion_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_diffusion_meridional_zero_grad) — status `active`
- simulation: [`latlon_diffusion_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/latlon_diffusion_meridional_zero_grad) — status `active`

#### Convergence — case [`latlon_diffusion_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_diffusion_meridional_zero_grad)

Error norms read from the committed golden [`tests/conformance/convergence/latlon_diffusion_meridional_zero_grad/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/latlon_diffusion_meridional_zero_grad/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.546e-03 | — | 1.424e-03 | — |
| 32 | 3.932e-04 | 1.98 | 3.597e-04 | 1.99 |
| 64 | 9.897e-05 | 1.99 | 9.016e-05 | 2.00 |
| 128 | 2.482e-05 | 2.00 | 2.255e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of latlon_diffusion_meridional_zero_grad (log-log)](/plots/convergence/latlon_diffusion_meridional_zero_grad.png)

### `central_D2_lon_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:lon,lat</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/latlon/rules/central_D2_lon_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/central_D2_lon_periodic.esm)

**Rewrites** `D(D(f)/Dlon)/Dlon` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "lon"}], "wrt": "lon"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLON - 1, 1:NLAT] = central_D2_lon_interior&lt;f=f&gt;, [1:1, 1:NLAT] = sum[i, j] ((f[2, j] + -2 * f[1, j] + f[NLON, j]) / (dlon_deg * dlon_deg)) where {i in 1:1, j in lat}, [NLON:NLON, 1:NLAT] = sum[i, j] ((f[1, j] + -2 * f[NLON, j] + f[NLON - 1, j]) / (dlon_deg * dlon_deg)) where {i in NLON:NLON, j in lat})</div>

Centered second derivative (per-degree diffusion / 1-D Laplacian) along longitude with periodic (wrap) boundary on the latlon grid — the physically natural rule for a closed zonal circle (global recipe: dlon_deg = 360/NLON; see grids/latlon/grid.esm). Compound match on D(D(f, wrt: lon), wrt: lon) at priority 10, so it fires on the whole second derivative before any plain-D rule could lower the inner longitudinal derivative (esm-spec 9.6.3 compound precedence); a PER-DEGREE operator (the physical eastward diffusion composes the (180/pi)^2 / (R_sphere coslat)^2 metric factors in the consuming equation — see grids/latlon/grid.esm; the full metric spherical Laplacian with cos(lat) weighting is a documented FOLLOW-UP, not this rule). Interior columns [2, NLON-1]: the imported 3-point centered stencil (f[i+1,j] - 2 f[i,j] + f[i-1,j]) / dlon_deg^2. The two face columns wrap explicitly, following the central_D_lon_periodic / cartesian upwind1_D_periodic precedent: i=1 gathers (f[2,j] - 2 f[1,j] + f[NLON,j]) / dlon_deg^2 and i=NLON gathers (f[1,j] - 2 f[NLON,j] + f[NLON-1,j]) / dlon_deg^2, so the scheme is the same O(dlon_deg^2) centered second difference at every column and conserves sum(f) exactly (each column of the operator telescopes to zero net weight). The face values are FULL-RANK aggregates (singleton output range on the lon axis) rather than the esm-spec 9.6.8 worked example's rank-reduced form, so every binding's simulator fills the (1 x NLAT) regions without rank promotion. `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (consumer-supplied geometry contract); the rule itself never assumes the global extent, but composing it with a non-closing regional grid is a modeling error — use a zero-gradient lon second-derivative rule there. Match scoping: the pattern is the compound `wrt: lon` twice plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1), so it fires only on the second longitudinal derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(D(., wrt: lon), wrt: lon) reusing the axis name, and, under import-edge renaming (esm-spec 9.7.7), the wrt literals and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains) and Chapter 2 (the [1, -2, 1] second-difference operator).

**Conformance fixtures**

- ast: [`central_D2_lon_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2_lon_periodic) — status `active`
- convergence: [`latlon_diffusion_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_diffusion_zonal_periodic) — status `active`
- simulation: [`latlon_diffusion_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/latlon_diffusion_zonal_periodic) — status `active`

#### Convergence — case [`latlon_diffusion_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_diffusion_zonal_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/latlon_diffusion_zonal_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/latlon_diffusion_zonal_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 7.792e-04 | — | 7.191e-04 | — |
| 32 | 1.955e-04 | 1.99 | 1.831e-04 | 1.97 |
| 64 | 4.892e-05 | 2.00 | 4.597e-05 | 1.99 |
| 128 | 1.223e-05 | 2.00 | 1.151e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of latlon_diffusion_zonal_periodic (log-log)](/plots/convergence/latlon_diffusion_zonal_periodic.png)

### `central_D_lat_zero_grad_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:lon,lat</code> <code>stencil width 3</code></p>

Source: [`grids/latlon/rules/central_D_lat_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/central_D_lat_zero_grad_bc.esm)

**Rewrites** `D(f)/Dlat` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "lat"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLON, 2:NLAT - 1] = central_D_lat_interior&lt;f=f&gt;, [1:NLON, 1:1] = sum[i, j] ((f[i, 2] - f[i, 1]) / dlat_deg) where {i in lon, j in 1:1}, [1:NLON, NLAT:NLAT] = sum[i, j] ((f[i, NLAT] - f[i, NLAT - 1]) / dlat_deg) where {i in lon, j in NLAT:NLAT})</div>

Centered first derivative along latitude with zero-gradient boundary rows on the latlon grid. Matches D(f, wrt: lat), a PER-DEGREE derivative (see grids/latlon/grid.esm for the physical-units metric recipe). Interior rows [2, NLAT-1]: the imported centered-difference stencil. Boundary rows j=1 and j=NLAT: one-sided differences (f[i,2] - f[i,1]) / dlat_deg and (f[i,NLAT] - f[i,NLAT-1]) / dlat_deg encoding the zero-gradient condition; the three regions tile the axes, so the discretized derivative is fully defined with its BC (esm-spec 9.6.8). On the GLOBAL recipe (lat0_deg = -90, dlat_deg = 180/(NLAT-1)) the boundary rows are the poles themselves, and zero-gradient there is a deliberately CRUDE pole treatment — standard for exemplar transport, adequate for fields whose meridional gradient vanishes at the poles (any regular scalar field on the sphere has d f/d lat = 0 at the poles along every meridian), but it neither couples the pole rows across longitude nor conserves anything through the pole; proper pole handling (averaged polar caps, or an unstructured cap) is a different rule/grid. On a REGIONAL grid the same rule is an ordinary zero-gradient (homogeneous Neumann) frame boundary. The face values are FULL-RANK aggregates (singleton output range on the lat axis) so every binding's simulator fills the (NLON x 1) regions without rank promotion — the rank-reduced face form of the esm-spec 9.6.8 worked example only fills leading-axis faces portably. `dlat_deg` is a free name resolving to the consuming model's meridional spacing variable (consumer-supplied geometry contract). The one-sided faces are O(dlat_deg) locally, which preserves global 2nd-order convergence for boundary-compatible solutions (see the MMS problem problems/latlon_advection_meridional_zero_grad.esm). Match scoping: the pattern is `wrt: lat` plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1), so it fires only on the meridional derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(., wrt: lat) reusing the axis name, and, under import-edge renaming (esm-spec 9.7.7), the wrt literal and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided boundary treatment).
- Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992). A standard test set for numerical approximations to the shallow water equations in spherical geometry. Journal of Computational Physics, 102(1), 211-224. (Context for pole treatments on latitude-longitude grids.)

**Conformance fixtures**

- ast: [`central_D_lat_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D_lat_zero_grad_bc) — status `active`
- convergence: [`latlon_advection_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_meridional_zero_grad) — status `active`
- simulation: [`latlon_advection_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/latlon_advection_meridional_zero_grad) — status `active`

#### Convergence — case [`latlon_advection_meridional_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_meridional_zero_grad)

Error norms read from the committed golden [`tests/conformance/convergence/latlon_advection_meridional_zero_grad/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/latlon_advection_meridional_zero_grad/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 7.907e-03 | — | 7.164e-03 | — |
| 32 | 2.010e-03 | 1.98 | 1.819e-03 | 1.98 |
| 64 | 5.048e-04 | 1.99 | 4.568e-04 | 1.99 |
| 128 | 1.263e-04 | 2.00 | 1.143e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of latlon_advection_meridional_zero_grad (log-log)](/plots/convergence/latlon_advection_meridional_zero_grad.png)

### `central_D_lon_periodic`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:lon,lat</code> <code>stencil width 3</code></p>

Source: [`grids/latlon/rules/central_D_lon_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/central_D_lon_periodic.esm)

**Rewrites** `D(f)/Dlon` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "lon"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLON - 1, 1:NLAT] = central_D_lon_interior&lt;f=f&gt;, [1:1, 1:NLAT] = sum[i, j] ((f[2, j] - f[NLON, j]) / (2 * dlon_deg)) where {i in 1:1, j in lat}, [NLON:NLON, 1:NLAT] = sum[i, j] ((f[1, j] - f[NLON - 1, j]) / (2 * dlon_deg)) where {i in NLON:NLON, j in lat})</div>

Centered first derivative along longitude with periodic (wrap) boundary on the latlon grid — the physically natural rule for a closed zonal circle (global recipe: dlon_deg = 360/NLON; see grids/latlon/grid.esm). Matches D(f, wrt: lon), a PER-DEGREE derivative. Interior columns [2, NLON-1]: the imported centered-difference stencil. The two face columns wrap explicitly, following the cartesian upwind1_D_periodic / ESS pde_simulation fixture precedent for periodic boundaries: i=1 gathers (f[2,j] - f[NLON,j]) / (2 dlon_deg) and i=NLON gathers (f[1,j] - f[NLON-1,j]) / (2 dlon_deg), so the scheme is the same O(dlon_deg^2) centered difference at every column and conserves sum(f) exactly (the wrap telescopes). The face values are FULL-RANK aggregates (singleton output range on the lon axis) rather than the esm-spec 9.6.8 worked example's rank-reduced form, so every binding's simulator fills the (1 x NLAT) regions without rank promotion. `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (consumer-supplied geometry contract); the rule itself never assumes the global extent, but composing it with a non-closing regional grid is a modeling error — use central_D_lon_zero_grad_bc there. Match scoping: the pattern is `wrt: lon` plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1), so it fires only on the longitudinal derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(., wrt: lon) reusing the axis name, and, under import-edge renaming (esm-spec 9.7.7), the wrt literal and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D_lon_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D_lon_periodic) — status `active`
- convergence: [`latlon_advection_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_zonal_periodic) — status `active`
- simulation: [`latlon_advection_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/latlon_advection_zonal_periodic) — status `active`

#### Convergence — case [`latlon_advection_zonal_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_zonal_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/latlon_advection_zonal_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/latlon_advection_zonal_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.003e-02 | — | 1.968e-02 | — |
| 32 | 5.037e-03 | 1.99 | 5.014e-03 | 1.97 |
| 64 | 1.261e-03 | 2.00 | 1.260e-03 | 1.99 |
| 128 | 3.154e-04 | 2.00 | 3.153e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of latlon_advection_zonal_periodic (log-log)](/plots/convergence/latlon_advection_zonal_periodic.png)

### `central_D_lon_zero_grad_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:lon,lat</code> <code>stencil width 3</code></p>

Source: [`grids/latlon/rules/central_D_lon_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/central_D_lon_zero_grad_bc.esm)

**Rewrites** `D(f)/Dlon` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "lon"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLON - 1, 1:NLAT] = central_D_lon_interior&lt;f=f&gt;, [1:1, 1:NLAT] = sum[j] ((f[2, j] - f[1, j]) / dlon_deg) where {j in lat}, [NLON:NLON, 1:NLAT] = sum[j] ((f[NLON, j] - f[NLON - 1, j]) / dlon_deg) where {j in lat})</div>

Centered first derivative along longitude with zero-gradient boundary conditions on the latlon grid — the esm-spec 9.6.8 worked example, generalized over NLON/NLAT metaparameters and the consumer-supplied spacing free name. Matches D(f, wrt: lon), a PER-DEGREE derivative (see grids/latlon/grid.esm for the physical-units metric recipe). Interior columns [2, NLON-1]: the imported centered-difference stencil. Face columns i=1 and i=NLON: one-sided differences (f[2,j]-f[1,j])/dlon_deg and (f[NLON,j]-f[NLON-1,j])/dlon_deg encoding the zero-gradient condition. The three regions tile the axes, so the discretized derivative is fully defined with its BC (there is nowhere else a boundary condition could live). This is the REGIONAL boundary choice — a zero-gradient frame for a limited-area lon extent; on the global recipe the zonal circle closes and central_D_lon_periodic is the physical rule (a zero gradient on a closed circle is an artificial seam). `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (consumer-supplied geometry contract) — the zero-gradient faces themselves are spacing-scale-invariant, so the rule serves any extent unchanged. The face values keep the rank-reduced aggregate form of the spec's worked example (output over the lat axis only, filling the (1 x NLAT) leading-axis face regions); the lat-axis rules use full-rank faces instead because rank-reduced values only fill leading-axis faces portably across bindings. Match scoping: the pattern is `wrt: lon` plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1), so it fires only on the longitudinal derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(., wrt: lon) reusing the axis name, and, under import-edge renaming (esm-spec 9.7.7), the wrt literal and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided boundary treatment).

**Conformance fixtures**

- ast: [`central_D_lon_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D_lon_zero_grad_bc) — status `active`

### `laplace_beltrami_lon_periodic_lat_zero_grad`

<p class="rule-meta"><code>op:laplacian</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:lon,lat</code></p>

Source: [`grids/latlon/rules/laplace_beltrami_lon_periodic_lat_zero_grad.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/laplace_beltrami_lon_periodic_lat_zero_grad.esm)

**Rewrites** `laplacian(u)` — match pattern:

```json
{"op": "laplacian", "args": ["u"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLON - 1, 2:NLAT - 1] = laplace_beltrami_interior&lt;u=u&gt;, [1:1, 2:NLAT - 1] = sum[i, j] ((u[2, j] + -2 * u[1, j] + u[NLON, j]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + (cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + 0.5) * dlat_deg)) * (u[1, j + 1] - u[1, j]) - cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + -0.5) * dlat_deg)) * (u[1, j] - u[1, j - 1])) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 1:1, j in 2:NLAT - 1}, [NLON:NLON, 2:NLAT - 1] = sum[i, j] ((u[1, j] + -2 * u[NLON, j] + u[NLON - 1, j]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + (cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + 0.5) * dlat_deg)) * (u[NLON, j + 1] - u[NLON, j]) - cos(3.141592653589793 / 180 * (lat0_deg + (j - 1 + -0.5) * dlat_deg)) * (u[NLON, j] - u[NLON, j - 1])) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (j - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in NLON:NLON, j in 2:NLAT - 1}, [2:NLON - 1, 1:1] = sum[i, j] ((u[i + 1, 1] + -2 * u[i, 1] + u[i - 1, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1 + 0.5) * dlat_deg)) * (u[i, 2] - u[i, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 2:NLON - 1, j in 1:1}, [2:NLON - 1, NLAT:NLAT] = sum[i, j] ((u[i + 1, NLAT] + -2 * u[i, NLAT] + u[i - 1, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1 + -0.5) * dlat_deg)) * (u[i, NLAT - 1] - u[i, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 2:NLON - 1, j in NLAT:NLAT}, [1:1, 1:1] = sum[i, j] ((u[2, 1] + -2 * u[1, 1] + u[NLON, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1 + 0.5) * dlat_deg)) * (u[1, 2] - u[1, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 1:1, j in 1:1}, [NLON:NLON, 1:1] = sum[i, j] ((u[1, 1] + -2 * u[NLON, 1] + u[NLON - 1, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1 + 0.5) * dlat_deg)) * (u[NLON, 2] - u[NLON, 1]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (1 - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in NLON:NLON, j in 1:1}, [1:1, NLAT:NLAT] = sum[i, j] ((u[2, NLAT] + -2 * u[1, NLAT] + u[NLON, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1 + -0.5) * dlat_deg)) * (u[1, NLAT - 1] - u[1, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in 1:1, j in NLAT:NLAT}, [NLON:NLON, NLAT:NLAT] = sum[i, j] ((u[1, NLAT] + -2 * u[NLON, NLAT] + u[NLON - 1, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg))^2 * (3.141592653589793 / 180 * dlon_deg)^2) + cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1 + -0.5) * dlat_deg)) * (u[NLON, NLAT - 1] - u[NLON, NLAT]) / (R_sphere^2 * cos(3.141592653589793 / 180 * (lat0_deg + (NLAT - 1) * dlat_deg)) * (3.141592653589793 / 180 * dlat_deg)^2)) where {i in NLON:NLON, j in NLAT:NLAT})</div>

METRIC spherical Laplacian (Laplace-Beltrami) of a scalar field on a REGIONAL latlon patch: matches the bare open-namespace rewrite-target {op: laplacian, args: [u]} (esm-spec 9.6.8/4.2 author sugar with no evaluable core) constrained by where {u:{shape:[lon,lat]}}, exactly mirroring grids/mpas/rules/fv_laplacian_cell -- a rule MAY match laplacian directly and emit a one-shot scheme rather than route through the nested-D identity. Emits the metric-weighted makearray for lap u = (1/(R_sphere^2 cos^2 phi)) d2u/dlambda^2 + (1/(R_sphere^2 cos phi)) d/dphi(cos phi du/dphi) (lambda=lon, phi=lat, RADIANS: dlambda=(pi/180)dlon_deg, dphi=(pi/180)dlat_deg, pi=3.141592653589793). The doubly-interior region [2,NLON-1]x[2,NLAT-1] imports stencils/laplace_beltrami_interior; the eight boundary regions tile the frame with two independent, composable closures. LONGITUDE = PERIODIC (the patch spans the full zonal circle, NLON dlon_deg=360): the i=1 and i=NLON columns wrap the zonal second difference (i-1->NLON at i=1, i+1->1 at i=NLON), following the central_D2_lon_periodic precedent. LATITUDE = ZERO-GRADIENT (homogeneous Neumann band edges): the mimetic meridional flux through the outer face is dropped, giving one-sided rows j=1: cos(phi_{3/2})(u[i,2]-u[i,1])/(R_sphere^2 cos phi_1 dphi^2) and j=NLAT: cos(phi_{NLAT-1/2})(u[i,NLAT-1]-u[i,NLAT])/(R_sphere^2 cos phi_NLAT dphi^2), the cos-weighted analogue of central_D2_lat_zero_grad_bc's ghost-elimination. cos(phi) at cell centers and faces is computed in the RUNTIME rule body from phi_j=(pi/180)(lat0_deg+(j-1)dlat_deg) and phi_{j+-1/2}=(pi/180)(lat0_deg+(j-1+-1/2)dlat_deg); the free names R_sphere, dlon_deg, dlat_deg, lat0_deg resolve in the consuming model's scope (grids/latlon/grid.esm keyed-factor contract). POLE SINGULARITY (the crux, deliberately avoided here): the zonal weight 1/cos^2 phi diverges at phi=+-90, so this rule is DEFINED ONLY on a regional patch with lat in (lat0,lat1) strictly inside (-90,90) (a Laplace-Beltrami on a spherical band/cap where cos phi is bounded away from 0). The GLOBAL pole-cap Laplacian -- an averaged polar-cap closure or an unstructured cap that removes the 1/cos^2 singularity -- is a documented FOLLOW-UP, not this rule. ACCURACY: 2nd order in both L2 and Linf on a regional band with a pole-supercompatible (boundary-vanishing) field; annihilates constants (every second difference is zero) and the meridional flux form is discretely conservative through the band interior. Match scoping mirrors fv_laplacian_cell: the where shape constraint fires only on the laplacian of a bare field over THIS grid's lon x lat index sets and, under esm-spec 9.7.7 import renaming, the shape follows the renamed axes so two latlon patches coexist; a consumer taking laplacian of a compound inline expression must first bind it to a declared lon x lat observed.

**References**

- Williamson, D. L., Drake, J. B., Hack, J. J., Jakob, R., & Swarztrauber, P. N. (1992). A standard test set for numerical approximations to the shallow water equations in spherical geometry. Journal of Computational Physics, 102(1), 211-224.
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (Neumann boundary by dropping the boundary flux) and Chapter 2 (the flux-form finite-volume Laplacian).

**Conformance fixtures**

- ast: [`laplace_beltrami_lon_periodic_lat_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/laplace_beltrami_lon_periodic_lat_zero_grad) — status `active`
- convergence: [`laplace_beltrami_band_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/laplace_beltrami_band_mms) — status `active`
- simulation: [`laplace_beltrami_band_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/laplace_beltrami_band_mms) — status `active`

#### Convergence — case [`laplace_beltrami_band_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/laplace_beltrami_band_mms)

Error norms read from the committed golden [`tests/conformance/convergence/laplace_beltrami_band_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/laplace_beltrami_band_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.642e-02 | — | 2.381e-02 | — |
| 32 | 6.653e-03 | 1.99 | 6.341e-03 | 1.91 |
| 64 | 1.671e-03 | 1.99 | 1.615e-03 | 1.97 |
| 128 | 4.188e-04 | 2.00 | 4.061e-04 | 1.99 |

Expected order: **2** (± 0.2).

![Convergence of laplace_beltrami_band_mms (log-log)](/plots/convergence/laplace_beltrami_band_mms.png)

### `upwind1_D_lon_periodic`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:periodic</code> <code>axes:lon,lat</code> <code>stencil width 2</code></p>

Source: [`grids/latlon/rules/upwind1_D_lon_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/latlon/rules/upwind1_D_lon_periodic.esm)

**Rewrites** `D(f)/Dlon` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "lon"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lon, lat]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLON, 1:NLAT] = upwind1_D_lon_interior&lt;f=f&gt;, [1:1, 1:NLAT] = sum[i, j] ((f[1, j] - f[NLON, j]) / dlon_deg) where {i in 1:1, j in lat})</div>

First-order upwind (backward) first derivative along longitude with periodic (wrap) boundary on the latlon grid — the dissipative companion to grids/latlon/rules/central_D_lon_periodic for a closed zonal circle (global recipe: dlon_deg = 360/NLON; see grids/latlon/grid.esm), for zonal advection with non-negative angular velocity (flow toward increasing longitude). Matches D(f, wrt: lon), a PER-DEGREE derivative. Interior columns [2, NLON]: the imported backward-difference stencil (f[i,j] - f[i-1,j]) / dlon_deg. The single i=1 face column wraps explicitly to (f[1,j] - f[NLON,j]) / dlon_deg, following the cartesian upwind1_D_periodic / ESS pde_simulation fixture precedent for periodic boundaries, so the scheme is the same O(dlon_deg) backward difference at every column and conserves sum(f) exactly (the periodic wrap telescopes). Only one face column is special here — upwind needs a neighbor upstream (toward smaller lon), which every column but i=1 has — so the makearray tiles [lon, lat] with just two regions rather than the centered rule's three. The face value is a FULL-RANK aggregate (singleton output range on the lon axis, output_idx [i, j] over the lat axis) rather than the esm-spec 9.6.8 worked example's rank-reduced form, mirroring central_D_lon_periodic so every binding's simulator fills the (1 x NLAT) region without rank promotion. `dlon_deg` is a free name resolving to the consuming model's zonal spacing variable (consumer-supplied geometry contract); the rule itself never assumes the global extent, but composing it with a non-closing regional grid is a modeling error — use central_D_lon_zero_grad_bc there. Match scoping: the pattern is `wrt: lon` plus the constraint `where: {f: {shape: [lon, lat]}}` (esm-spec 9.6.1) at DEFAULT priority, so it fires only on the longitudinal derivative of a bare field declared over this grid's lon x lat index sets — not on any unrelated D(., wrt: lon) reusing the axis name — and, under import-edge renaming (esm-spec 9.7.7), the wrt literal and the where shape follow the renamed axis together, so two latlon instances coexist scoped to their own meshes. The constraint requires a bare shaped variable; a consumer differentiating a compound inline expression must bind it to a declared lon x lat observed first. Sign convention: upwind for flow toward +lon; a westward variant belongs in a separate rule.

**References**

- LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge University Press. Section 4.8 (upwind methods on periodic domains).

**Conformance fixtures**

- ast: [`upwind1_D_lon_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/upwind1_D_lon_periodic) — status `active`
- convergence: [`latlon_advection_zonal_upwind`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_zonal_upwind) — status `active`
- simulation: [`latlon_advection_zonal_upwind`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/latlon_advection_zonal_upwind) — status `active`

#### Convergence — case [`latlon_advection_zonal_upwind`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/latlon_advection_zonal_upwind)

Error norms read from the committed golden [`tests/conformance/convergence/latlon_advection_zonal_upwind/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/latlon_advection_zonal_upwind/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.424e-01 | — | 1.420e-01 | — |
| 32 | 7.414e-02 | 0.94 | 7.409e-02 | 0.94 |
| 64 | 3.781e-02 | 0.97 | 3.780e-02 | 0.97 |
| 128 | 1.909e-02 | 0.99 | 1.909e-02 | 0.99 |

Expected order: **1** (± 0.15).

![Convergence of latlon_advection_zonal_upwind (log-log)](/plots/convergence/latlon_advection_zonal_upwind.png)

