---
title: "cartesian_uniform_2d"
description: "Uniform 2-D Cartesian grid with consumer-supplied extent: NX by NY cell-centered points, cell centers (x_i, y_j) = (x0 + (i - 1/2) dx, y0 + (j - 1/2) dy)."
esd_kind: "grid"
families: ["cartesian", "finite_difference"]
ops: ["D", "D2", "Dxy", "bc", "blend", "grad_norm", "limiter", "reconstruct", "selector", "varcoeff_laplacian"]
bcs: ["dirichlet", "extrapolate_linear", "neumann", "periodic", "zero_flux", "zero_gradient"]
tags: ["esd:grid", "family:cartesian", "ndim:2", "spacing:uniform", "axis:x", "axis:y", "esd:stencil", "family:finite_difference", "grid:cartesian_uniform_2d", "op:D", "axes:x", "order:4", "stencil_width:5", "order:2", "stencil_width:3", "axes:y", "op:D2", "stencil_width:6", "order:6", "stencil_width:8", "stencil_width:7", "op:Dxy", "op:limiter", "order:1", "op:selector", "axes:x,y", "op:bc", "stencil_width:2", "op:grad_norm", "op:varcoeff_laplacian", "op:reconstruct", "order:3", "stencil_width:1", "op:blend", "order:5", "esd:rule", "bc:dirichlet", "spec:9.6.8", "bc:periodic", "bc:zero_gradient", "bc:neumann", "bc:extrapolate_linear", "bc:zero_flux"]
source: "grids/cartesian_uniform_2d/grid.esm"
---

Uniform 2-D Cartesian grid with consumer-supplied extent: NX by NY cell-centered points, cell centers (x_i, y_j) = (x0 + (i - 1/2) dx, y0 + (j - 1/2) dy). The `x` and `y` index sets are the iteration axes; NX and NY (the cell counts) remain the esm-spec 9.7 metaparameters, bound at the import edge or the loader API. The real-valued geometry is NOT baked into this file: `x0`/`dx` (x-axis) and `y0`/`dy` (y-axis) are free names in the coordinate template bodies (and in this grid's stencils/rules), resolving in the consuming component's scope at evaluation — the same keyed-factor contract grids/cartesian_uniform_1d uses for x0/dx and grids/mpas uses for areaCell/dvEdge. A consuming model MUST define all four as ordinary real-valued model variables: for a box [ax, bx] x [ay, by], x0 = ax, dx = (bx - ax)/NX, y0 = ay, dy = (by - ay)/NY. Spell each spacing as an observed dividing by the METAPARAMETER NAME (dx = {op: /, args: [bx - ax, NX]}), which esm-spec 9.7.6 substitutes as an integer literal at load, so a loader-API rebinding of NX/NY — a convergence sweep — keeps the spacing consistent automatically; a model that binds NX/NY at the import edge closes them before they reach the model's own scope and must spell the matching literals instead. BUILD-TIME SCOPE CAVEAT (AGENTS.md): the x_coord/y_coord free-name templates resolve in RUNTIME positions (rule bodies, equation RHS) but NOT in ic equations or §6.6.5 test references, which are evaluated build-time and cellwise with the consumer's parameters out of scope — there, coordinates must fold from literals plus the metaparameter names (x_i = ax + (i - 1/2)((bx - ax)/NX)).

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NX` | integer | `64` | Number of cells along x. |
| `NY` | integer | `64` | Number of cells along y. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>bdy_eno1_cells</code> <code>dx</code> <code>dy</code> <code>gB</code> <code>gL</code> <code>gR</code> <code>gT</code> <code>pi_blend</code> <code>qXL</code> <code>qXR</code> <code>qYL</code> <code>qYR</code> <code>upwind_sign</code> <code>weno_band_cells</code> <code>wind_u</code> <code>wind_v</code> <code>x0</code> <code>y0</code></p>

The geometry is consumer-supplied so the same files serve any **domain extent**, not just any resolution. For a domain [a, b] the consuming model defines `x0 = a` (a parameter) and `dx = (b − a)/N` (an observed whose expression divides by the metaparameter name `N`, so a loader-API rebinding of `N` — a convergence sweep — keeps `dx` consistent). A model that instead closes `N` at the import edge spells the matching literal (`{op: /, args: [1, 8]}` for N = 8 on the unit interval). See the extent proof `problems/heat_1d_zero_grad_nonunit.esm` (x ∈ [−1.5, 2.5], observed order 2.00).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `x` | interval | `NX` |
| `y` | interval | `NY` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `x_coord` | — | `sum[i] (x0 + (i - 0.5) * dx) where {i in x}` | Cell-center x-coordinates x_i = x0 + (i - 1/2) dx, shaped over the x index set. x0 and dx are free names supplied by the consuming model (see metadata.description for the contract). |
| `y_coord` | — | `sum[j] (y0 + (j - 0.5) * dy) where {j in y}` | Cell-center y-coordinates y_j = y0 + (j - 1/2) dy, shaped over the y index set. y0 and dy are free names supplied by the consuming model (see metadata.description for the contract). |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `central_D1_4_point`

<p class="rule-meta"><code>op:D</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D1_4_point.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D1_4_point.esm)

Parameters: `h`, `m2`, `m1`, `p1`, `p2`

<div class="math-block">(8 * (p1 - m1) - (p2 - m2)) / (12 * h)</div>

Pointwise fourth-order centered first derivative along one axis from the four neighbours m2, m1, p1, p2 (cells i-2, i-1, i+1, i+2) and the spacing h: (8 (p1 - m1) - (p2 - m2))/(12 h). Algebraically identical to fire_behavior's Select_4th, which forms it as the difference of two fourth-order face fluxes (7 p1 + 7 c - p2 - m1)/(12 h) - (7 c + 7 m1 - p1 - m2)/(12 h) (the cell value cancels). Used by the hybrid level-set stencils only to build the upwind SELECTOR of the WENO5 derivative (the wind-relative propagation sign), never as the derivative itself. Match-less helper.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `central_D1x_4_point`

<p class="rule-meta"><code>op:D</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D1x_4_point.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D1x_4_point.esm)

Parameters: `u`, `i`, `j`

<div class="math-block">central_D1_4_point&lt;h=dx, m1=u[i - 1, j], m2=u[i - 2, j], p1=u[i + 1, j], p2=u[i + 2, j]&gt;</div>

Pointwise fourth-order central derivative d u/dx at the cell (i, j) of cartesian_uniform_2d: stencils/central_D1_4_point bound to u at offsets -2, -1, +1, +2 along x and spacing dx. fire_behavior's Select_4th (dx, lfn(i, j), ...). u is the field, i and j the cell indices, bound by the caller's faq. It is not a derivative the level-set rules return: it only builds their upwind selectors (stencils/grad4_dot_wind, and the per-axis factor of signo_x). It reaches two cells either side along x.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `central_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D1x_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i + 1, j] - f[i - 1, j]) / (2 * dx)) where {i in 2:NX - 1, j in y}</div>

Interior-only 3-point centered first derivative along x on cartesian_uniform_2d: (f[i+1,j] - f[i-1,j]) / (2 dx) over i in [2, NX-1], all j. Match-less named template; the boundary-condition layer (rules/central_D1x_periodic.esm, rules/central_D1x_zero_grad_bc.esm) wraps it with the x-edge faces. `dx` is a free name resolving to the consuming model's x-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dx^2). The per-axis first derivative D(f,x) is lowered by its own rule; a consumer needing D(f,y) imports central_D1y_interior's rule instead.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D1y_4_point`

<p class="rule-meta"><code>op:D</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D1y_4_point.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D1y_4_point.esm)

Parameters: `u`, `i`, `j`

<div class="math-block">central_D1_4_point&lt;h=dy, m1=u[i, j - 1], m2=u[i, j - 2], p1=u[i, j + 1], p2=u[i, j + 2]&gt;</div>

Pointwise fourth-order central derivative d u/dy at the cell (i, j) of cartesian_uniform_2d: stencils/central_D1_4_point bound to u at offsets -2, -1, +1, +2 along y and spacing dy. fire_behavior's Select_4th (dy, lfn(i, j), ...). u is the field, i and j the cell indices, bound by the caller's faq. It is not a derivative the level-set rules return: it only builds their upwind selectors (stencils/grad4_dot_wind, and the per-axis factor of signo_y). It reaches two cells either side along y.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `central_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D1y_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i, j + 1] - f[i, j - 1]) / (2 * dy)) where {i in x, j in 2:NY - 1}</div>

Interior-only 3-point centered first derivative along y on cartesian_uniform_2d: (f[i,j+1] - f[i,j-1]) / (2 dy) over j in [2, NY-1], all i. Match-less named template; the boundary-condition layer (rules/central_D1y_periodic.esm, rules/central_D1y_zero_grad_bc.esm) wraps it with the y-edge faces. `dy` is a free name resolving to the consuming model's y-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dy^2). The per-axis first derivative D(f,y) is lowered by its own rule; a consumer needing D(f,x) imports central_D1x_interior's rule instead.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_4_dirichlet_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_4_dirichlet_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_4_dirichlet_face1.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(2560 * g + -4095 * f1 + 1848 * f2 + -378 * f3 + 72 * f4 + -7 * f5) / (756 * dx * dx)</div>

Fourth-order one-sided second-derivative closure along x for the cell nearest a Dirichlet (fixed-value) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_4_dirichlet_bc) applies it at the first interior x-column, tiled over all j. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_dirichlet_face1 (same Fornberg 1988 offset-node weights), rewritten for the x-axis so the free name is the x-spacing `dx`. On this cell-centered grid the x-wall sits half a cell outside the first center, so the wall value g and the five nearest cells along x sit at offsets {0, 1/2, 3/2, 5/2, 7/2, 9/2}*dx from the wall. The unique degree-5 polynomial through those six data, differentiated twice at the first center (offset 1/2 dx), gives D2 = (2560 g - 4095 f1 + 1848 f2 - 378 f3 + 72 f4 - 7 f5) / (756 dx^2), where f1..f5 are the five cells (at fixed j) in order of increasing distance from the x-wall and g is the wall value. Exact for polynomials up to degree 5, so the local truncation error is O(dx^4) — matching the interior 5-point stencil central_D2x_4_interior. The six weights sum to zero, so a steady linear-in-x profile carrying the wall value is discretely exact. Reused unchanged at the right x-wall by reflection (the second derivative is even under x -> -x), with f1..f5 mapped to the cells in increasing distance from the right wall and g bound to the right wall value. `dx` is a free name resolving to the consuming model's x-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_4_dirichlet_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_4_dirichlet_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_4_dirichlet_face2.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(-256 * g + 1323 * f1 + -2100 * f2 + 1134 * f3 + -108 * f4 + 7 * f5) / (756 * dx * dx)</div>

Fourth-order one-sided second-derivative closure along x for the SECOND cell from a Dirichlet (fixed-value) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_4_dirichlet_bc) applies it at the second interior x-column, tiled over all j, whose centered 5-point stencil would reach one cell outside the domain. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_dirichlet_face2, rewritten for the x-axis so the free name is `dx`. The same six data as central_D2x_4_dirichlet_face1 (the wall value g at offset 0 and the five nearest cells f1..f5 at offsets {1/2,3/2,5/2,7/2,9/2}*dx) fit a degree-5 polynomial; differentiating it twice at the SECOND center (offset 3/2 dx) gives D2 = (-256 g + 1323 f1 - 2100 f2 + 1134 f3 - 108 f4 + 7 f5) / (756 dx^2). Exact for polynomials up to degree 5, so the local truncation error is O(dx^4); the six weights sum to zero, so a steady linear-in-x profile is discretely exact. Reused unchanged at the right x-wall's second cell by reflection, with f1..f5 mapped to the cells in increasing distance from the right wall and g bound to the right wall value. `dx` is a free name resolving to the consuming model's x-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_4_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_4_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_4_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((-1 * f[i - 2, j] + 16 * f[i - 1, j] + -30 * f[i, j] + 16 * f[i + 1, j] + -1 * f[i + 2, j]) / (12 * dx * dx)) where {i in 3:NX - 2, j in y}</div>

Interior-only 5-point fourth-order centered second derivative along x on cartesian_uniform_2d: (-f[i-2,j] + 16 f[i-1,j] - 30 f[i,j] + 16 f[i+1,j] - f[i+2,j]) / (12 dx^2) over i in [3, NX-2], all j. The per-axis lift of the 1-D operator grids/cartesian_uniform_1d/stencils/central_D2_4_interior with the y axis a passive free index (the same coefficients -1,16,-30,16,-1 / 12, applied along x at fixed j). Match-less named template; the boundary-condition layer (rules/central_D2x_4_dirichlet_bc.esm and siblings) wraps it with the x-edge fourth-order face closures. `dx` is a free name resolving to the consuming model's x-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dx^4). The 2-D fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis fourth-order rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_4_neumann_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_4_neumann_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_4_neumann_face1.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(-1.92e4 * q * dx + -1.5285e4 * f1 + 7544 * f2 + 1.2546e4 * f3 + -5784 * f4 + 979 * f5) / (2.0268e4 * dx * dx)</div>

Fourth-order one-sided second-derivative closure along x for the cell nearest a Neumann (fixed-flux) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_4_neumann_bc) applies it at the first interior x-column, tiled over all j, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_neumann_face1, rewritten for the x-axis so the free name is `dx`. The datum is the wall derivative q = df/dx at the x-wall (offset 0) together with the five nearest cells f1..f5 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2}*dx. The unique degree-5 Hermite-type polynomial matching those six data differentiated twice at the first center (offset 1/2 dx) gives D2 = (-19200 q dx - 15285 f1 + 7544 f2 + 12546 f3 - 5784 f4 + 979 f5) / (20268 dx^2). Exact for polynomials up to degree 5, so the local truncation error is O(dx^4) — matching the interior 5-point stencil central_D2x_4_interior. The five value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-x with the matching wall flux is discretely exact. Written for a LEFT (low-index) x-wall where q = +df/dx points into the domain; reused unchanged at the right x-wall by reflection (the second derivative is even under x -> -x while df/dx is odd), so the rule binds q -> -qXR there with f1..f5 mapped to the cells in increasing distance from the right wall. `dx` is a free name resolving to the consuming model's x-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_4_neumann_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_4_neumann_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_4_neumann_face2.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(1920 * q * dx + 2.6019e4 * f1 + -5.21e4 * f2 + 2.8134e4 * f3 + -2124 * f4 + 71 * f5) / (2.0268e4 * dx * dx)</div>

Fourth-order one-sided second-derivative closure along x for the SECOND cell from a Neumann (fixed-flux) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_4_neumann_bc) applies it at the second interior x-column, tiled over all j, whose centered 5-point stencil would reach one cell outside the domain. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_neumann_face2, rewritten for the x-axis so the free name is `dx`. The same six data as central_D2x_4_neumann_face1 (the wall derivative q = df/dx at offset 0 and the five nearest cells f1..f5 at offsets {1/2,3/2,5/2,7/2,9/2}*dx) determine the degree-5 Hermite-type polynomial; differentiating it twice at the SECOND center (offset 3/2 dx) gives D2 = (1920 q dx + 26019 f1 - 52100 f2 + 28134 f3 - 2124 f4 + 71 f5) / (20268 dx^2). Exact for polynomials up to degree 5, so the local truncation error is O(dx^4); the five value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-x with the matching wall flux is discretely exact. Written for a LEFT (low-index) x-wall where q = +df/dx points into the domain; reused unchanged at the right x-wall's second cell by reflection (the second derivative even, df/dx odd), so the rule binds q -> -qXR there with f1..f5 mapped to the cells in increasing distance from the right wall. `dx` is a free name resolving to the consuming model's x-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_dirichlet_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face1.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(8.6016e4 * g + -1.36136e5 * f1 + 5.6628e4 * f2 + -3861 * f3 + -5720 * f4 + 4290 * f5 + -1404 * f6 + 187 * f7) / (2.574e4 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the FIRST cell from a Dirichlet (fixed-value) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_dirichlet_bc) applies it at the first interior column, tiled over all j. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face1 (same Fornberg 1988 offset-node weights), rewritten for the x-axis so the free name is `dx`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx fit a degree-7 polynomial; differentiating it twice at the first center (offset 1/2 dx) gives D2 = (86016 g - 136136 f1 + 56628 f2 - 3861 f3 - 5720 f4 + 4290 f5 - 1404 f6 + 187 f7) / (25740 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The eight weights sum to zero, so a steady linear-in-x profile carrying the wall value is discretely exact. Reused unchanged at the far x-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_dirichlet_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face2.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-2048 * g + 1.1349e4 * f1 + -1.7745e4 * f2 + 8073 * f3 + 1170 * f4 + -1105 * f5 + 351 * f6 + -45 * f7) / (7020 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the SECOND cell from a Dirichlet (fixed-value) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_dirichlet_bc) applies it at the second interior column, tiled over all j. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face2 (same Fornberg 1988 offset-node weights), rewritten for the x-axis so the free name is `dx`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx fit a degree-7 polynomial; differentiating it twice at the second center (offset 3/2 dx) gives D2 = (-2048 g + 11349 f1 - 17745 f2 + 8073 f3 + 1170 f4 - 1105 f5 + 351 f6 - 45 f7) / (7020 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The eight weights sum to zero, so a steady linear-in-x profile carrying the wall value is discretely exact. Reused unchanged at the far x-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_dirichlet_face3`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_dirichlet_face3.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(4096 * g + -1.7589e4 * f1 + 1.21836e5 * f2 + -2.16216e5 * f3 + 1.2012e5 * f4 + -1.3585e4 * f5 + 1404 * f6 + -66 * f7) / (7.722e4 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the THIRD cell from a Dirichlet (fixed-value) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_dirichlet_bc) applies it at the third interior column, tiled over all j. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face3 (same Fornberg 1988 offset-node weights), rewritten for the x-axis so the free name is `dx`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx fit a degree-7 polynomial; differentiating it twice at the third center (offset 5/2 dx) gives D2 = (4096 g - 17589 f1 + 121836 f2 - 216216 f3 + 120120 f4 - 13585 f5 + 1404 f6 - 66 f7) / (77220 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The eight weights sum to zero, so a steady linear-in-x profile carrying the wall value is discretely exact. Reused unchanged at the far x-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((2 * f[i - 3, j] + -27 * f[i - 2, j] + 270 * f[i - 1, j] + -490 * f[i, j] + 270 * f[i + 1, j] + -27 * f[i + 2, j] + 2 * f[i + 3, j]) / (180 * dx * dx)) where {i in 4:NX - 3, j in y}</div>

Interior-only 7-point sixth-order centered second derivative along x on cartesian_uniform_2d: (2 f[i-3,j] - 27 f[..] + 270 f[..] - 490 f[i,j] + 270 f[..] - 27 f[..] + 2 f[..]) / (180 dx^2) over x in [4, NX-3], all y. The per-axis lift of the 1-D operator grids/cartesian_uniform_1d/stencils/central_D2_6_interior with the y axis a passive free index (the same coefficients 2,-27,270,-490,270,-27,2 / 180, applied along x at fixed y). Match-less named template; the boundary-condition layer (rules/central_D2x_6_dirichlet_bc.esm and siblings) wraps it with the x-edge sixth-order face closures. `dx` is a free name resolving to the consuming model's x-spacing variable (grids/cartesian_uniform_2d/grid.esm). O(dx^6). The 2-D sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis sixth-order rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_neumann_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face1.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-1.354752e7 * q * dx + -4.382308e6 * f1 + -1.8097596e7 * f2 + 4.5297765e7 * f3 + -3.595516e7 * f4 + 1.735677e7 * f5 + -4.804812e6 * f6 + 5.85341e5 * f7) / (1.585242e7 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the FIRST cell from a Neumann (fixed-flux) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_neumann_bc) applies it at the first interior column, tiled over all j, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face1, rewritten for the x-axis so the free name is `dx`. The datum is the wall derivative q = df/dx at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the first center (offset 1/2 dx) gives D2 = (-13547520 q dx - 4382308 f1 - 18097596 f2 + 45297765 f3 - 35955160 f4 + 17356770 f5 - 4804812 f6 + 585341 f7) / (15852420 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-x with the matching wall flux is discretely exact. Written for a LEFT (low-index) x-wall where q = +df/dx points into the domain; reused unchanged at the far x-wall by reflection (the second derivative even, df/dx odd), so the rule binds q -> -qXR there with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_neumann_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face2.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(1.18272e6 * q * dx + 1.8691149e7 * f1 + -3.5446775e7 * f2 + 1.4068125e7 * f3 + 5.47347e6 * f4 + -3.779905e6 * f5 + 1.136601e6 * f6 + -1.42665e5 * f7) / (1.585242e7 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the SECOND cell from a Neumann (fixed-flux) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_neumann_bc) applies it at the second interior column, tiled over all j, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face2, rewritten for the x-axis so the free name is `dx`. The datum is the wall derivative q = df/dx at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the second center (offset 3/2 dx) gives D2 = (1182720 q dx + 18691149 f1 - 35446775 f2 + 14068125 f3 + 5473470 f4 - 3779905 f5 + 1136601 f6 - 142665 f7) / (15852420 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-x with the matching wall flux is discretely exact. Written for a LEFT (low-index) x-wall where q = +df/dx points into the domain; reused unchanged at the far x-wall by reflection (the second derivative even, df/dx odd), so the rule binds q -> -qXR there with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_6_neumann_face3`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_6_neumann_face3.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-2.1504e5 * q * dx + -2.349569e6 * f1 + 2.4170756e7 * f2 + -4.363002e7 * f3 + 2.414452e7 * f4 + -2.555285e6 * f5 + 2.25684e5 * f6 + -6086 * f7) / (1.585242e7 * dx * dx)</div>

Sixth-order one-sided second-derivative closure along x for the THIRD cell from a Neumann (fixed-flux) x-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2x_6_neumann_bc) applies it at the third interior column, tiled over all j, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face3, rewritten for the x-axis so the free name is `dx`. The datum is the wall derivative q = df/dx at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed j) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dx. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the third center (offset 5/2 dx) gives D2 = (-215040 q dx - 2349569 f1 + 24170756 f2 - 43630020 f3 + 24144520 f4 - 2555285 f5 + 225684 f6 - 6086 f7) / (15852420 dx^2). Exact for polynomials up to degree 7, so the local truncation error is O(dx^6) — matching the interior 7-point stencil central_D2x_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-x with the matching wall flux is discretely exact. Written for a LEFT (low-index) x-wall where q = +df/dx points into the domain; reused unchanged at the far x-wall by reflection (the second derivative even, df/dx odd), so the rule binds q -> -qXR there with f1..f7 mapped to the cells in increasing distance from that wall. `dx` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2x_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2x_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i + 1, j] + -2 * f[i, j] + f[i - 1, j]) / (dx * dx)) where {i in 2:NX - 1, j in y}</div>

Interior-only 3-point centered second derivative along x on cartesian_uniform_2d: (f[i+1,j] - 2 f[i,j] + f[i-1,j]) / dx^2 over i in [2, NX-1], all j. Match-less named template; the boundary-condition layer (rules/central_D2x_neumann_bc.esm) wraps it with the x-edge faces. `dx` is a free name resolving to the consuming model's x-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dx^2). The 2-D Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_4_dirichlet_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_4_dirichlet_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_4_dirichlet_face1.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(2560 * g + -4095 * f1 + 1848 * f2 + -378 * f3 + 72 * f4 + -7 * f5) / (756 * dy * dy)</div>

Fourth-order one-sided second-derivative closure along y for the cell nearest a Dirichlet (fixed-value) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_4_dirichlet_bc) applies it at the first interior y-row, tiled over all i. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_dirichlet_face1 (same Fornberg 1988 offset-node weights), rewritten for the y-axis so the free name is the y-spacing `dy`. On this cell-centered grid the y-wall sits half a cell outside the first center, so the wall value g and the five nearest cells along y sit at offsets {0, 1/2, 3/2, 5/2, 7/2, 9/2}*dy from the wall. The unique degree-5 polynomial through those six data, differentiated twice at the first center (offset 1/2 dy), gives D2 = (2560 g - 4095 f1 + 1848 f2 - 378 f3 + 72 f4 - 7 f5) / (756 dy^2), where f1..f5 are the five cells (at fixed i) in order of increasing distance from the y-wall and g is the wall value. Exact for polynomials up to degree 5, so the local truncation error is O(dy^4) — matching the interior 5-point stencil central_D2y_4_interior. The six weights sum to zero, so a steady linear-in-y profile carrying the wall value is discretely exact. Reused unchanged at the top y-wall by reflection (the second derivative is even under y -> -y), with f1..f5 mapped to the cells in increasing distance from the top wall and g bound to the top wall value. `dy` is a free name resolving to the consuming model's y-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_4_dirichlet_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_4_dirichlet_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_4_dirichlet_face2.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(-256 * g + 1323 * f1 + -2100 * f2 + 1134 * f3 + -108 * f4 + 7 * f5) / (756 * dy * dy)</div>

Fourth-order one-sided second-derivative closure along y for the SECOND cell from a Dirichlet (fixed-value) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_4_dirichlet_bc) applies it at the second interior y-row, tiled over all i, whose centered 5-point stencil would reach one cell outside the domain. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_dirichlet_face2, rewritten for the y-axis so the free name is `dy`. The same six data as central_D2y_4_dirichlet_face1 (the wall value g at offset 0 and the five nearest cells f1..f5 at offsets {1/2,3/2,5/2,7/2,9/2}*dy) fit a degree-5 polynomial; differentiating it twice at the SECOND center (offset 3/2 dy) gives D2 = (-256 g + 1323 f1 - 2100 f2 + 1134 f3 - 108 f4 + 7 f5) / (756 dy^2). Exact for polynomials up to degree 5, so the local truncation error is O(dy^4); the six weights sum to zero, so a steady linear-in-y profile is discretely exact. Reused unchanged at the top y-wall's second cell by reflection, with f1..f5 mapped to the cells in increasing distance from the top wall and g bound to the top wall value. `dy` is a free name resolving to the consuming model's y-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_4_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_4_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_4_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((-1 * f[i, j - 2] + 16 * f[i, j - 1] + -30 * f[i, j] + 16 * f[i, j + 1] + -1 * f[i, j + 2]) / (12 * dy * dy)) where {i in x, j in 3:NY - 2}</div>

Interior-only 5-point fourth-order centered second derivative along y on cartesian_uniform_2d: (-f[i,j-2] + 16 f[i,j-1] - 30 f[i,j] + 16 f[i,j+1] - f[i,j+2]) / (12 dy^2) over j in [3, NY-2], all i. The per-axis lift of the 1-D operator grids/cartesian_uniform_1d/stencils/central_D2_4_interior with the x axis a passive free index (the same coefficients -1,16,-30,16,-1 / 12, applied along y at fixed i). Match-less named template; the boundary-condition layer (rules/central_D2y_4_dirichlet_bc.esm and siblings) wraps it with the y-edge fourth-order face closures. `dy` is a free name resolving to the consuming model's y-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dy^4). The 2-D fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis fourth-order rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_4_neumann_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_4_neumann_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_4_neumann_face1.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(-1.92e4 * q * dy + -1.5285e4 * f1 + 7544 * f2 + 1.2546e4 * f3 + -5784 * f4 + 979 * f5) / (2.0268e4 * dy * dy)</div>

Fourth-order one-sided second-derivative closure along y for the cell nearest a Neumann (fixed-flux) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_4_neumann_bc) applies it at the first interior y-row, tiled over all i, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_neumann_face1, rewritten for the y-axis so the free name is `dy`. The datum is the wall derivative q = df/dy at the y-wall (offset 0) together with the five nearest cells f1..f5 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2}*dy. The unique degree-5 Hermite-type polynomial matching those six data differentiated twice at the first center (offset 1/2 dy) gives D2 = (-19200 q dy - 15285 f1 + 7544 f2 + 12546 f3 - 5784 f4 + 979 f5) / (20268 dy^2). Exact for polynomials up to degree 5, so the local truncation error is O(dy^4) — matching the interior 5-point stencil central_D2y_4_interior. The five value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-y with the matching wall flux is discretely exact. Written for a BOTTOM (low-index) y-wall where q = +df/dy points into the domain; reused unchanged at the top y-wall by reflection (the second derivative is even under y -> -y while df/dy is odd), so the rule binds q -> -qYR there with f1..f5 mapped to the cells in increasing distance from the top wall. `dy` is a free name resolving to the consuming model's y-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_4_neumann_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>stencil width 6</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_4_neumann_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_4_neumann_face2.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`

<div class="math-block">(1920 * q * dy + 2.6019e4 * f1 + -5.21e4 * f2 + 2.8134e4 * f3 + -2124 * f4 + 71 * f5) / (2.0268e4 * dy * dy)</div>

Fourth-order one-sided second-derivative closure along y for the SECOND cell from a Neumann (fixed-flux) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_4_neumann_bc) applies it at the second interior y-row, tiled over all i, whose centered 5-point stencil would reach one cell outside the domain. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_4_neumann_face2, rewritten for the y-axis so the free name is `dy`. The same six data as central_D2y_4_neumann_face1 (the wall derivative q = df/dy at offset 0 and the five nearest cells f1..f5 at offsets {1/2,3/2,5/2,7/2,9/2}*dy) determine the degree-5 Hermite-type polynomial; differentiating it twice at the SECOND center (offset 3/2 dy) gives D2 = (1920 q dy + 26019 f1 - 52100 f2 + 28134 f3 - 2124 f4 + 71 f5) / (20268 dy^2). Exact for polynomials up to degree 5, so the local truncation error is O(dy^4); the five value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-y with the matching wall flux is discretely exact. Written for a BOTTOM (low-index) y-wall where q = +df/dy points into the domain; reused unchanged at the top y-wall's second cell by reflection (the second derivative even, df/dy odd), so the rule binds q -> -qYR there with f1..f5 mapped to the cells in increasing distance from the top wall. `dy` is a free name resolving to the consuming model's y-spacing variable (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_dirichlet_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face1.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(8.6016e4 * g + -1.36136e5 * f1 + 5.6628e4 * f2 + -3861 * f3 + -5720 * f4 + 4290 * f5 + -1404 * f6 + 187 * f7) / (2.574e4 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the FIRST cell from a Dirichlet (fixed-value) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_dirichlet_bc) applies it at the first interior row, tiled over all i. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face1 (same Fornberg 1988 offset-node weights), rewritten for the y-axis so the free name is `dy`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy fit a degree-7 polynomial; differentiating it twice at the first center (offset 1/2 dy) gives D2 = (86016 g - 136136 f1 + 56628 f2 - 3861 f3 - 5720 f4 + 4290 f5 - 1404 f6 + 187 f7) / (25740 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The eight weights sum to zero, so a steady linear-in-y profile carrying the wall value is discretely exact. Reused unchanged at the far y-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_dirichlet_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face2.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-2048 * g + 1.1349e4 * f1 + -1.7745e4 * f2 + 8073 * f3 + 1170 * f4 + -1105 * f5 + 351 * f6 + -45 * f7) / (7020 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the SECOND cell from a Dirichlet (fixed-value) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_dirichlet_bc) applies it at the second interior row, tiled over all i. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face2 (same Fornberg 1988 offset-node weights), rewritten for the y-axis so the free name is `dy`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy fit a degree-7 polynomial; differentiating it twice at the second center (offset 3/2 dy) gives D2 = (-2048 g + 11349 f1 - 17745 f2 + 8073 f3 + 1170 f4 - 1105 f5 + 351 f6 - 45 f7) / (7020 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The eight weights sum to zero, so a steady linear-in-y profile carrying the wall value is discretely exact. Reused unchanged at the far y-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_dirichlet_face3`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_dirichlet_face3.esm)

Parameters: `g`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(4096 * g + -1.7589e4 * f1 + 1.21836e5 * f2 + -2.16216e5 * f3 + 1.2012e5 * f4 + -1.3585e4 * f5 + 1404 * f6 + -66 * f7) / (7.722e4 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the THIRD cell from a Dirichlet (fixed-value) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_dirichlet_bc) applies it at the third interior row, tiled over all i. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_dirichlet_face3 (same Fornberg 1988 offset-node weights), rewritten for the y-axis so the free name is `dy`. The wall value g at offset 0 and the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy fit a degree-7 polynomial; differentiating it twice at the third center (offset 5/2 dy) gives D2 = (4096 g - 17589 f1 + 121836 f2 - 216216 f3 + 120120 f4 - 13585 f5 + 1404 f6 - 66 f7) / (77220 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The eight weights sum to zero, so a steady linear-in-y profile carrying the wall value is discretely exact. Reused unchanged at the far y-wall by reflection (the second derivative is even), with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((2 * f[i, j - 3] + -27 * f[i, j - 2] + 270 * f[i, j - 1] + -490 * f[i, j] + 270 * f[i, j + 1] + -27 * f[i, j + 2] + 2 * f[i, j + 3]) / (180 * dy * dy)) where {i in x, j in 4:NY - 3}</div>

Interior-only 7-point sixth-order centered second derivative along y on cartesian_uniform_2d: (2 f[i,j-3] - 27 f[..] + 270 f[..] - 490 f[i,j] + 270 f[..] - 27 f[..] + 2 f[..]) / (180 dy^2) over y in [4, NY-3], all x. The per-axis lift of the 1-D operator grids/cartesian_uniform_1d/stencils/central_D2_6_interior with the x axis a passive free index (the same coefficients 2,-27,270,-490,270,-27,2 / 180, applied along y at fixed x). Match-less named template; the boundary-condition layer (rules/central_D2y_6_dirichlet_bc.esm and siblings) wraps it with the y-edge sixth-order face closures. `dy` is a free name resolving to the consuming model's y-spacing variable (grids/cartesian_uniform_2d/grid.esm). O(dy^6). The 2-D sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis sixth-order rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_neumann_face1`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face1.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-1.354752e7 * q * dy + -4.382308e6 * f1 + -1.8097596e7 * f2 + 4.5297765e7 * f3 + -3.595516e7 * f4 + 1.735677e7 * f5 + -4.804812e6 * f6 + 5.85341e5 * f7) / (1.585242e7 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the FIRST cell from a Neumann (fixed-flux) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_neumann_bc) applies it at the first interior row, tiled over all i, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face1, rewritten for the y-axis so the free name is `dy`. The datum is the wall derivative q = df/dy at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the first center (offset 1/2 dy) gives D2 = (-13547520 q dy - 4382308 f1 - 18097596 f2 + 45297765 f3 - 35955160 f4 + 17356770 f5 - 4804812 f6 + 585341 f7) / (15852420 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-y with the matching wall flux is discretely exact. Written for a BOTTOM (low-index) y-wall where q = +df/dy points into the domain; reused unchanged at the far y-wall by reflection (the second derivative even, df/dy odd), so the rule binds q -> -qYR there with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_neumann_face2`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face2.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(1.18272e6 * q * dy + 1.8691149e7 * f1 + -3.5446775e7 * f2 + 1.4068125e7 * f3 + 5.47347e6 * f4 + -3.779905e6 * f5 + 1.136601e6 * f6 + -1.42665e5 * f7) / (1.585242e7 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the SECOND cell from a Neumann (fixed-flux) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_neumann_bc) applies it at the second interior row, tiled over all i, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face2, rewritten for the y-axis so the free name is `dy`. The datum is the wall derivative q = df/dy at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the second center (offset 3/2 dy) gives D2 = (1182720 q dy + 18691149 f1 - 35446775 f2 + 14068125 f3 + 5473470 f4 - 3779905 f5 + 1136601 f6 - 142665 f7) / (15852420 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-y with the matching wall flux is discretely exact. Written for a BOTTOM (low-index) y-wall where q = +df/dy points into the domain; reused unchanged at the far y-wall by reflection (the second derivative even, df/dy odd), so the rule binds q -> -qYR there with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_6_neumann_face3`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>stencil width 8</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_6_neumann_face3.esm)

Parameters: `q`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`

<div class="math-block">(-2.1504e5 * q * dy + -2.349569e6 * f1 + 2.4170756e7 * f2 + -4.363002e7 * f3 + 2.414452e7 * f4 + -2.555285e6 * f5 + 2.25684e5 * f6 + -6086 * f7) / (1.585242e7 * dy * dy)</div>

Sixth-order one-sided second-derivative closure along y for the THIRD cell from a Neumann (fixed-flux) y-wall on cartesian_uniform_2d, cell-centered. Match-less named scalar template; the boundary-condition rule (rules/central_D2y_6_neumann_bc) applies it at the third interior row, tiled over all i, against a prescribed wall gradient. Identical algebra to the 1-D closure grids/cartesian_uniform_1d/stencils/central_D2_6_neumann_face3, rewritten for the y-axis so the free name is `dy`. The datum is the wall derivative q = df/dy at the wall (offset 0) together with the seven nearest cells f1..f7 (fixed i) at offsets {1/2,3/2,5/2,7/2,9/2,11/2,13/2}*dy. The degree-7 Hermite-type polynomial matching those eight data differentiated twice at the third center (offset 5/2 dy) gives D2 = (-215040 q dy - 2349569 f1 + 24170756 f2 - 43630020 f3 + 24144520 f4 - 2555285 f5 + 225684 f6 - 6086 f7) / (15852420 dy^2). Exact for polynomials up to degree 7, so the local truncation error is O(dy^6) — matching the interior 7-point stencil central_D2y_6_interior. The seven value-weights sum to zero and the derivative term supplies the exact wall slope, so a steady ramp-in-y with the matching wall flux is discretely exact. Written for a BOTTOM (low-index) y-wall where q = +df/dy points into the domain; reused unchanged at the far y-wall by reflection (the second derivative even, df/dy odd), so the rule binds q -> -qYR there with f1..f7 mapped to the cells in increasing distance from that wall. `dy` is a free name (grids/cartesian_uniform_2d/grid.esm).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_D2y_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_D2y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_D2y_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i, j + 1] + -2 * f[i, j] + f[i, j - 1]) / (dy * dy)) where {i in x, j in 2:NY - 1}</div>

Interior-only 3-point centered second derivative along y on cartesian_uniform_2d: (f[i,j+1] - 2 f[i,j] + f[i,j-1]) / dy^2 over j in [2, NY-1], all i. Match-less named template; the boundary-condition layer (rules/central_D2y_neumann_bc.esm) wraps it with the y-edge faces. `dy` is a free name resolving to the consuming model's y-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dy^2). The 2-D Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y), each half lowered by its own per-axis rule.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

### `central_Dxy_interior`

<p class="rule-meta"><code>op:Dxy</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/central_Dxy_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/central_Dxy_interior.esm)

Parameters: `f`

<div class="math-block">sum[i, j] ((f[i + 1, j + 1] - f[i + 1, j - 1] - (f[i - 1, j + 1] - f[i - 1, j - 1])) / (4 * dx * dy)) where {i in 2:NX - 1, j in 2:NY - 1}</div>

Interior-only 4-point centered MIXED second derivative d^2f/dx dy along the x and y axes on cartesian_uniform_2d: (f[i+1,j+1] - f[i+1,j-1] - f[i-1,j+1] + f[i-1,j-1]) / (4 dx dy) over i in [2, NX-1] and j in [2, NY-1] (the doubly-interior block, both axes one cell in from either wall so all four corner samples f[i+-1,j+-1] exist). Match-less named template; the boundary-condition layer (rules/central_Dxy_periodic.esm) wraps it with the x- and y-edge and corner faces. This scheme is the tensor product of the two centered first differences (central_D1x_interior applied to central_D1y_interior), so it inherits their O(dx^2 + dy^2) accuracy; it is the off-diagonal term of a full anisotropic diffusion tensor K = [[Kxx, Kxy], [Kxy, Kyy]], the 2 Kxy d^2f/dx dy contribution that the diagonal central_D2x/central_D2y stencils do not supply. `dx` and `dy` are free names resolving to the consuming model's x- and y-spacing variables (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). O(dx^2 + dy^2).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (mixed partial derivatives as products of one-dimensional centered differences).

### `eno1_minmod`

<p class="rule-meta"><code>op:limiter</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/eno1_minmod.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/eno1_minmod.esm)

Parameters: `dl`, `dr`

<div class="math-block">max(0, min(dl, dr)) + min(0, max(dl, dr))</div>

Pointwise first-order ENO (minmod) selection of a one-sided derivative from the backward difference dl = (u[i]-u[i-1])/h and the forward difference dr = (u[i+1]-u[i])/h: the smaller-magnitude one when both have the same sign, zero when they disagree in sign (the derivative crosses zero inside the cell) or either vanishes. Written without sign() or a branch as max(0, min(dl, dr)) + min(0, max(dl, dr)): for dl, dr > 0 this is min(dl, dr); for dl, dr < 0 it is max(dl, dr) (the less negative); for mixed signs or a zero operand both terms vanish -- identical to fire_behavior's Select_eno (if dl*dr > 0 then sign(dl) min(|dl|,|dr|) else 0). First order O(h) on smooth data (a one-sided difference), exact on a linear field (dl = dr = a gives a), total-variation-diminishing through kinks. Match-less helper bound by the weno5_eno1_* and weno3_eno1_* stencils and rules of this grid.

**References**

- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `eno1_minmod_D1x`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/eno1_minmod_D1x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/eno1_minmod_D1x.esm)

Parameters: `u`, `i`, `j`

<div class="math-block">eno1_minmod&lt;dl=(u[i, j] - u[i - 1, j]) / dx, dr=(u[i + 1, j] - u[i, j]) / dx&gt;</div>

Pointwise first-order ENO derivative d u/dx at the cell (i, j) of cartesian_uniform_2d: stencils/eno1_minmod of the backward difference dl = (u[i,j] - u[i-1,j])/dx and the forward difference dr = (u[i+1,j] - u[i,j])/dx. fire_behavior's Select_eno (diffLx, diffRx) (respectively diffLy, diffRy). u is the field, i and j the cell indices, bound by the caller's faq (a pass-through of its loop symbols, as in stencils/eno1_wall_band); dx is the grid's free-name spacing. It reaches one cell either side along x, so a caller evaluates it only where both neighbours exist -- the wall lines use a ghost instead (rules/weno5_eno1_D1x_extrapolate_frame). The one spelling of the ENO1 branch shared by the weno5_eno1_*D1x_interior stencils and the ENO1 strips of the D1x frame.

**References**

- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `eno1_minmod_D1y`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/eno1_minmod_D1y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/eno1_minmod_D1y.esm)

Parameters: `u`, `i`, `j`

<div class="math-block">eno1_minmod&lt;dl=(u[i, j] - u[i, j - 1]) / dy, dr=(u[i, j + 1] - u[i, j]) / dy&gt;</div>

Pointwise first-order ENO derivative d u/dy at the cell (i, j) of cartesian_uniform_2d: stencils/eno1_minmod of the backward difference dl = (u[i,j] - u[i,j-1])/dy and the forward difference dr = (u[i,j+1] - u[i,j])/dy. fire_behavior's Select_eno (diffLx, diffRx) (respectively diffLy, diffRy). u is the field, i and j the cell indices, bound by the caller's faq (a pass-through of its loop symbols, as in stencils/eno1_wall_band); dy is the grid's free-name spacing. It reaches one cell either side along y, so a caller evaluates it only where both neighbours exist -- the wall lines use a ghost instead (rules/weno5_eno1_D1y_extrapolate_frame). The one spelling of the ENO1 branch shared by the weno5_eno1_*D1y_interior stencils and the ENO1 strips of the D1y frame.

**References**

- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `eno1_wall_band`

<p class="rule-meta"><code>op:selector</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/eno1_wall_band.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/eno1_wall_band.esm)

Parameters: `i`, `j`

<div class="math-block">i &lt; 1 + bdy_eno1_cells or i &gt; NX - bdy_eno1_cells or j &lt; 1 + bdy_eno1_cells or j &gt; NY - bdy_eno1_cells</div>

Pointwise predicate: true when the cell (i, j) of cartesian_uniform_2d lies within bdy_eno1_cells cells of any wall, i < 1 + bdy_eno1_cells or i > NX - bdy_eno1_cells or j < 1 + bdy_eno1_cells or j > NY - bdy_eno1_cells. This is the wall band of NCAR fire_behavior's Calc_tend_ls, where the level-set gradient is ALWAYS the first-order ENO1 minmod (stencils/eno1_minmod), whatever the fire_upwinding option: the test precedes the `select case (fire_upwinding)`, so every upwinding scheme is only an interior scheme there. bdy_eno1_cells is a free name resolved in the consuming model's scope (fire_behavior's BDY_ENO1 = 10, a compile-time constant of level_set_mod.F90); i and j are the cell's indices, bound by the caller's faq. The wall band depends on BOTH axes, so a per-axis derivative along x switches to ENO1 near the y walls too. Match-less helper bound by the godunov_eno1_D1, weno5_eno1_* and weno3_eno1_* stencils of this grid.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (the `i < ids + BDY_ENO1 .or. i > ide - BDY_ENO1 .or. j < jds + BDY_ENO1 .or. j > jde - BDY_ENO1` switch). The transcribed reference implementation.

### `extrapolate_linear_ghost`

<p class="rule-meta"><code>op:bc</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/extrapolate_linear_ghost.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/extrapolate_linear_ghost.esm)

Parameters: `a`, `b`

<div class="math-block">max((2 * a - b), a, b)</div>

Pointwise one-cell halo value outside a domain wall from the two innermost cells a (the wall cell) and b (its inward neighbour): the linear extrapolation 2 a - b, guarded from below by max(2 a - b, a, b) -- fire_behavior's Extrapol_var_at_bdys with BIAS = 1 ((1 - BIAS)(2a - b) + BIAS max(2a - b, a, b)). For a level-set (signed-distance) field increasing toward the wall (a > b, the fire inside the domain) the guard is inactive and the ghost is the pure linear extrapolation; when the field decreases toward the wall the ghost is clamped to the larger of the two inner values, so the one-sided difference across the wall is never steeper than the field's own slope and never changes sign spuriously. Exact on a linear field with a >= b. Match-less helper bound by every *_extrapolate_bc rule of this grid (the bc:extrapolate_linear closure).

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `godunov_eno1_D1`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/godunov_eno1_D1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/godunov_eno1_D1.esm)

Parameters: `dl`, `dr`, `i`, `j`

<div class="math-block">ifelse(eno1_wall_band&lt;i=i, j=j&gt;, eno1_minmod&lt;dl=dl, dr=dr&gt;, godunov_upwind_D1&lt;dl=dl, dr=dr&gt;)</div>

Pointwise per-axis quantity diff2 of NCAR fire_behavior's classical Godunov level-set option (fire_upwinding = 4) at cell (i, j), from the backward difference dl = D- u and the forward difference dr = D+ u: inside the wall band (stencils/eno1_wall_band, bdy_eno1_cells of any wall) it is the ENO1 minmod of the two (stencils/eno1_minmod), elsewhere the Godunov-upwind quantity max(dl, 0) - min(dr, 0) (stencils/godunov_upwind_D1) -- Calc_tend_ls's `use eno1 near domain boundaries` test, which precedes and overrides case (4). Both ifelse branches are three-point and evaluated on every cell (esm ifelse is not lazy); the caller keeps dl, dr in range. SIGN: the two branches disagree for a field decreasing along the axis -- ENO1 returns the signed (negative) derivative, Godunov the magnitude -D+ > 0 (fire_behavior note N44) -- exactly as the reference code does; their squares agree wherever the data are monotone or convex along the axis. Match-less helper bound by the godunov_upwind_D1x/D1y_interior stencils of this grid.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `godunov_norm_D2_interior`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/godunov_norm_D2_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/godunov_norm_D2_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (sqrt(max(((u[i, j] - u[i - 1, j]) / dx), 0)^2 + min(((u[i + 1, j] - u[i, j]) / dx), 0)^2 + max(((u[i, j] - u[i, j - 1]) / dy), 0)^2 + min(((u[i, j + 1] - u[i, j]) / dy), 0)^2)) where {i in 2:NX - 1, j in 2:NY - 1}</div>

Interior-only first-order Godunov (Rouy-Tourin / Osher-Sethian) upwind gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) on cartesian_uniform_2d over the doubly-interior block [2, NX-1] x [2, NY-1] where the 5-point stencil u[i,j], u[i-1,j], u[i+1,j], u[i,j-1], u[i,j+1] fits without reaching a boundary. Per-axis backward D-_x=(u[i,j]-u[i-1,j])/dx, D-_y=(u[i,j]-u[i,j-1])/dy and forward D+_x=(u[i+1,j]-u[i,j])/dx, D+_y=(u[i,j+1]-u[i,j])/dy; the norm is sqrt( max(D-_x,0)^2 + min(D+_x,0)^2 + max(D-_y,0)^2 + min(D+_y,0)^2 ), the per-axis entropy (upwind) fix of Osher-Sethian (1988) / Rouy-Tourin (1992). Match-less named template extracting the doubly-interior aggregate of grids/cartesian_uniform_2d/rules/godunov_norm_D2_periodic verbatim so a boundary-condition layer (rules/godunov_norm_D2_dirichlet_bc.esm) wraps it with prescribed-value edge and corner faces via apply_expression_template, exactly as central_D2x_dirichlet_bc reuses central_D2x_interior. EXACT on a linear field u = a x + b y + c (D-_x=D+_x=a, D-_y=D+_y=b, norm = sqrt(a^2+b^2) = |grad u|). First order O(dx + dy) on smooth data, uniform through the |grad| kinks. dx, dy are free names resolving to the consuming model's spacing variables (grids/cartesian_uniform_2d/grid.esm consumer-supplied geometry contract).

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. The max/min upwind (entropy) construction of \|grad phi\| for the level-set Hamilton-Jacobi term. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- Rouy, E., & Tourin, A. (1992). A viscosity solutions approach to shape-from-shading. SIAM Journal on Numerical Analysis, 29(3), 867-884. The Godunov Hamiltonian sqrt(max(D-,0)^2 + min(D+,0)^2) as the monotone, first-order-convergent discretization of the eikonal \|grad u\|. [`doi:10.1137/0729053`](https://doi.org/10.1137/0729053)

### `godunov_upwind_D1`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/godunov_upwind_D1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/godunov_upwind_D1.esm)

Parameters: `dl`, `dr`

<div class="math-block">max(dl, 0) - min(dr, 0)</div>

Pointwise Godunov-upwind one-sided derivative combination max(dl, 0) - min(dr, 0) from the backward difference dl and the forward difference dr along one axis: fire_behavior's diff2 of the classical Godunov option (Calc_tend_ls case 4, fire_upwinding = 4), used there as the per-axis 'derivative' that the front normal nvx = diff2x/sqrt(|grad psi|^2 + eps) is built from. NOTE (fire_behavior N44 in the EqWeFiC bug list): for a monotone field this equals |u_x| to first order -- a MAGNITUDE, not a signed derivative (dl > 0, dr > 0 gives dl; dl < 0, dr < 0 gives -dr > 0) -- so a consumer that wants the signed normal must supply the sign separately; the rules of this grid expose the quantity exactly as the reference code computes it. Exact on a linear field with positive slope; first order otherwise. Match-less helper of stencils/godunov_eno1_D1, which applies it only outside fire_behavior's ENO1 wall band.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `godunov_upwind_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/godunov_upwind_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/godunov_upwind_D1x_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (godunov_eno1_D1&lt;dl=(u[i, j] - u[i - 1, j]) / dx, dr=(u[i + 1, j] - u[i, j]) / dx, i=i, j=j&gt;) where {i in 2:NX - 1, j in y}</div>

Interior-only per-axis quantity diff2x of fire_behavior's classical Godunov level-set option (fire_upwinding = 4) on cartesian_uniform_2d over [2, NX-1] x [1, NY], where the three-point stencil along x fits: from the backward difference D- = (u[i,j] - u[i-1,j])/dx and the forward difference D+ = (u[i+1,j] - u[i,j])/dx (index shifted along x), each cell takes stencils/godunov_eno1_D1 -- the ENO1 minmod of D-, D+ within bdy_eno1_cells of any wall (fire_behavior's BDY_ENO1 band, which overrides every upwinding option), and the Godunov-upwind quantity max(D-, 0) - min(D+, 0) elsewhere. The Godunov quantity is the per-axis building block whose squares sum to the Godunov Hamiltonian only up to the cross term (max(D-,0) - min(D+,0))^2 = max(D-,0)^2 + min(D+,0)^2 - 2 max(D-,0) min(D+,0) >= the Rouy-Tourin Hamiltonian's per-axis term; fire_behavior uses it for the front NORMAL (nvx = diff2x / sqrt(grad^2 + eps)), while |grad psi| itself comes from godunov_norm_D2_*. First order; exact on a linear field with positive slope along x (both branches return the slope). See stencils/godunov_eno1_D1 for the sign caveat (outside the wall band, a magnitude for negative slopes). Free name bdy_eno1_cells (>= 1), resolved in the consuming model's scope like the grid geometry dx, dy.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `godunov_upwind_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/godunov_upwind_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/godunov_upwind_D1y_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (godunov_eno1_D1&lt;dl=(u[i, j] - u[i, j - 1]) / dy, dr=(u[i, j + 1] - u[i, j]) / dy, i=i, j=j&gt;) where {i in x, j in 2:NY - 1}</div>

Interior-only per-axis quantity diff2y of fire_behavior's classical Godunov level-set option (fire_upwinding = 4) on cartesian_uniform_2d over [1, NX] x [2, NY-1], where the three-point stencil along y fits: from the backward difference D- = (u[i,j] - u[i-1,j])/dy and the forward difference D+ = (u[i+1,j] - u[i,j])/dy (index shifted along y), each cell takes stencils/godunov_eno1_D1 -- the ENO1 minmod of D-, D+ within bdy_eno1_cells of any wall (fire_behavior's BDY_ENO1 band, which overrides every upwinding option), and the Godunov-upwind quantity max(D-, 0) - min(D+, 0) elsewhere. The Godunov quantity is the per-axis building block whose squares sum to the Godunov Hamiltonian only up to the cross term (max(D-,0) - min(D+,0))^2 = max(D-,0)^2 + min(D+,0)^2 - 2 max(D-,0) min(D+,0) >= the Rouy-Tourin Hamiltonian's per-axis term; fire_behavior uses it for the front NORMAL (nvx = diff2x / sqrt(grad^2 + eps)), while |grad psi| itself comes from godunov_norm_D2_*. First order; exact on a linear field with positive slope along y (both branches return the slope). See stencils/godunov_eno1_D1 for the sign caveat (outside the wall band, a magnitude for negative slopes). Free name bdy_eno1_cells (>= 1), resolved in the consuming model's scope like the grid geometry dx, dy.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `grad4_dot_wind`

<p class="rule-meta"><code>op:selector</code> <code>order:4</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/grad4_dot_wind.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/grad4_dot_wind.esm)

Parameters: `u`, `i`, `j`

<div class="math-block">central_D1x_4_point&lt;i=i, j=j, u=u&gt; * wind_u[i, j] + central_D1y_4_point&lt;i=i, j=j, u=u&gt; * wind_v[i, j]</div>

Pointwise advection term grad4 u . (wind_u, wind_v) at the cell (i, j) of cartesian_uniform_2d, with grad4 the fourth-order central gradient (stencils/central_D1x_4_point, stencils/central_D1y_4_point): fire_behavior's a_valor = Select_4th(dx, ...) uf(i, j) + Select_4th(dy, ...) vf(i, j) in Calc_tend_ls. Multiplied by the central derivative along an axis it gives that axis's WENO5 upwind selector signo = a_valor d4 u/d(axis), which is how the fire_upwinding = 9 and 10 schemes pick the upwind side from the wind. u is the field, i and j the cell indices, bound by the caller's faq. FREE NAMES (consumer-supplied, like the grid geometry dx, dy): wind_u and wind_v, the [x, y] fields of the advecting wind (fire_behavior's mid-flame wind uf, vf). It reaches two cells either side along both axes.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `varcoeff_laplacian_x_interior`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/varcoeff_laplacian_x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/varcoeff_laplacian_x_interior.esm)

Parameters: `u`, `k`

<div class="math-block">sum[i, j] ((0.5 * (k[i, j] + k[i + 1, j]) * (u[i + 1, j] - u[i, j]) - 0.5 * (k[i - 1, j] + k[i, j]) * (u[i, j] - u[i - 1, j])) / (dx * dx)) where {i in 2:NX - 1, j in y}</div>

Interior-only conservative finite-volume x-half of the variable-coefficient / nonlinear Laplacian div(k grad u) on cartesian_uniform_2d: the x-flux-divergence Lx u[i,j] = ( kf[i+1/2,j] (u[i+1,j]-u[i,j]) - kf[i-1/2,j] (u[i,j]-u[i-1,j]) ) / dx^2 over i in [2, NX-1], all j, the difference of the two x-face-normal fluxes divided by the uniform cell measure dx^2. This is the exact 2-D per-axis lift of grids/cartesian_uniform_1d/stencils/varcoeff_laplacian_interior with the y index a free (full-rank) spectator: the same arithmetic-mean face algebra runs independently on every y-row j. The x-face conductivity is the ARITHMETIC MEAN of the two adjacent cell values along x, kf[i+1/2,j] = (k[i,j]+k[i+1,j])/2; this is the standard second-order harmonic-vs-arithmetic choice for a smooth coefficient (both are O(dx^2)-accurate at the face for a differentiable k, and the arithmetic mean keeps the operator affine in k so it serves the nonlinear case k=k(u) with no reformulation). Both the state `u` and the coefficient `k` are consumer-supplied fields of shape [x, y] passed as template parameters (bound by the matching rule from the compound D(k*D(u,x),x) the consumer writes); the SAME rule serves the variable-coefficient case (k=k(x,y) a given field) and the nonlinear case (k declared as an observed k=f(u)) because the body only reads k by index at the x-faces. This interior aggregate loops over `i` and `j`; a coefficient k supplied as an aggregate observed must therefore loop over DIFFERENT indices (the library convention is `gi`/`gj`, as the grid geometry observeds use) so its inlining into this aggregate is not captured — see the rule (rules/varcoeff_laplacian_x_zero_flux_bc) authoring contract. Match-less named template; the boundary-condition layer (rules/varcoeff_laplacian_x_zero_flux_bc et al.) wraps it into the complete rewrite rule and adds the two x-edge faces. `dx` is a free name resolving to the consuming model's x-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y), each half lowered by its own per-axis rule; this stencil is the interior of the first half. Reduces to the constant-coefficient 3-point (k/dx^2)(u[i+1,j]-2u[i,j]+u[i-1,j]) when k is spatially constant.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient / conservative form) and Chapter 3 (multidimensional operators by dimension-wise assembly).

### `varcoeff_laplacian_y_interior`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/varcoeff_laplacian_y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/varcoeff_laplacian_y_interior.esm)

Parameters: `u`, `k`

<div class="math-block">sum[i, j] ((0.5 * (k[i, j] + k[i, j + 1]) * (u[i, j + 1] - u[i, j]) - 0.5 * (k[i, j - 1] + k[i, j]) * (u[i, j] - u[i, j - 1])) / (dy * dy)) where {i in x, j in 2:NY - 1}</div>

Interior-only conservative finite-volume y-half of the variable-coefficient / nonlinear Laplacian div(k grad u) on cartesian_uniform_2d: the y-flux-divergence Ly u[i,j] = ( kf[i,j+1/2] (u[i,j+1]-u[i,j]) - kf[i,j-1/2] (u[i,j]-u[i,j-1]) ) / dy^2 over j in [2, NY-1], all i, the difference of the two y-face-normal fluxes divided by the uniform cell measure dy^2. This is the exact 2-D per-axis lift of grids/cartesian_uniform_1d/stencils/varcoeff_laplacian_interior with the x index a free (full-rank) spectator: the same arithmetic-mean face algebra runs independently on every x-column i. The y-face conductivity is the ARITHMETIC MEAN of the two adjacent cell values along y, kf[i,j+1/2] = (k[i,j]+k[i,j+1])/2; this is the standard second-order harmonic-vs-arithmetic choice for a smooth coefficient (both are O(dy^2)-accurate at the face for a differentiable k, and the arithmetic mean keeps the operator affine in k so it serves the nonlinear case k=k(u) with no reformulation). Both the state `u` and the coefficient `k` are consumer-supplied fields of shape [x, y] passed as template parameters (bound by the matching rule from the compound D(k*D(u,y),y) the consumer writes); the SAME rule serves the variable-coefficient case (k=k(x,y) a given field) and the nonlinear case (k declared as an observed k=f(u)) because the body only reads k by index at the y-faces. This interior aggregate loops over `i` and `j`; a coefficient k supplied as an aggregate observed must therefore loop over DIFFERENT indices (the library convention is `gi`/`gj`, as the grid geometry observeds use) so its inlining into this aggregate is not captured — see the rule (rules/varcoeff_laplacian_y_zero_flux_bc) authoring contract. Match-less named template; the boundary-condition layer (rules/varcoeff_laplacian_y_zero_flux_bc et al.) wraps it into the complete rewrite rule and adds the two y-edge faces. `dy` is a free name resolving to the consuming model's y-spacing variable (the grid's consumer-supplied geometry contract; see grids/cartesian_uniform_2d/grid.esm). The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y), each half lowered by its own per-axis rule; this stencil is the interior of the second half. Reduces to the constant-coefficient 3-point (k/dy^2)(u[i,j+1]-2u[i,j]+u[i,j-1]) when k is spatially constant.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient / conservative form) and Chapter 3 (multidimensional operators by dimension-wise assembly).

### `weno3_eno1_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_eno1_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_eno1_D1x_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or abs(u[i, j]) &gt;= weno_band_cells * dx, eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;, weno3_js_upwind_D1x&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1x_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 4:NX - 3, j in 3:NY - 2}</div>

Interior-only hybrid WENO3/ENO1 upwind first derivative d u/dx of NCAR fire_behavior's fire_upwinding = 8 level-set scheme on cartesian_uniform_2d over the block [4, NX-3] x [3, NY-2]: stencils/weno5_eno1_D1x_interior (fire_upwinding = 9) with the WENO branch stencils/weno3_js_upwind_D1x in place of stencils/weno5_js_upwind_D1x, which is the only difference between Calc_tend_ls cases 8 and 9. Per cell the switch is the same: ENO1 minmod (stencils/eno1_minmod_D1x) within bdy_eno1_cells of any wall (stencils/eno1_wall_band) or where |u[i,j]| >= weno_band_cells * dx (outside the band around the front), otherwise the third-order WENO3 derivative with the upwind side picked by signo_x = (grad4 u . U) (d4 u/dx) (stencils/grad4_dot_wind times stencils/central_D1x_4_point). The block is the one the WENO5 sibling uses, so this stencil fits the shared frame rules/weno5_eno1_D1x_extrapolate_frame unchanged: with bdy_eno1_cells >= 3 (fire_behavior's BDY_ENO1 = 10 is a compile-time constant) every cell outside the block is ENO1 in the reference code too, so the shorter WENO3 reach needs no other tiling. FREE NAMES: wind_u, wind_v, bdy_eno1_cells (>= 3), weno_band_cells, and the grid's dx, dy. Both ifelse branches are evaluated on every cell (esm ifelse is not lazy). Order: third in smooth parts of the WENO band, first in the ENO1 regions; exact on a linear field.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno3_eno1_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_eno1_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_eno1_D1y_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or abs(u[i, j]) &gt;= weno_band_cells * dx, eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;, weno3_js_upwind_D1y&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1y_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 3:NX - 2, j in 4:NY - 3}</div>

Interior-only hybrid WENO3/ENO1 upwind first derivative d u/dy of NCAR fire_behavior's fire_upwinding = 8 level-set scheme on cartesian_uniform_2d over the block [3, NX-2] x [4, NY-3]: stencils/weno5_eno1_D1y_interior (fire_upwinding = 9) with the WENO branch stencils/weno3_js_upwind_D1y in place of stencils/weno5_js_upwind_D1y, which is the only difference between Calc_tend_ls cases 8 and 9. Per cell the switch is the same: ENO1 minmod (stencils/eno1_minmod_D1y) within bdy_eno1_cells of any wall (stencils/eno1_wall_band) or where |u[i,j]| >= weno_band_cells * dx (outside the band around the front), otherwise the third-order WENO3 derivative with the upwind side picked by signo_y = (grad4 u . U) (d4 u/dy) (stencils/grad4_dot_wind times stencils/central_D1y_4_point). The block is the one the WENO5 sibling uses, so this stencil fits the shared frame rules/weno5_eno1_D1y_extrapolate_frame unchanged: with bdy_eno1_cells >= 3 (fire_behavior's BDY_ENO1 = 10 is a compile-time constant) every cell outside the block is ENO1 in the reference code too, so the shorter WENO3 reach needs no other tiling. FREE NAMES: wind_u, wind_v, bdy_eno1_cells (>= 3), weno_band_cells, and the grid's dx, dy. Both ifelse branches are evaluated on every cell (esm ifelse is not lazy). Order: third in smooth parts of the WENO band, first in the ENO1 regions; exact on a linear field.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno3_js_face_recon`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:3</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_js_face_recon.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_js_face_recon.esm)

Parameters: `v0`, `v1`, `v2`

<div class="math-block">(weno_js_nonlinear_weight&lt;beta=(v1 - v0)^2, gamma=1 / 3&gt; * ((-1 * v0 + 3 * v1) / 2) + weno_js_nonlinear_weight&lt;beta=(v2 - v1)^2, gamma=2 / 3&gt; * ((v1 + v2) / 2)) / (weno_js_nonlinear_weight&lt;beta=(v1 - v0)^2, gamma=1 / 3&gt; + weno_js_nonlinear_weight&lt;beta=(v2 - v1)^2, gamma=2 / 3&gt;)</div>

Pointwise classical (Jiang-Shu 1996) WENO3 reconstruction of a point field at the face between the second and third cells of a three-cell window v0, v1, v2 (cells i-1, i, i+1 for the face i+1/2, upwind-biased to the left): the two second-order candidates p0 = (-v0 + 3 v1)/2 on {i-1, i} and p1 = (v1 + v2)/2 on {i, i+1}, linear weights gamma0 = 1/3 and gamma1 = 2/3 (whose linear combination is the third-order upwind-biased value (-v0 + 5 v1 + 2 v2)/6), smoothness indicators beta0 = (v1 - v0)^2 and beta1 = (v2 - v1)^2, and the classical nonlinear weights alpha_k = gamma_k/(beta_k + eps)^2 of stencils/weno_js_nonlinear_weight, normalised as (alpha0 p0 + alpha1 p1)/(alpha0 + alpha1) -- fire_behavior's flux_p / flux_m in Select_weno3. It is the third-order sibling of stencils/weno5_js_face_recon, with the same weight construction, eps and orientation, and NOT the same scheme as grids/cartesian_uniform_1d/stencils/weno3_recon (WENO-Z, Borges 2008: alpha_k = gamma_k (1 + (tau3/(eps + beta_k))^2), tau3 = |beta0 - beta1|): the two agree where the field is smooth and weight the candidates differently across a kink, so they must not be substituted for one another. Third order in smooth regions, degrading to the single-candidate value across a discontinuity. The mirrored (right-biased) face is obtained by binding the window in reverse order. Match-less helper of stencils/weno3_js_upwind_D1.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno3_js_upwind_D1`

<p class="rule-meta"><code>op:D</code> <code>order:3</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1.esm)

Parameters: `h`, `m2`, `m1`, `c`, `p1`, `p2`, `sgn`

<div class="math-block">ifelse(sgn &gt;= 0, (weno3_js_face_recon&lt;v0=m1, v1=c, v2=p1&gt; - weno3_js_face_recon&lt;v0=m2, v1=m1, v2=c&gt;) / h, (weno3_js_face_recon&lt;v0=p2, v1=p1, v2=c&gt; - weno3_js_face_recon&lt;v0=p1, v1=c, v2=m1&gt;) / h)</div>

Pointwise third-order upwind first derivative along one axis from the five cell values m2, m1 (two cells before), c (the cell), p1, p2 (two cells after) and the spacing h, as the difference of two classical Jiang-Shu WENO3 face reconstructions (stencils/weno3_js_face_recon): for a non-negative selector sgn the left-biased form (R(m1,c,p1) - R(m2,m1,c))/h (faces i+1/2 and i-1/2 reconstructed from the upwind side), for sgn < 0 the mirror image (R(p2,p1,c) - R(p1,c,m1))/h -- exactly fire_behavior's Select_weno3 with its uf argument as the selector: that code swaps the window (im1 <-> ip1, im2 <- ip2) and returns (flux_m - flux_p)/dx, which is the same expression. Each branch reads four of the five values (the non-negative branch never reads p2, the negative one never m2), so the reach is two cells upwind and one downwind. The selector is the caller's business (see stencils/weno3_js_upwind_D1x). The third-order sibling of stencils/weno5_js_upwind_D1; bound per axis by stencils/weno3_js_upwind_D1x and stencils/weno3_js_upwind_D1y.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno3_js_upwind_D1x`

<p class="rule-meta"><code>op:D</code> <code>order:3</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1x.esm)

Parameters: `u`, `i`, `j`, `sgn`

<div class="math-block">weno3_js_upwind_D1&lt;c=u[i, j], h=dx, m1=u[i - 1, j], m2=u[i - 2, j], p1=u[i + 1, j], p2=u[i + 2, j], sgn=sgn&gt;</div>

Pointwise third-order Jiang-Shu WENO3 upwind derivative d u/dx at the cell (i, j) of cartesian_uniform_2d: stencils/weno3_js_upwind_D1 bound to the five cells u at offsets -2..+2 along x and spacing dx, with the upwind side chosen by the sign of the supplied selector sgn (sgn >= 0 takes the left-biased reconstruction, sgn < 0 its mirror image, as in fire_behavior's Select_weno3). fire_behavior's Select_weno3 (dx, lfn(i, j), ..., signo_x). u is the field, i and j the cell indices, bound by the caller's faq; sgn is an expression the caller forms at (i, j) -- the wind-relative sign (stencils/grad4_dot_wind times stencils/central_D1x_4_point) in fire_upwinding = 8. It reaches two cells either side along x, so a caller evaluates it only on the block where that fits. The one spelling of the WENO3 branch along x; the third-order sibling of stencils/weno5_js_upwind_D1x.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno3_js_upwind_D1y`

<p class="rule-meta"><code>op:D</code> <code>order:3</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno3_js_upwind_D1y.esm)

Parameters: `u`, `i`, `j`, `sgn`

<div class="math-block">weno3_js_upwind_D1&lt;c=u[i, j], h=dy, m1=u[i, j - 1], m2=u[i, j - 2], p1=u[i, j + 1], p2=u[i, j + 2], sgn=sgn&gt;</div>

Pointwise third-order Jiang-Shu WENO3 upwind derivative d u/dy at the cell (i, j) of cartesian_uniform_2d: stencils/weno3_js_upwind_D1 bound to the five cells u at offsets -2..+2 along y and spacing dy, with the upwind side chosen by the sign of the supplied selector sgn (sgn >= 0 takes the left-biased reconstruction, sgn < 0 its mirror image, as in fire_behavior's Select_weno3). fire_behavior's Select_weno3 (dy, lfn(i, j), ..., signo_y). u is the field, i and j the cell indices, bound by the caller's faq; sgn is an expression the caller forms at (i, j) -- the wind-relative sign (stencils/grad4_dot_wind times stencils/central_D1y_4_point) in fire_upwinding = 8. It reaches two cells either side along y, so a caller evaluates it only on the block where that fits. The one spelling of the WENO3 branch along y; the third-order sibling of stencils/weno5_js_upwind_D1y.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_D1x_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or abs(u[i, j]) &gt;= weno_band_cells * dx, eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;, weno5_js_upwind_D1x&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1x_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 4:NX - 3, j in 3:NY - 2}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dx of NCAR fire_behavior's default level-set scheme (fire_upwinding = 9; Munoz-Esparza et al. 2018) on cartesian_uniform_2d over the block [4, NX-3] x [3, NY-2], where the seven-point WENO stencil along x and the five-point fourth-order selector stencil along BOTH axes fit without reaching a wall. Per cell the SWITCH is written on the field itself and on the cell's distance from the walls, as the reference code does: the cell takes the first-order ENO1 minmod of its two one-sided differences (stencils/eno1_minmod_D1x) when it lies within bdy_eno1_cells of any wall (stencils/eno1_wall_band: i < 1 + bdy_eno1_cells, i > NX - bdy_eno1_cells, likewise j) OR when |u[i,j]| >= weno_band_cells * dx (outside the band around the zero level set, the front); otherwise -- inside the band, away from the walls -- it takes the fifth-order Jiang-Shu WENO5 upwind derivative (stencils/weno5_js_upwind_D1x) whose upwind side is selected by the sign of signo_x = (grad4 u . U) (d4 u/dx), the product of the fourth-order central advection term grad4 u . (wind_u, wind_v) with the fourth-order central derivative along this axis (stencils/grad4_dot_wind times stencils/central_D1x_4_point) -- fire_behavior's a_valor * Select_4th. FREE NAMES (consumer-supplied, the same contract as the grid geometry dx, dy): wind_u and wind_v, the [x, y] fields of the advecting wind that biases the WENO stencil (fire_behavior's mid-flame wind uf, vf; a level set with no wind may bind zero fields, which selects the left-biased stencil everywhere); bdy_eno1_cells, the ENO1 wall band in cells (fire_behavior's BDY_ENO1 = 10, a compile-time constant there; it MUST be >= 3 so that the WENO7-point reach never crosses a wall, which this block guarantees structurally only for the outermost three cells -- the rule regions -- while the band itself is a value comparison on the index); weno_band_cells, the half-width of the WENO band in cells (fire_behavior's fire_lsm_band_ngp, namelist default 4). Both branches of the ifelse are evaluated on every cell (esm ifelse is not lazy), which is why this template is restricted to the block where the WENO stencil is in range; the boundary strips are ENO1-only regions of the rule's frame (rules/weno5_eno1_D1x_extrapolate_frame). Order: fifth in smooth parts of the WENO band away from critical points (Jiang-Shu weights), first in the ENO1 regions -- the scheme is formally first order; the rule is tagged order:1. Exact on a linear field (every candidate reconstruction and the minmod reproduce the slope).

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_D1y_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or abs(u[i, j]) &gt;= weno_band_cells * dx, eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;, weno5_js_upwind_D1y&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1y_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 3:NX - 2, j in 4:NY - 3}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dy of NCAR fire_behavior's default level-set scheme (fire_upwinding = 9; Munoz-Esparza et al. 2018) on cartesian_uniform_2d over the block [3, NX-2] x [4, NY-3], where the seven-point WENO stencil along y and the five-point fourth-order selector stencil along BOTH axes fit without reaching a wall. Per cell the SWITCH is written on the field itself and on the cell's distance from the walls, as the reference code does: the cell takes the first-order ENO1 minmod of its two one-sided differences (stencils/eno1_minmod_D1y) when it lies within bdy_eno1_cells of any wall (stencils/eno1_wall_band: i < 1 + bdy_eno1_cells, i > NX - bdy_eno1_cells, likewise j) OR when |u[i,j]| >= weno_band_cells * dx (outside the band around the zero level set, the front); otherwise -- inside the band, away from the walls -- it takes the fifth-order Jiang-Shu WENO5 upwind derivative (stencils/weno5_js_upwind_D1y) whose upwind side is selected by the sign of signo_y = (grad4 u . U) (d4 u/dy), the product of the fourth-order central advection term grad4 u . (wind_u, wind_v) with the fourth-order central derivative along this axis (stencils/grad4_dot_wind times stencils/central_D1y_4_point) -- fire_behavior's a_valor * Select_4th. FREE NAMES (consumer-supplied, the same contract as the grid geometry dx, dy): wind_u and wind_v, the [x, y] fields of the advecting wind that biases the WENO stencil (fire_behavior's mid-flame wind uf, vf; a level set with no wind may bind zero fields, which selects the left-biased stencil everywhere); bdy_eno1_cells, the ENO1 wall band in cells (fire_behavior's BDY_ENO1 = 10, a compile-time constant there; it MUST be >= 3 so that the WENO7-point reach never crosses a wall, which this block guarantees structurally only for the outermost three cells -- the rule regions -- while the band itself is a value comparison on the index); weno_band_cells, the half-width of the WENO band in cells (fire_behavior's fire_lsm_band_ngp, namelist default 4). Both branches of the ifelse are evaluated on every cell (esm ifelse is not lazy), which is why this template is restricted to the block where the WENO stencil is in range; the boundary strips are ENO1-only regions of the rule's frame (rules/weno5_eno1_D1y_extrapolate_frame). Order: fifth in smooth parts of the WENO band away from critical points (Jiang-Shu weights), first in the ENO1 regions -- the scheme is formally first order; the rule is tagged order:1. Exact on a linear field (every candidate reconstruction and the minmod reproduce the slope).

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_blend_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_blend_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_blend_D1x_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (weno5_eno1_masked_blend&lt;e=eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;, m=ifelse(eno1_wall_band&lt;i=i, j=j&gt;, 0, weno5_eno1_cosine_mask&lt;a=abs(u[i, j]), b=weno_band_cells * dx, p=pi_blend&gt;), w=weno5_js_upwind_D1x&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1x_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;&gt;) where {i in 4:NX - 3, j in 3:NY - 2}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dx of NCAR fire_behavior's fire_upwinding = 10 -- the same scheme as fire_upwinding = 9 (stencils/weno5_eno1_D1x_interior) but with a BLENDING ZONE instead of a hard band edge -- on cartesian_uniform_2d over the block [4, NX-3] x [3, NY-2]. Per cell it is stencils/weno5_eno1_masked_blend of the WENO5 derivative (stencils/weno5_js_upwind_D1x, upwind side by the wind-relative sign stencils/grad4_dot_wind times stencils/central_D1x_4_point, exactly as in the fire_upwinding = 9 sibling) and the ENO1 derivative (stencils/eno1_minmod_D1x), with mask m = stencils/weno5_eno1_cosine_mask(|u[i,j]|, weno_band_cells * dx, pi_blend) -- 1 in the band, a cosine taper over B < |u| <= 2B, 0 beyond -- forced to 0 on the ENO1 wall band (stencils/eno1_wall_band), because case 10 sits inside the same wall `else` as cases 8 and 9. Where m = 0 (the wall band and the far field) the blend SELECTS the ENO1 value, so the result there is bitwise the ENO1 derivative, as in the reference, which skips WENO5 when mask = 0; the WENO5 value is still evaluated there (esm ifelse is not lazy) but cannot leak in. ENO1 appears once: the wall band enters through the mask, not through a second ifelse. Free names: wind_u, wind_v, bdy_eno1_cells, weno_band_cells, pi_blend, dx, dy. pi_blend is pi as the consumer's constants spell it (fire_behavior's PI = 3.1415926); it is a free name rather than a literal so the consumer supplies it.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls case (10) (band_width = 2 threshold_hlu, transition_width = threshold_hlu, the three-branch cosine mask and the `if (mask > 0.0)` blend), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_blend_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_blend_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_blend_D1y_interior.esm)

Parameters: `u`

<div class="math-block">sum[i, j] (weno5_eno1_masked_blend&lt;e=eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;, m=ifelse(eno1_wall_band&lt;i=i, j=j&gt;, 0, weno5_eno1_cosine_mask&lt;a=abs(u[i, j]), b=weno_band_cells * dx, p=pi_blend&gt;), w=weno5_js_upwind_D1y&lt;i=i, j=j, sgn=grad4_dot_wind&lt;i=i, j=j, u=u&gt; * central_D1y_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;&gt;) where {i in 3:NX - 2, j in 4:NY - 3}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dy of NCAR fire_behavior's fire_upwinding = 10 -- the same scheme as fire_upwinding = 9 (stencils/weno5_eno1_D1y_interior) but with a BLENDING ZONE instead of a hard band edge -- on cartesian_uniform_2d over the block [3, NX-2] x [4, NY-3]. Per cell it is stencils/weno5_eno1_masked_blend of the WENO5 derivative (stencils/weno5_js_upwind_D1y, upwind side by the wind-relative sign stencils/grad4_dot_wind times stencils/central_D1y_4_point, exactly as in the fire_upwinding = 9 sibling) and the ENO1 derivative (stencils/eno1_minmod_D1y), with mask m = stencils/weno5_eno1_cosine_mask(|u[i,j]|, weno_band_cells * dx, pi_blend) -- 1 in the band, a cosine taper over B < |u| <= 2B, 0 beyond -- forced to 0 on the ENO1 wall band (stencils/eno1_wall_band), because case 10 sits inside the same wall `else` as cases 8 and 9. Where m = 0 (the wall band and the far field) the blend SELECTS the ENO1 value, so the result there is bitwise the ENO1 derivative, as in the reference, which skips WENO5 when mask = 0; the WENO5 value is still evaluated there (esm ifelse is not lazy) but cannot leak in. ENO1 appears once: the wall band enters through the mask, not through a second ifelse. Free names: wind_u, wind_v, bdy_eno1_cells, weno_band_cells, pi_blend, dx, dy. pi_blend is pi as the consumer's constants spell it (fire_behavior's PI = 3.1415926); it is a free name rather than a literal so the consumer supplies it.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls case (10) (band_width = 2 threshold_hlu, transition_width = threshold_hlu, the three-branch cosine mask and the `if (mask > 0.0)` blend), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_cosine_mask`

<p class="rule-meta"><code>op:selector</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_cosine_mask.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_cosine_mask.esm)

Parameters: `a`, `b`, `p`

<div class="math-block">ifelse(a &lt;= b, 1, ifelse(a &lt;= 2 * b, 0.5 * (1 + cos(p * (a - b) / b)), 0))</div>

Pointwise cosine blending mask of NCAR fire_behavior's fire_upwinding = 10 (physics/level_set_mod.F90 Calc_tend_ls case 10), the weight of the WENO5 derivative in the blend of stencils/weno5_eno1_masked_blend. With a = |u[i,j]| the cell's distance from the front and b the band half-width (threshold_hlu = weno_band_cells dx in the reference, where band_width = 2 b and transition_width = b), the mask is the reference's three branches written as they are: 1 for a <= b (inside the band: pure WENO5), 0.5 (1 + cos(p (a - b) / b)) for b < a <= 2 b (the transition: a cosine taper from 1 to 0), and 0 for a > 2 b (outside: pure ENO1). The two ends are therefore EXACTLY 1 and 0 whatever p is -- in particular with fire_behavior's truncated PI = 3.1415926, for which the taper itself only reaches 0.5 (1 + cos(3.1415926)) = 7e-16 at a = 2 b, as it does in the reference. The taper spells the reference's own operation order, PI * (a - (band_width - transition_width)) / transition_width, with band_width - transition_width = 2 b - b = b exactly in floating point. p is pi, bound by the caller (the blend stencils bind the free name pi_blend) so the consumer supplies the value its own constants use. Params: a (>= 0), b (> 0), p.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls case (10) (band_width = 2 threshold_hlu, transition_width = threshold_hlu, the three-branch cosine mask and the `if (mask > 0.0)` blend), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_masked_blend`

<p class="rule-meta"><code>op:blend</code> <code>order:1</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_masked_blend.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_masked_blend.esm)

Parameters: `m`, `w`, `e`

<div class="math-block">ifelse(m &gt; 0, m * w + (1 - m) * e, e)</div>

Pointwise masked convex combination of a high-order and a low-order derivative, the blending step of NCAR fire_behavior's fire_upwinding = 10 (physics/level_set_mod.F90 Calc_tend_ls case 10): m w + (1 - m) e where the mask m > 0, and e itself where m = 0 -- the reference's `if (mask > 0.0) then diff2 = mask diff2_weno + (1 - mask) diff2_eno else diff2 = diff2_eno`. The else branch is SELECTED, not multiplied by zero, so where m = 0 the result is exactly e and w (the WENO5 value, which is still evaluated because esm ifelse is not lazy) cannot leak into it, even if it is non-finite. m is the caller's mask (stencils/weno5_eno1_cosine_mask, forced to 0 in the ENO1 wall band by the blend stencils), w the WENO5 derivative, e the ENO1 derivative. No grid dependence.

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls case (10) (band_width = 2 threshold_hlu, transition_width = threshold_hlu, the three-branch cosine mask and the `if (mask > 0.0)` blend), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_signed_D1x_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_signed_D1x_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_signed_D1x_interior.esm)

Parameters: `u`, `s`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or u[i, j] &gt;= weno_band_cells * dx, eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;, weno5_js_upwind_D1x&lt;i=i, j=j, sgn=s[i, j] * central_D1x_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 4:NX - 3, j in 3:NY - 2}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dx for the level-set REINITIALISATION equation psi_tau = S(psi_0)(1 - |grad psi|) of NCAR fire_behavior (physics/level_set_mod.F90 Advance_ls_reinit, fire_upwinding_reinit = 4), on cartesian_uniform_2d over the block [4, NX-3] x [3, NY-2]. It is stencils/weno5_eno1_D1x_interior with the two changes the reinitialisation makes, and nothing else. (1) THE UPWIND SIDE IS CHOSEN BY A SUPPLIED SIGN FIELD, not by the wind: the propagation scheme selects WENO5's side by signo = a_valor (d4 u/dx) with a_valor = grad4 u . (wind_u, wind_v) (stencils/grad4_dot_wind); the reinitialisation by s[i,j] (d4 u/dx) (stencils/central_D1x_4_point), where s is the param bound to fire_behavior's lfn_s0 -- the smoothed sign S(psi_0) = psi_0 / sqrt(psi_0^2 + dx^2), frozen at the start of the pseudo-time iteration. This stencil reads no wind. (2) THE BAND TEST IS SIGNED, NOT ABSOLUTE: ENO1 is taken when u[i,j] >= weno_band_cells * dx, whereas the propagation stencil takes it when |u[i,j]| >= weno_band_cells * dx. That is faithful to the reference code, which writes `if (lfn_curr(i, j) < threshold_hlu)` in every WENO case of Advance_ls_reinit while Calc_tend_ls writes `if (abs (lfn(i, j)) < threshold_hlu)`. The missing abs is almost certainly an oversight -- it puts the whole BURNED INTERIOR, where the level set is large and negative, on the fifth-order branch meant for the band around the front -- but it is what the model computes, and a transcription that silently inserted the abs would not reproduce it. Everything else is the propagation sibling's own pieces: the wall band (stencils/eno1_wall_band), the ENO1 branch (stencils/eno1_minmod_D1x), the WENO5 branch (stencils/weno5_js_upwind_D1x), the block bounds. Params: u, the level set, and s, the [x, y] sign field (indexed at the cell). Free names: bdy_eno1_cells, weno_band_cells, dx, dy.

**References**

- Sussman, M., Smereka, P., & Osher, S. (1994). A level set approach for computing solutions to incompressible two-phase flow. Journal of Computational Physics, 114(1), 146-159. The reinitialisation equation psi_tau = S(psi_0)(1 - \|grad psi\|) with the smoothed sign S(psi_0) = psi_0/sqrt(psi_0^2 + dx^2). [`doi:10.1006/jcph.1994.1155`](https://doi.org/10.1006/jcph.1994.1155)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Advance_ls_reinit (fire_upwinding_reinit = 4), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_eno1_signed_D1y_interior`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_eno1_signed_D1y_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_eno1_signed_D1y_interior.esm)

Parameters: `u`, `s`

<div class="math-block">sum[i, j] (ifelse(eno1_wall_band&lt;i=i, j=j&gt; or u[i, j] &gt;= weno_band_cells * dx, eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;, weno5_js_upwind_D1y&lt;i=i, j=j, sgn=s[i, j] * central_D1y_4_point&lt;i=i, j=j, u=u&gt;, u=u&gt;)) where {i in 3:NX - 2, j in 4:NY - 3}</div>

Interior-only hybrid WENO5/ENO1 upwind first derivative d u/dy for the level-set REINITIALISATION equation psi_tau = S(psi_0)(1 - |grad psi|) of NCAR fire_behavior (physics/level_set_mod.F90 Advance_ls_reinit, fire_upwinding_reinit = 4), on cartesian_uniform_2d over the block [3, NX-2] x [4, NY-3]. It is stencils/weno5_eno1_D1y_interior with the two changes the reinitialisation makes, and nothing else. (1) THE UPWIND SIDE IS CHOSEN BY A SUPPLIED SIGN FIELD, not by the wind: the propagation scheme selects WENO5's side by signo = a_valor (d4 u/dy) with a_valor = grad4 u . (wind_u, wind_v) (stencils/grad4_dot_wind); the reinitialisation by s[i,j] (d4 u/dy) (stencils/central_D1y_4_point), where s is the param bound to fire_behavior's lfn_s0 -- the smoothed sign S(psi_0) = psi_0 / sqrt(psi_0^2 + dx^2), frozen at the start of the pseudo-time iteration. This stencil reads no wind. (2) THE BAND TEST IS SIGNED, NOT ABSOLUTE: ENO1 is taken when u[i,j] >= weno_band_cells * dx, whereas the propagation stencil takes it when |u[i,j]| >= weno_band_cells * dx. That is faithful to the reference code, which writes `if (lfn_curr(i, j) < threshold_hlu)` in every WENO case of Advance_ls_reinit while Calc_tend_ls writes `if (abs (lfn(i, j)) < threshold_hlu)`. The missing abs is almost certainly an oversight -- it puts the whole BURNED INTERIOR, where the level set is large and negative, on the fifth-order branch meant for the band around the front -- but it is what the model computes, and a transcription that silently inserted the abs would not reproduce it. Everything else is the propagation sibling's own pieces: the wall band (stencils/eno1_wall_band), the ENO1 branch (stencils/eno1_minmod_D1y), the WENO5 branch (stencils/weno5_js_upwind_D1y), the block bounds. Params: u, the level set, and s, the [x, y] sign field (indexed at the cell). Free names: bdy_eno1_cells, weno_band_cells, dx, dy.

**References**

- Sussman, M., Smereka, P., & Osher, S. (1994). A level set approach for computing solutions to incompressible two-phase flow. Journal of Computational Physics, 114(1), 146-159. The reinitialisation equation psi_tau = S(psi_0)(1 - \|grad psi\|) with the smoothed sign S(psi_0) = psi_0/sqrt(psi_0^2 + dx^2). [`doi:10.1006/jcph.1994.1155`](https://doi.org/10.1006/jcph.1994.1155)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Advance_ls_reinit (fire_upwinding_reinit = 4), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_js_face_recon`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:5</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_js_face_recon.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_js_face_recon.esm)

Parameters: `v0`, `v1`, `v2`, `v3`, `v4`

<div class="math-block">(1 / 10 / (13 / 12 * (v0 + -2 * v1 + v2)^2 + 1 / 4 * (v0 + -4 * v1 + 3 * v2)^2 + 1.0e-6)^2 * ((2 * v0 + -7 * v1 + 11 * v2) / 6) + 3 / 5 / (13 / 12 * (v1 + -2 * v2 + v3)^2 + 1 / 4 * (v1 - v3)^2 + 1.0e-6)^2 * ((-1 * v1 + 5 * v2 + 2 * v3) / 6) + 3 / 10 / (13 / 12 * (v2 + -2 * v3 + v4)^2 + 1 / 4 * (3 * v2 + -4 * v3 + v4)^2 + 1.0e-6)^2 * ((2 * v2 + 5 * v3 + -1 * v4) / 6)) / (1 / 10 / (13 / 12 * (v0 + -2 * v1 + v2)^2 + 1 / 4 * (v0 + -4 * v1 + 3 * v2)^2 + 1.0e-6)^2 + 3 / 5 / (13 / 12 * (v1 + -2 * v2 + v3)^2 + 1 / 4 * (v1 - v3)^2 + 1.0e-6)^2 + 3 / 10 / (13 / 12 * (v2 + -2 * v3 + v4)^2 + 1 / 4 * (3 * v2 + -4 * v3 + v4)^2 + 1.0e-6)^2)</div>

Pointwise classical (Jiang-Shu 1996) WENO5 reconstruction of a cell-averaged/point field at the face between the third and fourth cells of a five-cell window v0..v4 (cells i-2, i-1, i, i+1, i+2 for the face i+1/2, upwind-biased to the left): the three third-order candidate stencils f1 = (2 v0 - 7 v1 + 11 v2)/6, f2 = (-v1 + 5 v2 + 2 v3)/6, f3 = (2 v2 + 5 v3 - v4)/6 with linear weights 1/10, 3/5, 3/10, smoothness indicators beta_k = 13/12 (second difference)^2 + 1/4 (first-difference combination)^2, and the classical nonlinear weights w_k = gamma_k/(beta_k + eps)^2 with eps = 1e-6 and power 2 -- NOT the WENO-Z weights of stencils/weno5_recon.esm on cartesian_uniform_1d, and NOT the Hamilton-Jacobi divided-difference form of hjweno_recon_*. eps = 1e-6 is the literal fire_behavior uses (TOL in Select_weno5), on the field's own value scale; it is not dimensionally scaled. Fifth order in smooth regions away from critical points, degrading to third order where a first derivative vanishes (the known Jiang-Shu limitation that WENO-Z repairs). Match-less helper of the weno5_eno1_* stencils of this grid; the mirrored (right-biased) face is obtained by binding the window in reverse order.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_js_upwind_D1`

<p class="rule-meta"><code>op:D</code> <code>order:5</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1.esm)

Parameters: `h`, `m3`, `m2`, `m1`, `c`, `p1`, `p2`, `p3`, `sgn`

<div class="math-block">ifelse(sgn &gt;= 0, (weno5_js_face_recon&lt;v0=m2, v1=m1, v2=c, v3=p1, v4=p2&gt; - weno5_js_face_recon&lt;v0=m3, v1=m2, v2=m1, v3=c, v4=p1&gt;) / h, (weno5_js_face_recon&lt;v0=p3, v1=p2, v2=p1, v3=c, v4=m1&gt; - weno5_js_face_recon&lt;v0=p2, v1=p1, v2=c, v3=m1, v4=m2&gt;) / h)</div>

Pointwise fifth-order upwind first derivative along one axis from the seven cell values m3, m2, m1 (three cells before), c (the cell), p1, p2, p3 (three cells after) and the spacing h, as the difference of two Jiang-Shu WENO5 face reconstructions: for a non-negative selector sgn the left-biased form (R(m2,m1,c,p1,p2) - R(m3,m2,m1,c,p1))/h (faces i+1/2 and i-1/2 reconstructed from the upwind side), for sgn < 0 the mirror image (R(p3,p2,p1,c,m1) - R(p2,p1,c,m1,m2))/h -- exactly fire_behavior's Select_weno5 with its uf argument as the selector: the code swaps the window (im1<->ip1, im2<->ip2, im3<-ip3) and returns (flux_m - flux_p)/dx, which is the same expression. The selector is the caller's business: fire_behavior's level set passes signo = (grad4 psi . U) (d4 psi/dx), the product of the fourth-order central advection term with the fourth-order central derivative along this axis, so the upwind side is chosen by the wind-relative propagation direction rather than by the sign of the derivative alone. Match-less helper of weno5_eno1_D1x_interior / weno5_eno1_D1y_interior.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_js_upwind_D1x`

<p class="rule-meta"><code>op:D</code> <code>order:5</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1x.esm)

Parameters: `u`, `i`, `j`, `sgn`

<div class="math-block">weno5_js_upwind_D1&lt;c=u[i, j], h=dx, m1=u[i - 1, j], m2=u[i - 2, j], m3=u[i - 3, j], p1=u[i + 1, j], p2=u[i + 2, j], p3=u[i + 3, j], sgn=sgn&gt;</div>

Pointwise fifth-order Jiang-Shu WENO5 upwind derivative d u/dx at the cell (i, j) of cartesian_uniform_2d: stencils/weno5_js_upwind_D1 bound to the seven cells u at offsets -3..+3 along x and spacing dx, with the upwind side chosen by the sign of the supplied selector sgn (sgn >= 0 takes the left-biased reconstruction, sgn < 0 its mirror image, as in fire_behavior's Select_weno5). fire_behavior's Select_weno5 (dx, lfn(i, j), ..., signo_x). u is the field, i and j the cell indices, bound by the caller's faq; sgn is an expression the caller forms at (i, j) -- the wind-relative sign (stencils/grad4_dot_wind times stencils/central_D1x_4_point) in the propagation schemes, the frozen smoothed sign times the same central derivative in the reinitialisation scheme. It reaches three cells either side along x, so a caller evaluates it only on the block where that fits. The one spelling of the WENO5 branch shared by the weno5_eno1_*D1x_interior stencils.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno5_js_upwind_D1y`

<p class="rule-meta"><code>op:D</code> <code>order:5</code> <code>stencil width 7</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno5_js_upwind_D1y.esm)

Parameters: `u`, `i`, `j`, `sgn`

<div class="math-block">weno5_js_upwind_D1&lt;c=u[i, j], h=dy, m1=u[i, j - 1], m2=u[i, j - 2], m3=u[i, j - 3], p1=u[i, j + 1], p2=u[i, j + 2], p3=u[i, j + 3], sgn=sgn&gt;</div>

Pointwise fifth-order Jiang-Shu WENO5 upwind derivative d u/dy at the cell (i, j) of cartesian_uniform_2d: stencils/weno5_js_upwind_D1 bound to the seven cells u at offsets -3..+3 along y and spacing dy, with the upwind side chosen by the sign of the supplied selector sgn (sgn >= 0 takes the left-biased reconstruction, sgn < 0 its mirror image, as in fire_behavior's Select_weno5). fire_behavior's Select_weno5 (dy, lfn(i, j), ..., signo_y). u is the field, i and j the cell indices, bound by the caller's faq; sgn is an expression the caller forms at (i, j) -- the wind-relative sign (stencils/grad4_dot_wind times stencils/central_D1y_4_point) in the propagation schemes, the frozen smoothed sign times the same central derivative in the reinitialisation scheme. It reaches three cells either side along y, so a caller evaluates it only on the block where that fits. The one spelling of the WENO5 branch shared by the weno5_eno1_*D1y_interior stencils.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

### `weno_js_nonlinear_weight`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_uniform_2d/stencils/weno_js_nonlinear_weight.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/stencils/weno_js_nonlinear_weight.esm)

Parameters: `gamma`, `beta`

<div class="math-block">gamma / (beta + 1.0e-6)^2</div>

Pointwise classical (Jiang-Shu 1996) unnormalised WENO weight alpha = gamma/(beta + eps)^2 of one candidate stencil, from its linear weight gamma and its smoothness indicator beta, with eps = 1e-6 and power 2 -- the w1t/w2t of fire_behavior's Select_weno3 (TOL = 1e-6, on the field's own value scale; not dimensionally scaled). A face reconstruction normalises these over its candidates, sum_k alpha_k p_k / sum_k alpha_k, so the one spelling here is what both the numerator and the denominator of stencils/weno3_js_face_recon read, and eps or the power cannot drift between them. NOT the WENO-Z weight of grids/cartesian_uniform_1d/stencils/weno3_recon. Match-less helper.

**References**

- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `central_D1x_dirichlet_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1x_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1x_dirichlet_bc.esm)

**Rewrites** `D(f)/Dx` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D1x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((f[2, j] + f[1, j] - 2 * gL) / (2 * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((2 * gR - f[NX, j] - f[NX - 1, j]) / (2 * dx)) where {i in NX:NX, j in y})</div>

Second-order centered first derivative along x with inhomogeneous DIRICHLET (fixed-value) conditions on the two x-walls, over the full 2-D field f[x,y] on cartesian_uniform_2d, cell-centered — the first-derivative companion of central_D2x_dirichlet_bc (same mirror ghost, different stencil). Interior i in [2, NX-1] (all j): the imported 3-point centered x-stencil (f[i+1,j]-f[i-1,j])/(2 dx). Faces i=1 and i=NX: mirror-ghost elimination against a prescribed wall value. `gL` and `gR` are consumer-supplied free names (the same free-name contract as the grid geometry `dx`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the left x-wall x=x0 and the right x-wall x=x0+NX*dx, constant along that wall (all j). On this cell-centered grid the wall sits half a cell outside the first center, so linear reconstruction to the wall gives the ghost f[0,j]=2*gL-f[1,j] (and f[NX+1,j]=2*gR-f[NX,j]); substituting into the centered stencil (f[2,j]-f[0,j])/(2 dx) and (f[NX+1,j]-f[NX-1,j])/(2 dx) yields D1x[1,j] = (f[2,j]+f[1,j]-2*gL)/(2 dx) and D1x[NX,j] = (2*gR-f[NX,j]-f[NX-1,j])/(2 dx). A consuming model declares gL/gR as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec §6: default, or 0), so gL=gR=0 is the homogeneous-Dirichlet x-wall and both wall terms drop out. The gL/gR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). `dx` is the grid's consumer-supplied x-spacing free name. The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_dirichlet_bc and central_D1x_zero_grad_bc, so every binding's simulator fills the (1 x NY) regions without rank promotion. FACE TRUNCATION: D1x[1,j] - u_x(x_1,y_j) = (dx/8) u_xx(x0,y_j) + O(dx^2), so the mirror-ghost face is O(dx) locally for a generic field (SUPERCONVERGENT to O(dx^2) when u_xx=0 at the wall). WELL-POSED USE: a fixed value pins the field at the wall, the DIFFUSION operator's natural boundary condition, so this gradient closure belongs in a well-posed advection-DIFFUSION setting (or flow into the wall), NOT a pure first-order advection with a fixed value at an OUTFLOW wall, where the pinned value reflects (a Kreiss/GKS-unstable mode); the MMS driver problems/advection_2d_x_dirichlet_central.esm exercises it as the advective x-term of an advection-diffusion decay problem whose physical x-diffusion supplies the boundary dissipation, confirming 2nd order. Matches plain D(f, wrt: x) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the first x-derivative of a bare field declared over this 2-D grid's [x, y] index sets (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec §9.7.7) the wrt literal and the where shape follow the renamed axes together. A model imports exactly one first-derivative rule for x (this Dirichlet rule, the zero-gradient rule, or the periodic variant), never more than one.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (cell-centered Dirichlet via mirror ghost) and Chapter 10 (well-posed boundary conditions for advection-diffusion).
- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

**Conformance fixtures**

- ast: [`central_D1x_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1x_dirichlet_bc) — status `active`
- convergence: [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_dirichlet_central) — status `active`
- simulation: [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_x_dirichlet_central) — status `active`

#### Convergence — case [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_dirichlet_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_x_dirichlet_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_x_dirichlet_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.690e-04 | — | 6.559e-03 | — |
| 32 | 1.442e-04 | 1.98 | 1.829e-03 | 1.84 |
| 64 | 3.617e-05 | 2.00 | 4.809e-04 | 1.93 |
| 128 | 9.050e-06 | 2.00 | 1.232e-04 | 1.96 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_x_dirichlet_central (log-log)](/plots/convergence/advection_2d_x_dirichlet_central.png)

### `central_D1x_periodic`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1x_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1x_periodic.esm)

**Rewrites** `D(f)/Dx` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D1x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((f[2, j] - f[NX, j]) / (2 * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((f[1, j] - f[NX - 1, j]) / (2 * dx)) where {i in NX:NX, j in y})</div>

Second-order centered first derivative along x with periodic (wrap) boundary on cartesian_uniform_2d, over the full 2-D field f[x,y]. Interior i in [2, NX-1] (all j): the imported x-stencil (f[i+1,j]-f[i-1,j])/(2 dx). The two x-face columns wrap periodically: i=1 uses the ghost f[0,j]=f[NX,j], giving (f[2,j]-f[NX,j])/(2 dx); i=NX uses f[NX+1,j]=f[1,j], giving (f[1,j]-f[NX-1,j])/(2 dx). The scheme is the same O(dx^2) centered difference at every column and conserves sum_x(f) exactly (the wrap telescopes). The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) rather than the esm-spec 9.6.8 worked example's rank-reduced form, mirroring grids/latlon/rules/central_D_lon_periodic so every binding's simulator fills the (1 x NY) regions without rank promotion. `dx` is the grid's consumer-supplied x-spacing free name; the periodic wrap is dx-scale-invariant, so the rule serves any x-extent unchanged. Matches plain D(f, wrt: x) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the first x-derivative of a bare field declared over this grid's [x, y] index sets (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. O(dx^2). A model imports exactly one first-derivative rule for x (this centered periodic rule or the zero-gradient variant), never both.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D1x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1x_periodic) — status `active`
- convergence: [`advection_2d_x_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_periodic_central) — status `active`
- simulation: [`advection_2d_x_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_x_periodic_central) — status `active`

#### Convergence — case [`advection_2d_x_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_periodic_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_x_periodic_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_x_periodic_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.578e-02 | — | 1.568e-02 | — |
| 32 | 3.964e-03 | 1.99 | 3.964e-03 | 1.98 |
| 64 | 9.924e-04 | 2.00 | 9.919e-04 | 2.00 |
| 128 | 2.482e-04 | 2.00 | 2.482e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_x_periodic_central (log-log)](/plots/convergence/advection_2d_x_periodic_central.png)

### `central_D1x_zero_grad_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1x_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1x_zero_grad_bc.esm)

**Rewrites** `D(f)/Dx` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D1x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((-3 * f[1, j] + 4 * f[2, j] - f[3, j]) / (2 * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((f[NX - 2, j] - 4 * f[NX - 1, j] + 3 * f[NX, j]) / (2 * dx)) where {i in NX:NX, j in y})</div>

Second-order centered first derivative along x with one-sided (zero-gradient frame) boundary columns on cartesian_uniform_2d, over the full 2-D field f[x,y]. Interior i in [2, NX-1] (all j): the imported x-stencil (f[i+1,j]-f[i-1,j])/(2 dx). The two x-face columns use 2nd-order one-sided differences that reach only into the interior (no ghost, no periodic wrap): i=1 gives (-3 f[1,j]+4 f[2,j]-f[3,j])/(2 dx) and i=NX gives (f[NX-2,j]-4 f[NX-1,j]+3 f[NX,j])/(2 dx), so the discretized derivative is fully defined with its BC and stays O(dx^2) at every column (esm-spec 9.6.8). The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring grids/cartesian_uniform_2d/rules/central_D2y_neumann_bc so every binding's simulator fills the (1 x NY) regions without rank promotion. `dx` is the grid's consumer-supplied x-spacing free name. The one-sided faces are O(dx^2) locally and preserve global 2nd-order convergence for boundary-compatible solutions — those whose odd x-derivatives vanish at the walls, so the discrete derivative matches the zero-gradient closure superconvergently (the MMS driver problems/advection_2d_x_zero_grad_central.esm uses cos(pi x), whose x-derivative is zero on both walls). Matches plain D(f, wrt: x) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the first x-derivative of a bare field declared over this grid's [x, y] index sets (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axes together. A model imports exactly one first-derivative rule for x (this zero-gradient rule or the periodic variant), never both.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided boundary treatment).
- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

**Conformance fixtures**

- ast: [`central_D1x_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1x_zero_grad_bc) — status `active`
- convergence: [`advection_2d_x_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_zero_grad_central) — status `active`
- simulation: [`advection_2d_x_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_x_zero_grad_central) — status `active`

#### Convergence — case [`advection_2d_x_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_zero_grad_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_x_zero_grad_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_x_zero_grad_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 6.276e-04 | — | 6.034e-04 | — |
| 32 | 1.527e-04 | 2.04 | 1.513e-04 | 2.00 |
| 64 | 3.799e-05 | 2.01 | 3.782e-05 | 2.00 |
| 128 | 9.494e-06 | 2.00 | 9.458e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_x_zero_grad_central (log-log)](/plots/convergence/advection_2d_x_zero_grad_central.png)

### `central_D1y_dirichlet_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1y_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1y_dirichlet_bc.esm)

**Rewrites** `D(f)/Dy` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D1y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((f[i, 2] + f[i, 1] - 2 * gB) / (2 * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((2 * gT - f[i, NY] - f[i, NY - 1]) / (2 * dy)) where {i in x, j in NY:NY})</div>

Second-order centered first derivative along y with inhomogeneous DIRICHLET (fixed-value) conditions on the two y-walls, over the full 2-D field f[x,y] on cartesian_uniform_2d, cell-centered — the first-derivative companion of central_D2y_dirichlet_bc (same mirror ghost, different stencil). Interior j in [2, NY-1] (all i): the imported 3-point centered y-stencil (f[i,j+1]-f[i,j-1])/(2 dy). Faces j=1 and j=NY: mirror-ghost elimination against a prescribed wall value. `gB` and `gT` are consumer-supplied free names (the same free-name contract as the grid geometry `dy`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the BOTTOM y-wall y=y0 (the j=1 side) and the TOP y-wall y=y0+NY*dy (the j=NY side), constant along that wall (all i). On this cell-centered grid the wall sits half a cell outside the first center, so linear reconstruction to the wall gives the ghost f[i,0]=2*gB-f[i,1] (and f[i,NY+1]=2*gT-f[i,NY]); substituting into the centered stencil (f[i,2]-f[i,0])/(2 dy) and (f[i,NY+1]-f[i,NY-1])/(2 dy) yields D1y[i,1] = (f[i,2]+f[i,1]-2*gB)/(2 dy) and D1y[i,NY] = (2*gT-f[i,NY]-f[i,NY-1])/(2 dy). A consuming model declares gB/gT as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec §6: default, or 0), so gB=gT=0 is the homogeneous-Dirichlet y-wall and both wall terms drop out. The gB/gT terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). `dy` is the grid's consumer-supplied y-spacing free name. The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_dirichlet_bc and central_D1y_zero_grad_bc, so every binding's simulator fills the (NX x 1) regions without rank promotion. FACE TRUNCATION: D1y[i,1] - u_y(x_i,y_1) = (dy/8) u_yy(x_i,y0) + O(dy^2), so the mirror-ghost face is O(dy) locally for a generic field (SUPERCONVERGENT to O(dy^2) when u_yy=0 at the wall). WELL-POSED USE: a fixed value pins the field at the wall, the DIFFUSION operator's natural boundary condition, so this gradient closure belongs in a well-posed advection-DIFFUSION setting (or flow into the wall), NOT a pure first-order advection with a fixed value at an OUTFLOW wall, where the pinned value reflects (a Kreiss/GKS-unstable mode); the MMS driver problems/advection_2d_y_dirichlet_central.esm exercises it as the advective y-term of an advection-diffusion decay problem whose physical y-diffusion supplies the boundary dissipation, confirming 2nd order. Matches plain D(f, wrt: y) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the first y-derivative of a bare field declared over this 2-D grid's [x, y] index sets; under import-edge renaming (esm-spec §9.7.7) the wrt literal and the where shape follow the renamed axes together. A model imports exactly one first-derivative rule for y (this Dirichlet rule, the zero-gradient rule, or the periodic variant), never more than one.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (cell-centered Dirichlet via mirror ghost) and Chapter 10 (well-posed boundary conditions for advection-diffusion).
- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

**Conformance fixtures**

- ast: [`central_D1y_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1y_dirichlet_bc) — status `active`
- convergence: [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_dirichlet_central) — status `active`
- simulation: [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_y_dirichlet_central) — status `active`

#### Convergence — case [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_dirichlet_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_y_dirichlet_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_y_dirichlet_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.690e-04 | — | 6.559e-03 | — |
| 32 | 1.442e-04 | 1.98 | 1.829e-03 | 1.84 |
| 64 | 3.617e-05 | 2.00 | 4.809e-04 | 1.93 |
| 128 | 9.050e-06 | 2.00 | 1.232e-04 | 1.96 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_y_dirichlet_central (log-log)](/plots/convergence/advection_2d_y_dirichlet_central.png)

### `central_D1y_periodic`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1y_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1y_periodic.esm)

**Rewrites** `D(f)/Dy` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D1y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((f[i, 2] - f[i, NY]) / (2 * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((f[i, 1] - f[i, NY - 1]) / (2 * dy)) where {i in x, j in NY:NY})</div>

Second-order centered first derivative along y with periodic (wrap) boundary on cartesian_uniform_2d, over the full 2-D field f[x,y]. Interior j in [2, NY-1] (all i): the imported y-stencil (f[i,j+1]-f[i,j-1])/(2 dy). The two y-face rows wrap periodically: j=1 uses the ghost f[i,0]=f[i,NY], giving (f[i,2]-f[i,NY])/(2 dy); j=NY uses f[i,NY+1]=f[i,1], giving (f[i,1]-f[i,NY-1])/(2 dy). The scheme is the same O(dy^2) centered difference at every row and conserves sum_y(f) exactly (the wrap telescopes). The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring grids/cartesian_uniform_2d/rules/central_D2y_neumann_bc so every binding's simulator fills the (NX x 1) regions without rank promotion — the rank-reduced face form of the esm-spec 9.6.8 worked example only fills leading-axis faces portably. `dy` is the grid's consumer-supplied y-spacing free name; the periodic wrap is dy-scale-invariant, so the rule serves any y-extent unchanged. Matches plain D(f, wrt: y) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the first y-derivative of a bare field declared over this grid's [x, y] index sets; under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. O(dy^2). A model imports exactly one first-derivative rule for y (this centered periodic rule or the zero-gradient variant), never both.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D1y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1y_periodic) — status `active`
- convergence: [`advection_2d_y_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_periodic_central) — status `active`
- simulation: [`advection_2d_y_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_y_periodic_central) — status `active`

#### Convergence — case [`advection_2d_y_periodic_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_periodic_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_y_periodic_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_y_periodic_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.578e-02 | — | 1.568e-02 | — |
| 32 | 3.964e-03 | 1.99 | 3.964e-03 | 1.98 |
| 64 | 9.924e-04 | 2.00 | 9.919e-04 | 2.00 |
| 128 | 2.482e-04 | 2.00 | 2.482e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_y_periodic_central (log-log)](/plots/convergence/advection_2d_y_periodic_central.png)

### `central_D1y_zero_grad_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:x,y</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D1y_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D1y_zero_grad_bc.esm)

**Rewrites** `D(f)/Dy` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D1y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((-3 * f[i, 1] + 4 * f[i, 2] - f[i, 3]) / (2 * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((f[i, NY - 2] - 4 * f[i, NY - 1] + 3 * f[i, NY]) / (2 * dy)) where {i in x, j in NY:NY})</div>

Second-order centered first derivative along y with one-sided (zero-gradient frame) boundary rows on cartesian_uniform_2d, over the full 2-D field f[x,y]. Interior j in [2, NY-1] (all i): the imported y-stencil (f[i,j+1]-f[i,j-1])/(2 dy). The two y-face rows use 2nd-order one-sided differences that reach only into the interior (no ghost, no periodic wrap): j=1 gives (-3 f[i,1]+4 f[i,2]-f[i,3])/(2 dy) and j=NY gives (f[i,NY-2]-4 f[i,NY-1]+3 f[i,NY])/(2 dy), so the discretized derivative is fully defined with its BC and stays O(dy^2) at every row (esm-spec 9.6.8). The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring grids/cartesian_uniform_2d/rules/central_D2y_neumann_bc so every binding's simulator fills the (NX x 1) regions without rank promotion — the rank-reduced face form of the esm-spec 9.6.8 worked example only fills leading-axis faces portably. `dy` is the grid's consumer-supplied y-spacing free name. The one-sided faces are O(dy^2) locally and preserve global 2nd-order convergence for boundary-compatible solutions — those whose odd y-derivatives vanish at the walls, so the discrete derivative matches the zero-gradient closure superconvergently (the MMS driver problems/advection_2d_y_zero_grad_central.esm uses cos(pi y), whose y-derivative is zero on both walls). Matches plain D(f, wrt: y) at default priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the first y-derivative of a bare field declared over this grid's [x, y] index sets; under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axes together. A model imports exactly one first-derivative rule for y (this zero-gradient rule or the periodic variant), never both.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided boundary treatment).
- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)

**Conformance fixtures**

- ast: [`central_D1y_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1y_zero_grad_bc) — status `active`
- convergence: [`advection_2d_y_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_zero_grad_central) — status `active`
- simulation: [`advection_2d_y_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_y_zero_grad_central) — status `active`

#### Convergence — case [`advection_2d_y_zero_grad_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_zero_grad_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_y_zero_grad_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_y_zero_grad_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 6.276e-04 | — | 6.034e-04 | — |
| 32 | 1.527e-04 | 2.04 | 1.513e-04 | 2.00 |
| 64 | 3.799e-05 | 2.01 | 3.782e-05 | 2.00 |
| 128 | 9.494e-06 | 2.00 | 9.458e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_y_zero_grad_central (log-log)](/plots/convergence/advection_2d_y_zero_grad_central.png)

### `central_D2x_4_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_4_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_4_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([3:NX - 2, 1:NY] = central_D2x_4_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] (central_D2x_4_dirichlet_face1&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], g=gL&gt;) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] (central_D2x_4_dirichlet_face2&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], g=gL&gt;) where {i in 2:2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] (central_D2x_4_dirichlet_face2&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], g=gR&gt;) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] (central_D2x_4_dirichlet_face1&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], g=gR&gt;) where {i in NX:NX, j in y})</div>

The x-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two x-walls, over the full 2-D field f[x,y], cell-centered. Interior i in [3, NX-2] (all j): the imported 5-point fourth-order x-stencil central_D2x_4_interior, (-f[i-2,j] + 16 f[i-1,j] - 30 f[i,j] + 16 f[i+1,j] - f[i+2,j]) / (12 dx^2), O(dx^4). The four near-x-wall columns i in {1, 2, NX-1, NX}, whose centered 5-point stencils would reach outside the domain, are closed with the imported one-sided fourth-order face templates central_D2x_4_dirichlet_face1 (nearest-wall column) and central_D2x_4_dirichlet_face2 (second column), each applied as a FULL-RANK aggregate over all j (singleton output range on x): the unique degree-5 polynomial through the prescribed x-wall value and the five nearest cells (at fixed j), differentiated twice at that cell. The wall values `gL` (left wall x=x0) and `gR` (right wall x=x0+NX*dx) are consumer-supplied free names bound into the face templates (the g parameter); same free-name contract as the grid geometry `dx` (grids/cartesian_uniform_2d/grid.esm). A consuming model declares gL/gR as ordinary real parameters, constant along the wall (all j); an un-overridden parameter falls back to its default (esm-spec 6: default, or 0), so gL=gR=0 is the homogeneous-Dirichlet case. gL/gR appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The right x-wall reuses the same face weights by reflection (the second derivative is even under x -> -x), so f1..f5 there map to cells NX, NX-1, NX-2, NX-3, NX-4 (increasing distance from the right wall) and g binds to gR. Because BOTH the interior stencil and all four boundary closures are O(dx^4), the x-operator is uniformly fourth order with no boundary reduction; a steady linear-in-x profile carrying the wall values is discretely exact at every column (all closure weight-sums are zero). The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NX >= 5. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2x_4_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_4_dirichlet_bc) — status `active`
- convergence: [`heat_2d_x_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet_o4) — status `active`
- simulation: [`heat_2d_x_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_x_dirichlet_o4) — status `active`

#### Convergence — case [`heat_2d_x_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_dirichlet_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_dirichlet_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 9.756e-07 | — | 2.500e-06 | — |
| 16 | 6.667e-08 | 3.87 | 1.601e-07 | 3.96 |
| 32 | 4.225e-09 | 3.98 | 1.007e-08 | 3.99 |
| 64 | 2.647e-10 | 4.00 | 6.301e-10 | 4.00 |

Expected order: **4** (± 0.2).

![Convergence of heat_2d_x_dirichlet_o4 (log-log)](/plots/convergence/heat_2d_x_dirichlet_o4.png)

### `central_D2x_4_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_4_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_4_neumann_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([3:NX - 2, 1:NY] = central_D2x_4_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] (central_D2x_4_neumann_face1&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], q=qXL&gt;) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] (central_D2x_4_neumann_face2&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], q=qXL&gt;) where {i in 2:2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] (central_D2x_4_neumann_face2&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], q=-1 * qXR&gt;) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] (central_D2x_4_neumann_face1&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], q=-1 * qXR&gt;) where {i in NX:NX, j in y})</div>

The x-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous NEUMANN (fixed-flux) conditions on the two x-walls, over the full 2-D field f[x,y], cell-centered. Interior i in [3, NX-2] (all j): the imported 5-point fourth-order x-stencil central_D2x_4_interior, (-f[i-2,j] + 16 f[i-1,j] - 30 f[i,j] + 16 f[i+1,j] - f[i+2,j]) / (12 dx^2), O(dx^4). The four near-x-wall columns i in {1, 2, NX-1, NX}, whose centered 5-point stencils would reach outside the domain, are closed with the imported one-sided fourth-order face templates central_D2x_4_neumann_face1 (nearest-wall column) and central_D2x_4_neumann_face2 (second column), each applied as a FULL-RANK aggregate over all j (singleton output range on x): the unique degree-5 Hermite-type polynomial matching the prescribed x-wall derivative and the five nearest cells (at fixed j), differentiated twice at that cell. The wall gradients `qXL` (df/dx at the left wall x=x0, +x direction) and `qXR` (df/dx at the right wall x=x0+NX*dx, +x direction) are consumer-supplied free names; same free-name contract as the grid geometry `dx` (grids/cartesian_uniform_2d/grid.esm). A consuming model declares qXL/qXR as ordinary real parameters, constant along the wall (all j); an un-overridden parameter falls back to its default (esm-spec 6: default, or 0), so qXL=qXR=0 is the homogeneous (zero-flux) case and the rule reduces exactly to a zero-gradient fourth-order closure. qXL/qXR appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The left face templates take the wall derivative datum q as their +x-into-domain slope, so the left columns bind q -> qXL. The right x-wall reuses the SAME face weights by reflection (the second derivative is even under x -> -x while df/dx is odd), so the right columns bind q -> -qXR and map f1..f5 to cells NX, NX-1, NX-2, NX-3, NX-4 (increasing distance from the right wall). Because BOTH the interior stencil and all four boundary closures are O(dx^4), the x-operator is uniformly fourth order with no boundary reduction; a steady ramp-in-x carrying the wall fluxes is discretely exact at every column (the value-weight sums are zero and the derivative term supplies the exact wall slope). The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term and central_D2y_4_neumann_bc the second, and the two per-axis makearrays superpose to the correct 9-point interior and one-sided edges/corners. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NX >= 5. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (Neumann one-sided high-order closures) and Chapter 3 (multidimensional heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2x_4_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_4_neumann_bc) — status `active`
- convergence: [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o4) — status `active`
- simulation: [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_neumann_flux_o4) — status `active`

#### Convergence — case [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 4.366e-06 | — | 1.413e-05 | — |
| 16 | 2.077e-07 | 4.39 | 1.006e-06 | 3.81 |
| 32 | 9.348e-09 | 4.47 | 5.036e-08 | 4.32 |
| 64 | 5.271e-10 | 4.15 | 1.688e-09 | 4.90 |

Expected order: **4** (± 0.5).

![Convergence of heat_2d_neumann_flux_o4 (log-log)](/plots/convergence/heat_2d_neumann_flux_o4.png)

### `central_D2x_4_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_4_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_4_periodic.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([3:NX - 2, 1:NY] = central_D2x_4_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((-1 * f[NX - 1, j] + 16 * f[NX, j] + -30 * f[1, j] + 16 * f[2, j] + -1 * f[3, j]) / (12 * dx * dx)) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] ((-1 * f[NX, j] + 16 * f[1, j] + -30 * f[2, j] + 16 * f[3, j] + -1 * f[4, j]) / (12 * dx * dx)) where {i in 2:2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] ((-1 * f[NX - 3, j] + 16 * f[NX - 2, j] + -30 * f[NX - 1, j] + 16 * f[NX, j] + -1 * f[1, j]) / (12 * dx * dx)) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((-1 * f[NX - 2, j] + 16 * f[NX - 1, j] + -30 * f[NX, j] + 16 * f[1, j] + -1 * f[2, j]) / (12 * dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with a PERIODIC x-boundary: fourth-order centered second derivative along x over the full 2-D field f[x,y], the x direction wrapped as a circle. Interior i in [3, NX-2] (all j): the imported 5-point x-stencil central_D2x_4_interior, (-f[i-2,j] + 16 f[i-1,j] - 30 f[i,j] + 16 f[i+1,j] - f[i+2,j]) / (12 dx^2). The four x-face columns i in {1, 2, NX-1, NX} carry the SAME 5-point weights with the out-of-range x-neighbors wrapped periodically (i-2, i-1 wrap past 1 to NX-1, NX; i+1, i+2 wrap past NX to 1, 2), so the whole x-operator is the exact circulant fourth-order Laplacian along x — O(dx^4) with no boundary order reduction, the 2-D per-axis analogue of grids/cartesian_uniform_1d/rules/central_D2_4_periodic. The four face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_periodic so every binding's simulator fills the (1 x NY) columns without rank promotion. `dx` is the grid's consumer-supplied x-spacing free name; the periodic wrap is dx-scale-invariant. Requires NX >= 5 for the 5-point stencil to have distinct points. The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D2x_4_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_4_periodic) — status `active`
- convergence: [`heat_2d_x_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic_o4) — status `active`
- simulation: [`heat_2d_x_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_x_periodic_o4) — status `active`

#### Convergence — case [`heat_2d_x_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_periodic_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_periodic_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.580e-04 | — | 1.403e-04 | — |
| 16 | 1.029e-05 | 3.94 | 9.701e-06 | 3.85 |
| 32 | 6.497e-07 | 3.99 | 6.216e-07 | 3.96 |
| 64 | 4.071e-08 | 4.00 | 3.909e-08 | 3.99 |

Expected order: **4** (± 0.2).

![Convergence of heat_2d_x_periodic_o4 (log-log)](/plots/convergence/heat_2d_x_periodic_o4.png)

### `central_D2x_6_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_6_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_6_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([4:NX - 3, 1:NY] = central_D2x_6_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face1&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], g=gL&gt;) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face2&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], g=gL&gt;) where {i in 2:2, j in y}, [3:3, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face3&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], g=gL&gt;) where {i in 3:3, j in y}, [NX - 2:NX - 2, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face3&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], g=gR&gt;) where {i in NX - 2:NX - 2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face2&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], g=gR&gt;) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] (central_D2x_6_dirichlet_face1&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], g=gR&gt;) where {i in NX:NX, j in y})</div>

The x-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET conditions on the two x-walls, over the full 2-D field f[x,y], cell-centered. Interior x in [4, NX-3] (all y): the imported 7-point sixth-order x-stencil central_D2x_6_interior, O(dx^6). The six near-x-wall lines x in {1,2,3,NX-2,NX-1,NX}, whose centered 7-point stencils would reach outside the domain, are closed with the imported one-sided sixth-order face templates central_D2x_6_dirichlet_face1/face2/face3 (nearest, second, third column), each applied as a FULL-RANK aggregate over all y (singleton output range on x): the unique degree-7 polynomial matching the prescribed wall datum and the seven nearest cells (at fixed y), differentiated twice at that cell. The wall values `gL` and `gR` are consumer-supplied free names bound into the face templates (the g parameter); an un-overridden parameter defaults to 0, so gL=gR=0 is the homogeneous-Dirichlet case. The far wall reuses the same face weights by reflection (the second derivative is even), so f1..f7 there map to cells NX, NX-1, ..., NX-6. The wall free names appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Because BOTH the interior stencil and all six boundary closures are O(dx^6), the x-operator is uniformly sixth order with no boundary reduction; a steady linear-in-x profile carrying the wall values is discretely exact at every line (all closure weight-sums are zero). The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the x-term. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NX >= 7. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2x_6_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_6_dirichlet_bc) — status `active`
- convergence: [`heat_2d_x_dirichlet_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet_o6) — status `active`

#### Convergence — case [`heat_2d_x_dirichlet_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_dirichlet_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_dirichlet_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 7.259e-08 | — | 1.873e-07 | — |
| 16 | 4.874e-10 | 7.22 | 9.892e-10 | 7.56 |
| 32 | 6.520e-12 | 6.22 | 1.558e-11 | 5.99 |
| 64 | 1.422e-13 | 5.52 | 1.044e-12 | 3.90 |

Expected order: **6** (± 0.6).

![Convergence of heat_2d_x_dirichlet_o6 (log-log)](/plots/convergence/heat_2d_x_dirichlet_o6.png)

### `central_D2x_6_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_6_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_6_neumann_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([4:NX - 3, 1:NY] = central_D2x_6_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] (central_D2x_6_neumann_face1&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], q=qXL&gt;) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] (central_D2x_6_neumann_face2&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], q=qXL&gt;) where {i in 2:2, j in y}, [3:3, 1:NY] = sum[i, j] (central_D2x_6_neumann_face3&lt;f1=f[1, j], f2=f[2, j], f3=f[3, j], f4=f[4, j], f5=f[5, j], f6=f[6, j], f7=f[7, j], q=qXL&gt;) where {i in 3:3, j in y}, [NX - 2:NX - 2, 1:NY] = sum[i, j] (central_D2x_6_neumann_face3&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], q=-1 * qXR&gt;) where {i in NX - 2:NX - 2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] (central_D2x_6_neumann_face2&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], q=-1 * qXR&gt;) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] (central_D2x_6_neumann_face1&lt;f1=f[NX, j], f2=f[NX - 1, j], f3=f[NX - 2, j], f4=f[NX - 3, j], f5=f[NX - 4, j], f6=f[NX - 5, j], f7=f[NX - 6, j], q=-1 * qXR&gt;) where {i in NX:NX, j in y})</div>

The x-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous NEUMANN conditions on the two x-walls, over the full 2-D field f[x,y], cell-centered. Interior x in [4, NX-3] (all y): the imported 7-point sixth-order x-stencil central_D2x_6_interior, O(dx^6). The six near-x-wall lines x in {1,2,3,NX-2,NX-1,NX}, whose centered 7-point stencils would reach outside the domain, are closed with the imported one-sided sixth-order face templates central_D2x_6_neumann_face1/face2/face3 (nearest, second, third column), each applied as a FULL-RANK aggregate over all y (singleton output range on x): the unique degree-7 polynomial matching the prescribed wall datum and the seven nearest cells (at fixed y), differentiated twice at that cell. The wall gradients `qXL` (df/dx at the low wall) and `qXR` (df/dx at the high wall) are consumer-supplied free names; an un-overridden parameter defaults to 0, so qXL=qXR=0 is the homogeneous (zero-flux) case and the rule reduces exactly to a zero-gradient sixth-order closure. The low face templates take q as the +x-into-domain slope, so the low rows bind q -> qXL. The far wall reuses the SAME face weights by reflection (the second derivative even, df/dx odd), so the far rows bind q -> -qXR and map f1..f7 to cells NX, NX-1, ..., NX-6. The wall free names appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Because BOTH the interior stencil and all six boundary closures are O(dx^6), the x-operator is uniformly sixth order with no boundary reduction; a steady ramp-in-x carrying the wall fluxes is discretely exact at every line (all closure value-weight sums are zero). The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the x-term. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NX >= 7. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2x_6_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_6_neumann_bc) — status `active`
- convergence: [`heat_2d_neumann_flux_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o6) — status `active`

#### Convergence — case [`heat_2d_neumann_flux_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 3.229e-07 | — | 1.019e-06 | — |
| 16 | 5.963e-09 | 5.76 | 2.827e-08 | 5.17 |
| 32 | 6.581e-11 | 6.50 | 4.193e-10 | 6.08 |
| 64 | 6.494e-13 | 6.66 | 4.322e-12 | 6.60 |

Expected order: **6** (± 0.6).

![Convergence of heat_2d_neumann_flux_o6 (log-log)](/plots/convergence/heat_2d_neumann_flux_o6.png)

### `central_D2x_6_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_6_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_6_periodic.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([4:NX - 3, 1:NY] = central_D2x_6_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((2 * f[NX - 2, j] + -27 * f[NX - 1, j] + 270 * f[NX, j] + -490 * f[1, j] + 270 * f[2, j] + -27 * f[3, j] + 2 * f[4, j]) / (180 * dx * dx)) where {i in 1:1, j in y}, [2:2, 1:NY] = sum[i, j] ((2 * f[NX - 1, j] + -27 * f[NX, j] + 270 * f[1, j] + -490 * f[2, j] + 270 * f[3, j] + -27 * f[4, j] + 2 * f[5, j]) / (180 * dx * dx)) where {i in 2:2, j in y}, [3:3, 1:NY] = sum[i, j] ((2 * f[NX, j] + -27 * f[1, j] + 270 * f[2, j] + -490 * f[3, j] + 270 * f[4, j] + -27 * f[5, j] + 2 * f[6, j]) / (180 * dx * dx)) where {i in 3:3, j in y}, [NX - 2:NX - 2, 1:NY] = sum[i, j] ((2 * f[NX - 5, j] + -27 * f[NX - 4, j] + 270 * f[NX - 3, j] + -490 * f[NX - 2, j] + 270 * f[NX - 1, j] + -27 * f[NX, j] + 2 * f[1, j]) / (180 * dx * dx)) where {i in NX - 2:NX - 2, j in y}, [NX - 1:NX - 1, 1:NY] = sum[i, j] ((2 * f[NX - 4, j] + -27 * f[NX - 3, j] + 270 * f[NX - 2, j] + -490 * f[NX - 1, j] + 270 * f[NX, j] + -27 * f[1, j] + 2 * f[2, j]) / (180 * dx * dx)) where {i in NX - 1:NX - 1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((2 * f[NX - 3, j] + -27 * f[NX - 2, j] + 270 * f[NX - 1, j] + -490 * f[NX, j] + 270 * f[1, j] + -27 * f[2, j] + 2 * f[3, j]) / (180 * dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with a PERIODIC x-boundary: sixth-order centered second derivative along x over the full 2-D field f[x,y], the x direction wrapped as a circle. Interior x in [4, NX-3] (all y): the imported 7-point x-stencil central_D2x_6_interior, (2 f[..-3] - 27 f[..-2] + 270 f[..-1] - 490 f[..] + 270 f[..+1] - 27 f[..+2] + 2 f[..+3]) / (180 dx^2). The six x-face lines x in {1,2,3,NX-2,NX-1,NX} carry the SAME 7-point weights with out-of-range x-neighbors wrapped periodically, so the whole x-operator is the exact circulant sixth-order Laplacian along x — O(dx^6) with no boundary order reduction, the 2-D per-axis analogue of grids/cartesian_uniform_1d/rules/central_D2_6_periodic. The six face values are FULL-RANK aggregates (singleton output range on the x axis, all y) so every binding's simulator fills the lines without rank promotion. Requires NX >= 7 for the 7-point stencil to have distinct points. `dx` is the grid's consumer-supplied x-spacing free name; the periodic wrap is dx-scale-invariant. The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the x-term. Compound match on D(D(f,x),x) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. A model imports exactly one second-derivative rule for x.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D2x_6_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_6_periodic) — status `active`
- convergence: [`heat_2d_x_periodic_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic_o6) — status `active`

#### Convergence — case [`heat_2d_x_periodic_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_periodic_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_periodic_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.503e-05 | — | 1.335e-05 | — |
| 16 | 2.524e-07 | 5.90 | 2.380e-07 | 5.81 |
| 32 | 4.016e-09 | 5.97 | 3.842e-09 | 5.95 |
| 64 | 6.303e-11 | 5.99 | 6.052e-11 | 5.99 |

Expected order: **6** (± 0.3).

![Convergence of heat_2d_x_periodic_o6 (log-log)](/plots/convergence/heat_2d_x_periodic_o6.png)

### `central_D2x_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D2x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((f[2, j] - 3 * f[1, j] + 2 * gL) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((f[NX - 1, j] - 3 * f[NX, j] + 2 * gR) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two x-walls, over the full 2-D field f[x,y], cell-centered. Interior i in [2, NX-1] (all j): the imported x-stencil (f[i+1,j] - 2 f[i,j] + f[i-1,j]) / dx^2. Faces i=1 and i=NX: mirror-ghost elimination against a prescribed wall value. `gL` and `gR` are consumer-supplied free names (the same free-name contract as the grid geometry `dx`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the left x-wall x=x0 and the right x-wall x=x0+NX*dx, constant along that wall (all j). On this cell-centered grid the wall sits half a cell outside the first center, so linear reconstruction to the wall gives the ghost f[0,j]=2*gL-f[1,j] (and f[NX+1,j]=2*gR-f[NX,j]); substituting into the 3-point stencil yields D2x[1,j] = (f[2,j]-3*f[1,j]+2*gL)/dx^2 and D2x[NX,j] = (f[NX-1,j]-3*f[NX,j]+2*gR)/dx^2. A consuming model declares gL/gR as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec §6: default, or 0), so gL=gR=0 is the homogeneous-Dirichlet x-wall. The gL/gR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). `dx` is the grid's consumer-supplied x-spacing free name. The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_neumann_bc. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term. Compound match on D(D(f,x),x) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together. Interior O(dx^2); the one-sided mirror-ghost faces preserve global 2nd-order convergence for Dirichlet-compatible solutions.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (boundary conditions; cell-centered Dirichlet via mirror ghost) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2x_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_dirichlet_bc) — status `active`
- convergence: [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_dirichlet_central) — status `active`
- convergence: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet) — status `active`
- convergence: [`heat_2d_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet) — status `active`
- simulation: [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_x_dirichlet_central) — status `active`
- simulation: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_dirichlet) — status `active`
- simulation: [`heat_2d_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_x_dirichlet) — status `active`

#### Convergence — case [`advection_2d_x_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_x_dirichlet_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_x_dirichlet_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_x_dirichlet_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.690e-04 | — | 6.559e-03 | — |
| 32 | 1.442e-04 | 1.98 | 1.829e-03 | 1.84 |
| 64 | 3.617e-05 | 2.00 | 4.809e-04 | 1.93 |
| 128 | 9.050e-06 | 2.00 | 1.232e-04 | 1.96 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_x_dirichlet_central (log-log)](/plots/convergence/advection_2d_x_dirichlet_central.png)

#### Convergence — case [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.657e-03 | — | 2.377e-03 | — |
| 32 | 6.589e-04 | 2.01 | 5.932e-04 | 2.00 |
| 64 | 1.644e-04 | 2.00 | 1.481e-04 | 2.00 |
| 128 | 4.107e-05 | 2.00 | 3.703e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_dirichlet (log-log)](/plots/convergence/anisotropic_diffusion_2d_dirichlet.png)

#### Convergence — case [`heat_2d_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 3.167e-05 | — | 3.121e-05 | — |
| 32 | 7.925e-06 | 2.00 | 7.837e-06 | 1.99 |
| 64 | 1.982e-06 | 2.00 | 1.962e-06 | 2.00 |
| 128 | 4.954e-07 | 2.00 | 4.905e-07 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_x_dirichlet (log-log)](/plots/convergence/heat_2d_x_dirichlet.png)

### `central_D2x_extrapolate_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_extrapolate_bc.esm)

**Rewrites** `D(D(u)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["u"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D2x_interior&lt;f=u&gt;, [1:1, 1:NY] = sum[i, j] ((u[2, j] + -2 * u[1, j] + extrapolate_linear_ghost&lt;a=u[1, j], b=u[2, j]&gt;) / dx^2) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((u[NX - 1, j] + -2 * u[NX, j] + extrapolate_linear_ghost&lt;a=u[NX, j], b=u[NX - 1, j]&gt;) / dx^2) where {i in NX:NX, j in y})</div>

Second-order centred second derivative d2u/dx2 over the full 2-D field u[x, y] on cartesian_uniform_2d with the linear-extrapolation halo closure of NCAR fire_behavior's level-set solver -- the closure its artificial-viscosity term nu |R| ((diffrx - difflx) + ...) sees at the walls. Interior [2, NX-1] x [1, NY]: the imported three-point stencil stencils/central_D2x_interior, (u[i+1,j] - 2 u[i,j] + u[i-1,j])/dx^2 (respectively along y). Wall lines: the same three-point second difference with the out-of-range neighbour replaced by the one-cell guarded linear-extrapolation ghost of stencils/extrapolate_linear_ghost, g = max(2 a - b, a, b) with a the wall cell and b its inward neighbour (fire_behavior's Extrapol_var_at_bdys, BIAS = 1). CONSEQUENCE, stated plainly: with the guard inactive (a >= b, the field increasing toward the wall) the ghost is exactly 2a - b and the wall-cell second difference is (b - 2a + 2a - b)/dx^2 = 0 -- the closure imposes ZERO curvature on the wall cells, which is what the reference code computes and what makes its wall viscosity vanish; it is not a consistent approximation of d2u/dx2 there for a curved field. When the guard is active (a < b) the ghost is b and the wall cell reads (b - 2a + b)/dx^2 = 2 (b - a)/dx^2. This rule is therefore tagged order:2 for its interior and makes no convergence claim for the wall cells (no convergence case is shipped: a whole-field error would be dominated by the wall cells' zeroth-order closure by construction). Match D(D(u, wrt: x), wrt: x) at priority 10, scoped with where {u: {shape: [x, y]}}, exactly like central_D2x_dirichlet_bc; consumers pick the closure by importing this rule instead of the Dirichlet/Neumann/periodic siblings. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; instead each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j], u[i,0], u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the ENO1 (three-point) stencil is ever evaluated at a wall cell, so one ghost suffices, exactly as in the reference code, whose halo is one cell for the same reason. Corner cells use the two edge ghosts of their walls (the reference code's diagonal corner ghost is never read by ENO1). dx, dy are free names (grid geometry contract). Requires NX >= 3 (respectively NY >= 3).

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)

**Conformance fixtures**

- ast: [`central_D2x_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_extrapolate_bc) — status `active`

### `central_D2x_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_neumann_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D2x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((f[i + 1, j] - f[i, j] - qXL * dx) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((f[i - 1, j] - f[i, j] + qXR * dx) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D Laplacian on cartesian_uniform_2d: second-order centered second derivative along x with inhomogeneous Neumann (fixed-flux) conditions on the two x-walls, over the full 2-D field f[x,y]. Interior i in [2, NX-1] (all j): the imported x-stencil. x-left face i=1: ghost elimination (f[i,j]-f[0,j])/dx = qXL gives D2x[1,j] = (f[2,j]-f[1,j]-qXL*dx)/dx^2. x-right face i=NX: (f[NX+1,j]-f[NX,j])/dx = qXR gives D2x[NX,j] = (f[NX-1,j]-f[NX,j]+qXR*dx)/dx^2. `qXL`/`qXR` are consumer-supplied free names (df/dx at the left/right walls, +x direction) declared as ordinary real parameters; an un-overridden parameter falls back to its default (esm-spec §6: default, or 0), so qXL=qXR=0 is the homogeneous (zero-gradient) x-wall. `dx` is the grid's consumer-supplied spacing free name. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term and central_D2y_neumann_bc the second, and the two per-axis makearrays superpose to the correct 5-point interior and one-sided edges/corners. Compound match on D(D(f,x),x) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together. Interior O(dx^2); the one-sided faces preserve global 2nd-order convergence for Neumann-compatible solutions.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 and Chapter 3 (Neumann conditions; multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`laplacian_2d_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/laplacian_2d_neumann_bc) — status `active`
- convergence: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann) — status `active`
- convergence: [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux) — status `active`
- simulation: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_neumann) — status `active`
- simulation: [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_neumann_flux) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 9.382e-03 | — | 6.067e-03 | — |
| 32 | 2.347e-03 | 2.00 | 1.520e-03 | 2.00 |
| 64 | 5.868e-04 | 2.00 | 3.802e-04 | 2.00 |
| 128 | 1.467e-04 | 2.00 | 9.507e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_neumann (log-log)](/plots/convergence/anisotropic_diffusion_2d_neumann.png)

#### Convergence — case [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.300e-04 | — | 2.380e-04 | — |
| 16 | 3.261e-05 | 2.00 | 6.150e-05 | 1.95 |
| 32 | 8.159e-06 | 2.00 | 1.550e-05 | 1.99 |
| 64 | 2.040e-06 | 2.00 | 3.883e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_neumann_flux (log-log)](/plots/convergence/heat_2d_neumann_flux.png)

### `central_D2x_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2x_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2x_periodic.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = central_D2x_interior&lt;f=f&gt;, [1:1, 1:NY] = sum[i, j] ((f[2, j] + -2 * f[1, j] + f[NX, j]) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((f[1, j] + -2 * f[NX, j] + f[NX - 1, j]) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D Laplacian on cartesian_uniform_2d with a PERIODIC x-boundary: second-order centered second derivative along x over the full 2-D field f[x,y], the x direction wrapped as a circle. Interior i in [2, NX-1] (all j): the imported x-stencil (f[i+1,j] - 2 f[i,j] + f[i-1,j]) / dx^2. The two x-face columns wrap periodically: i=1 uses the ghost f[0,j]=f[NX,j], giving (f[2,j] - 2 f[1,j] + f[NX,j]) / dx^2; i=NX uses f[NX+1,j]=f[1,j], giving (f[1,j] - 2 f[NX,j] + f[NX-1,j]) / dx^2. The scheme is the same O(dx^2) centered second difference at every column, so the whole x-operator is the exact circulant symmetric Laplacian along x — no boundary order reduction, the natural second-order analogue of central_D1x_periodic. The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_neumann_bc so every binding's simulator fills the (1 x NY) columns without rank promotion. `dx` is the grid's consumer-supplied x-spacing free name; the periodic wrap is dx-scale-invariant, so the rule serves any x-extent unchanged. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the first term. Compound match on D(D(f,x),x) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the xx-second-derivative of a bare field declared over this 2-D grid's [x, y] axes (distinct from the 1-D shape-[x] rules); under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. O(dx^2). A model imports exactly one second-derivative rule for x (this periodic rule or a Neumann/Dirichlet variant), never more than one.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains) and Section 2.12 (second differences).

**Conformance fixtures**

- ast: [`central_D2x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2x_periodic) — status `active`
- convergence: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic) — status `active`
- convergence: [`heat_2d_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic) — status `active`
- simulation: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_periodic) — status `active`
- simulation: [`heat_2d_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_x_periodic) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.489e-02 | — | 2.396e-02 | — |
| 32 | 6.326e-03 | 1.98 | 6.266e-03 | 1.94 |
| 64 | 1.589e-03 | 1.99 | 1.585e-03 | 1.98 |
| 128 | 3.976e-04 | 2.00 | 3.973e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_periodic (log-log)](/plots/convergence/anisotropic_diffusion_2d_periodic.png)

#### Convergence — case [`heat_2d_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_x_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_x_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_x_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.049e-04 | — | 4.760e-04 | — |
| 32 | 1.267e-04 | 1.99 | 1.212e-04 | 1.97 |
| 64 | 3.170e-05 | 2.00 | 3.044e-05 | 1.99 |
| 128 | 7.927e-06 | 2.00 | 7.617e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_x_periodic (log-log)](/plots/convergence/heat_2d_x_periodic.png)

### `central_D2y_4_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_4_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_4_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 3:NY - 2] = central_D2y_4_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] (central_D2y_4_dirichlet_face1&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], g=gB&gt;) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] (central_D2y_4_dirichlet_face2&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], g=gB&gt;) where {i in x, j in 2:2}, [1:NX, NY - 1:NY - 1] = sum[i, j] (central_D2y_4_dirichlet_face2&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], g=gT&gt;) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] (central_D2y_4_dirichlet_face1&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], g=gT&gt;) where {i in x, j in NY:NY})</div>

The y-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two y-walls, over the full 2-D field f[x,y], cell-centered. Interior j in [3, NY-2] (all i): the imported 5-point fourth-order y-stencil central_D2y_4_interior, (-f[i,j-2] + 16 f[i,j-1] - 30 f[i,j] + 16 f[i,j+1] - f[i,j+2]) / (12 dy^2), O(dy^4). The four near-y-wall rows j in {1, 2, NY-1, NY}, whose centered 5-point stencils would reach outside the domain, are closed with the imported one-sided fourth-order face templates central_D2y_4_dirichlet_face1 (nearest-wall row) and central_D2y_4_dirichlet_face2 (second row), each applied as a FULL-RANK aggregate over all i (singleton output range on y): the unique degree-5 polynomial through the prescribed y-wall value and the five nearest cells (at fixed i), differentiated twice at that cell. The wall values `gB` (bottom wall y=y0) and `gT` (top wall y=y0+NY*dy) are consumer-supplied free names bound into the face templates (the g parameter); same free-name contract as the grid geometry `dy` (grids/cartesian_uniform_2d/grid.esm). A consuming model declares gB/gT as ordinary real parameters, constant along the wall (all i); an un-overridden parameter falls back to its default (esm-spec 6: default, or 0), so gB=gT=0 is the homogeneous-Dirichlet case. gB/gT appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The top y-wall reuses the same face weights by reflection (the second derivative is even under y -> -y), so f1..f5 there map to cells NY, NY-1, NY-2, NY-3, NY-4 (increasing distance from the top wall) and g binds to gT. Because BOTH the interior stencil and all four boundary closures are O(dy^4), the y-operator is uniformly fourth order with no boundary reduction; a steady linear-in-y profile carrying the wall values is discretely exact at every row (all closure weight-sums are zero). The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the second term. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NY >= 5. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2y_4_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_4_dirichlet_bc) — status `active`
- convergence: [`heat_2d_y_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet_o4) — status `active`
- simulation: [`heat_2d_y_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_y_dirichlet_o4) — status `active`

#### Convergence — case [`heat_2d_y_dirichlet_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_dirichlet_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_dirichlet_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 9.756e-07 | — | 2.500e-06 | — |
| 16 | 6.667e-08 | 3.87 | 1.601e-07 | 3.96 |
| 32 | 4.225e-09 | 3.98 | 1.007e-08 | 3.99 |
| 64 | 2.647e-10 | 4.00 | 6.301e-10 | 4.00 |

Expected order: **4** (± 0.2).

![Convergence of heat_2d_y_dirichlet_o4 (log-log)](/plots/convergence/heat_2d_y_dirichlet_o4.png)

### `central_D2y_4_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_4_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_4_neumann_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 3:NY - 2] = central_D2y_4_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] (central_D2y_4_neumann_face1&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], q=qYL&gt;) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] (central_D2y_4_neumann_face2&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], q=qYL&gt;) where {i in x, j in 2:2}, [1:NX, NY - 1:NY - 1] = sum[i, j] (central_D2y_4_neumann_face2&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], q=-1 * qYR&gt;) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] (central_D2y_4_neumann_face1&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], q=-1 * qYR&gt;) where {i in x, j in NY:NY})</div>

The y-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous NEUMANN (fixed-flux) conditions on the two y-walls, over the full 2-D field f[x,y], cell-centered. Interior j in [3, NY-2] (all i): the imported 5-point fourth-order y-stencil central_D2y_4_interior, (-f[i,j-2] + 16 f[i,j-1] - 30 f[i,j] + 16 f[i,j+1] - f[i,j+2]) / (12 dy^2), O(dy^4). The four near-y-wall rows j in {1, 2, NY-1, NY}, whose centered 5-point stencils would reach outside the domain, are closed with the imported one-sided fourth-order face templates central_D2y_4_neumann_face1 (nearest-wall row) and central_D2y_4_neumann_face2 (second row), each applied as a FULL-RANK aggregate over all i (singleton output range on y): the unique degree-5 Hermite-type polynomial matching the prescribed y-wall derivative and the five nearest cells (at fixed i), differentiated twice at that cell. The wall gradients `qYL` (df/dy at the bottom wall y=y0, +y direction) and `qYR` (df/dy at the top wall y=y0+NY*dy, +y direction) are consumer-supplied free names; same free-name contract as the grid geometry `dy` (grids/cartesian_uniform_2d/grid.esm). A consuming model declares qYL/qYR as ordinary real parameters, constant along the wall (all i); an un-overridden parameter falls back to its default (esm-spec 6: default, or 0), so qYL=qYR=0 is the homogeneous (zero-flux) case and the rule reduces exactly to a zero-gradient fourth-order closure. qYL/qYR appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The bottom face templates take the wall derivative datum q as their +y-into-domain slope, so the bottom rows bind q -> qYL. The top y-wall reuses the SAME face weights by reflection (the second derivative is even under y -> -y while df/dy is odd), so the top rows bind q -> -qYR and map f1..f5 to cells NY, NY-1, NY-2, NY-3, NY-4 (increasing distance from the top wall). Because BOTH the interior stencil and all four boundary closures are O(dy^4), the y-operator is uniformly fourth order with no boundary reduction; a steady ramp-in-y carrying the wall fluxes is discretely exact at every row (the value-weight sums are zero and the derivative term supplies the exact wall slope). The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): central_D2x_4_neumann_bc lowers the first term and this rule the second, and the two per-axis makearrays superpose to the correct 9-point interior and one-sided edges/corners. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NY >= 5. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (Neumann one-sided high-order closures) and Chapter 3 (multidimensional heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2y_4_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_4_neumann_bc) — status `active`
- convergence: [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o4) — status `active`
- simulation: [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_neumann_flux_o4) — status `active`

#### Convergence — case [`heat_2d_neumann_flux_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 4.366e-06 | — | 1.413e-05 | — |
| 16 | 2.077e-07 | 4.39 | 1.006e-06 | 3.81 |
| 32 | 9.348e-09 | 4.47 | 5.036e-08 | 4.32 |
| 64 | 5.271e-10 | 4.15 | 1.688e-09 | 4.90 |

Expected order: **4** (± 0.5).

![Convergence of heat_2d_neumann_flux_o4 (log-log)](/plots/convergence/heat_2d_neumann_flux_o4.png)

### `central_D2y_4_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:4</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_4_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_4_periodic.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 3:NY - 2] = central_D2y_4_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((-1 * f[i, NY - 1] + 16 * f[i, NY] + -30 * f[i, 1] + 16 * f[i, 2] + -1 * f[i, 3]) / (12 * dy * dy)) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] ((-1 * f[i, NY] + 16 * f[i, 1] + -30 * f[i, 2] + 16 * f[i, 3] + -1 * f[i, 4]) / (12 * dy * dy)) where {i in x, j in 2:2}, [1:NX, NY - 1:NY - 1] = sum[i, j] ((-1 * f[i, NY - 3] + 16 * f[i, NY - 2] + -30 * f[i, NY - 1] + 16 * f[i, NY] + -1 * f[i, 1]) / (12 * dy * dy)) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] ((-1 * f[i, NY - 2] + 16 * f[i, NY - 1] + -30 * f[i, NY] + 16 * f[i, 1] + -1 * f[i, 2]) / (12 * dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a fourth-order 2-D Laplacian on cartesian_uniform_2d with a PERIODIC y-boundary: fourth-order centered second derivative along y over the full 2-D field f[x,y], the y direction wrapped as a circle. Interior j in [3, NY-2] (all i): the imported 5-point y-stencil central_D2y_4_interior, (-f[i,j-2] + 16 f[i,j-1] - 30 f[i,j] + 16 f[i,j+1] - f[i,j+2]) / (12 dy^2). The four y-face rows j in {1, 2, NY-1, NY} carry the SAME 5-point weights with the out-of-range y-neighbors wrapped periodically (j-2, j-1 wrap past 1 to NY-1, NY; j+1, j+2 wrap past NY to 1, 2), so the whole y-operator is the exact circulant fourth-order Laplacian along y — O(dy^4) with no boundary order reduction, the 2-D per-axis analogue of grids/cartesian_uniform_1d/rules/central_D2_4_periodic. The four face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_periodic so every binding's simulator fills the (NX x 1) rows without rank promotion. `dy` is the grid's consumer-supplied y-spacing free name; the periodic wrap is dy-scale-invariant. Requires NY >= 5 for the 5-point stencil to have distinct points. The full fourth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the second term. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D2y_4_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_4_periodic) — status `active`
- convergence: [`heat_2d_y_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic_o4) — status `active`
- simulation: [`heat_2d_y_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_y_periodic_o4) — status `active`

#### Convergence — case [`heat_2d_y_periodic_o4`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic_o4)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_periodic_o4/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_periodic_o4/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.580e-04 | — | 1.403e-04 | — |
| 16 | 1.029e-05 | 3.94 | 9.701e-06 | 3.85 |
| 32 | 6.497e-07 | 3.99 | 6.216e-07 | 3.96 |
| 64 | 4.071e-08 | 4.00 | 3.909e-08 | 3.99 |

Expected order: **4** (± 0.2).

![Convergence of heat_2d_y_periodic_o4 (log-log)](/plots/convergence/heat_2d_y_periodic_o4.png)

### `central_D2y_6_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_6_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_6_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 4:NY - 3] = central_D2y_6_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] (central_D2y_6_dirichlet_face1&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], g=gB&gt;) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] (central_D2y_6_dirichlet_face2&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], g=gB&gt;) where {i in x, j in 2:2}, [1:NX, 3:3] = sum[i, j] (central_D2y_6_dirichlet_face3&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], g=gB&gt;) where {i in x, j in 3:3}, [1:NX, NY - 2:NY - 2] = sum[i, j] (central_D2y_6_dirichlet_face3&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], g=gT&gt;) where {i in x, j in NY - 2:NY - 2}, [1:NX, NY - 1:NY - 1] = sum[i, j] (central_D2y_6_dirichlet_face2&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], g=gT&gt;) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] (central_D2y_6_dirichlet_face1&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], g=gT&gt;) where {i in x, j in NY:NY})</div>

The y-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET conditions on the two y-walls, over the full 2-D field f[x,y], cell-centered. Interior y in [4, NY-3] (all x): the imported 7-point sixth-order y-stencil central_D2y_6_interior, O(dy^6). The six near-y-wall lines y in {1,2,3,NY-2,NY-1,NY}, whose centered 7-point stencils would reach outside the domain, are closed with the imported one-sided sixth-order face templates central_D2y_6_dirichlet_face1/face2/face3 (nearest, second, third row), each applied as a FULL-RANK aggregate over all x (singleton output range on y): the unique degree-7 polynomial matching the prescribed wall datum and the seven nearest cells (at fixed x), differentiated twice at that cell. The wall values `gB` and `gT` are consumer-supplied free names bound into the face templates (the g parameter); an un-overridden parameter defaults to 0, so gB=gT=0 is the homogeneous-Dirichlet case. The far wall reuses the same face weights by reflection (the second derivative is even), so f1..f7 there map to cells NY, NY-1, ..., NY-6. The wall free names appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Because BOTH the interior stencil and all six boundary closures are O(dy^6), the y-operator is uniformly sixth order with no boundary reduction; a steady linear-in-y profile carrying the wall values is discretely exact at every line (all closure weight-sums are zero). The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the y-term. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NY >= 7. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2y_6_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_6_dirichlet_bc) — status `active`
- convergence: [`heat_2d_y_dirichlet_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet_o6) — status `active`

#### Convergence — case [`heat_2d_y_dirichlet_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_dirichlet_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_dirichlet_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 7.259e-08 | — | 1.873e-07 | — |
| 16 | 4.874e-10 | 7.22 | 9.892e-10 | 7.56 |
| 32 | 6.520e-12 | 6.22 | 1.558e-11 | 5.99 |
| 64 | 1.422e-13 | 5.52 | 1.044e-12 | 3.90 |

Expected order: **6** (± 0.6).

![Convergence of heat_2d_y_dirichlet_o6 (log-log)](/plots/convergence/heat_2d_y_dirichlet_o6.png)

### `central_D2y_6_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_6_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_6_neumann_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 4:NY - 3] = central_D2y_6_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] (central_D2y_6_neumann_face1&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], q=qYL&gt;) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] (central_D2y_6_neumann_face2&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], q=qYL&gt;) where {i in x, j in 2:2}, [1:NX, 3:3] = sum[i, j] (central_D2y_6_neumann_face3&lt;f1=f[i, 1], f2=f[i, 2], f3=f[i, 3], f4=f[i, 4], f5=f[i, 5], f6=f[i, 6], f7=f[i, 7], q=qYL&gt;) where {i in x, j in 3:3}, [1:NX, NY - 2:NY - 2] = sum[i, j] (central_D2y_6_neumann_face3&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], q=-1 * qYR&gt;) where {i in x, j in NY - 2:NY - 2}, [1:NX, NY - 1:NY - 1] = sum[i, j] (central_D2y_6_neumann_face2&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], q=-1 * qYR&gt;) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] (central_D2y_6_neumann_face1&lt;f1=f[i, NY], f2=f[i, NY - 1], f3=f[i, NY - 2], f4=f[i, NY - 3], f5=f[i, NY - 4], f6=f[i, NY - 5], f7=f[i, NY - 6], q=-1 * qYR&gt;) where {i in x, j in NY:NY})</div>

The y-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with inhomogeneous NEUMANN conditions on the two y-walls, over the full 2-D field f[x,y], cell-centered. Interior y in [4, NY-3] (all x): the imported 7-point sixth-order y-stencil central_D2y_6_interior, O(dy^6). The six near-y-wall lines y in {1,2,3,NY-2,NY-1,NY}, whose centered 7-point stencils would reach outside the domain, are closed with the imported one-sided sixth-order face templates central_D2y_6_neumann_face1/face2/face3 (nearest, second, third row), each applied as a FULL-RANK aggregate over all x (singleton output range on y): the unique degree-7 polynomial matching the prescribed wall datum and the seven nearest cells (at fixed x), differentiated twice at that cell. The wall gradients `qYL` (df/dy at the low wall) and `qYR` (df/dy at the high wall) are consumer-supplied free names; an un-overridden parameter defaults to 0, so qYL=qYR=0 is the homogeneous (zero-flux) case and the rule reduces exactly to a zero-gradient sixth-order closure. The low face templates take q as the +y-into-domain slope, so the low rows bind q -> qYL. The far wall reuses the SAME face weights by reflection (the second derivative even, df/dy odd), so the far rows bind q -> -qYR and map f1..f7 to cells NY, NY-1, ..., NY-6. The wall free names appear only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Because BOTH the interior stencil and all six boundary closures are O(dy^6), the y-operator is uniformly sixth order with no boundary reduction; a steady ramp-in-y carrying the wall fluxes is discretely exact at every line (all closure value-weight sums are zero). The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the y-term. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. Requires NY >= 7. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (one-sided high-order closures) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2y_6_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_6_neumann_bc) — status `active`
- convergence: [`heat_2d_neumann_flux_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o6) — status `active`

#### Convergence — case [`heat_2d_neumann_flux_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 3.229e-07 | — | 1.019e-06 | — |
| 16 | 5.963e-09 | 5.76 | 2.827e-08 | 5.17 |
| 32 | 6.581e-11 | 6.50 | 4.193e-10 | 6.08 |
| 64 | 6.494e-13 | 6.66 | 4.322e-12 | 6.60 |

Expected order: **6** (± 0.6).

![Convergence of heat_2d_neumann_flux_o6 (log-log)](/plots/convergence/heat_2d_neumann_flux_o6.png)

### `central_D2y_6_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:6</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_6_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_6_periodic.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 4:NY - 3] = central_D2y_6_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((2 * f[i, NY - 2] + -27 * f[i, NY - 1] + 270 * f[i, NY] + -490 * f[i, 1] + 270 * f[i, 2] + -27 * f[i, 3] + 2 * f[i, 4]) / (180 * dy * dy)) where {i in x, j in 1:1}, [1:NX, 2:2] = sum[i, j] ((2 * f[i, NY - 1] + -27 * f[i, NY] + 270 * f[i, 1] + -490 * f[i, 2] + 270 * f[i, 3] + -27 * f[i, 4] + 2 * f[i, 5]) / (180 * dy * dy)) where {i in x, j in 2:2}, [1:NX, 3:3] = sum[i, j] ((2 * f[i, NY] + -27 * f[i, 1] + 270 * f[i, 2] + -490 * f[i, 3] + 270 * f[i, 4] + -27 * f[i, 5] + 2 * f[i, 6]) / (180 * dy * dy)) where {i in x, j in 3:3}, [1:NX, NY - 2:NY - 2] = sum[i, j] ((2 * f[i, NY - 5] + -27 * f[i, NY - 4] + 270 * f[i, NY - 3] + -490 * f[i, NY - 2] + 270 * f[i, NY - 1] + -27 * f[i, NY] + 2 * f[i, 1]) / (180 * dy * dy)) where {i in x, j in NY - 2:NY - 2}, [1:NX, NY - 1:NY - 1] = sum[i, j] ((2 * f[i, NY - 4] + -27 * f[i, NY - 3] + 270 * f[i, NY - 2] + -490 * f[i, NY - 1] + 270 * f[i, NY] + -27 * f[i, 1] + 2 * f[i, 2]) / (180 * dy * dy)) where {i in x, j in NY - 1:NY - 1}, [1:NX, NY:NY] = sum[i, j] ((2 * f[i, NY - 3] + -27 * f[i, NY - 2] + 270 * f[i, NY - 1] + -490 * f[i, NY] + 270 * f[i, 1] + -27 * f[i, 2] + 2 * f[i, 3]) / (180 * dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a sixth-order 2-D Laplacian on cartesian_uniform_2d with a PERIODIC y-boundary: sixth-order centered second derivative along y over the full 2-D field f[x,y], the y direction wrapped as a circle. Interior y in [4, NY-3] (all x): the imported 7-point y-stencil central_D2y_6_interior, (2 f[..-3] - 27 f[..-2] + 270 f[..-1] - 490 f[..] + 270 f[..+1] - 27 f[..+2] + 2 f[..+3]) / (180 dy^2). The six y-face lines y in {1,2,3,NY-2,NY-1,NY} carry the SAME 7-point weights with out-of-range y-neighbors wrapped periodically, so the whole y-operator is the exact circulant sixth-order Laplacian along y — O(dy^6) with no boundary order reduction, the 2-D per-axis analogue of grids/cartesian_uniform_1d/rules/central_D2_6_periodic. The six face values are FULL-RANK aggregates (singleton output range on the y axis, all x) so every binding's simulator fills the lines without rank promotion. Requires NY >= 7 for the 7-point stencil to have distinct points. `dy` is the grid's consumer-supplied y-spacing free name; the periodic wrap is dy-scale-invariant. The full sixth-order Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the y-term. Compound match on D(D(f,y),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. A model imports exactly one second-derivative rule for y.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`central_D2y_6_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_6_periodic) — status `active`
- convergence: [`heat_2d_y_periodic_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic_o6) — status `active`

#### Convergence — case [`heat_2d_y_periodic_o6`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic_o6)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_periodic_o6/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_periodic_o6/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.503e-05 | — | 1.335e-05 | — |
| 16 | 2.524e-07 | 5.90 | 2.380e-07 | 5.81 |
| 32 | 4.016e-09 | 5.97 | 3.842e-09 | 5.95 |
| 64 | 6.303e-11 | 5.99 | 6.052e-11 | 5.99 |

Expected order: **6** (± 0.3).

![Convergence of heat_2d_y_periodic_o6 (log-log)](/plots/convergence/heat_2d_y_periodic_o6.png)

### `central_D2y_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D2y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((f[i, 2] - 3 * f[i, 1] + 2 * gB) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((f[i, NY - 1] - 3 * f[i, NY] + 2 * gT) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D Laplacian on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two y-walls, over the full 2-D field f[x,y], cell-centered. Interior j in [2, NY-1] (all i): the imported y-stencil (f[i,j+1] - 2 f[i,j] + f[i,j-1]) / dy^2. Faces j=1 and j=NY: mirror-ghost elimination against a prescribed wall value. `gB` and `gT` are consumer-supplied free names (the same free-name contract as the grid geometry `dy`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the BOTTOM y-wall y=y0 (the j=1 side) and the TOP y-wall y=y0+NY*dy (the j=NY side), constant along that wall (all i). The names deliberately read gB/gT (bottom/top) so the y-wall assignment is unambiguous, the y-analogue of gL/gR on the x-walls. On this cell-centered grid the wall sits half a cell outside the first center, so linear reconstruction to the wall gives the ghost f[i,0]=2*gB-f[i,1] (and f[i,NY+1]=2*gT-f[i,NY]); substituting into the 3-point stencil yields D2y[i,1] = (f[i,2]-3*f[i,1]+2*gB)/dy^2 and D2y[i,NY] = (f[i,NY-1]-3*f[i,NY]+2*gT)/dy^2. A consuming model declares gB/gT as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec §6: default, or 0), so gB=gT=0 is the homogeneous-Dirichlet y-wall. The gB/gT terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). `dy` is the grid's consumer-supplied y-spacing free name. The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_neumann_bc. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the second term. Compound match on D(D(f,y),y) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together. Interior O(dy^2); the one-sided mirror-ghost faces preserve global 2nd-order convergence for Dirichlet-compatible solutions.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (boundary conditions; cell-centered Dirichlet via mirror ghost) and Chapter 3 (multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- ast: [`central_D2y_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_dirichlet_bc) — status `active`
- convergence: [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_dirichlet_central) — status `active`
- convergence: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet) — status `active`
- convergence: [`heat_2d_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet) — status `active`
- simulation: [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_2d_y_dirichlet_central) — status `active`
- simulation: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_dirichlet) — status `active`
- simulation: [`heat_2d_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_y_dirichlet) — status `active`

#### Convergence — case [`advection_2d_y_dirichlet_central`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_2d_y_dirichlet_central)

Error norms read from the committed golden [`tests/conformance/convergence/advection_2d_y_dirichlet_central/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_2d_y_dirichlet_central/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.690e-04 | — | 6.559e-03 | — |
| 32 | 1.442e-04 | 1.98 | 1.829e-03 | 1.84 |
| 64 | 3.617e-05 | 2.00 | 4.809e-04 | 1.93 |
| 128 | 9.050e-06 | 2.00 | 1.232e-04 | 1.96 |

Expected order: **2** (± 0.2).

![Convergence of advection_2d_y_dirichlet_central (log-log)](/plots/convergence/advection_2d_y_dirichlet_central.png)

#### Convergence — case [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.657e-03 | — | 2.377e-03 | — |
| 32 | 6.589e-04 | 2.01 | 5.932e-04 | 2.00 |
| 64 | 1.644e-04 | 2.00 | 1.481e-04 | 2.00 |
| 128 | 4.107e-05 | 2.00 | 3.703e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_dirichlet (log-log)](/plots/convergence/anisotropic_diffusion_2d_dirichlet.png)

#### Convergence — case [`heat_2d_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 3.167e-05 | — | 3.121e-05 | — |
| 32 | 7.925e-06 | 2.00 | 7.837e-06 | 1.99 |
| 64 | 1.982e-06 | 2.00 | 1.962e-06 | 2.00 |
| 128 | 4.954e-07 | 2.00 | 4.905e-07 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_y_dirichlet (log-log)](/plots/convergence/heat_2d_y_dirichlet.png)

### `central_D2y_extrapolate_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_extrapolate_bc.esm)

**Rewrites** `D(D(u)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["u"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D2y_interior&lt;f=u&gt;, [1:NX, 1:1] = sum[i, j] ((u[i, 2] + -2 * u[i, 1] + extrapolate_linear_ghost&lt;a=u[i, 1], b=u[i, 2]&gt;) / dy^2) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((u[i, NY - 1] + -2 * u[i, NY] + extrapolate_linear_ghost&lt;a=u[i, NY], b=u[i, NY - 1]&gt;) / dy^2) where {i in x, j in NY:NY})</div>

Second-order centred second derivative d2u/dy2 over the full 2-D field u[x, y] on cartesian_uniform_2d with the linear-extrapolation halo closure of NCAR fire_behavior's level-set solver -- the closure its artificial-viscosity term nu |R| ((diffry - diffly) + ...) sees at the walls. Interior [1, NX] x [2, NY-1]: the imported three-point stencil stencils/central_D2y_interior, (u[i+1,j] - 2 u[i,j] + u[i-1,j])/dx^2 (respectively along y). Wall lines: the same three-point second difference with the out-of-range neighbour replaced by the one-cell guarded linear-extrapolation ghost of stencils/extrapolate_linear_ghost, g = max(2 a - b, a, b) with a the wall cell and b its inward neighbour (fire_behavior's Extrapol_var_at_bdys, BIAS = 1). CONSEQUENCE, stated plainly: with the guard inactive (a >= b, the field increasing toward the wall) the ghost is exactly 2a - b and the wall-cell second difference is (b - 2a + 2a - b)/dx^2 = 0 -- the closure imposes ZERO curvature on the wall cells, which is what the reference code computes and what makes its wall viscosity vanish; it is not a consistent approximation of d2u/dy2 there for a curved field. When the guard is active (a < b) the ghost is b and the wall cell reads (b - 2a + b)/dx^2 = 2 (b - a)/dx^2. This rule is therefore tagged order:2 for its interior and makes no convergence claim for the wall cells (no convergence case is shipped: a whole-field error would be dominated by the wall cells' zeroth-order closure by construction). Match D(D(u, wrt: y), wrt: y) at priority 10, scoped with where {u: {shape: [x, y]}}, exactly like central_D2y_dirichlet_bc; consumers pick the closure by importing this rule instead of the Dirichlet/Neumann/periodic siblings. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; instead each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j], u[i,0], u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the ENO1 (three-point) stencil is ever evaluated at a wall cell, so one ghost suffices, exactly as in the reference code, whose halo is one cell for the same reason. Corner cells use the two edge ghosts of their walls (the reference code's diagonal corner ghost is never read by ENO1). dx, dy are free names (grid geometry contract). Requires NX >= 3 (respectively NY >= 3).

**References**

- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)

**Conformance fixtures**

- ast: [`central_D2y_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_extrapolate_bc) — status `active`

### `central_D2y_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_neumann_bc.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D2y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((f[i, j + 1] - f[i, j] - qYL * dy) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((f[i, j - 1] - f[i, j] + qYR * dy) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D Laplacian on cartesian_uniform_2d: second-order centered second derivative along y with inhomogeneous Neumann (fixed-flux) conditions on the two y-walls, over the full 2-D field f[x,y]. Interior j in [2, NY-1] (all i): the imported y-stencil. y-bottom face j=1: ghost elimination (f[i,j]-f[i,0])/dy = qYL gives D2y[i,1] = (f[i,2]-f[i,1]-qYL*dy)/dy^2. y-top face j=NY: (f[i,NY+1]-f[i,NY])/dy = qYR gives D2y[i,NY] = (f[i,NY-1]-f[i,NY]+qYR*dy)/dy^2. `qYL`/`qYR` are consumer-supplied free names (df/dy at the bottom/top walls, +y direction) declared as ordinary real parameters; an un-overridden parameter falls back to its default (esm-spec §6: default, or 0), so qYL=qYR=0 is the homogeneous (zero-gradient) y-wall. `dy` is the grid's consumer-supplied spacing free name. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): central_D2x_neumann_bc lowers the first term and this rule the second, and the two per-axis makearrays superpose to the correct 5-point interior and one-sided edges/corners. Compound match on D(D(f,y),y) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together. Interior O(dy^2); the one-sided faces preserve global 2nd-order convergence for Neumann-compatible solutions.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 and Chapter 3 (Neumann conditions; multidimensional Poisson/heat by dimension-wise second differences).

**Conformance fixtures**

- convergence: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann) — status `active`
- convergence: [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux) — status `active`
- simulation: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_neumann) — status `active`
- simulation: [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_neumann_flux) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 9.382e-03 | — | 6.067e-03 | — |
| 32 | 2.347e-03 | 2.00 | 1.520e-03 | 2.00 |
| 64 | 5.868e-04 | 2.00 | 3.802e-04 | 2.00 |
| 128 | 1.467e-04 | 2.00 | 9.507e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_neumann (log-log)](/plots/convergence/anisotropic_diffusion_2d_neumann.png)

#### Convergence — case [`heat_2d_neumann_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_neumann_flux)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_neumann_flux/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_neumann_flux/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 8 | 1.300e-04 | — | 2.380e-04 | — |
| 16 | 3.261e-05 | 2.00 | 6.150e-05 | 1.95 |
| 32 | 8.159e-06 | 2.00 | 1.550e-05 | 1.99 |
| 64 | 2.040e-06 | 2.00 | 3.883e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_neumann_flux (log-log)](/plots/convergence/heat_2d_neumann_flux.png)

### `central_D2y_periodic`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_D2y_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_D2y_periodic.esm)

**Rewrites** `D(D(f)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "y"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = central_D2y_interior&lt;f=f&gt;, [1:NX, 1:1] = sum[i, j] ((f[i, 2] + -2 * f[i, 1] + f[i, NY]) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((f[i, 1] + -2 * f[i, NY] + f[i, NY - 1]) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D Laplacian on cartesian_uniform_2d with a PERIODIC y-boundary: second-order centered second derivative along y over the full 2-D field f[x,y], the y direction wrapped as a circle. Interior j in [2, NY-1] (all i): the imported y-stencil (f[i,j+1] - 2 f[i,j] + f[i,j-1]) / dy^2. The two y-face rows wrap periodically: j=1 uses the ghost f[i,0]=f[i,NY], giving (f[i,2] - 2 f[i,1] + f[i,NY]) / dy^2; j=NY uses f[i,NY+1]=f[i,1], giving (f[i,1] - 2 f[i,NY] + f[i,NY-1]) / dy^2. The scheme is the same O(dy^2) centered second difference at every row, so the whole y-operator is the exact circulant symmetric Laplacian along y — no boundary order reduction, the natural second-order analogue of central_D1y_periodic. The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_neumann_bc so every binding's simulator fills the (NX x 1) rows without rank promotion. `dy` is the grid's consumer-supplied y-spacing free name; the periodic wrap is dy-scale-invariant, so the rule serves any y-extent unchanged. The full Laplacian is assembled in the consuming equation as D(D(f,x),x) + D(D(f,y),y): this rule lowers the second term. Compound match on D(D(f,y),y) at priority 10 (esm-spec §9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec §9.6.1) so it fires only on the yy-second-derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec §9.7.7) the wrt literals and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. O(dy^2). A model imports exactly one second-derivative rule for y (this periodic rule or a Neumann/Dirichlet variant), never more than one.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains) and Section 2.12 (second differences).

**Conformance fixtures**

- ast: [`central_D2y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2y_periodic) — status `active`
- convergence: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic) — status `active`
- convergence: [`heat_2d_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic) — status `active`
- simulation: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_periodic) — status `active`
- simulation: [`heat_2d_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_y_periodic) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.489e-02 | — | 2.396e-02 | — |
| 32 | 6.326e-03 | 1.98 | 6.266e-03 | 1.94 |
| 64 | 1.589e-03 | 1.99 | 1.585e-03 | 1.98 |
| 128 | 3.976e-04 | 2.00 | 3.973e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_periodic (log-log)](/plots/convergence/anisotropic_diffusion_2d_periodic.png)

#### Convergence — case [`heat_2d_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_y_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_y_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_y_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 5.049e-04 | — | 4.760e-04 | — |
| 32 | 1.267e-04 | 1.99 | 1.212e-04 | 1.97 |
| 64 | 3.170e-05 | 2.00 | 3.044e-05 | 1.99 |
| 128 | 7.927e-06 | 2.00 | 7.617e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_y_periodic (log-log)](/plots/convergence/heat_2d_y_periodic.png)

### `central_Dxy_dirichlet_bc`

<p class="rule-meta"><code>op:Dxy</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_Dxy_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_Dxy_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 2:NY - 1] = central_Dxy_interior&lt;f=f&gt;, [1:1, 2:NY - 1] = sum[i, j] ((f[2, j + 1] + f[1, j + 1] - (f[2, j - 1] + f[1, j - 1])) / (4 * dx * dy)) where {i in 1:1, j in 2:NY - 1}, [NX:NX, 2:NY - 1] = sum[i, j] ((f[NX, j - 1] + f[NX - 1, j - 1] - (f[NX, j + 1] + f[NX - 1, j + 1])) / (4 * dx * dy)) where {i in NX:NX, j in 2:NY - 1}, [2:NX - 1, 1:1] = sum[i, j] ((f[i + 1, 2] + f[i + 1, 1] - (f[i - 1, 2] + f[i - 1, 1])) / (4 * dx * dy)) where {i in 2:NX - 1, j in 1:1}, [2:NX - 1, NY:NY] = sum[i, j] ((f[i - 1, NY] + f[i - 1, NY - 1] - (f[i + 1, NY] + f[i + 1, NY - 1])) / (4 * dx * dy)) where {i in 2:NX - 1, j in NY:NY}, [1:1, 1:1] = sum[i, j] ((f[2, 2] + f[2, 1] + f[1, 2] + f[1, 1]) / (4 * dx * dy)) where {i in 1:1, j in 1:1}, [1:1, NY:NY] = sum[i, j] ((0 - (f[2, NY] + f[2, NY - 1] + f[1, NY] + f[1, NY - 1])) / (4 * dx * dy)) where {i in 1:1, j in NY:NY}, [NX:NX, 1:1] = sum[i, j] ((0 - (f[NX, 2] + f[NX, 1] + f[NX - 1, 2] + f[NX - 1, 1])) / (4 * dx * dy)) where {i in NX:NX, j in 1:1}, [NX:NX, NY:NY] = sum[i, j] ((f[NX, NY] + f[NX, NY - 1] + f[NX - 1, NY] + f[NX - 1, NY - 1]) / (4 * dx * dy)) where {i in NX:NX, j in NY:NY})</div>

Second-order centered MIXED second derivative d^2f/dx dy on cartesian_uniform_2d with DIRICHLET (fixed-value) walls on all four sides, over the full 2-D field f[x,y], cell-centered. This is the off-diagonal 2 Kxy d^2f/dx dy term of a full anisotropic diffusion tensor K = [[Kxx, Kxy], [Kxy, Kyy]] that the diagonal central_D2x/central_D2y rules do not supply; the non-periodic wall analogue of central_Dxy_periodic. INTERIOR block [2,NX-1]x[2,NY-1]: the imported central_Dxy_interior 4-point stencil (f[i+1,j+1]-f[i+1,j-1]-f[i-1,j+1]+f[i-1,j-1])/(4 dx dy), which reaches diagonally to all four corner neighbours, so the wall closure needs ghosts across whichever axis is out of the domain (one axis on the four edge strips, two axes plus the diagonal corner ghost on the four corner cells). GHOST RECONSTRUCTION: on this cell-centered grid the wall sits half a cell outside the first center, so the Dirichlet value g_wall gives the mirror (odd-reflection) ghost f_ghost = 2 g_wall - f_interior (left f[0,j]=2 gL-f[1,j], right f[NX+1,j]=2 gR-f[NX,j], bottom f[i,0]=2 gB-f[i,1], top f[i,NY+1]=2 gT-f[i,NY]); the diagonal corner ghost is that reconstruction applied in both directions (for the homogeneous-compatible case f[0,0]=f[1,1], etc.). CANCELLATION OF THE WALL VALUE (why no free names appear): the mixed derivative differences ALONG the wall (the tangential direction), and the prescribed value is constant along each wall (scalar gL/gR/gB/gT, matching central_D2x_dirichlet_bc / central_D2y_dirichlet_bc). At the left x-wall edge (i=1) the stencil's antisymmetric y-difference of the ghost column is -f[0,j+1]+f[0,j-1] = -(2 gL-f[1,j+1])+(2 gL-f[1,j-1]) = f[1,j+1]-f[1,j-1]: the 2 gL drops IDENTICALLY. The same tangential-difference cancellation removes gR, gB, gT on the other three edges and every corner (the diagonal-ghost constant terms cancel against the two edge-ghost constant terms). So this Dirichlet cross closure carries NO wall free names - it encodes only the structural fact that df/dy=0 along the x-walls and df/dx=0 along the y-walls (both implied by f being constant along the wall), which is exactly the boundary information the mixed derivative needs; the specific wall values are irrelevant to d^2f/dx dy. (A wall datum that VARIES along the wall would add a tangential-derivative term this scalar-wall closure omits, consistent with the scalar-wall D2 rules; out of scope.) RESULTING FACE FORMULAS (all over 4 dx dy): left edge i=1 -> (f[2,j+1]-f[2,j-1])+(f[1,j+1]-f[1,j-1]); right edge i=NX -> (f[NX,j-1]-f[NX,j+1])+(f[NX-1,j-1]-f[NX-1,j+1]); bottom edge j=1 -> (f[i+1,2]-f[i-1,2])+(f[i+1,1]-f[i-1,1]); top edge j=NY -> (f[i-1,NY]-f[i+1,NY])+(f[i-1,NY-1]-f[i+1,NY-1]); corner (1,1) -> f[2,2]+f[2,1]+f[1,2]+f[1,1]; corner (1,NY) -> -(f[2,NY]+f[2,NY-1]+f[1,NY]+f[1,NY-1]); corner (NX,1) -> -(f[NX,2]+f[NX,1]+f[NX-1,2]+f[NX-1,1]); corner (NX,NY) -> f[NX,NY]+f[NX,NY-1]+f[NX-1,NY]+f[NX-1,NY-1]. Each edge/corner is the mirror-ghost centered cross difference, O(dx^2+dy^2) like the interior, and is exactly the odd-reflected form of the interior stencil (the odd extension is exact for a solution that is odd about the wall, e.g. sin(pi x) about x=0). The makearray tiles [x,y] as the 3x3 block partition {1}u[2,NX-1]u{NX} x {1}u[2,NY-1]u{NY} (interior + four edge strips + four corner cells), mirroring central_Dxy_periodic and central_D2x/D2y_dirichlet_bc so every binding's simulator fills the strips and corners without rank promotion. `dx`/`dy` are the grid's consumer-supplied spacing free names (grids/cartesian_uniform_2d/grid.esm). WELL-POSEDNESS: the isolated cross operator is symmetric indefinite (eigenvalues span -N^2..+N^2), not a stable time-evolution operator on its own; the accompanying MMS (problems/anisotropic_diffusion_2d_dirichlet.esm) embeds it in the positive-definite tensor K=[[1,1/2],[1/2,1]] so L = D2x+Dxy+D2y is dissipative. With homogeneous Dirichlet walls the discrete energy estimate's boundary flux vanishes (the test function is zero on the wall), so the combined operator is provably negative definite and the parabolic evolution is well-posed - unlike the Neumann cross closure, whose walls leave a residual conormal boundary flux. Compound match on D(D(f,x),y) at priority 10 (esm-spec 9.6.3); the consumer writes d^2u/dx dy as D(D(u,x),y) (inner wrt x, outer wrt y), and the priority-10 compound fires before any first-derivative rule could rewrite the inner D. Match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the mixed second derivative of a bare field declared over this 2-D grid's [x,y] axes; under import-edge renaming (esm-spec 9.7.7) the two wrt literals and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. O(dx^2+dy^2).

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (mixed partial derivatives; cell-centered ghost reconstruction for Dirichlet/Neumann walls) and Chapter 3 (multidimensional operators by tensor products of one-dimensional differences).
- Morton, K. W. and Mayers, D. F. (2005). Numerical Solution of Partial Differential Equations, 2nd ed. Cambridge University Press. Section 3 (boundary closures by ghost cells; tangential vs normal derivatives at walls).

**Conformance fixtures**

- ast: [`central_Dxy_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_Dxy_dirichlet_bc) — status `active`
- convergence: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet) — status `active`
- simulation: [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_dirichlet) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.657e-03 | — | 2.377e-03 | — |
| 32 | 6.589e-04 | 2.01 | 5.932e-04 | 2.00 |
| 64 | 1.644e-04 | 2.00 | 1.481e-04 | 2.00 |
| 128 | 4.107e-05 | 2.00 | 3.703e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_dirichlet (log-log)](/plots/convergence/anisotropic_diffusion_2d_dirichlet.png)

### `central_Dxy_neumann_bc`

<p class="rule-meta"><code>op:Dxy</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_Dxy_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_Dxy_neumann_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 2:NY - 1] = central_Dxy_interior&lt;f=f&gt;, [1:1, 2:NY - 1] = sum[i, j] ((f[2, j + 1] + f[1, j - 1] - (f[2, j - 1] + f[1, j + 1])) / (4 * dx * dy)) where {i in 1:1, j in 2:NY - 1}, [NX:NX, 2:NY - 1] = sum[i, j] ((f[NX, j + 1] + f[NX - 1, j - 1] - (f[NX, j - 1] + f[NX - 1, j + 1])) / (4 * dx * dy)) where {i in NX:NX, j in 2:NY - 1}, [2:NX - 1, 1:1] = sum[i, j] ((f[i + 1, 2] + f[i - 1, 1] - (f[i + 1, 1] + f[i - 1, 2])) / (4 * dx * dy)) where {i in 2:NX - 1, j in 1:1}, [2:NX - 1, NY:NY] = sum[i, j] ((f[i + 1, NY] + f[i - 1, NY - 1] - (f[i + 1, NY - 1] + f[i - 1, NY])) / (4 * dx * dy)) where {i in 2:NX - 1, j in NY:NY}, [1:1, 1:1] = sum[i, j] ((f[2, 2] + f[1, 1] - (f[2, 1] + f[1, 2])) / (4 * dx * dy)) where {i in 1:1, j in 1:1}, [1:1, NY:NY] = sum[i, j] ((f[2, NY] + f[1, NY - 1] - (f[2, NY - 1] + f[1, NY])) / (4 * dx * dy)) where {i in 1:1, j in NY:NY}, [NX:NX, 1:1] = sum[i, j] ((f[NX, 2] + f[NX - 1, 1] - (f[NX, 1] + f[NX - 1, 2])) / (4 * dx * dy)) where {i in NX:NX, j in 1:1}, [NX:NX, NY:NY] = sum[i, j] ((f[NX, NY] + f[NX - 1, NY - 1] - (f[NX, NY - 1] + f[NX - 1, NY])) / (4 * dx * dy)) where {i in NX:NX, j in NY:NY})</div>

Second-order centered MIXED second derivative d^2f/dx dy on cartesian_uniform_2d with NEUMANN (fixed-flux) walls on all four sides, over the full 2-D field f[x,y], cell-centered. This is the off-diagonal 2 Kxy d^2f/dx dy term of a full anisotropic diffusion tensor K = [[Kxx, Kxy], [Kxy, Kyy]] that the diagonal central_D2x/central_D2y rules do not supply; the fixed-flux wall analogue of central_Dxy_periodic and the Neumann counterpart of central_Dxy_dirichlet_bc. INTERIOR block [2,NX-1]x[2,NY-1]: the imported central_Dxy_interior 4-point stencil (f[i+1,j+1]-f[i+1,j-1]-f[i-1,j+1]+f[i-1,j-1])/(4 dx dy), which reaches diagonally, so the wall closure needs ghosts across whichever axis is out of the domain (one axis on the edges, two axes plus the diagonal corner ghost on the corners). GHOST RECONSTRUCTION: the prescribed normal derivative q_wall gives the even-reflection ghost from the one-sided flux relation (left (f[1,j]-f[0,j])/dx=qXL -> f[0,j]=f[1,j]-qXL dx; right f[NX+1,j]=f[NX,j]+qXR dx; bottom f[i,0]=f[i,1]-qYL dy; top f[i,NY+1]=f[i,NY]+qYR dy); the diagonal corner ghost applies the flux relation in both directions and is UNAMBIGUOUS (both orders give f[0,0]=f[1,1]-qXL dx-qYL dy). CANCELLATION OF THE WALL FLUX (why no free names appear): the mixed derivative differences ALONG the wall, and q_wall is constant along each wall (scalar qXL/qXR/qYL/qYR, matching central_D2x_neumann_bc / central_D2y_neumann_bc). At the left x-wall edge (i=1) the antisymmetric y-difference of the ghost column is -f[0,j+1]+f[0,j-1] = -(f[1,j+1]-qXL dx)+(f[1,j-1]-qXL dx) = -f[1,j+1]+f[1,j-1]: the qXL dx drops IDENTICALLY. The same tangential-difference cancellation removes qXR, qYL, qYR on the other edges and every corner (the diagonal-ghost -qXL dx-qYL dy terms cancel against the two edge-ghost flux terms). So this Neumann cross closure carries NO wall free names either; it encodes only the EVEN reflection of the field across each wall (the fixed-flux continuation), the sign-flipped counterpart of the odd-reflection Dirichlet closure. RESULTING FACE FORMULAS (all over 4 dx dy), the interior stencil with each out-of-range sample replaced by its EVEN mirror (reflected index, same sign): left edge i=1 -> (f[2,j+1]-f[2,j-1])-(f[1,j+1]-f[1,j-1]); right edge i=NX -> (f[NX,j+1]-f[NX,j-1])-(f[NX-1,j+1]-f[NX-1,j-1]); bottom edge j=1 -> (f[i+1,2]-f[i-1,2])-(f[i+1,1]-f[i-1,1]); top edge j=NY -> (f[i+1,NY]-f[i-1,NY])-(f[i+1,NY-1]-f[i-1,NY-1]); corner (1,1) -> (f[2,2]+f[1,1])-(f[2,1]+f[1,2]); corner (1,NY) -> (f[2,NY]+f[1,NY-1])-(f[2,NY-1]+f[1,NY]); corner (NX,1) -> (f[NX,2]+f[NX-1,1])-(f[NX,1]+f[NX-1,2]); corner (NX,NY) -> (f[NX,NY]+f[NX-1,NY-1])-(f[NX,NY-1]+f[NX-1,NY]). Each face is the even-reflected centered cross difference, O(dx^2+dy^2) like the interior, exact for a solution even about the wall (e.g. cos(pi x) about x=0). The makearray tiles [x,y] as the 3x3 block partition {1}u[2,NX-1]u{NX} x {1}u[2,NY-1]u{NY}, mirroring central_Dxy_periodic and central_D2x/D2y_neumann_bc. `dx`/`dy` are the grid's consumer-supplied spacing free names (grids/cartesian_uniform_2d/grid.esm). WELL-POSEDNESS CAVEAT: the isolated cross operator is symmetric indefinite; embedding it in the positive-definite tensor K=[[1,1/2],[1/2,1]] (L = D2x+Dxy+D2y) makes the DIRICHLET-walled operator dissipative, but with per-axis Neumann walls (df/dx=0 on x-walls, df/dy=0 on y-walls) the discrete energy estimate retains a residual conormal boundary flux n.K grad u = Kxy u_tangential on each wall that the imposed normal condition does not kill, so the combined operator is NOT guaranteed negative definite and the pure time-evolution MMS is only conditionally well-posed. See problems/anisotropic_diffusion_2d_neumann.esm for the manufactured-solution setup and the observed stability at the tested resolutions. Compound match on D(D(f,x),y) at priority 10 (esm-spec 9.6.3), match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the mixed second derivative of a bare field declared over this 2-D grid's [x,y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axes together. O(dx^2+dy^2).

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (mixed partial derivatives; cell-centered ghost reconstruction for Dirichlet/Neumann walls) and Chapter 3 (multidimensional operators by tensor products of one-dimensional differences).
- Morton, K. W. and Mayers, D. F. (2005). Numerical Solution of Partial Differential Equations, 2nd ed. Cambridge University Press. Section 3 (boundary closures by ghost cells; tangential vs normal derivatives at walls).

**Conformance fixtures**

- ast: [`central_Dxy_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_Dxy_neumann_bc) — status `active`
- convergence: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann) — status `active`
- simulation: [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_neumann) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 9.382e-03 | — | 6.067e-03 | — |
| 32 | 2.347e-03 | 2.00 | 1.520e-03 | 2.00 |
| 64 | 5.868e-04 | 2.00 | 3.802e-04 | 2.00 |
| 128 | 1.467e-04 | 2.00 | 9.507e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_neumann (log-log)](/plots/convergence/anisotropic_diffusion_2d_neumann.png)

### `central_Dxy_periodic`

<p class="rule-meta"><code>op:Dxy</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/central_Dxy_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/central_Dxy_periodic.esm)

**Rewrites** `D(D(f)/Dx)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 2:NY - 1] = central_Dxy_interior&lt;f=f&gt;, [1:1, 2:NY - 1] = sum[i, j] ((f[2, j + 1] - f[2, j - 1] - (f[NX, j + 1] - f[NX, j - 1])) / (4 * dx * dy)) where {i in 1:1, j in 2:NY - 1}, [NX:NX, 2:NY - 1] = sum[i, j] ((f[1, j + 1] - f[1, j - 1] - (f[NX - 1, j + 1] - f[NX - 1, j - 1])) / (4 * dx * dy)) where {i in NX:NX, j in 2:NY - 1}, [2:NX - 1, 1:1] = sum[i, j] ((f[i + 1, 2] - f[i + 1, NY] - (f[i - 1, 2] - f[i - 1, NY])) / (4 * dx * dy)) where {i in 2:NX - 1, j in 1:1}, [2:NX - 1, NY:NY] = sum[i, j] ((f[i + 1, 1] - f[i + 1, NY - 1] - (f[i - 1, 1] - f[i - 1, NY - 1])) / (4 * dx * dy)) where {i in 2:NX - 1, j in NY:NY}, [1:1, 1:1] = sum[i, j] ((f[2, 2] - f[2, NY] - (f[NX, 2] - f[NX, NY])) / (4 * dx * dy)) where {i in 1:1, j in 1:1}, [1:1, NY:NY] = sum[i, j] ((f[2, 1] - f[2, NY - 1] - (f[NX, 1] - f[NX, NY - 1])) / (4 * dx * dy)) where {i in 1:1, j in NY:NY}, [NX:NX, 1:1] = sum[i, j] ((f[1, 2] - f[1, NY] - (f[NX - 1, 2] - f[NX - 1, NY])) / (4 * dx * dy)) where {i in NX:NX, j in 1:1}, [NX:NX, NY:NY] = sum[i, j] ((f[1, 1] - f[1, NY - 1] - (f[NX - 1, 1] - f[NX - 1, NY - 1])) / (4 * dx * dy)) where {i in NX:NX, j in NY:NY})</div>

Second-order centered MIXED second derivative d^2f/dx dy on cartesian_uniform_2d with a DOUBLY-PERIODIC boundary (both x and y wrapped as circles), over the full 2-D field f[x,y]. This is the off-diagonal term of a full anisotropic diffusion tensor K = [[Kxx, Kxy], [Kxy, Kyy]] — the 2 Kxy d^2f/dx dy contribution the diagonal central_D2x/central_D2y rules do not supply. The 4-point stencil (f[i+1,j+1] - f[i+1,j-1] - f[i-1,j+1] + f[i-1,j-1]) / (4 dx dy) reaches diagonally to all four corner neighbours, so the periodic closure needs BOTH x-walls AND both y-walls to wrap, plus the four corners where the two wraps meet. The makearray therefore tiles [x,y] as a 3x3 block partition {1}u[2,NX-1]u{NX} x {1}u[2,NY-1]u{NY}: the doubly-interior block [2,NX-1]x[2,NY-1] holds the imported central_Dxy_interior stencil; the four edge strips wrap the single out-of-range axis (i=1 uses the ghost column NX for i-1, i=NX uses column 1 for i+1; j=1 uses row NY for j-1, j=NY uses row 1 for j+1); the four corner cells wrap both axes. Every cell uses the identical O(dx^2 + dy^2) centered cross difference (the wrap is exact on the periodic torus, no boundary order reduction), so the whole operator is the circulant tensor product D1x (x) D1y — the periodic cross-derivative analogue of central_D1x_periodic/central_D1y_periodic. The edge and corner faces are FULL-RANK aggregates (singleton output range on the wrapped axis) mirroring central_D2x_periodic/central_D2y_periodic so every binding's simulator fills the strips and corner cells without rank promotion. `dx`/`dy` are the grid's consumer-supplied spacing free names; the periodic wrap is scale-invariant, so the rule serves any x/y extent unchanged. Compound match on D(D(f,x),y) at priority 10 (esm-spec 9.6.3) — the consumer writes d^2u/dx dy as D(D(u,x),y) (inner wrt x, outer wrt y); the priority-10 compound fires before any first-derivative rule could rewrite the inner D. Match-scoped with `where: {f: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the mixed second derivative of a bare field declared over this 2-D grid's [x, y] axes; under import-edge renaming (esm-spec 9.7.7) the two wrt literals and the where shape follow the renamed axes together, so two instances coexist scoped to their own meshes. A one-sided (non-periodic) closure of the cross derivative is a follow-up; periodic is the cleanest self-contained boundary for a diagonal-reaching cross stencil (no ghost-gradient bookkeeping on four walls plus four corners), which is why the accompanying MMS (problems/anisotropic_diffusion_2d_periodic.esm) is doubly periodic. NOTE the isolated cross operator is symmetric indefinite (its discrete eigenvalues span -N^2..+N^2), so it is not a well-posed time-evolution operator on its own; the MMS embeds it in a positive-definite diffusion tensor, whose combined operator is dissipative — see the problem file for the derivation.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 10.4 (centered differences on periodic domains) and Section 2.12 (mixed partial derivatives).

**Conformance fixtures**

- ast: [`central_Dxy_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_Dxy_periodic) — status `active`
- convergence: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic) — status `active`
- simulation: [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/anisotropic_diffusion_2d_periodic) — status `active`

#### Convergence — case [`anisotropic_diffusion_2d_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/anisotropic_diffusion_2d_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.489e-02 | — | 2.396e-02 | — |
| 32 | 6.326e-03 | 1.98 | 6.266e-03 | 1.94 |
| 64 | 1.589e-03 | 1.99 | 1.585e-03 | 1.98 |
| 128 | 3.976e-04 | 2.00 | 3.973e-04 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of anisotropic_diffusion_2d_periodic (log-log)](/plots/convergence/anisotropic_diffusion_2d_periodic.png)

### `godunov_norm_D2_dirichlet_bc`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/godunov_norm_D2_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/godunov_norm_D2_dirichlet_bc.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 10 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 2:NY - 1] = godunov_norm_D2_interior&lt;u=u&gt;, [1:1, 2:NY - 1] = sum[i, j] (sqrt(max(((u[1, j] - (2 * gL - u[1, j])) / dx), 0)^2 + min(((u[2, j] - u[1, j]) / dx), 0)^2 + max(((u[1, j] - u[1, j - 1]) / dy), 0)^2 + min(((u[1, j + 1] - u[1, j]) / dy), 0)^2)) where {i in 1:1, j in 2:NY - 1}, [NX:NX, 2:NY - 1] = sum[i, j] (sqrt(max(((u[NX, j] - u[NX - 1, j]) / dx), 0)^2 + min(((2 * gR - u[NX, j] - u[NX, j]) / dx), 0)^2 + max(((u[NX, j] - u[NX, j - 1]) / dy), 0)^2 + min(((u[NX, j + 1] - u[NX, j]) / dy), 0)^2)) where {i in NX:NX, j in 2:NY - 1}, [2:NX - 1, 1:1] = sum[i, j] (sqrt(max(((u[i, 1] - u[i - 1, 1]) / dx), 0)^2 + min(((u[i + 1, 1] - u[i, 1]) / dx), 0)^2 + max(((u[i, 1] - (2 * gB - u[i, 1])) / dy), 0)^2 + min(((u[i, 2] - u[i, 1]) / dy), 0)^2)) where {i in 2:NX - 1, j in 1:1}, [2:NX - 1, NY:NY] = sum[i, j] (sqrt(max(((u[i, NY] - u[i - 1, NY]) / dx), 0)^2 + min(((u[i + 1, NY] - u[i, NY]) / dx), 0)^2 + max(((u[i, NY] - u[i, NY - 1]) / dy), 0)^2 + min(((2 * gT - u[i, NY] - u[i, NY]) / dy), 0)^2)) where {i in 2:NX - 1, j in NY:NY}, [1:1, 1:1] = sum[i, j] (sqrt(max(((u[1, 1] - (2 * gL - u[1, 1])) / dx), 0)^2 + min(((u[2, 1] - u[1, 1]) / dx), 0)^2 + max(((u[1, 1] - (2 * gB - u[1, 1])) / dy), 0)^2 + min(((u[1, 2] - u[1, 1]) / dy), 0)^2)) where {i in 1:1, j in 1:1}, [1:1, NY:NY] = sum[i, j] (sqrt(max(((u[1, NY] - (2 * gL - u[1, NY])) / dx), 0)^2 + min(((u[2, NY] - u[1, NY]) / dx), 0)^2 + max(((u[1, NY] - u[1, NY - 1]) / dy), 0)^2 + min(((2 * gT - u[1, NY] - u[1, NY]) / dy), 0)^2)) where {i in 1:1, j in NY:NY}, [NX:NX, 1:1] = sum[i, j] (sqrt(max(((u[NX, 1] - u[NX - 1, 1]) / dx), 0)^2 + min(((2 * gR - u[NX, 1] - u[NX, 1]) / dx), 0)^2 + max(((u[NX, 1] - (2 * gB - u[NX, 1])) / dy), 0)^2 + min(((u[NX, 2] - u[NX, 1]) / dy), 0)^2)) where {i in NX:NX, j in 1:1}, [NX:NX, NY:NY] = sum[i, j] (sqrt(max(((u[NX, NY] - u[NX - 1, NY]) / dx), 0)^2 + min(((2 * gR - u[NX, NY] - u[NX, NY]) / dx), 0)^2 + max(((u[NX, NY] - u[NX, NY - 1]) / dy), 0)^2 + min(((2 * gT - u[NX, NY] - u[NX, NY]) / dy), 0)^2)) where {i in NX:NX, j in NY:NY})</div>

First-order Godunov (Rouy-Tourin / Osher-Sethian) upwind discretization of the gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) for the non-negative-speed eikonal term psi_t = -S|grad psi| (S >= 0) on cartesian_uniform_2d, with an inhomogeneous DIRICHLET (state-constraint) boundary u = g prescribed on all four rectangle walls -- the natural boundary condition for eikonal / level-set problems |grad u| = f with u given on the boundary. It is the Dirichlet-closed sibling of grids/cartesian_uniform_2d/rules/godunov_norm_D2_periodic and the 2-D companion of godunov_norm_D1_dirichlet_bc: the interior Godunov Hamiltonian and its per-axis entropy (max/min) upwind selection are UNCHANGED -- only the out-of-range neighbor on each boundary strip changes from a periodic wrap to a prescribed-value (mirror-ghost) source. COMPOUND MATCH (esm-spec 9.6.3, the canonical 2-D gradient-magnitude compound-precedence example): matches the whole compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) rooted at the outer sqrt at priority 10, IDENTICAL to godunov_norm_D2_periodic's match, firing before any plain first-derivative rule could lower an inner D. Match-scoped with where {u: {shape: [x, y]}} (esm-spec 9.6.1). THE SCHEME: the makearray tiles [x, y] as the same 3x3 block partition {1} u [2,NX-1] u {NX} x {1} u [2,NY-1] u {NY} as central_D2x/y and the periodic sibling. The doubly-interior block [2,NX-1] x [2,NY-1] is a single apply_expression_template of the imported stencils/godunov_norm_D2_interior aggregate. The four edge strips wrap their single out-of-range axis with a mirror ghost against the prescribed wall value: on cell-centered walls half a cell outside the first/last center, the ghost is u[0,j]=2*gL-u[1,j] (left), u[NX+1,j]=2*gR-u[NX,j] (right), u[i,0]=2*gB-u[i,1] (bottom), u[i,NY+1]=2*gT-u[i,NY] (top); the four corner cells wrap BOTH axes with the corresponding pair of ghosts. Every cell evaluates the identical per-axis max/min entropy norm sqrt(max(D-_x,0)^2 + min(D+_x,0)^2 + max(D-_y,0)^2 + min(D+_y,0)^2); the edge and corner faces are FULL-RANK aggregates (singleton output range on the wrapped axis) mirroring central_D2x_dirichlet_bc so every simulator fills the strips and corners without rank promotion. gL, gR are consumer-supplied free names giving the value on the left x-wall x=x0 and right x-wall x=x0+NX*dx (constant along that wall, all j); gB, gT give the value on the bottom y-wall y=y0 and top y-wall y=y0+NY*dy (constant along that wall, all i); a parameter not overridden defaults to 0 (esm-spec 6), so gL=gR=gB=gT=0 is the homogeneous case. The g* terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). EXACTNESS: the scheme ANNIHILATES a constant field exactly (every difference is 0 to the bit, the null-space property |grad(const)| = 0) and reproduces |grad u| exactly for any field whose wall traces are captured by the scalar wall values (e.g. u = c constant, or the manufactured g = c + sin(pi x) sin(pi y) whose four walls all equal c). A general non-constant affine field u = a x + b y + c has wall traces that VARY along each wall (the left wall trace c + b*y depends on y), which a scalar wall value cannot represent exactly; such array-valued wall data is out of scope for the constant-along-wall Dirichlet free-name convention shared with central_D2x/y_dirichlet_bc, so 2-D exactness is stated for constant / constant-wall-trace fields. First order O(dx + dy) in space on smooth data; the mirror ghost is 2nd-order-consistent so the boundary retains the interior order for Dirichlet-compatible solutions. dx, dy are free names (consumer-supplied geometry contract). A model imports exactly one gradient-magnitude rule for the [x, y] pair. IMPLEMENTATION NOTE: the interior Godunov Hamiltonian's sqrt(sum of squares) / min-max is non-smooth and the 2-D lowered RHS is large; the eikonal MMS driver problems/godunov_norm_2d_eikonal_dirichlet_mms.esm verifies this rule in Julia (Tsit5) and Python (LSODA); the Rust diffsol integrator lacks a step-size fail-fast guard and wedges on this RHS (see the case manifests' blocked_upstream_bindings, by analogy to the periodic sibling).

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. The max/min upwind (entropy) construction of \|grad phi\| for the level-set Hamilton-Jacobi term. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- Rouy, E., & Tourin, A. (1992). A viscosity solutions approach to shape-from-shading. SIAM Journal on Numerical Analysis, 29(3), 867-884. The Godunov Hamiltonian sqrt(max(D-,0)^2 + min(D+,0)^2) as the monotone, first-order-convergent discretization of the eikonal \|grad u\|. [`doi:10.1137/0729053`](https://doi.org/10.1137/0729053)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.12 (cell-centered Dirichlet via mirror ghost).

**Conformance fixtures**

- ast: [`godunov_norm_D2_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/godunov_norm_D2_dirichlet_bc) — status `active`
- convergence: [`godunov_norm_2d_eikonal_dirichlet_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_norm_2d_eikonal_dirichlet_mms) — status `active`
- simulation: [`godunov_norm_2d_eikonal_dirichlet_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/godunov_norm_2d_eikonal_dirichlet_mms) — status `active`

#### Convergence — case [`godunov_norm_2d_eikonal_dirichlet_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_norm_2d_eikonal_dirichlet_mms)

Error norms read from the committed golden [`tests/conformance/convergence/godunov_norm_2d_eikonal_dirichlet_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/godunov_norm_2d_eikonal_dirichlet_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.604e-02 | — | 6.653e-02 | — |
| 32 | 1.344e-02 | 0.95 | 3.502e-02 | 0.93 |
| 64 | 6.838e-03 | 0.98 | 1.819e-02 | 0.95 |
| 128 | 3.449e-03 | 0.99 | 9.366e-03 | 0.96 |

Expected order: **1** (± 0.2).

![Convergence of godunov_norm_2d_eikonal_dirichlet_mms (log-log)](/plots/convergence/godunov_norm_2d_eikonal_dirichlet_mms.png)

### `godunov_norm_D2_periodic`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/godunov_norm_D2_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/godunov_norm_D2_periodic.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 10 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 2:NY - 1] = sum[i, j] (sqrt(max(((u[i, j] - u[i - 1, j]) / dx), 0)^2 + min(((u[i + 1, j] - u[i, j]) / dx), 0)^2 + max(((u[i, j] - u[i, j - 1]) / dy), 0)^2 + min(((u[i, j + 1] - u[i, j]) / dy), 0)^2)) where {i in 2:NX - 1, j in 2:NY - 1}, [1:1, 2:NY - 1] = sum[i, j] (sqrt(max(((u[1, j] - u[NX, j]) / dx), 0)^2 + min(((u[2, j] - u[1, j]) / dx), 0)^2 + max(((u[1, j] - u[1, j - 1]) / dy), 0)^2 + min(((u[1, j + 1] - u[1, j]) / dy), 0)^2)) where {i in 1:1, j in 2:NY - 1}, [NX:NX, 2:NY - 1] = sum[i, j] (sqrt(max(((u[NX, j] - u[NX - 1, j]) / dx), 0)^2 + min(((u[1, j] - u[NX, j]) / dx), 0)^2 + max(((u[NX, j] - u[NX, j - 1]) / dy), 0)^2 + min(((u[NX, j + 1] - u[NX, j]) / dy), 0)^2)) where {i in NX:NX, j in 2:NY - 1}, [2:NX - 1, 1:1] = sum[i, j] (sqrt(max(((u[i, 1] - u[i - 1, 1]) / dx), 0)^2 + min(((u[i + 1, 1] - u[i, 1]) / dx), 0)^2 + max(((u[i, 1] - u[i, NY]) / dy), 0)^2 + min(((u[i, 2] - u[i, 1]) / dy), 0)^2)) where {i in 2:NX - 1, j in 1:1}, [2:NX - 1, NY:NY] = sum[i, j] (sqrt(max(((u[i, NY] - u[i - 1, NY]) / dx), 0)^2 + min(((u[i + 1, NY] - u[i, NY]) / dx), 0)^2 + max(((u[i, NY] - u[i, NY - 1]) / dy), 0)^2 + min(((u[i, 1] - u[i, NY]) / dy), 0)^2)) where {i in 2:NX - 1, j in NY:NY}, [1:1, 1:1] = sum[i, j] (sqrt(max(((u[1, 1] - u[NX, 1]) / dx), 0)^2 + min(((u[2, 1] - u[1, 1]) / dx), 0)^2 + max(((u[1, 1] - u[1, NY]) / dy), 0)^2 + min(((u[1, 2] - u[1, 1]) / dy), 0)^2)) where {i in 1:1, j in 1:1}, [1:1, NY:NY] = sum[i, j] (sqrt(max(((u[1, NY] - u[NX, NY]) / dx), 0)^2 + min(((u[2, NY] - u[1, NY]) / dx), 0)^2 + max(((u[1, NY] - u[1, NY - 1]) / dy), 0)^2 + min(((u[1, 1] - u[1, NY]) / dy), 0)^2)) where {i in 1:1, j in NY:NY}, [NX:NX, 1:1] = sum[i, j] (sqrt(max(((u[NX, 1] - u[NX - 1, 1]) / dx), 0)^2 + min(((u[1, 1] - u[NX, 1]) / dx), 0)^2 + max(((u[NX, 1] - u[NX, NY]) / dy), 0)^2 + min(((u[NX, 2] - u[NX, 1]) / dy), 0)^2)) where {i in NX:NX, j in 1:1}, [NX:NX, NY:NY] = sum[i, j] (sqrt(max(((u[NX, NY] - u[NX - 1, NY]) / dx), 0)^2 + min(((u[1, NY] - u[NX, NY]) / dx), 0)^2 + max(((u[NX, NY] - u[NX, NY - 1]) / dy), 0)^2 + min(((u[NX, 1] - u[NX, NY]) / dy), 0)^2)) where {i in NX:NX, j in NY:NY})</div>

First-order Godunov (Rouy-Tourin / Osher-Sethian) upwind discretization of the gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) for a non-negative-speed Hamilton-Jacobi term (psi_t = -S|grad psi|, S >= 0; the level-set / eikonal norm) on cartesian_uniform_2d, with a DOUBLY-PERIODIC boundary (both x and y wrapped as circles). This is the full 2-D scheme of archive/discretizations/finite_difference/godunov_norm_1st_uniform_cartesian_2d.json expressed in the D convention that grids/cartesian_uniform_2d already uses (central_D1x_*, central_D2x_*), and the 2-D companion of grids/cartesian_uniform_1d/rules/godunov_norm_D1_periodic (the 1-D reduction keeps only the x term; this rule sums the x AND y axes). COMPOUND MATCH (esm-spec 9.6.3 compound precedence): the rule matches the whole compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) rooted at the outer sqrt — cartesian's plain-D first-derivative idiom squared and summed under a square root, the exact structural analogue of the archive's applies_to sqrt(grad(u,x)^2 + grad(u,y)^2). This 2-D gradient magnitude is precisely the canonical compound-precedence example esm-spec 9.6.3 (and 9.6.8) names. It carries priority 10 (>= the D1/D2 rules' 10) so that, under the section-9.6.3 outermost-first fixpoint, it fires on the entire sqrt(D_x^2 + D_y^2) term BEFORE any plain first-derivative rule (central_D1x_periodic / central_D1y_periodic) could lower the inner D(u, wrt: x) or D(u, wrt: y); the whole compound is replaced by this one makearray in a single step and the engine does not descend into the fresh body. Match-scoped with `where: {u: {shape: [x, y]}}` (esm-spec 9.6.1) so it fires only on the gradient magnitude of a bare field declared over this grid's [x, y] index sets; under import-edge renaming (esm-spec 9.7.7) the two `wrt` literals and the where shape follow the renamed axes together, so two cartesian instances coexist scoped to their own grids. Because the constraint requires a bare shaped variable, a consumer forming the gradient magnitude of a COMPOUND inner expression must first bind it to a declared [x,y]-shaped observed. THE SCHEME: with per-axis backward differences D-_x = (u[i,j] - u[i-1,j]) / dx, D-_y = (u[i,j] - u[i,j-1]) / dy and forward differences D+_x = (u[i+1,j] - u[i,j]) / dx, D+_y = (u[i,j+1] - u[i,j]) / dy, the discrete norm is sqrt( max(D-_x,0)^2 + min(D+_x,0)^2 + max(D-_y,0)^2 + min(D+_y,0)^2 ). The per-axis max/min is the entropy (upwind) fix of Osher & Sethian (1988) and Rouy & Tourin (1992): for the non-negative-speed term it selects the one-sided difference along the incoming characteristic per axis per point (a centered difference is unstable for this hyperbolic operator). It is EXACT on a linear field: for u = a x + b y + c the interior gives D-_x = D+_x = a and D-_y = D+_y = b, so sqrt(max(a,0)^2 + min(a,0)^2 + max(b,0)^2 + min(b,0)^2) = sqrt(a^2 + b^2) = |grad u| to rounding, reproducing the continuous gradient magnitude with no truncation error (the low-order analogue of the null-space property that annihilates constants, a = b = 0 -> 0). First order in space, O(dx + dy), on smooth data; it retains first-order accuracy uniformly through the |grad| kinks (points where a gradient component changes sign, or where grad u vanishes) that the entropy fix is designed to resolve. The makearray tiles [x, y] as a 3x3 block partition {1}u[2,NX-1]u{NX} x {1}u[2,NY-1]u{NY} — a 9-region closure identical in shape to grids/cartesian_uniform_2d/rules/central_Dxy_periodic, because the Godunov 5-point stencil reaches i+-1 AND j+-1: the doubly-interior block [2,NX-1]x[2,NY-1] holds the Godunov-norm aggregate over the full 5-point stencil (u[i,j], u[i-1,j], u[i+1,j], u[i,j-1], u[i,j+1]); the four edge strips wrap their single out-of-range axis (i=1 uses the ghost column NX for D-_x, i=NX uses column 1 for D+_x; j=1 uses row NY for D-_y, j=NY uses row 1 for D+_y) with the in-range axis interior; the four corner cells wrap BOTH axes. Every cell evaluates the identical max/min entropy norm; the periodic wrap is exact on the torus, so the operator is self-contained and testable on a closed 2-torus. The edge and corner faces are FULL-RANK aggregates (singleton output range on the wrapped axis) mirroring central_Dxy_periodic so every binding's simulator fills the strips and corner cells without rank promotion. A doubly-periodic linear ramp is a contradiction (a non-zero ramp cannot wrap), so the wrap faces carry the seam of a ramp rather than sqrt(a^2+b^2); the exactness property is the interior statement, and the periodic gate that ships with this rule uses the constant/interior form (problems/godunov_norm_2d_periodic_exact.esm). `dx`/`dy` are free names resolving to the consuming model's spacing variables (grids/cartesian_uniform_2d/grid.esm consumer-supplied geometry contract); the periodic wrap is dx/dy-scale-invariant, so the rule serves any x/y extent unchanged. A model imports exactly one gradient-magnitude rule for the [x, y] pair.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. The max/min upwind (entropy) construction of \|grad phi\| for the level-set Hamilton-Jacobi term. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- Rouy, E., & Tourin, A. (1992). A viscosity solutions approach to shape-from-shading. SIAM Journal on Numerical Analysis, 29(3), 867-884. The Godunov Hamiltonian sqrt(sum_d max(D-_d,0)^2 + min(D+_d,0)^2) as the monotone, first-order-convergent discretization of the eikonal \|grad u\|. [`doi:10.1137/0729053`](https://doi.org/10.1137/0729053)

**Conformance fixtures**

- ast: [`godunov_norm_D2_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/godunov_norm_D2_periodic) — status `active`
- convergence: [`godunov_norm_2d_eikonal_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_norm_2d_eikonal_mms) — status `active`
- simulation: [`godunov_norm_2d_eikonal_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/godunov_norm_2d_eikonal_mms) — status `active`
- simulation: [`godunov_norm_2d_periodic_exact`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/godunov_norm_2d_periodic_exact) — status `active`

#### Convergence — case [`godunov_norm_2d_eikonal_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_norm_2d_eikonal_mms)

Error norms read from the committed golden [`tests/conformance/convergence/godunov_norm_2d_eikonal_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/godunov_norm_2d_eikonal_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.558e-01 | — | 3.467e-01 | — |
| 32 | 7.807e-02 | 1.00 | 1.742e-01 | 0.99 |
| 64 | 3.940e-02 | 0.99 | 8.720e-02 | 1.00 |
| 128 | 1.986e-02 | 0.99 | 4.361e-02 | 1.00 |

Expected order: **1** (± 0.15).

![Convergence of godunov_norm_2d_eikonal_mms (log-log)](/plots/convergence/godunov_norm_2d_eikonal_mms.png)

### `godunov_upwind_D1x_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/godunov_upwind_D1x_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/godunov_upwind_D1x_extrapolate_bc.esm)

**Rewrites** `D(u)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = godunov_upwind_D1x_interior&lt;u=u&gt;, [1:1, 1:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[1, j] - extrapolate_linear_ghost&lt;a=u[1, j], b=u[2, j]&gt;) / dx, dr=(u[1 + 1, j] - u[1, j]) / dx&gt;) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[NX, j] - u[NX - 1, j]) / dx, dr=(extrapolate_linear_ghost&lt;a=u[NX, j], b=u[NX - 1, j]&gt; - u[NX, j]) / dx&gt;) where {i in NX:NX, j in y})</div>

Per-axis quantity diff2x of NCAR fire_behavior's classical Godunov level-set option (fire_upwinding = 4) on cartesian_uniform_2d with the linear-extrapolation halo closure, from which fire_behavior builds the front normal nvx = diff2x/sqrt(|grad psi|^2 + eps) while taking |grad psi| itself from the Rouy-Tourin Hamiltonian (rules/godunov_norm_D2_*). Match D(u, wrt: x) at priority 10, scoped with where {u: {shape: [x, y]}}. THE SCHEME, as Calc_tend_ls computes it: within bdy_eno1_cells of ANY wall (either axis) the ENO1 minmod of the two one-sided differences (fire_behavior's BDY_ENO1 band, whose test precedes the `select case (fire_upwinding)` and so overrides case (4)), and the Godunov-upwind quantity max(D-_x u, 0) - min(D+_x u, 0) everywhere else. Three regions -- the interior [2, NX-1] x [1, NY] is a single apply of stencils/godunov_upwind_D1x_interior (the per-cell wall-band switch, stencils/godunov_eno1_D1); the two wall lines i=1 and i=NX lie in the band for every bdy_eno1_cells >= 1 and are ENO1-only (stencils/eno1_minmod), with the out-of-range neighbour replaced by the guarded linear-extrapolation ghost. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; instead each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j], u[i,0], u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the three-point ENO1 stencil is evaluated at a wall cell, so one ghost suffices, as in the reference code. SIGN CAVEAT (fire_behavior note N44): for a field decreasing along x the Godunov branch returns -D+ > 0, a magnitude, while the ENO1 branch in the wall band returns the signed (negative) derivative -- the reference code's own inconsistency, reproduced; a consumer building a signed normal must carry the sign separately. EXACTNESS: on a linear field with positive slope along x the result is the slope on every cell EXCEPT the i=1 wall, where the field decreases toward the wall, the guard clamps the ghost to the neighbour (u[0] = u[2]), and the ENO1 minmod of D- = -slope, D+ = +slope is 0; first order otherwise. FREE NAMES (consumer-supplied, resolved in the consuming model's scope like the grid geometry dx): bdy_eno1_cells, the ENO1 wall band in cells (>= 1; fire_behavior's BDY_ENO1 = 10). Requires NX >= 3.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`godunov_upwind_D1x_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/godunov_upwind_D1x_extrapolate_bc) — status `active`
- convergence: [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms) — status `active`
- simulation: [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/godunov_upwind_D1_2d_paraboloid_mms) — status `active`

#### Convergence — case [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms)

Error norms read from the committed golden [`tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.246e-01 | — | 1.768e-02 | — |
| 32 | 6.215e-02 | 1.00 | 8.839e-03 | 1.00 |
| 64 | 3.114e-02 | 1.00 | 4.419e-03 | 1.00 |
| 128 | 1.560e-02 | 1.00 | 2.210e-03 | 1.00 |

Expected order: **1** (± 0.25).

![Convergence of godunov_upwind_D1_2d_paraboloid_mms (log-log)](/plots/convergence/godunov_upwind_D1_2d_paraboloid_mms.png)

### `godunov_upwind_D1y_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/godunov_upwind_D1y_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/godunov_upwind_D1y_extrapolate_bc.esm)

**Rewrites** `D(u)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = godunov_upwind_D1y_interior&lt;u=u&gt;, [1:NX, 1:1] = sum[i, j] (eno1_minmod&lt;dl=(u[i, 1] - extrapolate_linear_ghost&lt;a=u[i, 1], b=u[i, 2]&gt;) / dy, dr=(u[i, 1 + 1] - u[i, 1]) / dy&gt;) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[i, NY] - u[i, NY - 1]) / dy, dr=(extrapolate_linear_ghost&lt;a=u[i, NY], b=u[i, NY - 1]&gt; - u[i, NY]) / dy&gt;) where {i in x, j in NY:NY})</div>

Per-axis quantity diff2y of NCAR fire_behavior's classical Godunov level-set option (fire_upwinding = 4) on cartesian_uniform_2d with the linear-extrapolation halo closure, from which fire_behavior builds the front normal nvy = diff2y/sqrt(|grad psi|^2 + eps) while taking |grad psi| itself from the Rouy-Tourin Hamiltonian (rules/godunov_norm_D2_*). Match D(u, wrt: y) at priority 10, scoped with where {u: {shape: [x, y]}}. THE SCHEME, as Calc_tend_ls computes it: within bdy_eno1_cells of ANY wall (either axis) the ENO1 minmod of the two one-sided differences (fire_behavior's BDY_ENO1 band, whose test precedes the `select case (fire_upwinding)` and so overrides case (4)), and the Godunov-upwind quantity max(D-_y u, 0) - min(D+_y u, 0) everywhere else. Three regions -- the interior [1, NX] x [2, NY-1] is a single apply of stencils/godunov_upwind_D1y_interior (the per-cell wall-band switch, stencils/godunov_eno1_D1); the two wall lines j=1 and j=NY lie in the band for every bdy_eno1_cells >= 1 and are ENO1-only (stencils/eno1_minmod), with the out-of-range neighbour replaced by the guarded linear-extrapolation ghost. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; instead each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j], u[i,0], u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the three-point ENO1 stencil is evaluated at a wall cell, so one ghost suffices, as in the reference code. SIGN CAVEAT (fire_behavior note N44): for a field decreasing along y the Godunov branch returns -D+ > 0, a magnitude, while the ENO1 branch in the wall band returns the signed (negative) derivative -- the reference code's own inconsistency, reproduced; a consumer building a signed normal must carry the sign separately. EXACTNESS: on a linear field with positive slope along y the result is the slope on every cell EXCEPT the j=1 wall, where the field decreases toward the wall, the guard clamps the ghost to the neighbour (u[0] = u[2]), and the ENO1 minmod of D- = -slope, D+ = +slope is 0; first order otherwise. FREE NAMES (consumer-supplied, resolved in the consuming model's scope like the grid geometry dy): bdy_eno1_cells, the ENO1 wall band in cells (>= 1; fire_behavior's BDY_ENO1 = 10). Requires NY >= 3.

**References**

- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`godunov_upwind_D1y_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/godunov_upwind_D1y_extrapolate_bc) — status `active`
- convergence: [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms) — status `active`
- simulation: [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/godunov_upwind_D1_2d_paraboloid_mms) — status `active`

#### Convergence — case [`godunov_upwind_D1_2d_paraboloid_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms)

Error norms read from the committed golden [`tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/godunov_upwind_D1_2d_paraboloid_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.246e-01 | — | 1.768e-02 | — |
| 32 | 6.215e-02 | 1.00 | 8.839e-03 | 1.00 |
| 64 | 3.114e-02 | 1.00 | 4.419e-03 | 1.00 |
| 128 | 1.560e-02 | 1.00 | 2.210e-03 | 1.00 |

Expected order: **1** (± 0.25).

![Convergence of godunov_upwind_D1_2d_paraboloid_mms (log-log)](/plots/convergence/godunov_upwind_D1_2d_paraboloid_mms.png)

### `varcoeff_laplacian_x_dirichlet_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_dirichlet_bc.esm)

**Rewrites** `D(k * D(u)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "x"}]}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = varcoeff_laplacian_x_interior&lt;k=k, u=u&gt;, [1:1, 1:NY] = sum[i, j] ((0.5 * (k[1, j] + k[2, j]) * (u[2, j] - u[1, j]) - 2 * k[1, j] * (u[1, j] - gL)) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((0.5 * (k[NX - 1, j] + k[NX, j]) * (u[NX - 1, j] - u[NX, j]) + 2 * k[NX, j] * (gR - u[NX, j])) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two x-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:x)]}, wrt:x ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative. This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_dirichlet_bc: the operated axis x gets an interior band plus two mirror-ghost x-faces, while the y axis is full-rank (all rows j treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:x). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k (identical to varcoeff_laplacian_x_zero_flux_bc): k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [2,NX-1]x[1,NY] aggregate and a shared loop-index name is captured. Face conductivity is the ARITHMETIC MEAN kf[i+1/2,j]=(k[i,j]+k[i+1,j])/2 (O(dx^2) at the face for a differentiable coefficient, affine in k so it composes with the nonlinear dual). Interior i in [2,NX-1] (all j): the imported conservative x-flux-difference stencil. Faces i=1 and i=NX: mirror-ghost elimination against a prescribed wall value combined with the boundary-face conductivity. On this cell-centered grid the wall sits half a cell outside the first/last x-center, so linear reconstruction to the wall gives u[0,j]=2*gL-u[1,j] (and u[NX+1,j]=2*gR-u[NX,j]); the left x-face flux F[1/2,j]=kf[1/2,j]*(u[1,j]-u[0,j])/dx uses the wall conductivity k[1,j] (the arithmetic-mean face value degenerates to the first cell value under the even ghost k[0,j]=k[1,j]), so F[1/2,j]=k[1,j]*2*(u[1,j]-gL)/dx and Lx u[1,j] = ( (k[1,j]+k[2,j])/2 (u[2,j]-u[1,j]) - 2 k[1,j] (u[1,j]-gL) ) / dx^2; symmetrically at the right x-wall F[NX+1/2,j]=k[NX,j]*2*(gR-u[NX,j])/dx and Lx u[NX,j] = ( (k[NX-1,j]+k[NX,j])/2 (u[NX-1,j]-u[NX,j]) + 2 k[NX,j] (gR-u[NX,j]) ) / dx^2. `gL` and `gR` are consumer-supplied free names (the same free-name contract as the grid geometry `dx`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the left x-wall x=x0 and the right x-wall x=x0+NX*dx, constant along that wall (all j); a consuming model declares them as ordinary real parameters, and a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so gL=gR=0 is the homogeneous-Dirichlet x-wall. The gL/gR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_dirichlet_bc. `dx` is the grid's consumer-supplied x-spacing free name. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the first term. Reduces to central_D2x_dirichlet_bc when k is spatially constant (k=1 gives the mirror-ghost residuals (u[2,j]-3u[1,j]+2gL)/dx^2 and (u[NX-1,j]-3u[NX,j]+2gR)/dx^2). Interior O(dx^2); the one-sided mirror-ghost faces are O(dx) locally, preserving global 2nd-order L2 convergence for Dirichlet-compatible solutions. A model imports exactly one second-derivative/Laplacian rule for x.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (cell-centered Dirichlet via mirror ghost).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_x_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_x_dirichlet_bc) — status `active`
- convergence: [`heat_2d_varcoeff_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_dirichlet) — status `active`
- simulation: [`heat_2d_varcoeff_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_x_dirichlet) — status `active`

#### Convergence — case [`heat_2d_varcoeff_x_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_x_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_x_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.366e-04 | — | 2.108e-04 | — |
| 32 | 5.687e-05 | 2.06 | 5.193e-05 | 2.02 |
| 64 | 1.408e-05 | 2.01 | 1.290e-05 | 2.01 |
| 128 | 3.513e-06 | 2.00 | 3.218e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_x_dirichlet (log-log)](/plots/convergence/heat_2d_varcoeff_x_dirichlet.png)

### `varcoeff_laplacian_x_neumann_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_neumann_bc.esm)

**Rewrites** `D(k * D(u)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "x"}]}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = varcoeff_laplacian_x_interior&lt;k=k, u=u&gt;, [1:1, 1:NY] = sum[i, j] (0.5 * (k[1, j] + k[2, j]) * (u[2, j] - u[1, j]) / (dx * dx) - qXL / dx) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] (0.5 * (k[NX - 1, j] + k[NX, j]) * (u[NX - 1, j] - u[NX, j]) / (dx * dx) + qXR / dx) where {i in NX:NX, j in y})</div>

The x-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with inhomogeneous NEUMANN (fixed-flux) conditions on the two x-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:x)]}, wrt:x ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative. This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_neumann_bc: the operated axis x gets an interior band plus two prescribed-flux x-faces, while the y axis is full-rank (all rows j treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:x). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k (identical to varcoeff_laplacian_x_zero_flux_bc): k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [2,NX-1]x[1,NY] aggregate and a shared loop-index name is captured. Face conductivity is the ARITHMETIC MEAN kf[i+1/2,j]=(k[i,j]+k[i+1,j])/2. Interior i in [2,NX-1] (all j): the imported conservative x-flux-difference stencil. Faces i=1 and i=NX: the domain-boundary x-face flux is PRESCRIBED DIRECTLY -- the natural Neumann treatment for the finite-volume form, requiring no ghost geometry and no boundary conductivity. `qXL` and `qXR` are consumer-supplied free names (the same free-name contract as the grid geometry `dx`; see grids/cartesian_uniform_2d/grid.esm) giving the PHYSICAL FLUX (k du/dx, in the +x direction) at the left x-wall x=x0 and the right x-wall x=x0+NX*dx respectively, constant along that wall (all j). The left cell balance Lx u[1,j] = (F[3/2,j]-F[1/2,j])/dx with F[1/2,j]=qXL gives Lx u[1,j] = (k[1,j]+k[2,j])/2 (u[2,j]-u[1,j]) / dx^2 - qXL/dx; the right cell balance Lx u[NX,j] = (F[NX+1/2,j]-F[NX-1/2,j])/dx with F[NX+1/2,j]=qXR gives Lx u[NX,j] = (k[NX-1,j]+k[NX,j])/2 (u[NX-1,j]-u[NX,j]) / dx^2 + qXR/dx. A consuming model declares qXL/qXR as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so the homogeneous case is the zero default and this rule REDUCES EXACTLY to varcoeff_laplacian_x_zero_flux_bc when qXL=qXR=0 -- the prescribed-flux terms qXL/dx and qXR/dx vanish and the two face values become byte-identical to that rule's dropped-flux faces (kf[3/2,j\](u[2,j]-u[1,j])/dx^2 at i=1, kf[NX-1/2,j\](u[NX-1,j]-u[NX,j])/dx^2 at i=NX), the flux at each x-wall dropped to zero. The qXL/qXR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_neumann_bc. `dx` is the grid's consumer-supplied x-spacing free name. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the first term. Interior O(dx^2); the prescribed-flux faces impose the exact boundary flux, preserving global 2nd-order L2 convergence for Neumann-compatible solutions. A model imports exactly one second-derivative/Laplacian rule for x.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (Neumann boundary conditions by prescribing the domain-face flux).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_x_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_x_neumann_bc) — status `active`
- convergence: [`heat_2d_varcoeff_x_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_neumann) — status `active`
- simulation: [`heat_2d_varcoeff_x_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_x_neumann) — status `active`

#### Convergence — case [`heat_2d_varcoeff_x_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_x_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_x_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.293e-03 | — | 8.851e-04 | — |
| 32 | 3.147e-04 | 2.04 | 2.165e-04 | 2.03 |
| 64 | 7.814e-05 | 2.01 | 5.354e-05 | 2.02 |
| 128 | 1.950e-05 | 2.00 | 1.338e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_x_neumann (log-log)](/plots/convergence/heat_2d_varcoeff_x_neumann.png)

### `varcoeff_laplacian_x_periodic`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_periodic.esm)

**Rewrites** `D(k * D(u)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "x"}]}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = varcoeff_laplacian_x_interior&lt;k=k, u=u&gt;, [1:1, 1:NY] = sum[i, j] ((0.5 * (k[1, j] + k[2, j]) * (u[2, j] - u[1, j]) - 0.5 * (k[NX, j] + k[1, j]) * (u[1, j] - u[NX, j])) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] ((0.5 * (k[NX, j] + k[1, j]) * (u[1, j] - u[NX, j]) - 0.5 * (k[NX - 1, j] + k[NX, j]) * (u[NX, j] - u[NX - 1, j])) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with a PERIODIC x-boundary, over the full 2-D field u[x,y]. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:x)]}, wrt:x ) -- outer-D of (k times the inner x-gradient) -- at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative (exactly as central_D2x_periodic fires on D(D(u,x),x)). This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_periodic: the operated axis x gets an interior band plus two wrapped x-faces, while the y axis is full-rank (all rows j treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state, inner D's argument) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write the coefficient on the LEFT: k * D(u, wrt:x). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on the xx-conservative-compound of bare fields declared over this 2-D grid's [x,y] axes (distinct from the 1-D shape-[x] varcoeff rules); under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together, so two instances coexist scoped to their own meshes. ONE rule serves two duals: the VARIABLE-COEFFICIENT case (k=k(x,y), a given periodic field) and the NONLINEAR case (k declared as an observed k=f(u)); the body only reads k by index at the faces, so a state-dependent k just makes the expanded stencil nonlinear in u with no reformulation (the keyed-factor field contract, like the grid's dx). AUTHORING CONTRACT for k: it is a declared [x,y]-shaped field; when supplied as an aggregate observed it MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because the coefficient observed is inlined into the stencil's [2,NX-1]x[1,NY] aggregate at evaluation and a shared loop-index name is captured (silently wrong interior values). Face conductivity is the ARITHMETIC MEAN kf[i+1/2,j]=(k[i,j]+k[i+1,j])/2 (O(dx^2) at the face for a differentiable coefficient, and affine in k so it composes with the nonlinear dual). Interior i in [2,NX-1] (all j): the imported conservative x-flux-difference stencil. The two x-face columns i=1 and i=NX carry the same conservative x-flux difference with the out-of-range x-neighbor AND its coefficient wrapped periodically: at i=1 the left x-face flux uses the ghost u[0,j]=u[NX,j] and k[0,j]=k[NX,j], giving kf[1/2,j]=(k[NX,j]+k[1,j])/2 and Lx u[1,j] = ( (k[1,j]+k[2,j])/2 (u[2,j]-u[1,j]) - (k[NX,j]+k[1,j])/2 (u[1,j]-u[NX,j]) ) / dx^2; at i=NX the right x-face flux uses u[NX+1,j]=u[1,j] and k[NX+1,j]=k[1,j], giving kf[NX+1/2,j]=(k[NX,j]+k[1,j])/2 and Lx u[NX,j] = ( (k[NX,j]+k[1,j])/2 (u[1,j]-u[NX,j]) - (k[NX-1,j]+k[NX,j])/2 (u[NX,j]-u[NX-1,j]) ) / dx^2. The SAME wrapped x-face flux (k[NX,j]+k[1,j])/2 (u[1,j]-u[NX,j])/dx is the left face of column 1 and the right face of column NX, so the discrete x-operator is exactly conservative on each row (the circulant variable-coefficient x-Laplacian) and its row sums vanish. The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_neumann_bc so every binding's simulator fills the (1 x NY) columns without rank promotion. `dx` is the grid's consumer-supplied x-spacing free name; the periodic wrap is dx-scale-invariant, so the rule serves any x-extent unchanged. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the first term. Reduces to central_D2x_periodic when k is spatially constant. Interior and periodic faces are all O(dx^2), so the whole x-operator is 2nd order with no boundary reduction. A model imports exactly one second-derivative/Laplacian rule for x.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_x_periodic) — status `active`
- convergence: [`heat_2d_varcoeff_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_periodic) — status `active`
- simulation: [`heat_2d_varcoeff_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_x_periodic) — status `active`

#### Convergence — case [`heat_2d_varcoeff_x_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_x_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_x_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 4.046e-02 | — | 8.114e-04 | — |
| 32 | 9.975e-03 | 2.02 | 2.028e-04 | 2.00 |
| 64 | 2.485e-03 | 2.01 | 5.069e-05 | 2.00 |
| 128 | 6.207e-04 | 2.00 | 1.267e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_x_periodic (log-log)](/plots/convergence/heat_2d_varcoeff_x_periodic.png)

### `varcoeff_laplacian_x_zero_flux_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:zero_flux</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_zero_flux_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_x_zero_flux_bc.esm)

**Rewrites** `D(k * D(u)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "x"}]}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NX - 1, 1:NY] = varcoeff_laplacian_x_interior&lt;k=k, u=u&gt;, [1:1, 1:NY] = sum[i, j] (0.5 * (k[1, j] + k[2, j]) * (u[2, j] - u[1, j]) / (dx * dx)) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] (0.5 * (k[NX - 1, j] + k[NX, j]) * (u[NX - 1, j] - u[NX, j]) / (dx * dx)) where {i in NX:NX, j in y})</div>

The x-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with zero-flux (homogeneous Neumann) conditions on the two x-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:x)]}, wrt:x ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative (exactly as central_D2x fires on D(D(u,x),x)). This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_zero_flux_bc: the operated axis x gets an interior band plus two dropped-flux x-faces, while the y axis is full-rank (all rows j treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:x). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k: k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [2,NX-1]x[1,NY] aggregate and a shared loop-index name is captured (silently wrong interior values). Face conductivity is the ARITHMETIC MEAN kf[i+1/2,j]=(k[i,j]+k[i+1,j])/2 (O(dx^2) at the face for a differentiable coefficient, affine in k so it composes with the nonlinear dual). Interior i in [2,NX-1] (all j): the imported conservative x-flux-difference stencil. The x-boundary faces impose zero flux exactly (natural for the finite-volume form -- no ghost geometry): at i=1 the left x-domain-face flux is dropped, giving Lx u[1,j] = kf[3/2,j\](u[2,j]-u[1,j])/dx^2 = (k[1,j]+k[2,j])/2 (u[2,j]-u[1,j])/dx^2; at i=NX the right x-domain-face flux is dropped, giving Lx u[NX,j] = -kf[NX-1/2,j\](u[NX,j]-u[NX-1,j])/dx^2 = (k[NX-1,j]+k[NX,j])/2 (u[NX-1,j]-u[NX,j])/dx^2. This is exactly varcoeff_laplacian_x_neumann_bc at qXL=qXR=0. The face values are FULL-RANK aggregates (singleton output range on the x axis, all j) mirroring central_D2x_neumann_bc. `dx` is the grid's consumer-supplied x-spacing free name; the zero-flux closure is dx-scale-invariant, so the rule serves any x-extent unchanged. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the first term. Reduces to the constant-coefficient zero-gradient x-half when k is spatially constant. Interior O(dx^2); the one-sided zero-flux faces are O(dx) locally, which preserves global 2nd-order L2 convergence for zero-flux-compatible solutions (flux k u_x vanishing at both x-walls). A model imports exactly one second-derivative/Laplacian rule for x.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (Neumann boundary conditions by dropping the domain-face flux).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_x_zero_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_x_zero_flux_bc) — status `active`
- convergence: [`heat_2d_varcoeff_x_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_zero_flux) — status `active`
- simulation: [`heat_2d_varcoeff_x_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_x_zero_flux) — status `active`

#### Convergence — case [`heat_2d_varcoeff_x_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_x_zero_flux)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_x_zero_flux/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_x_zero_flux/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 6.726e-03 | — | 2.584e-03 | — |
| 32 | 1.692e-03 | 1.99 | 6.522e-04 | 1.99 |
| 64 | 4.236e-04 | 2.00 | 1.630e-04 | 2.00 |
| 128 | 1.060e-04 | 2.00 | 4.078e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_x_zero_flux (log-log)](/plots/convergence/heat_2d_varcoeff_x_zero_flux.png)

### `varcoeff_laplacian_y_dirichlet_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_dirichlet_bc.esm)

**Rewrites** `D(k * D(u)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "y"}]}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = varcoeff_laplacian_y_interior&lt;k=k, u=u&gt;, [1:NX, 1:1] = sum[i, j] ((0.5 * (k[i, 1] + k[i, 2]) * (u[i, 2] - u[i, 1]) - 2 * k[i, 1] * (u[i, 1] - gB)) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((0.5 * (k[i, NY - 1] + k[i, NY]) * (u[i, NY - 1] - u[i, NY]) + 2 * k[i, NY] * (gT - u[i, NY])) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with inhomogeneous DIRICHLET (fixed-value) conditions on the two y-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:y)]}, wrt:y ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative. This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_dirichlet_bc: the operated axis y gets an interior band plus two mirror-ghost y-faces, while the x axis is full-rank (all columns i treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:y). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k (identical to varcoeff_laplacian_y_zero_flux_bc): k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [1,NX]x[2,NY-1] aggregate and a shared loop-index name is captured. Face conductivity is the ARITHMETIC MEAN kf[i,j+1/2]=(k[i,j]+k[i,j+1])/2 (O(dy^2) at the face for a differentiable coefficient, affine in k so it composes with the nonlinear dual). Interior j in [2,NY-1] (all i): the imported conservative y-flux-difference stencil. Faces j=1 and j=NY: mirror-ghost elimination against a prescribed wall value combined with the boundary-face conductivity. On this cell-centered grid the wall sits half a cell outside the first/last y-center, so linear reconstruction to the wall gives u[i,0]=2*gB-u[i,1] (and u[i,NY+1]=2*gT-u[i,NY]); the bottom y-face flux F[i,1/2]=kf[i,1/2]*(u[i,1]-u[i,0])/dy uses the wall conductivity k[i,1] (the arithmetic-mean face value degenerates to the first cell value under the even ghost k[i,0]=k[i,1]), so F[i,1/2]=k[i,1]*2*(u[i,1]-gB)/dy and Ly u[i,1] = ( (k[i,1]+k[i,2])/2 (u[i,2]-u[i,1]) - 2 k[i,1] (u[i,1]-gB) ) / dy^2; symmetrically at the top y-wall F[i,NY+1/2]=k[i,NY]*2*(gT-u[i,NY])/dy and Ly u[i,NY] = ( (k[i,NY-1]+k[i,NY])/2 (u[i,NY-1]-u[i,NY]) + 2 k[i,NY] (gT-u[i,NY]) ) / dy^2. `gB` and `gT` are consumer-supplied free names (the same free-name contract as the grid geometry `dy`; see grids/cartesian_uniform_2d/grid.esm) giving the field value on the BOTTOM y-wall y=y0 (the j=1 side) and the TOP y-wall y=y0+NY*dy (the j=NY side), constant along that wall (all i). The names deliberately read gB/gT (bottom/top) so the y-wall assignment is unambiguous, the y-analogue of gL/gR on the x-walls. A consuming model declares them as ordinary real parameters, and a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so gB=gT=0 is the homogeneous-Dirichlet y-wall. The gB/gT terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_dirichlet_bc. `dy` is the grid's consumer-supplied y-spacing free name. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the second term. Reduces to central_D2y_dirichlet_bc when k is spatially constant (k=1 gives the mirror-ghost residuals (u[i,2]-3u[i,1]+2gB)/dy^2 and (u[i,NY-1]-3u[i,NY]+2gT)/dy^2). Interior O(dy^2); the one-sided mirror-ghost faces are O(dy) locally, preserving global 2nd-order L2 convergence for Dirichlet-compatible solutions. A model imports exactly one second-derivative/Laplacian rule for y.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (cell-centered Dirichlet via mirror ghost).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_y_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_y_dirichlet_bc) — status `active`
- convergence: [`heat_2d_varcoeff_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_dirichlet) — status `active`
- simulation: [`heat_2d_varcoeff_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_y_dirichlet) — status `active`

#### Convergence — case [`heat_2d_varcoeff_y_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_y_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_y_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 2.366e-04 | — | 2.108e-04 | — |
| 32 | 5.687e-05 | 2.06 | 5.193e-05 | 2.02 |
| 64 | 1.408e-05 | 2.01 | 1.290e-05 | 2.01 |
| 128 | 3.513e-06 | 2.00 | 3.218e-06 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_y_dirichlet (log-log)](/plots/convergence/heat_2d_varcoeff_y_dirichlet.png)

### `varcoeff_laplacian_y_neumann_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_neumann_bc.esm)

**Rewrites** `D(k * D(u)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "y"}]}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = varcoeff_laplacian_y_interior&lt;k=k, u=u&gt;, [1:NX, 1:1] = sum[i, j] (0.5 * (k[i, 1] + k[i, 2]) * (u[i, 2] - u[i, 1]) / (dy * dy) - qYL / dy) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] (0.5 * (k[i, NY - 1] + k[i, NY]) * (u[i, NY - 1] - u[i, NY]) / (dy * dy) + qYR / dy) where {i in x, j in NY:NY})</div>

The y-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with inhomogeneous NEUMANN (fixed-flux) conditions on the two y-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:y)]}, wrt:y ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative. This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_neumann_bc: the operated axis y gets an interior band plus two prescribed-flux y-faces, while the x axis is full-rank (all columns i treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:y). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k (identical to varcoeff_laplacian_y_zero_flux_bc): k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [1,NX]x[2,NY-1] aggregate and a shared loop-index name is captured. Face conductivity is the ARITHMETIC MEAN kf[i,j+1/2]=(k[i,j]+k[i,j+1])/2. Interior j in [2,NY-1] (all i): the imported conservative y-flux-difference stencil. Faces j=1 and j=NY: the domain-boundary y-face flux is PRESCRIBED DIRECTLY -- the natural Neumann treatment for the finite-volume form, requiring no ghost geometry and no boundary conductivity. `qYL` and `qYR` are consumer-supplied free names (the same free-name contract as the grid geometry `dy`; see grids/cartesian_uniform_2d/grid.esm) giving the PHYSICAL FLUX (k du/dy, in the +y direction) at the bottom y-wall y=y0 and the top y-wall y=y0+NY*dy respectively, constant along that wall (all i). The bottom cell balance Ly u[i,1] = (F[i,3/2]-F[i,1/2])/dy with F[i,1/2]=qYL gives Ly u[i,1] = (k[i,1]+k[i,2])/2 (u[i,2]-u[i,1]) / dy^2 - qYL/dy; the top cell balance Ly u[i,NY] = (F[i,NY+1/2]-F[i,NY-1/2])/dy with F[i,NY+1/2]=qYR gives Ly u[i,NY] = (k[i,NY-1]+k[i,NY])/2 (u[i,NY-1]-u[i,NY]) / dy^2 + qYR/dy. A consuming model declares qYL/qYR as ordinary real parameters; a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so the homogeneous case is the zero default and this rule REDUCES EXACTLY to varcoeff_laplacian_y_zero_flux_bc when qYL=qYR=0 -- the prescribed-flux terms qYL/dy and qYR/dy vanish and the two face values become byte-identical to that rule's dropped-flux faces (kf[i,3/2\](u[i,2]-u[i,1])/dy^2 at j=1, kf[i,NY-1/2\](u[i,NY-1]-u[i,NY])/dy^2 at j=NY), the flux at each y-wall dropped to zero. The qYL/qYR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_neumann_bc. `dy` is the grid's consumer-supplied y-spacing free name. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the second term. Interior O(dy^2); the prescribed-flux faces impose the exact boundary flux, preserving global 2nd-order L2 convergence for Neumann-compatible solutions. A model imports exactly one second-derivative/Laplacian rule for y.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (Neumann boundary conditions by prescribing the domain-face flux).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_y_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_y_neumann_bc) — status `active`
- convergence: [`heat_2d_varcoeff_y_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_neumann) — status `active`
- simulation: [`heat_2d_varcoeff_y_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_y_neumann) — status `active`

#### Convergence — case [`heat_2d_varcoeff_y_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_y_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_y_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 1.293e-03 | — | 8.851e-04 | — |
| 32 | 3.147e-04 | 2.04 | 2.165e-04 | 2.03 |
| 64 | 7.814e-05 | 2.01 | 5.354e-05 | 2.02 |
| 128 | 1.950e-05 | 2.00 | 1.338e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_y_neumann (log-log)](/plots/convergence/heat_2d_varcoeff_y_neumann.png)

### `varcoeff_laplacian_y_periodic`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:periodic</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_periodic.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_periodic.esm)

**Rewrites** `D(k * D(u)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "y"}]}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = varcoeff_laplacian_y_interior&lt;k=k, u=u&gt;, [1:NX, 1:1] = sum[i, j] ((0.5 * (k[i, 1] + k[i, 2]) * (u[i, 2] - u[i, 1]) - 0.5 * (k[i, NY] + k[i, 1]) * (u[i, 1] - u[i, NY])) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] ((0.5 * (k[i, NY] + k[i, 1]) * (u[i, 1] - u[i, NY]) - 0.5 * (k[i, NY - 1] + k[i, NY]) * (u[i, NY] - u[i, NY - 1])) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with a PERIODIC y-boundary, over the full 2-D field u[x,y]. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:y)]}, wrt:y ) -- outer-D of (k times the inner y-gradient) -- at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative (exactly as central_D2y_periodic fires on D(D(u,y),y)). This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_periodic: the operated axis y gets an interior band plus two wrapped y-faces, while the x axis is full-rank (all columns i treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state, inner D's argument) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write the coefficient on the LEFT: k * D(u, wrt:y). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on the yy-conservative-compound of bare fields declared over this 2-D grid's [x,y] axes; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together, so two instances coexist scoped to their own meshes. ONE rule serves two duals: the VARIABLE-COEFFICIENT case (k=k(x,y), a given periodic field) and the NONLINEAR case (k declared as an observed k=f(u)); the body only reads k by index at the faces, so a state-dependent k just makes the expanded stencil nonlinear in u with no reformulation (the keyed-factor field contract, like the grid's dy). AUTHORING CONTRACT for k: it is a declared [x,y]-shaped field; when supplied as an aggregate observed it MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because the coefficient observed is inlined into the stencil's [1,NX]x[2,NY-1] aggregate at evaluation and a shared loop-index name is captured (silently wrong interior values). Face conductivity is the ARITHMETIC MEAN kf[i,j+1/2]=(k[i,j]+k[i,j+1])/2 (O(dy^2) at the face for a differentiable coefficient, and affine in k so it composes with the nonlinear dual). Interior j in [2,NY-1] (all i): the imported conservative y-flux-difference stencil. The two y-face rows j=1 and j=NY carry the same conservative y-flux difference with the out-of-range y-neighbor AND its coefficient wrapped periodically: at j=1 the bottom y-face flux uses the ghost u[i,0]=u[i,NY] and k[i,0]=k[i,NY], giving kf[i,1/2]=(k[i,NY]+k[i,1])/2 and Ly u[i,1] = ( (k[i,1]+k[i,2])/2 (u[i,2]-u[i,1]) - (k[i,NY]+k[i,1])/2 (u[i,1]-u[i,NY]) ) / dy^2; at j=NY the top y-face flux uses u[i,NY+1]=u[i,1] and k[i,NY+1]=k[i,1], giving kf[i,NY+1/2]=(k[i,NY]+k[i,1])/2 and Ly u[i,NY] = ( (k[i,NY]+k[i,1])/2 (u[i,1]-u[i,NY]) - (k[i,NY-1]+k[i,NY])/2 (u[i,NY]-u[i,NY-1]) ) / dy^2. The SAME wrapped y-face flux (k[i,NY]+k[i,1])/2 (u[i,1]-u[i,NY])/dy is the bottom face of row 1 and the top face of row NY, so the discrete y-operator is exactly conservative on each column (the circulant variable-coefficient y-Laplacian) and its column sums vanish. The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_neumann_bc so every binding's simulator fills the (NX x 1) rows without rank promotion. `dy` is the grid's consumer-supplied y-spacing free name; the periodic wrap is dy-scale-invariant, so the rule serves any y-extent unchanged. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the second term. Reduces to central_D2y_periodic when k is spatially constant. Interior and periodic faces are all O(dy^2), so the whole y-operator is 2nd order with no boundary reduction. A model imports exactly one second-derivative/Laplacian rule for y.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 10.4 (centered differences on periodic domains).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_y_periodic) — status `active`
- convergence: [`heat_2d_varcoeff_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_periodic) — status `active`
- simulation: [`heat_2d_varcoeff_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_y_periodic) — status `active`

#### Convergence — case [`heat_2d_varcoeff_y_periodic`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_periodic)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_y_periodic/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_y_periodic/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 4.046e-02 | — | 8.114e-04 | — |
| 32 | 9.975e-03 | 2.02 | 2.028e-04 | 2.00 |
| 64 | 2.485e-03 | 2.01 | 5.069e-05 | 2.00 |
| 128 | 6.207e-04 | 2.00 | 1.267e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_y_periodic (log-log)](/plots/convergence/heat_2d_varcoeff_y_periodic.png)

### `varcoeff_laplacian_y_zero_flux_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:zero_flux</code> <code>axes:x,y</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_zero_flux_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/varcoeff_laplacian_y_zero_flux_bc.esm)

**Rewrites** `D(k * D(u)/Dy)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["k", {"op": "D", "args": ["u"], "wrt": "y"}]}], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>; <code>k</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NX, 2:NY - 1] = varcoeff_laplacian_y_interior&lt;k=k, u=u&gt;, [1:NX, 1:1] = sum[i, j] (0.5 * (k[i, 1] + k[i, 2]) * (u[i, 2] - u[i, 1]) / (dy * dy)) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] (0.5 * (k[i, NY - 1] + k[i, NY]) * (u[i, NY - 1] - u[i, NY]) / (dy * dy)) where {i in x, j in NY:NY})</div>

The y-half of a 2-D conservative finite-volume VARIABLE-COEFFICIENT / NONLINEAR Laplacian div(k grad u) on cartesian_uniform_2d with zero-flux (homogeneous Neumann) conditions on the two y-walls, over the full 2-D field u[x,y], cell-centered. The operator carries a coefficient, so a bare laplacian match cannot express it; the rule matches the CONSERVATIVE COMPOUND the consumer writes, D( {op:*, args:[k, D(u,wrt:y)]}, wrt:y ) at priority 10, so under esm-spec 9.6.3 outermost-first compound precedence it fires on the whole compound before any plain-D rule could lower the inner derivative (exactly as central_D2y fires on D(D(u,y),y)). This is the 2-D per-axis lift of grids/cartesian_uniform_1d/rules/varcoeff_laplacian_zero_flux_bc: the operated axis y gets an interior band plus two dropped-flux y-faces, while the x axis is full-rank (all columns i treated identically). BOTH operands are match wildcards bound to consumer fields: `u` (the state) and `k` (the coefficient, the LEFT factor of the product); matching is structural and non-commutative, so the consumer MUST write k * D(u, wrt:y). Match-scoped with where {u:{shape:[x,y]}, k:{shape:[x,y]}} (esm-spec 9.6.1) so it fires only on this 2-D grid's declared [x,y]-shaped fields; under import-edge renaming (esm-spec 9.7.7) the wrt literals and both where shapes follow the renamed axes together. ONE rule serves the VARIABLE-COEFFICIENT case (k=k(x,y)) and the NONLINEAR case (k=f(u)); the body only reads k by index at the faces. AUTHORING CONTRACT for k: k is a declared [x,y]-shaped field; an aggregate-observed coefficient MUST loop over indices DISTINCT from the stencil's interior indices `i`,`j` (the library convention is `gi`/`gj`), because it is inlined into the stencil's [1,NX]x[2,NY-1] aggregate and a shared loop-index name is captured (silently wrong interior values). Face conductivity is the ARITHMETIC MEAN kf[i,j+1/2]=(k[i,j]+k[i,j+1])/2 (O(dy^2) at the face for a differentiable coefficient, affine in k so it composes with the nonlinear dual). Interior j in [2,NY-1] (all i): the imported conservative y-flux-difference stencil. The y-boundary faces impose zero flux exactly (natural for the finite-volume form -- no ghost geometry): at j=1 the bottom y-domain-face flux is dropped, giving Ly u[i,1] = kf[i,3/2\](u[i,2]-u[i,1])/dy^2 = (k[i,1]+k[i,2])/2 (u[i,2]-u[i,1])/dy^2; at j=NY the top y-domain-face flux is dropped, giving Ly u[i,NY] = -kf[i,NY-1/2\](u[i,NY]-u[i,NY-1])/dy^2 = (k[i,NY-1]+k[i,NY])/2 (u[i,NY-1]-u[i,NY])/dy^2. This is exactly varcoeff_laplacian_y_neumann_bc at qYL=qYR=0. The face values are FULL-RANK aggregates (singleton output range on the y axis, all i) mirroring central_D2y_neumann_bc. `dy` is the grid's consumer-supplied y-spacing free name; the zero-flux closure is dy-scale-invariant, so the rule serves any y-extent unchanged. The full variable-coefficient 2-D Laplacian is assembled in the consuming equation as D(k*D(u,x),x) + D(k*D(u,y),y): this rule lowers the second term. Reduces to the constant-coefficient zero-gradient y-half when k is spatially constant. Interior O(dy^2); the one-sided zero-flux faces are O(dy) locally, which preserves global 2nd-order L2 convergence for zero-flux-compatible solutions (flux k u_y vanishing at both y-walls). A model imports exactly one second-derivative/Laplacian rule for y.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.15 (variable-coefficient conservative form) and Section 2.12 (Neumann boundary conditions by dropping the domain-face flux).

**Conformance fixtures**

- ast: [`varcoeff_laplacian_y_zero_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_laplacian_y_zero_flux_bc) — status `active`
- convergence: [`heat_2d_varcoeff_y_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_zero_flux) — status `active`
- simulation: [`heat_2d_varcoeff_y_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_2d_varcoeff_y_zero_flux) — status `active`

#### Convergence — case [`heat_2d_varcoeff_y_zero_flux`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_2d_varcoeff_y_zero_flux)

Error norms read from the committed golden [`tests/conformance/convergence/heat_2d_varcoeff_y_zero_flux/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_2d_varcoeff_y_zero_flux/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order | Linf_error | observed order |
|---|---|---|---|---|
| 16 | 6.726e-03 | — | 2.584e-03 | — |
| 32 | 1.692e-03 | 1.99 | 6.522e-04 | 1.99 |
| 64 | 4.236e-04 | 2.00 | 1.630e-04 | 2.00 |
| 128 | 1.060e-04 | 2.00 | 4.078e-05 | 2.00 |

Expected order: **2** (± 0.2).

![Convergence of heat_2d_varcoeff_y_zero_flux (log-log)](/plots/convergence/heat_2d_varcoeff_y_zero_flux.png)

### `weno3_eno1_D1x_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno3_eno1_D1x_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno3_eno1_D1x_extrapolate_bc.esm)

**Rewrites** `D(u)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">weno5_eno1_D1x_extrapolate_frame&lt;interior=weno3_eno1_D1x_interior&lt;u=u&gt;, u=u&gt;</div>

Hybrid WENO3/ENO1 upwind first derivative d u/dx on cartesian_uniform_2d with the linear-extrapolation halo closure -- the per-axis derivative diff2x of NCAR fire_behavior's fire_upwinding = 8 level-set scheme: third-order Jiang-Shu WENO3 in a band of weno_band_cells cells around the zero level set, first-order ENO1 (minmod) elsewhere and within bdy_eno1_cells of the walls, upwind side chosen by the wind-relative sign signo = (grad4 u . (wind_u, wind_v)) d4 u/dx. Calc_tend_ls case 8 differs from the default case 9 only in calling Select_weno3 instead of Select_weno5, and this rule differs from rules/weno5_eno1_D1x_extrapolate_bc in the same single place: its body applies the same frame, rules/weno5_eno1_D1x_extrapolate_frame (ENO1 strips, guarded-extrapolation wall lines; lint-tiled in its own file), around stencils/weno3_eno1_D1x_interior instead of the WENO5 interior. Match D(u, wrt: x) at priority 10, scoped with where {u: {shape: [x, y]}} (esm-spec 9.6.1). It matches the same pattern at the same priority as the WENO5 sibling, so a document imports one of the two (esm-spec 9.6.8: choosing a scheme = choosing a rule). FREE NAMES: wind_u, wind_v ([x, y] fields), bdy_eno1_cells (>= 3, as for the WENO5 sibling: the shared frame's ENO1 strips are the outermost three cells; fire_behavior's BDY_ENO1 = 10), weno_band_cells (fire_lsm_band_ngp), and the grid's dx, dy. Requires NX, NY >= 7. ORDER: formally first (ENO1 strips and far field); third in the smooth part of the WENO band. Exact on a linear field except at the wall the field decreases toward, where the guarded ghost clamps and ENO1 returns 0 -- fire_behavior's own wall behaviour.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno3_eno1_D1x_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno3_eno1_D1x_extrapolate_bc) — status `active`

### `weno3_eno1_D1y_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno3_eno1_D1y_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno3_eno1_D1y_extrapolate_bc.esm)

**Rewrites** `D(u)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">weno5_eno1_D1y_extrapolate_frame&lt;interior=weno3_eno1_D1y_interior&lt;u=u&gt;, u=u&gt;</div>

Hybrid WENO3/ENO1 upwind first derivative d u/dy on cartesian_uniform_2d with the linear-extrapolation halo closure -- the per-axis derivative diff2y of NCAR fire_behavior's fire_upwinding = 8 level-set scheme: third-order Jiang-Shu WENO3 in a band of weno_band_cells cells around the zero level set, first-order ENO1 (minmod) elsewhere and within bdy_eno1_cells of the walls, upwind side chosen by the wind-relative sign signo = (grad4 u . (wind_u, wind_v)) d4 u/dy. Calc_tend_ls case 8 differs from the default case 9 only in calling Select_weno3 instead of Select_weno5, and this rule differs from rules/weno5_eno1_D1y_extrapolate_bc in the same single place: its body applies the same frame, rules/weno5_eno1_D1y_extrapolate_frame (ENO1 strips, guarded-extrapolation wall lines; lint-tiled in its own file), around stencils/weno3_eno1_D1y_interior instead of the WENO5 interior. Match D(u, wrt: y) at priority 10, scoped with where {u: {shape: [x, y]}} (esm-spec 9.6.1). It matches the same pattern at the same priority as the WENO5 sibling, so a document imports one of the two (esm-spec 9.6.8: choosing a scheme = choosing a rule). FREE NAMES: wind_u, wind_v ([x, y] fields), bdy_eno1_cells (>= 3, as for the WENO5 sibling: the shared frame's ENO1 strips are the outermost three cells; fire_behavior's BDY_ENO1 = 10), weno_band_cells (fire_lsm_band_ngp), and the grid's dx, dy. Requires NX, NY >= 7. ORDER: formally first (ENO1 strips and far field); third in the smooth part of the WENO band. Exact on a linear field except at the wall the field decreases toward, where the guarded ghost clamps and ENO1 returns 0 -- fire_behavior's own wall behaviour.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno3_eno1_D1y_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno3_eno1_D1y_extrapolate_bc) — status `active`

### `weno3_eno1_norm_D2_extrapolate_bc`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 5</code> <code>priority 20</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno3_eno1_norm_D2_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno3_eno1_norm_D2_extrapolate_bc.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 20 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sqrt(weno5_eno1_D1x_extrapolate_frame&lt;interior=weno3_eno1_D1x_interior&lt;u=u&gt;, u=u&gt;^2 + weno5_eno1_D1y_extrapolate_frame&lt;interior=weno3_eno1_D1y_interior&lt;u=u&gt;, u=u&gt;^2)</div>

Hybrid WENO3/ENO1 gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) for the level-set term psi_t = -S|grad psi| on cartesian_uniform_2d with the linear-extrapolation halo closure -- NCAR fire_behavior's fire_upwinding = 8: grad = sqrt(diff2x^2 + diff2y^2) with diff2x, diff2y the per-axis derivatives of rules/weno3_eno1_D1x_extrapolate_bc and rules/weno3_eno1_D1y_extrapolate_bc. COMPOUND MATCH (esm-spec 9.6.3): the whole compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) at priority 20, above the per-axis rules' 10, scoped with where {u: {shape: [x, y]}}. The body is sqrt(Mx^2 + My^2) with Mx, My exactly the per-axis bodies: the shared frames rules/weno5_eno1_D1x_extrapolate_frame and rules/weno5_eno1_D1y_extrapolate_frame around stencils/weno3_eno1_D1x_interior and stencils/weno3_eno1_D1y_interior. It is rules/weno5_eno1_norm_D2_extrapolate_bc (fire_upwinding = 9) with the two interior stencils swapped, Select_weno3 for Select_weno5, the only difference between Calc_tend_ls cases 8 and 9. SIBLINGS: rules/weno5_eno1_norm_D2_extrapolate_bc, rules/weno5_eno1_blend_norm_D2_extrapolate_bc and rules/weno5_eno1_signed_norm_D2_extrapolate_bc match the same compound at the same priority; a component imports exactly one (if it imported two, the first declared would win every match by the 9.6.3 declaration-order tie-break). The front-distance switch |u[i,j]| < weno_band_cells * dx is on the field itself, so a signed-distance level set is required. FREE NAMES: wind_u, wind_v ([x, y] fields), bdy_eno1_cells (>= 3; fire_behavior's BDY_ENO1 = 10), weno_band_cells, and the grid's dx, dy. Requires NX, NY >= 7. ORDER: formally first; third in the smooth WENO band. BOUNDARY CLOSURE bc:extrapolate_linear: one-cell guarded linear-extrapolation ghosts at the walls (stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys); exact on a linear field except on the wall lines the field decreases toward, where ENO1 returns 0 along that axis -- fire_behavior's own wall behaviour.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO/ENO1 level-set schemes (fire_upwinding = 8, 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls (case 8), Select_eno, Select_weno3, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno3_eno1_norm_D2_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno3_eno1_norm_D2_extrapolate_bc) — status `active`

### `weno5_eno1_D1x_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_D1x_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_D1x_extrapolate_bc.esm)

**Rewrites** `D(u)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">weno5_eno1_D1x_extrapolate_frame&lt;interior=weno5_eno1_D1x_interior&lt;u=u&gt;, u=u&gt;</div>

Hybrid WENO5/ENO1 upwind first derivative d u/dx on cartesian_uniform_2d with the linear-extrapolation halo closure -- the per-axis derivative diff2x of NCAR fire_behavior's default level-set scheme (fire_upwinding = 9; Munoz-Esparza et al. 2018): fifth-order Jiang-Shu WENO5 in a band of weno_band_cells cells around the zero level set, first-order ENO1 (minmod) elsewhere and within bdy_eno1_cells of the walls, upwind side of the WENO reconstruction chosen by the wind-relative propagation sign signo = (grad4 u . (wind_u, wind_v)) d4 u/dx. This is the quantity fire_behavior feeds both to |grad psi| = sqrt(diff2x^2 + diff2y^2) and to the front normal nvx = diff2x/sqrt(|grad psi|^2 + eps); rules/weno5_eno1_norm_D2_extrapolate_bc packages the norm as a compound match so a consumer may write the norm as sqrt((D(u,x))^2 + (D(u,y))^2) and the normal's ingredients as bare D(u, wrt: x) / D(u, wrt: y) with the same scheme. Match D(u, wrt: x) at priority 10 (before the priority-less centered central_D1x_* rules if a document imports both), scoped with where {u: {shape: [x, y]}} (esm-spec 9.6.1). THE SCHEME: the body applies rules/weno5_eno1_D1x_extrapolate_frame -- the seven-region makearray shared by this family, with its ENO1 strips and guarded-extrapolation wall lines -- around the imported stencils/weno5_eno1_D1x_interior on the doubly-interior block [4, NX-3] x [3, NY-2] (the hybrid switch on the field and the index, see that stencil). The frame's regions are lint-tiled (L006) in the frame's own file. BOUNDARY CLOSURE bc:extrapolate_linear: one-cell guarded linear-extrapolation ghosts at the x walls -- fire_behavior's Extrapol_var_at_bdys; see the frame. FREE NAMES (consumer-supplied, resolved in the consuming model's scope like the grid geometry dx, dy): wind_u, wind_v ([x, y] fields, the advecting wind), bdy_eno1_cells (>= 3; fire_behavior's BDY_ENO1 = 10), weno_band_cells (fire_behavior's fire_lsm_band_ngp = 4). Both ifelse branches are evaluated on every interior cell (esm ifelse is not lazy); the region structure keeps every evaluated index in range. Requires NX >= 7, NY >= 5 (respectively NY >= 7, NX >= 5 for the y rule). ORDER: formally first (ENO1 strips and far field); fifth in the smooth part of the WENO band away from critical points. EXACT on a linear field u = a x + b y + c on every cell except the wall toward which the field decreases along x (i = 1 for a > 0, i = NX for a < 0): there the guarded extrapolation clamps the ghost to the inner neighbour, so the two one-sided differences have opposite signs and the ENO1 minmod returns 0 instead of a -- fire_behavior's own wall behaviour, reproduced (the ghost is the plain linear one, and the result exact, only at a wall the field increases toward). Tagged order:1. The consumer's units checker sees a derivative of u per unit length.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno5_eno1_D1x_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno5_eno1_D1x_extrapolate_bc) — status `active`

### `weno5_eno1_D1x_extrapolate_frame`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_D1x_extrapolate_frame.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_D1x_extrapolate_frame.esm)

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([4:NX - 3, 3:NY - 2] = interior, [4:NX - 3, 1:2] = sum[i, j] (eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;) where {i in 4:NX - 3, j in 1:2}, [4:NX - 3, NY - 1:NY] = sum[i, j] (eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;) where {i in 4:NX - 3, j in NY - 1:NY}, [2:3, 1:NY] = sum[i, j] (eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;) where {i in 2:3, j in y}, [NX - 2:NX - 1, 1:NY] = sum[i, j] (eno1_minmod_D1x&lt;i=i, j=j, u=u&gt;) where {i in NX - 2:NX - 1, j in y}, [1:1, 1:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[1, j] - extrapolate_linear_ghost&lt;a=u[1, j], b=u[2, j]&gt;) / dx, dr=(u[1 + 1, j] - u[1, j]) / dx&gt;) where {i in 1:1, j in y}, [NX:NX, 1:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[NX, j] - u[NX - 1, j]) / dx, dr=(extrapolate_linear_ghost&lt;a=u[NX, j], b=u[NX - 1, j]&gt; - u[NX, j]) / dx&gt;) where {i in NX:NX, j in y})</div>

The boundary frame shared by every hybrid WENO5/ENO1 derivative d u/dx on cartesian_uniform_2d with the linear-extrapolation halo closure: a seven-region makearray tiling [x, y] whose doubly-interior block [4, NX-3] x [3, NY-2] is the caller-supplied `interior` (an interior-only stencil over exactly that block: stencils/weno5_eno1_D1x_interior for fire_upwinding = 9, stencils/weno5_eno1_blend_D1x_interior for 10, stencils/weno5_eno1_signed_D1x_interior for the reinitialisation, stencils/weno3_eno1_D1x_interior for fire_upwinding = 8, whose shorter WENO3 reach fits the same block), and whose other six regions are the same for all of them because they are ENO1 whatever the scheme -- the two-cell strips [4, NX-3] x [1, 2] and [4, NX-3] x [NY-1, NY] where the fourth-order selector along y would leave the domain, and the strips [2, 3] x [1, NY] and [NX-2, NX-1] x [1, NY] where the seven-point WENO stencil along x would (stencils/eno1_minmod_D1x; the reference code is ENO1 there by its bdy_eno1_cells rule, which must be >= 3 for the WENO5 stencil never to reach a wall; fire_behavior fixes it at BDY_ENO1 = 10, so these strips are ENO1 in the reference code for every scheme, WENO3 included), and the two wall lines i = 1 and i = NX. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the ENO1 (three-point) stencil is ever evaluated at a wall cell, so one ghost suffices, exactly as in the reference code, whose halo is one cell for the same reason. Corner cells use the edge ghost of their x wall (the reference code's diagonal corner ghost is never read by ENO1). At the wall the field decreases toward, the guarded extrapolation clamps the ghost to the inner neighbour, the two one-sided differences have opposite signs, and the ENO1 minmod returns 0 -- fire_behavior's own wall behaviour, reproduced. A match-less rule applied by the D1x and norm rules of this family (rules/weno5_eno1_D1x_extrapolate_bc, rules/weno5_eno1_norm_D2_extrapolate_bc, rules/weno5_eno1_blend_norm_D2_extrapolate_bc, rules/weno5_eno1_signed_norm_D2_extrapolate_bc, rules/weno3_eno1_D1x_extrapolate_bc, rules/weno3_eno1_norm_D2_extrapolate_bc), so the boundary treatment has one spelling. Params: u, the [x, y] field; interior, its interior-block derivative. Requires NX >= 7, NY >= 5.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- none yet

### `weno5_eno1_D1y_extrapolate_bc`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 10</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_D1y_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_D1y_extrapolate_bc.esm)

**Rewrites** `D(u)/Dy` at priority 10 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">weno5_eno1_D1y_extrapolate_frame&lt;interior=weno5_eno1_D1y_interior&lt;u=u&gt;, u=u&gt;</div>

Hybrid WENO5/ENO1 upwind first derivative d u/dy on cartesian_uniform_2d with the linear-extrapolation halo closure -- the per-axis derivative diff2y of NCAR fire_behavior's default level-set scheme (fire_upwinding = 9; Munoz-Esparza et al. 2018): fifth-order Jiang-Shu WENO5 in a band of weno_band_cells cells around the zero level set, first-order ENO1 (minmod) elsewhere and within bdy_eno1_cells of the walls, upwind side of the WENO reconstruction chosen by the wind-relative propagation sign signo = (grad4 u . (wind_u, wind_v)) d4 u/dy. This is the quantity fire_behavior feeds both to |grad psi| = sqrt(diff2x^2 + diff2y^2) and to the front normal nvy = diff2y/sqrt(|grad psi|^2 + eps); rules/weno5_eno1_norm_D2_extrapolate_bc packages the norm as a compound match so a consumer may write the norm as sqrt((D(u,x))^2 + (D(u,y))^2) and the normal's ingredients as bare D(u, wrt: x) / D(u, wrt: y) with the same scheme. Match D(u, wrt: y) at priority 10 (before the priority-less centered central_D1y_* rules if a document imports both), scoped with where {u: {shape: [x, y]}} (esm-spec 9.6.1). THE SCHEME: the body applies rules/weno5_eno1_D1y_extrapolate_frame -- the seven-region makearray shared by this family, with its ENO1 strips and guarded-extrapolation wall lines -- around the imported stencils/weno5_eno1_D1y_interior on the doubly-interior block [3, NX-2] x [4, NY-3] (the hybrid switch on the field and the index, see that stencil). The frame's regions are lint-tiled (L006) in the frame's own file. BOUNDARY CLOSURE bc:extrapolate_linear: one-cell guarded linear-extrapolation ghosts at the y walls -- fire_behavior's Extrapol_var_at_bdys; see the frame. FREE NAMES (consumer-supplied, resolved in the consuming model's scope like the grid geometry dx, dy): wind_u, wind_v ([x, y] fields, the advecting wind), bdy_eno1_cells (>= 3; fire_behavior's BDY_ENO1 = 10), weno_band_cells (fire_behavior's fire_lsm_band_ngp = 4). Both ifelse branches are evaluated on every interior cell (esm ifelse is not lazy); the region structure keeps every evaluated index in range. Requires NX >= 7, NY >= 5 (respectively NY >= 7, NX >= 5 for the y rule). ORDER: formally first (ENO1 strips and far field); fifth in the smooth part of the WENO band away from critical points. EXACT on a linear field u = a x + b y + c on every cell except the wall toward which the field decreases along y (j = 1 for b > 0, j = NY for b < 0): there the guarded extrapolation clamps the ghost to the inner neighbour, so the two one-sided differences have opposite signs and the ENO1 minmod returns 0 instead of b -- fire_behavior's own wall behaviour, reproduced (the ghost is the plain linear one, and the result exact, only at a wall the field increases toward). Tagged order:1. The consumer's units checker sees a derivative of u per unit length.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno5_eno1_D1y_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno5_eno1_D1y_extrapolate_bc) — status `active`

### `weno5_eno1_D1y_extrapolate_frame`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_D1y_extrapolate_frame.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_D1y_extrapolate_frame.esm)

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([3:NX - 2, 4:NY - 3] = interior, [1:2, 4:NY - 3] = sum[i, j] (eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;) where {i in 1:2, j in 4:NY - 3}, [NX - 1:NX, 4:NY - 3] = sum[i, j] (eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;) where {i in NX - 1:NX, j in 4:NY - 3}, [1:NX, 2:3] = sum[i, j] (eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;) where {i in x, j in 2:3}, [1:NX, NY - 2:NY - 1] = sum[i, j] (eno1_minmod_D1y&lt;i=i, j=j, u=u&gt;) where {i in x, j in NY - 2:NY - 1}, [1:NX, 1:1] = sum[i, j] (eno1_minmod&lt;dl=(u[i, 1] - extrapolate_linear_ghost&lt;a=u[i, 1], b=u[i, 2]&gt;) / dy, dr=(u[i, 1 + 1] - u[i, 1]) / dy&gt;) where {i in x, j in 1:1}, [1:NX, NY:NY] = sum[i, j] (eno1_minmod&lt;dl=(u[i, NY] - u[i, NY - 1]) / dy, dr=(extrapolate_linear_ghost&lt;a=u[i, NY], b=u[i, NY - 1]&gt; - u[i, NY]) / dy&gt;) where {i in x, j in NY:NY})</div>

The boundary frame shared by every hybrid WENO5/ENO1 derivative d u/dy on cartesian_uniform_2d with the linear-extrapolation halo closure: a seven-region makearray tiling [x, y] whose doubly-interior block [3, NX-2] x [4, NY-3] is the caller-supplied `interior` (an interior-only stencil over exactly that block: stencils/weno5_eno1_D1y_interior for fire_upwinding = 9, stencils/weno5_eno1_blend_D1y_interior for 10, stencils/weno5_eno1_signed_D1y_interior for the reinitialisation, stencils/weno3_eno1_D1y_interior for fire_upwinding = 8, whose shorter WENO3 reach fits the same block), and whose other six regions are the same for all of them because they are ENO1 whatever the scheme -- the two-cell strips [1, 2] x [4, NY-3] and [NX-1, NX] x [4, NY-3] where the fourth-order selector along x would leave the domain, and the strips [1, NX] x [2, 3] and [1, NX] x [NY-2, NY-1] where the seven-point WENO stencil along y would (stencils/eno1_minmod_D1y; the reference code is ENO1 there by its bdy_eno1_cells rule, which must be >= 3 for the WENO5 stencil never to reach a wall; fire_behavior fixes it at BDY_ENO1 = 10, so these strips are ENO1 in the reference code for every scheme, WENO3 included), and the two wall lines j = 1 and j = NY. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[i,0] = max(2 u[i,1] - u[i,2], u[i,1], u[i,2]) (and symmetrically u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the ENO1 (three-point) stencil is ever evaluated at a wall cell, so one ghost suffices, exactly as in the reference code, whose halo is one cell for the same reason. Corner cells use the edge ghost of their y wall (the reference code's diagonal corner ghost is never read by ENO1). At the wall the field decreases toward, the guarded extrapolation clamps the ghost to the inner neighbour, the two one-sided differences have opposite signs, and the ENO1 minmod returns 0 -- fire_behavior's own wall behaviour, reproduced. A match-less rule applied by the D1y and norm rules of this family (rules/weno5_eno1_D1y_extrapolate_bc, rules/weno5_eno1_norm_D2_extrapolate_bc, rules/weno5_eno1_blend_norm_D2_extrapolate_bc, rules/weno5_eno1_signed_norm_D2_extrapolate_bc, rules/weno3_eno1_D1y_extrapolate_bc, rules/weno3_eno1_norm_D2_extrapolate_bc), so the boundary treatment has one spelling. Params: u, the [x, y] field; interior, its interior-block derivative. Requires NY >= 7, NX >= 5.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- none yet

### `weno5_eno1_blend_norm_D2_extrapolate_bc`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 20</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_blend_norm_D2_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_blend_norm_D2_extrapolate_bc.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 20 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sqrt(weno5_eno1_D1x_extrapolate_frame&lt;interior=weno5_eno1_blend_D1x_interior&lt;u=u&gt;, u=u&gt;^2 + weno5_eno1_D1y_extrapolate_frame&lt;interior=weno5_eno1_blend_D1y_interior&lt;u=u&gt;, u=u&gt;^2)</div>

Hybrid WENO5/ENO1 gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) of NCAR fire_behavior's fire_upwinding = 10 (physics/level_set_mod.F90 Calc_tend_ls case 10), the blending-zone variant of the default fire_upwinding = 9 scheme, on cartesian_uniform_2d with the linear-extrapolation halo closure. It differs from rules/weno5_eno1_norm_D2_extrapolate_bc only in its interior stencils, stencils/weno5_eno1_blend_D1x_interior and stencils/weno5_eno1_blend_D1y_interior, which replace the hard band edge with the reference's cosine mask (stencils/weno5_eno1_cosine_mask): pure WENO5 for |u| <= B, a cosine taper over B < |u| <= 2B, pure ENO1 beyond (B = weno_band_cells dx), with the ENO1 value SELECTED wherever the mask is 0 (stencils/weno5_eno1_masked_blend), so the far field and the wall band are bitwise ENO1, as in the reference. COMPOUND MATCH (esm-spec 9.6.3): the whole compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) at priority 20, scoped with where {u: {shape: [x, y]}} -- the same pattern and priority as rules/weno5_eno1_norm_D2_extrapolate_bc, because it is another discretization of the same operator |grad u|. Which one lowers a component's norm is chosen by which rule that component imports (esm-spec 9.6.8: choosing a scheme = choosing a rule); a component imports one of them. Importing two is deterministic but pointless: the 9.6.3 declaration-order tie-break gives every match to the first declared, and the other never fires. Propagation and reinitialisation are separate components (a tendency model and a pseudo-time model), each importing its own. The body is sqrt(Mx^2 + My^2) with Mx, My the family's per-axis frames (rules/weno5_eno1_D1x_extrapolate_frame, rules/weno5_eno1_D1y_extrapolate_frame: the ENO1 strips and one-cell guarded linear-extrapolation ghosts of fire_behavior's Extrapol_var_at_bdys, lint-tiled in their own files) around this rule's own interior stencils. ORDER: tagged order:1 for the same reason as the sibling, and no convergence case is shipped. Free names: wind_u, wind_v, bdy_eno1_cells, weno_band_cells, pi_blend, dx, dy. Requires NX, NY >= 7.

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls case (10) (band_width = 2 threshold_hlu, transition_width = threshold_hlu, the three-branch cosine mask and the `if (mask > 0.0)` blend), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno5_eno1_blend_norm_D2_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno5_eno1_blend_norm_D2_extrapolate_bc) — status `active`

### `weno5_eno1_norm_D2_extrapolate_bc`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 20</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_norm_D2_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_norm_D2_extrapolate_bc.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 20 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sqrt(weno5_eno1_D1x_extrapolate_frame&lt;interior=weno5_eno1_D1x_interior&lt;u=u&gt;, u=u&gt;^2 + weno5_eno1_D1y_extrapolate_frame&lt;interior=weno5_eno1_D1y_interior&lt;u=u&gt;, u=u&gt;^2)</div>

Hybrid WENO5/ENO1 gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) for the level-set Hamilton-Jacobi term psi_t = -S|grad psi| on cartesian_uniform_2d with the linear-extrapolation halo closure -- NCAR fire_behavior's default level-set scheme (fire_upwinding = 9; Munoz-Esparza et al. 2018): grad = sqrt(diff2x^2 + diff2y^2) with diff2x, diff2y the per-axis hybrid upwind derivatives of rules/weno5_eno1_D1x_extrapolate_bc and rules/weno5_eno1_D1y_extrapolate_bc (fifth-order Jiang-Shu WENO5 within weno_band_cells cells of the zero level set with the upwind side chosen by the wind-relative propagation sign, first-order ENO1 minmod elsewhere and within bdy_eno1_cells of the walls, one-cell guarded linear-extrapolation ghosts at the walls). It is the fire_behavior-default sibling of rules/godunov_norm_D2_dirichlet_bc (fire_upwinding = 4): unlike the Godunov Hamiltonian, this norm is NOT a compound of max/min one-sided pieces -- it is exactly the root-sum-square of the two per-axis derivatives -- so the compound spelling and the per-axis spelling agree, and this rule exists so that a consumer who writes the norm as the compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) gets the scheme without importing the per-axis rules (and so that importing both cannot lower the compound differently). COMPOUND MATCH (esm-spec 9.6.3): the whole compound rooted at the outer sqrt at priority 20 -- above the per-axis rules' 10 -- so it fires first when a document imports both; scoped with where {u: {shape: [x, y]}}. The body is sqrt(Mx^2 + My^2) where Mx and My are exactly the two per-axis rules' bodies: rules/weno5_eno1_D1x_extrapolate_frame and rules/weno5_eno1_D1y_extrapolate_frame (seven regions each: the hybrid interior block, the ENO1 strips, the two ghost wall lines) around stencils/weno5_eno1_D1x_interior and stencils/weno5_eno1_D1y_interior. It is therefore not itself a makearray-rooted body; the two frames are lint-tiled (L006) in their own files. SIBLINGS: rules/weno5_eno1_blend_norm_D2_extrapolate_bc (fire_upwinding = 10) and rules/weno5_eno1_signed_norm_D2_extrapolate_bc (the reinitialisation) are the same frames around different interior stencils, and match the same compound at the same priority. The scheme is chosen by which one a component imports (esm-spec 9.6.8: choosing a scheme = choosing a rule); if a component imported two, the first declared would win every match by the 9.6.3 declaration-order tie-break and the other would never fire. THE FRONT-DISTANCE SWITCH IS ON THE FIELD ITSELF: |u[i,j]| < weno_band_cells * dx selects WENO5, so a signed-distance level set is required for the band to mean 'near the front' (fire_behavior reinitialises for this reason). FREE NAMES: wind_u, wind_v ([x, y] fields), bdy_eno1_cells (>= 3), weno_band_cells, and the grid's dx, dy. Requires NX, NY >= 7. ORDER: formally first; fifth in the smooth WENO band. EXACT on a linear field u = a x + b y + c (norm sqrt(a^2 + b^2)) on every cell except the wall lines toward which the field decreases (i = 1 for a > 0, i = NX for a < 0, likewise j with b): there the guarded extrapolation clamps the ghost, the ENO1 minmod along that axis returns 0, and the norm is O(1) wrong -- fire_behavior's own wall behaviour, reproduced. BOUNDARY CLOSURE bc:extrapolate_linear: the walls carry no prescribed value; instead each wall cell sees a one-cell halo ghost extrapolated from the two innermost cells, u[0,j] = max(2 u[1,j] - u[2,j], u[1,j], u[2,j]) (and symmetrically u[NX+1,j], u[i,0], u[i,NY+1]) -- stencils/extrapolate_linear_ghost, fire_behavior's Extrapol_var_at_bdys (HALO = 1, BIAS = 1). Only the ENO1 (three-point) stencil is ever evaluated at a wall cell, so one ghost suffices, exactly as in the reference code, whose halo is one cell for the same reason. Corner cells use the two edge ghosts of their walls (the reference code's diagonal corner ghost is never read by ENO1).

**References**

- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- Osher, S., & Sethian, J. A. (1988). Fronts propagating with curvature-dependent speed: Algorithms based on Hamilton-Jacobi formulations. Journal of Computational Physics, 79(1), 12-49. [`doi:10.1016/0021-9991(88)90002-2`](https://doi.org/10.1016/0021-9991%2888%2990002-2)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Calc_tend_ls, Select_eno, Select_weno5, Select_4th, Select_godunov (case 4), Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno5_eno1_norm_D2_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno5_eno1_norm_D2_extrapolate_bc) — status `active`
- convergence: [`levelset_weno5_eno1_2d_bump_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/levelset_weno5_eno1_2d_bump_mms) — status `active`
- simulation: [`levelset_weno5_eno1_2d_bump_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/levelset_weno5_eno1_2d_bump_mms) — status `active`

#### Convergence — case [`levelset_weno5_eno1_2d_bump_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/levelset_weno5_eno1_2d_bump_mms)

Error norms read from the committed golden [`tests/conformance/convergence/levelset_weno5_eno1_2d_bump_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/levelset_weno5_eno1_2d_bump_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.000e-01 | — |
| 32 | 3.675e-02 | 1.44 |
| 64 | 1.325e-02 | 1.47 |
| 128 | 4.743e-03 | 1.48 |

Expected order: **1.5** (± 0.35).

![Convergence of levelset_weno5_eno1_2d_bump_mms (log-log)](/plots/convergence/levelset_weno5_eno1_2d_bump_mms.png)

### `weno5_eno1_signed_norm_D2_extrapolate_bc`

<p class="rule-meta"><code>op:grad_norm</code> <code>order:1</code> <code>bc:extrapolate_linear</code> <code>axes:x,y</code> <code>stencil width 7</code> <code>priority 20</code></p>

Source: [`grids/cartesian_uniform_2d/rules/weno5_eno1_signed_norm_D2_extrapolate_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_uniform_2d/rules/weno5_eno1_signed_norm_D2_extrapolate_bc.esm)

**Rewrites** `sqrt(D(u)/Dx^2 + D(u)/Dy^2)` at priority 20 — match pattern:

```json
{"op": "sqrt", "args": [{"op": "+", "args": [{"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "x"}, 2]}, {"op": "^", "args": [{"op": "D", "args": ["u"], "wrt": "y"}, 2]}]}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[x, y]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">sqrt(weno5_eno1_D1x_extrapolate_frame&lt;interior=weno5_eno1_signed_D1x_interior&lt;s=upwind_sign, u=u&gt;, u=u&gt;^2 + weno5_eno1_D1y_extrapolate_frame&lt;interior=weno5_eno1_signed_D1y_interior&lt;s=upwind_sign, u=u&gt;, u=u&gt;^2)</div>

Hybrid WENO5/ENO1 gradient magnitude |grad u| = sqrt(u_x^2 + u_y^2) for the level-set REINITIALISATION equation psi_tau = S(psi_0)(1 - |grad psi|) of NCAR fire_behavior (physics/level_set_mod.F90 Advance_ls_reinit, fire_upwinding_reinit = 4; Sussman, Smereka & Osher 1994 as adapted by Munoz-Esparza et al. 2018), on cartesian_uniform_2d with the linear-extrapolation halo closure. It differs from rules/weno5_eno1_norm_D2_extrapolate_bc only in its interior stencils, stencils/weno5_eno1_signed_D1x_interior and stencils/weno5_eno1_signed_D1y_interior: the upwind side is chosen by the frozen sign field instead of the wind, and the band test that selects WENO5 over ENO1 is SIGNED rather than absolute, as the reference writes it (see the stencils, including why the missing abs() is reproduced rather than corrected). The sign field is the free name `upwind_sign` ([x, y]; fire_behavior's lfn_s0, S(psi_0) = psi_0 / sqrt(psi_0^2 + dx^2) frozen at the start of the pseudo-time iteration), bound to the stencils' param s; the rule references no wind_u/wind_v. COMPOUND MATCH (esm-spec 9.6.3): the whole compound sqrt((D(u, wrt: x))^2 + (D(u, wrt: y))^2) at priority 20, scoped with where {u: {shape: [x, y]}} -- the same pattern and priority as rules/weno5_eno1_norm_D2_extrapolate_bc, because it is another discretization of the same operator |grad u|. Which one lowers a component's norm is chosen by which rule that component imports (esm-spec 9.6.8: choosing a scheme = choosing a rule); a component imports one of them. Importing two is deterministic but pointless: the 9.6.3 declaration-order tie-break gives every match to the first declared, and the other never fires. Propagation and reinitialisation are separate components (a tendency model and a pseudo-time model), each importing its own. The body is sqrt(Mx^2 + My^2) with Mx, My the family's per-axis frames (rules/weno5_eno1_D1x_extrapolate_frame, rules/weno5_eno1_D1y_extrapolate_frame: the ENO1 strips and one-cell guarded linear-extrapolation ghosts of fire_behavior's Extrapol_var_at_bdys, lint-tiled in their own files) around this rule's own interior stencils. ORDER: tagged order:1 for the same reason as its sibling -- the ENO1 branch and the wall closure are first-order, and the fifth-order branch applies only inside the band, so no whole-field convergence claim is made and no convergence case is shipped. Free names: upwind_sign, bdy_eno1_cells, weno_band_cells, dx, dy. Requires NX, NY >= 7.

**References**

- Sussman, M., Smereka, P., & Osher, S. (1994). A level set approach for computing solutions to incompressible two-phase flow. Journal of Computational Physics, 114(1), 146-159. The reinitialisation equation psi_tau = S(psi_0)(1 - \|grad psi\|) with the smoothed sign S(psi_0) = psi_0/sqrt(psi_0^2 + dx^2). [`doi:10.1006/jcph.1994.1155`](https://doi.org/10.1006/jcph.1994.1155)
- Munoz-Esparza, D., Kosovic, B., Jimenez, P. A., & Coen, J. L. (2018). An accurate fire-spread algorithm in the Weather Research and Forecasting model using the level-set method. Journal of Advances in Modeling Earth Systems, 10(4), 908-926. The hybrid WENO5/ENO1 level-set scheme (fire_upwinding = 9) of WRF-Fire / NCAR fire_behavior. [`doi:10.1002/2017MS001108`](https://doi.org/10.1002/2017MS001108)
- Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228. The classical WENO5 smoothness indicators beta_k and nonlinear weights w_k = gamma_k/(beta_k + eps)^2. [`doi:10.1006/jcph.1996.0130`](https://doi.org/10.1006/jcph.1996.0130)
- Harten, A., Engquist, B., Osher, S., & Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. Journal of Computational Physics, 71(2), 231-303. ENO stencil selection by smallest divided difference (first-order form: the minmod limiter). [`doi:10.1016/0021-9991(87)90031-3`](https://doi.org/10.1016/0021-9991%2887%2990031-3)
- NCAR Community Fire Behavior Model (fire_behavior), develop @ eb77580, physics/level_set_mod.F90: Advance_ls_reinit (fire_upwinding_reinit = 4), Select_eno, Select_weno5, Select_4th, Extrapol_var_at_bdys. The transcribed reference implementation.

**Conformance fixtures**

- ast: [`weno5_eno1_signed_norm_D2_extrapolate_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/weno5_eno1_signed_norm_D2_extrapolate_bc) — status `active`

