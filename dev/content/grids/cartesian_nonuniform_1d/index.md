---
title: "cartesian_nonuniform_1d"
description: "Non-uniform 1-D Cartesian grid with consumer-supplied edge geometry: N cell-centered points whose spacing varies along x."
esd_kind: "grid"
families: ["cartesian", "finite_difference"]
ops: ["D", "D2"]
bcs: ["dirichlet", "neumann", "zero_flux", "zero_gradient"]
tags: ["esd:grid", "family:cartesian", "ndim:1", "spacing:nonuniform", "axis:x", "esd:stencil", "family:finite_difference", "grid:cartesian_nonuniform_1d", "op:D", "order:2", "stencil_width:3", "op:D2", "esd:rule", "axes:x", "bc:dirichlet", "spec:9.6.8", "bc:zero_gradient", "bc:neumann", "bc:zero_flux"]
source: "grids/cartesian_nonuniform_1d/grid.esm"
---

Non-uniform 1-D Cartesian grid with consumer-supplied edge geometry: N cell-centered points whose spacing varies along x. Two iteration axes: `x` (the N cells) and `x_nodes` (the N+1 cell edges/nodes). N (the cell count) is the esm-spec 9.7 metaparameter, bound at the import edge or the loader API; the node axis is sized N+1 by the metaparameter expression {op:+, args:[N, 1]} (esm-schema IndexSet.size MAY be a metaparameter expression folded at load, esm-spec 9.7.6). The geometry is NOT baked in: a consuming model supplies a SINGLE data array `xe` -- the node (edge) coordinates, shape [x_nodes], with xe[1] and xe[N+1] the domain ends -- and the GRID derives everything else from it through two match-less templates, so the mesh derivation lives in the grid (like grids/mpas derives cell rings from vertex arrays) rather than being duplicated per consumer: `nonuniform_cell_center` gives the cell centers xc[i] = (xe[i] + xe[i+1])/2 and `nonuniform_cell_width` gives the finite-volume control widths dxw[i] = xe[i+1] - xe[i]. A consuming model exposes those as its x-shaped `xc` and `dxw` observeds by APPLYING the templates to its `xe`, and the finite-volume rules then read xc/dxw by bare name (the keyed-factor contract grids/mpas uses for areaCell/dvEdge). Because centers and widths both descend from ONE edge array they are automatically mutually consistent -- there is no way to supply a center that disagrees with a width. The map itself is the consumer's: for a smooth stretch x(xi) on uniform computational nodes xi_k = (k-1)/N (k = 1..N+1) the consumer sets xe[k] = x((k-1)/N). On a smoothly stretched mesh the conservative finite-volume Laplacian built on this geometry is second-order in the L2 (supraconvergence) norm even though the pointwise truncation error is lower order near the walls (LeVeque 2007 section 2.13). A loader-API rebinding of N re-sizes both axes and the consumer's `xe` aggregate together, so a convergence sweep re-derives the whole mesh consistently; the reals-are-consumer-supplied build-time scope caveat (AGENTS.md) still applies to any ic/6.6.5 reference, which must fold the cell coordinate from literals and the metaparameter name N inline.

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `N` | integer | `64` | Number of cells along x (the node axis x_nodes has N+1 members). |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>dxw</code> <code>gL</code> <code>gR</code> <code>qL</code> <code>qR</code> <code>xc</code></p>

The geometry is consumer-supplied so the same files serve any **domain extent**, not just any resolution. For a domain [a, b] the consuming model defines `x0 = a` (a parameter) and `dx = (b − a)/N` (an observed whose expression divides by the metaparameter name `N`, so a loader-API rebinding of `N` — a convergence sweep — keeps `dx` consistent). A model that instead closes `N` at the import edge spells the matching literal (`{op: /, args: [1, 8]}` for N = 8 on the unit interval). See the extent proof `problems/heat_1d_zero_grad_nonunit.esm` (x ∈ [−1.5, 2.5], observed order 2.00).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `x` | interval | `N` |
| `x_nodes` | interval | `{'op': '+', 'args': ['N', 1]}` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `nonuniform_cell_center` | `xe` | `sum[gi] (0.5 * (xe[gi] + xe[gi + 1])) where {gi in x}` | Cell centers from the node array: xc[gi] = (xe[gi] + xe[gi+1]) / 2 over the x (cell) axis, gathering the N+1-node edge array xe (shape [x_nodes]). A consuming model applies this to its `xe` to define its `xc` observed. The aggregation index is named `gi`, NOT `i`, so the derived observed does not capture a consuming finite-volume rule's own aggregate index `i` when it is inlined into the rule body (the rule reads xc by variable index xc[i-1], xc[i], xc[i+1]). |
| `nonuniform_cell_width` | `xe` | `sum[gi] (xe[gi + 1] - xe[gi]) where {gi in x}` | Finite-volume control widths from the node array: dxw[gi] = xe[gi+1] - xe[gi] over the x (cell) axis. A consuming model applies this to its `xe` to define its `dxw` observed. Aggregation index `gi` (distinct from a rule's `i`) for the same hygiene reason as nonuniform_cell_center. |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `central_D1_nonuniform_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_nonuniform_1d/stencils/central_D1_nonuniform_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/stencils/central_D1_nonuniform_interior.esm)

Parameters: `f`

<div class="math-block">sum[i] ((-(xc[i + 1] - xc[i])) / ((xc[i] - xc[i - 1]) * (xc[i] - xc[i - 1] + (xc[i + 1] - xc[i]))) * f[i - 1] + (xc[i + 1] - xc[i] - (xc[i] - xc[i - 1])) / ((xc[i] - xc[i - 1]) * (xc[i + 1] - xc[i])) * f[i] + (xc[i] - xc[i - 1]) / ((xc[i + 1] - xc[i]) * (xc[i] - xc[i - 1] + (xc[i + 1] - xc[i]))) * f[i + 1]) where {i in 2:N - 1}</div>

Interior-only 3-point centered first derivative on cartesian_nonuniform_1d: D f[i] = -b/(a(a+b)) f[i-1] + (b-a)/(a b) f[i] + a/(b(a+b)) f[i+1] over i in [2, N-1], with a = xc[i]-xc[i-1] and b = xc[i+1]-xc[i] the two center-to-center spacings -- the exact Fornberg (1988) weights for the 1st derivative on the three unequally spaced cell centers xc[i-1], xc[i], xc[i+1]. Match-less named template; boundary-condition layers (rules/) wrap it into complete rewrite rules. `xc` (cell centers) is a free-name array resolving to the consuming model's geometry observed (the grid's consumer-supplied geometry contract; see grids/cartesian_nonuniform_1d/grid.esm) -- a first derivative needs only the centers, not the control widths dxw. On a smoothly stretched mesh this is second-order accurate in L2 (supraconvergent) even though the local truncation error is O(h) where the spacing varies; it reduces to the 3-point (f[i+1]-f[i-1])/(2 dx) when the mesh is uniform (a = b = dx).

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids) and the supraconvergence discussion.

### `central_D2_nonuniform_interior`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cartesian_nonuniform_1d/stencils/central_D2_nonuniform_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/stencils/central_D2_nonuniform_interior.esm)

Parameters: `f`

<div class="math-block">sum[i] (((f[i + 1] - f[i]) / (xc[i + 1] - xc[i]) - (f[i] - f[i - 1]) / (xc[i] - xc[i - 1])) / dxw[i]) where {i in 2:N - 1}</div>

Interior-only conservative finite-volume second derivative on cartesian_nonuniform_1d: L f[i] = ( (f[i+1]-f[i])/(xc[i+1]-xc[i]) - (f[i]-f[i-1])/(xc[i]-xc[i-1]) ) / dxw[i] over i in [2, N-1] -- the flux difference of the two face gradients divided by the control-cell width. Match-less named template; boundary-condition layers (rules/) wrap it into complete rewrite rules. `xc` (cell centers) and `dxw` (cell widths) are free-name arrays resolving to the consuming model's geometry observeds (the grid's consumer-supplied geometry contract; see grids/cartesian_nonuniform_1d/grid.esm). On a smoothly stretched mesh this is second-order accurate in L2 (supraconvergent) even though the local truncation error is O(h) near the walls; it reduces to the 3-point (f[i+1]-2f[i]+f[i-1])/dx^2 when the mesh is uniform (xc[i+1]-xc[i] = xc[i]-xc[i-1] = dxw[i] = dx).

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `central_D1_nonuniform_dirichlet_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_nonuniform_1d/rules/central_D1_nonuniform_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/rules/central_D1_nonuniform_dirichlet_bc.esm)

**Rewrites** `D(f)/Dx` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N - 1] = central_D1_nonuniform_interior&lt;f=f&gt;, [1:1] = (-(xc[2] - xc[1])) / (dxw[1] / 2 * (dxw[1] / 2 + (xc[2] - xc[1]))) * gL + (xc[2] - xc[1] - dxw[1] / 2) / (dxw[1] / 2 * (xc[2] - xc[1])) * f[1] + dxw[1] / 2 / ((xc[2] - xc[1]) * (dxw[1] / 2 + (xc[2] - xc[1]))) * f[2], [N:N] = (-dxw[N] / 2) / ((xc[N] - xc[N - 1]) * (xc[N] - xc[N - 1] + dxw[N] / 2)) * f[N - 1] + (dxw[N] / 2 - (xc[N] - xc[N - 1])) / ((xc[N] - xc[N - 1]) * (dxw[N] / 2)) * f[N] + (xc[N] - xc[N - 1]) / (dxw[N] / 2 * (xc[N] - xc[N - 1] + dxw[N] / 2)) * gR)</div>

Second-order centered first derivative with inhomogeneous DIRICHLET (fixed-value wall) boundary closures on cartesian_nonuniform_1d, cell-centered. Matches plain D(f, wrt: x) at DEFAULT priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x]}}` (esm-spec 9.6.1) so it fires only on the first derivative of a bare field over this grid's `x` index set; under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axis together. Interior [2, N-1]: the imported 3-point centered nonuniform stencil central_D1_nonuniform_interior. The boundary faces use the exact second-order 3-point (Fornberg 1988) first derivative evaluated at the boundary CELL CENTER, built on the prescribed WALL VALUE plus the two nearest cell centers. At i=1 the stencil is on the three points {xe[1] carrying value gL, xc[1] carrying f[1], xc[2] carrying f[2]} with left spacing p = xc[1]-xe[1] = dxw[1]/2 (the distance from the boundary cell center to the wall node, exactly half the boundary control width) and right spacing q = xc[2]-xc[1], giving D f[1] = -q/(p(p+q)) gL + (q-p)/(pq) f[1] + p/(q(p+q)) f[2]. Symmetrically at i=N on {xc[N-1] carrying f[N-1], xc[N] carrying f[N], xe[N+1] carrying gR} with left spacing p = xc[N]-xc[N-1] and right spacing q = xe[N+1]-xc[N] = dxw[N]/2, giving D f[N] = -q/(p(p+q)) f[N-1] + (q-p)/(pq) f[N] + p/(q(p+q)) gR. Because the boundary center xc[1] (resp. xc[N]) is the MIDDLE node of {xe[1], xc[1], xc[2]} (resp. {xc[N-1], xc[N], xe[N+1]}) the wall closure is a genuinely centered 3-point formula (exact for quadratics), hence second-order POINTWISE at the wall -- stronger than the mirror-ghost uniform closure grids/cartesian_uniform_1d/rules/central_D1_dirichlet_bc, which is only O(dx) locally. As the mesh uniformizes (dxw[i] -> dx, xc[i+1]-xc[i] -> dx, so p -> dx/2, q -> dx) each wall term reduces to (-4 gL + 3 f[1] + f[2])/(3 dx), the standard second-order first derivative at a point half a cell from the boundary given the wall value. `xc` (cell centers) and `dxw` (cell widths) are free-name arrays resolving to the consuming model's geometry observeds (grids/cartesian_nonuniform_1d/grid.esm); the first-derivative Dirichlet closure reads the half control-width dxw[1]/2, dxw[N]/2 as the center-to-wall distance (unlike the zero-gradient sibling central_D1_nonuniform_zero_grad_bc, which needs only centers). `gL`, `gR` are consumer-supplied free names giving the field value at the left wall x=xe[1] and the right wall x=xe[N+1]; a consuming model declares them as ordinary real parameters and a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so gL=gR=0 is the homogeneous-Dirichlet case. The gL/gR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Second-order in L2 on smoothly stretched meshes (supraconvergence); a Dirichlet-compatible manufactured mode -- e.g. sin(pi x), whose value vanishes at both walls (gL=gR=0) while its wall gradient +/- pi is nonzero, genuinely exercising the fixed-value closure rather than a zero-gradient one -- matches the walls. A model imports exactly one first-derivative rule for `x`.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Sections 2.12 (boundary conditions; cell-centered Dirichlet via one-sided wall reconstruction) and 2.13 (nonuniform grids).

**Conformance fixtures**

- ast: [`central_D1_nonuniform_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1_nonuniform_dirichlet_bc) — status `active`
- convergence: [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet) — status `active`
- simulation: [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_1d_nonuniform_dirichlet) — status `active`

#### Convergence — case [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/advection_1d_nonuniform_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 8.853e-04 | — |
| 32 | 2.223e-04 | 1.99 |
| 64 | 5.480e-05 | 2.02 |
| 128 | 1.355e-05 | 2.02 |

Expected order: **2** (± 0.25).

![Convergence of advection_1d_nonuniform_dirichlet (log-log)](/plots/convergence/advection_1d_nonuniform_dirichlet.png)

### `central_D1_nonuniform_zero_grad_bc`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:zero_gradient</code> <code>axes:x</code> <code>stencil width 3</code></p>

Source: [`grids/cartesian_nonuniform_1d/rules/central_D1_nonuniform_zero_grad_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/rules/central_D1_nonuniform_zero_grad_bc.esm)

**Rewrites** `D(f)/Dx` — match pattern:

```json
{"op": "D", "args": ["f"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N - 1] = central_D1_nonuniform_interior&lt;f=f&gt;, [1:1] = (-(2 * (xc[2] - xc[1]) + (xc[3] - xc[2]))) / ((xc[2] - xc[1]) * (xc[2] - xc[1] + (xc[3] - xc[2]))) * f[1] + (xc[2] - xc[1] + (xc[3] - xc[2])) / ((xc[2] - xc[1]) * (xc[3] - xc[2])) * f[2] + (-(xc[2] - xc[1])) / ((xc[3] - xc[2]) * (xc[2] - xc[1] + (xc[3] - xc[2]))) * f[3], [N:N] = (xc[N] - xc[N - 1]) / ((xc[N - 1] - xc[N - 2]) * (xc[N - 1] - xc[N - 2] + (xc[N] - xc[N - 1]))) * f[N - 2] + (-(xc[N - 1] - xc[N - 2] + (xc[N] - xc[N - 1]))) / ((xc[N - 1] - xc[N - 2]) * (xc[N] - xc[N - 1])) * f[N - 1] + (xc[N - 1] - xc[N - 2] + 2 * (xc[N] - xc[N - 1])) / ((xc[N - 1] - xc[N - 2] + (xc[N] - xc[N - 1])) * (xc[N] - xc[N - 1])) * f[N])</div>

Second-order centered first derivative with second-order one-sided (zero-gradient-family) boundary closures on cartesian_nonuniform_1d. Matches plain D(f, wrt: x) at DEFAULT priority, so compound rules (e.g. a second-derivative rule at priority 10) fire first, match-scoped with `where: {f: {shape: [x]}}` (esm-spec 9.6.1) so it fires only on the first derivative of a bare field over this grid's `x` index set; under import-edge renaming (esm-spec 9.7.7) the wrt literal and the where shape follow the renamed axis together. Interior [2, N-1]: the imported 3-point centered nonuniform stencil central_D1_nonuniform_interior. The boundary faces use the exact second-order ONE-SIDED nonuniform first derivative (Fornberg 1988) built on the three nearest cell centers -- at i=1 the forward stencil on {xc[1], xc[2], xc[3]}, at i=N the backward stencil on {xc[N-2], xc[N-1], xc[N]} -- so the operator is second-order at the walls without any ghost cell, the natural zero-gradient-compatible closure for the manufactured cos(pi x) mode (whose gradient vanishes at both domain ends). `xc` (cell centers) is a free-name array resolving to the consuming model's geometry observed (grids/cartesian_nonuniform_1d/grid.esm); the first derivative needs only centers, not widths. Second-order in L2 on smoothly stretched meshes (supraconvergence). A model imports exactly one first-derivative rule for `x`.

**References**

- Fornberg, B. (1988). Generation of finite difference formulas on arbitrarily spaced grids. Mathematics of Computation 51(184), 699-706. [`doi:10.1090/S0025-5718-1988-0935077-0`](https://doi.org/10.1090/S0025-5718-1988-0935077-0)
- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Sections 2.13 (nonuniform grids) and 2.12 (one-sided boundary closures).

**Conformance fixtures**

- ast: [`central_D1_nonuniform_zero_grad_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D1_nonuniform_zero_grad_bc) — status `active`
- convergence: [`advection_1d_nonuniform_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_zero_grad) — status `active`
- simulation: [`advection_1d_nonuniform_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_1d_nonuniform_zero_grad) — status `active`

#### Convergence — case [`advection_1d_nonuniform_zero_grad`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_zero_grad)

Error norms read from the committed golden [`tests/conformance/convergence/advection_1d_nonuniform_zero_grad/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_1d_nonuniform_zero_grad/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 2.250e-03 | — |
| 32 | 5.621e-04 | 2.00 |
| 64 | 1.443e-04 | 1.96 |
| 128 | 3.925e-05 | 1.88 |

Expected order: **2** (± 0.25).

![Convergence of advection_1d_nonuniform_zero_grad (log-log)](/plots/convergence/advection_1d_nonuniform_zero_grad.png)

### `central_D2_nonuniform_dirichlet_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:dirichlet</code> <code>axes:x</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_dirichlet_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_dirichlet_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N - 1] = central_D2_nonuniform_interior&lt;f=f&gt;, [1:1] = ((f[2] - f[1]) / (xc[2] - xc[1]) - (f[1] - gL) / (dxw[1] / 2)) / dxw[1], [N:N] = ((gR - f[N]) / (dxw[N] / 2) - (f[N] - f[N - 1]) / (xc[N] - xc[N - 1])) / dxw[N])</div>

Conservative finite-volume second derivative with inhomogeneous Dirichlet (fixed-value wall) boundaries on cartesian_nonuniform_1d, cell-centered. Compound match on D(D(f,x),x) at priority 10, so it fires on the whole second derivative before any plain-D rule lowers the inner derivative (esm-spec 9.6.3 compound precedence), match-scoped with `where: {f: {shape: [x]}}` (esm-spec 9.6.1); under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axis together. Interior [2, N-1]: the imported conservative flux-difference stencil. The boundary faces impose a prescribed value at the domain wall (the finite-volume Dirichlet closure -- no ghost geometry): the left wall sits at node xe[1], a half control-width dxw[1]/2 below the first center xc[1], and the right wall sits at node xe[N+1], a half control-width dxw[N]/2 above the last center xc[N]. Replacing the missing outside face-gradient by the one-sided gradient to the wall value gives L f[1] = ( (f[2]-f[1])/(xc[2]-xc[1]) - (f[1]-gL)/(dxw[1]/2) ) / dxw[1] at i=1 and L f[N] = ( (gR-f[N])/(dxw[N]/2) - (f[N]-f[N-1])/(xc[N]-xc[N-1]) ) / dxw[N] at i=N. As the mesh uniformizes (xc[i+1]-xc[i] -> dx, dxw[i] -> dx) the wall term (f[1]-gL)/(dx/2) is the linear reconstruction to the half-cell-outside wall, recovering the cell-centered mirror-ghost Dirichlet face (grids/cartesian_uniform_1d/rules/central_D2_dirichlet_bc); at gL=gR=0 it is the homogeneous-Dirichlet closure. `xc` (cell centers) and `dxw` (cell widths) are free-name arrays resolving to the consuming model's geometry observeds (grids/cartesian_nonuniform_1d/grid.esm), and `gL`, `gR` are consumer-supplied free names giving the field value at the left wall x=xe[1] and the right wall x=xe[N+1]; a consuming model declares gL/gR as ordinary real parameters and a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so gL=gR=0 is the homogeneous case. The gL/gR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Second-order in L2 on smoothly stretched meshes (supraconvergence); a Dirichlet-compatible manufactured mode (e.g. sin(pi x), zero at both walls) matches the closure so no forcing is needed. A model imports exactly one second-derivative rule for `x`.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Sections 2.12 (boundary conditions; cell-centered Dirichlet via one-sided wall reconstruction) and 2.13 (nonuniform grids).

**Conformance fixtures**

- ast: [`central_D2_nonuniform_dirichlet_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2_nonuniform_dirichlet_bc) — status `active`
- convergence: [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet) — status `active`
- convergence: [`heat_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform_dirichlet) — status `active`
- simulation: [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/advection_1d_nonuniform_dirichlet) — status `active`
- simulation: [`heat_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_1d_nonuniform_dirichlet) — status `active`

#### Convergence — case [`advection_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/advection_1d_nonuniform_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/advection_1d_nonuniform_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 8.853e-04 | — |
| 32 | 2.223e-04 | 1.99 |
| 64 | 5.480e-05 | 2.02 |
| 128 | 1.355e-05 | 2.02 |

Expected order: **2** (± 0.25).

![Convergence of advection_1d_nonuniform_dirichlet (log-log)](/plots/convergence/advection_1d_nonuniform_dirichlet.png)

#### Convergence — case [`heat_1d_nonuniform_dirichlet`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform_dirichlet)

Error norms read from the committed golden [`tests/conformance/convergence/heat_1d_nonuniform_dirichlet/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_1d_nonuniform_dirichlet/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 5.049e-05 | — |
| 32 | 1.265e-05 | 2.00 |
| 64 | 3.165e-06 | 2.00 |
| 128 | 7.913e-07 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of heat_1d_nonuniform_dirichlet (log-log)](/plots/convergence/heat_1d_nonuniform_dirichlet.png)

### `central_D2_nonuniform_neumann_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:neumann</code> <code>axes:x</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_neumann_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_neumann_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N - 1] = central_D2_nonuniform_interior&lt;f=f&gt;, [1:1] = ((f[2] - f[1]) / (xc[2] - xc[1]) - qL) / dxw[1], [N:N] = (qR - (f[N] - f[N - 1]) / (xc[N] - xc[N - 1])) / dxw[N])</div>

Conservative finite-volume second derivative with inhomogeneous NEUMANN (fixed-flux wall) boundaries on cartesian_nonuniform_1d, cell-centered. Compound match on D(D(f,x),x) at priority 10, so it fires on the whole second derivative before any plain-D rule lowers the inner derivative (esm-spec 9.6.3 compound precedence), match-scoped with `where: {f: {shape: [x]}}` (esm-spec 9.6.1); under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axis together. Interior [2, N-1]: the imported conservative flux-difference stencil central_D2_nonuniform_interior, L f[i] = ( (f[i+1]-f[i])/(xc[i+1]-xc[i]) - (f[i]-f[i-1])/(xc[i]-xc[i-1]) ) / dxw[i], the difference of the two face gradients divided by the control-cell width. The boundary faces prescribe the PHYSICAL WALL FLUX directly: q = df/dx (in the +x direction) is fixed at the domain wall, so the missing outside face-gradient is REPLACED by the prescribed value rather than dropped. In the conservative form L f[i] = (F_right - F_left)/dxw[i] the wall face flux F is set to the prescribed q -- at i=1 the left face sits at node xe[1] and its flux is qL, giving L f[1] = ( (f[2]-f[1])/(xc[2]-xc[1]) - qL ) / dxw[1]; at i=N the right face sits at node xe[N+1] and its flux is qR, giving L f[N] = ( qR - (f[N]-f[N-1])/(xc[N]-xc[N-1]) ) / dxw[N]. Because q enters as a bare additive/subtractive flux over the boundary control width dxw (NOT scaled by any interior center-to-center distance), this rule reduces EXACTLY to grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_zero_flux_bc when qL=qR=0: the qL term drops leaving ((f[2]-f[1])/(xc[2]-xc[1]))/dxw[1], and the qR term drops leaving -( (f[N]-f[N-1])/(xc[N]-xc[N-1]) )/dxw[N] = ( (f[N-1]-f[N])/(xc[N]-xc[N-1]) )/dxw[N], the zero-flux faces byte-for-byte in value. As the mesh uniformizes (xc[i+1]-xc[i] -> dx, dxw[i] -> dx) L f[1] -> (f[2]-f[1]-qL*dx)/dx^2 and L f[N] -> (f[N-1]-f[N]+qR*dx)/dx^2, recovering the uniform ghost-elimination Neumann face grids/cartesian_uniform_1d/rules/central_D2_neumann_bc. `xc` (cell centers) and `dxw` (cell widths) are free-name arrays resolving to the consuming model's geometry observeds (grids/cartesian_nonuniform_1d/grid.esm), and `qL`, `qR` are consumer-supplied free names giving the value of df/dx (in the +x direction) at the left wall x=xe[1] and the right wall x=xe[N+1]; a consuming model declares qL/qR as ordinary real parameters and a parameter not overridden falls back to its default (esm-spec 6: default, or 0), so qL=qR=0 is the homogeneous (zero-flux) case. The qL/qR terms live only in the runtime rule body, never in build-time ic/reference positions (AGENTS.md build-time scope caveat). Second-order in L2 on smoothly stretched meshes (supraconvergence). A Neumann-compatible manufactured mode with CONSTANT-in-time wall flux -- e.g. u(x,t) = qL*x + exp(-alpha pi^2 t) cos(pi x), whose gradient at both walls is the constant qL=qR while the decaying cos part carries zero wall gradient -- exercises the fixed-flux closure exactly (the linear ramp qL*x is annihilated to machine precision by the discrete operator because its discrete face gradients equal qL at every face including the two walls, so the wall term cancels it, leaving only the zero-flux evolution of the cos mode). A model imports exactly one second-derivative rule for `x`.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Sections 2.12 (Neumann boundary conditions by prescribing the wall flux) and 2.13 (nonuniform grids).

**Conformance fixtures**

- ast: [`central_D2_nonuniform_neumann_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2_nonuniform_neumann_bc) — status `active`
- convergence: [`heat_1d_nonuniform_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform_neumann) — status `active`
- simulation: [`heat_1d_nonuniform_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_1d_nonuniform_neumann) — status `active`

#### Convergence — case [`heat_1d_nonuniform_neumann`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform_neumann)

Error norms read from the committed golden [`tests/conformance/convergence/heat_1d_nonuniform_neumann/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_1d_nonuniform_neumann/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 7.811e-05 | — |
| 32 | 2.104e-05 | 1.89 |
| 64 | 5.403e-06 | 1.96 |
| 128 | 1.361e-06 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of heat_1d_nonuniform_neumann (log-log)](/plots/convergence/heat_1d_nonuniform_neumann.png)

### `central_D2_nonuniform_zero_flux_bc`

<p class="rule-meta"><code>op:D2</code> <code>order:2</code> <code>bc:zero_flux</code> <code>axes:x</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_zero_flux_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cartesian_nonuniform_1d/rules/central_D2_nonuniform_zero_flux_bc.esm)

**Rewrites** `D(D(f)/Dx)/Dx` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "D", "args": ["f"], "wrt": "x"}], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[x]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N - 1] = central_D2_nonuniform_interior&lt;f=f&gt;, [1:1] = (f[2] - f[1]) / (xc[2] - xc[1]) / dxw[1], [N:N] = (f[N - 1] - f[N]) / (xc[N] - xc[N - 1]) / dxw[N])</div>

Conservative finite-volume second derivative with zero-flux (homogeneous Neumann) boundaries on cartesian_nonuniform_1d. Compound match on D(D(f,x),x) at priority 10, so it fires on the whole second derivative before any plain-D rule lowers the inner derivative (esm-spec 9.6.3 compound precedence), match-scoped with `where: {f: {shape: [x]}}` (esm-spec 9.6.1); under import-edge renaming (esm-spec 9.7.7) the wrt literals and the where shape follow the renamed axis together. Interior [2, N-1]: the imported conservative flux-difference stencil. The boundary faces impose zero flux exactly (natural for the finite-volume form -- no ghost geometry): at i=1 the left domain-face flux is dropped, giving L f[1] = ((f[2]-f[1])/(xc[2]-xc[1])) / dxw[1]; at i=N the right domain-face flux is dropped, giving L f[N] = ((f[N-1]-f[N])/(xc[N]-xc[N-1])) / dxw[N]. `xc` (cell centers) and `dxw` (cell widths) are free-name arrays resolving to the consuming model's geometry observeds (grids/cartesian_nonuniform_1d/grid.esm). Second-order in L2 on smoothly stretched meshes (supraconvergence); the manufactured cos(pi x) mode has zero gradient at both domain ends, matching the zero-flux closure so no forcing is needed. A model imports exactly one second-derivative rule for `x`.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Sections 2.12 (Neumann boundary conditions) and 2.13 (nonuniform grids).

**Conformance fixtures**

- ast: [`central_D2_nonuniform_zero_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/central_D2_nonuniform_zero_flux_bc) — status `active`
- convergence: [`heat_1d_nonuniform`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform) — status `active`
- simulation: [`heat_1d_nonuniform`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_1d_nonuniform) — status `active`

#### Convergence — case [`heat_1d_nonuniform`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_1d_nonuniform)

Error norms read from the committed golden [`tests/conformance/convergence/heat_1d_nonuniform/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_1d_nonuniform/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 8 | 2.681e-04 | — |
| 16 | 7.863e-05 | 1.77 |
| 32 | 2.118e-05 | 1.89 |
| 64 | 5.439e-06 | 1.96 |
| 128 | 1.370e-06 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of heat_1d_nonuniform (log-log)](/plots/convergence/heat_1d_nonuniform.png)

