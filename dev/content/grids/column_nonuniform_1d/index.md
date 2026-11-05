---
title: "column_nonuniform_1d"
description: "Single atmospheric (or oceanic/soil) COLUMN with consumer-supplied edge geometry: NLEV cell-centered layers stacked along the vertical index `lev`, plus the NLEV+1 layer interfaces `lev_nodes`."
esd_kind: "grid"
families: ["column", "finite_volume"]
ops: ["D", "advection", "flux", "integral", "interface_flux", "interpolate", "limiter", "reconstruct", "varcoeff_laplacian"]
bcs: ["no-inflow", "none", "prescribed_flux", "supplied_flux"]
tags: ["esd:grid", "family:column", "ndim:1", "spacing:nonuniform", "axis:lev", "esd:stencil", "family:finite_volume", "grid:column_nonuniform_1d", "op:D", "order:2", "stencil_width:2", "op:integral", "stencil_width:1", "op:interface_flux", "op:advection", "order:1", "scheme:upwind", "op:flux", "axes:lev", "stencil_width:3", "op:interpolate", "stencil_width:4", "op:reconstruct", "order:3", "stencil_width:5", "scheme:semi_lagrangian_plm_limit", "op:limiter", "op:varcoeff_laplacian", "esd:rule", "bc:supplied_flux", "spec:9.6.8", "bc:none", "axes:lev_nodes", "axes:none", "bc:no-inflow", "form:flux", "bc:prescribed_flux"]
source: "grids/column_nonuniform_1d/grid.esm"
---

Single atmospheric (or oceanic/soil) COLUMN with consumer-supplied edge geometry: NLEV cell-centered layers stacked along the vertical index `lev`, plus the NLEV+1 layer interfaces `lev_nodes`. Layer 1 is the BOTTOM (surface) layer and `lev` increases upward, following the WRF/CCPP convention (k=1 lowest model level); a model whose vertical index runs top-down reverses its arrays at the boundary of the component. NLEV (the layer count) is the esm-spec 9.7 metaparameter, bound at the import edge or the loader API; the interface axis is sized by the metaparameter expression {op:+, args:[NLEV, 1]} (esm-spec 9.7.6). The geometry is NOT baked in: a consuming model supplies a SINGLE data array `ze` -- the interface coordinates, shape [lev_nodes], monotone increasing with the lev index, ze[1] the surface and ze[NLEV+1] the model top -- and the GRID derives everything else from it. The two geometry templates are not re-implemented here: this grid imports grids/cartesian_nonuniform_1d/grid.esm at an esm-spec 9.7.7 renaming edge (x -> lev, x_nodes -> lev_nodes, with the cartesian cell count N closed by the symbolic 9.7.6 edge binding N = NLEV, so the renamed index sets are deep-equal to the ones declared here and dedupe), and `column_cell_center` / `column_cell_thickness` are one-line applications of the imported `nonuniform_cell_center` / `nonuniform_cell_width` to the interface array. The layer mid-points zc[k] = (ze[k] + ze[k+1])/2 and the layer thicknesses dz[k] = ze[k+1] - ze[k] therefore have exactly one definition in the library. A consuming model exposes those as its lev-shaped `zc` and `dz` observeds by APPLYING the two templates to its `ze`, and every rule on this grid then reads `zc`, `dz` (and, for the integral rules, only `dz`) by bare name -- the keyed-factor contract grids/mpas uses for areaCell/dvEdge. The one-sided and whole-column integral rules additionally name the two DOMAIN-FACE coordinates by the free names `z_sfc` and `z_top` (rules/integral_lev_whole): a consumer declares them as the observeds ze[1] and ze[NLEV+1] (or as parameters) or rebinds them at its import edge (esm-spec 9.7.7). Because centers and thicknesses both descend from ONE interface array they are automatically mutually consistent. The coordinate need not be geometric height: a pressure- or mass-based column supplies its own monotone interface coordinate as `ze` (or rebinds `dz` at the import edge to a mass-per-area layer measure such as dp/g), and the finite-volume operators and Riemann-sum integrals are then taken with respect to that coordinate. This grid exists separately from grids/cartesian_nonuniform_1d -- rather than as an import-edge rename of it in every consumer -- because column physics is the dominant use of a 1-D non-uniform mesh in Earth-system models and reads naturally in the lev/zc/dz/surface/top vocabulary shared with grids/latlon3d: the interface-coefficient diffusion rule (a PBL scheme's exchange coefficients), its prescribed-surface-flux closure, the flux-form vertical divergence and the surface-to-interface / interface-to-top cumulative integrals (optical depths, burdens) are column operators with no cartesian counterpart to rename. (Since EarthSciAST#174 esm-spec 9.7.7 does rewrite an `integral` node's `var` and bare-axis-name bounds under a rename; the geometry import above is that same mechanism.) The aggregation index of the imported geometry templates is `gi`, NOT `k`, so a derived zc/dz observed does not capture a consuming rule's own aggregate index `k` when it is inlined (the non-hygienic-inlining footgun documented in grids/cartesian_nonuniform_1d). The minimum extent is NLEV = 2: the interior region [2, NLEV-1] of the three-point diffusion rule folds to the empty [2, 1] and the two face closures cover the column; at NLEV = 1 it folds to the inverted [2, 0] and the load is rejected with makearray_region_inverted (esm-spec 4.3.2). On a smoothly stretched column the conservative finite-volume operators built on this geometry are second-order in L2 (supraconvergence, LeVeque 2007 section 2.13), and the midpoint-rule integrals are second-order. The build-time scope caveat (AGENTS.md) applies to any ic / 6.6.5 reference, which must fold the layer coordinate from literals and the metaparameter name NLEV inline. The quadratic-stretch test column shared by this grid's MMS problems and AST fixtures is the fixture library mesh/mms_column.esm.

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NLEV` | integer | `64` | Number of layers (cells) along lev; the interface axis lev_nodes has NLEV+1 members. Minimum 2 (see the grid description). |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>dz</code> <code>kdudz_bot</code> <code>kdudz_top</code> <code>zc</code></p>

These names are supplied by the consuming model's scope at evaluation (the keyed-factor / free-name geometry contract); they are neither metaparameters (load-time integers) nor template parameters (bound per apply site).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `lev` | interval | `NLEV` |
| `lev_nodes` | interval | `{'op': '+', 'args': ['NLEV', 1]}` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `column_cell_center` | `ze` | `nonuniform_cell_center<xe=ze>` | Layer mid-points from the interface array: zc[gi] = (ze[gi] + ze[gi+1]) / 2 over the lev (cell) axis, gathering the NLEV+1-node interface array ze (shape [lev_nodes]). The cartesian_nonuniform_1d `nonuniform_cell_center` template applied on the renamed lev axis -- one definition in the library. A consuming model applies this to its `ze` to define its `zc` observed. Aggregation index `gi` (not `k`) so the derived observed never captures a consuming rule's own aggregate index when inlined. |
| `column_cell_thickness` | `ze` | `nonuniform_cell_width<xe=ze>` | Layer thicknesses from the interface array: dz[gi] = ze[gi+1] - ze[gi] over the lev (cell) axis (positive when ze increases with the lev index, the required orientation). The cartesian_nonuniform_1d `nonuniform_cell_width` template applied on the renamed lev axis. A consuming model applies this to its `ze` to define its `dz` observed, the finite-volume control measure read by every rule on this grid. |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `face_flux_D_lev_interior`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/face_flux_D_lev_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/face_flux_D_lev_interior.esm)

Parameters: `F`

<div class="math-block">sum[k] ((F[k + 1] - F[k]) / dz[k]) where {k in 1:NLEV}</div>

Finite-volume divergence of an INTERFACE (lev_nodes-shaped) flux on column_nonuniform_1d: (d/dz F)[k] = (F[k+1] - F[k]) / dz[k] over every layer k in [1, NLEV]. Because the flux is supplied at all NLEV+1 interfaces -- including the surface F[1] and the model top F[NLEV+1] -- the operator is complete on the whole axis and has no separate boundary layer; the consuming model's flux field carries its own boundary values (a rigid lid is F[NLEV+1] = 0, a surface flux is F[1]). Match-less named template; rules/face_flux_D_lev_supplied_faces wraps it into the complete rewrite rule. `dz` is the free-name layer-thickness array resolving to the consuming model's geometry observed (grids/column_nonuniform_1d/grid.esm). Exact telescoping: the sum over k of dz[k] (d/dz F)[k] equals F[NLEV+1] - F[1], so column totals are conserved to round-off whatever F is.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

### `integral_lev_prefix_exclusive`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/integral_lev_prefix_exclusive.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/integral_lev_prefix_exclusive.esm)

Parameters: `f`

<div class="math-block">sum[k] (f[j] * dz[j]) where {j in lev, k in lev} if j &lt; k</div>

EXCLUSIVE prefix (surface-upward) Riemann sum of a layer field on column_nonuniform_1d: P[k] = sum over j < k of f[j] dz[j], a lev-shaped field with P[1] = 0 (the empty prefix). P[k] is the midpoint-rule integral from the surface ze[1] to the interface ze[k] BELOW layer k -- whole layers only, so it is exact for the layer sum. It is spelled in the esm-spec 4.3.1 prefix-reduction shape -- both ranges drawn from the ONE index set lev, monotone filter j < k -- which is the form the spec licenses a binding to evaluate as a single O(NLEV) running accumulator (a forward scan is bit-identical to the triangular double loop). rules/integral_lev_cumulative_from_bottom adds the half-layer f[k] dz[k]/2 to reach the layer centre; rules/integral_lev_nodes_cumulative_from_bottom places P at interfaces 1..NLEV and the column total at node NLEV+1. Match-less named template. `dz` is the consumer's free-name layer-thickness observed (grids/column_nonuniform_1d/grid.esm).

**References**

- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

### `integral_lev_suffix_inclusive`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/integral_lev_suffix_inclusive.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/integral_lev_suffix_inclusive.esm)

Parameters: `f`

<div class="math-block">sum[k] (f[j] * dz[j]) where {j in lev, k in lev} if j &gt;= k</div>

INCLUSIVE suffix (top-downward) Riemann sum of a layer field on column_nonuniform_1d: S[k] = sum over j >= k of f[j] dz[j], a lev-shaped field with S[NLEV] = f[NLEV] dz[NLEV]. S[k] is the midpoint-rule integral from the interface ze[k] BELOW layer k to the model top ze[NLEV+1] -- whole layers only, exact for the layer sum; this is the optical depth measured down from the top of the atmosphere to the base of layer k. Spelled in the esm-spec 4.3.1 prefix-reduction shape (both ranges from the one index set lev, monotone filter j >= k); note that 4.3.1 licenses the O(N) running-accumulator evaluation for FORWARD scans only, so a binding evaluates this reverse scan as the triangular double loop (each suffix folded from its own low end, which is what keeps the result bit-identical across bindings). rules/integral_lev_cumulative_to_top adds the half-layer f[k] dz[k]/2 to start from the layer centre; rules/integral_lev_nodes_cumulative_to_top places S at interfaces 1..NLEV and 0 at the top node NLEV+1. Match-less named template. `dz` is the consumer's free-name layer-thickness observed (grids/column_nonuniform_1d/grid.esm).

**References**

- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

### `integral_lev_total`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/integral_lev_total.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/integral_lev_total.esm)

Parameters: `f`

<div class="math-block">sum[] (f[j] * dz[j]) where {j in lev}</div>

Whole-column composite midpoint sum of a layer field on column_nonuniform_1d: sum over j in lev of f[j] dz[j], a SCALAR (a full reduction, output_idx []). The one definition of the column total in the library: rules/integral_lev_whole wraps it as the whole-column `integral` lowering, and rules/integral_lev_nodes_cumulative_from_bottom places it at the top interface (node NLEV+1) of the surface-to-interface cumulative integral. Match-less named template. `dz` is the consumer's free-name layer-thickness observed (grids/column_nonuniform_1d/grid.esm); because the layer values are midpoint samples the sum is second-order accurate on a smoothly stretched column.

**References**

- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

### `interface_flux_K_lev`

<p class="rule-meta"><code>op:interface_flux</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/interface_flux_K_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/interface_flux_K_lev.esm)

Parameters: `u`, `K`, `n`

<div class="math-block">K[n] * ((u[n] - u[n - 1]) / (zc[n] - zc[n - 1]))</div>

Diffusive flux K du/dz through ONE interior interface n (2 <= n <= NLEV) of column_nonuniform_1d: K[n] (u[n] - u[n-1]) / (zc[n] - zc[n-1]), the interface-coefficient gradient between the two layers the interface separates (node n lies between layers n-1 and n; node 1 is the surface and node NLEV+1 the model top, neither of which has a computable gradient). The interface index `n` is a template PARAMETER bound to an index expression at each call site -- the enclosing aggregate's loop symbol k or k+1 in stencils/varcoeff_face_laplacian_lev_interior, the literal 2 and the metaparameter NLEV in the two face closures of rules/varcoeff_face_laplacian_lev_flux_bc -- so the flux formula has ONE definition in the library and the telescoping identity sum_k dz[k] L u[k] = kdudz_top - kdudz_bot is structural: every interior flux enters the column sum exactly once with each sign. Positive in the +lev (upward) direction. Match-less named template. `zc` is the consumer's free-name layer-centre observed (grids/column_nonuniform_1d/grid.esm).

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

### `interface_total_flux_lev`

<p class="rule-meta"><code>op:interface_flux</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/interface_total_flux_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/interface_total_flux_lev.esm)

Parameters: `u`, `K`, `F`, `n`

<div class="math-block">interface_flux_K_lev&lt;K=K, n=n, u=u&gt; - F[n]</div>

Total flux K du/dz - F through ONE interior interface n (2 <= n <= NLEV) of column_nonuniform_1d: the diffusive flux of stencils/interface_flux_K_lev at n minus the prescribed interface flux F[n], with K and F both on the interfaces (shape [lev_nodes]; node n lies between layers n-1 and n, node 1 is the surface and node NLEV+1 the model top). The interface index `n` is a template PARAMETER bound to an index expression at each call site -- k+1 and k in stencils/varcoeff_face_flux_laplacian_lev_interior, the literal 2 and the metaparameter NLEV in the two face closures of rules/varcoeff_face_flux_laplacian_lev -- so the total-flux formula has ONE definition and the telescoping identity sum_k dz[k] L u[k] = F[1] - F[NLEV+1] is structural. F is the prescribed UPWARD (+lev) flux and -K du/dz the upward turbulent flux, so this quantity is minus the total upward flux through the interface. Match-less named template.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

### `sedimentation_donor_flux_lev`

<p class="rule-meta"><code>op:advection</code> <code>order:1</code> <code>stencil width 1</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_donor_flux_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_donor_flux_lev.esm)

Parameters: `q`, `W`, `n`

<div class="math-block">W[n] * q[n]</div>

Donor-cell sedimentation flux through the LOWER face of ONE layer n (1 <= n <= NLEV) of column_nonuniform_1d: F[n] = W[n] q[n], the layer's own falling density q carried at the layer's own signed velocity W along +lev (W = -v_t <= 0 for a terminal fall speed v_t; both shape [lev]). Face n is the lower edge of layer n; because the flow is downward everywhere the upstream (donor) cell of that face is layer n itself, so the upwind flux and the microphysics notion that mass leaves a layer at that layer's mass-weighted terminal velocity coincide. The face index `n` is a template PARAMETER so the same flux serves the interior difference (n = k+1 and n = k) and the model-top closure (n = NLEV) of rules/sedimentation_upwind1_flux_D_lev, with every face flux one reference to this template and the divergence telescoping exactly (adjacent layers reuse the identical face expression). Face NLEV+1 (the model top) has no flux: nothing falls in from above the domain. Valid only for W <= 0: a sign-changing velocity moves the donor cell to the other side of the face and needs a face-velocity upwind rule instead.

**References**

- LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge University Press. Section 4.1 (the upwind / donor-cell flux).
- Hong, S.-Y., & Lim, J.-O. J. (2006). The WRF single-moment 6-class microphysics scheme (WSM6). Journal of the Korean Meteorological Society, 42(2), 129-151.

### `sedimentation_jh2010_end_face_flux`

<p class="rule-meta"><code>op:flux</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_end_face_flux.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_end_face_flux.esm)

Parameters: `Wm1`, `W0`, `qm1`, `q0`, `qp1`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">sedimentation_jh2010_end_face_velocity&lt;W0=W0, Wm1=Wm1&gt; * sedimentation_jh2010_face_value&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;</div>

FLUX through an END face of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit: F = W_face q_face through face 2 (the lower face of layer 2) or face NLEV (the lower face of the top layer), the two faces whose velocity is the two-point sedimentation_jh2010_end_face_velocity(Wm1, W0). Because that velocity can never turn upward, the donor is always the layer ABOVE the face and q_face is its lower-face value sedimentation_jh2010_face_value(qm1, q0, qp1, dzm1, dz0, dzp1), with q0 the donor layer, qm1 the layer below the face and qp1 the layer above the donor. Applied at face 2 by the rule's layer-1 and layer-2 regions (with q[1], q[2], q[3]) and, through sedimentation_jh2010_lid_face_flux, at face NLEV where the layer above the donor is the lid. The interior faces 3 <= j <= NLEV-1 use sedimentation_jh2010_face_flux instead, which carries the upwind-side switch.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)

### `sedimentation_jh2010_end_face_velocity`

<p class="rule-meta"><code>op:interpolate</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_end_face_velocity.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_end_face_velocity.esm)

Parameters: `Wm1`, `W0`

<div class="math-block">ifelse(W0 == 0, Wm1, 0.5 * (W0 + Wm1))</div>

END-FACE fall velocity of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, at the two faces where the four-point interior interpolation (sedimentation_jh2010_face_velocity) has no room: face 2 (the lower face of layer 2) and face NLEV (the lower face of the top layer). There WRF uses the TWO-POINT average, wi(2) = 0.5 (ww(2) + ww(1)) and wi(km) = 0.5 (ww(km) + ww(km-1)), so W_face = 0.5 (W0 + Wm1) with W0 the layer above the face and Wm1 the layer below, for the layer-centred signed velocity W = -v_t <= 0 along +lev. WRF's 'terminate of top of raingroup' override applies here as at every face above the first: where the layer above the face carries no falling speed at all, W0 == 0, the face takes the speed of the layer below it, W_face = Wm1. That branch is part of the scheme and is kept as a branch. Because both branches are averages or copies of non-positive speeds this face velocity can never turn upward, so the end-face fluxes need no upwind-side switch (sedimentation_jh2010_end_face_flux). Face 1 takes the bottom layer's own velocity and is written in the rule with the donor-flux template.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: wi(2) = 0.5*(ww(2)+ww(1)), wi(km) = 0.5*(ww(km)+ww(km-1)) and the 'terminate of top of raingroup' override.

### `sedimentation_jh2010_face_flux`

<p class="rule-meta"><code>op:flux</code> <code>order:2</code> <code>stencil width 4</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_flux.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_flux.esm)

Parameters: `Wm2`, `Wm1`, `W0`, `Wp1`, `qm2`, `qm1`, `q0`, `qp1`, `dzm2`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">sedimentation_jh2010_face_velocity&lt;W0=W0, Wm1=Wm1, Wm2=Wm2, Wp1=Wp1&gt; * ifelse(sedimentation_jh2010_face_velocity&lt;W0=W0, Wm1=Wm1, Wm2=Wm2, Wp1=Wp1&gt; &gt; 0, sedimentation_jh2010_top_face_value&lt;dz0=dzm1, dzm1=dzm2, dzp1=dz0, q0=qm1, qm1=qm2, qp1=q0&gt;, sedimentation_jh2010_face_value&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;)</div>

FLUX through an INTERIOR face of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit: F = W_face q_face through the lower face of layer k, for 3 <= k <= NLEV-1, the faces whose velocity is sedimentation_jh2010_face_velocity. The face value is taken from the UPWIND side of the face. For W_face <= 0 (downward, the normal case) that is the LOWER-face value of layer k, sedimentation_jh2010_face_value(q[k-1], q[k], q[k+1]); for W_face > 0 it is the UPPER-face value of layer k-1, sedimentation_jh2010_top_face_value(q[k-2], q[k-1], q[k]). The second branch is not hypothetical: although every layer-centred speed is downward (W <= 0), the third-order interpolation 9/16 (W0 + Wm1) - 1/16 (Wp1 + Wm2) overshoots to an UPWARD face speed wherever a slow layer sits under a fast one -- at the edge of a falling layer, 3-12 faces per column in EqWeFiC's WSM6 ice columns. In WRF the departure point of such a face, za(k) = zi(k) - wi(k) dt, lies ABOVE the face, so the remap draws the mass crossing it from the top of the layer below (qpi(k-1)), and the dt -> 0 flux is wi(k) qpi(k-1). Written without this switch the rule leaves 5.2e-4 of the peak tendency against WRF on a synthetic column with such a face (4.4e-10 with it), and 5.6e-7 to 6.4e-6 of the peak ice tendency on the four EqWeFiC WSM6 sedimentation columns (1.0e-8 to 3.9e-8 with it). The other interior-face ingredients -- the velocity's W0 == 0 override and the two limiter switches behind both face values (sedimentation_jh2010_half_jump) -- live in the stencils this one applies. The face velocity is applied twice, once for the sign test and once as the factor, because a template body has no local binding. Only interior faces need the switch: face 1 carries W[1] <= 0, and faces 2 and NLEV are two-point averages of non-positive speeds, which cannot turn upward. q[k-2] exists because k >= 3.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: the departure points za(k) = zi(k) - wi(k) dt and the 'interpolation to regular point' remap.

### `sedimentation_jh2010_face_value`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_value.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_value.esm)

Parameters: `qm1`, `q0`, `qp1`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">q0 - sedimentation_jh2010_half_jump&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;</div>

LOWER-FACE VALUE of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit: the value at the LOWER face of a layer of falling density q0, reconstructed from the layer and its two neighbours (qm1 below, qp1 above) on a non-uniform column (thicknesses dzm1, dz0, dzp1). It is q0 - h with h the limited half-jump sedimentation_jh2010_half_jump (the average of the two one-sided slopes carried half a layer, zeroed at an extremum and dropped where a face value would go negative), exactly WRF's qmi(k) = 2 qa(k) - qpi(k). For downward flow the lower face of a layer is fed by that layer, so this is the DONOR value of the face below it, and it is the value every face flux uses whenever the face speed is downward (sedimentation_jh2010_face_flux, sedimentation_jh2010_end_face_flux). Its twin sedimentation_jh2010_top_face_value = q0 + h is the upper-face value (WRF's qpi).

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: qmi(k) = 2.0*qa(k) - qpi(k).

### `sedimentation_jh2010_face_velocity`

<p class="rule-meta"><code>op:interpolate</code> <code>order:3</code> <code>stencil width 4</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_velocity.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_face_velocity.esm)

Parameters: `Wm2`, `Wm1`, `W0`, `Wp1`

<div class="math-block">ifelse(W0 == 0, Wm1, 0.5625 * (W0 + Wm1) - 0.0625 * (Wp1 + Wm2))</div>

INTERFACE fall velocity of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, at an INTERIOR face: the lower face of layer k, for 3 <= k <= NLEV-1. The layer-centred velocity W (signed along +lev, W = -v_t <= 0) is carried to the face by the THIRD-ORDER four-point interpolation W_face = 9/16 (W0 + Wm1) - 1/16 (Wp1 + Wm2), with W0 = W[k] (the layer above the face), Wm1 = W[k-1] (the layer below), Wp1 = W[k+1] and Wm2 = W[k-2]. The interpolation weights are those of a uniform grid and WRF applies them unchanged on its stretched eta column, so they are NOT thickness-weighted here either. One exception, WRF's "terminate of top of raingroup": where the layer above the face carries no falling speed at all, W0 == 0, the face takes the speed of the layer BELOW it, W_face = Wm1. That branch is part of the scheme, not a limiter added to it, and it is kept as a branch. Because the interpolation is linear, writing it for W = -v_t is the same as writing WRF's for v_t with the sign flipped. The two end faces (face 2 and face NLEV) use WRF's two-point averages and face 1 the bottom layer's own velocity; they are written in the rule's boundary regions, not here.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: the third-order wi interpolation and the 'terminate of top of raingroup' override.

### `sedimentation_jh2010_flux_D_lev_interior`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>stencil width 5</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_flux_D_lev_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_flux_D_lev_interior.esm)

Parameters: `q`, `W`

<div class="math-block">sum[k] ((sedimentation_jh2010_face_flux&lt;W0=W[k + 1], Wm1=W[k], Wm2=W[k - 1], Wp1=W[k + 2], dz0=dz[k + 1], dzm1=dz[k], dzm2=dz[k - 1], dzp1=dz[k + 2], q0=q[k + 1], qm1=q[k], qm2=q[k - 1], qp1=q[k + 2]&gt; - sedimentation_jh2010_face_flux&lt;W0=W[k], Wm1=W[k - 1], Wm2=W[k - 2], Wp1=W[k + 1], dz0=dz[k], dzm1=dz[k - 1], dzm2=dz[k - 2], dzp1=dz[k + 1], q0=q[k], qm1=q[k - 1], qm2=q[k - 2], qp1=q[k + 1]&gt;) / dz[k]) where {k in 3:NLEV - 2}</div>

INTERIOR layers k = 3 .. NLEV-2 of the Juang & Hong (2010) sedimentation divergence on column_nonuniform_1d, in the dt -> 0 limit of WRF WSM6's forward semi-Lagrangian fallout (nislfv_rain_plm): (d/dz (W q))[k] = (F[k+1] - F[k]) / dz[k], where F[j] is the flux through the LOWER face of layer j, F[j] = W_face[j] q_face[j], given by sedimentation_jh2010_face_flux. The face velocity is sedimentation_jh2010_face_velocity (third-order, four layers W[j-2..j+1]); the face value is taken from the UPWIND side, which for a downward face speed is the lower-face value of layer j (sedimentation_jh2010_face_value, the monotone average-slope reconstruction over q[j-1..j+1]) and, where the third-order speed overshoots upward, the upper-face value of layer j-1 (over q[j-2..j]). Layer k therefore reads W, q and dz over k-2 .. k+2 (five layers); both faces of every layer in this range are interior faces (3 <= j <= NLEV-1), which is what bounds the range. Same signs and the same velocity-inside-the-operator convention as sedimentation_upwind1_flux_D_lev_interior: q is the layer-centred falling density, W = -v_t <= 0 the layer-centred signed velocity along +lev, and a consuming model writes dq/dt = -D(W*q, wrt: lev).

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)

### `sedimentation_jh2010_half_jump`

<p class="rule-meta"><code>op:limiter</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_half_jump.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_half_jump.esm)

Parameters: `qm1`, `q0`, `qp1`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">ifelse(sedimentation_jh2010_one_sided_slope&lt;dza=dzp1, dzb=dz0, qa=qp1, qb=q0&gt; * sedimentation_jh2010_one_sided_slope&lt;dza=dz0, dzb=dzm1, qa=q0, qb=qm1&gt; &gt; 0 and abs(sedimentation_jh2010_unlimited_half_jump&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;) &lt;= q0, sedimentation_jh2010_unlimited_half_jump&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;, 0)</div>

LIMITED half-jump of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit: the rise h of the monotone piecewise-linear reconstruction from the centre of a layer of value q0 to either of its faces, so that the lower-face value is q0 - h (sedimentation_jh2010_face_value, WRF's qmi) and the upper-face value q0 + h (sedimentation_jh2010_top_face_value, WRF's qpi). h is the unlimited average-slope half-jump sedimentation_jh2010_unlimited_half_jump, or 0 where either of the scheme's two switches fires. Both switches are kept as branches because they ARE the scheme: (1) at an extremum of q, where the two one-sided slopes dip and dim (sedimentation_jh2010_one_sided_slope) satisfy dip dim <= 0, the slope is zeroed (WRF: qmi(k) = qpi(k) = qa(k)); (2) if either face value q0 - h or q0 + h would be negative the slope is also dropped (WRF: if qpi(k) < 0 or qmi(k) < 0). Branch (2) is written as |h| <= q0, which is the same test exactly: q0 + h and q0 - h are both >= 0 iff |h| <= q0, and the subtraction is exact near the boundary by Sterbenz's lemma. The template keeps the slope only when NEITHER switch fires, so the condition is dip dim > 0 AND |h| <= q0. THIS IS NOT A TEXTBOOK LIMITER: the slope is the average of the one-sided differences, not their minmod and not the monotonized-central limit, and on a measured WRF column a minmod reconstruction in its place leaves 31.6 % of the difference to WRF that this closes. Do not 'simplify' it to a library limiter. The bottom layer takes no reconstruction in WRF (qmi(1) = qpi(1) = qa(1)) and the rule never applies this template there; at the model top the rule binds qp1 = 0 and dzp1 = 0 (WRF's qa(km+1) = 0, dza(km+1) = wi(km+1) dt -> 0).

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: 'estimate values at arrival cell interface with monotone' (the dip*dim <= 0 and qpi/qmi < 0 branches).

### `sedimentation_jh2010_lid_face_flux`

<p class="rule-meta"><code>op:flux</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_lid_face_flux.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_lid_face_flux.esm)

Parameters: `Wm1`, `W0`, `qm1`, `q0`, `dzm1`, `dz0`

<div class="math-block">sedimentation_jh2010_end_face_flux&lt;W0=W0, Wm1=Wm1, dz0=dz0, dzm1=dzm1, dzp1=0, q0=q0, qm1=qm1, qp1=0&gt;</div>

FLUX through face NLEV, the lower face of the TOP layer, of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit. It is the end-face flux sedimentation_jh2010_end_face_flux with the LID CLOSURE bound in one place: the layer above the donor (the top layer, q0 = q[NLEV], dz0 = dz[NLEV]) is the model lid, where WRF sets qa(km+1) = 0 and whose arrival thickness dza(km+1) = wi(km+1) dt vanishes as dt -> 0, so qp1 = 0 and dzp1 = 0. Wm1 = W[NLEV-1] and W0 = W[NLEV] give the two-point face velocity, qm1 = q[NLEV-1] and dzm1 = dz[NLEV-1] the layer below the face. The rule applies it twice, as the inflow of layer NLEV-1 and as the only flux of layer NLEV (there is no inflow through the lid, F[NLEV+1] = 0), so the two regions share one spelling of the closure.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: qa(km+1) = 0.0 and dza(km+1) = zi(km+1) - za(km+1).

### `sedimentation_jh2010_one_sided_slope`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:1</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_one_sided_slope.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_one_sided_slope.esm)

Parameters: `qa`, `qb`, `dza`, `dzb`

<div class="math-block">(qa - qb) / (dza + dzb)</div>

ONE-SIDED SLOPE of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d: the difference of two adjacent layer values over the SUM of their thicknesses, (qa - qb) / (dza + dzb), with a the layer above and b the layer below. WRF forms both one-sided differences of its monotone reconstruction this way, dip = (qa(k+1) - qa(k)) / (dza(k+1) + dza(k)) and dim = (qa(k) - qa(k-1)) / (dza(k-1) + dza(k)), and in the dt -> 0 limit the arrival thicknesses dza are the layer thicknesses. Note the denominator is the sum of the two thicknesses, not the distance between the two layer centres (half that sum): the factor 1/2 is carried by sedimentation_jh2010_unlimited_half_jump, which multiplies the sum of the two slopes by dz0 / 2. Applied twice per layer, once for each neighbour.

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: 'estimate values at arrival cell interface with monotone' (dip, dim).

### `sedimentation_jh2010_top_face_value`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_top_face_value.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_top_face_value.esm)

Parameters: `qm1`, `q0`, `qp1`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">q0 + sedimentation_jh2010_half_jump&lt;dz0=dz0, dzm1=dzm1, dzp1=dzp1, q0=q0, qm1=qm1, qp1=qp1&gt;</div>

UPPER-FACE VALUE of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d, in its dt -> 0 limit: the value at the UPPER face of a layer of falling density q0, q0 + h with h the limited half-jump sedimentation_jh2010_half_jump over the layer and its two neighbours (qm1 below, qp1 above; thicknesses dzm1, dz0, dzp1), exactly WRF's qpi(k) = qa(k) + 0.5 (dip + dim) dza(k) after its switches. The same reconstruction as sedimentation_jh2010_face_value = q0 - h, read at the other face. It is needed because the donor of a face is its UPWIND side: wherever WRF's third-order interface speed overshoots to an upward face speed, the mass crossing the face comes up from the top of the layer BELOW the face, and the flux is W_face times this value of that layer (sedimentation_jh2010_face_flux).

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: qpi(k) = qa(k) + 0.5*(dip+dim)*dza(k).

### `sedimentation_jh2010_unlimited_half_jump`

<p class="rule-meta"><code>op:reconstruct</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_jh2010_unlimited_half_jump.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_jh2010_unlimited_half_jump.esm)

Parameters: `qm1`, `q0`, `qp1`, `dzm1`, `dz0`, `dzp1`

<div class="math-block">0.5 * dz0 * (sedimentation_jh2010_one_sided_slope&lt;dza=dzp1, dzb=dz0, qa=qp1, qb=q0&gt; + sedimentation_jh2010_one_sided_slope&lt;dza=dz0, dzb=dzm1, qa=q0, qb=qm1&gt;)</div>

UNLIMITED half-jump of the Juang & Hong (2010) forward semi-Lagrangian fallout (WRF WSM6 nislfv_rain_plm / nislfv_rain_plm6) on column_nonuniform_1d: the rise of the piecewise-linear reconstruction from the centre of a layer of value q0 to either of its faces, before the scheme's two switches are applied. It is h = (dip + dim) dz0 / 2, the AVERAGE of the two one-sided slopes (sedimentation_jh2010_one_sided_slope, dip towards the layer above with value qp1 and thickness dzp1, dim towards the layer below with qm1 and dzm1) carried half a layer, exactly WRF's qpi(k) - qa(k) = 0.5 (dip + dim) dza(k). The scheme's switches (zero the slope at an extremum, drop it if either face value would go negative) are applied to this quantity by sedimentation_jh2010_half_jump; this template is the value they choose between and is not a limiter in itself. THIS IS NOT A TEXTBOOK SLOPE: the average of the one-sided differences over non-uniform thicknesses is neither the minmod nor the monotonized-central choice, and replacing it with either leaves a measurable difference to WRF (see sedimentation_jh2010_flux_D_lev).

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutine nislfv_rain_plm: qpi(k) = qa(k) + 0.5 (dip + dim) dza(k).

### `sedimentation_upwind1_flux_D_lev_interior`

<p class="rule-meta"><code>op:advection</code> <code>order:1</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/sedimentation_upwind1_flux_D_lev_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/sedimentation_upwind1_flux_D_lev_interior.esm)

Parameters: `q`, `W`

<div class="math-block">sum[k] ((sedimentation_donor_flux_lev&lt;W=W, n=k + 1, q=q&gt; - sedimentation_donor_flux_lev&lt;W=W, n=k, q=q&gt;) / dz[k]) where {k in 1:NLEV - 1}</div>

INTERIOR layers k = 1 .. NLEV-1 of the conservative first-order DONOR-CELL divergence of a downward hydrometeor flux on column_nonuniform_1d: (d/dz (W q))[k] = ( F[k+1] - F[k] ) / dz[k], where F[n] = W[n] q[n] is the donor flux of stencils/sedimentation_donor_flux_lev through the lower face of layer n, applied here at the face above (n = k+1) and below (n = k) layer k. Both the transported density q (rho times a mixing ratio, or any layer quantity) and the velocity W are LAYER-centred ([lev]); W is the signed velocity along +lev (upward), so a fall speed v_t enters as W = -v_t <= 0. Because the flow is downward everywhere, the upstream (donor) cell of face k -- the lower edge of layer k -- is layer k itself, so the mass leaving a layer falls at that layer's OWN speed: exactly the structure of a sedimentation scheme, in which the mass leaving a layer falls at its mass-weighted terminal velocity. Adjacent layers reuse the identical face expression, so the divergence TELESCOPES: sum_k dz[k] (d/dz (W q))[k] = F[NLEV+1] - F[1] = -F[1], the surface outflow, and the column mass budget closes to round-off (the surface precipitation rate is -F[1] = v_t[1] q[1]). Layer NLEV, whose upper face NLEV+1 carries no flux, is closed by rules/sedimentation_upwind1_flux_D_lev itself. Valid only for W <= 0: for a velocity that changes sign the donor cell changes side and a face-velocity upwind rule (latlon3d upwind1_flux_D_lev_noflux_bc) is needed instead. `dz` is the free-name layer thickness of grids/column_nonuniform_1d/grid.esm. First order (O(dz)), monotone and positivity-preserving.

**References**

- LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge University Press. Section 4.1 (the upwind / donor-cell flux).
- Hong, S.-Y., & Lim, J.-O. J. (2006). The WRF single-moment 6-class microphysics scheme (WSM6). Journal of the Korean Meteorological Society, 42(2), 129-151.
- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)

### `varcoeff_face_flux_laplacian_lev_interior`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/varcoeff_face_flux_laplacian_lev_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/varcoeff_face_flux_laplacian_lev_interior.esm)

Parameters: `u`, `K`, `F`

<div class="math-block">sum[k] ((interface_total_flux_lev&lt;F=F, K=K, n=k + 1, u=u&gt; - interface_total_flux_lev&lt;F=F, K=K, n=k, u=u&gt;) / dz[k]) where {k in 2:NLEV - 1}</div>

Interior-only conservative finite-volume divergence of the TOTAL interface flux K du/dz - F on column_nonuniform_1d: L u[k] = ( H[k+1] - H[k] ) / dz[k] over k in [2, NLEV-1], where H[n] = K[n] (u[n]-u[n-1])/(zc[n]-zc[n-1]) - F[n] is the total flux of stencils/interface_total_flux_lev (itself stencils/interface_flux_K_lev minus F[n]), applied at the interface above (n = k+1) and below (n = k) each layer, with the coefficient K and the prescribed flux F both supplied AT THE INTERFACES (shape [lev_nodes]; node 1 is the surface, node k separates layers k-1 and k, node NLEV+1 is the model top). It is stencils/varcoeff_face_laplacian_lev_interior with F subtracted inside the same difference, so the diffusive and the prescribed parts of the flux telescope together. Match-less named template; rules/varcoeff_face_flux_laplacian_lev wraps it with the two boundary layers, which close the same total-flux template against the prescribed domain-face fluxes. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). Second-order in L2 on a smoothly stretched column (supraconvergent).

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

### `varcoeff_face_laplacian_lev_interior`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d/stencils/varcoeff_face_laplacian_lev_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/stencils/varcoeff_face_laplacian_lev_interior.esm)

Parameters: `u`, `K`

<div class="math-block">sum[k] ((interface_flux_K_lev&lt;K=K, n=k + 1, u=u&gt; - interface_flux_K_lev&lt;K=K, n=k, u=u&gt;) / dz[k]) where {k in 2:NLEV - 1}</div>

Interior-only conservative finite-volume variable-coefficient Laplacian d/dz(K du/dz) on column_nonuniform_1d with the coefficient K supplied AT THE INTERFACES (shape [lev_nodes]): L u[k] = ( G[k+1] - G[k] ) / dz[k] over k in [2, NLEV-1], where G[n] = K[n] (u[n]-u[n-1])/(zc[n]-zc[n-1]) is the interface flux of stencils/interface_flux_K_lev, applied here at the interface above (n = k+1) and below (n = k) each layer. Interface numbering: K[k] lives on the interface BELOW layer k (node 1 is the surface, node k separates layers k-1 and k, node NLEV+1 is the model top), so the flux through the interface above layer k uses K[k+1]. This is exactly how atmospheric PBL schemes define their exchange coefficients -- e.g. WRF/CCPP YSU's xkzh(k), 'at the interface between k and k+1', is K[k+1] here -- so no face averaging is performed or needed (contrast grids/cartesian_uniform_1d/stencils/varcoeff_laplacian_interior, which averages a CELL-centred coefficient to the faces). Match-less named template; rules/varcoeff_face_laplacian_lev_flux_bc wraps it with the two boundary-flux layers, which close the same interface flux against the prescribed domain-face fluxes. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). Second-order in L2 on a smoothly stretched column (supraconvergent); reduces to the classical three-point form on a uniform column with uniform K.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `face_flux_D_lev_supplied_faces`

<p class="rule-meta"><code>op:D</code> <code>order:2</code> <code>bc:supplied_flux</code> <code>axes:lev</code> <code>stencil width 2</code> <code>priority 0</code></p>

Source: [`grids/column_nonuniform_1d/rules/face_flux_D_lev_supplied_faces.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/face_flux_D_lev_supplied_faces.esm)

**Rewrites** `D(F)/Dlev` at priority 0 — match pattern:

```json
{"op": "D", "args": ["F"], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>F</code> is a bare field shaped <code>[lev_nodes]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLEV] = face_flux_D_lev_interior&lt;F=F&gt;)</div>

Vertical flux divergence D(F, wrt: lev) of an INTERFACE flux field F (shape [lev_nodes]) on column_nonuniform_1d, lowered to the conservative finite-volume difference (F[k+1] - F[k]) / dz[k] over all NLEV layers. This is the primitive that every flux-form column process reduces to: a consuming model writes its conservation law as du/dt = -D(F, wrt: lev) + ... with F its upward (+lev) flux at the NLEV+1 interfaces -- turbulent flux, sedimentation flux, radiative net flux, precipitation flux -- and supplies the two domain-face values inside F itself (F[1] the surface exchange, F[NLEV+1] = 0 for a rigid model top), so the rule carries no boundary regions of its own (bc:supplied_flux). Match-scoped with where {F: {shape: [lev_nodes]}} (esm-spec 9.6.1) so it fires only on a bare interface-shaped field of THIS grid and never on a layer-shaped field or a compound: a consumer differentiating a compound flux (e.g. w*q, K*D(u)) first binds it to a declared [lev_nodes] observed, or uses a compound-matching rule such as varcoeff_face_laplacian_lev_flux_bc. Priority 0 (a plain-D rule); compound rules on this grid fire first at priority 10 under esm-spec 9.6.3 outermost-first precedence. `dz` is the consumer's free-name geometry observed. Exactly telescoping, second-order on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

**Conformance fixtures**

- ast: [`face_flux_D_lev_supplied_faces`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/face_flux_D_lev_supplied_faces) — status `active`

### `integral_lev_cumulative_from_bottom`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>bc:none</code> <code>axes:lev</code> <code>stencil width 1</code></p>

Source: [`grids/column_nonuniform_1d/rules/integral_lev_cumulative_from_bottom.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/integral_lev_cumulative_from_bottom.esm)

**Rewrites** `integral(f, lev, z_sfc, lev)` — match pattern:

```json
{"op": "integral", "args": ["f"], "var": "lev", "lower": "z_sfc", "upper": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">integral_lev_prefix_exclusive&lt;f=f&gt; + 0.5 * f * dz</div>

Cumulative integral from the SURFACE to the centre of the current layer: lowers integral(f, var: lev, lower: z_sfc, upper: lev) -- the upper bound literally the layer axis name, esm-spec 4.2 cumulative mode -- to the layer field F[k] = P[k] + f[k] dz[k]/2, where P is the exclusive prefix sum of stencils/integral_lev_prefix_exclusive (the whole layers below, an O(NLEV) forward scan in every binding) and the half layer reaches the layer centre zc[k] from the interface below it: the midpoint-rule integral from ze[1] to zc[k]. F[1] = f[1] dz[1]/2 (the empty prefix is zero). BOUNDS ARE GROUND, never wildcards (esm-spec 9.6.1: a non-parameter string in a match binds only that exact bare name): `lower` and `upper` must each be spelled as one of the surface free name `z_sfc`, the top free name `z_top`, the layer axis name `lev` (the current layer centre, esm-spec 4.2 cumulative mode) or the interface axis name `lev_nodes` (the current interface). The five integral rules on this grid partition those spellings -- (z_sfc, z_top) the whole column (integral_lev_whole); (z_sfc, lev) and (lev, z_top) from the surface to, and from, each layer centre (integral_lev_cumulative_from_bottom / _to_top); (z_sfc, lev_nodes) and (lev_nodes, z_top) from the surface to, and from, each interface (integral_lev_nodes_cumulative_from_bottom / _to_top) -- so at most one rule matches any integral node and no priority is needed. An integral over any other bounds (a boundary-layer height, an interior interface, a numeric literal, or two axis names) matches nothing, survives lowering, and is reported by the pre-evaluation unlowered_operator gate instead of being lowered silently to a column sum. The consumer declares z_sfc = ze[1] and z_top = ze[NLEV+1] as observeds (or as parameters), or rebinds the two names at its import edge (esm-spec 9.7.7). esm-spec 4.2 leaves the bounds measure-free, so the bounds carry no numeric weight; the measure is the consumer's free-name layer-thickness observed `dz` (grids/column_nonuniform_1d/grid.esm) -- rebind it at the import edge to integrate against another monotone layer measure such as dp/g. The integrand `f` must be a bare layer-shaped (lev) field of THIS grid (where {f: {shape: [lev]}}); bind a compound integrand to a declared [lev] observed first. Because the layer values are midpoint samples every form is the composite midpoint rule, second-order accurate on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

**Conformance fixtures**

- ast: [`integral_lev_cumulative_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/integral_lev_cumulative_from_bottom) — status `active`
- convergence: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom) — status `active`
- convergence: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top) — status `active`
- convergence: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom) — status `active`
- convergence: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_from_bottom) — status `active`
- simulation: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_to_top) — status `active`
- simulation: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_from_bottom) — status `active`
- simulation: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_whole) — status `active`

#### Convergence — case [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.511e-03 | — |
| 32 | 3.775e-04 | 2.00 |
| 64 | 9.435e-05 | 2.00 |
| 128 | 2.359e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_from_bottom (log-log)](/plots/convergence/integral_column_center_from_bottom.png)

#### Convergence — case [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.386e-03 | — |
| 32 | 3.457e-04 | 2.00 |
| 64 | 8.637e-05 | 2.00 |
| 128 | 2.159e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_to_top (log-log)](/plots/convergence/integral_column_center_to_top.png)

#### Convergence — case [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.418e-04 | — |
| 32 | 8.288e-05 | 2.04 |
| 64 | 2.043e-05 | 2.02 |
| 128 | 5.073e-06 | 2.01 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_from_bottom (log-log)](/plots/convergence/integral_column_nodes_from_bottom.png)

#### Convergence — case [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.298e-03 | — |
| 32 | 3.331e-04 | 1.96 |
| 64 | 8.438e-05 | 1.98 |
| 128 | 2.124e-05 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_to_top (log-log)](/plots/convergence/integral_column_nodes_to_top.png)

### `integral_lev_cumulative_to_top`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>bc:none</code> <code>axes:lev</code> <code>stencil width 1</code></p>

Source: [`grids/column_nonuniform_1d/rules/integral_lev_cumulative_to_top.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/integral_lev_cumulative_to_top.esm)

**Rewrites** `integral(f, lev, lev, z_top)` — match pattern:

```json
{"op": "integral", "args": ["f"], "var": "lev", "lower": "lev", "upper": "z_top"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">integral_lev_suffix_inclusive&lt;f=f&gt; - 0.5 * f * dz</div>

Cumulative integral from the centre of the current layer to the MODEL TOP: lowers integral(f, var: lev, lower: lev, upper: z_top) -- the lower bound literally the layer axis name -- to the layer field G[k] = S[k] - f[k] dz[k]/2, where S is the inclusive suffix sum of stencils/integral_lev_suffix_inclusive (the whole layers from the interface below layer k to the top) and the half layer below the centre is removed; written out, G[k] = f[k] dz[k]/2 + sum_{j>k} f[j] dz[j], the midpoint-rule integral from zc[k] to ze[NLEV+1]. The inclusive suffix is shared with the interface form (integral_lev_nodes_cumulative_to_top) so the grid carries one reverse scan. G[NLEV] = f[NLEV] dz[NLEV]/2. This is the form of an optical depth measured downward from the top of the atmosphere to a layer centre. BOUNDS ARE GROUND, never wildcards (esm-spec 9.6.1: a non-parameter string in a match binds only that exact bare name): `lower` and `upper` must each be spelled as one of the surface free name `z_sfc`, the top free name `z_top`, the layer axis name `lev` (the current layer centre, esm-spec 4.2 cumulative mode) or the interface axis name `lev_nodes` (the current interface). The five integral rules on this grid partition those spellings -- (z_sfc, z_top) the whole column (integral_lev_whole); (z_sfc, lev) and (lev, z_top) from the surface to, and from, each layer centre (integral_lev_cumulative_from_bottom / _to_top); (z_sfc, lev_nodes) and (lev_nodes, z_top) from the surface to, and from, each interface (integral_lev_nodes_cumulative_from_bottom / _to_top) -- so at most one rule matches any integral node and no priority is needed. An integral over any other bounds (a boundary-layer height, an interior interface, a numeric literal, or two axis names) matches nothing, survives lowering, and is reported by the pre-evaluation unlowered_operator gate instead of being lowered silently to a column sum. The consumer declares z_sfc = ze[1] and z_top = ze[NLEV+1] as observeds (or as parameters), or rebinds the two names at its import edge (esm-spec 9.7.7). esm-spec 4.2 leaves the bounds measure-free, so the bounds carry no numeric weight; the measure is the consumer's free-name layer-thickness observed `dz` (grids/column_nonuniform_1d/grid.esm) -- rebind it at the import edge to integrate against another monotone layer measure such as dp/g. The integrand `f` must be a bare layer-shaped (lev) field of THIS grid (where {f: {shape: [lev]}}); bind a compound integrand to a declared [lev] observed first. Because the layer values are midpoint samples every form is the composite midpoint rule, second-order accurate on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

**Conformance fixtures**

- ast: [`integral_lev_cumulative_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/integral_lev_cumulative_to_top) — status `active`
- convergence: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom) — status `active`
- convergence: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top) — status `active`
- convergence: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom) — status `active`
- convergence: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_from_bottom) — status `active`
- simulation: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_to_top) — status `active`
- simulation: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_from_bottom) — status `active`
- simulation: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_whole) — status `active`

#### Convergence — case [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.511e-03 | — |
| 32 | 3.775e-04 | 2.00 |
| 64 | 9.435e-05 | 2.00 |
| 128 | 2.359e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_from_bottom (log-log)](/plots/convergence/integral_column_center_from_bottom.png)

#### Convergence — case [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.386e-03 | — |
| 32 | 3.457e-04 | 2.00 |
| 64 | 8.637e-05 | 2.00 |
| 128 | 2.159e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_to_top (log-log)](/plots/convergence/integral_column_center_to_top.png)

#### Convergence — case [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.418e-04 | — |
| 32 | 8.288e-05 | 2.04 |
| 64 | 2.043e-05 | 2.02 |
| 128 | 5.073e-06 | 2.01 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_from_bottom (log-log)](/plots/convergence/integral_column_nodes_from_bottom.png)

#### Convergence — case [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.298e-03 | — |
| 32 | 3.331e-04 | 1.96 |
| 64 | 8.438e-05 | 1.98 |
| 128 | 2.124e-05 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_to_top (log-log)](/plots/convergence/integral_column_nodes_to_top.png)

### `integral_lev_nodes_cumulative_from_bottom`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>bc:none</code> <code>axes:lev_nodes</code> <code>stencil width 1</code></p>

Source: [`grids/column_nonuniform_1d/rules/integral_lev_nodes_cumulative_from_bottom.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/integral_lev_nodes_cumulative_from_bottom.esm)

**Rewrites** `integral(f, lev, z_sfc, lev_nodes)` — match pattern:

```json
{"op": "integral", "args": ["f"], "var": "lev", "lower": "z_sfc", "upper": "lev_nodes"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLEV] = integral_lev_prefix_exclusive&lt;f=f&gt;, [NLEV + 1:NLEV + 1] = integral_lev_total&lt;f=f&gt;)</div>

Cumulative integral from the SURFACE to each INTERFACE: lowers integral(f, var: lev, lower: z_sfc, upper: lev_nodes) -- the upper bound literally the interface axis name, meaning 'up to the coordinate of the current interface' -- to the interface field F[kn] = sum_{j<kn} f[j] dz[j] over kn in lev_nodes (shape [lev_nodes]): a makearray placing the exclusive prefix scan of stencils/integral_lev_prefix_exclusive on interfaces 1..NLEV (P[k] is exactly the integral up to the interface below layer k; P[1] = 0 at the surface) and the column total of stencils/integral_lev_total at the top interface NLEV+1 -- the esm-spec 4.3.1 O(NLEV) forward-scan shape plus one scalar reduction, never a quadratic double sum over two different axes. Whole layers only, so every layer is integrated exactly by the midpoint rule and the sum telescopes with face_flux_D_lev_supplied_faces: D(F, wrt: lev) recovers f exactly. Use it for interface-valued burdens and for a flux that is the running integral of a layer source (a precipitation flux from layer conversion rates, a radiative flux from layer heating). BOUNDS ARE GROUND, never wildcards (esm-spec 9.6.1: a non-parameter string in a match binds only that exact bare name): `lower` and `upper` must each be spelled as one of the surface free name `z_sfc`, the top free name `z_top`, the layer axis name `lev` (the current layer centre, esm-spec 4.2 cumulative mode) or the interface axis name `lev_nodes` (the current interface). The five integral rules on this grid partition those spellings -- (z_sfc, z_top) the whole column (integral_lev_whole); (z_sfc, lev) and (lev, z_top) from the surface to, and from, each layer centre (integral_lev_cumulative_from_bottom / _to_top); (z_sfc, lev_nodes) and (lev_nodes, z_top) from the surface to, and from, each interface (integral_lev_nodes_cumulative_from_bottom / _to_top) -- so at most one rule matches any integral node and no priority is needed. An integral over any other bounds (a boundary-layer height, an interior interface, a numeric literal, or two axis names) matches nothing, survives lowering, and is reported by the pre-evaluation unlowered_operator gate instead of being lowered silently to a column sum. The consumer declares z_sfc = ze[1] and z_top = ze[NLEV+1] as observeds (or as parameters), or rebinds the two names at its import edge (esm-spec 9.7.7). esm-spec 4.2 leaves the bounds measure-free, so the bounds carry no numeric weight; the measure is the consumer's free-name layer-thickness observed `dz` (grids/column_nonuniform_1d/grid.esm) -- rebind it at the import edge to integrate against another monotone layer measure such as dp/g. The integrand `f` must be a bare layer-shaped (lev) field of THIS grid (where {f: {shape: [lev]}}); bind a compound integrand to a declared [lev] observed first. Because the layer values are midpoint samples every form is the composite midpoint rule, second-order accurate on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

**Conformance fixtures**

- ast: [`integral_lev_nodes_cumulative_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/integral_lev_nodes_cumulative_from_bottom) — status `active`
- convergence: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom) — status `active`
- convergence: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top) — status `active`
- convergence: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom) — status `active`
- convergence: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_from_bottom) — status `active`
- simulation: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_to_top) — status `active`
- simulation: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_from_bottom) — status `active`
- simulation: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_whole) — status `active`

#### Convergence — case [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.511e-03 | — |
| 32 | 3.775e-04 | 2.00 |
| 64 | 9.435e-05 | 2.00 |
| 128 | 2.359e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_from_bottom (log-log)](/plots/convergence/integral_column_center_from_bottom.png)

#### Convergence — case [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.386e-03 | — |
| 32 | 3.457e-04 | 2.00 |
| 64 | 8.637e-05 | 2.00 |
| 128 | 2.159e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_to_top (log-log)](/plots/convergence/integral_column_center_to_top.png)

#### Convergence — case [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.418e-04 | — |
| 32 | 8.288e-05 | 2.04 |
| 64 | 2.043e-05 | 2.02 |
| 128 | 5.073e-06 | 2.01 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_from_bottom (log-log)](/plots/convergence/integral_column_nodes_from_bottom.png)

#### Convergence — case [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.298e-03 | — |
| 32 | 3.331e-04 | 1.96 |
| 64 | 8.438e-05 | 1.98 |
| 128 | 2.124e-05 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_to_top (log-log)](/plots/convergence/integral_column_nodes_to_top.png)

### `integral_lev_nodes_cumulative_to_top`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>bc:none</code> <code>axes:lev_nodes</code> <code>stencil width 1</code></p>

Source: [`grids/column_nonuniform_1d/rules/integral_lev_nodes_cumulative_to_top.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/integral_lev_nodes_cumulative_to_top.esm)

**Rewrites** `integral(f, lev, lev_nodes, z_top)` — match pattern:

```json
{"op": "integral", "args": ["f"], "var": "lev", "lower": "lev_nodes", "upper": "z_top"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLEV] = integral_lev_suffix_inclusive&lt;f=f&gt;, [NLEV + 1:NLEV + 1] = 0)</div>

Cumulative integral from each INTERFACE to the MODEL TOP: lowers integral(f, var: lev, lower: lev_nodes, upper: z_top) -- the lower bound literally the interface axis name -- to the interface field G[kn] = sum_{j>=kn} f[j] dz[j] over kn in lev_nodes (shape [lev_nodes]): a makearray placing the inclusive suffix scan of stencils/integral_lev_suffix_inclusive on interfaces 1..NLEV (S[k] is exactly the integral from the interface below layer k to the top; S[1] is the whole column) and the literal 0 at the top interface NLEV+1. Whole layers only, exact telescoping with face_flux_D_lev_supplied_faces (D(G, wrt: lev) = -f). This is the optical depth from the top of the atmosphere down to each interface, the quantity a two-stream or Beer-Lambert shortwave scheme (e.g. WRF's Dudhia scheme) needs at the level interfaces where it evaluates fluxes. BOUNDS ARE GROUND, never wildcards (esm-spec 9.6.1: a non-parameter string in a match binds only that exact bare name): `lower` and `upper` must each be spelled as one of the surface free name `z_sfc`, the top free name `z_top`, the layer axis name `lev` (the current layer centre, esm-spec 4.2 cumulative mode) or the interface axis name `lev_nodes` (the current interface). The five integral rules on this grid partition those spellings -- (z_sfc, z_top) the whole column (integral_lev_whole); (z_sfc, lev) and (lev, z_top) from the surface to, and from, each layer centre (integral_lev_cumulative_from_bottom / _to_top); (z_sfc, lev_nodes) and (lev_nodes, z_top) from the surface to, and from, each interface (integral_lev_nodes_cumulative_from_bottom / _to_top) -- so at most one rule matches any integral node and no priority is needed. An integral over any other bounds (a boundary-layer height, an interior interface, a numeric literal, or two axis names) matches nothing, survives lowering, and is reported by the pre-evaluation unlowered_operator gate instead of being lowered silently to a column sum. The consumer declares z_sfc = ze[1] and z_top = ze[NLEV+1] as observeds (or as parameters), or rebinds the two names at its import edge (esm-spec 9.7.7). esm-spec 4.2 leaves the bounds measure-free, so the bounds carry no numeric weight; the measure is the consumer's free-name layer-thickness observed `dz` (grids/column_nonuniform_1d/grid.esm) -- rebind it at the import edge to integrate against another monotone layer measure such as dp/g. The integrand `f` must be a bare layer-shaped (lev) field of THIS grid (where {f: {shape: [lev]}}); bind a compound integrand to a declared [lev] observed first. Because the layer values are midpoint samples every form is the composite midpoint rule, second-order accurate on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

**Conformance fixtures**

- ast: [`integral_lev_nodes_cumulative_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/integral_lev_nodes_cumulative_to_top) — status `active`
- convergence: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom) — status `active`
- convergence: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top) — status `active`
- convergence: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom) — status `active`
- convergence: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_from_bottom) — status `active`
- simulation: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_to_top) — status `active`
- simulation: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_from_bottom) — status `active`
- simulation: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_whole) — status `active`

#### Convergence — case [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.511e-03 | — |
| 32 | 3.775e-04 | 2.00 |
| 64 | 9.435e-05 | 2.00 |
| 128 | 2.359e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_from_bottom (log-log)](/plots/convergence/integral_column_center_from_bottom.png)

#### Convergence — case [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.386e-03 | — |
| 32 | 3.457e-04 | 2.00 |
| 64 | 8.637e-05 | 2.00 |
| 128 | 2.159e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_to_top (log-log)](/plots/convergence/integral_column_center_to_top.png)

#### Convergence — case [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.418e-04 | — |
| 32 | 8.288e-05 | 2.04 |
| 64 | 2.043e-05 | 2.02 |
| 128 | 5.073e-06 | 2.01 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_from_bottom (log-log)](/plots/convergence/integral_column_nodes_from_bottom.png)

#### Convergence — case [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.298e-03 | — |
| 32 | 3.331e-04 | 1.96 |
| 64 | 8.438e-05 | 1.98 |
| 128 | 2.124e-05 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_to_top (log-log)](/plots/convergence/integral_column_nodes_to_top.png)

### `integral_lev_whole`

<p class="rule-meta"><code>op:integral</code> <code>order:2</code> <code>bc:none</code> <code>axes:none</code> <code>stencil width 1</code></p>

Source: [`grids/column_nonuniform_1d/rules/integral_lev_whole.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/integral_lev_whole.esm)

**Rewrites** `integral(f, lev, z_sfc, z_top)` — match pattern:

```json
{"op": "integral", "args": ["f"], "var": "lev", "lower": "z_sfc", "upper": "z_top"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>f</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">integral_lev_total&lt;f=f&gt;</div>

Whole-column integral: lowers integral(f, var: lev, lower: z_sfc, upper: z_top) -- the two bounds spelled as the surface and model-top free names -- to the composite midpoint sum over all NLEV layers, sum_j f[j] dz[j], a SCALAR (the stencil integral_lev_total; tag axes:none, the rank-0 output form of the L006 tiling lint). Typical uses: column burden, column-integrated heating, precipitable water, the total optical depth. BOUNDS ARE GROUND, never wildcards (esm-spec 9.6.1: a non-parameter string in a match binds only that exact bare name): `lower` and `upper` must each be spelled as one of the surface free name `z_sfc`, the top free name `z_top`, the layer axis name `lev` (the current layer centre, esm-spec 4.2 cumulative mode) or the interface axis name `lev_nodes` (the current interface). The five integral rules on this grid partition those spellings -- (z_sfc, z_top) the whole column (integral_lev_whole); (z_sfc, lev) and (lev, z_top) from the surface to, and from, each layer centre (integral_lev_cumulative_from_bottom / _to_top); (z_sfc, lev_nodes) and (lev_nodes, z_top) from the surface to, and from, each interface (integral_lev_nodes_cumulative_from_bottom / _to_top) -- so at most one rule matches any integral node and no priority is needed. An integral over any other bounds (a boundary-layer height, an interior interface, a numeric literal, or two axis names) matches nothing, survives lowering, and is reported by the pre-evaluation unlowered_operator gate instead of being lowered silently to a column sum. The consumer declares z_sfc = ze[1] and z_top = ze[NLEV+1] as observeds (or as parameters), or rebinds the two names at its import edge (esm-spec 9.7.7). esm-spec 4.2 leaves the bounds measure-free, so the bounds carry no numeric weight; the measure is the consumer's free-name layer-thickness observed `dz` (grids/column_nonuniform_1d/grid.esm) -- rebind it at the import edge to integrate against another monotone layer measure such as dp/g. The integrand `f` must be a bare layer-shaped (lev) field of THIS grid (where {f: {shape: [lev]}}); bind a compound integrand to a declared [lev] observed first. Because the layer values are midpoint samples every form is the composite midpoint rule, second-order accurate on a smoothly stretched column.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Davis, P. J. and Rabinowitz, P. (1984). Methods of Numerical Integration, 2nd ed. Academic Press. Section 2.1 (the midpoint rule and its error).

**Conformance fixtures**

- ast: [`integral_lev_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/integral_lev_whole) — status `active`
- convergence: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom) — status `active`
- convergence: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top) — status `active`
- convergence: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom) — status `active`
- convergence: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_from_bottom) — status `active`
- simulation: [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_center_to_top) — status `active`
- simulation: [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_from_bottom) — status `active`
- simulation: [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_nodes_to_top) — status `active`
- simulation: [`integral_column_whole`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/integral_column_whole) — status `active`

#### Convergence — case [`integral_column_center_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.511e-03 | — |
| 32 | 3.775e-04 | 2.00 |
| 64 | 9.435e-05 | 2.00 |
| 128 | 2.359e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_from_bottom (log-log)](/plots/convergence/integral_column_center_from_bottom.png)

#### Convergence — case [`integral_column_center_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_center_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_center_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_center_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.386e-03 | — |
| 32 | 3.457e-04 | 2.00 |
| 64 | 8.637e-05 | 2.00 |
| 128 | 2.159e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_center_to_top (log-log)](/plots/convergence/integral_column_center_to_top.png)

#### Convergence — case [`integral_column_nodes_from_bottom`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_from_bottom)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_from_bottom/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.418e-04 | — |
| 32 | 8.288e-05 | 2.04 |
| 64 | 2.043e-05 | 2.02 |
| 128 | 5.073e-06 | 2.01 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_from_bottom (log-log)](/plots/convergence/integral_column_nodes_from_bottom.png)

#### Convergence — case [`integral_column_nodes_to_top`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/integral_column_nodes_to_top)

Error norms read from the committed golden [`tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/integral_column_nodes_to_top/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 1.298e-03 | — |
| 32 | 3.331e-04 | 1.96 |
| 64 | 8.438e-05 | 1.98 |
| 128 | 2.124e-05 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of integral_column_nodes_to_top (log-log)](/plots/convergence/integral_column_nodes_to_top.png)

### `sedimentation_jh2010_flux_D_lev`

<p class="rule-meta"><code>op:advection</code> <code>order:1</code> <code>bc:no-inflow</code> <code>axes:lev</code> <code>stencil width 5</code> <code>priority 10</code></p>

Source: [`grids/column_nonuniform_1d/rules/sedimentation_jh2010_flux_D_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/sedimentation_jh2010_flux_D_lev.esm)

**Rewrites** `D(W * q)/Dlev` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["W", "q"]}], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>q</code> is a bare field shaped <code>[lev]</code>; <code>W</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:1] = sum[k] ((sedimentation_jh2010_end_face_flux&lt;W0=W[k + 1], Wm1=W[k], dz0=dz[k + 1], dzm1=dz[k], dzp1=dz[k + 2], q0=q[k + 1], qm1=q[k], qp1=q[k + 2]&gt; - sedimentation_donor_flux_lev&lt;W=W, n=k, q=q&gt;) / dz[k]) where {k in 1:1}, [2:2] = sum[k] ((sedimentation_jh2010_face_flux&lt;W0=W[k + 1], Wm1=W[k], Wm2=W[k - 1], Wp1=W[k + 2], dz0=dz[k + 1], dzm1=dz[k], dzm2=dz[k - 1], dzp1=dz[k + 2], q0=q[k + 1], qm1=q[k], qm2=q[k - 1], qp1=q[k + 2]&gt; - sedimentation_jh2010_end_face_flux&lt;W0=W[k], Wm1=W[k - 1], dz0=dz[k], dzm1=dz[k - 1], dzp1=dz[k + 1], q0=q[k], qm1=q[k - 1], qp1=q[k + 1]&gt;) / dz[k]) where {k in 2:2}, [3:NLEV - 2] = sedimentation_jh2010_flux_D_lev_interior&lt;W=W, q=q&gt;, [NLEV - 1:NLEV - 1] = sum[k] ((sedimentation_jh2010_lid_face_flux&lt;W0=W[k + 1], Wm1=W[k], dz0=dz[k + 1], dzm1=dz[k], q0=q[k + 1], qm1=q[k]&gt; - sedimentation_jh2010_face_flux&lt;W0=W[k], Wm1=W[k - 1], Wm2=W[k - 2], Wp1=W[k + 1], dz0=dz[k], dzm1=dz[k - 1], dzm2=dz[k - 2], dzp1=dz[k + 1], q0=q[k], qm1=q[k - 1], qm2=q[k - 2], qp1=q[k + 1]&gt;) / dz[k]) where {k in NLEV - 1:NLEV - 1}, [NLEV:NLEV] = sum[k] ((-sedimentation_jh2010_lid_face_flux&lt;W0=W[k], Wm1=W[k - 1], dz0=dz[k], dzm1=dz[k - 1], q0=q[k], qm1=q[k - 1]&gt;) / dz[k]) where {k in NLEV:NLEV})</div>

Conservative lowering of the sedimentation (fallout) flux divergence D(W*q, wrt: lev) of a DOWNWARD-falling quantity on column_nonuniform_1d that reproduces the dt -> 0 LIMIT of the Juang & Hong (2010) forward semi-Lagrangian scheme -- the fallout WRF's WSM6 microphysics applies to cloud ice (nislfv_rain_plm), rain (nislfv_rain_plm) and snow and graupel (nislfv_rain_plm6, the same code run once per species with a shared velocity). Same contract as sedimentation_upwind1_flux_D_lev, and MUTUALLY EXCLUSIVE with it (identical match, identical priority): a consumer imports one. q is the layer-centred falling density (rho times a mixing ratio, shape [lev]), W = -v_t <= 0 the layer-centred signed velocity along +lev (shape [lev]), and the consumer writes dq/dt = -D(W*q, wrt: lev)/rho + ... . The divergence is (F[k+1] - F[k])/dz[k] with F[j] = W_face[j] q_face[j] the flux through the lower face of layer j and F[NLEV+1] = 0 (no inflow through the model top). WHY NOT DONOR-CELL. WSM6's fallout increment converges as dt -> 0, but not to the donor-cell divergence of sedimentation_upwind1_flux_D_lev: on a WRF column (em_scm_xy cirrus, WRF step 60, real64 kernel replay at dtcld = 0.01/0.02/0.04 s) the two differ by 58 % Linf layer by layer, with ratios from 0.26 to 2.07. The converged operator has three ingredients donor-cell lacks, all transcribed from nislfv_rain_plm into one stencil template each: (1) a THIRD-ORDER interface velocity, W_face[j] = 9/16 (W[j] + W[j-1]) - 1/16 (W[j+1] + W[j-2]) on interior faces 3 <= j <= NLEV-1 (sedimentation_jh2010_face_velocity), two-point averages on faces 2 and NLEV (sedimentation_jh2010_end_face_velocity), W[1] on face 1, and WRF's override W_face[j] = W[j-1] wherever W[j] == 0 ('terminate of top of raingroup'); (2) a monotone piecewise-linear FACE VALUE whose slope is the AVERAGE of the two one-sided differences over the non-uniform thicknesses (sedimentation_jh2010_one_sided_slope, sedimentation_jh2010_unlimited_half_jump), zeroed at extrema and dropped if either face value would go negative (sedimentation_jh2010_half_jump), read at the lower face as q0 - h (sedimentation_jh2010_face_value) and at the upper face as q0 + h (sedimentation_jh2010_top_face_value), with no reconstruction in the bottom layer (WRF's qmi(1) = qa(1)) and, above the top layer, q = 0 and zero thickness beyond the lid (WRF's qa(km+1) = 0, dza(km+1) = wi(km+1) dt -> 0, bound once in sedimentation_jh2010_lid_face_flux); (3) an UPWIND choice of donor side at every interior face (sedimentation_jh2010_face_flux): the third-order velocity of (1) overshoots to an UPWARD face speed, W_face > 0, wherever a slow layer sits under a fast one (the edge of a falling layer), and there WRF's departure point lies above the face and the mass crossing it comes from the TOP face value of the layer BELOW (q[j-1] + h[j-1], WRF's qpi(j-1)), so F[j] = W_face[j] qpi[j-1] for W_face[j] > 0 and W_face[j] qmi[j] otherwise. Omitting (3) leaves 5.2e-4 of the peak tendency on a column with such a face, and 6e-7 to 6e-6 of the peak ice tendency on EqWeFiC's four WRF WSM6 ice columns (3-12 upward faces each); with it, 4.4e-10 and 1e-8 to 4e-8. Faces 1, 2 and NLEV cannot turn upward (a layer speed and two-point averages of non-positive speeds), so (3) lives only in the interior-face stencil and the end faces use sedimentation_jh2010_end_face_flux. THIS IS NOT A GENERIC PLM OR PPM LIMITER. On the WRF column a textbook minmod reconstruction in place of (2) leaves 31.6 % Linf of the difference to WRF, and the interface velocity alone leaves 37.2 %; only the scheme's own slope rule closes it. Do not 'simplify' the face value to a library limiter or the velocity to a donor value. EXCLUDED, deliberately: WRF's finite-dt velocity clamp ('diffusivity of wi', decfl = (wi[k+1] - wi[k]) dt/dz[k] > 0.05) is a property of the time step, vanishes as dt -> 0, and is not part of this operator; likewise the mean departure/arrival velocity iteration (iter >= 1), whose limit is the departure velocity used here. MEASURED. Transcribed in float64 and evaluated against WRF's own nislfv_rain_plm, extracted verbatim from mp_wsm6.F90 and Richardson-extrapolated to dt -> 0 on six stretched 12-layer columns (an ice-like bump aloft, a rain-like profile reaching the surface, a top-loaded profile, a random profile with many extrema, a layer with W == 0 but q > 0 that exercises the override, and a slow layer under a fast one that turns two interior face speeds upward and exercises (3)): relative Linf 2.4e-9, 3.8e-10, 9.4e-9, 7.8e-10, 1.1e-7 and 4.4e-10 against 42 %, 37 %, 19 %, 49 %, 278 % and 96 % for donor-cell, the fifth case limited by the extrapolation itself (7e-7), not the formula. On the WRF cirrus column the transcription matches the kernel at every layer below the top two to five significant figures. In this repository the rule is exercised by problems/sedimentation_column_jh2010_forced (the manufactured solution exp(-t) sin(pi z) with W = -(1 + z) on the shared stretched test column), whose recorded L2 errors fall with observed order 1.5: second order on the smooth monotone flanks and first order across the extremum, where switch (1) of the half-jump zeroes the slope. WHERE WRF IS NOT A REFERENCE: the kernel empties the TOP TWO LAYERS of the column on every call, by an amount independent of dt, and the mass does not reach the surface (EqWeFiC FORTRAN_BUGS B16); this rule keeps those layers conservative, so do not assert against a WSM6 fallout increment in the top two layers. SCOPE. Valid only for W <= 0 everywhere (a sign-changing velocity needs latlon3d's face-velocity upwind rule), and only for NLEV >= 5, so that the five layer regions below are distinct (a smaller NLEV is rejected at load, makearray_region_inverted). Exactly telescoping: sum_k dz[k] D[k] = -F[1] = -W[1] q[1], the surface precipitation flux, which a consumer reads back as for the donor-cell rule. Guaranteed order 1 (the extremum and positivity switches and the two first-order end faces), second order at smooth monotone interior points. THE REGIONS, each face spelled by exactly one template so that shared faces cancel exactly: [1,1] the bottom layer, F[2] (sedimentation_jh2010_end_face_flux over q[1..3]) minus F[1] = W[1] q[1] (the donor-cell rule's sedimentation_donor_flux_lev, WRF's qmi(1) = qa(1)); [2,2] F[3] (sedimentation_jh2010_face_flux, the first interior face) minus F[2]; [3, NLEV-2] the interior stencil sedimentation_jh2010_flux_D_lev_interior; [NLEV-1, NLEV-1] F[NLEV] (sedimentation_jh2010_lid_face_flux) minus F[NLEV-1] (the last interior face); [NLEV, NLEV] the top layer, 0 minus F[NLEV].

**References**

- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)
- Hong, S.-Y., & Lim, J.-O. J. (2006). The WRF single-moment 6-class microphysics scheme (WSM6). Journal of the Korean Meteorological Society, 42(2), 129-151.
- WRF v4.8.0 / NCAR MMM-physics 550b5b4, phys/physics_mmm/mp_wsm6.F90, subroutines nislfv_rain_plm (lines 1879-2122) and nislfv_rain_plm6 (2125-2400).

**Conformance fixtures**

- ast: [`sedimentation_jh2010_flux_D_lev`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/sedimentation_jh2010_flux_D_lev) — status `active`
- convergence: [`sedimentation_column_jh2010_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/sedimentation_column_jh2010_forced) — status `active`
- simulation: [`sedimentation_column_jh2010_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/sedimentation_column_jh2010_forced) — status `active`

#### Convergence — case [`sedimentation_column_jh2010_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/sedimentation_column_jh2010_forced)

Error norms read from the committed golden [`tests/conformance/convergence/sedimentation_column_jh2010_forced/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/sedimentation_column_jh2010_forced/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 2.291e-02 | — |
| 32 | 8.068e-03 | 1.51 |
| 64 | 2.794e-03 | 1.53 |
| 128 | 9.789e-04 | 1.51 |

Expected order: **1.5** (± 0.25).

![Convergence of sedimentation_column_jh2010_forced (log-log)](/plots/convergence/sedimentation_column_jh2010_forced.png)

### `sedimentation_upwind1_flux_D_lev`

<p class="rule-meta"><code>op:advection</code> <code>order:1</code> <code>bc:no-inflow</code> <code>axes:lev</code> <code>stencil width 2</code> <code>priority 10</code></p>

Source: [`grids/column_nonuniform_1d/rules/sedimentation_upwind1_flux_D_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/sedimentation_upwind1_flux_D_lev.esm)

**Rewrites** `D(W * q)/Dlev` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["W", "q"]}], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>q</code> is a bare field shaped <code>[lev]</code>; <code>W</code> is a bare field shaped <code>[lev]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NLEV - 1] = sedimentation_upwind1_flux_D_lev_interior&lt;W=W, q=q&gt;, [NLEV:NLEV] = (-sedimentation_donor_flux_lev&lt;W=W, n=NLEV, q=q&gt;) / dz[NLEV])</div>

Conservative first-order DONOR-CELL lowering of the flux divergence D(W*q, wrt: lev) of a DOWNWARD-falling quantity on column_nonuniform_1d -- the sedimentation (fallout) operator of a microphysics scheme. The consuming model writes its continuity law as dq/dt = -D(W*q, wrt: lev)/rho + ... with q the layer-centred falling density (rho times the hydrometeor mixing ratio, shape [lev]) and W the layer-centred signed velocity along +lev, W = -v_t <= 0 for a terminal fall speed v_t (shape [lev]); the rule lowers the whole compound to the two-region makearray [1, NLEV-1] (the interior donor-cell difference of stencils/sedimentation_upwind1_flux_D_lev_interior) and [NLEV, NLEV] (the top layer, closed here with NO INFLOW through the model top: (0 - F[NLEV]) / dz[NLEV], the flux through face NLEV+1 being literally absent rather than formed from a value beyond the column). Every face flux, in both regions, is one reference to stencils/sedimentation_donor_flux_lev, F[n] = W[n] q[n], bound at n = k+1, k and NLEV. THE VELOCITY IS INSIDE THE OPERATOR AND LAYER-CENTRED: unlike the face-staggered wind of an advection scheme, a terminal velocity is a property of the falling particles in a layer, so the mass leaving layer k through its lower face falls at layer k's own speed -- which is also the upwind (donor-cell) choice for a downward flow, since the upstream cell of the lower face of layer k is layer k. The two coincide only because the flow direction is fixed; the rule is therefore valid ONLY for W <= 0 everywhere and a consumer must not feed it a sign-changing velocity (use latlon3d's face-velocity upwind rule for that). The surface (face 1) needs no boundary treatment: the layer-1 outflow F[1] = W[1] q[1] is an ordinary donor flux that the consumer exposes as its precipitation rate -- no free names, no boundary values to supply. Exactly telescoping: sum_k dz[k] D[k] = -F[1], so the column loses exactly what precipitates. COMPOUND MATCH: pattern D(W*q, wrt: lev) at priority 10 with where {q: [lev], W: [lev]}; the velocity is written on the LEFT of the product (structural, non-commutative matching). It cannot fire on this grid's diffusion compounds D(K*D(u)) (K is [lev_nodes] and D(u) is not a bare variable, so both shape constraints fail, esm-spec 9.6.1) nor on a bare interface flux D(F) (F is [lev_nodes], plain-D rule at priority 0). ACCURACY AND SCOPE: O(dz), monotone, positivity-preserving; it is the METHOD-OF-LINES (instantaneous derivative) form of fallout. Operational schemes such as WSM6's forward semi-Lagrangian piecewise-linear fallout (nislfv_rain_plm, Juang & Hong 2010) are FULLY DISCRETE in time -- the Courant number v_t dt/dz lives inside the operator -- and their per-step increment does not equal dt times this derivative at a finite step; at t = 0 the exact divergence d(rho v_t q)/dz is scheme-independent, which is what derivative tests need. A higher-order (PLM/PPM) reconstruction of the donor value is a possible later rule on the same face-flux structure. Requires NLEV >= 2 (the grid's minimum; at NLEV = 1 the interior region [1, 0] is rejected at load as makearray_region_inverted).

**References**

- LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge University Press. Section 4.1 (the upwind / donor-cell flux).
- Hong, S.-Y., & Lim, J.-O. J. (2006). The WRF single-moment 6-class microphysics scheme (WSM6). Journal of the Korean Meteorological Society, 42(2), 129-151.
- Juang, H.-M. H., & Hong, S.-Y. (2010). Forward semi-Lagrangian advection with mass conservation and positive definiteness for falling hydrometeors. Monthly Weather Review, 138(5), 1778-1791. [`doi:10.1175/2009MWR3109.1`](https://doi.org/10.1175/2009MWR3109.1)

**Conformance fixtures**

- ast: [`sedimentation_upwind1_flux_D_lev`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/sedimentation_upwind1_flux_D_lev) — status `active`
- convergence: [`sedimentation_column_upwind_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/sedimentation_column_upwind_forced) — status `active`
- simulation: [`sedimentation_column_upwind_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/sedimentation_column_upwind_forced) — status `active`

#### Convergence — case [`sedimentation_column_upwind_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/sedimentation_column_upwind_forced)

Error norms read from the committed golden [`tests/conformance/convergence/sedimentation_column_upwind_forced/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/sedimentation_column_upwind_forced/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 7.229e-02 | — |
| 32 | 3.958e-02 | 0.87 |
| 64 | 2.094e-02 | 0.92 |
| 128 | 1.085e-02 | 0.95 |

Expected order: **1** (± 0.25).

![Convergence of sedimentation_column_upwind_forced (log-log)](/plots/convergence/sedimentation_column_upwind_forced.png)

### `varcoeff_face_flux_laplacian_lev`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:supplied_flux</code> <code>axes:lev</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/column_nonuniform_1d/rules/varcoeff_face_flux_laplacian_lev.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/varcoeff_face_flux_laplacian_lev.esm)

**Rewrites** `D(K * D(u)/Dlev - F)/Dlev` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "-", "args": [{"op": "*", "args": ["K", {"op": "D", "args": ["u"], "wrt": "lev"}]}, "F"]}], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[lev]</code>; <code>K</code> is a bare field shaped <code>[lev_nodes]</code>; <code>F</code> is a bare field shaped <code>[lev_nodes]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLEV - 1] = varcoeff_face_flux_laplacian_lev_interior&lt;F=F, K=K, u=u&gt;, [1:1] = (interface_total_flux_lev&lt;F=F, K=K, n=2, u=u&gt; - (-F[1])) / dz[1], [NLEV:NLEV] = (-F[NLEV + 1] - interface_total_flux_lev&lt;F=F, K=K, n=NLEV, u=u&gt;) / dz[NLEV])</div>

Conservative finite-volume lowering of the flux-form column diffusion operator d/dz( K du/dz - F ) on column_nonuniform_1d, with the diffusivity K and a PRESCRIBED interface flux F both supplied on the NLEV+1 interfaces (shape [lev_nodes]) and NO free names. The consumer writes the compound D( {op:-, args:[ {op:*, args:[K, D(u, wrt: lev)]}, F ]}, wrt: lev ) -- diffusive term on the LEFT, coefficient on the LEFT of its product (matching is structural and non-commutative) -- which is the natural spelling of a column conservation law du/dt = d/dz(K du/dz - F): one closed-form flux, diffusive minus prescribed, differentiated once. The rule fires on the whole compound at priority 10 (esm-spec 9.6.3 outermost-first compound precedence) and is match-scoped with where {u: {shape: [lev]}, K: {shape: [lev_nodes]}, F: {shape: [lev_nodes]}} (esm-spec 9.6.1), so it fires only on a bare layer field, a bare interface coefficient and a bare interface flux of THIS grid. Every flux in the operator is the ONE total interface flux H[n] = K[n] (u[n]-u[n-1])/(zc[n]-zc[n-1]) - F[n] of stencils/interface_total_flux_lev (built on stencils/interface_flux_K_lev, the flux varcoeff_face_laplacian_lev_flux_bc uses). Interior [2, NLEV-1]: the imported stencil, (H[k+1] - H[k]) / dz[k]. Boundary layers: the DIFFUSIVE part of the flux through a DOMAIN face is not computed from a ghost value; it is dropped, and the total flux through the surface interface (node 1) and the top interface (node NLEV+1) is the prescribed part alone, -F[1] and -F[NLEV+1] -- the places varcoeff_face_laplacian_lev_flux_bc reads its free names kdudz_bot / kdudz_top -- so L u[1] = ( H[2] + F[1] ) / dz[1] and L u[NLEV] = ( -F[NLEV+1] - H[NLEV] ) / dz[NLEV]. F is the UPWARD (+lev) flux, the same sign convention as the upward turbulent flux -K du/dz. So a WRF surface sensible-heat flux hfx (W m^-2, positive upward) into a potential-temperature equation is F[1] = hfx / (rho cp), a moisture flux qfx is F[1] = qfx / rho, and a surface stress on a wind component -- a DOWNWARD momentum flux, drag removing momentum from the column -- is F[1] = -(u_*^2) u/|U| (all in the layer-1 air density); a rigid model top is F[NLEV+1] = 0, and a PBL scheme's nonlocal (counter-gradient) and entrainment fluxes are the interior entries of F. Each diffused field (theta, u, v, qv, ...) therefore carries its own F and the one rule serves them all with no per-field free names. The match is structurally disjoint from those of varcoeff_face_laplacian_lev_flux_bc (D(K*D(u, lev), lev), no '- F') and face_flux_D_lev_supplied_faces (D(F, lev) on a bare interface field), so a model may import this rule beside either or both of them: each compound is lowered by the one rule whose spelling it uses. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). The operator telescopes exactly: sum_k dz[k] L u[k] = F[1] - F[NLEV+1], so the column integral of u changes only through the two prescribed domain-face fluxes. Minimum extent NLEV = 2: the interior region folds to the empty [2, 1] and the two face closures, H[2] = H[NLEV], cover both layers; at NLEV = 1 the interior region folds to the inverted [2, 0] and the load is rejected with makearray_region_inverted (esm-spec 4.3.2) rather than reading u[0] or writing layer 1 twice. Second-order in L2 on a smoothly stretched column (supraconvergent); MMS convergence is pinned by problems/heat_column_varcoeff_faceflux_forced, whose prescribed flux is nonzero at the surface and curved, so both face closures and the interior F difference are exercised.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Hong, S.-Y., Noh, Y., and Dudhia, J. (2006). A new vertical diffusion package with an explicit treatment of entrainment processes. Mon. Wea. Rev. 134, 2318-2341. (Surface, counter-gradient and entrainment fluxes of the YSU scheme as interface fluxes.) [`doi:10.1175/MWR3199.1`](https://doi.org/10.1175/MWR3199.1)

**Conformance fixtures**

- ast: [`varcoeff_face_flux_laplacian_lev`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_face_flux_laplacian_lev) — status `active`
- convergence: [`heat_column_varcoeff_faceflux_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_faceflux_forced) — status `active`
- simulation: [`heat_column_varcoeff_faceflux_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_column_varcoeff_faceflux_forced) — status `active`

#### Convergence — case [`heat_column_varcoeff_faceflux_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_faceflux_forced)

Error norms read from the committed golden [`tests/conformance/convergence/heat_column_varcoeff_faceflux_forced/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_column_varcoeff_faceflux_forced/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.052e-03 | — |
| 32 | 7.636e-04 | 2.00 |
| 64 | 1.909e-04 | 2.00 |
| 128 | 4.773e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of heat_column_varcoeff_faceflux_forced (log-log)](/plots/convergence/heat_column_varcoeff_faceflux_forced.png)

### `varcoeff_face_laplacian_lev_flux_bc`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:prescribed_flux</code> <code>axes:lev</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/column_nonuniform_1d/rules/varcoeff_face_laplacian_lev_flux_bc.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d/rules/varcoeff_face_laplacian_lev_flux_bc.esm)

**Rewrites** `D(K * D(u)/Dlev)/Dlev` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["K", {"op": "D", "args": ["u"], "wrt": "lev"}]}], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>u</code> is a bare field shaped <code>[lev]</code>; <code>K</code> is a bare field shaped <code>[lev_nodes]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:NLEV - 1] = varcoeff_face_laplacian_lev_interior&lt;K=K, u=u&gt;, [1:1] = (interface_flux_K_lev&lt;K=K, n=2, u=u&gt; - kdudz_bot) / dz[1], [NLEV:NLEV] = (kdudz_top - interface_flux_K_lev&lt;K=K, n=NLEV, u=u&gt;) / dz[NLEV])</div>

Conservative finite-volume variable-coefficient vertical diffusion d/dz(K du/dz) on column_nonuniform_1d with the coefficient K on the INTERFACES and PRESCRIBED diffusive fluxes at the surface and the model top. The consumer writes the compound D( {op:*, args:[K, D(u, wrt: lev)]}, wrt: lev ) -- coefficient on the LEFT of the product (matching is structural and non-commutative) -- and the rule fires on the whole compound at priority 10 (esm-spec 9.6.3 outermost-first compound precedence) before any plain-D rule can lower the inner derivative; it is match-scoped with where {u: {shape: [lev]}, K: {shape: [lev_nodes]}} (esm-spec 9.6.1), so it fires only when the state is a bare layer field and the coefficient a bare interface field of THIS grid (a cell-centred coefficient does not satisfy the K constraint and falls through to whichever rule the consumer imports for that case). Every flux in the operator is the ONE interface flux G[n] = K[n] (u[n]-u[n-1])/(zc[n]-zc[n-1]) of stencils/interface_flux_K_lev (no face averaging: K[k] is the value on the interface below layer k, node 1 = surface, so a PBL scheme's exchange coefficient 'between k and k+1' is K[k+1]). Interior [2, NLEV-1]: the imported stencil, (G[k+1] - G[k]) / dz[k]. Boundary layers: the flux through a DOMAIN face is not computed from a ghost value but taken from the free names `kdudz_bot` and `kdudz_top`, the consumer-supplied values of K du/dz ON the surface interface (z = ze[1]) and the top interface (z = ze[NLEV+1]), positive in the +lev (upward) direction, so L u[1] = ( G[2] - kdudz_bot ) / dz[1] and L u[NLEV] = ( kdudz_top - G[NLEV] ) / dz[NLEV], the same interface-flux template applied at n = 2 and n = NLEV. Because every interior flux enters the column sum once with each sign, the operator telescopes exactly: sum_k dz[k] L u[k] = kdudz_top - kdudz_bot, so the column integral of u changes only through the two prescribed fluxes. The consumer declares the two free names as ordinary parameters or observeds: zero for an insulated face (the zero-flux/rigid-lid case, kdudz_top = 0 in every atmospheric column), and for a surface exchange the negative of the kinematic surface flux, since the upward turbulent flux is -K du/dz -- a WRF surface sensible-heat flux hfx (W m^-2) into the potential-temperature equation is kdudz_bot = -hfx / (rho cp), a moisture flux qfx (kg m^-2 s^-1) is kdudz_bot = -qfx / rho, and a surface stress on a wind component is kdudz_bot = -(u_* ^2) u/|U| (all in the layer-1 air density) -- which is precisely how the WRF YSU and MYNN schemes inject their surface fluxes into layer 1. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). Minimum extent NLEV = 2: the interior region folds to the empty [2, 1] and the two face closures, G[2] = G[NLEV], cover both layers; at NLEV = 1 the interior region folds to the inverted [2, 0] and the load is rejected with makearray_region_inverted (esm-spec 4.3.2) rather than reading u[0]. Second-order in L2 on a smoothly stretched column (supraconvergent); MMS convergence is pinned by problems/heat_column_varcoeff_flux_bc and problems/heat_column_varcoeff_forced. Its match is structurally disjoint from that of rules/varcoeff_face_flux_laplacian_lev (D(K*D(u, lev) - F, lev), which carries every prescribed flux in an interface field F and has no free names), so a model may import both and spell each field's diffusion in whichever form suits it.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Hong, S.-Y., Noh, Y., and Dudhia, J. (2006). A new vertical diffusion package with an explicit treatment of entrainment processes. Mon. Wea. Rev. 134, 2318-2341. (Surface-flux boundary condition of the YSU scheme.) [`doi:10.1175/MWR3199.1`](https://doi.org/10.1175/MWR3199.1)

**Conformance fixtures**

- ast: [`varcoeff_face_laplacian_lev_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_face_laplacian_lev_flux_bc) — status `active`
- convergence: [`heat_column_varcoeff_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc) — status `active`
- convergence: [`heat_column_varcoeff_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_forced) — status `active`
- simulation: [`heat_column_varcoeff_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_column_varcoeff_flux_bc) — status `active`
- simulation: [`heat_column_varcoeff_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_column_varcoeff_forced) — status `active`

#### Convergence — case [`heat_column_varcoeff_flux_bc`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc)

Error norms read from the committed golden [`tests/conformance/convergence/heat_column_varcoeff_flux_bc/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 7.811e-05 | — |
| 32 | 2.104e-05 | 1.89 |
| 64 | 5.403e-06 | 1.96 |
| 128 | 1.361e-06 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of heat_column_varcoeff_flux_bc (log-log)](/plots/convergence/heat_column_varcoeff_flux_bc.png)

#### Convergence — case [`heat_column_varcoeff_forced`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_forced)

Error norms read from the committed golden [`tests/conformance/convergence/heat_column_varcoeff_forced/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_column_varcoeff_forced/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 3.050e-03 | — |
| 32 | 7.631e-04 | 2.00 |
| 64 | 1.908e-04 | 2.00 |
| 128 | 4.770e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of heat_column_varcoeff_forced (log-log)](/plots/convergence/heat_column_varcoeff_forced.png)

