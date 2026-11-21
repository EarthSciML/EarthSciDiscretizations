---
title: "mom6_cgrid"
description: "The Arakawa C-grid of the MOM6 ocean model: the cell axes `i` (x) and `j` (y) of the N_I x N_J thickness (T) cells, plus the layer axis `k` of the N_K layers that carry the per-layer state (the layer thickness h and the layer face velocities)."
esd_kind: "grid"
families: ["mom6_cgrid"]
ops: ["advection", "avg", "clamp", "diff", "div", "flux", "grad", "limit", "mask", "pvf", "reconstruct", "sum", "vort"]
bcs: ["wall"]
tags: ["esd:grid", "family:mom6_cgrid", "ndim:2", "stagger:c", "axis:i", "axis:j", "axis:k", "esd:stencil", "grid:mom6_cgrid", "op:avg", "op:clamp", "op:limit", "op:reconstruct", "op:flux", "op:advection", "op:sum", "op:pvf", "esd:rule", "order:0", "bc:wall", "axes:i,j,k", "op:diff", "order:1", "op:div", "spec:9.6.8", "order:3", "op:grad", "op:mask", "axes:i,j", "op:vort"]
source: "grids/mom6_cgrid/grid.esm"
---

The Arakawa C-grid of the MOM6 ocean model: the cell axes `i` (x) and `j` (y) of the N_I x N_J thickness (T) cells, plus the layer axis `k` of the N_K layers that carry the per-layer state (the layer thickness h and the layer face velocities). The C-grid staggering: the zonal velocity u lives on the EAST faces of the T cells and the meridional velocity v on the NORTH faces; a face flux of a cell is the face velocity times the face length times the reconstructed upwind edge thickness. N_I, N_J, N_K are the esm-spec 9.7 metaparameters, bound at an import edge or the loader API; the consuming EqESM stage-2 model (stage2/mom6_dynamics.esm) binds N_I = N_J = 8, N_K = 2 for the traced 8x8 subdomain (global T cells i = -1..6, j = 74..81 of the stage-1 continuity dump) and the pinned dump layers k = 1 and 37. FREE-NAME CONTRACT: a consuming model supplies the real-valued geometry and state as its own variables with the names u_E/v_N (face velocities, shape [i,j,k]), dy_Cu_E/dx_Cv_N (east/north face lengths, [i,j]), mask (0/1 ocean-cell indicator, [i,j]), IdxT/IdxT_e and IdyT/IdyT_n (the 1/dxT of the cell and its east neighbor, the 1/dyT of the cell and its north neighbor, [i,j]), IareaT (1/cell area, [i,j]), h_min_recon (the PPM_limit_pos thickness floor) and dt (the tracer time step) as scalars — exactly the columns of the stage-1 continuity dump (mom6_continuity.csv) — together with the layer thickness field h [i,j,k] the rules advect; the rules below read them by bare name, so a consumer MAY either bind them from its own state (the stage-2 pin tests override them from the dump) or derive them through this file's templates. The single geometry template is the exact MOM6 reciprocal bridge (Adcroft_reciprocal, MOM_grid_initialize.F90:1173-1179, applied to the grid metrics at MOM_grid.F90:420-424): a consumer holding dxT/dyT/areaT (or dx_Cu/dy_Cv for the face metrics of a future pressure-gradient component) derives IdxT/IdyT/IareaT/IdxCu/IdyCv through it. The stage-2 model binds the dumped reciprocals directly so its pins compare the rules against the Fortran's own values bit-for-bit. The rules are the PPM continuity pair (stage-1 dump description; MOM_continuity_PPM.F90): mom6_cgrid_flux_x/y lower the per-face PPM layer volume flux uh = dy*u*int(parabola) (zonal_flux_layer F90:903-979, merid_flux_layer F90:1801-1885) and mom6_cgrid_div_x/y lower the C-grid flux-form divergence D(uh,x) = (uh(i,j) - uh(i-1,j))*IareaT, together the thickness continuity dh/dt = -div(uh) - div(vh). The reconstruction stencils transcribe PPM_reconstruction_x/y (F90:2327-2562) and PPM_limit_pos (F90:2636-2674) of the reference run (FIRST_DIRECTION=0, MONOTONIC_CONTINUITY=False, SIMPLE_2ND_PPM=False, UPWIND_1ST=False, VOL_CFL=False), including the masked-neighbor edge replacement and the zero slope on masked cells (G%mask2dT = 0 over land), as a DAG of one-template laws: the center slope mom6_ppm_slope (over mom6_ppm_slope_raw/_mag), the raw edges mom6_ppm_edge_w/_e (over mom6_ppm_edge), the PPM_limit_pos scale mom6_ppm_cell_scale (over mom6_ppm_limit_scale and its mom6_ppm_curv/_dh/_limit_scale_num/_limit_scale_den sub-laws), the limited edge mom6_ppm_limited_edge, the parabola integral mom6_ppm_face_integral, and the upwind-cell face fluxes mom6_ppm_upwind_flux_pos/_neg that the flux rules apply, one per flow direction, to the upwind cell's 5-cell stencil.

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `N_I` | integer | `8` | Number of T cells along the C-grid x axis (cell-based axis), [-]. |
| `N_J` | integer | `8` | Number of T cells along the C-grid y axis (cell-based axis), [-]. |
| `N_K` | integer | `2` | Number of layers carried by the component (the traced continuity dump pins layers 1 and 37 of the reference run), [-]. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>IareaT</code> <code>KE</code> <code>dxCu</code> <code>dyCv</code></p>

These names are supplied by the consuming model's scope at evaluation (the keyed-factor / free-name geometry contract); they are neither metaparameters (load-time integers) nor template parameters (bound per apply site).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `i` | interval | \(N_{\mathrm{I}}\) |
| `j` | interval | \(N_{\mathrm{J}}\) |
| `k` | interval | \(N_{\mathrm{K}}\) |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `mom6_adcroft_reciprocal` | `x` | \(\begin{cases}0 & \text{if } x = 0 \\ \frac{1}{x} & \text{otherwise}\end{cases}\) | The guarded reciprocal of the MOM6 grid metrics (Adcroft_reciprocal, MOM_grid_initialize.F90:1173-1179): 1/x, or exactly 0.0 when x == 0.0. The grid's metrics are derived with it at MOM_grid.F90:420-424 (IdxT = 1/dxT, IdyT = 1/dyT, IareaT = 1/areaT); a consumer holding the positive metric lengths applies this template to each to obtain the free-name reciprocals (the dump supplies the reciprocals directly, so the traced stage-2 model binds them rather than deriving). |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `mom6_cell_face_mean`

<p class="rule-meta"><code>op:avg</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_cell_face_mean.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_cell_face_mean.esm)

Parameters: `x_e`, `x_w`, `x_n`, `x_s`

<div class="math-block">$$0.25 \left(\left(x_{\mathrm{e}} + x_{\mathrm{w}}\right) + \left(x_{\mathrm{n}} + x_{\mathrm{s}}\right)\right)$$</div>

The four-face mean of a T cell (MOM_CoriolisAdv.F90 gradKE :1290-1297, the KE_ARAKAWA sum without its IareaT factor): mean = 0.25*((x_e + x_w) + (x_n + x_s)), the cell's east and west u-face values plus its north and south v-face values, in the Fortran's grouping. mom6_cgrid_avg_cell binds it at every cell; a closed wall face contributes 0 (MOM6's zero wall velocity). Pure scalar algebra.

**References**

- Sadourny, R. (1975). The dynamics of finite-difference models of the shallow-water equations. Journal of the Atmospheric Sciences, 32(4), 680-689. [`doi:10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2`](https://doi.org/10.1175/1520-0469%281975%29032%3C0680%3ATDOFDM%3E2.0.CO%3B2)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

### `mom6_coriolis_clamp4`

<p class="rule-meta"><code>op:clamp</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_coriolis_clamp4.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_coriolis_clamp4.esm)

Parameters: `x`, `f1`, `f2`, `f3`, `f4`

<div class="math-block">$$\min(\max(x, \min(\mathit{f1}, \min(\mathit{f2}, \min(\mathit{f3}, \mathit{f4})))), \max(\mathit{f1}, \max(\mathit{f2}, \max(\mathit{f3}, \mathit{f4}))))$$</div>

The bound_Coriolis clamp of MOM6 (MOM_CoriolisAdv.F90 CorAdCalc :897-908 for CAu, :1139-1150 for CAv): clamp4 = min(max(x, min(f1, f2, f3, f4)), max(f1, f2, f3, f4)), the raw Coriolis term x limited to the range of the four products f1..f4 of the absolute vorticity and the transverse velocity at the q-points around the face (fv1..fv4 = abs_vort*v for CAu, fu1..fu4 = -abs_vort*u for CAv). The Fortran applies CAu = min(CAu, max_fv) then CAu = max(CAu, min_fv); min(max(x, lo), hi) gives the same bits whenever lo <= hi, which always holds here. Each product is bound once, so the min and the max read the same four values. At a solid wall the rules bind a product whose transverse velocity sits on a closed wall face (or on a land face beyond the domain) to 0, so 0 enters the min and the max exactly as MOM6's zero wall velocity makes it. Pure scalar algebra.

**References**

- Sadourny, R. (1975). The dynamics of finite-difference models of the shallow-water equations. Journal of the Atmospheric Sciences, 32(4), 680-689. [`doi:10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2`](https://doi.org/10.1175/1520-0469%281975%29032%3C0680%3ATDOFDM%3E2.0.CO%3B2)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

### `mom6_ppm_cell_scale`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_cell_scale.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_cell_scale.esm)

Parameters: `h_m2`, `h_m1`, `h_0`, `h_p1`, `h_p2`, `m_m2`, `m_m1`, `m_0`, `m_p1`, `m_p2`, `h_min_recon`

<div class="math-block">$$\operatorname{mom6\_ppm\_limit\_scale}(\mathit{hL\_raw} = \operatorname{mom6\_ppm\_edge\_w}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}), \mathit{hR\_raw} = \operatorname{mom6\_ppm\_edge\_e}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}}), h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}})$$</div>

The PPM_limit_pos scale of one cell from its 5-cell stencil (PPM_reconstruction_x then PPM_limit_pos, MOM_continuity_PPM.F90:2430-2433, 2665-2677): mom6_ppm_limit_scale with hL_raw = mom6_ppm_edge_w and hR_raw = mom6_ppm_edge_e of the cell, h_0 the cell's thickness and h_min_recon the floor. Both limited edges of the cell (mom6_ppm_upwind_flux_pos/neg) scale toward h_0 by this one number. Positional stencil params h_m2..h_p2 / m_m2..m_p2 (two cells left of the cell to two cells right of it, along the reconstruction axis). Pure scalar algebra, dimensionless.

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_curv`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_curv.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_curv.esm)

Parameters: `hL_raw`, `hR_raw`, `h_0`

<div class="math-block">$$3 \left(\mathit{hL\_raw} + \mathit{hR\_raw} - 2 h_{\mathrm{0}}\right)$$</div>

The grid-normalized curvature of one cell's PPM parabola (PPM_limit_pos, MOM_continuity_PPM.F90:2665): curv = 3*((hL_raw + hR_raw) - 2*h_0), positive exactly when the parabola through the raw edges and the cell mean is convex (a minimum, the only case PPM_limit_pos limits). Read by mom6_ppm_limit_scale and its numerator/denominator sub-laws mom6_ppm_limit_scale_num/den. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The positive-definiteness argument behind PPM_limit_pos.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_dh`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_dh.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_dh.esm)

Parameters: `hL_raw`, `hR_raw`

<div class="math-block">$$\mathit{hR\_raw} - \mathit{hL\_raw}$$</div>

The edge difference of one cell's PPM parabola (PPM_limit_pos, MOM_continuity_PPM.F90:2667): dh = hR_raw - hL_raw. With mom6_ppm_curv it locates the parabola's extremum (|dh| < curv puts the minimum inside the cell) and sizes its depth (curv^2 + 3*dh^2, mom6_ppm_limit_scale_den). Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The positive-definiteness argument behind PPM_limit_pos.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_edge`

<p class="rule-meta"><code>op:reconstruct</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_edge.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_edge.esm)

Parameters: `h_side`, `h_0`, `slp_side`, `slp_0`, `m_side`

<div class="math-block">$$0.5 \left(\left(m_{\mathrm{side}} h_{\mathrm{side}} + \left(1 - m_{\mathrm{side}}\right) h_{\mathrm{0}}\right) + h_{\mathrm{0}}\right) + \frac{1}{6} \left(\mathit{slp\_side} - \mathit{slp\_0}\right)$$</div>

One raw edge of the MOM6 PPM parabola (PPM_reconstruction_x, MOM_continuity_PPM.F90:2428-2449), shared by the left and right edges of a cell: edge = 0.5*(h_side_masked + h_0) + (slp_side - slp_0)/6 with h_side_masked = m_side*h_side + (1 - m_side)*h_0 (a masked neighbor is replaced by the cell's own thickness, the Fortran's h_im1/h_ip1 masked replacement at F90:2430-2431 and 2441-2442). The caller binds the two edges with swapped slope roles: the LEFT edge binds slp_side to the west neighbor's slope and slp_0 to the cell's own (hL = 0.5*(h_im1 + h_0) + oneSixth*(slp(i-1) - slp(i))), the RIGHT edge binds slp_side to the cell's own slope and slp_0 to the east neighbor's (hR = 0.5*(h_ip1 + h_0) + oneSixth*(slp(i) - slp(i+1))). One law, both edges; the un-limited (raw) edge feeds mom6_ppm_limit_scale. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (Eq. 1.6, the edge interpolation.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_edge_e`

<p class="rule-meta"><code>op:reconstruct</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_edge_e.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_edge_e.esm)

Parameters: `h_m1`, `h_0`, `h_p1`, `h_p2`, `m_m1`, `m_0`, `m_p1`, `m_p2`

<div class="math-block">$$\operatorname{mom6\_ppm\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{side}} = h_{\mathrm{p1}}, m_{\mathrm{side}} = m_{\mathrm{p1}}, \mathit{slp\_0} = \operatorname{mom6\_ppm\_slope}(h_{\mathrm{0}} = h_{\mathrm{p1}}, h_{\mathrm{m1}} = h_{\mathrm{0}}, h_{\mathrm{p1}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{p1}}, m_{\mathrm{m1}} = m_{\mathrm{0}}, m_{\mathrm{p1}} = m_{\mathrm{p2}}), \mathit{slp\_side} = \operatorname{mom6\_ppm\_slope}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}))$$</div>

The raw (un-limited) RIGHT (east/north) edge of one cell's MOM6 PPM parabola from the cell's 4-cell stencil (PPM_reconstruction_x, MOM_continuity_PPM.F90:2431, 2433): h_E = mom6_ppm_edge with h_side = h_p1, m_side = m_p1, slp_side = the cell's own slope mom6_ppm_slope(h_m1, h_0, h_p1) and slp_0 = the right neighbour's slope mom6_ppm_slope(h_0, h_p1, h_p2), i.e. h_E = 0.5*(h_ip1 + h_0) + (slp(i) - slp(i+1))/6 with the masked-neighbour replacement h_ip1. The mirror of mom6_ppm_edge_w (the slope roles swap, which is why the two edges are two laws over the one mom6_ppm_edge); positional stencil params as there. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_edge_w`

<p class="rule-meta"><code>op:reconstruct</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_edge_w.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_edge_w.esm)

Parameters: `h_m2`, `h_m1`, `h_0`, `h_p1`, `m_m2`, `m_m1`, `m_0`, `m_p1`

<div class="math-block">$$\operatorname{mom6\_ppm\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{side}} = h_{\mathrm{m1}}, m_{\mathrm{side}} = m_{\mathrm{m1}}, \mathit{slp\_0} = \operatorname{mom6\_ppm\_slope}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}), \mathit{slp\_side} = \operatorname{mom6\_ppm\_slope}(h_{\mathrm{0}} = h_{\mathrm{m1}}, h_{\mathrm{m1}} = h_{\mathrm{m2}}, h_{\mathrm{p1}} = h_{\mathrm{0}}, m_{\mathrm{0}} = m_{\mathrm{m1}}, m_{\mathrm{m1}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{0}}))$$</div>

The raw (un-limited) LEFT (west/south) edge of one cell's MOM6 PPM parabola from the cell's 4-cell stencil (PPM_reconstruction_x, MOM_continuity_PPM.F90:2430, 2432): h_W = mom6_ppm_edge with h_side = h_m1, m_side = m_m1, slp_side = the left neighbour's slope mom6_ppm_slope(h_m2, h_m1, h_0) and slp_0 = the cell's own slope mom6_ppm_slope(h_m1, h_0, h_p1), i.e. h_W = 0.5*(h_im1 + h_0) + (slp(i-1) - slp(i))/6 with the masked-neighbour replacement h_im1. The stencil params are positional relative to the cell (m2 = two cells to the left, ..., p1 = one cell to the right) so the same law serves both axes and every cell the rules bind it to. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_face_integral`

<p class="rule-meta"><code>op:flux</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_face_integral.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_face_integral.esm)

Parameters: `u`, `dt`, `dy`, `I_dx`, `h_face`, `h_back`, `h_0`

<div class="math-block">$$\begin{cases}0 &amp; \text{if } u = 0 \\ \mathit{dy} \left(u \left(h_{\mathrm{face}} + \left(\left|u\right| \left(\mathit{dt} I_{\mathrm{dx}}\right)\right) \left(0.5 \left(h_{\mathrm{back}} - h_{\mathrm{face}}\right) + \left(h_{\mathrm{face}} + h_{\mathrm{back}} - 2 h_{\mathrm{0}}\right) \left(\left|u\right| \left(\mathit{dt} I_{\mathrm{dx}}\right) - 1.5\right)\right)\right)\right) &amp; \text{otherwise}\end{cases}$$</div>

The MOM6 PPM layer-flux integral at one cell face (zonal_flux_layer, MOM_continuity_PPM.F90:933-956, VOL_CFL=False), ONE law for both axes and both flow directions: flux = dy*u*(h_face + CFL*(0.5*(h_back - h_face) + curv_3*(CFL - 1.5))) with CFL = |u|*dt*I_dx (the upwind-cell Courant number) and curv_3 = (h_face + h_back) - 2*h_0, and flux = 0 for u == 0 (the Fortran's else branch, F90:953-955). The caller selects the UPWIND cell: for u > 0 the upwind cell is the cell itself (h_face = its east edge, h_back = its west edge, h_0 = its thickness, I_dx = its 1/dxT, the F90:933-939 branch: uh = dy*u*(h_E(i) + CFL*(0.5*(h_W(i) - h_E(i)) + curv_3*(CFL - 1.5)))); for u < 0 the upwind cell is the east/north neighbor (h_face = its west/south edge, h_back = its east/north edge, h_0 = its thickness, I_dx = its 1/dxT, the F90:940-949 branch: uh = dy*u*(h_W(i+1) + CFL*(0.5*(h_E(i+1) - h_W(i+1)) + curv_3*(CFL - 1.5)))); the two Fortran branches are the same parabola integral with the two edge roles swapped, which is why one law serves both. The edges h_face/h_back are the PPM_limit_pos-limited edges of the upwind cell (mom6_ppm_limited_edge), the face length dy is dy_Cu (zonal) or dx_Cv (meridional), and the flux is the layer volume flux [m3 s-1]. Pure scalar algebra; no branch inside the law beyond the zero-velocity floor, so the sign selection lives at the rule call site.

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The parabolic face integral.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_limit_scale`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_limit_scale.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_limit_scale.esm)

Parameters: `hL_raw`, `hR_raw`, `h_0`, `h_min_recon`

<div class="math-block">$$\begin{cases}\begin{cases}\begin{cases}0 &amp; \text{if } h_{\mathrm{0}} \le h_{\mathrm{min\_recon}} \\ \frac{\operatorname{mom6\_ppm\_limit\_scale\_num}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}})}{\operatorname{mom6\_ppm\_limit\_scale\_den}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}})} &amp; \text{if } \operatorname{mom6\_ppm\_limit\_scale\_num}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}) &lt; \operatorname{mom6\_ppm\_limit\_scale\_den}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) \\ 1 &amp; \text{otherwise}\end{cases} &amp; \text{if } \left|\operatorname{mom6\_ppm\_dh}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw})\right| &lt; \operatorname{mom6\_ppm\_curv}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) \\ 1 &amp; \text{otherwise}\end{cases} &amp; \text{if } \operatorname{mom6\_ppm\_curv}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) &gt; 0 \\ 1 &amp; \text{otherwise}\end{cases}$$</div>

The positive-definiteness limiter scale of the MOM6 PPM reconstruction (PPM_limit_pos, MOM_continuity_PPM.F90:2665-2677): with curv = mom6_ppm_curv = 3*((hL_raw + hR_raw) - 2*h_0), dh = mom6_ppm_dh = hR_raw - hL_raw, num = mom6_ppm_limit_scale_num = 12*curv*(h_0 - h_min_recon) and den = mom6_ppm_limit_scale_den = curv^2 + 3*dh^2, scale = num/den if curv > 0 and |dh| < curv and h_0 > h_min_recon and num < den; scale = 0 (the constant h_0 reconstruction) if h_0 <= h_min_recon; scale = 1 otherwise (the raw edges are kept). The limited edges are h_L = h_0 + scale*(h_L_raw - h_0) and h_R = h_0 + scale*(h_R_raw - h_0), applied by the sibling mom6_ppm_limited_edge. The h_min_recon floor is 2*GV%Angstrom_H (the reference run: h_min_recon = 2e-10 m). Pure scalar algebra, dimensionless.

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The positive-definiteness argument behind PPM_limit_pos.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_limit_scale_den`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_limit_scale_den.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_limit_scale_den.esm)

Parameters: `hL_raw`, `hR_raw`, `h_0`

<div class="math-block">$$\operatorname{mom6\_ppm\_curv}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) \operatorname{mom6\_ppm\_curv}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) + 3 \left(\operatorname{mom6\_ppm\_dh}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}) \operatorname{mom6\_ppm\_dh}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw})\right)$$</div>

The denominator of the PPM_limit_pos scale (MOM_continuity_PPM.F90:2671, 2674): den = curv^2 + 3*dh^2 with curv = mom6_ppm_curv and dh = mom6_ppm_dh (the Fortran's curv**2 + 3.0*dh**2, spelled as products). den/(12*curv) is the depth of the parabola's minimum below the cell mean. Pure scalar algebra, [m^2].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The positive-definiteness argument behind PPM_limit_pos.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_limit_scale_num`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_limit_scale_num.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_limit_scale_num.esm)

Parameters: `hL_raw`, `hR_raw`, `h_0`, `h_min_recon`

<div class="math-block">$$12 \left(\operatorname{mom6\_ppm\_curv}(\mathit{hL\_raw} = \mathit{hL\_raw}, \mathit{hR\_raw} = \mathit{hR\_raw}, h_{\mathrm{0}} = h_{\mathrm{0}}) \left(h_{\mathrm{0}} - h_{\mathrm{min\_recon}}\right)\right)$$</div>

The numerator of the PPM_limit_pos scale (MOM_continuity_PPM.F90:2671, 2674): num = 12*curv*(h_0 - h_min_recon) with curv = mom6_ppm_curv. The parabola's minimum is h_0 - den/(12*curv) (den = mom6_ppm_limit_scale_den), so num < den is exactly the test that the minimum undershoots the floor h_min_recon, and num/den is the scale that lifts it back to the floor. Pure scalar algebra, [m^2].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. (The positive-definiteness argument behind PPM_limit_pos.) [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_limited_edge`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_limited_edge.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_limited_edge.esm)

Parameters: `h_0`, `scale`, `raw`

<div class="math-block">$$h_{\mathrm{0}} + \mathit{scale} \left(\mathit{raw} - h_{\mathrm{0}}\right)$$</div>

One PPM_limit_pos-limited edge of a cell (MOM_continuity_PPM.F90:2663-2667): the raw edge scaled toward the cell mean, edge = h_0 + scale*(raw - h_0), where scale is the mom6_ppm_limit_scale of the cell (1 = the raw edge kept, 0 = the flat h_0 reconstruction, (0, 1) = partially limited). The single law replaces the four inline copies of the Fortran's h_L = h_in + scale*(h_L - h_in) / h_R = h_in + scale*(h_R - h_in) applications (F90:2665-2666); the caller binds raw to either the left or the right raw edge of the cell. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_slope`

<p class="rule-meta"><code>op:reconstruct</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_slope.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_slope.esm)

Parameters: `h_m1`, `h_0`, `h_p1`, `m_m1`, `m_0`, `m_p1`

<div class="math-block">$$\begin{cases}\begin{cases}\operatorname{mom6\_ppm\_slope\_mag}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}) &amp; \text{if } \operatorname{mom6\_ppm\_slope\_raw}(h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}) &gt; 0 \\ -\operatorname{mom6\_ppm\_slope\_mag}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}) &amp; \text{if } \operatorname{mom6\_ppm\_slope\_raw}(h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}) &lt; 0 \\ 0 &amp; \text{otherwise}\end{cases} &amp; \text{if } m_{\mathrm{m1}} \left(m_{\mathrm{0}} m_{\mathrm{p1}}\right) &gt; 0.5 \\ 0 &amp; \text{otherwise}\end{cases}$$</div>

The monotonicity-limited slope of the MOM6 PPM reconstruction at a cell center (PPM_reconstruction_x, MOM_continuity_PPM.F90:2397-2406, the monotonic=False branch of the reference run): the sign of the raw slope s = mom6_ppm_slope_raw = 0.5*(h_p1 - h_m1) times the Lin (1994) Eq. B2 limited magnitude mom6_ppm_slope_mag = min(|s|, 2*min(dMx, dMn)) (the Fortran's sign(1., s)*min(abs(s), 2*min(dMx, dMn))), and s = 0 if the cell is masked (the 0/1 ocean-cell indicator product m_m1*m_0*m_p1 == 0, G%mask2dT = 0 over land) or if s == 0 exactly. Pure scalar algebra (one cell's five stencil values), the center-slope piece of the C-grid PPM reconstruction; grid-agnostic. The layer thicknesses h_m1/h_0/h_p1 are of one layer, [m]; the masks are dimensionless.

**References**

- Lin, S.-J. (1994). A nearly conservative flux-reconstruction scheme. Monthly Weather Review, 122(5), 994-1011. (Eq. B2, the monotonic slope constraint.)
- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_slope_mag`

<p class="rule-meta"><code>op:limit</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_slope_mag.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_slope_mag.esm)

Parameters: `h_m1`, `h_0`, `h_p1`

<div class="math-block">$$\min(\left|\operatorname{mom6\_ppm\_slope\_raw}(h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}})\right|, 2 \min(\max(h_{\mathrm{p1}}, \max(h_{\mathrm{m1}}, h_{\mathrm{0}})) - h_{\mathrm{0}}, h_{\mathrm{0}} - \min(h_{\mathrm{p1}}, \min(h_{\mathrm{m1}}, h_{\mathrm{0}}))))$$</div>

The Lin (1994) Eq. B2 limited magnitude of the MOM6 PPM center slope (PPM_reconstruction_x, MOM_continuity_PPM.F90:2403-2405): |s| = min(|mom6_ppm_slope_raw|, 2*min(dMx, dMn)) with dMx = max(h_p1, h_m1, h_0) - h_0 and dMn = h_0 - min(h_p1, h_m1, h_0), so the slope never pushes an edge outside the range of the cell and its two neighbours. mom6_ppm_slope restores the sign. Pure scalar algebra, [m].

**References**

- Lin, S.-J. (1994). A nearly conservative flux-reconstruction scheme. Monthly Weather Review, 122(5), 994-1011. (Eq. B2, the monotonic slope constraint.)

### `mom6_ppm_slope_raw`

<p class="rule-meta"><code>op:reconstruct</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_slope_raw.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_slope_raw.esm)

Parameters: `h_m1`, `h_p1`

<div class="math-block">$$0.5 \left(h_{\mathrm{p1}} - h_{\mathrm{m1}}\right)$$</div>

The unlimited centered slope of one cell of the MOM6 PPM reconstruction (PPM_reconstruction_x, MOM_continuity_PPM.F90:2401, slp = 0.5*(h_in(i+1) - h_in(i-1))): s = 0.5*(h_p1 - h_m1). The sign test and the magnitude of mom6_ppm_slope both read it; mom6_ppm_slope_mag limits its magnitude. Pure scalar algebra, [m].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_upwind_flux_neg`

<p class="rule-meta"><code>op:advection</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_upwind_flux_neg.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_upwind_flux_neg.esm)

Parameters: `u`, `dt`, `dy`, `I_dx`, `h_m2`, `h_m1`, `h_0`, `h_p1`, `h_p2`, `m_m2`, `m_m1`, `m_0`, `m_p1`, `m_p2`, `h_min_recon`

<div class="math-block">$$\operatorname{mom6\_ppm\_face\_integral}(I_{\mathrm{dx}} = I_{\mathrm{dx}}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dy}, h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{back}} = \operatorname{mom6\_ppm\_limited\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, \mathit{raw} = \operatorname{mom6\_ppm\_edge\_e}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}}), \mathit{scale} = \operatorname{mom6\_ppm\_cell\_scale}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}})), h_{\mathrm{face}} = \operatorname{mom6\_ppm\_limited\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, \mathit{raw} = \operatorname{mom6\_ppm\_edge\_w}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}), \mathit{scale} = \operatorname{mom6\_ppm\_cell\_scale}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}})), u = u)$$</div>

The MOM6 PPM layer flux through a face whose UPWIND cell is given by its 5-cell stencil, for u < 0: the upwind cell lies on the plus side of the face (the east/north neighbour for the east/north face of a cell), its LEFT (west/south) limited edge is the face edge h_face and its RIGHT edge the back edge h_back (zonal_flux_layer F90:958-964: uh = dy*u*(h_W(i+1) + CFL*(0.5*(h_E(i+1) - h_W(i+1)) + curv_3*(CFL - 1.5)))). The law is mom6_ppm_face_integral of the upwind cell (CFL = |u|*dt*I_dx, I_dx the upwind cell's 1/dx) with both edges the PPM_limit_pos-limited edges mom6_ppm_limited_edge of the cell: h_face = h_0 + scale*(mom6_ppm_edge_w - h_0), h_back = h_0 + scale*(mom6_ppm_edge_e - h_0), scale = mom6_ppm_cell_scale of the cell. The sibling mom6_ppm_upwind_flux_pos is the same law with the face/back edge roles swapped; the rules bind this one to the upwind cell's stencil, so the u > 0 and u < 0 branches of every flux rule are one reconstruction applied to two one-cell-shifted stencils. The stencil params h_m2..h_p2 / m_m2..m_p2 are positional relative to the upwind cell along the flux axis; dy is the face length, [m3 s-1].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_ppm_upwind_flux_pos`

<p class="rule-meta"><code>op:advection</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_ppm_upwind_flux_pos.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_ppm_upwind_flux_pos.esm)

Parameters: `u`, `dt`, `dy`, `I_dx`, `h_m2`, `h_m1`, `h_0`, `h_p1`, `h_p2`, `m_m2`, `m_m1`, `m_0`, `m_p1`, `m_p2`, `h_min_recon`

<div class="math-block">$$\operatorname{mom6\_ppm\_face\_integral}(I_{\mathrm{dx}} = I_{\mathrm{dx}}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dy}, h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{back}} = \operatorname{mom6\_ppm\_limited\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, \mathit{raw} = \operatorname{mom6\_ppm\_edge\_w}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}), \mathit{scale} = \operatorname{mom6\_ppm\_cell\_scale}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}})), h_{\mathrm{face}} = \operatorname{mom6\_ppm\_limited\_edge}(h_{\mathrm{0}} = h_{\mathrm{0}}, \mathit{raw} = \operatorname{mom6\_ppm\_edge\_e}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}}), \mathit{scale} = \operatorname{mom6\_ppm\_cell\_scale}(h_{\mathrm{0}} = h_{\mathrm{0}}, h_{\mathrm{m1}} = h_{\mathrm{m1}}, h_{\mathrm{m2}} = h_{\mathrm{m2}}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{\mathrm{p1}}, h_{\mathrm{p2}} = h_{\mathrm{p2}}, m_{\mathrm{0}} = m_{\mathrm{0}}, m_{\mathrm{m1}} = m_{\mathrm{m1}}, m_{\mathrm{m2}} = m_{\mathrm{m2}}, m_{\mathrm{p1}} = m_{\mathrm{p1}}, m_{\mathrm{p2}} = m_{\mathrm{p2}})), u = u)$$</div>

The MOM6 PPM layer flux through a face whose UPWIND cell is given by its 5-cell stencil, for u > 0: the upwind cell lies on the minus side of the face (the cell itself for the east/north face of a cell), its RIGHT (east/north) limited edge is the face edge h_face and its LEFT edge the back edge h_back (zonal_flux_layer F90:952-957: uh = dy*u*(h_E(i) + CFL*(0.5*(h_W(i) - h_E(i)) + curv_3*(CFL - 1.5)))). The law is mom6_ppm_face_integral of the upwind cell (CFL = |u|*dt*I_dx, I_dx the upwind cell's 1/dx) with both edges the PPM_limit_pos-limited edges mom6_ppm_limited_edge of the cell: h_face = h_0 + scale*(mom6_ppm_edge_e - h_0), h_back = h_0 + scale*(mom6_ppm_edge_w - h_0), scale = mom6_ppm_cell_scale of the cell. The sibling mom6_ppm_upwind_flux_neg is the same law with the face/back edge roles swapped; the rules bind this one to the upwind cell's stencil, so the u > 0 and u < 0 branches of every flux rule are one reconstruction applied to two one-cell-shifted stencils. The stencil params h_m2..h_p2 / m_m2..m_p2 are positional relative to the upwind cell along the flux axis; dy is the face length, [m3 s-1].

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)

### `mom6_q_corner_sum`

<p class="rule-meta"><code>op:sum</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_q_corner_sum.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_q_corner_sum.esm)

Parameters: `x_s`, `x_n`, `x_w`, `x_e`

<div class="math-block">$$\left(x_{\mathrm{s}} + x_{\mathrm{n}}\right) + \left(x_{\mathrm{w}} + x_{\mathrm{e}}\right)$$</div>

The q-cell corner sum of MOM6 (MOM_CoriolisAdv.F90 CorAdCalc :542, hArea_q = (hArea_u(I,j) + hArea_u(I,j+1)) + (hArea_v(i,J) + hArea_v(i+1,J))): sum = (x_s + x_n) + (x_w + x_e), the two u-face values south and north of the q-point plus the two v-face values west and east of it, in the Fortran's grouping. mom6_cgrid_sum_q binds it at every q-point; a face beyond the domain is a land face whose value is 0 (Area_h = 0 on land makes hArea_u/hArea_v vanish there). Pure scalar algebra.

**References**

- Sadourny, R. (1975). The dynamics of finite-difference models of the shallow-water equations. Journal of the Atmospheric Sciences, 32(4), 680-689. [`doi:10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2`](https://doi.org/10.1175/1520-0469%281975%29032%3C0680%3ATDOFDM%3E2.0.CO%3B2)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

### `mom6_sadourny_pvf`

<p class="rule-meta"><code>op:pvf</code> <code>interior-only (match-less)</code></p>

Source: [`grids/mom6_cgrid/stencils/mom6_sadourny_pvf.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/stencils/mom6_sadourny_pvf.esm)

Parameters: `q_a`, `f_a1`, `f_a2`, `q_b`, `f_b1`, `f_b2`, `I_d`

<div class="math-block">$$0.25 \left(\left(q_{\mathrm{a}} \left(f_{\mathrm{a1}} + f_{\mathrm{a2}}\right) + q_{\mathrm{b}} \left(f_{\mathrm{b1}} + f_{\mathrm{b2}}\right)\right) I_{\mathrm{d}}\right)$$</div>

The SADOURNY75_ENERGY potential-vorticity flux at one C-grid velocity face (MOM_CoriolisAdv.F90 CorAdCalc :724-726 for CAu, :951-953 for CAv): pvf = 0.25*(((q_a*(f_a1 + f_a2)) + (q_b*(f_b1 + f_b2)))*I_d), the two q-points a/b at the ends of the face, each times the sum of the two transverse layer volume fluxes adjacent to it, scaled by the face's 1/length I_d. ONE law for both faces: mom6_cgrid_pvf_u binds q_a = q(I,J), f_a = vh(i+1,J), vh(i,J), q_b = q(I,J-1), f_b = vh(i,J-1), vh(i+1,J-1), I_d = IdxCu(I,j); mom6_cgrid_pvf_v binds q_a = q(I-1,J), f_a = uh(I-1,j), uh(I-1,j+1), q_b = q(I,J), f_b = uh(I,j), uh(I,j+1), I_d = IdyCv(i,J) and negates the result. At a solid wall the rules bind the fluxes through the closed wall faces (and through land faces beyond the domain) to 0, so a q-point term whose two fluxes are wall fluxes vanishes, as MOM6's zero wall fluxes make it. Pure scalar algebra, [L T-2].

**References**

- Sadourny, R. (1975). The dynamics of finite-difference models of the shallow-water equations. Journal of the Atmospheric Sciences, 32(4), 680-689. [`doi:10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2`](https://doi.org/10.1175/1520-0469%281975%29032%3C0680%3ATDOFDM%3E2.0.CO%3B2)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `mom6_cgrid_avg_cell`

<p class="rule-meta"><code>op:avg</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_avg_cell.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_avg_cell.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_avg\_cell}(x_{\mathrm{u}}, x_{\mathrm{v}})\) — match pattern:

```json
{"op": "mom6_cgrid_avg_cell", "args": ["x_u", "x_v"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>x_u</code> is a bare field shaped <code>[i, j, k]</code>; <code>x_v</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits\operatorname{mom6\_cell\_face\_mean}(x_{\mathrm{e}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{n}} = x_{\mathrm{v}}[i, j, k], x_{\mathrm{s}} = x_{\mathrm{v}}[i, j - 1, k], x_{\mathrm{w}} = x_{\mathrm{u}}[i - 1, j, k])\right]_{2 \le i \le N_{\mathrm{I}}, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits\operatorname{mom6\_cell\_face\_mean}(x_{\mathrm{e}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{n}} = x_{\mathrm{v}}[i, j, k], x_{\mathrm{s}} = x_{\mathrm{v}}[i, j - 1, k], x_{\mathrm{w}} = 0)\right]_{1 \le i \le 1, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits\operatorname{mom6\_cell\_face\_mean}(x_{\mathrm{e}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{n}} = x_{\mathrm{v}}[i, j, k], x_{\mathrm{s}} = 0, x_{\mathrm{w}} = x_{\mathrm{u}}[i - 1, j, k])\right]_{2 \le i \le N_{\mathrm{I}}, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [1, 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits\operatorname{mom6\_cell\_face\_mean}(x_{\mathrm{e}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{n}} = x_{\mathrm{v}}[i, j, k], x_{\mathrm{s}} = 0, x_{\mathrm{w}} = 0)\right]_{1 \le i \le 1, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [1, 1] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The face-to-cell average on the mom6_cgrid: avg_cell(i,j,k) = 0.25*((x_u(i,j,k) + x_u(i-1,j,k)) + (x_v(i,j,k) + x_v(i,j-1,k))), the quarter-sum of the four face values of a cell, with x_u [i,j,k] the east-face field (x_u(i-1,j) is the cell's west face) and x_v [i,j,k] the north-face field (x_v(i,j-1) is its south face). Every region applies the one law mom6_cell_face_mean with region bindings. KE ASSUMPTION: the average is UNWEIGHTED. It reproduces the KE_ARAKAWA kinetic energy of MOM_CoriolisAdv.F90 gradKE :1290-1297, KE = (((areaCu(I,j)*u(I,j)^2) + (areaCu(I-1,j)*u(I-1,j)^2)) + ((areaCv(i,J)*v(i,J)^2) + (areaCv(i,J-1)*v(i,J-1)^2)))*0.25*IareaT, only if the consumer passes the area-weighted squares x_u = areaCu*u^2 and x_v = areaCv*v^2 and multiplies the result by IareaT itself (KE = IareaT*avg_cell(eu2, nv2) in the momentum fixture); with plain u^2/v^2 it is the KE_SIMPLE-style unweighted mean, not KE_ARAKAWA. SOLID-WALL BOUNDARIES: the west face of cell i = 1 and the south face of cell j = 1 are the closed domain walls, where the velocity (hence any face value built from it) is 0, the same zero-wall convention as mom6_cgrid_div_x/y, so those slots are bound to 0: avg_cell(1,j,k) = 0.25*((x_u(1,j,k) + 0) + (x_v(1,j,k) + x_v(1,j-1,k))), avg_cell(i,1,k) = 0.25*((x_u(i,1,k) + x_u(i-1,1,k)) + (x_v(i,1,k) + 0)), both at the (1,1) corner, which is the Fortran's grouping with u(0,j) = v(i,0) = 0. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_avg_cell(x_u, x_v); both params must be declared [i,j,k] fields (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 gradKE :1290-1297, the KE_ARAKAWA cell sum; the face conventions of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_avg_u`

<p class="rule-meta"><code>op:avg</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_avg_u.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_avg_u.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_avg\_u}(X)\) — match pattern:

```json
{"op": "mom6_cgrid_avg_u", "args": ["X"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>X</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits0.5 \left(X_{i, j, k} + X_{i + 1, j, k}\right)\right]_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits X_{N_{\mathrm{I}}, j, k}\right]_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal cell-to-face average on the mom6_cgrid: the arithmetic mean of a cell field X [i,j,k] across an east face, avg_u(i,j,k) = 0.5*(X(i,j,k) + X(i+1,j,k)) (the MOM6 u-cell average, e.g. the hArea_u = 0.5*(Area_h(i,j)*h(i,j) + Area_h(i+1,j)*h(i+1,j)) of MOM_CoriolisAdv.F90 CorAdCalc and the h_c + h_e face-thickness sums of MOM_PressureForce_FV.F90). SOLID-WALL BOUNDARY: the east face at the domain edge (i = N_I) is a wall with a zero-gradient halo (no normal variation of the cell fields through the wall), so the i = N_I region is the EXPLICIT mirror wall avg_u(N_I,j,k) = 0.5*(X(N_I,j,k) + X(N_I,j,k)) = X(N_I,j,k) — the mirrored neighbor is the cell itself, NOT a clamped read of a nonexistent cell. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_avg_u(X); X is a param whose shape must be declared [i,j,k] (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc hArea_u/hArea_v :348-357; the east-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_avg_v`

<p class="rule-meta"><code>op:avg</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_avg_v.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_avg_v.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_avg\_v}(X)\) — match pattern:

```json
{"op": "mom6_cgrid_avg_v", "args": ["X"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>X</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits0.5 \left(X_{i, j, k} + X_{i, j + 1, k}\right)\right]_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits X_{i, N_{\mathrm{J}}, k}\right]_{1 \le i \le N_{\mathrm{I}}, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional cell-to-face average on the mom6_cgrid, the y mirror of mom6_cgrid_avg_u: the arithmetic mean of a cell field X [i,j,k] across a north face, avg_v(i,j,k) = 0.5*(X(i,j,k) + X(i,j+1,k)) (the MOM6 v-cell average, e.g. the hArea_v = 0.5*(Area_h(i,j)*h(i,j) + Area_h(i,j+1)*h(i,j+1)) of MOM_CoriolisAdv.F90 CorAdCalc :348-357 and the h_c + h_n face-thickness sums of MOM_PressureForce_FV.F90). SOLID-WALL BOUNDARY: the north face at the domain edge (j = N_J) is a wall with a zero-gradient halo, so the j = N_J region is the EXPLICIT mirror wall avg_v(i,N_J,k) = 0.5*(X(i,N_J,k) + X(i,N_J,k)) = X(i,N_J,k) — the mirrored neighbor is the cell itself, NOT a clamped read of a nonexistent cell. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_avg_v(X); X is a param whose shape must be declared [i,j,k] (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc hArea_u/hArea_v :348-357; the north-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_bound_coriolis_u`

<p class="rule-meta"><code>op:clamp</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_bound_coriolis_u.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_bound_coriolis_u.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_bound\_coriolis\_u}(\mathit{pvf\_u}, \mathit{av}, v, \mathit{KEx})\) — match pattern:

```json
{"op": "mom6_cgrid_bound_coriolis_u", "args": ["pvf_u", "av", "v", "KEx"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>KEx</code> is a bare field shaped <code>[i, j, k]</code>; <code>av</code> is a bare field shaped <code>[i, j, k]</code>; <code>pvf_u</code> is a bare field shaped <code>[i, j, k]</code>; <code>v</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = \mathit{av}_{i, j, k} v_{i + 1, j, k}, \mathit{f2} = \mathit{av}_{i, j, k} v_{i, j, k}, \mathit{f3} = \mathit{av}_{i, j - 1, k} v_{i + 1, j - 1, k}, \mathit{f4} = \mathit{av}_{i, j - 1, k} v_{i, j - 1, k}, x = \mathit{pvf\_u}_{i, j, k}) - \mathit{KEx}_{i, j, k}\right)\right]^{\{\mathit{pvf\_u}, \mathit{av}, v, \mathit{KEx}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = 0, \mathit{f2} = \mathit{av}_{i, j, k} v_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = \mathit{av}_{i, j - 1, k} v_{i, j - 1, k}, x = \mathit{pvf\_u}_{i, j, k}) - \mathit{KEx}_{i, j, k}\right)\right]^{\{\mathit{pvf\_u}, \mathit{av}, v, \mathit{KEx}\}}_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = \mathit{av}_{i, j, k} v_{i + 1, j, k}, \mathit{f2} = \mathit{av}_{i, j, k} v_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = 0, x = \mathit{pvf\_u}_{i, j, k}) - \mathit{KEx}_{i, j, k}\right)\right]^{\{\mathit{pvf\_u}, \mathit{av}, v, \mathit{KEx}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = 0, \mathit{f2} = \mathit{av}_{i, j, k} v_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = 0, x = \mathit{pvf\_u}_{i, j, k}) - \mathit{KEx}_{i, j, k}\right)\right]^{\{\mathit{pvf\_u}, \mathit{av}, v, \mathit{KEx}\}}_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, 1] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The bound_coriolis east-face Coriolis + momentum-advection acceleration of the mom6_cgrid: CAu(i,j,k) = mom6_coriolis_clamp4(pvf_u(i,j,k), fv1, fv2, fv3, fv4) - KEx(i,j,k), the MOM_CoriolisAdv.F90 CorAdCalc bound_Coriolis law (the fv1..fv4 bounds and the min/max clamp at :897-908, then the -d(KE)/dx term at :911-914): fv1 = av(i,j,k)*v(i+1,j,k), fv2 = av(i,j,k)*v(i,j,k), fv3 = av(i,j-1,k)*v(i+1,j-1,k), fv4 = av(i,j-1,k)*v(i,j-1,k), the av*v products at the two q-points at the ends of the east face (av the absolute vorticity, v the north-face velocities of the cells on each side of the face), with the raw potential-vorticity flux pvf_u (the sibling mom6_cgrid_pvf_u) clamped to their range and the zonal KE gradient KEx subtracted. Every region applies the one clamp template, each product written once. SOLID-WALL BOUNDARIES (a basin closed by land, MOM6's zero wall velocities): at j = 1 the velocities v(i,0)/v(i+1,0) sit on the closed south wall, so fv3 = fv4 = 0 (bound to 0) and the clamp range of the real south-row faces includes 0, as MOM6's min/max over fv1..fv4 does. At i = N_I (the east wall face itself) v(N_I+1,*) are velocities of the land cell beyond, so fv1 = fv3 = 0 (the literal MOM6 value at a wall face, whose velocity is 0). The SE corner applies both. (Before this revision the wall regions copied the edge cell's own av/v into the out-of-domain slots, dropping 0 from the range on the whole j = 1 row of real faces.) The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_bound_coriolis_u(pvf_u, av, v, KEx); all four params must be declared [i,j,k] fields (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc bound_Coriolis fv1..fv4 :897-908 and the -d(KE)/dx term :911-914; the face conventions of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_bound_coriolis_v`

<p class="rule-meta"><code>op:clamp</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_bound_coriolis_v.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_bound_coriolis_v.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_bound\_coriolis\_v}(\mathit{pvf\_v}, \mathit{av}, u, \mathit{KEy})\) — match pattern:

```json
{"op": "mom6_cgrid_bound_coriolis_v", "args": ["pvf_v", "av", "u", "KEy"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>KEy</code> is a bare field shaped <code>[i, j, k]</code>; <code>av</code> is a bare field shaped <code>[i, j, k]</code>; <code>pvf_v</code> is a bare field shaped <code>[i, j, k]</code>; <code>u</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = -\mathit{av}_{i, j, k} u_{i, j + 1, k}, \mathit{f2} = -\mathit{av}_{i, j, k} u_{i, j, k}, \mathit{f3} = -\mathit{av}_{i - 1, j, k} u_{i - 1, j + 1, k}, \mathit{f4} = -\mathit{av}_{i - 1, j, k} u_{i - 1, j, k}, x = \mathit{pvf\_v}_{i, j, k}) - \mathit{KEy}_{i, j, k}\right)\right]^{\{\mathit{pvf\_v}, \mathit{av}, u, \mathit{KEy}\}}_{2 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = -\mathit{av}_{i, j, k} u_{i, j + 1, k}, \mathit{f2} = -\mathit{av}_{i, j, k} u_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = 0, x = \mathit{pvf\_v}_{i, j, k}) - \mathit{KEy}_{i, j, k}\right)\right]^{\{\mathit{pvf\_v}, \mathit{av}, u, \mathit{KEy}\}}_{1 \le i \le 1, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = 0, \mathit{f2} = -\mathit{av}_{i, j, k} u_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = -\mathit{av}_{i - 1, j, k} u_{i - 1, j, k}, x = \mathit{pvf\_v}_{i, j, k}) - \mathit{KEy}_{i, j, k}\right)\right]^{\{\mathit{pvf\_v}, \mathit{av}, u, \mathit{KEy}\}}_{2 \le i \le N_{\mathrm{I}}, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\left(\operatorname{mom6\_coriolis\_clamp4}(\mathit{f1} = 0, \mathit{f2} = -\mathit{av}_{i, j, k} u_{i, j, k}, \mathit{f3} = 0, \mathit{f4} = 0, x = \mathit{pvf\_v}_{i, j, k}) - \mathit{KEy}_{i, j, k}\right)\right]^{\{\mathit{pvf\_v}, \mathit{av}, u, \mathit{KEy}\}}_{1 \le i \le 1, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The bound_coriolis north-face Coriolis + momentum-advection acceleration of the mom6_cgrid: CAv(i,j,k) = mom6_coriolis_clamp4(pvf_v(i,j,k), fu1, fu2, fu3, fu4) - KEy(i,j,k), the MOM_CoriolisAdv.F90 CorAdCalc bound_Coriolis law (the fu1..fu4 bounds and the min/max clamp at :1139-1150, then the -d(KE)/dy term at :1153-1156): fu1 = -av(i,j,k)*u(i,j+1,k), fu2 = -av(i,j,k)*u(i,j,k), fu3 = -av(i-1,j,k)*u(i-1,j+1,k), fu4 = -av(i-1,j,k)*u(i-1,j,k), the -av*u products at the two q-points at the ends of the north face (av the absolute vorticity, u the east-face velocities of the cells on each side of the face), with the raw potential-vorticity flux pvf_v (the sibling mom6_cgrid_pvf_v) clamped to their range and the meridional KE gradient KEy subtracted. Every region applies the one clamp template (the same as mom6_cgrid_bound_coriolis_u), each product written once. SOLID-WALL BOUNDARIES (a basin closed by land, MOM6's zero wall velocities): at i = 1 the velocities u(0,j)/u(0,j+1) sit on the closed west wall, so fu3 = fu4 = 0 (bound to 0) and the clamp range of the real west-column faces includes 0, as MOM6's min/max over fu1..fu4 does. At j = N_J (the north wall face itself) u(*,N_J+1) are velocities of the land row beyond, so fu1 = fu3 = 0 (the literal MOM6 value at a wall face, whose velocity is 0). The NW corner applies both. (Before this revision the wall regions copied the edge cell's own av/u into the out-of-domain slots, dropping 0 from the range on the whole i = 1 column of real faces.) The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_bound_coriolis_v(pvf_v, av, u, KEy); all four params must be declared [i,j,k] fields (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc bound_Coriolis fu1..fu4 :1139-1150 and the -d(KE)/dy term :1153-1156; the face conventions of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_diff_x`

<p class="rule-meta"><code>op:diff</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_diff_x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_diff_x.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_diff\_x}(X)\) — match pattern:

```json
{"op": "mom6_cgrid_diff_x", "args": ["X"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>X</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits\left(X_{i + 1, j, k} - X_{i, j, k}\right)\right]_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal face forward difference on the mom6_cgrid: diff_x(i,j,k) = X(i+1,j,k) - X(i,j,k), the east-neighbor-minus-cell difference of a cell field X [i,j,k] lowered onto the east faces. It is the difference-building block of the finite-volume pressure-gradient force of MOM_PressureForce_FV.F90 PressureForce_FV_Bouss (:1797-1805): the consumer assembles the east-face pressure-difference group dP_x = -diff_x(P_int) + diff_x(h)*intx_pa - diff_x(e_bot)*intx_dpa*Z_to_H from the P_int/h/e_bot cell fields and the intx_pa/intx_dpa east-face interface integrals, keeping the Fortran's exact g1 + (g2f - g2s) grouping (verified bit-exact against the dump), and the face-thickness denominator h_c + h_e of the PFu acceleration from 2*mom6_cgrid_avg_u(h) (exact since the 0.5 factor is a power of two). SOLID-WALL BOUNDARY: the east face at the domain edge (i = N_I) is a wall with a zero-gradient halo (no normal variation of the cell fields through the wall), so the i = N_I region is the EXPLICIT mirror wall diff_x(N_I,j,k) = X(N_I,j,k) - X(N_I,j,k) = 0 — the mirrored neighbor is the cell itself, NOT a clamped read of a nonexistent cell. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_diff_x(X); X is a param whose shape must be declared [i,j,k] (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_PressureForce_FV.F90 PressureForce_FV_Bouss g1+g2 at :1797-1805; the east-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_diff_y`

<p class="rule-meta"><code>op:diff</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_diff_y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_diff_y.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_diff\_y}(X)\) — match pattern:

```json
{"op": "mom6_cgrid_diff_y", "args": ["X"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>X</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits\left(X_{i, j + 1, k} - X_{i, j, k}\right)\right]_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [1, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional face forward difference on the mom6_cgrid, the y mirror of mom6_cgrid_diff_x: diff_y(i,j,k) = X(i,j+1,k) - X(i,j,k), the north-neighbor-minus-cell difference of a cell field X [i,j,k] lowered onto the north faces. It is the difference-building block of the finite-volume pressure-gradient force of MOM_PressureForce_FV.F90 PressureForce_FV_Bouss (:1797-1805): the consumer assembles the north-face pressure-difference group dP_y = -diff_y(P_int) + diff_y(h)*inty_pa - diff_y(e_bot)*inty_dpa*Z_to_H from the P_int/h/e_bot cell fields and the inty_pa/inty_dpa north-face interface integrals, keeping the Fortran's exact g1 + (g2f - g2s) grouping (verified bit-exact against the dump), and the face-thickness denominator h_c + h_n of the PFv acceleration from 2*mom6_cgrid_avg_v(h) (exact since the 0.5 factor is a power of two). SOLID-WALL BOUNDARY: the north face at the domain edge (j = N_J) is a wall with a zero-gradient halo, so the j = N_J region is the EXPLICIT mirror wall diff_y(i,N_J,k) = X(i,N_J,k) - X(i,N_J,k) = 0 (the y mirror of the diff_x explicit-wall statement). The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_diff_y(X); X is a param whose shape must be declared [i,j,k] (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_PressureForce_FV.F90 PressureForce_FV_Bouss g1+g2 at :1797-1805; the north-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_div_x`

<p class="rule-meta"><code>op:div</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 10</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_div_x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_div_x.esm)

**Rewrites** \(\frac{\partial\mathit{uh}}{\partial x}\) at priority 10 — match pattern:

```json
{"op": "D", "args": ["uh"], "wrt": "x"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>uh</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\mathit{IareaT}_{i, j} \left(\mathit{uh}_{i, j, k} - \mathit{uh}_{i - 1, j, k}\right)\right]^{\{\mathit{uh}, \mathit{IareaT}\}}_{2 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\mathit{IareaT}_{i, j} \mathit{uh}_{i, j, k}\right]^{\{\mathit{uh}, \mathit{IareaT}\}}_{1 \le i \le 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal C-grid flux-form divergence rewrite target D(uh, wrt: x) on the mom6_cgrid: lowers the spatial derivative of the east-face layer volume flux field uh [i,j,k] to the cell divergence div_x(i,j,k) = (uh(i,j,k) - uh(i-1,j,k)) * IareaT(i,j) (the MOM6 continuity convergence (uh(i,j) - uh(i-1,j))*IareaT of continuity_zonal_convergence, MOM_continuity_PPM.F90:1085-1094, with uh the east-face flux of each cell and uh(i-1,j) the east-face flux of the west neighbor, i.e. the west face of cell (i,j)). SOLID-WALL BOUNDARIES: the west edge of cell i = 1 is a wall whose face flux is identically zero, so the i = 1 region is the EXPLICIT wall divergence div_x(1,j,k) = uh(1,j,k) * IareaT(1,j) — NOT the clamped difference (uh(1) - uh(1))*IareaT = 0 of a naive halo closure (the bug it replaces: a clamped west boundary silently zeroes the divergence of the first cell column). The east wall is carried by the flux rule, which emits the zero face flux at the east domain edge (i = N_I), so the interior region formula extends to i = N_I unchanged: div_x(N_I,j,k) = (0 - uh(N_I-1,j,k))*IareaT(N_I,j). The rule matches any D(uh, x) whose operand is a declared [i,j,k] field (esm-spec 9.6.1); IareaT is a free name resolved in the consuming model (the free-name contract of grids/mom6_cgrid/grid.esm). Pure telescoping difference of the face fluxes; the PPM reconstruction lives in the sibling mom6_cgrid_flux_x rule.

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

**Conformance fixtures**

- simulation: [`mom6_cgrid_flux_ppm`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_flux_ppm) — status `active`

### `mom6_cgrid_div_y`

<p class="rule-meta"><code>op:div</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 10</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_div_y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_div_y.esm)

**Rewrites** \(\frac{\partial\mathit{vh}}{\partial y}\) at priority 10 — match pattern:

```json
{"op": "D", "args": ["vh"], "wrt": "y"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>vh</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\mathit{IareaT}_{i, j} \left(\mathit{vh}_{i, j, k} - \mathit{vh}_{i, j - 1, k}\right)\right]^{\{\mathit{vh}, \mathit{IareaT}\}}_{1 \le i \le N_{\mathrm{I}}, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\mathit{IareaT}_{i, j} \mathit{vh}_{i, j, k}\right]^{\{\mathit{vh}, \mathit{IareaT}\}}_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, 1] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional C-grid flux-form divergence rewrite target D(vh, wrt: y) on the mom6_cgrid, the y mirror of mom6_cgrid_div_x: lowers the spatial derivative of the north-face layer volume flux field vh [i,j,k] to the cell divergence div_y(i,j,k) = (vh(i,j,k) - vh(i,j-1,k)) * IareaT(i,j) (the MOM6 continuity convergence (vh(i,j) - vh(i,j-1))*IareaT of continuity_merid_convergence, MOM_continuity_PPM.F90:1564-1573, with vh(i,j) the north-face flux of each cell and vh(i,j-1) the north-face flux of the south neighbor, i.e. the south face of cell (i,j)). SOLID-WALL BOUNDARIES: the south edge of cell j = 1 is a wall whose face flux is identically zero, so the j = 1 region is the EXPLICIT wall divergence div_y(i,1,k) = vh(i,1,k) * IareaT(i,1) — NOT the clamped difference (vh(1) - vh(1))*IareaT = 0 (the bug it replaces, the y mirror of the div_x fix). The north wall is carried by the flux rule, which emits the zero face flux at the north domain edge (j = N_J), so the interior region formula extends to j = N_J unchanged: div_y(i,N_J,k) = (0 - vh(i,N_J-1,k))*IareaT(i,N_J). The rule matches any D(vh, y) whose operand is a declared [i,j,k] field; IareaT is a free name resolved in the consuming model. Pure telescoping difference of the face fluxes; the PPM reconstruction lives in the sibling mom6_cgrid_flux_y rule.

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

**Conformance fixtures**

- simulation: [`mom6_cgrid_flux_ppm`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_flux_ppm) — status `active`

### `mom6_cgrid_flux_x`

<p class="rule-meta"><code>op:advection</code> <code>order:3</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 10</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_flux_x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_flux_x.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_flux\_x}(h, u_{\mathrm{E}}, \mathit{dy\_Cu\_E}, \mathit{mask}, h_{\mathrm{min\_recon}}, \mathit{dt}, \mathit{IdxT}, \mathit{IdxT\_e})\) at priority 10 — match pattern:

```json
{"op": "mom6_cgrid_flux_x", "args": ["h", "u_E", "dy_Cu_E", "mask", "h_min_recon", "dt", "IdxT", "IdxT_e"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>IdxT</code> is a bare field shaped <code>[i, j]</code>; <code>IdxT_e</code> is a bare field shaped <code>[i, j]</code>; <code>dy_Cu_E</code> is a bare field shaped <code>[i, j]</code>; <code>h</code> is a bare field shaped <code>[i, j, k]</code>; <code>mask</code> is a bare field shaped <code>[i, j]</code>; <code>u_E</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\begin{cases}\operatorname{mom6\_ppm\_upwind\_flux\_pos}(I_{\mathrm{dx}} = \mathit{IdxT}_{i, j}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dy\_Cu\_E}_{i, j}, h_{\mathrm{0}} = h_{i, j, k}, h_{\mathrm{m1}} = h_{\max(1, i - 1), j, k}, h_{\mathrm{m2}} = h_{\max(1, i - 2), j, k}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{i + 1, j, k}, h_{\mathrm{p2}} = h_{\min(N_{\mathrm{I}}, i + 2), j, k}, m_{\mathrm{0}} = \mathit{mask}_{i, j}, m_{\mathrm{m1}} = \mathit{mask}_{\max(1, i - 1), j}, m_{\mathrm{m2}} = \mathit{mask}_{\max(1, i - 2), j}, m_{\mathrm{p1}} = \mathit{mask}_{i + 1, j}, m_{\mathrm{p2}} = \mathit{mask}_{\min(N_{\mathrm{I}}, i + 2), j}, u = u_{\mathrm{E}}[i, j, k]) &amp; \text{if } u_{\mathrm{E}}[i, j, k] &gt; 0 \\ \operatorname{mom6\_ppm\_upwind\_flux\_neg}(I_{\mathrm{dx}} = \mathit{IdxT\_e}_{i, j}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dy\_Cu\_E}_{i, j}, h_{\mathrm{0}} = h_{i + 1, j, k}, h_{\mathrm{m1}} = h_{i, j, k}, h_{\mathrm{m2}} = h_{\max(1, i - 1), j, k}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{\min(N_{\mathrm{I}}, i + 2), j, k}, h_{\mathrm{p2}} = h_{\min(N_{\mathrm{I}}, i + 3), j, k}, m_{\mathrm{0}} = \mathit{mask}_{i + 1, j}, m_{\mathrm{m1}} = \mathit{mask}_{i, j}, m_{\mathrm{m2}} = \mathit{mask}_{\max(1, i - 1), j}, m_{\mathrm{p1}} = \mathit{mask}_{\min(N_{\mathrm{I}}, i + 2), j}, m_{\mathrm{p2}} = \mathit{mask}_{\min(N_{\mathrm{I}}, i + 3), j}, u = u_{\mathrm{E}}[i, j, k]) &amp; \text{if } u_{\mathrm{E}}[i, j, k] &lt; 0 \\ 0 &amp; \text{otherwise}\end{cases}\right]^{\{h, u_{\mathrm{E}}, \mathit{dy\_Cu\_E}, \mathit{mask}, h_{\mathrm{min\_recon}}, \mathit{dt}, \mathit{IdxT}, \mathit{IdxT\_e}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal PPM layer-flux lowering rule of the mom6_cgrid: matches the user rewrite target mom6_cgrid_flux_x(h, u_E, dy_Cu_E, mask, h_min_recon, dt, IdxT, IdxT_e) (the east-face layer volume flux of every cell) and emits the [i,j,k] makearray of the MOM6 zonal flux law (zonal_flux_layer, MOM_continuity_PPM.F90:950-966, VOL_CFL=False, por_face_areaU = 1): uh(i,j,k) = dy_Cu_E(i,j) * u_E(i,j,k) * (the PPM parabola integral of the UPWIND cell's PPM_limit_pos-limited edges). The sign of the face velocity selects the upwind cell, and the two branches apply the one per-cell reconstruction to two one-cell-shifted stencils: u_E > 0 upwinds the cell itself, mom6_ppm_upwind_flux_pos bound to the stencil i-2..i+2 with I_dx = IdxT(i,j) (F90:952-957: CFL = u*dt*IdxT(i,j), uh = dy*u*(h_E(i) + CFL*(0.5*(h_W(i) - h_E(i)) + curv_3*(CFL - 1.5)))); u_E < 0 upwinds the east neighbour, mom6_ppm_upwind_flux_neg bound to the stencil i-1..i+3 with I_dx = IdxT_e(i,j) (F90:958-964: CFL = -u*dt*IdxT(i+1,j), uh = dy*u*(h_W(i+1) + CFL*(0.5*(h_E(i+1) - h_W(i+1)) + curv_3*(CFL - 1.5)))); u_E == 0 gives the zero flux. The reconstruction (center slopes, raw edges, PPM_limit_pos scale, limited edges, parabola integral) lives in the stencil laws, so the meridional rule mom6_cgrid_flux_y applies the same two templates along j. LAND CLOSURE: the stencil index expressions CLAMP at the domain edges (max(1, i-2), max(1, i-1) low; min(N_I, i+2), min(N_I, i+3) high), and the clamp reproduces MOM6's land-masked closure exactly, as if the cells outside the domain were land (mask2dT = 0): a clamped thickness equals the masked-neighbour replacement h_im1 = h(i) of PPM_reconstruction_x (F90:2430-2431), and the clamped center slope of an edge cell and of a ghost cell is identically zero (the Lin bound 2*min(dMx, dMn) vanishes when two of the three stencil values coincide), which is the masked zero slope of F90:2397-2398. The faces next to the walls (i = 1, i = N_I-1) are therefore the MOM6 faces of a basin closed by land, pinned numerically by problems/mom6_cgrid_flux_ppm.esm. The EAST DOMAIN WALL is an explicit region: the face at i = N_I (the east edge of the last cell column, between the last ocean cell and the land beyond) carries the zero flux of a solid wall whatever u_E holds there (MOM6 sets u = 0 on a land face); the convergence (uh(i,j) - uh(i-1,j))*IareaT then telescopes to -uh(N_I-1)*IareaT in the sibling mom6_cgrid_div_x. The consuming model spells uh = mom6_cgrid_flux_x(h, u_E, dy_Cu_E, mask, h_min_recon, dt, IdxT, IdxT_e) with h [i,j,k], u_E [i,j,k], the horizontal metrics dy_Cu_E/IdxT/IdxT_e [i,j], mask [i,j] and the scalars h_min_recon/dt (esm-spec 9.6.1: non-param strings in the match are literal per-variable selectors, so the rule fires only on THIS operator with exactly these argument variables).

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)
- Lin, S.-J. (1994). A nearly conservative flux-reconstruction scheme. Monthly Weather Review, 122(5), 994-1011. (The monotonic slope constraint, Eq. B2.)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

**Conformance fixtures**

- simulation: [`mom6_cgrid_flux_ppm`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_flux_ppm) — status `active`

### `mom6_cgrid_flux_y`

<p class="rule-meta"><code>op:advection</code> <code>order:3</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 10</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_flux_y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_flux_y.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_flux\_y}(h, v_{\mathrm{N}}, \mathit{dx\_Cv\_N}, \mathit{mask}, h_{\mathrm{min\_recon}}, \mathit{dt}, \mathit{IdyT}, \mathit{IdyT\_n})\) at priority 10 — match pattern:

```json
{"op": "mom6_cgrid_flux_y", "args": ["h", "v_N", "dx_Cv_N", "mask", "h_min_recon", "dt", "IdyT", "IdyT_n"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>IdyT</code> is a bare field shaped <code>[i, j]</code>; <code>IdyT_n</code> is a bare field shaped <code>[i, j]</code>; <code>dx_Cv_N</code> is a bare field shaped <code>[i, j]</code>; <code>h</code> is a bare field shaped <code>[i, j, k]</code>; <code>mask</code> is a bare field shaped <code>[i, j]</code>; <code>v_N</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\begin{cases}\operatorname{mom6\_ppm\_upwind\_flux\_pos}(I_{\mathrm{dx}} = \mathit{IdyT}_{i, j}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dx\_Cv\_N}_{i, j}, h_{\mathrm{0}} = h_{i, j, k}, h_{\mathrm{m1}} = h_{i, \max(1, j - 1), k}, h_{\mathrm{m2}} = h_{i, \max(1, j - 2), k}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{i, j + 1, k}, h_{\mathrm{p2}} = h_{i, \min(N_{\mathrm{J}}, j + 2), k}, m_{\mathrm{0}} = \mathit{mask}_{i, j}, m_{\mathrm{m1}} = \mathit{mask}_{i, \max(1, j - 1)}, m_{\mathrm{m2}} = \mathit{mask}_{i, \max(1, j - 2)}, m_{\mathrm{p1}} = \mathit{mask}_{i, j + 1}, m_{\mathrm{p2}} = \mathit{mask}_{i, \min(N_{\mathrm{J}}, j + 2)}, u = v_{\mathrm{N}}[i, j, k]) &amp; \text{if } v_{\mathrm{N}}[i, j, k] &gt; 0 \\ \operatorname{mom6\_ppm\_upwind\_flux\_neg}(I_{\mathrm{dx}} = \mathit{IdyT\_n}_{i, j}, \mathit{dt} = \mathit{dt}, \mathit{dy} = \mathit{dx\_Cv\_N}_{i, j}, h_{\mathrm{0}} = h_{i, j + 1, k}, h_{\mathrm{m1}} = h_{i, j, k}, h_{\mathrm{m2}} = h_{i, \max(1, j - 1), k}, h_{\mathrm{min\_recon}} = h_{\mathrm{min\_recon}}, h_{\mathrm{p1}} = h_{i, \min(N_{\mathrm{J}}, j + 2), k}, h_{\mathrm{p2}} = h_{i, \min(N_{\mathrm{J}}, j + 3), k}, m_{\mathrm{0}} = \mathit{mask}_{i, j + 1}, m_{\mathrm{m1}} = \mathit{mask}_{i, j}, m_{\mathrm{m2}} = \mathit{mask}_{i, \max(1, j - 1)}, m_{\mathrm{p1}} = \mathit{mask}_{i, \min(N_{\mathrm{J}}, j + 2)}, m_{\mathrm{p2}} = \mathit{mask}_{i, \min(N_{\mathrm{J}}, j + 3)}, u = v_{\mathrm{N}}[i, j, k]) &amp; \text{if } v_{\mathrm{N}}[i, j, k] &lt; 0 \\ 0 &amp; \text{otherwise}\end{cases}\right]^{\{h, v_{\mathrm{N}}, \mathit{dx\_Cv\_N}, \mathit{mask}, h_{\mathrm{min\_recon}}, \mathit{dt}, \mathit{IdyT}, \mathit{IdyT\_n}\}}_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [1, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional PPM layer-flux lowering rule of the mom6_cgrid, the y mirror of mom6_cgrid_flux_x: matches the user rewrite target mom6_cgrid_flux_y(h, v_N, dx_Cv_N, mask, h_min_recon, dt, IdyT, IdyT_n) (the north-face layer volume flux of every cell) and emits the [i,j,k] makearray of the MOM6 meridional flux law (merid_flux_layer, MOM_continuity_PPM.F90:1852-1872, VOL_CFL=False, por_face_areaV = 1): vh(i,j,k) = dx_Cv_N(i,j) * v_N(i,j,k) * (the PPM parabola integral of the UPWIND cell's PPM_limit_pos-limited edges), applying the same two per-cell templates as the zonal rule along j (PPM_reconstruction_y is PPM_reconstruction_x with the axes swapped). The sign of the face velocity selects the upwind cell: v_N > 0 upwinds the cell itself, mom6_ppm_upwind_flux_pos bound to the stencil j-2..j+2 with I_dx = IdyT(i,j) (F90:1854-1861: CFL = v*dt*IdyT(i,j), vh = dx*v*(h_N(i,j) + CFL*(0.5*(h_S(i,j) - h_N(i,j)) + curv_3*(CFL - 1.5)))); v_N < 0 upwinds the north neighbour, mom6_ppm_upwind_flux_neg bound to the stencil j-1..j+3 with I_dx = IdyT_n(i,j) (F90:1862-1869: CFL = -v*dt*IdyT(i,j+1), vh = dx*v*(h_S(i,j+1) + CFL*(0.5*(h_N(i,j+1) - h_S(i,j+1)) + curv_3*(CFL - 1.5)))); v_N == 0 gives the zero flux. LAND CLOSURE: the stencil index expressions CLAMP at the domain edges (max(1, j-2), max(1, j-1) low; min(N_J, j+2), min(N_J, j+3) high), which reproduces MOM6's land-masked closure exactly, as if the cells outside the domain were land (the clamped thickness is the masked-neighbour replacement and the clamped slope of an edge or ghost cell is identically zero; see mom6_cgrid_flux_x). The NORTH DOMAIN WALL is an explicit region: the face at j = N_J (the north edge of the last cell row) carries the zero flux of a solid wall whatever v_N holds there; the convergence (vh(i,j) - vh(i,j-1))*IareaT then telescopes to -vh(i,N_J-1)*IareaT in the sibling mom6_cgrid_div_y. The consuming model spells vh = mom6_cgrid_flux_y(h, v_N, dx_Cv_N, mask, h_min_recon, dt, IdyT, IdyT_n) with h [i,j,k], v_N [i,j,k], the horizontal metrics dx_Cv_N/IdyT/IdyT_n [i,j], mask [i,j] and the scalars h_min_recon/dt (esm-spec 9.6.1: non-param strings in the match are literal per-variable selectors, so the rule fires only on THIS operator with exactly these argument variables).

**References**

- Colella, P., & Woodward, P. R. (1984). The Piecewise Parabolic Method (PPM) for gas-dynamical simulations. Journal of Computational Physics, 54(1), 174-201. [`doi:10.1016/0021-9991(84)90143-8`](https://doi.org/10.1016/0021-9991%2884%2990143-8)
- Lin, S.-J. (1994). A nearly conservative flux-reconstruction scheme. Monthly Weather Review, 122(5), 994-1011. (The monotonic slope constraint, Eq. B2.)
- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)

**Conformance fixtures**

- simulation: [`mom6_cgrid_flux_ppm`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_flux_ppm) — status `active`

### `mom6_cgrid_grad_x`

<p class="rule-meta"><code>op:grad</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 20</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_grad_x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_grad_x.esm)

**Rewrites** \(\frac{\partial\mathit{KE}}{\partial x}\) at priority 20 — match pattern:

```json
{"op": "D", "args": ["KE"], "wrt": "x"}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\mathit{IdxCu}_{i, j} \left(\mathit{KE}_{i + 1, j, k} - \mathit{KE}_{i, j, k}\right)\right]^{\{\mathit{KE}, \mathit{IdxCu}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal C-grid kinetic-energy-gradient rewrite target D(KE, wrt: x) on the mom6_cgrid: lowers the spatial derivative of the cell kinetic energy KE [i,j,k] to the east-face zonal gradient KEx(i,j,k) = (KE(i+1,j,k) - KE(i,j,k)) * IdxCu(i,j) (the KE_ARAKAWA gradient term -d(KE)/dx of the Coriolis + momentum-advection acceleration, MOM_CoriolisAdv.F90 gradKE :1249-1397 with the face difference (KE(i+1,j) - KE(i,j))*G%IdxCu(I,j) of the KEx loop :1300-1320; IdxCu(i,j) = 1/dxCu is the east-face 1/length so the gradient lives on the east faces). SOLID-WALL BOUNDARY: the east face at the domain edge (i = N_I) is a wall whose normal kinetic-energy gradient is identically zero, so the i = N_I region is the EXPLICIT free-slip wall value KEx(N_I,j,k) = 0 — NOT the clamped difference (KE(N_I) - KE(N_I))*IdxCu = 0 of a naive halo closure (the clamped form happens to evaluate to zero here, but it hides the wall as if the cell duplicated; this rule states the wall explicitly). The rule is a GROUND per-variable selector (esm-spec 9.6.1): KE is NOT a parameter, so the match arg "KE" is a literal that fires only on the derivative of the declared variable named KE, with priority 20 to outrank the wildcard mom6_cgrid_div_x (priority 10) when both are imported. SHAPE GUARD: esm-spec 9.6.1 admits `where` keys only for declared params (a non-param key is apply_expression_template_invalid_declaration), so a ground selector cannot carry a where shape constraint; making KE a param to constrain it would turn the rule into a wildcard D(*, x) over [i,j,k] fields that outranks mom6_cgrid_div_x on every operand. The guard is therefore the literal name: the rule fires only on D(KE, x) of a variable named KE, and a consumer importing it must declare KE as the [i,j,k] field this rule assumes (the momentum fixtures and problems/mom6_cgrid_momentum.esm do). Unlike the vort rules, D(KE, x) IS a derivative: the face difference is scaled by IdxCu (KEx of gradKE, MOM_CoriolisAdv.F90:1404). KE and IdxCu are free names resolved in the consuming model (the free-name contract of grids/mom6_cgrid/grid.esm); the body's KE is the matched operand.

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc gradKE :1249-1397, KEx at :1300-1320; the east-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_grad_y`

<p class="rule-meta"><code>op:grad</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 20</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_grad_y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_grad_y.esm)

**Rewrites** \(\frac{\partial\mathit{KE}}{\partial y}\) at priority 20 — match pattern:

```json
{"op": "D", "args": ["KE"], "wrt": "y"}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\mathit{IdyCv}_{i, j} \left(\mathit{KE}_{i, j + 1, k} - \mathit{KE}_{i, j, k}\right)\right]^{\{\mathit{KE}, \mathit{IdyCv}\}}_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [1, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional C-grid kinetic-energy-gradient rewrite target D(KE, wrt: y) on the mom6_cgrid, the y mirror of mom6_cgrid_grad_x: lowers the spatial derivative of the cell kinetic energy KE [i,j,k] to the north-face meridional gradient KEy(i,j,k) = (KE(i,j+1,k) - KE(i,j,k)) * IdyCv(i,j) (the KE_ARAKAWA gradient term -d(KE)/dy of the Coriolis + momentum-advection acceleration, MOM_CoriolisAdv.F90 gradKE :1249-1397 with the face difference (KE(i,j+1) - KE(i,j))*G%IdyCv(i,J) of the KEy loop; IdyCv(i,j) = 1/dyCv is the north-face 1/length so the gradient lives on the north faces). SOLID-WALL BOUNDARY: the north face at the domain edge (j = N_J) is a wall whose normal kinetic-energy gradient is identically zero, so the j = N_J region is the EXPLICIT free-slip wall value KEy(i,N_J,k) = 0 (the y mirror of the grad_x explicit-wall statement). The rule is a GROUND per-variable selector (esm-spec 9.6.1): KE is NOT a parameter, so the match arg "KE" is a literal that fires only on the derivative of the declared variable named KE, with priority 20 to outrank the wildcard mom6_cgrid_div_y (priority 10) when both are imported. SHAPE GUARD: esm-spec 9.6.1 admits `where` keys only for declared params (a non-param key is apply_expression_template_invalid_declaration), so a ground selector cannot carry a where shape constraint; making KE a param to constrain it would turn the rule into a wildcard D(*, y) over [i,j,k] fields that outranks mom6_cgrid_div_y on every operand. The guard is therefore the literal name: the rule fires only on D(KE, y) of a variable named KE, and a consumer importing it must declare KE as the [i,j,k] field this rule assumes (the momentum fixtures and problems/mom6_cgrid_momentum.esm do). Unlike the vort rules, D(KE, y) IS a derivative: the face difference is scaled by IdyCv (KEy of gradKE, MOM_CoriolisAdv.F90:1409). KE and IdyCv are free names resolved in the consuming model (the free-name contract of grids/mom6_cgrid/grid.esm); the body's KE is the matched operand.

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc gradKE :1249-1397, KEy at :1322-1342; the north-face convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_mask2dBu`

<p class="rule-meta"><code>op:mask</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_mask2dBu.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_mask2dBu.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_mask2dBu}(\mathit{mask})\) — match pattern:

```json
{"op": "mom6_cgrid_mask2dBu", "args": ["mask"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>mask</code> is a bare field shaped <code>[i, j]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\limits\left(\mathit{mask}_{i, j} \mathit{mask}_{i + 1, j}\right) \left(\mathit{mask}_{i, j + 1} \mathit{mask}_{i + 1, j + 1}\right)\right]_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}} - 1} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}} - 1] \\ 0 &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \\ 0 &amp; [1, N_{\mathrm{I}} - 1] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \end{array}\right]$$</div>

The q-point mask of the mom6_cgrid: the product of the masks of the four T cells around a q-point, mask2dBu(i,j) = mask(i,j) * mask(i+1,j) * mask(i,j+1) * mask(i+1,j+1), where mask [i,j] is the 0/1 ocean-cell indicator and q-point (i,j) sits at the NE corner of cell (i,j). This is the MOM6 vertex mask of MOM_grid_initialize.F90:1294 (G%mask2dBu(I,J) = (G%mask2dCu(I,j) * G%mask2dCu(I,j+1)) * (G%mask2dCv(i,J) * G%mask2dCv(i+1,J)), the product of the four surrounding face masks, which equals the product of the four surrounding cell masks on the interior), used by CorAdCalc to zero the relative vorticity at land q-points. SOLID-WALL BOUNDARIES: at the east domain edge (i = N_I) or the north domain edge (j = N_J) there is no cell beyond the wall, so mask2dBu is EXPLICITLY zero there (the q-point touches land — a cell mask of 0) — the i = N_I and j = N_J regions are the wall value 0.0, NOT clamped reads of nonexistent cells. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_mask2dBu(mask); the param mask must be a declared [i,j] field (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_grid_initialize.F90:1294 (mask2dBu from the face masks); the q-point convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_pvf_u`

<p class="rule-meta"><code>op:pvf</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_pvf_u.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_pvf_u.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_pvf\_u}(q, \mathit{vh}, \mathit{IdxCu})\) — match pattern:

```json
{"op": "mom6_cgrid_pvf_u", "args": ["q", "vh", "IdxCu"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>IdxCu</code> is a bare field shaped <code>[i, j]</code>; <code>q</code> is a bare field shaped <code>[i, j, k]</code>; <code>vh</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdxCu}_{i, j}, f_{\mathrm{a1}} = \mathit{vh}_{i + 1, j, k}, f_{\mathrm{a2}} = \mathit{vh}_{i, j, k}, f_{\mathrm{b1}} = \mathit{vh}_{i, j - 1, k}, f_{\mathrm{b2}} = \mathit{vh}_{i + 1, j - 1, k}, q_{\mathrm{a}} = q_{i, j, k}, q_{\mathrm{b}} = q_{i, j - 1, k})\right]^{\{q, \mathit{vh}, \mathit{IdxCu}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdxCu}_{i, j}, f_{\mathrm{a1}} = 0, f_{\mathrm{a2}} = \mathit{vh}_{i, j, k}, f_{\mathrm{b1}} = \mathit{vh}_{i, j - 1, k}, f_{\mathrm{b2}} = 0, q_{\mathrm{a}} = q_{i, j, k}, q_{\mathrm{b}} = q_{i, j - 1, k})\right]^{\{q, \mathit{vh}, \mathit{IdxCu}\}}_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 2 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [2, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdxCu}_{i, j}, f_{\mathrm{a1}} = \mathit{vh}_{i + 1, j, k}, f_{\mathrm{a2}} = \mathit{vh}_{i, j, k}, f_{\mathrm{b1}} = 0, f_{\mathrm{b2}} = 0, q_{\mathrm{a}} = q_{i, j, k}, q_{\mathrm{b}} = 0)\right]^{\{q, \mathit{vh}, \mathit{IdxCu}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdxCu}_{i, j}, f_{\mathrm{a1}} = 0, f_{\mathrm{a2}} = \mathit{vh}_{i, j, k}, f_{\mathrm{b1}} = 0, f_{\mathrm{b2}} = 0, q_{\mathrm{a}} = q_{i, j, k}, q_{\mathrm{b}} = 0)\right]^{\{q, \mathit{vh}, \mathit{IdxCu}\}}_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 1 \le j \le 1, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, 1] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The raw east-face potential-vorticity flux of the SADOURNY75_ENERGY Coriolis scheme on the mom6_cgrid: pvf_u(i,j,k) = 0.25 * ((q(i,j,k)*(vh(i+1,j,k) + vh(i,j,k))) + (q(i,j-1,k)*(vh(i,j-1,k) + vh(i+1,j-1,k)))) * IdxCu(i,j), the CAu of MOM_CoriolisAdv.F90 CorAdCalc :724-726 before the bound_coriolis clamp and the -KEx term, with q(i,j) the layer potential vorticity at the q-point at the NE corner of cell (i,j), vh the north-face layer volume fluxes of the two cells on each side of the east face (vh(i,j)/vh(i+1,j) the north faces of the cell and its east neighbour, vh(i,j-1)/vh(i+1,j-1) their south faces) and IdxCu(i,j) the east-face 1/length. Every region applies the one law mom6_sadourny_pvf with region bindings. SOLID-WALL BOUNDARIES (a basin closed by land, MOM6's zero wall fluxes): at j = 1 the south faces vh(i,0)/vh(i+1,0) are the closed south wall, whose flux is 0, so the q(i,0) term vanishes (q_b, f_b1, f_b2 are bound to 0; the wall q-point is not stored and is multiplied only by zero fluxes) and pvf_u(i,1) = 0.25*((q(i,1)*(vh(i+1,1) + vh(i,1))) + 0)*IdxCu(i,1), an interior east face whose value MOM6 computes from the north-side term alone. At i = N_I the east face is the closed east wall itself and vh(N_I+1,*) are north faces of the land cell beyond it, flux 0, so f_a1 and f_b2 are bound to 0 (the literal MOM6 value at the wall face; the wall face's velocity is 0, so it never accelerates anything). The SE corner applies both. (Before this revision the wall regions copied the edge cell's own q/vh into the out-of-domain slots, double-counting the interior term on the whole j = 1 row of real faces.) The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_pvf_u(q, vh, IdxCu); the where clause requires q and vh to be [i,j,k] fields and IdxCu an [i,j] field (esm-spec 9.6.1).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc SADOURNY75_ENERGY CAu at :722-727; the face conventions of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_pvf_v`

<p class="rule-meta"><code>op:pvf</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_pvf_v.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_pvf_v.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_pvf\_v}(q, \mathit{uh}, \mathit{IdyCv})\) — match pattern:

```json
{"op": "mom6_cgrid_pvf_v", "args": ["q", "uh", "IdyCv"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>IdyCv</code> is a bare field shaped <code>[i, j]</code>; <code>q</code> is a bare field shaped <code>[i, j, k]</code>; <code>uh</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum-\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdyCv}_{i, j}, f_{\mathrm{a1}} = \mathit{uh}_{i - 1, j, k}, f_{\mathrm{a2}} = \mathit{uh}_{i - 1, j + 1, k}, f_{\mathrm{b1}} = \mathit{uh}_{i, j, k}, f_{\mathrm{b2}} = \mathit{uh}_{i, j + 1, k}, q_{\mathrm{a}} = q_{i - 1, j, k}, q_{\mathrm{b}} = q_{i, j, k})\right]^{\{q, \mathit{uh}, \mathit{IdyCv}\}}_{2 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum-\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdyCv}_{i, j}, f_{\mathrm{a1}} = 0, f_{\mathrm{a2}} = 0, f_{\mathrm{b1}} = \mathit{uh}_{i, j, k}, f_{\mathrm{b2}} = \mathit{uh}_{i, j + 1, k}, q_{\mathrm{a}} = 0, q_{\mathrm{b}} = q_{i, j, k})\right]^{\{q, \mathit{uh}, \mathit{IdyCv}\}}_{1 \le i \le 1, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum-\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdyCv}_{i, j}, f_{\mathrm{a1}} = \mathit{uh}_{i - 1, j, k}, f_{\mathrm{a2}} = 0, f_{\mathrm{b1}} = \mathit{uh}_{i, j, k}, f_{\mathrm{b2}} = 0, q_{\mathrm{a}} = q_{i - 1, j, k}, q_{\mathrm{b}} = q_{i, j, k})\right]^{\{q, \mathit{uh}, \mathit{IdyCv}\}}_{2 \le i \le N_{\mathrm{I}}, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [2, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum-\operatorname{mom6\_sadourny\_pvf}(I_{\mathrm{d}} = \mathit{IdyCv}_{i, j}, f_{\mathrm{a1}} = 0, f_{\mathrm{a2}} = 0, f_{\mathrm{b1}} = \mathit{uh}_{i, j, k}, f_{\mathrm{b2}} = 0, q_{\mathrm{a}} = 0, q_{\mathrm{b}} = q_{i, j, k})\right]^{\{q, \mathit{uh}, \mathit{IdyCv}\}}_{1 \le i \le 1, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, 1] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The raw north-face potential-vorticity flux of the SADOURNY75_ENERGY Coriolis scheme on the mom6_cgrid: pvf_v(i,j,k) = -0.25 * ((q(i-1,j,k)*(uh(i-1,j,k) + uh(i-1,j+1,k))) + (q(i,j,k)*(uh(i,j,k) + uh(i,j+1,k)))) * IdyCv(i,j), the CAv of MOM_CoriolisAdv.F90 CorAdCalc :951-953 before the bound_coriolis clamp and the -KEy term, with q the layer potential vorticity at the q-points at the two ends of the north face (q(i-1,j) west, q(i,j) east), uh the east-face layer volume fluxes adjacent to them (uh(i-1,j)/uh(i-1,j+1) the west faces of the cell and its north neighbour, uh(i,j)/uh(i,j+1) their east faces) and IdyCv(i,j) the north-face 1/length. Every region applies the one law mom6_sadourny_pvf (the same template as mom6_cgrid_pvf_u) with region bindings, negated. SOLID-WALL BOUNDARIES (a basin closed by land, MOM6's zero wall fluxes): at i = 1 the west faces uh(0,j)/uh(0,j+1) are the closed west wall, flux 0, so the q(0,j) term vanishes (q_a, f_a1, f_a2 are bound to 0; the wall q-point is not stored and is multiplied only by zero fluxes) and pvf_v(1,j) = -0.25*(0 + (q(1,j)*(uh(1,j) + uh(1,j+1))))*IdyCv(1,j), an interior north face whose value MOM6 computes from the east-side term alone. At j = N_J the north face is the closed north wall itself and uh(*,N_J+1) are east faces of the land row beyond it, flux 0, so f_a2 and f_b2 are bound to 0 (the literal MOM6 value at the wall face, whose velocity is 0). The NW corner applies both. (Before this revision the wall regions copied the edge cell's own q/uh into the out-of-domain slots, double-counting the interior term on the whole i = 1 column of real faces.) The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_pvf_v(q, uh, IdyCv); the where clause requires q and uh to be [i,j,k] fields and IdyCv an [i,j] field (esm-spec 9.6.1).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc SADOURNY75_ENERGY CAv at :949-954; the face conventions of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_sum_q`

<p class="rule-meta"><code>op:sum</code> <code>order:0</code> <code>bc:wall</code> <code>axes:i,j,k</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_sum_q.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_sum_q.esm)

**Rewrites** \(\operatorname{mom6\_cgrid\_sum\_q}(x_{\mathrm{u}}, x_{\mathrm{v}})\) — match pattern:

```json
{"op": "mom6_cgrid_sum_q", "args": ["x_u", "x_v"]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>x_u</code> is a bare field shaped <code>[i, j, k]</code>; <code>x_v</code> is a bare field shaped <code>[i, j, k]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\operatorname{mom6\_q\_corner\_sum}(x_{\mathrm{e}} = x_{\mathrm{v}}[i + 1, j, k], x_{\mathrm{n}} = x_{\mathrm{u}}[i, j + 1, k], x_{\mathrm{s}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{w}} = x_{\mathrm{v}}[i, j, k])\right]^{\{x_{\mathrm{u}}, x_{\mathrm{v}}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits\operatorname{mom6\_q\_corner\_sum}(x_{\mathrm{e}} = 0, x_{\mathrm{n}} = x_{\mathrm{u}}[i, j + 1, k], x_{\mathrm{s}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{w}} = x_{\mathrm{v}}[i, j, k])\right]_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ \left[\sum\operatorname{mom6\_q\_corner\_sum}(x_{\mathrm{e}} = x_{\mathrm{v}}[i + 1, j, k], x_{\mathrm{n}} = 0, x_{\mathrm{s}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{w}} = x_{\mathrm{v}}[i, j, k])\right]^{\{x_{\mathrm{u}}, x_{\mathrm{v}}\}}_{1 \le i \le N_{\mathrm{I}} - 1, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ \left[\sum\limits\operatorname{mom6\_q\_corner\_sum}(x_{\mathrm{e}} = 0, x_{\mathrm{n}} = 0, x_{\mathrm{s}} = x_{\mathrm{u}}[i, j, k], x_{\mathrm{w}} = x_{\mathrm{v}}[i, j, k])\right]_{N_{\mathrm{I}} \le i \le N_{\mathrm{I}}, N_{\mathrm{J}} \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The q-point corner sum on the mom6_cgrid: sum_q(i,j,k) = (x_u(i,j,k) + x_u(i,j+1,k)) + (x_v(i,j,k) + x_v(i+1,j,k)), the two east-face values south and north of the q-point at the NE corner of cell (i,j) plus the two north-face values west and east of it, the MOM6 q-cell layer volume hArea_q = (hArea_u(I,j) + hArea_u(I,j+1)) + (hArea_v(i,J) + hArea_v(i+1,J)) of MOM_CoriolisAdv.F90 CorAdCalc :542 (with hArea_u/hArea_v of :361-366). Every region applies the one law mom6_q_corner_sum with region bindings. WALL q-POINTS (a basin closed by land): the q-points at i = N_I and j = N_J lie on the east/north walls, where the face beyond the domain (x_v(N_I+1,j), x_u(i,N_J+1)) belongs to a land cell whose Area_h = 0 makes hArea_v/hArea_u vanish, so those slots are bound to 0, MOM6's value. (Before this revision they copied the edge face's own value.) A wall q-point never reaches an interior face: it multiplies only the zero wall-face fluxes in mom6_cgrid_pvf_u/v. The rule is a named-op rule (esm-spec 9.6.8), applied once per consumer equation as mom6_cgrid_sum_q(x_u, x_v); both params must be declared [i,j,k] fields (esm-spec 9.6.1 where clause).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc hArea_q :541-545; the q-point convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_vort_x`

<p class="rule-meta"><code>op:vort</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 20</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_vort_x.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_vort_x.esm)

**Rewrites** \(\frac{\partial v}{\partial x}\) at priority 20 — match pattern:

```json
{"op": "D", "args": ["v"], "wrt": "x"}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\left(\mathit{dyCv}_{i + 1, j} v_{i + 1, j, k} - \mathit{dyCv}_{i, j} v_{i, j, k}\right)\right]^{\{v, \mathit{dyCv}\}}_{1 \le i \le N_{\mathrm{I}} - 1, 1 \le j \le N_{\mathrm{J}}, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}} - 1] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [N_{\mathrm{I}}, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The zonal C-grid vorticity rewrite target D(v, wrt: x) on the mom6_cgrid: lowers the spatial derivative of the meridional face velocity field v [i,j,k] (the north-face velocity of each T cell) to the q-point zonal derivative of the meridional volume transport dvdx(i,j,k) = v(i+1,j,k)*dyCv(i+1,j) - v(i,j,k)*dyCv(i,j) (the forward difference dvdx(I,J) = v(i+1,J,k)*G%dyCv(i+1,J) - v(i,J,k)*G%dyCv(i,J) of MOM_CoriolisAdv.F90:322-323, with v(i+1,j) the north-face velocity of the east neighbor cell, i.e. the v at the east side of the q-point, and dyCv(i+1,j) its face length; the vorticity is second-order accurate with free-slip boundaries). SOLID-WALL BOUNDARY: the q-points at the east domain edge (i = N_I) lie on the east wall. MOM6 still forms dvdx there (= -v(N_I,j)*dyCv(N_I,j), the land v beyond being 0), but its free-slip relative vorticity multiplies dvdx by mask2dBu = 0 on wall q-points (MOM_CoriolisAdv.F90:526), so no wall value ever reaches rel_vort; this rule states the free-slip result explicitly, dvdx(N_I,j,k) = 0. UNDIVIDED: D(v, x) lowers to the undivided circulation difference v(i+1,j)*dyCv(i+1,j) - v(i,j)*dyCv(i,j) [m2 s-1], NOT a derivative dv/dx: no 1/dx factor is applied, because MOM6 divides the q-cell circulation (dvdx - dudy) by the q-cell area once, through IareaBu (:526). A consumer that wants dv/dx must not import this rule. The rule is a GROUND per-variable selector (esm-spec 9.6.1): v is NOT a parameter, so the match arg "v" is a literal that fires only on the derivative of the declared variable named v — mixed-scheme per-variable selection on one axis — and carries priority 20 so it outranks the wildcard sibling mom6_cgrid_div_x (priority 10, whose uh IS a parameter), which matches any [i,j,k] operand, whenever both are imported. SHAPE GUARD: esm-spec 9.6.1 admits `where` keys only for declared params (a non-param key is apply_expression_template_invalid_declaration), so a ground selector cannot carry a where shape constraint; making v a param to constrain it would turn the rule into a wildcard D(*, x) over [i,j,k] fields that outranks mom6_cgrid_div_x on every operand. The guard is therefore the literal name: the rule fires only on D(v, x) of a variable named v, and a consumer importing it must declare v as the [i,j,k] field this rule assumes (the momentum fixtures and problems/mom6_cgrid_momentum.esm do). The body is a pure forward difference of the face transport — no reconstruction (the momentum PV advection is centered, PV_ADV_CENTERED). v and dyCv are free names resolved in the consuming model (the free-name contract of grids/mom6_cgrid/grid.esm); the body's v is the matched operand, so a consumer must bind the derivative target's operand to a declared field named v (esm-spec 9.6.1, the same bare-variable discipline as a where shape constraint).

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc, dvdx at :322-323 (Stokes_VF = .false. branch); the q-point convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

### `mom6_cgrid_vort_y`

<p class="rule-meta"><code>op:vort</code> <code>order:1</code> <code>bc:wall</code> <code>axes:i,j,k</code> <code>priority 20</code></p>

Source: [`grids/mom6_cgrid/rules/mom6_cgrid_vort_y.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/mom6_cgrid/rules/mom6_cgrid_vort_y.esm)

**Rewrites** \(\frac{\partial u}{\partial y}\) at priority 20 — match pattern:

```json
{"op": "D", "args": ["u"], "wrt": "y"}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \left[\sum\left(\mathit{dxCu}_{i, j + 1} u_{i, j + 1, k} - \mathit{dxCu}_{i, j} u_{i, j, k}\right)\right]^{\{u, \mathit{dxCu}\}}_{1 \le i \le N_{\mathrm{I}}, 1 \le j \le N_{\mathrm{J}} - 1, 1 \le k \le N_{\mathrm{K}}} &amp; [1, N_{\mathrm{I}}] \times [1, N_{\mathrm{J}} - 1] \times [1, N_{\mathrm{K}}] \\ 0 &amp; [1, N_{\mathrm{I}}] \times [N_{\mathrm{J}}, N_{\mathrm{J}}] \times [1, N_{\mathrm{K}}] \end{array}\right]$$</div>

The meridional C-grid vorticity rewrite target D(u, wrt: y) on the mom6_cgrid, the y mirror of mom6_cgrid_vort_x: lowers the spatial derivative of the zonal face velocity field u [i,j,k] (the east-face velocity of each T cell) to the q-point meridional derivative of the zonal volume transport dudy(i,j,k) = u(i,j+1,k)*dxCu(i,j+1) - u(i,j,k)*dxCu(i,j) (the forward difference dudy(I,J) = u(I,j+1,k)*G%dxCu(I,j+1) - u(I,j,k)*G%dxCu(I,j) of MOM_CoriolisAdv.F90:325-326, with u(i,j+1) the east-face velocity of the north neighbor cell, i.e. the u at the north side of the q-point, and dxCu(i,j+1) its face length). SOLID-WALL BOUNDARY: the q-points at the north domain edge (j = N_J) lie on the north wall; MOM6's free-slip relative vorticity multiplies dudy by mask2dBu = 0 there (MOM_CoriolisAdv.F90:526), so this rule states the free-slip result explicitly, dudy(i,N_J,k) = 0 (the y mirror of vort_x). UNDIVIDED: D(u, y) lowers to the undivided circulation difference u(i,j+1)*dxCu(i,j+1) - u(i,j)*dxCu(i,j) [m2 s-1], NOT a derivative du/dy; MOM6 divides the q-cell circulation by the q-cell area once, through IareaBu. The rule is a GROUND per-variable selector (esm-spec 9.6.1): u is NOT a parameter, so the match arg "u" is a literal that fires only on the derivative of the declared variable named u, with priority 20 to outrank the wildcard mom6_cgrid_div_y (priority 10) when both are imported; SHAPE GUARD: esm-spec 9.6.1 admits `where` keys only for declared params (a non-param key is apply_expression_template_invalid_declaration), so a ground selector cannot carry a where shape constraint; making u a param to constrain it would turn the rule into a wildcard D(*, y) over [i,j,k] fields that outranks mom6_cgrid_div_y on every operand. The guard is therefore the literal name: the rule fires only on D(u, y) of a variable named u, and a consumer importing it must declare u as the [i,j,k] field this rule assumes (the momentum fixtures and problems/mom6_cgrid_momentum.esm do); the body is a pure forward difference of the face transport (PV_ADV_CENTERED momentum PV advection). u and dxCu are free names resolved in the consuming model (the free-name contract of grids/mom6_cgrid/grid.esm); the body's u is the matched operand.

**References**

- Adcroft, A., et al. (2019). The GFDL global ocean and sea ice model OM4.0. Journal of Advances in Modeling Earth Systems, 11(10), 3167-3211. [`doi:10.1029/2019MS001726`](https://doi.org/10.1029/2019MS001726)
- MOM6 source: MOM_CoriolisAdv.F90 CorAdCalc, dudy at :325-326 (Stokes_VF = .false. branch); the q-point convention of the mom6_cgrid grid.esm.

**Conformance fixtures**

- simulation: [`mom6_cgrid_momentum`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/mom6_cgrid_momentum) — status `active`

