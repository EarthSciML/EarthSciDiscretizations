---
title: "cice_bgrid"
description: "Two-dimensional B-grid of the CICE sea-ice model (Hunke et al., Los Alamos Sea Ice Model), the horizontal mesh of the momentum/stress/transport components: `i` (x) and `j` (y) are the cell axes of N_I x N_J T cells (the grid cells that carry the ice state), `corner` is the 4-member axis of the cell-corner (staggered U-point) stress values (corner 1 = NE = the U point U(i,j) at the cell's north-east corner, 2 = NW, 3 = SW, 4 = SE), and `ncat` is the thickness-category axis of NCAT categories."
esd_kind: "grid"
families: ["cice_bgrid"]
ops: ["D", "advection", "strain-rate", "stress"]
bcs: ["halo-copy"]
tags: ["esd:grid", "family:cice_bgrid", "ndim:2", "stagger:b", "axis:i", "axis:j", "esd:stencil", "grid:cice_bgrid", "op:stress", "op:strain-rate", "op:advection", "esd:rule", "op:D", "axes:i,j,corner", "order:1", "bc:halo-copy", "spec:9.6.8", "axes:i,j", "axes:i,j,ncat", "form:flux"]
source: "grids/cice_bgrid/grid.esm"
---

Two-dimensional B-grid of the CICE sea-ice model (Hunke et al., Los Alamos Sea Ice Model), the horizontal mesh of the momentum/stress/transport components: `i` (x) and `j` (y) are the cell axes of N_I x N_J T cells (the grid cells that carry the ice state), `corner` is the 4-member axis of the cell-corner (staggered U-point) stress values (corner 1 = NE = the U point U(i,j) at the cell's north-east corner, 2 = NW, 3 = SW, 4 = SE), and `ncat` is the thickness-category axis of NCAT categories. N_I, N_J, NCAT are the esm-spec 9.7 metaparameters, bound at an import edge or the loader API; the consumer (EqESM stage2 cice_dynamics/cice_transport) binds N_I = 3, N_J = 5, NCAT = 5 for the traced 3x5 subdomain of the stage-1 dump. FREE-NAME CONTRACT: a consuming model supplies the real-valued horizontal geometry as its own variables with the names dxT, dyT (T-cell E-W and N-S dimensions), dxhy, dyhx (the 1/2 east-edge-length and north-edge-length differences), cxp, cyp, cxm, cym (the strain-rate metric coefficients), HTE, HTN (the east- and north-edge lengths of the T cells), uarear (1/U-cell area), and tarear (1/T-cell area); the templates below derive each metric from the minimal set (HTE, HTN, dxT, dyT, tarea) with the Fortran's exact algebra, and a consumer MAY either bind the metrics directly from its own state (the stage-2 pin tests override them from the dump) or derive them through these templates. METRIC DERIVATIONS (ice_dyn_shared.F90 init_dyn_shared): dxhy(i,j) = 0.5*(HTE(i,j) - HTE(i-1,j)) (F90:413-414), dyhx(i,j) = 0.5*(HTN(i,j) - HTN(i,j-1)) (F90:415-416), cyp(i,j) = 1.5*HTE(i,j) - 0.5*HTE(i-1,j) (F90:434), cxp(i,j) = 1.5*HTN(i,j) - 0.5*HTN(i,j-1) (F90:435), cym(i,j) = -(1.5*HTE(i-1,j) - 0.5*HTE(i,j)) (F90:436), cxm(i,j) = -(1.5*HTN(i,j-1) - 0.5*HTN(i,j)) (F90:437); and (ice_grid.F90 primary_grid_lengths_HTN/HTE): dxT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1)) (F90:1342), dyT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j)) (F90:1380), tarea = dxT*dyT (F90:682) with tarear = 1/tarea (F90:698), and uarea(i,j) = 0.25*(tarea(i,j) + tarea(i+1,j) + tarea(i,j+1) + tarea(i+1,j+1)) (F90:1643-1645) with uarear = 1/uarea (F90:703). The derivations read the west (i-1) and south (j-1) neighbors of the T cell; the Fortran halos HTE/HTN before deriving (ice_HaloUpdate with fillValue = 1.0, F90:418-428), and a consumer whose HTE/HTN carry a halo-copy closure at the i = 1 / j = 1 boundary reproduces the halo values exactly (the stage-2 traced subdomain supplies real values on all four sides of the pin cells, so the closures do not enter the pin tests). The B-grid spatial rules that lower the D operators of the momentum/stress/transport equations live in rules/ (strain_rates_divu/tension/shear for D(uvel,x)+D(vvel,y) etc., stress_divergence_strintx/strinty for D(sigma11,x)+D(sigma12,y) etc., transport_flux_divergence for D(c*u,x)+D(c*v,y)) with their scalar atoms in stencils/, following the latlon3d flux-rule precedent of the velocity inside the compound operator.

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `N_I` | integer | `1` | Number of T cells along the B-grid x axis (cell-based axis), [-]. |
| `N_J` | integer | `1` | Number of T cells along the B-grid y axis (cell-based axis), [-]. |
| `NCAT` | integer | `5` | Number of thickness categories of the transported ice state (namelist ncat = 5 default), [-]. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>sigma11</code> <code>sigma12</code> <code>sigma22</code> <code>uarear</code></p>

These names are supplied by the consuming model's scope at evaluation (the keyed-factor / free-name geometry contract); they are neither metaparameters (load-time integers) nor template parameters (bound per apply site).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `i` | interval | `N_I` |
| `j` | interval | `N_J` |
| `corner` | interval | `4` |
| `ncat` | interval | `NCAT` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `cice_cyp` | `HTE` | `sum[gi, gj] (1.5 * HTE[gi, gj] - 0.5 * HTE[gi - 1, gj]) where {gi in i, gj in j}` | Strain-rate metric coefficient cyp of the NE corner combination (ice_dyn_shared.F90:434): cyp(i,j) = 1.5*HTE(i,j) - 0.5*HTE(i-1,j), gathered over the [i,j] cell axes (aggregation indices gi/gj for inlining hygiene). A consuming model applies this to its HTE to derive its [i,j] cyp, or supplies cyp directly under the free-name contract. |
| `cice_cym` | `HTE` | `sum[gi, gj] (-(1.5 * HTE[gi - 1, gj] - 0.5 * HTE[gi, gj])) where {gi in i, gj in j}` | Strain-rate metric coefficient cym of the NW/SW corner combinations (ice_dyn_shared.F90:436): cym(i,j) = -(1.5*HTE(i-1,j) - 0.5*HTE(i,j)), gathered over the [i,j] cell axes. |
| `cice_cxp` | `HTN` | `sum[gi, gj] (1.5 * HTN[gi, gj] - 0.5 * HTN[gi, gj - 1]) where {gi in i, gj in j}` | Strain-rate metric coefficient cxp of the NE/NW corner combinations (ice_dyn_shared.F90:435): cxp(i,j) = 1.5*HTN(i,j) - 0.5*HTN(i,j-1), gathered over the [i,j] cell axes. |
| `cice_cxm` | `HTN` | `sum[gi, gj] (-(1.5 * HTN[gi, gj - 1] - 0.5 * HTN[gi, gj])) where {gi in i, gj in j}` | Strain-rate metric coefficient cxm of the SW/SE corner combinations (ice_dyn_shared.F90:437): cxm(i,j) = -(1.5*HTN(i,j-1) - 0.5*HTN(i,j)), gathered over the [i,j] cell axes. |
| `cice_dxhy` | `HTE` | `sum[gi, gj] (0.5 * (HTE[gi, gj] - HTE[gi - 1, gj])) where {gi in i, gj in j}` | Half the east-edge-length difference of the str combinations (ice_dyn_shared.F90:413-414): dxhy(i,j) = 0.5*(HTE(i,j) - HTE(i-1,j)), gathered over the [i,j] cell axes. |
| `cice_dyhx` | `HTN` | `sum[gi, gj] (0.5 * (HTN[gi, gj] - HTN[gi, gj - 1])) where {gi in i, gj in j}` | Half the north-edge-length difference of the str combinations (ice_dyn_shared.F90:415-416): dyhx(i,j) = 0.5*(HTN(i,j) - HTN(i,j-1)), gathered over the [i,j] cell axes. |
| `cice_dxT` | `HTN` | `sum[gi, gj] (0.5 * (HTN[gi, gj] + HTN[gi, gj - 1])) where {gi in i, gj in j}` | E-W dimension of the T cell (ice_grid.F90:1342): dxT(i,j) = 0.5*(HTN(i,j) + HTN(i,j-1)), gathered over the [i,j] cell axes. |
| `cice_dyT` | `HTE` | `sum[gi, gj] (0.5 * (HTE[gi, gj] + HTE[gi - 1, gj])) where {gi in i, gj in j}` | N-S dimension of the T cell (ice_grid.F90:1380): dyT(i,j) = 0.5*(HTE(i,j) + HTE(i-1,j)), gathered over the [i,j] cell axes. |
| `cice_tarea` | `dxT`, `dyT` | `sum[gi, gj] (dxT[gi, gj] * dyT[gi, gj]) where {gi in i, gj in j}` | T-cell area (ice_grid.F90:682): tarea(i,j) = dxT(i,j)*dyT(i,j), gathered over the [i,j] cell axes; tarear = 1/tarea (F90:698) is the reciprocal the consumer derives from it. |
| `cice_uarea` | `tarea` | `sum[gi, gj] (0.25 * (tarea[gi, gj] + tarea[gi + 1, gj] + tarea[gi, gj + 1] + tarea[gi + 1, gj + 1])) where {gi in i, gj in j}` | U-cell (corner cell) area from the four surrounding T-cell areas (ice_grid.F90:1643-1645): uarea(i,j) = 0.25*(tarea(i,j) + tarea(i+1,j) + tarea(i,j+1) + tarea(i+1,j+1)), gathered over the [i,j] U-point axes; uarear = 1/uarea (F90:703) is the reciprocal the consumer derives from it. This is the standard-B-grid form the Fortran uses away from the tripole fold (the special tripole form of F90:1641-1642 is out of scope for the traced subdomain). |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `str_csig`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_csig.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_csig.esm)

Parameters: `a1`, `a2`, `a3`, `a4`, `c1`, `c2`, `c3`

<div class="math-block">c1 * a1 + c2 * (a2 + a4) + c3 * a3</div>

Scalar atom of the CICE stress-combination algebra (ice_dyn_evp.F90:1667-1686, the csigp*/csigm*/csig12* corner averages): csig = c1*a1 + c2*(a2 + a4) + c3*a3. Each csig term of the Fortran is p111 times one corner value plus p055/p111 times the opposite-corner pair plus p027/p055 times the third corner; e.g. csigpne = p111*p1 + p055*(p2 + p4) + p027*p3 (F90:1667) binds a1 = p1, a2 = p2, a3 = p3, a4 = p4, c1 = p111, c2 = p055, c3 = p027, and csig12ne = p222*s1 + p111*(s2 + s4) + p055*s3 (F90:1683) binds c1 = p222, c2 = p111, c3 = p055. Pure combination algebra: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dx_cell_k1`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dx_cell_k1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dx_cell_k1.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-(str_pairsum_lin&lt;a1=p1, a2=p2, a3=p3, a4=p4, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=m1, a2=m2, a3=m3, a4=m4, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=s1, a2=s4, a3=s2, a4=s3, c=0.5, metric=dxT&gt;) + dxhy * (-str_csig&lt;a1=p1, a2=p2, a3=p3, a4=p4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m1, a2=m2, a3=m3, a4=m4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;) + dyhx * str_csig&lt;a1=s1, a2=s2, a3=s3, a4=s4, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,1), the NE-corner contribution to the x divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1727-1728): str1 = -strp_tmp - strm_tmp - str12ew + dxhy*(-csigpne + csigmne) + dyhx*csig12ne, with strp_tmp/strm_tmp = p25*dyT*(p333*ssigpn + p166*ssigps) etc. (F90:1723-1724, the p/m north-south pair sums), str12ew = p5*dxT*(p333*ssig12e + p166*ssig12w) (F90:1715), csigpne/csigmne = p111*p1 + p055*(p2+p4) + p027*p3 (F90:1696, 1701) and csig12ne = p222*s1 + p111*(s2+s4) + p055*s3 (F90:1706). The p/m/s params are the cell's four corner stressp/stressm/stress12 values (p_k = sigma11+sigma22, m_k = sigma11-sigma22, s_k = sigma12 at corner k), bound by the stress_divergence_strintx rule to the cells (i,j), (i+1,j), (i,j+1), (i+1,j+1) that surround the U point, and the dxT/dyT/dxhy/dyhx params are the metrics of the cell where the str value lives. Pure combination algebra built on the str_csig and str_pairsum_lin atoms: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dx_cell_k2`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dx_cell_k2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dx_cell_k2.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">str_pairsum_lin&lt;a1=p1, a2=p2, a3=p3, a4=p4, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=m1, a2=m2, a3=m3, a4=m4, c=0.25, metric=dyT&gt; + (-str_pairsum_lin&lt;a1=s2, a2=s3, a3=s1, a4=s4, c=0.5, metric=dxT&gt;) + dxhy * (-str_csig&lt;a1=p2, a2=p1, a3=p4, a4=p3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m2, a2=m1, a3=m4, a4=m3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;) + dyhx * str_csig&lt;a1=s2, a2=s1, a3=s4, a4=s3, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,2), the NW-corner contribution to the x divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1731-1732): str2 = strp_tmp + strm_tmp - str12we + dxhy*(-csigpnw + csigmnw) + dyhx*csig12nw, with strp_tmp/strm_tmp = p25*dyT*(p333*ssigpn + p166*ssigps) (F90:1723-1724), str12we = p5*dxT*(p333*ssig12w + p166*ssig12e) (F90:1716), csigpnw/csigmnw = p111*p2 + p055*(p1+p3) + p027*p4 (F90:1697, 1702) and csig12nw = p222*s2 + p111*(s1+s3) + p055*s4 (F90:1708). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the NW cell (i+1,j) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dx_cell_k3`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dx_cell_k3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dx_cell_k3.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-(str_pairsum_lin&lt;a1=p3, a2=p4, a3=p1, a4=p2, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=m3, a2=m4, a3=m1, a4=m2, c=0.25, metric=dyT&gt;) + str_pairsum_lin&lt;a1=s1, a2=s4, a3=s2, a4=s3, c=0.5, metric=dxT&gt; + dxhy * (-str_csig&lt;a1=p4, a2=p1, a3=p2, a4=p3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m4, a2=m1, a3=m2, a4=m3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;) + dyhx * str_csig&lt;a1=s4, a2=s1, a3=s2, a4=s3, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,3), the SE-corner contribution to the x divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1738-1739): str3 = -strp_tmp - strm_tmp + str12ew + dxhy*(-csigpse + csigmse) + dyhx*csig12se, with strp_tmp/strm_tmp = p25*dyT*(p333*ssigps + p166*ssigpn) (F90:1734-1735), str12ew = p5*dxT*(p333*ssig12e + p166*ssig12w) (F90:1715), csigpse/csigmse = p111*p4 + p055*(p1+p3) + p027*p2 (F90:1699, 1704) and csig12se = p222*s4 + p111*(s1+s3) + p055*s2 (F90:1712). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the SE cell (i,j+1) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dx_cell_k4`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dx_cell_k4.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dx_cell_k4.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">str_pairsum_lin&lt;a1=p3, a2=p4, a3=p1, a4=p2, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=m3, a2=m4, a3=m1, a4=m2, c=0.25, metric=dyT&gt; + str_pairsum_lin&lt;a1=s2, a2=s3, a3=s1, a4=s4, c=0.5, metric=dxT&gt; + dxhy * (-str_csig&lt;a1=p3, a2=p2, a3=p1, a4=p4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m3, a2=m2, a3=m1, a4=m4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;) + dyhx * str_csig&lt;a1=s3, a2=s2, a3=s1, a4=s4, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,4), the SW-corner contribution to the x divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1742-1743): str4 = strp_tmp + strm_tmp + str12we + dxhy*(-csigpsw + csigmsw) + dyhx*csig12sw, with strp_tmp/strm_tmp = p25*dyT*(p333*ssigps + p166*ssigpn) (F90:1734-1735), str12we = p5*dxT*(p333*ssig12w + p166*ssig12e) (F90:1716), csigpsw/csigmsw = p111*p3 + p055*(p2+p4) + p027*p1 (F90:1698, 1703) and csig12sw = p222*s3 + p111*(s2+s4) + p055*s1 (F90:1710). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the SW cell (i+1,j+1) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dy_cell_k5`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dy_cell_k5.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dy_cell_k5.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-(str_pairsum_lin&lt;a1=p1, a2=p4, a3=p2, a4=p3, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=s1, a2=s2, a3=s3, a4=s4, c=0.5, metric=dyT&gt;) + str_pairsum_lin&lt;a1=m1, a2=m4, a3=m2, a4=m3, c=0.25, metric=dxT&gt; + (-dyhx * (str_csig&lt;a1=p1, a2=p2, a3=p3, a4=p4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m1, a2=m2, a3=m3, a4=m4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;)) + dxhy * str_csig&lt;a1=s1, a2=s2, a3=s3, a4=s4, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,5), the NE-corner contribution to the y divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1752-1753): str5 = -strp_tmp + strm_tmp - str12ns - dyhx*(csigpne + csigmne) + dxhy*csig12ne, with strp_tmp/strm_tmp = p25*dxT*(p333*ssigpe + p166*ssigpw) (F90:1748-1749, the p/m east-west pair sums), str12ns = p5*dyT*(p333*ssig12n + p166*ssig12s) (F90:1717), csigpne/csigmne = p111*p1 + p055*(p2+p4) + p027*p3 (F90:1696, 1701) and csig12ne = p222*s1 + p111*(s2+s4) + p055*s3 (F90:1706). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the NE cell (i,j) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dy_cell_k6`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dy_cell_k6.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dy_cell_k6.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-(str_pairsum_lin&lt;a1=m1, a2=m4, a3=m2, a4=m3, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=s3, a2=s4, a3=s1, a4=s2, c=0.5, metric=dyT&gt;) + str_pairsum_lin&lt;a1=p1, a2=p4, a3=p2, a4=p3, c=0.25, metric=dxT&gt; + (-dyhx * (str_csig&lt;a1=p4, a2=p1, a3=p2, a4=p3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m4, a2=m1, a3=m2, a4=m3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;)) + dxhy * str_csig&lt;a1=s4, a2=s1, a3=s2, a4=s3, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,6), the SE-corner contribution to the y divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1756-1757): str6 = strp_tmp - strm_tmp - str12sn - dyhx*(csigpse + csigmse) + dxhy*csig12se, with strp_tmp/strm_tmp = p25*dxT*(p333*ssigpe + p166*ssigpw) (F90:1748-1749), str12sn = p5*dyT*(p333*ssig12s + p166*ssig12n) (F90:1718), csigpse/csigmse = p111*p4 + p055*(p1+p3) + p027*p2 (F90:1699, 1704) and csig12se = p222*s4 + p111*(s1+s3) + p055*s2 (F90:1712). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the SE cell (i,j+1) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dy_cell_k7`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dy_cell_k7.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dy_cell_k7.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-str_pairsum_lin&lt;a1=p2, a2=p3, a3=p1, a4=p4, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=m2, a2=m3, a3=m1, a4=m4, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=s1, a2=s2, a3=s3, a4=s4, c=0.5, metric=dyT&gt; + (-dyhx * (str_csig&lt;a1=p2, a2=p1, a3=p4, a4=p3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m2, a2=m1, a3=m4, a4=m3, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;)) + dxhy * str_csig&lt;a1=s2, a2=s1, a3=s4, a4=s3, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,7), the NW-corner contribution to the y divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1760-1761): str7 = -strp_tmp + strm_tmp + str12ns - dyhx*(csigpnw + csigmnw) + dxhy*csig12nw, with strp_tmp/strm_tmp = p25*dxT*(p333*ssigpw + p166*ssigpe) (F90:1754-1755), str12ns = p5*dyT*(p333*ssig12n + p166*ssig12s) (F90:1717), csigpnw/csigmnw = p111*p2 + p055*(p1+p3) + p027*p4 (F90:1697, 1702) and csig12nw = p222*s2 + p111*(s1+s3) + p055*s4 (F90:1708). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the NW cell (i+1,j) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_dy_cell_k8`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_dy_cell_k8.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_dy_cell_k8.esm)

Parameters: `p1`, `p2`, `p3`, `p4`, `m1`, `m2`, `m3`, `m4`, `s1`, `s2`, `s3`, `s4`, `dxT`, `dyT`, `dxhy`, `dyhx`

<div class="math-block">-str_pairsum_lin&lt;a1=m2, a2=m3, a3=m1, a4=m4, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=p2, a2=p3, a3=p1, a4=p4, c=0.25, metric=dxT&gt; + str_pairsum_lin&lt;a1=s3, a2=s4, a3=s1, a4=s2, c=0.5, metric=dyT&gt; + (-dyhx * (str_csig&lt;a1=p3, a2=p2, a3=p1, a4=p4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt; + str_csig&lt;a1=m3, a2=m2, a3=m1, a4=m4, c1=1 / 9, c2=1 / 18, c3=1 / 36&gt;)) + dxhy * str_csig&lt;a1=s3, a2=s2, a3=s1, a4=s4, c1=2 / 9, c2=1 / 9, c3=1 / 18&gt;</div>

Scalar atom of the CICE stress-combination algebra: str(i,j,8), the SW-corner contribution to the y divergence of the internal stress at the U point U(i,j) (stress, ice_dyn_evp.F90:1764-1765): str8 = strp_tmp - strm_tmp + str12sn - dyhx*(csigpsw + csigmsw) + dxhy*csig12sw, with strp_tmp/strm_tmp = p25*dxT*(p333*ssigpw + p166*ssigpe) (F90:1754-1755), str12sn = p5*dyT*(p333*ssig12s + p166*ssig12n) (F90:1718), csigpsw/csigmsw = p111*p3 + p055*(p2+p4) + p027*p1 (F90:1698, 1703) and csig12sw = p222*s3 + p111*(s2+s4) + p055*s1 (F90:1710). Params and structure as on str_dx_cell_k1; the rule binds the p/m/s values of the SW cell (i+1,j+1) around the U point.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `str_pairsum_lin`

<p class="rule-meta"><code>op:stress</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/str_pairsum_lin.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/str_pairsum_lin.esm)

Parameters: `a1`, `a2`, `a3`, `a4`, `c`, `metric`

<div class="math-block">c * metric * (1 / 3 * (a1 + a2) + 1 / 6 * (a3 + a4))</div>

Scalar atom of the CICE stress-combination algebra (ice_dyn_evp.F90:1687-1696, the str12* edge transports, and F90:1697-1700, the strp_tmp/strm_tmp pair sums): lin = c*metric*(p333*(a1 + a2) + p166*(a3 + a4)). The str12 terms bind a1/a2 = the two corners of one side pair and a3/a4 = the other side pair with c = p5 and metric = dxT (str12ew, F90:1687-1688) or dyT (str12ns, F90:1691-1692); the strp_tmp/strm_tmp terms bind a1/a2 = the north pair and a3/a4 = the south pair with c = p25 and metric = dyT (F90:1697-1698) or dxT (F90:1701-1702). Pure combination algebra: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the stress() corner combinations.

### `strain_divu_corner1`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_divu_corner1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_divu_corner1.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">cyp * uvel[oi, oj] - dyT * uvel[wim1, oj] + (cxp * vvel[oi, oj] - dxT * vvel[oi, wjm1])</div>

Scalar atom of the B-grid strain-rate divergence at corner 1 (the NE corner, the U point U(i,j)) of T cell (i,j): divu_1 = cyp*u(i,j) - dyT*u(i-1,j) + cxp*v(i,j) - dxT*v(i,j-1) (strain_rates, ice_dyn_shared.F90:2148-2151; the area-weighted form of divu = D(u,x) + D(v,y) on the B grid). The four U-point coordinates enter through the index params: oi/oj are the cell's own corner indices (oi = i, oj = j) and wim1/wjm1 are the west/south neighbor indices (i-1, j-1); the strain_rates_divu rule binds wim1 := 1 / wjm1 := 1 on the i = 1 / j = 1 halo strips, reproducing the Fortran's halo-copy of the velocity at the domain edge. The metric coefficients cxp/cyp/cxm/cym and cell dimensions dxT/dyT are scalars of the cell (ice_dyn_shared.F90:410-413, 433-437), bound by the rule to the cell's own values. Pure combination algebra: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), section on the EVP dynamics strain_rates.

### `strain_divu_corner2`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_divu_corner2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_divu_corner2.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">cym * uvel[wim1, oj] + dyT * uvel[oi, oj] + (cxp * vvel[wim1, oj] - dxT * vvel[wim1, wjm1])</div>

Scalar atom of the B-grid strain-rate divergence at corner 2 (the NW corner, U(i-1,j)) of T cell (i,j): divu_2 = cym*u(i-1,j) + dyT*u(i,j) + cxp*v(i-1,j) - dxT*v(i-1,j-1) (strain_rates, ice_dyn_shared.F90:2152-2155). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_divu rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_divu_corner3`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_divu_corner3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_divu_corner3.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">cym * uvel[wim1, wjm1] + dyT * uvel[oi, wjm1] + cxm * vvel[wim1, wjm1] + dxT * vvel[wim1, oj]</div>

Scalar atom of the B-grid strain-rate divergence at corner 3 (the SW corner, U(i-1,j-1)) of T cell (i,j): divu_3 = cym*u(i-1,j-1) + dyT*u(i,j-1) + cxm*v(i-1,j-1) + dxT*v(i-1,j) (strain_rates, ice_dyn_shared.F90:2156-2159). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_divu rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_divu_corner4`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_divu_corner4.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_divu_corner4.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">cyp * uvel[oi, wjm1] - dyT * uvel[wim1, wjm1] + cxm * vvel[oi, wjm1] + dxT * vvel[oi, oj]</div>

Scalar atom of the B-grid strain-rate divergence at corner 4 (the SE corner, U(i,j-1)) of T cell (i,j): divu_4 = cyp*u(i,j-1) - dyT*u(i-1,j-1) + cxm*v(i,j-1) + dxT*v(i,j) (strain_rates, ice_dyn_shared.F90:2160-2163). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_divu rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_shear_corner1`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_shear_corner1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_shear_corner1.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-(cym * vvel[oi, oj] + dyT * vvel[wim1, oj] + cxm * uvel[oi, oj] + dxT * uvel[oi, wjm1])</div>

Scalar atom of the B-grid shearing strain rate at corner 1 (NE, U(i,j)) of T cell (i,j): shear_1 = -cym*v(i,j) - dyT*v(i-1,j) - cxm*u(i,j) - dxT*u(i,j-1) (strain_rates, ice_dyn_shared.F90:2168-2169; the area-weighted form of shear = D(u,y) + D(v,x) = 2*e12 on the B grid). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_shear rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_shear_corner2`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_shear_corner2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_shear_corner2.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-cyp * vvel[wim1, oj] + dyT * vvel[oi, oj] - (cxm * uvel[wim1, oj] + dxT * uvel[wim1, wjm1])</div>

Scalar atom of the B-grid shearing strain rate at corner 2 (NW, U(i-1,j)) of T cell (i,j): shear_2 = -cyp*v(i-1,j) + dyT*v(i,j) - cxm*u(i-1,j) - dxT*u(i-1,j-1) (strain_rates, ice_dyn_shared.F90:2170-2171). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_shear rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_shear_corner3`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_shear_corner3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_shear_corner3.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-cyp * vvel[wim1, wjm1] + dyT * vvel[oi, wjm1] + -cxp * uvel[wim1, wjm1] + dxT * uvel[wim1, oj]</div>

Scalar atom of the B-grid shearing strain rate at corner 3 (SW, U(i-1,j-1)) of T cell (i,j): shear_3 = -cyp*v(i-1,j-1) + dyT*v(i,j-1) - cxp*u(i-1,j-1) + dxT*u(i-1,j) (strain_rates, ice_dyn_shared.F90:2172-2173). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_shear rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_shear_corner4`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_shear_corner4.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_shear_corner4.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-(cym * vvel[oi, wjm1] + dyT * vvel[wim1, wjm1]) + -cxp * uvel[oi, wjm1] + dxT * uvel[oi, oj]</div>

Scalar atom of the B-grid shearing strain rate at corner 4 (SE, U(i,j-1)) of T cell (i,j): shear_4 = -cym*v(i,j-1) - dyT*v(i-1,j-1) - cxp*u(i,j-1) + dxT*u(i,j) (strain_rates, ice_dyn_shared.F90:2174-2175). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_shear rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_tension_corner1`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_tension_corner1.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_tension_corner1.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-(cym * uvel[oi, oj] + dyT * uvel[wim1, oj]) + cxm * vvel[oi, oj] + dxT * vvel[oi, wjm1]</div>

Scalar atom of the B-grid tension strain rate at corner 1 (NE, U(i,j)) of T cell (i,j): tension_1 = -cym*u(i,j) - dyT*u(i-1,j) + cxm*v(i,j) + dxT*v(i,j-1) (strain_rates, ice_dyn_shared.F90:2158-2159; the area-weighted form of tension = D(u,x) - D(v,y) on the B grid). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_tension rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_tension_corner2`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_tension_corner2.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_tension_corner2.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-cyp * uvel[wim1, oj] + dyT * uvel[oi, oj] + cxm * vvel[wim1, oj] + dxT * vvel[wim1, wjm1]</div>

Scalar atom of the B-grid tension strain rate at corner 2 (NW, U(i-1,j)) of T cell (i,j): tension_2 = -cyp*u(i-1,j) + dyT*u(i,j) + cxm*v(i-1,j) + dxT*v(i-1,j-1) (strain_rates, ice_dyn_shared.F90:2160-2161). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_tension rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_tension_corner3`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_tension_corner3.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_tension_corner3.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-cyp * uvel[wim1, wjm1] + dyT * uvel[oi, wjm1] + (cxp * vvel[wim1, wjm1] - dxT * vvel[wim1, oj])</div>

Scalar atom of the B-grid tension strain rate at corner 3 (SW, U(i-1,j-1)) of T cell (i,j): tension_3 = -cyp*u(i-1,j-1) + dyT*u(i,j-1) + cxp*v(i-1,j-1) - dxT*v(i-1,j) (strain_rates, ice_dyn_shared.F90:2162-2163). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_tension rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `strain_tension_corner4`

<p class="rule-meta"><code>op:strain-rate</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/strain_tension_corner4.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/strain_tension_corner4.esm)

Parameters: `uvel`, `vvel`, `oi`, `oj`, `wim1`, `wjm1`, `cxp`, `cyp`, `cxm`, `cym`, `dxT`, `dyT`

<div class="math-block">-(cym * uvel[oi, wjm1] + dyT * uvel[wim1, wjm1]) + (cxp * vvel[oi, wjm1] - dxT * vvel[oi, oj])</div>

Scalar atom of the B-grid tension strain rate at corner 4 (SE, U(i,j-1)) of T cell (i,j): tension_4 = -cym*u(i,j-1) - dyT*u(i-1,j-1) + cxp*v(i,j-1) - dxT*v(i,j) (strain_rates, ice_dyn_shared.F90:2164-2165). Index params as on strain_divu_corner1 (oi/oj = i/j, wim1/wjm1 = the west/south neighbors, bound to 1 on the i = 1 / j = 1 halo strips by the strain_rates_tension rule); metric coefficients and cell dimensions are scalars of the cell.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851).

### `transport_flux_e`

<p class="rule-meta"><code>op:advection</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/transport_flux_e.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/transport_flux_e.esm)

Parameters: `cL`, `cR`, `uS`, `uN`, `HTE`

<div class="math-block">HTE * (uS + uN) / 2 * ((cL + cR) / 2)</div>

Scalar atom of the B-grid east-edge advective mass flux of the conservative transport: F_e = HTE * ((uS + uN)/2) * ((cL + cR)/2), the centered corner-average reconstruction of the CICE transport edge algebra. The east edge of T cell (i,j) runs between the two corner U points U(i,j-1) (south endpoint, velocity uS) and U(i,j) (north endpoint, velocity uN); the edge velocity is the average of the two corner velocities times the edge length HTE(i,j), exactly the Fortran's edge-area reconstruction edgearea_e(i,j) = (uvel(i,j) + uvel(i,j-1))*p5*HTE(i,j)*dt (ice_transport_remap.F90:658-661, the linear branch), and cL/cR are the tracer values in the two cells sharing the edge (cell (i,j) and cell (i+1,j)), averaged to the face. The rule transport_flux_divergence binds the six values from the surrounding cells and U points and differences the four edge fluxes around the cell; the incremental-remap scheme of the traced Fortran run computes its mflxe by departure-triangle integrals over the limited linear reconstruction (transport_integrals, ice_transport_remap.F90:3188-3472) rather than by this face average, so the rule implements the centered flux-divergence discretization of D(c*u,x) + D(c*v,y) (the mass update mm -= (mflxe(i,j) - mflxe(i-1,j) + mflxn(i,j) - mflxn(i,j-1))*tarear, F90:3586-3590), documented in the stage-2 transport pins. Pure combination algebra: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the horizontal transport edge reconstruction.

### `transport_flux_n`

<p class="rule-meta"><code>op:advection</code> <code>interior-only (match-less)</code></p>

Source: [`grids/cice_bgrid/stencils/transport_flux_n.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/stencils/transport_flux_n.esm)

Parameters: `cS`, `cN`, `vW`, `vE`, `HTN`

<div class="math-block">HTN * (vW + vE) / 2 * ((cS + cN) / 2)</div>

Scalar atom of the B-grid north-edge advective mass flux of the conservative transport: F_n = HTN * ((vW + vE)/2) * ((cS + cN)/2), the centered corner-average reconstruction of the CICE transport edge algebra. The north edge of T cell (i,j) runs between the two corner U points U(i-1,j) (west endpoint, velocity vW) and U(i,j) (east endpoint, velocity vE); the edge velocity is the average of the two corner velocities times the edge length HTN(i,j), exactly the Fortran's edge-area reconstruction edgearea_n(i,j) = (vvel(i,j) + vvel(i-1,j))*p5*HTN(i,j)*dt (ice_transport_remap.F90:662-666), and cS/cN are the tracer values in the two cells sharing the edge (cell (i,j) and cell (i,j+1)), averaged to the face. The rule transport_flux_divergence binds the six values from the surrounding cells and U points and differences the four edge fluxes around the cell; the scheme caveat of transport_flux_e (incremental remap vs centered flux divergence, mass update at F90:3586-3590) applies identically. Pure combination algebra: no faq, no index sets of its own.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the horizontal transport edge reconstruction.

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `strain_rates_divu`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j,corner</code></p>

Source: [`grids/cice_bgrid/rules/strain_rates_divu.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/strain_rates_divu.esm)

**Rewrites** `D(uvel)/Dx + D(vvel)/Dy` — match pattern:

```json
{"op": "+", "args": [{"op": "D", "args": ["uvel"], "wrt": "x"}, {"op": "D", "args": ["vvel"], "wrt": "y"}]}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N_I, 2:N_J, 1:1] = sum[i, j, k] (strain_divu_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 1:1}, [2:N_I, 2:N_J, 2:2] = sum[i, j, k] (strain_divu_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 2:2}, [2:N_I, 2:N_J, 3:3] = sum[i, j, k] (strain_divu_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 3:3}, [2:N_I, 2:N_J, 4:4] = sum[i, j, k] (strain_divu_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 4:4}, [1:1, 2:N_J, 1:1] = sum[i, j, k] (strain_divu_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 1:1}, [1:1, 2:N_J, 2:2] = sum[i, j, k] (strain_divu_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 2:2}, [1:1, 2:N_J, 3:3] = sum[i, j, k] (strain_divu_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 3:3}, [1:1, 2:N_J, 4:4] = sum[i, j, k] (strain_divu_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 4:4}, [2:N_I, 1:1, 1:1] = sum[i, j, k] (strain_divu_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 1:1}, [2:N_I, 1:1, 2:2] = sum[i, j, k] (strain_divu_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 2:2}, [2:N_I, 1:1, 3:3] = sum[i, j, k] (strain_divu_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 3:3}, [2:N_I, 1:1, 4:4] = sum[i, j, k] (strain_divu_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 4:4}, [1:1, 1:1, 1:1] = sum[i, j, k] (strain_divu_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 1:1}, [1:1, 1:1, 2:2] = sum[i, j, k] (strain_divu_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 2:2}, [1:1, 1:1, 3:3] = sum[i, j, k] (strain_divu_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 3:3}, [1:1, 1:1, 4:4] = sum[i, j, k] (strain_divu_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 4:4})</div>

B-grid lowering rule of the strain-rate divergence rewrite target divu = D(uvel, x) + D(vvel, y) (strain_rates, ice_dyn_shared.F90:2094-2105; the corner formulas 2148-2176): the rule matches the two-term sum by GROUND PATTERN (the arguments are literal variable names uvel/vvel, not params, so no where-shape constraint is possible or needed, esm-spec 9.6.1/9.6.8) and emits the [i,j,corner] makearray of the four corner combinations, one region per corner slice (corner 1 = NE = U(i,j), 2 = NW = U(i-1,j), 3 = SW = U(i-1,j-1), 4 = SE = U(i,j-1)) tiled over four i/j position classes (interior [2..N_I]x[2..N_J], the i = 1 strip, the j = 1 strip, and the (1,1) corner) with halo-copy closures: on the i = 1 strip the west-neighbor index wim1 is bound to 1 and on the j = 1 strip the south-neighbor index wjm1 to 1 (the velocity of the ghost cell copies the edge value, as the Fortran's ice_HaloUpdate of uvel/vvel does before strain_rates). Each region value is a faq over its box that applies the scalar corner atom (strain_divu_corner1..4, ice_dyn_shared.F90:2148-2176) with the index params oi = i, oj = j, wim1, wjm1 and the cell's own metric coefficients cxp/cyp/cxm/cym and dimensions dxT/dyT (F90:410-413, 433-437), e.g. corner 1: divu = cyp*u(i,j) - dyT*u(i-1,j) + cxp*v(i,j) - dxT*v(i,j-1). Requires N_I >= 2 and N_J >= 2.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the EVP dynamics strain_rates.

**Conformance fixtures**

- ast: [`cice_bgrid_strain_divu`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_strain_divu) — status `active`

### `strain_rates_shear`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j,corner</code></p>

Source: [`grids/cice_bgrid/rules/strain_rates_shear.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/strain_rates_shear.esm)

**Rewrites** `D(uvel)/Dy + D(vvel)/Dx` — match pattern:

```json
{"op": "+", "args": [{"op": "D", "args": ["uvel"], "wrt": "y"}, {"op": "D", "args": ["vvel"], "wrt": "x"}]}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N_I, 2:N_J, 1:1] = sum[i, j, k] (strain_shear_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 1:1}, [2:N_I, 2:N_J, 2:2] = sum[i, j, k] (strain_shear_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 2:2}, [2:N_I, 2:N_J, 3:3] = sum[i, j, k] (strain_shear_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 3:3}, [2:N_I, 2:N_J, 4:4] = sum[i, j, k] (strain_shear_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 4:4}, [1:1, 2:N_J, 1:1] = sum[i, j, k] (strain_shear_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 1:1}, [1:1, 2:N_J, 2:2] = sum[i, j, k] (strain_shear_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 2:2}, [1:1, 2:N_J, 3:3] = sum[i, j, k] (strain_shear_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 3:3}, [1:1, 2:N_J, 4:4] = sum[i, j, k] (strain_shear_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 4:4}, [2:N_I, 1:1, 1:1] = sum[i, j, k] (strain_shear_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 1:1}, [2:N_I, 1:1, 2:2] = sum[i, j, k] (strain_shear_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 2:2}, [2:N_I, 1:1, 3:3] = sum[i, j, k] (strain_shear_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 3:3}, [2:N_I, 1:1, 4:4] = sum[i, j, k] (strain_shear_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 4:4}, [1:1, 1:1, 1:1] = sum[i, j, k] (strain_shear_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 1:1}, [1:1, 1:1, 2:2] = sum[i, j, k] (strain_shear_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 2:2}, [1:1, 1:1, 3:3] = sum[i, j, k] (strain_shear_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 3:3}, [1:1, 1:1, 4:4] = sum[i, j, k] (strain_shear_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 4:4})</div>

B-grid lowering rule of the strain-rate shearing rewrite target shear = D(uvel, y) + D(vvel, x) (strain_rates, ice_dyn_shared.F90:2094-2105; the corner formulas 2168-2175): the rule matches the two-term sum by GROUND PATTERN (the arguments are literal variable names uvel/vvel, not params, so no where-shape constraint is possible or needed, esm-spec 9.6.1/9.6.8) and emits the [i,j,corner] makearray of the four corner combinations, one region per corner slice (corner 1 = NE = U(i,j), 2 = NW = U(i-1,j), 3 = SW = U(i-1,j-1), 4 = SE = U(i,j-1)) tiled over four i/j position classes (interior [2..N_I]x[2..N_J], the i = 1 strip, the j = 1 strip, and the (1,1) corner) with halo-copy closures: on the i = 1 strip the west-neighbor index wim1 is bound to 1 and on the j = 1 strip the south-neighbor index wjm1 to 1 (the velocity of the ghost cell copies the edge value, as the Fortran's ice_HaloUpdate of uvel/vvel does before strain_rates). Each region value is a faq over its box that applies the scalar corner atom (strain_shear_corner1..4, ice_dyn_shared.F90:2168-2175) with the index params oi = i, oj = j, wim1, wjm1 and the cell's own metric coefficients cxp/cyp/cxm/cym and dimensions dxT/dyT (F90:410-413, 433-437), e.g. corner 1: shear = -cym*v(i,j) - dyT*v(i-1,j) - cxm*u(i,j) - dxT*u(i,j-1). Requires N_I >= 2 and N_J >= 2.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the EVP dynamics strain_rates.

**Conformance fixtures**

- ast: [`cice_bgrid_strain_shear`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_strain_shear) — status `active`

### `strain_rates_tension`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j,corner</code></p>

Source: [`grids/cice_bgrid/rules/strain_rates_tension.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/strain_rates_tension.esm)

**Rewrites** `D(uvel)/Dx - D(vvel)/Dy` — match pattern:

```json
{"op": "-", "args": [{"op": "D", "args": ["uvel"], "wrt": "x"}, {"op": "D", "args": ["vvel"], "wrt": "y"}]}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N_I, 2:N_J, 1:1] = sum[i, j, k] (strain_tension_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 1:1}, [2:N_I, 2:N_J, 2:2] = sum[i, j, k] (strain_tension_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 2:2}, [2:N_I, 2:N_J, 3:3] = sum[i, j, k] (strain_tension_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 3:3}, [2:N_I, 2:N_J, 4:4] = sum[i, j, k] (strain_tension_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=j - 1&gt;) where {i in 2:N_I, j in 2:N_J, k in 4:4}, [1:1, 2:N_J, 1:1] = sum[i, j, k] (strain_tension_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 1:1}, [1:1, 2:N_J, 2:2] = sum[i, j, k] (strain_tension_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 2:2}, [1:1, 2:N_J, 3:3] = sum[i, j, k] (strain_tension_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 3:3}, [1:1, 2:N_J, 4:4] = sum[i, j, k] (strain_tension_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=j - 1&gt;) where {i in 1:1, j in 2:N_J, k in 4:4}, [2:N_I, 1:1, 1:1] = sum[i, j, k] (strain_tension_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 1:1}, [2:N_I, 1:1, 2:2] = sum[i, j, k] (strain_tension_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 2:2}, [2:N_I, 1:1, 3:3] = sum[i, j, k] (strain_tension_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 3:3}, [2:N_I, 1:1, 4:4] = sum[i, j, k] (strain_tension_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=i - 1, wjm1=1&gt;) where {i in 2:N_I, j in 1:1, k in 4:4}, [1:1, 1:1, 1:1] = sum[i, j, k] (strain_tension_corner1&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 1:1}, [1:1, 1:1, 2:2] = sum[i, j, k] (strain_tension_corner2&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 2:2}, [1:1, 1:1, 3:3] = sum[i, j, k] (strain_tension_corner3&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 3:3}, [1:1, 1:1, 4:4] = sum[i, j, k] (strain_tension_corner4&lt;cxm=cxm[i, j], cxp=cxp[i, j], cym=cym[i, j], cyp=cyp[i, j], dxT=dxT[i, j], dyT=dyT[i, j], oi=i, oj=j, uvel=uvel, vvel=vvel, wim1=1, wjm1=1&gt;) where {i in 1:1, j in 1:1, k in 4:4})</div>

B-grid lowering rule of the strain-rate tension rewrite target tension = D(uvel, x) - D(vvel, y) (strain_rates, ice_dyn_shared.F90:2094-2105; the corner formulas 2158-2165): the rule matches the two-term sum by GROUND PATTERN (the arguments are literal variable names uvel/vvel, not params, so no where-shape constraint is possible or needed, esm-spec 9.6.1/9.6.8) and emits the [i,j,corner] makearray of the four corner combinations, one region per corner slice (corner 1 = NE = U(i,j), 2 = NW = U(i-1,j), 3 = SW = U(i-1,j-1), 4 = SE = U(i,j-1)) tiled over four i/j position classes (interior [2..N_I]x[2..N_J], the i = 1 strip, the j = 1 strip, and the (1,1) corner) with halo-copy closures: on the i = 1 strip the west-neighbor index wim1 is bound to 1 and on the j = 1 strip the south-neighbor index wjm1 to 1 (the velocity of the ghost cell copies the edge value, as the Fortran's ice_HaloUpdate of uvel/vvel does before strain_rates). Each region value is a faq over its box that applies the scalar corner atom (strain_tension_corner1..4, ice_dyn_shared.F90:2158-2165) with the index params oi = i, oj = j, wim1, wjm1 and the cell's own metric coefficients cxp/cyp/cxm/cym and dimensions dxT/dyT (F90:410-413, 433-437), e.g. corner 1: tension = -cym*u(i,j) - dyT*u(i-1,j) + cxm*v(i,j) + dxT*v(i,j-1). Requires N_I >= 2 and N_J >= 2.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the EVP dynamics strain_rates.

**Conformance fixtures**

- ast: [`cice_bgrid_strain_tension`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_strain_tension) — status `active`

### `stress_divergence_strintx`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j</code></p>

Source: [`grids/cice_bgrid/rules/stress_divergence_strintx.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/stress_divergence_strintx.esm)

**Rewrites** `D(sigma11)/Dx + D(sigma12)/Dy` — match pattern:

```json
{"op": "+", "args": [{"op": "D", "args": ["sigma11"], "wrt": "x"}, {"op": "D", "args": ["sigma12"], "wrt": "y"}]}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:N_I - 1, 1:N_J - 1] = sum[i, j] (uarear[i, j] * (str_dx_cell_k1&lt;dxT=dxT[i, j], dxhy=dxhy[i, j], dyT=dyT[i, j], dyhx=dyhx[i, j], m1=sigma11[i, j, 1] - sigma22[i, j, 1], m2=sigma11[i, j, 2] - sigma22[i, j, 2], m3=sigma11[i, j, 3] - sigma22[i, j, 3], m4=sigma11[i, j, 4] - sigma22[i, j, 4], p1=sigma11[i, j, 1] + sigma22[i, j, 1], p2=sigma11[i, j, 2] + sigma22[i, j, 2], p3=sigma11[i, j, 3] + sigma22[i, j, 3], p4=sigma11[i, j, 4] + sigma22[i, j, 4], s1=sigma12[i, j, 1], s2=sigma12[i, j, 2], s3=sigma12[i, j, 3], s4=sigma12[i, j, 4]&gt; + str_dx_cell_k2&lt;dxT=dxT[i + 1, j], dxhy=dxhy[i + 1, j], dyT=dyT[i + 1, j], dyhx=dyhx[i + 1, j], m1=sigma11[i + 1, j, 1] - sigma22[i + 1, j, 1], m2=sigma11[i + 1, j, 2] - sigma22[i + 1, j, 2], m3=sigma11[i + 1, j, 3] - sigma22[i + 1, j, 3], m4=sigma11[i + 1, j, 4] - sigma22[i + 1, j, 4], p1=sigma11[i + 1, j, 1] + sigma22[i + 1, j, 1], p2=sigma11[i + 1, j, 2] + sigma22[i + 1, j, 2], p3=sigma11[i + 1, j, 3] + sigma22[i + 1, j, 3], p4=sigma11[i + 1, j, 4] + sigma22[i + 1, j, 4], s1=sigma12[i + 1, j, 1], s2=sigma12[i + 1, j, 2], s3=sigma12[i + 1, j, 3], s4=sigma12[i + 1, j, 4]&gt; + str_dx_cell_k3&lt;dxT=dxT[i, j + 1], dxhy=dxhy[i, j + 1], dyT=dyT[i, j + 1], dyhx=dyhx[i, j + 1], m1=sigma11[i, j + 1, 1] - sigma22[i, j + 1, 1], m2=sigma11[i, j + 1, 2] - sigma22[i, j + 1, 2], m3=sigma11[i, j + 1, 3] - sigma22[i, j + 1, 3], m4=sigma11[i, j + 1, 4] - sigma22[i, j + 1, 4], p1=sigma11[i, j + 1, 1] + sigma22[i, j + 1, 1], p2=sigma11[i, j + 1, 2] + sigma22[i, j + 1, 2], p3=sigma11[i, j + 1, 3] + sigma22[i, j + 1, 3], p4=sigma11[i, j + 1, 4] + sigma22[i, j + 1, 4], s1=sigma12[i, j + 1, 1], s2=sigma12[i, j + 1, 2], s3=sigma12[i, j + 1, 3], s4=sigma12[i, j + 1, 4]&gt; + str_dx_cell_k4&lt;dxT=dxT[i + 1, j + 1], dxhy=dxhy[i + 1, j + 1], dyT=dyT[i + 1, j + 1], dyhx=dyhx[i + 1, j + 1], m1=sigma11[i + 1, j + 1, 1] - sigma22[i + 1, j + 1, 1], m2=sigma11[i + 1, j + 1, 2] - sigma22[i + 1, j + 1, 2], m3=sigma11[i + 1, j + 1, 3] - sigma22[i + 1, j + 1, 3], m4=sigma11[i + 1, j + 1, 4] - sigma22[i + 1, j + 1, 4], p1=sigma11[i + 1, j + 1, 1] + sigma22[i + 1, j + 1, 1], p2=sigma11[i + 1, j + 1, 2] + sigma22[i + 1, j + 1, 2], p3=sigma11[i + 1, j + 1, 3] + sigma22[i + 1, j + 1, 3], p4=sigma11[i + 1, j + 1, 4] + sigma22[i + 1, j + 1, 4], s1=sigma12[i + 1, j + 1, 1], s2=sigma12[i + 1, j + 1, 2], s3=sigma12[i + 1, j + 1, 3], s4=sigma12[i + 1, j + 1, 4]&gt;)) where {i in 1:N_I - 1, j in 1:N_J - 1}, [N_I:N_I, 1:N_J - 1] = sum[i, j] (uarear[i, j] * (str_dx_cell_k1&lt;dxT=dxT[N_I, j], dxhy=dxhy[N_I, j], dyT=dyT[N_I, j], dyhx=dyhx[N_I, j], m1=sigma11[N_I, j, 1] - sigma22[N_I, j, 1], m2=sigma11[N_I, j, 2] - sigma22[N_I, j, 2], m3=sigma11[N_I, j, 3] - sigma22[N_I, j, 3], m4=sigma11[N_I, j, 4] - sigma22[N_I, j, 4], p1=sigma11[N_I, j, 1] + sigma22[N_I, j, 1], p2=sigma11[N_I, j, 2] + sigma22[N_I, j, 2], p3=sigma11[N_I, j, 3] + sigma22[N_I, j, 3], p4=sigma11[N_I, j, 4] + sigma22[N_I, j, 4], s1=sigma12[N_I, j, 1], s2=sigma12[N_I, j, 2], s3=sigma12[N_I, j, 3], s4=sigma12[N_I, j, 4]&gt; + str_dx_cell_k2&lt;dxT=dxT[N_I, j], dxhy=dxhy[N_I, j], dyT=dyT[N_I, j], dyhx=dyhx[N_I, j], m1=sigma11[N_I, j, 1] - sigma22[N_I, j, 1], m2=sigma11[N_I, j, 2] - sigma22[N_I, j, 2], m3=sigma11[N_I, j, 3] - sigma22[N_I, j, 3], m4=sigma11[N_I, j, 4] - sigma22[N_I, j, 4], p1=sigma11[N_I, j, 1] + sigma22[N_I, j, 1], p2=sigma11[N_I, j, 2] + sigma22[N_I, j, 2], p3=sigma11[N_I, j, 3] + sigma22[N_I, j, 3], p4=sigma11[N_I, j, 4] + sigma22[N_I, j, 4], s1=sigma12[N_I, j, 1], s2=sigma12[N_I, j, 2], s3=sigma12[N_I, j, 3], s4=sigma12[N_I, j, 4]&gt; + str_dx_cell_k3&lt;dxT=dxT[N_I, j + 1], dxhy=dxhy[N_I, j + 1], dyT=dyT[N_I, j + 1], dyhx=dyhx[N_I, j + 1], m1=sigma11[N_I, j + 1, 1] - sigma22[N_I, j + 1, 1], m2=sigma11[N_I, j + 1, 2] - sigma22[N_I, j + 1, 2], m3=sigma11[N_I, j + 1, 3] - sigma22[N_I, j + 1, 3], m4=sigma11[N_I, j + 1, 4] - sigma22[N_I, j + 1, 4], p1=sigma11[N_I, j + 1, 1] + sigma22[N_I, j + 1, 1], p2=sigma11[N_I, j + 1, 2] + sigma22[N_I, j + 1, 2], p3=sigma11[N_I, j + 1, 3] + sigma22[N_I, j + 1, 3], p4=sigma11[N_I, j + 1, 4] + sigma22[N_I, j + 1, 4], s1=sigma12[N_I, j + 1, 1], s2=sigma12[N_I, j + 1, 2], s3=sigma12[N_I, j + 1, 3], s4=sigma12[N_I, j + 1, 4]&gt; + str_dx_cell_k4&lt;dxT=dxT[N_I, j + 1], dxhy=dxhy[N_I, j + 1], dyT=dyT[N_I, j + 1], dyhx=dyhx[N_I, j + 1], m1=sigma11[N_I, j + 1, 1] - sigma22[N_I, j + 1, 1], m2=sigma11[N_I, j + 1, 2] - sigma22[N_I, j + 1, 2], m3=sigma11[N_I, j + 1, 3] - sigma22[N_I, j + 1, 3], m4=sigma11[N_I, j + 1, 4] - sigma22[N_I, j + 1, 4], p1=sigma11[N_I, j + 1, 1] + sigma22[N_I, j + 1, 1], p2=sigma11[N_I, j + 1, 2] + sigma22[N_I, j + 1, 2], p3=sigma11[N_I, j + 1, 3] + sigma22[N_I, j + 1, 3], p4=sigma11[N_I, j + 1, 4] + sigma22[N_I, j + 1, 4], s1=sigma12[N_I, j + 1, 1], s2=sigma12[N_I, j + 1, 2], s3=sigma12[N_I, j + 1, 3], s4=sigma12[N_I, j + 1, 4]&gt;)) where {i in N_I:N_I, j in 1:N_J - 1}, [1:N_I - 1, N_J:N_J] = sum[i, j] (uarear[i, j] * (str_dx_cell_k1&lt;dxT=dxT[i, N_J], dxhy=dxhy[i, N_J], dyT=dyT[i, N_J], dyhx=dyhx[i, N_J], m1=sigma11[i, N_J, 1] - sigma22[i, N_J, 1], m2=sigma11[i, N_J, 2] - sigma22[i, N_J, 2], m3=sigma11[i, N_J, 3] - sigma22[i, N_J, 3], m4=sigma11[i, N_J, 4] - sigma22[i, N_J, 4], p1=sigma11[i, N_J, 1] + sigma22[i, N_J, 1], p2=sigma11[i, N_J, 2] + sigma22[i, N_J, 2], p3=sigma11[i, N_J, 3] + sigma22[i, N_J, 3], p4=sigma11[i, N_J, 4] + sigma22[i, N_J, 4], s1=sigma12[i, N_J, 1], s2=sigma12[i, N_J, 2], s3=sigma12[i, N_J, 3], s4=sigma12[i, N_J, 4]&gt; + str_dx_cell_k2&lt;dxT=dxT[i + 1, N_J], dxhy=dxhy[i + 1, N_J], dyT=dyT[i + 1, N_J], dyhx=dyhx[i + 1, N_J], m1=sigma11[i + 1, N_J, 1] - sigma22[i + 1, N_J, 1], m2=sigma11[i + 1, N_J, 2] - sigma22[i + 1, N_J, 2], m3=sigma11[i + 1, N_J, 3] - sigma22[i + 1, N_J, 3], m4=sigma11[i + 1, N_J, 4] - sigma22[i + 1, N_J, 4], p1=sigma11[i + 1, N_J, 1] + sigma22[i + 1, N_J, 1], p2=sigma11[i + 1, N_J, 2] + sigma22[i + 1, N_J, 2], p3=sigma11[i + 1, N_J, 3] + sigma22[i + 1, N_J, 3], p4=sigma11[i + 1, N_J, 4] + sigma22[i + 1, N_J, 4], s1=sigma12[i + 1, N_J, 1], s2=sigma12[i + 1, N_J, 2], s3=sigma12[i + 1, N_J, 3], s4=sigma12[i + 1, N_J, 4]&gt; + str_dx_cell_k3&lt;dxT=dxT[i, N_J], dxhy=dxhy[i, N_J], dyT=dyT[i, N_J], dyhx=dyhx[i, N_J], m1=sigma11[i, N_J, 1] - sigma22[i, N_J, 1], m2=sigma11[i, N_J, 2] - sigma22[i, N_J, 2], m3=sigma11[i, N_J, 3] - sigma22[i, N_J, 3], m4=sigma11[i, N_J, 4] - sigma22[i, N_J, 4], p1=sigma11[i, N_J, 1] + sigma22[i, N_J, 1], p2=sigma11[i, N_J, 2] + sigma22[i, N_J, 2], p3=sigma11[i, N_J, 3] + sigma22[i, N_J, 3], p4=sigma11[i, N_J, 4] + sigma22[i, N_J, 4], s1=sigma12[i, N_J, 1], s2=sigma12[i, N_J, 2], s3=sigma12[i, N_J, 3], s4=sigma12[i, N_J, 4]&gt; + str_dx_cell_k4&lt;dxT=dxT[i + 1, N_J], dxhy=dxhy[i + 1, N_J], dyT=dyT[i + 1, N_J], dyhx=dyhx[i + 1, N_J], m1=sigma11[i + 1, N_J, 1] - sigma22[i + 1, N_J, 1], m2=sigma11[i + 1, N_J, 2] - sigma22[i + 1, N_J, 2], m3=sigma11[i + 1, N_J, 3] - sigma22[i + 1, N_J, 3], m4=sigma11[i + 1, N_J, 4] - sigma22[i + 1, N_J, 4], p1=sigma11[i + 1, N_J, 1] + sigma22[i + 1, N_J, 1], p2=sigma11[i + 1, N_J, 2] + sigma22[i + 1, N_J, 2], p3=sigma11[i + 1, N_J, 3] + sigma22[i + 1, N_J, 3], p4=sigma11[i + 1, N_J, 4] + sigma22[i + 1, N_J, 4], s1=sigma12[i + 1, N_J, 1], s2=sigma12[i + 1, N_J, 2], s3=sigma12[i + 1, N_J, 3], s4=sigma12[i + 1, N_J, 4]&gt;)) where {i in 1:N_I - 1, j in N_J:N_J}, [N_I:N_I, N_J:N_J] = sum[i, j] (uarear[i, j] * (str_dx_cell_k1&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dx_cell_k2&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dx_cell_k3&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dx_cell_k4&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt;)) where {i in N_I:N_I, j in N_J:N_J})</div>

B-grid lowering rule of the internal-stress x-divergence rewrite target strintx = D(sigma11, x) + D(sigma12, y), the momentum-equation divergence of the internal stress tensor: the Fortran assembles it at each U point as strintx(i,j) = uarear(i,j)*(str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4)) (stepu, ice_dyn_shared.F90:958-959) from the str corner combinations of the stress() routine (ice_dyn_evp.F90:1727-1743). The rule matches the two-term sum by GROUND PATTERN (sigma11/sigma12 literal, so no where-shape constraint is possible or needed, esm-spec 9.6.1/9.6.8) and emits the [i,j] makearray of the four surrounding-cell corner contributions tiled over four i/j position classes (interior [1..N_I-1]x[1..N_J-1], the i = N_I strip, the j = N_J strip, and the (N_I,N_J) corner) with halo-copy closures: on the i = N_I strip the cells (i+1,j) and (i+1,j+1) that supply the k2/k4 str atoms are replaced by (N_I,j)/(N_I,j+1), and on the j = N_J strip the cells (i,j+1) and (i+1,j+1) by (i,N_J)/(i+1,N_J) (the ghost-cell stress copies the edge value, as the Fortran's ice_HaloUpdate of the stress arrays does before stress()). Each region value is a faq over its box that applies the four scalar str atoms str_dx_cell_k1..k4, each bound to its own cell's corner stresses and metrics: the p_k = sigma11 + sigma22, m_k = sigma11 - sigma22, s_k = sigma12 corner combinations (the stressp/stressm/stress12 of the stress() routine, F90:1667-1686) and the cell's dxT/dyT/dxhy/dyhx (ice_dyn_evp.F90:410-416), the whole sum scaled by uarear(i,j), the U-cell inverse area (ice_grid.F90:703). Requires N_I >= 2 and N_J >= 2.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the EVP dynamics stress and stepu.

**Conformance fixtures**

- ast: [`cice_bgrid_strintx`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_strintx) — status `active`

### `stress_divergence_strinty`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j</code></p>

Source: [`grids/cice_bgrid/rules/stress_divergence_strinty.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/stress_divergence_strinty.esm)

**Rewrites** `D(sigma12)/Dx + D(sigma22)/Dy` — match pattern:

```json
{"op": "+", "args": [{"op": "D", "args": ["sigma12"], "wrt": "x"}, {"op": "D", "args": ["sigma22"], "wrt": "y"}]}
```

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:N_I - 1, 1:N_J - 1] = sum[i, j] (uarear[i, j] * (str_dy_cell_k5&lt;dxT=dxT[i, j], dxhy=dxhy[i, j], dyT=dyT[i, j], dyhx=dyhx[i, j], m1=sigma11[i, j, 1] - sigma22[i, j, 1], m2=sigma11[i, j, 2] - sigma22[i, j, 2], m3=sigma11[i, j, 3] - sigma22[i, j, 3], m4=sigma11[i, j, 4] - sigma22[i, j, 4], p1=sigma11[i, j, 1] + sigma22[i, j, 1], p2=sigma11[i, j, 2] + sigma22[i, j, 2], p3=sigma11[i, j, 3] + sigma22[i, j, 3], p4=sigma11[i, j, 4] + sigma22[i, j, 4], s1=sigma12[i, j, 1], s2=sigma12[i, j, 2], s3=sigma12[i, j, 3], s4=sigma12[i, j, 4]&gt; + str_dy_cell_k6&lt;dxT=dxT[i, j + 1], dxhy=dxhy[i, j + 1], dyT=dyT[i, j + 1], dyhx=dyhx[i, j + 1], m1=sigma11[i, j + 1, 1] - sigma22[i, j + 1, 1], m2=sigma11[i, j + 1, 2] - sigma22[i, j + 1, 2], m3=sigma11[i, j + 1, 3] - sigma22[i, j + 1, 3], m4=sigma11[i, j + 1, 4] - sigma22[i, j + 1, 4], p1=sigma11[i, j + 1, 1] + sigma22[i, j + 1, 1], p2=sigma11[i, j + 1, 2] + sigma22[i, j + 1, 2], p3=sigma11[i, j + 1, 3] + sigma22[i, j + 1, 3], p4=sigma11[i, j + 1, 4] + sigma22[i, j + 1, 4], s1=sigma12[i, j + 1, 1], s2=sigma12[i, j + 1, 2], s3=sigma12[i, j + 1, 3], s4=sigma12[i, j + 1, 4]&gt; + str_dy_cell_k7&lt;dxT=dxT[i + 1, j], dxhy=dxhy[i + 1, j], dyT=dyT[i + 1, j], dyhx=dyhx[i + 1, j], m1=sigma11[i + 1, j, 1] - sigma22[i + 1, j, 1], m2=sigma11[i + 1, j, 2] - sigma22[i + 1, j, 2], m3=sigma11[i + 1, j, 3] - sigma22[i + 1, j, 3], m4=sigma11[i + 1, j, 4] - sigma22[i + 1, j, 4], p1=sigma11[i + 1, j, 1] + sigma22[i + 1, j, 1], p2=sigma11[i + 1, j, 2] + sigma22[i + 1, j, 2], p3=sigma11[i + 1, j, 3] + sigma22[i + 1, j, 3], p4=sigma11[i + 1, j, 4] + sigma22[i + 1, j, 4], s1=sigma12[i + 1, j, 1], s2=sigma12[i + 1, j, 2], s3=sigma12[i + 1, j, 3], s4=sigma12[i + 1, j, 4]&gt; + str_dy_cell_k8&lt;dxT=dxT[i + 1, j + 1], dxhy=dxhy[i + 1, j + 1], dyT=dyT[i + 1, j + 1], dyhx=dyhx[i + 1, j + 1], m1=sigma11[i + 1, j + 1, 1] - sigma22[i + 1, j + 1, 1], m2=sigma11[i + 1, j + 1, 2] - sigma22[i + 1, j + 1, 2], m3=sigma11[i + 1, j + 1, 3] - sigma22[i + 1, j + 1, 3], m4=sigma11[i + 1, j + 1, 4] - sigma22[i + 1, j + 1, 4], p1=sigma11[i + 1, j + 1, 1] + sigma22[i + 1, j + 1, 1], p2=sigma11[i + 1, j + 1, 2] + sigma22[i + 1, j + 1, 2], p3=sigma11[i + 1, j + 1, 3] + sigma22[i + 1, j + 1, 3], p4=sigma11[i + 1, j + 1, 4] + sigma22[i + 1, j + 1, 4], s1=sigma12[i + 1, j + 1, 1], s2=sigma12[i + 1, j + 1, 2], s3=sigma12[i + 1, j + 1, 3], s4=sigma12[i + 1, j + 1, 4]&gt;)) where {i in 1:N_I - 1, j in 1:N_J - 1}, [N_I:N_I, 1:N_J - 1] = sum[i, j] (uarear[i, j] * (str_dy_cell_k5&lt;dxT=dxT[N_I, j], dxhy=dxhy[N_I, j], dyT=dyT[N_I, j], dyhx=dyhx[N_I, j], m1=sigma11[N_I, j, 1] - sigma22[N_I, j, 1], m2=sigma11[N_I, j, 2] - sigma22[N_I, j, 2], m3=sigma11[N_I, j, 3] - sigma22[N_I, j, 3], m4=sigma11[N_I, j, 4] - sigma22[N_I, j, 4], p1=sigma11[N_I, j, 1] + sigma22[N_I, j, 1], p2=sigma11[N_I, j, 2] + sigma22[N_I, j, 2], p3=sigma11[N_I, j, 3] + sigma22[N_I, j, 3], p4=sigma11[N_I, j, 4] + sigma22[N_I, j, 4], s1=sigma12[N_I, j, 1], s2=sigma12[N_I, j, 2], s3=sigma12[N_I, j, 3], s4=sigma12[N_I, j, 4]&gt; + str_dy_cell_k6&lt;dxT=dxT[N_I, j + 1], dxhy=dxhy[N_I, j + 1], dyT=dyT[N_I, j + 1], dyhx=dyhx[N_I, j + 1], m1=sigma11[N_I, j + 1, 1] - sigma22[N_I, j + 1, 1], m2=sigma11[N_I, j + 1, 2] - sigma22[N_I, j + 1, 2], m3=sigma11[N_I, j + 1, 3] - sigma22[N_I, j + 1, 3], m4=sigma11[N_I, j + 1, 4] - sigma22[N_I, j + 1, 4], p1=sigma11[N_I, j + 1, 1] + sigma22[N_I, j + 1, 1], p2=sigma11[N_I, j + 1, 2] + sigma22[N_I, j + 1, 2], p3=sigma11[N_I, j + 1, 3] + sigma22[N_I, j + 1, 3], p4=sigma11[N_I, j + 1, 4] + sigma22[N_I, j + 1, 4], s1=sigma12[N_I, j + 1, 1], s2=sigma12[N_I, j + 1, 2], s3=sigma12[N_I, j + 1, 3], s4=sigma12[N_I, j + 1, 4]&gt; + str_dy_cell_k7&lt;dxT=dxT[N_I, j], dxhy=dxhy[N_I, j], dyT=dyT[N_I, j], dyhx=dyhx[N_I, j], m1=sigma11[N_I, j, 1] - sigma22[N_I, j, 1], m2=sigma11[N_I, j, 2] - sigma22[N_I, j, 2], m3=sigma11[N_I, j, 3] - sigma22[N_I, j, 3], m4=sigma11[N_I, j, 4] - sigma22[N_I, j, 4], p1=sigma11[N_I, j, 1] + sigma22[N_I, j, 1], p2=sigma11[N_I, j, 2] + sigma22[N_I, j, 2], p3=sigma11[N_I, j, 3] + sigma22[N_I, j, 3], p4=sigma11[N_I, j, 4] + sigma22[N_I, j, 4], s1=sigma12[N_I, j, 1], s2=sigma12[N_I, j, 2], s3=sigma12[N_I, j, 3], s4=sigma12[N_I, j, 4]&gt; + str_dy_cell_k8&lt;dxT=dxT[N_I, j + 1], dxhy=dxhy[N_I, j + 1], dyT=dyT[N_I, j + 1], dyhx=dyhx[N_I, j + 1], m1=sigma11[N_I, j + 1, 1] - sigma22[N_I, j + 1, 1], m2=sigma11[N_I, j + 1, 2] - sigma22[N_I, j + 1, 2], m3=sigma11[N_I, j + 1, 3] - sigma22[N_I, j + 1, 3], m4=sigma11[N_I, j + 1, 4] - sigma22[N_I, j + 1, 4], p1=sigma11[N_I, j + 1, 1] + sigma22[N_I, j + 1, 1], p2=sigma11[N_I, j + 1, 2] + sigma22[N_I, j + 1, 2], p3=sigma11[N_I, j + 1, 3] + sigma22[N_I, j + 1, 3], p4=sigma11[N_I, j + 1, 4] + sigma22[N_I, j + 1, 4], s1=sigma12[N_I, j + 1, 1], s2=sigma12[N_I, j + 1, 2], s3=sigma12[N_I, j + 1, 3], s4=sigma12[N_I, j + 1, 4]&gt;)) where {i in N_I:N_I, j in 1:N_J - 1}, [1:N_I - 1, N_J:N_J] = sum[i, j] (uarear[i, j] * (str_dy_cell_k5&lt;dxT=dxT[i, N_J], dxhy=dxhy[i, N_J], dyT=dyT[i, N_J], dyhx=dyhx[i, N_J], m1=sigma11[i, N_J, 1] - sigma22[i, N_J, 1], m2=sigma11[i, N_J, 2] - sigma22[i, N_J, 2], m3=sigma11[i, N_J, 3] - sigma22[i, N_J, 3], m4=sigma11[i, N_J, 4] - sigma22[i, N_J, 4], p1=sigma11[i, N_J, 1] + sigma22[i, N_J, 1], p2=sigma11[i, N_J, 2] + sigma22[i, N_J, 2], p3=sigma11[i, N_J, 3] + sigma22[i, N_J, 3], p4=sigma11[i, N_J, 4] + sigma22[i, N_J, 4], s1=sigma12[i, N_J, 1], s2=sigma12[i, N_J, 2], s3=sigma12[i, N_J, 3], s4=sigma12[i, N_J, 4]&gt; + str_dy_cell_k6&lt;dxT=dxT[i, N_J], dxhy=dxhy[i, N_J], dyT=dyT[i, N_J], dyhx=dyhx[i, N_J], m1=sigma11[i, N_J, 1] - sigma22[i, N_J, 1], m2=sigma11[i, N_J, 2] - sigma22[i, N_J, 2], m3=sigma11[i, N_J, 3] - sigma22[i, N_J, 3], m4=sigma11[i, N_J, 4] - sigma22[i, N_J, 4], p1=sigma11[i, N_J, 1] + sigma22[i, N_J, 1], p2=sigma11[i, N_J, 2] + sigma22[i, N_J, 2], p3=sigma11[i, N_J, 3] + sigma22[i, N_J, 3], p4=sigma11[i, N_J, 4] + sigma22[i, N_J, 4], s1=sigma12[i, N_J, 1], s2=sigma12[i, N_J, 2], s3=sigma12[i, N_J, 3], s4=sigma12[i, N_J, 4]&gt; + str_dy_cell_k7&lt;dxT=dxT[i + 1, N_J], dxhy=dxhy[i + 1, N_J], dyT=dyT[i + 1, N_J], dyhx=dyhx[i + 1, N_J], m1=sigma11[i + 1, N_J, 1] - sigma22[i + 1, N_J, 1], m2=sigma11[i + 1, N_J, 2] - sigma22[i + 1, N_J, 2], m3=sigma11[i + 1, N_J, 3] - sigma22[i + 1, N_J, 3], m4=sigma11[i + 1, N_J, 4] - sigma22[i + 1, N_J, 4], p1=sigma11[i + 1, N_J, 1] + sigma22[i + 1, N_J, 1], p2=sigma11[i + 1, N_J, 2] + sigma22[i + 1, N_J, 2], p3=sigma11[i + 1, N_J, 3] + sigma22[i + 1, N_J, 3], p4=sigma11[i + 1, N_J, 4] + sigma22[i + 1, N_J, 4], s1=sigma12[i + 1, N_J, 1], s2=sigma12[i + 1, N_J, 2], s3=sigma12[i + 1, N_J, 3], s4=sigma12[i + 1, N_J, 4]&gt; + str_dy_cell_k8&lt;dxT=dxT[i + 1, N_J], dxhy=dxhy[i + 1, N_J], dyT=dyT[i + 1, N_J], dyhx=dyhx[i + 1, N_J], m1=sigma11[i + 1, N_J, 1] - sigma22[i + 1, N_J, 1], m2=sigma11[i + 1, N_J, 2] - sigma22[i + 1, N_J, 2], m3=sigma11[i + 1, N_J, 3] - sigma22[i + 1, N_J, 3], m4=sigma11[i + 1, N_J, 4] - sigma22[i + 1, N_J, 4], p1=sigma11[i + 1, N_J, 1] + sigma22[i + 1, N_J, 1], p2=sigma11[i + 1, N_J, 2] + sigma22[i + 1, N_J, 2], p3=sigma11[i + 1, N_J, 3] + sigma22[i + 1, N_J, 3], p4=sigma11[i + 1, N_J, 4] + sigma22[i + 1, N_J, 4], s1=sigma12[i + 1, N_J, 1], s2=sigma12[i + 1, N_J, 2], s3=sigma12[i + 1, N_J, 3], s4=sigma12[i + 1, N_J, 4]&gt;)) where {i in 1:N_I - 1, j in N_J:N_J}, [N_I:N_I, N_J:N_J] = sum[i, j] (uarear[i, j] * (str_dy_cell_k5&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dy_cell_k6&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dy_cell_k7&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt; + str_dy_cell_k8&lt;dxT=dxT[N_I, N_J], dxhy=dxhy[N_I, N_J], dyT=dyT[N_I, N_J], dyhx=dyhx[N_I, N_J], m1=sigma11[N_I, N_J, 1] - sigma22[N_I, N_J, 1], m2=sigma11[N_I, N_J, 2] - sigma22[N_I, N_J, 2], m3=sigma11[N_I, N_J, 3] - sigma22[N_I, N_J, 3], m4=sigma11[N_I, N_J, 4] - sigma22[N_I, N_J, 4], p1=sigma11[N_I, N_J, 1] + sigma22[N_I, N_J, 1], p2=sigma11[N_I, N_J, 2] + sigma22[N_I, N_J, 2], p3=sigma11[N_I, N_J, 3] + sigma22[N_I, N_J, 3], p4=sigma11[N_I, N_J, 4] + sigma22[N_I, N_J, 4], s1=sigma12[N_I, N_J, 1], s2=sigma12[N_I, N_J, 2], s3=sigma12[N_I, N_J, 3], s4=sigma12[N_I, N_J, 4]&gt;)) where {i in N_I:N_I, j in N_J:N_J})</div>

B-grid lowering rule of the internal-stress y-divergence rewrite target strinty = D(sigma12, x) + D(sigma22, y), the momentum-equation divergence of the internal stress tensor: the Fortran assembles it at each U point as strinty(i,j) = uarear(i,j)*(str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8)) (stepu, ice_dyn_shared.F90:960-961) from the str corner combinations of the stress() routine (ice_dyn_evp.F90:1752-1765). The rule matches the two-term sum by GROUND PATTERN (sigma12/sigma22 literal, so no where-shape constraint is possible or needed, esm-spec 9.6.1/9.6.8) and emits the [i,j] makearray of the four surrounding-cell corner contributions tiled over four i/j position classes (interior [1..N_I-1]x[1..N_J-1], the i = N_I strip, the j = N_J strip, and the (N_I,N_J) corner) with halo-copy closures: on the i = N_I strip the cells (i+1,j) and (i+1,j+1) that supply the k7/k8 str atoms are replaced by (N_I,j)/(N_I,j+1), and on the j = N_J strip the cells (i,j+1) and (i+1,j+1) by (i,N_J)/(i+1,N_J) (the ghost-cell stress copies the edge value, as the Fortran's ice_HaloUpdate of the stress arrays does before stress()). Each region value is a faq over its box that applies the four scalar str atoms str_dy_cell_k5..k8, each bound to its own cell's corner stresses and metrics: the p_k = sigma11 + sigma22, m_k = sigma11 - sigma22, s_k = sigma12 corner combinations (the stressp/stressm/stress12 of the stress() routine, F90:1667-1686) and the cell's dxT/dyT/dxhy/dyhx (ice_dyn_evp.F90:410-416), the whole sum scaled by uarear(i,j), the U-cell inverse area (ice_grid.F90:703). Requires N_I >= 2 and N_J >= 2.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the EVP dynamics stress and stepu.

**Conformance fixtures**

- ast: [`cice_bgrid_strinty`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_strinty) — status `active`

### `transport_flux_divergence`

<p class="rule-meta"><code>op:D</code> <code>order:1</code> <code>bc:halo-copy</code> <code>axes:i,j,ncat</code> <code>priority 10</code></p>

Source: [`grids/cice_bgrid/rules/transport_flux_divergence.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/cice_bgrid/rules/transport_flux_divergence.esm)

**Rewrites** `D(c * u)/Dx + D(c * v)/Dy` at priority 10 — match pattern:

```json
{"op": "+", "args": [{"op": "D", "args": [{"op": "*", "args": ["c", "u"]}], "wrt": "x"}, {"op": "D", "args": [{"op": "*", "args": ["c", "v"]}], "wrt": "y"}]}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>c</code> is a bare field shaped <code>[i, j, ncat]</code>; <code>u</code> is a bare field shaped <code>[i, j]</code>; <code>v</code> is a bare field shaped <code>[i, j]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([2:N_I - 1, 2:N_J - 1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[i - 1, j - 1]&gt;)) / tarea[i, j]) where {i in 2:N_I - 1, j in 2:N_J - 1, k in 1:NCAT}, [1:1, 2:N_J - 1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[1, j], cL=c[1, j, k], cR=c[i, j, k], uN=u[1, j], uS=u[1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[1, j - 1]&gt;)) / tarea[i, j]) where {i in 1:1, j in 2:N_J - 1, k in 1:NCAT}, [N_I:N_I, 2:N_J - 1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[N_I, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[i - 1, j - 1]&gt;)) / tarea[i, j]) where {i in N_I:N_I, j in 2:N_J - 1, k in 1:NCAT}, [2:N_I - 1, 1:1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, 1], cN=c[i, j, k], cS=c[i, 1, k], vE=v[i, 1], vW=v[i - 1, 1]&gt;)) / tarea[i, j]) where {i in 2:N_I - 1, j in 1:1, k in 1:NCAT}, [2:N_I - 1, N_J:N_J, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, N_J, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[i - 1, j - 1]&gt;)) / tarea[i, j]) where {i in 2:N_I - 1, j in N_J:N_J, k in 1:NCAT}, [1:1, 1:1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, 1]&gt; + (-transport_flux_e&lt;HTE=HTE[1, j], cL=c[1, j, k], cR=c[i, j, k], uN=u[1, j], uS=u[1, 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, 1], cN=c[i, j, k], cS=c[i, 1, k], vE=v[i, 1], vW=v[1, 1]&gt;)) / tarea[i, j]) where {i in 1:1, j in 1:1, k in 1:NCAT}, [1:1, N_J:N_J, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[i + 1, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[1, j], cL=c[1, j, k], cR=c[i, j, k], uN=u[1, j], uS=u[1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, N_J, k], cS=c[i, j, k], vE=v[i, j], vW=v[1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[1, j - 1]&gt;)) / tarea[i, j]) where {i in 1:1, j in N_J:N_J, k in 1:NCAT}, [N_I:N_I, 1:1, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[N_I, j, k], uN=u[i, j], uS=u[i, 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, j + 1, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, 1], cN=c[i, j, k], cS=c[i, 1, k], vE=v[i, 1], vW=v[i - 1, 1]&gt;)) / tarea[i, j]) where {i in N_I:N_I, j in 1:1, k in 1:NCAT}, [N_I:N_I, N_J:N_J, 1:NCAT] = sum[i, j, k] ((transport_flux_e&lt;HTE=HTE[i, j], cL=c[i, j, k], cR=c[N_I, j, k], uN=u[i, j], uS=u[i, j - 1]&gt; + (-transport_flux_e&lt;HTE=HTE[i - 1, j], cL=c[i - 1, j, k], cR=c[i, j, k], uN=u[i - 1, j], uS=u[i - 1, j - 1]&gt;) + transport_flux_n&lt;HTN=HTN[i, j], cN=c[i, N_J, k], cS=c[i, j, k], vE=v[i, j], vW=v[i - 1, j]&gt; + (-transport_flux_n&lt;HTN=HTN[i, j - 1], cN=c[i, j, k], cS=c[i, j - 1, k], vE=v[i, j - 1], vW=v[i - 1, j - 1]&gt;)) / tarea[i, j]) where {i in N_I:N_I, j in N_J:N_J, k in 1:NCAT})</div>

B-grid lowering rule of the conservative flux-form tracer divergence rewrite target D(c*u, x) + D(c*v, y) for a per-category tracer c on the T cells advected by the B-grid corner velocities uvel/vvel on the U points: the rule matches the compound two-term sum with params c/u/v as wildcards and where-shape constraints (c: [i,j,ncat], u: [i,j], v: [i,j], esm-spec 9.6.1) and emits the [i,j,ncat] makearray of the conservative cell divergence div = (F_e(i,j,k) - F_e(i-1,j,k) + F_n(i,j,k) - F_n(i,j-1,k)) / tarea(i,j) (the Fortran mass update mm = mm - (mflxe(i,j) - mflxe(i-1,j) + mflxn(i,j) - mflxn(i,j-1))*tarear, ice_transport_remap.F90:3586-3590), where the east/north edge fluxes are the centered corner-average reconstructions F_e = HTE*((uS + uN)/2)*((cL + cR)/2) and F_n = HTN*((vW + vE)/2)*((cS + cN)/2) of the scalar atoms transport_flux_e/n: the edge runs between the two corner U points and its velocity is the average of the two corner velocities times the edge length, exactly the Fortran's edge-area reconstruction (edgearea_e/edgearea_n, ice_transport_remap.F90:658-666, the linear branch), and cL/cR (cS/cN) are the tracer values in the two cells sharing the edge, averaged to the face. The makearray tiles nine i/j position classes (interior [2..N_I-1]x[2..N_J-1], the four strips i = 1, i = N_I, j = 1, j = N_J, and the four corners) with halo-copy closures: every neighbor index that leaves the domain (c, u/v, HTE/HTN at i-1, i+1, j-1, j+1) is bound to the nearest interior value, reproducing the Fortran's ice_HaloUpdate of the transported fields and edge lengths. NOTE: the traced Fortran run integrates the incremental-remap transport (departure-region integrals over the limited linear reconstruction, transport_integrals, ice_transport_remap.F90:3188-3472), whose per-cell flux differs from this centered face average by O(C) slope corrections; the scheme difference is documented with the stage-2 transport pins.

**References**

- Hunke, E. C., et al. (2023). CICE v6.3.22 documentation and user's guide (LA-UR-23-31851), the horizontal transport edge reconstruction.

**Conformance fixtures**

- ast: [`cice_bgrid_transport_flux_divergence`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/cice_bgrid_transport_flux_divergence) — status `active`

