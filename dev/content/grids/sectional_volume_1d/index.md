---
title: "sectional_volume_1d"
description: "SECTIONAL particle size distribution in single-particle VOLUME v, the size grid of a sectional aerosol model (MOSAIC, MADRID, CMAQ-sectional, ...): NBIN contiguous sections `bin` with the NBIN+1 section edges `bin_nodes`, carried independently in each of NCELL spatial cells `cell` (a column level, a grid box -- whatever the consumer's own cell axis is; rename it at the import edge, e.g."
esd_kind: "grid"
families: ["finite_volume", "sectional"]
ops: ["advection"]
bcs: ["closed"]
tags: ["esd:grid", "family:sectional", "ndim:2", "spacing:nonuniform", "axis:cell", "axis:bin", "axis:tracer", "esd:stencil", "family:finite_volume", "grid:sectional_volume_1d", "op:advection", "order:2", "scheme:linear-discrete", "esd:rule", "scheme:upwind", "bc:closed", "form:flux", "spec:9.6.8", "axes:cell,bin", "axes:cell,bin,tracer"]
source: "grids/sectional_volume_1d/grid.esm"
---

SECTIONAL particle size distribution in single-particle VOLUME v, the size grid of a sectional aerosol model (MOSAIC, MADRID, CMAQ-sectional, ...): NBIN contiguous sections `bin` with the NBIN+1 section edges `bin_nodes`, carried independently in each of NCELL spatial cells `cell` (a column level, a grid box -- whatever the consumer's own cell axis is; rename it at the import edge, e.g. {cell: lev}). Section j spans [vb[j], vb[j+1]], so consecutive sections share an edge and the sections tile [vb[1], vb[NBIN+1]]. Each section carries INTEGRATED moments of the number distribution n(v): the number N[c,j] = integral of n dv over the section and the volume V[c,j] = integral of v n dv, plus any number of volume-carried quantities (component masses, water) on the trailing axis `tracer` (NTRACER), which are assumed uniformly mixed over the section's particles. Because the moments are section integrals, not densities, a flux-form operator on this grid returns the INTEGRATED divergence: the net flux through a section's two edges, with no division by the section width. The geometry is consumer-supplied (the integers-are-metaparameters, reals-are-consumer-supplied convention): a consuming model supplies the edge volumes `vb`, shape [bin_nodes], strictly increasing with the bin index, as a free name resolved in its own scope. Nothing here assumes a particular spacing; MOSAIC's sections are geometric (vb[j+1]/vb[j] = 8, a doubling of diameter).

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NCELL` | integer | `1` | Number of spatial cells, each carrying its own size distribution. |
| `NBIN` | integer | `8` | Number of size sections; the edge axis bin_nodes has NBIN+1 members. Minimum 2. |
| `NTRACER` | integer | `1` | Number of volume-carried quantities on the trailing axis. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>vb</code></p>

These names are supplied by the consuming model's scope at evaluation (the keyed-factor / free-name geometry contract); they are neither metaparameters (load-time integers) nor template parameters (bound per apply site).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `cell` | interval | `NCELL` |
| `bin` | interval | `NBIN` |
| `bin_nodes` | interval | `{'op': '+', 'args': ['NBIN', 1]}` |
| `tracer` | interval | `NTRACER` |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `sectional_total` | `N` | `sum[c] (N[c, gj]) where {c in cell, gj in bin}` | Per-cell total sum_j N[c,j] over the size sections of a [cell, bin] section moment (number, or volume), shape [cell]: the quantity every flux-form operator on this grid conserves exactly, since each lowers to edge fluxes that cancel in pairs between neighbouring sections and vanish at the two closed ends of the size range. Aggregation index `gj`, so it never captures a consuming rule's own index when inlined. problems/growth_sectional_lr2m_mms asserts it (its Ntot observed) both with the distribution inside the size range and after it has piled up against both closed ends. |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `sectional_lr2m_density_bottom`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_density_bottom.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_density_bottom.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(sectional_lr2m_interior_mean&lt;N=N, V=V, c=c, j=j&gt; and sectional_lr2m_support_bottom&lt;N=N, V=V, c=c, j=j&gt; &lt;= vb[j], (1 - sectional_lr2m_slope&lt;N=N, V=V, c=c, j=j&gt; / 2) * N[c, j] / (sectional_lr2m_support_top&lt;N=N, V=V, c=c, j=j&gt; - sectional_lr2m_support_bottom&lt;N=N, V=V, c=c, j=j&gt;), 0)</div>

Number density per unit particle volume at the LOWER edge v = vb[j] of section j, N (1 - b/2)/(hi - lo): the reconstructed linear density at z = 0, zero when the support starts above the edge (an upper-third triangle) or the section has no interior mean.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_density_top`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_density_top.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_density_top.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(sectional_lr2m_interior_mean&lt;N=N, V=V, c=c, j=j&gt; and sectional_lr2m_support_top&lt;N=N, V=V, c=c, j=j&gt; &gt;= vb[j + 1], (1 + sectional_lr2m_slope&lt;N=N, V=V, c=c, j=j&gt; / 2) * N[c, j] / (sectional_lr2m_support_top&lt;N=N, V=V, c=c, j=j&gt; - sectional_lr2m_support_bottom&lt;N=N, V=V, c=c, j=j&gt;), 0)</div>

Number density per unit particle volume at the UPPER edge v = vb[j+1] of section j, N (1 + b/2)/(hi - lo): the reconstructed linear density at z = 1, zero when the support stops short of the edge (a lower-third triangle) or the section has no interior mean (sectional_lr2m_interior_mean).

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_growth_number_interior`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_growth_number_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_growth_number_interior.esm)

Parameters: `lam`, `N`, `V`

<div class="math-block">sum[c, j] (sectional_lr2m_number_flux_up&lt;N=N, V=V, c=c, j=j, lam=lam&gt; + sectional_lr2m_number_flux_down&lt;N=N, V=V, c=c, j=j, lam=lam&gt; - (sectional_lr2m_number_flux_up&lt;N=N, V=V, c=c, j=j - 1, lam=lam&gt; + sectional_lr2m_number_flux_down&lt;N=N, V=V, c=c, j=j + 1, lam=lam&gt;)) where {c in cell, j in 2:NBIN - 1}</div>

INTERIOR sections j = 2 .. NBIN-1 of the integrated growth divergence of the number moment: what leaves section j through its two edges (its own upward and downward fluxes) minus what enters from its neighbours (section j-1's upward flux and section j+1's downward flux). Every edge flux is one reference to sectional_lr2m_number_flux_up / _down, so adjacent sections reuse the identical edge expression and the total number is conserved exactly.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_growth_tracer_interior`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_growth_tracer_interior.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_growth_tracer_interior.esm)

Parameters: `lam`, `N`, `V`, `m`

<div class="math-block">sum[c, j, s] (sectional_lr2m_tracer_flux_up&lt;N=N, V=V, c=c, j=j, lam=lam, m=m, s=s&gt; + sectional_lr2m_tracer_flux_down&lt;N=N, V=V, c=c, j=j, lam=lam, m=m, s=s&gt; - (sectional_lr2m_tracer_flux_up&lt;N=N, V=V, c=c, j=j - 1, lam=lam, m=m, s=s&gt; + sectional_lr2m_tracer_flux_down&lt;N=N, V=V, c=c, j=j + 1, lam=lam, m=m, s=s&gt;)) where {c in cell, j in 2:NBIN - 1, s in tracer}</div>

INTERIOR sections j = 2 .. NBIN-1 of the integrated growth divergence of a volume-carried quantity m (shape [cell, bin, tracer]): outflow through both edges minus inflow from the neighbours, each edge flux one reference to sectional_lr2m_tracer_flux_up / _down.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_high_triangle`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_high_triangle.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_high_triangle.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">sectional_lr2m_position&lt;N=N, V=V, c=c, j=j&gt; &gt;= 2 / 3</div>

True when the mean volume of section j sits in the UPPER third of the section, x >= 2/3 with x = sectional_lr2m_position: the reconstruction is cut back to a triangle on [vb[j+1] - 3 (vb[j+1] - vbar), vb[j+1]] that vanishes at its bottom (sectional_lr2m_support_bottom). Exclusive of sectional_lr2m_low_triangle (x <= 1/3), so the two cut-backs never compete. MOSAIC's form, vbar/vb[j] >= 0.9999 + 2 gamma/3, carries an absolute ratio margin that does not shrink with the section width; see sectional_lr2m_low_triangle for why the position form is used and why it reproduces MOSAIC wherever MOSAIC's margins do not overlap.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_interior_mean`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_interior_mean.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_interior_mean.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">N[c, j] &gt; 0 and sectional_lr2m_position&lt;N=N, V=V, c=c, j=j&gt; &gt; 0.01 and sectional_lr2m_position&lt;N=N, V=V, c=c, j=j&gt; &lt; 0.99</div>

True when section j holds particles (N > 0) and their mean volume lies inside the section by more than 1 % of its width at each end, 0.01 < x < 0.99 with x = sectional_lr2m_position. Otherwise the section has no linear reconstruction and carries nothing across its edges at an infinitesimal growth step: MOSAIC moves such a section's whole content to the section that holds its mean ('moving center', iforce_movecenter 12/13), a finite relocation that is not part of this rate operator -- a mean outside its own section is a state this operator never produces, so a consumer that can create one (coagulation, emissions, nucleation) must project it back itself. An interior mean implies V = N vbar > 0, which is what lets sectional_lr2m_tracer_per_volume divide by the donor volume under this same predicate.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_low_triangle`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_low_triangle.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_low_triangle.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">sectional_lr2m_position&lt;N=N, V=V, c=c, j=j&gt; &lt;= 1 / 3</div>

True when the mean volume of section j sits in the LOWER third of the section, x <= 1/3 with x = sectional_lr2m_position the mean's fractional position: a linear number density on [vb[j], vb[j+1]] with that mean would go negative at the upper edge, so the reconstruction is cut back to a triangle on [vb[j], vb[j] + 3 (vbar - vb[j])] that vanishes at its top (sectional_lr2m_support_top). The test has no tolerance, and so is exclusive of sectional_lr2m_high_triangle: the reconstruction is continuous across x = 1/3 (the full-section profile there already vanishes at the top), so which branch an exact tie takes is immaterial. MOSAIC (move_sections_calc_masnumadvect) writes the same test as vbar/vb[j] <= 1.0001 + gamma/3 with gamma = vb[j+1]/vb[j] - 1, an absolute 1e-4 on the volume RATIO; in position units that is a margin of 1e-4/gamma, which is harmless for MOSAIC's own sections (gamma = 7) but swallows the whole upper third once gamma < ~6e-4 and then yields a negative density at the lower edge. For non-overlapping margins the two forms give the same reconstruction, since inside MOSAIC's margin its clamped support is the full section.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_number_flux_down`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_number_flux_down.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_number_flux_down.esm)

Parameters: `lam`, `N`, `V`, `c`, `j`

<div class="math-block">max((-lam[c, j]), 0) * vb[j] * sectional_lr2m_density_bottom&lt;N=N, V=V, c=c, j=j&gt;</div>

Number flux DOWN through the lower edge of section j carried by section j's own shrinkage: max(-lam[c,j], 0) vb[j] n_bottom. Upwind: a shrinking section (lam < 0) is the donor of its lower edge; a growing one sends nothing down.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_number_flux_up`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_number_flux_up.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_number_flux_up.esm)

Parameters: `lam`, `N`, `V`, `c`, `j`

<div class="math-block">max(lam[c, j], 0) * vb[j + 1] * sectional_lr2m_density_top&lt;N=N, V=V, c=c, j=j&gt;</div>

Number flux UP through the upper edge of section j carried by section j's own growth: G(v) n(v) at v = vb[j+1] with the growth law G(v) = lam v of the section, i.e. max(lam[c,j], 0) vb[j+1] n_top. Upwind: a growing section (lam > 0) is the donor of its upper edge; a shrinking one sends nothing up. lam is the section's relative volume growth rate (1/V dV/dt at fixed number), shape [cell, bin].

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

### `sectional_lr2m_position`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_position.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_position.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">(sectional_mean_volume&lt;N=N, V=V, c=c, j=j&gt; - vb[j]) / (vb[j + 1] - vb[j])</div>

Position of the mean volume of section j within the section, as a fraction of its width: x = (vbar - vb[j]) / (vb[j+1] - vb[j]), 0 at the lower edge and 1 at the upper one (vbar = sectional_mean_volume). Every branch test of the two-moment reconstruction is a test on x -- the lower and upper thirds (sectional_lr2m_low_triangle / _high_triangle), the 1 % edge band (sectional_lr2m_interior_mean) -- and the slope of the full-section profile is 12 (x - 1/2) (sectional_lr2m_slope). Measuring the position in units of the section width keeps those tests independent of how wide the section is relative to its volume, so they hold on arbitrarily fine size grids.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_slope`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_slope.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_slope.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(sectional_lr2m_low_triangle&lt;N=N, V=V, c=c, j=j&gt;, -2, ifelse(sectional_lr2m_high_triangle&lt;N=N, V=V, c=c, j=j&gt;, 2, 12 * (sectional_lr2m_position&lt;N=N, V=V, c=c, j=j&gt; - 1 / 2)))</div>

Slope b of the linear number density of section j on its support [lo, hi] (sectional_lr2m_support_bottom / _top), in the normalised coordinate z = (v - lo)/(hi - lo): n(z) = N (a + b z)/(hi - lo) with a = 1 - b/2, so that the density integrates to N, and b = 12 ((vbar - lo)/(hi - lo) - 1/2), so that its first moment is the section volume V = N vbar. On the full section (lo, hi) = (vb[j], vb[j+1]) the normalised mean is the position x itself (sectional_lr2m_position), so b = 12 (x - 1/2); on a cut-back support the normalised mean is exactly 1/3 (the lower-third triangle, b = -2, vanishing at its top) or 2/3 (the upper-third triangle, b = +2, vanishing at its bottom). The slope is written by branch so that it reads the position once and the supports not at all. MOSAIC's dumbb.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_support_bottom`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_support_bottom.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_support_bottom.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(sectional_lr2m_high_triangle&lt;N=N, V=V, c=c, j=j&gt;, vb[j + 1] + 3 * (sectional_mean_volume&lt;N=N, V=V, c=c, j=j&gt; - vb[j + 1]), vb[j])</div>

Lower end of the support of the reconstructed number density of section j: vb[j], or, for a mean in the upper third (sectional_lr2m_high_triangle), vb[j+1] + 3 (vbar - vb[j+1]) -- the bottom of the triangle that has the section's mean and vanishes there, which lies inside the section because x >= 2/3. A support that starts above vb[j] puts no particles at the lower edge, so evaporation carries nothing across it.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_support_top`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_support_top.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_support_top.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(sectional_lr2m_low_triangle&lt;N=N, V=V, c=c, j=j&gt;, vb[j] + 3 * (sectional_mean_volume&lt;N=N, V=V, c=c, j=j&gt; - vb[j]), vb[j + 1])</div>

Upper end of the support of the reconstructed number density of section j: vb[j+1], or, for a mean in the lower third (sectional_lr2m_low_triangle), vb[j] + 3 (vbar - vb[j]) -- the top of the triangle that has the section's mean and vanishes there, which lies inside the section because x <= 1/3. A support that ends below vb[j+1] puts no particles at the upper edge, so growth carries nothing across it.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_tracer_flux_down`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_flux_down.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_flux_down.esm)

Parameters: `lam`, `N`, `V`, `m`, `c`, `j`, `s`

<div class="math-block">vb[j] * sectional_lr2m_number_flux_down&lt;N=N, V=V, c=c, j=j, lam=lam&gt; * sectional_lr2m_tracer_per_volume&lt;N=N, V=V, c=c, j=j, m=m, s=s&gt;</div>

Flux of a volume-carried quantity m DOWN through the lower edge of section j: vb[j] F_N m[c,j,s]/V[c,j] with F_N = sectional_lr2m_number_flux_down and m/V = sectional_lr2m_tracer_per_volume. Zero whenever the number flux is, under the one predicate (sectional_lr2m_interior_mean) both read.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_tracer_flux_up`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_flux_up.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_flux_up.esm)

Parameters: `lam`, `N`, `V`, `m`, `c`, `j`, `s`

<div class="math-block">vb[j + 1] * sectional_lr2m_number_flux_up&lt;N=N, V=V, c=c, j=j, lam=lam&gt; * sectional_lr2m_tracer_per_volume&lt;N=N, V=V, c=c, j=j, m=m, s=s&gt;</div>

Flux of a volume-carried quantity m (shape [cell, bin, tracer]: a component mass, water, ...) UP through the upper edge of section j: every particle crossing carries the edge volume vb[j+1], and the quantity is uniformly mixed over the donor section's volume, so the flux is vb[j+1] F_N m[c,j,s]/V[c,j] with F_N = sectional_lr2m_number_flux_up and m/V = sectional_lr2m_tracer_per_volume. Zero whenever the number flux is, under the one predicate (sectional_lr2m_interior_mean) both read.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_lr2m_tracer_per_volume`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_per_volume.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_lr2m_tracer_per_volume.esm)

Parameters: `N`, `V`, `m`, `c`, `j`, `s`

<div class="math-block">m[c, j, s] / ifelse(sectional_lr2m_interior_mean&lt;N=N, V=V, c=c, j=j&gt;, V[c, j], 1)</div>

Amount of a volume-carried quantity m (shape [cell, bin, tracer]: a component mass, water, ...) per unit particle volume in section j, m[c,j,s] / V[c,j]: the quantity is uniformly mixed over the section's particles, so every particle volume crossing an edge out of section j carries this much of it (sectional_lr2m_tracer_flux_up / _down). It is guarded by the SAME predicate that gates the number flux it multiplies, sectional_lr2m_interior_mean, so the two zero conditions coincide: wherever the number flux can be nonzero the mean is interior and V = N vbar > 0, and everywhere else the denominator is 1 and the (zero) number flux makes the product zero. The denominator is therefore finite on every path, including the untaken branch an AD/taped binding differentiates.

**References**

- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).

### `sectional_mean_volume`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/sectional_volume_1d/stencils/sectional_mean_volume.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/stencils/sectional_mean_volume.esm)

Parameters: `N`, `V`, `c`, `j`

<div class="math-block">ifelse(N[c, j] &gt; 0, V[c, j] / ifelse(N[c, j] &gt; 0, N[c, j], 1), vb[j])</div>

Mean single-particle volume vbar = V[c,j] / N[c,j] of section j in cell c, from the section's integrated volume V and number N (both shape [cell, bin]); the one quantity the two-moment reconstruction of the sibling stencils is built from. An EMPTY section (N <= 0, e.g. every section outside a compactly supported distribution) has no mean, and gets vbar = vb[j], its lower edge, which sectional_lr2m_interior_mean rejects, so the section carries nothing. The guard is applied twice -- ifelse(N > 0, V / ifelse(N > 0, N, 1), vb[j]) -- so the division is finite on BOTH branches: a binding that differentiates through ifelse (a taped/AD Jacobian for an implicit integrator) multiplies the untaken branch's derivative by zero, and a 0/0 there would turn that zero into NaN and poison the Jacobian.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `growth_lr2m_D_bin_closed`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>bc:closed</code> <code>axes:cell,bin</code> <code>priority 10</code></p>

Source: [`grids/sectional_volume_1d/rules/growth_lr2m_D_bin_closed.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/rules/growth_lr2m_D_bin_closed.esm)

**Rewrites** `D(lam * N, V)` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["lam", "N"]}, "V"], "wrt": "bin"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>lam</code> is a bare field shaped <code>[cell, bin]</code>; <code>N</code> is a bare field shaped <code>[cell, bin]</code>; <code>V</code> is a bare field shaped <code>[cell, bin]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NCELL, 2:NBIN - 1] = sectional_lr2m_growth_number_interior&lt;N=N, V=V, lam=lam&gt;, [1:NCELL, 1:1] = sum[c, j] (sectional_lr2m_number_flux_up&lt;N=N, V=V, c=c, j=j, lam=lam&gt; - sectional_lr2m_number_flux_down&lt;N=N, V=V, c=c, j=j + 1, lam=lam&gt;) where {c in cell, j in 1:1}, [1:NCELL, NBIN:NBIN] = sum[c, j] (sectional_lr2m_number_flux_down&lt;N=N, V=V, c=c, j=j, lam=lam&gt; - sectional_lr2m_number_flux_up&lt;N=N, V=V, c=c, j=j - 1, lam=lam&gt;) where {c in cell, j in NBIN:NBIN})</div>

Integrated divergence, over each SECTION of a sectional size distribution, of the condensational-growth flux of particle NUMBER, int_section d(G n)/dv dv, for the growth law G(v) = lam v (each section's particles grow at that section's uniform relative rate lam, the volume-proportional growth of a sectional model in which the condensed volume of a section is shared over its particles in proportion to their volume). The consumer writes the number-moment equation as dN/dt = -D(lam*N, V, wrt: bin) + ..., where N [cell, bin] is the section NUMBER and the trailing operand V [cell, bin] (esm-spec 9.6.8 'trailing operands on D') is the section dry VOLUME moment the reconstruction reads; lam [cell, bin] is 1/V dV/dt of the section's own condensation (zero for a section the consumer considers inactive). The lowered value per section is outflow through both edges minus inflow from the neighbours, telescoping, so sum_j dN_j/dt = 0 to round-off. CLOSED size range: nothing grows out of the top section nor shrinks out of the bottom one (the first and last regions drop the outer edge fluxes), so particles growing past vb[NBIN+1] stay in section NBIN and shrinking ones below vb[1] stay in section 1, as in MOSAIC. RECONSTRUCTION AND EDGE FLUX. Within each section the number density is linear in v on a support [lo, hi] inside the section, fixed by the two moments (N, V) -- MOSAIC's 'linear-discrete' / mass-number-advection reconstruction (move_sections_calc_masnumadvect): the full section when the mean volume vbar = V/N lies in the middle third, a triangle that vanishes at the far end when it lies in the lower or upper third (sectional_lr2m_support_top / _bottom), and no reconstruction (zero edge flux) when vbar lies within 1 % of the section's width of either edge (sectional_lr2m_interior_mean). The edge flux is upwind: the growth flux G n at an edge is evaluated in the section that grows into it (lam > 0 sends particles up through the upper edge) or shrinks into it (lam < 0 sends them down through the lower edge), G(v_edge) n(v_edge) = |lam| v_edge n_edge (sectional_lr2m_number_flux_up / _down). WHY THIS IS THE RIGHT OPERATOR. MOSAIC's movesect is fully discrete in time: after a step's growth it maps every particle of section j from v to s v (s = vbar_after/vbar_before, the section's volume growth factor) and integrates the linear density over the pulled-back section edges. As the growth step goes to zero (s = 1 + lam dt), the fraction of number carried across the upper edge is (1 + b/2) N (vb[j+1] - vb[j+1]/s)/(hi - lo) = lam dt vb[j+1] n_top, and that limit is this flux; verified against a real64 replay of WRF's own move_sections on 16 recorded WRF-Chem column states (59 levels x 8 sections), with each step's growth scaled to fractions 1e-2 and 2e-2 and Richardson-extrapolated: all 98,176 section rates (number, ten components, two waters) agree within 1e-5 relative or 1e-6 of the largest rate of their quantity. MOSAIC's finite-step relocations ('moving center': a section whose mean has left it is moved whole; and the number conform that clips a section's mean back inside its edges) are projections, not rates, and are not part of this operator. ORDER. Second order on smooth data: the two-moment linear reconstruction is exact for a density linear within the section, so the upwind edge value carries an O(h^2) error, as in an unlimited MUSCL scheme; the triangle cut-backs and the 1 % edge band act as a limiter and degrade it locally where they engage (at the edges of a distribution's support). Measured on problems/growth_sectional_lr2m_mms (a compactly supported C2 bump on a geometric grid, growing and shrinking): pairwise observed orders 1.77, 2.13, 2.24 over NBIN = 32..256 (median 2.1; tests/conformance/convergence/growth_sectional_lr2m_mms). The same problem asserts the telescoping: each cell's total number stays at its initial value to round-off, also after the whole distribution has piled up against both closed ends. Requires NBIN >= 2.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

**Conformance fixtures**

- ast: [`growth_lr2m_D_bin_closed`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/growth_lr2m_D_bin_closed) — status `active`
- convergence: [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_mms) — status `active`
- convergence: [`growth_sectional_lr2m_tracer_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms) — status `active`
- simulation: [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/growth_sectional_lr2m_mms) — status `active`

#### Convergence — case [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_mms)

Error norms read from the committed golden [`tests/conformance/convergence/growth_sectional_lr2m_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/growth_sectional_lr2m_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 32 | 3.138e-02 | — |
| 64 | 9.226e-03 | 1.77 |
| 128 | 2.109e-03 | 2.13 |
| 256 | 4.452e-04 | 2.24 |

Expected order: **2** (± 0.25).

![Convergence of growth_sectional_lr2m_mms (log-log)](/plots/convergence/growth_sectional_lr2m_mms.png)

#### Convergence — case [`growth_sectional_lr2m_tracer_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms)

Error norms read from the committed golden [`tests/conformance/convergence/growth_sectional_lr2m_tracer_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 32 | 1.713e-02 | — |
| 64 | 4.926e-03 | 1.80 |
| 128 | 1.349e-03 | 1.87 |
| 256 | 2.727e-04 | 2.31 |

Expected order: **2** (± 0.25).

![Convergence of growth_sectional_lr2m_tracer_mms (log-log)](/plots/convergence/growth_sectional_lr2m_tracer_mms.png)

### `growth_lr2m_tracer_D_bin_closed`

<p class="rule-meta"><code>op:advection</code> <code>order:2</code> <code>bc:closed</code> <code>axes:cell,bin,tracer</code> <code>priority 10</code></p>

Source: [`grids/sectional_volume_1d/rules/growth_lr2m_tracer_D_bin_closed.esm`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/sectional_volume_1d/rules/growth_lr2m_tracer_D_bin_closed.esm)

**Rewrites** `D(lam * m, N, V)` at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["lam", "m"]}, "N", "V"], "wrt": "bin"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>lam</code> is a bare field shaped <code>[cell, bin]</code>; <code>m</code> is a bare field shaped <code>[cell, bin, tracer]</code>; <code>N</code> is a bare field shaped <code>[cell, bin]</code>; <code>V</code> is a bare field shaped <code>[cell, bin]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">makearray([1:NCELL, 2:NBIN - 1, 1:NTRACER] = sectional_lr2m_growth_tracer_interior&lt;N=N, V=V, lam=lam, m=m&gt;, [1:NCELL, 1:1, 1:NTRACER] = sum[c, j, s] (sectional_lr2m_tracer_flux_up&lt;N=N, V=V, c=c, j=j, lam=lam, m=m, s=s&gt; - sectional_lr2m_tracer_flux_down&lt;N=N, V=V, c=c, j=j + 1, lam=lam, m=m, s=s&gt;) where {c in cell, j in 1:1, s in tracer}, [1:NCELL, NBIN:NBIN, 1:NTRACER] = sum[c, j, s] (sectional_lr2m_tracer_flux_down&lt;N=N, V=V, c=c, j=j, lam=lam, m=m, s=s&gt; - sectional_lr2m_tracer_flux_up&lt;N=N, V=V, c=c, j=j - 1, lam=lam, m=m, s=s&gt;) where {c in cell, j in NBIN:NBIN, s in tracer})</div>

Companion of growth_lr2m_D_bin_closed for any quantity carried by the particles' VOLUME -- component masses, aerosol water, a dry-volume tracer -- on the trailing axis `tracer`: the consumer writes dm/dt = -D(lam*m, N, V, wrt: bin) + ... with m [cell, bin, tracer] and the section number N and dry volume V [cell, bin] as trailing operands (esm-spec 9.6.8). The particles crossing an edge each carry the edge volume, and m is uniformly mixed over its donor section's volume, so the flux of m through an edge is v_edge F_N m_donor/V_donor, with F_N the number flux of growth_lr2m_D_bin_closed (sectional_lr2m_tracer_flux_up / _down). That is MOSAIC's apply_moves, which moves every component and the water by the same volume fraction xferfracvol, in the dt -> 0 limit. Telescoping per tracer (sum_j dm_j/dt = 0), closed at both ends of the size range exactly as the number rule. The arity (three operands) keeps the two rules disjoint: the number rule's two-operand pattern cannot match this node, nor this pattern that one. Second order on smooth data, as the number rule whose number flux it carries: measured on its own on problems/growth_sectional_lr2m_mms (the volume-carried m at t = 0.1), pairwise observed orders 1.80, 1.87, 2.31 over NBIN = 32..256 (median 1.9; tests/conformance/convergence/growth_sectional_lr2m_tracer_mms). Requires NBIN >= 2.

**References**

- Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008). Model for Simulating Aerosol Interactions and Chemistry (MOSAIC). Journal of Geophysical Research, 113, D13204. [`doi:10.1029/2007JD008782`](https://doi.org/10.1029/2007JD008782)
- WRF-Chem v4.8.0 chem/module_mosaic_movesect.F, subroutines move_sections_calc_masnumadvect and move_sections_apply_moves (R. C. Easter), the 'linear-discrete' / mass-number-advection section-moving algorithm (msectional = 20).
- Gelbard, F., Tambour, Y., & Seinfeld, J. H. (1980). Sectional representations for simulating aerosol dynamics. Journal of Colloid and Interface Science, 76(2), 541-556. [`doi:10.1016/0021-9797(80)90394-X`](https://doi.org/10.1016/0021-9797%2880%2990394-X)

**Conformance fixtures**

- ast: [`growth_lr2m_tracer_D_bin_closed`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/growth_lr2m_tracer_D_bin_closed) — status `active`
- convergence: [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_mms) — status `active`
- convergence: [`growth_sectional_lr2m_tracer_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms) — status `active`
- simulation: [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/growth_sectional_lr2m_mms) — status `active`

#### Convergence — case [`growth_sectional_lr2m_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_mms)

Error norms read from the committed golden [`tests/conformance/convergence/growth_sectional_lr2m_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/growth_sectional_lr2m_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 32 | 3.138e-02 | — |
| 64 | 9.226e-03 | 1.77 |
| 128 | 2.109e-03 | 2.13 |
| 256 | 4.452e-04 | 2.24 |

Expected order: **2** (± 0.25).

![Convergence of growth_sectional_lr2m_mms (log-log)](/plots/convergence/growth_sectional_lr2m_mms.png)

#### Convergence — case [`growth_sectional_lr2m_tracer_mms`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms)

Error norms read from the committed golden [`tests/conformance/convergence/growth_sectional_lr2m_tracer_mms/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/growth_sectional_lr2m_tracer_mms/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 32 | 1.713e-02 | — |
| 64 | 4.926e-03 | 1.80 |
| 128 | 1.349e-03 | 1.87 |
| 256 | 2.727e-04 | 2.31 |

Expected order: **2** (± 0.25).

![Convergence of growth_sectional_lr2m_tracer_mms (log-log)](/plots/convergence/growth_sectional_lr2m_tracer_mms.png)

