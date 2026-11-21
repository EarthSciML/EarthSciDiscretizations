---
title: "column_nonuniform_1d_tracer"
description: "grids/column_nonuniform_1d with a trailing TRACER axis: fields shaped [lev, tracer] are NTRACER transported constituents (chemical species, aerosol sections) that share ONE column's geometry and coefficients."
esd_kind: "grid"
families: ["column", "finite_volume"]
ops: ["interface_flux", "varcoeff_laplacian"]
bcs: ["prescribed_flux", "supplied_face_flux"]
tags: ["esd:grid", "family:column", "ndim:2", "spacing:nonuniform", "axis:lev", "axis:tracer", "esd:stencil", "family:finite_volume", "grid:column_nonuniform_1d_tracer", "op:interface_flux", "order:2", "stencil_width:2", "op:varcoeff_laplacian", "stencil_width:3", "esd:rule", "axes:lev,tracer", "bc:prescribed_flux", "spec:9.6.8", "bc:supplied_face_flux"]
source: "grids/column_nonuniform_1d_tracer/grid.typ"
---

grids/column_nonuniform_1d with a trailing TRACER axis: fields shaped [lev, tracer] are NTRACER transported constituents (chemical species, aerosol sections) that share ONE column's geometry and coefficients. The vertical axes are exactly the column grid's -- `lev` (NLEV layers, layer 1 the surface) and `lev_nodes` (NLEV+1 interfaces) -- declared here -- with the NLEV metaparameter -- verbatim as on that grid, so they deep-equal and dedupe against grids/column_nonuniform_1d when a model imports from both grids; `tracer` (NTRACER members) is the only new axis. The geometry is the column grid's, unchanged: a consumer supplies the interface array `ze` [lev_nodes] and derives `zc` and `dz` [lev] by applying column_nonuniform_1d's `column_cell_center` / `column_cell_thickness` (import that grid file, or its mesh fixture, beside this grid's rules); this grid adds no geometry, so the per-tracer operators read the same `zc`/`dz` by bare name. Its one template, `tracer_column_total`, is the per-tracer column total sum_k u[k,s] dz[k] that the flux-form operators here conserve tracer by tracer. It is a separate grid rather than an extra index set on column_nonuniform_1d because the grid-first layout (AGENTS.md section 3) puts an operator that takes a different form on a different field shape under its own grid, and rule match-scoping (esm-spec 9.6.1 `where` shapes) already tells a [lev] field from a [lev, tracer] one by axis names, so the two grids' rules are disjoint and a model may import both. A consumer whose species axis has another name imports this grid's rules with `rename: {tracer: <its axis>}` (esm-spec 9.7.7) and binds NTRACER to that axis's size; NLEV is bound as on the column grid. Motivation: one mount of a column tracer-mixing component for every species of an atmospheric-chemistry mechanism (WRF-Chem's vertmx is called once per species with the same exchange coefficients).

## Metaparameters

| Name | Type | Default | Description |
|---|---|---|---|
| `NLEV` | integer | `64` | Number of layers (cells) along lev; the interface axis lev_nodes has NLEV+1 members. Minimum 2 (see the grid description). |
| `NTRACER` | integer | `1` | Number of tracers on the trailing axis. |

All sizes are load-time metaparameters (esm-spec §9.7.6): rebind them at the import edge (or through the loader API) and the same files serve every resolution.

## Consumer-supplied free names

<p class="rule-meta"><code>dz</code> <code>kdudz_bot</code> <code>kdudz_top</code> <code>zc</code></p>

These names are supplied by the consuming model's scope at evaluation (the keyed-factor / free-name geometry contract); they are neither metaparameters (load-time integers) nor template parameters (bound per apply site).

## Index sets

| Name | Kind | Size |
|---|---|---|
| `lev` | interval | \(\mathit{NLEV}\) |
| `lev_nodes` | interval | \(\mathit{NLEV} + 1\) |
| `tracer` | interval | \(\mathit{NTRACER}\) |

## Geometry templates

| Template | Parameters | Definition | Description |
|---|---|---|---|
| `tracer_column_total` | `u` | \(\operatorname{faq}(u; \mathrm{expr} = u_{\mathit{gk}, s} \mathit{dz}_{\mathit{gk}}, \mathrm{output\_idx} = [\text{"s"}], \mathrm{ranges} = {\text{"gk"}: {\text{"from"}: \text{"lev"}}, \text{"s"}: {\text{"from"}: \text{"tracer"}}}, \mathrm{reduce} = \text{"+"})\) | Per-tracer column total sum_k u[k,s] dz[k] of a [lev, tracer] field, shape [tracer]: the finite-volume measure the per-tracer flux-form operators on this grid conserve, so for the prescribed-flux diffusion rule d/dt of it equals kdudz_top[s] - kdudz_bot[s] exactly (telescoping). The same sum as grids/column_nonuniform_1d/stencils/integral_lev_total applied to each column u[:, s]; `dz` is read by bare name as on the column grid. Aggregation index `gk`, not `k`, so it never captures a consuming rule's own aggregate index when inlined. |

## Stencils

Interior-only, match-less templates importing `grid.esm`; the rules below wrap them with boundary-condition face regions into complete rewrite rules.

### `interface_flux_K_lev_tracer`

<p class="rule-meta"><code>op:interface_flux</code> <code>order:2</code> <code>stencil width 2</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d_tracer/stencils/interface_flux_K_lev_tracer.typ`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d_tracer/stencils/interface_flux_K_lev_tracer.typ)

Parameters: `u`, `K`, `n`, `s`

<div class="math-block">$$\operatorname{interface\_flux\_K\_two\_point}(K_{\mathrm{face}} = K_{n}, u_{\mathrm{hi}} = u_{n, s}, u_{\mathrm{lo}} = u_{n - 1, s}, z_{\mathrm{hi}} = \mathit{zc}_{n}, z_{\mathrm{lo}} = \mathit{zc}_{n - 1})$$</div>

The interface flux of grids/column_nonuniform_1d/stencils/interface_flux_K_lev for ONE member s of a trailing tracer axis: G[n,s] = K[n] (u[n,s]-u[n-1,s])/(zc[n]-zc[n-1]). The arithmetic is the column grid's one two-point flux law, grids/column_nonuniform_1d/stencils/interface_flux_K_two_point (imported across grids; it reads no grid names), applied with K_face = K[n], u_lo = u[n-1,s], u_hi = u[n,s], z_lo = zc[n-1], z_hi = zc[n]; this template only fixes which values an interior interface of tracer s reads. The coefficient K and the geometry are shared by every tracer; only u carries the second index. Column s is exactly interface_flux_K_lev applied to u[:,s]. The interface index n and the tracer index s are template parameters bound at each call site (the enclosing aggregate's loop symbols k or k+1 and s in stencils/varcoeff_face_laplacian_lev_interior_tracer, the literal 2 and the metaparameter NLEV in the two boundary layers of rules/varcoeff_face_laplacian_lev_tracer_frame). `zc` is the consumer's free-name layer-centre observed (grids/column_nonuniform_1d/grid.esm).

### `varcoeff_face_laplacian_lev_interior_tracer`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>stencil width 3</code> <code>interior-only (match-less)</code></p>

Source: [`grids/column_nonuniform_1d_tracer/stencils/varcoeff_face_laplacian_lev_interior_tracer.typ`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d_tracer/stencils/varcoeff_face_laplacian_lev_interior_tracer.typ)

Parameters: `u`, `K`

<div class="math-block">$$\left[\sum\frac{\operatorname{interface\_flux\_K\_lev\_tracer}(K = K, n = k + 1, s = s, u = u) - \operatorname{interface\_flux\_K\_lev\_tracer}(K = K, n = k, s = s, u = u)}{\mathit{dz}_{k}}\right]^{\{u, K\}}_{2 \le k \le \mathit{NLEV} - 1, s \in \mathit{tracer}}$$</div>

Interior [2, NLEV-1] of the variable-coefficient vertical diffusion for a [lev, tracer] field: (G[k+1,s] - G[k,s]) / dz[k] for every tracer s, G being stencils/interface_flux_K_lev_tracer. Column s is grids/column_nonuniform_1d/stencils/varcoeff_face_laplacian_lev_interior applied to u[:, s].

## Rules

Complete auto-applied rewrite rules on spatial `D`: imported stencil + boundary conditions in one `makearray` (esm-spec §9.6.8). Import a rule and every matching derivative in your model lowers through it.

### `varcoeff_face_laplacian_lev_flux_bc_tracer`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:prescribed_flux</code> <code>axes:lev,tracer</code> <code>stencil width 3</code> <code>priority 10</code></p>

Source: [`grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer.md`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer.md)

**Rewrites** \(\frac{\partial}{\partial\mathit{lev}} \left(K \frac{\partial u}{\partial\mathit{lev}}\right)\) at priority 10 — match pattern:

```json
{"op": "D", "args": [{"op": "*", "args": ["K", {"op": "D", "args": ["u"], "wrt": "lev"}]}], "wrt": "lev"}
```

<div class="callout callout-scope"><strong>Match scope</strong> (esm-spec §9.6.1 <code>where</code>): fires only when <code>K</code> is a bare field shaped <code>[lev_nodes]</code>; <code>u</code> is a bare field shaped <code>[lev, tracer]</code>. Under import-edge renaming (§9.7.7) the <code>wrt</code> literal and this shape follow the renamed axis together, so the rule can be imported more than once (each instance scoped to its own grid) without a first-declared-wins collision. A consumer differentiating a compound inline expression must bind it to a declared shaped observed first.</div>

**Definition:**

<div class="math-block">$$\operatorname{varcoeff\_face\_laplacian\_lev\_tracer\_frame}(G_{\mathrm{bot}} = \mathit{kdudz\_bot}, G_{\mathrm{top}} = \mathit{kdudz\_top}, K = K, u = u)$$</div>

grids/column_nonuniform_1d/rules/varcoeff_face_laplacian_lev_flux_bc applied to a field with a trailing TRACER axis, u shaped [lev, tracer]: the same compound match D( K*D(u, lev), lev ), the same stencil and the same prescribed-flux closure, column by column, with the coefficient K [lev_nodes] and the geometry shared by every tracer and the two prescribed face fluxes `kdudz_bot` / `kdudz_top` supplied per tracer ([tracer]). The body is this grid's three-region layout rules/varcoeff_face_laplacian_lev_tracer_frame with its face-flux parameters bound to those two free names, so L u[1,s] = ( G[2,s] - kdudz_bot[s] ) / dz[1] and L u[NLEV,s] = ( kdudz_top[s] - G[NLEV,s] ) / dz[NLEV], every interior flux G[n,s] being the column grid's two-point law (grids/column_nonuniform_1d/stencils/interface_flux_K_two_point) applied to column s. Column s of the result is exactly what the 1-D rule gives for u[:, s], so the telescoping identity sum_k dz[k] L u[k,s] = kdudz_top[s] - kdudz_bot[s] holds tracer by tracer. Match-scoped with where {u: {shape: [lev, tracer]}, K: {shape: [lev_nodes]}}, so it is disjoint from the 1-D rule (where u: [lev]) and a model may import both: each fires only on its own field shape (tests/conformance/ast/varcoeff_face_laplacian_lev_flux_bc_tracer imports both and lowers one field of each shape). The tracer axis comes from this grid (grids/column_nonuniform_1d_tracer/grid.typ); a consumer whose axis has another name imports this rule with `rename: {tracer: <name>}` and binds NTRACER. Motivation: one mount of a column tracer-mixing component for every chemical species of a WRF-Chem mechanism (WRF-Chem's vertmx is called once per species with the same exch_h). Second-order in L2 on a smoothly stretched column, as the 1-D rule; MMS convergence with a variable K and unequal per-tracer face fluxes is pinned by problems/heat_column_varcoeff_flux_bc_tracer.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).
- Hong, S.-Y., Noh, Y., and Dudhia, J. (2006). A new vertical diffusion package with an explicit treatment of entrainment processes. Mon. Wea. Rev. 134, 2318-2341. (Surface-flux boundary condition of the YSU scheme.) [`doi:10.1175/MWR3199.1`](https://doi.org/10.1175/MWR3199.1)

**Conformance fixtures**

- ast: [`varcoeff_face_laplacian_lev_flux_bc_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/ast/varcoeff_face_laplacian_lev_flux_bc_tracer) — status `active`
- convergence: [`heat_column_varcoeff_flux_bc_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc_tracer) — status `active`
- convergence: [`heat_column_varcoeff_forced_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_forced_tracer) — status `active`
- simulation: [`heat_column_varcoeff_flux_bc_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_column_varcoeff_flux_bc_tracer) — status `active`
- simulation: [`heat_column_varcoeff_forced_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/simulation/heat_column_varcoeff_forced_tracer) — status `active`

#### Convergence — case [`heat_column_varcoeff_flux_bc_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc_tracer)

Error norms read from the committed golden [`tests/conformance/convergence/heat_column_varcoeff_flux_bc_tracer/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_column_varcoeff_flux_bc_tracer/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 4.452e-05 | — |
| 32 | 1.199e-05 | 1.89 |
| 64 | 3.080e-06 | 1.96 |
| 128 | 7.755e-07 | 1.99 |

Expected order: **2** (± 0.25).

![Convergence of heat_column_varcoeff_flux_bc_tracer (log-log)](/plots/convergence/heat_column_varcoeff_flux_bc_tracer.png)

#### Convergence — case [`heat_column_varcoeff_forced_tracer`](https://github.com/EarthSciML/earthscidiscretizations/tree/main/tests/conformance/convergence/heat_column_varcoeff_forced_tracer)

Error norms read from the committed golden [`tests/conformance/convergence/heat_column_varcoeff_forced_tracer/golden/errors.json`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/tests/conformance/convergence/heat_column_varcoeff_forced_tracer/golden/errors.json) (binding: `julia`; nothing recomputed for display).

| n | L2_error | observed order |
|---|---|---|
| 16 | 2.115e-03 | — |
| 32 | 5.284e-04 | 2.00 |
| 64 | 1.321e-04 | 2.00 |
| 128 | 3.302e-05 | 2.00 |

Expected order: **2** (± 0.25).

![Convergence of heat_column_varcoeff_forced_tracer (log-log)](/plots/convergence/heat_column_varcoeff_forced_tracer.png)

### `varcoeff_face_laplacian_lev_tracer_frame`

<p class="rule-meta"><code>op:varcoeff_laplacian</code> <code>order:2</code> <code>bc:supplied_face_flux</code> <code>axes:lev,tracer</code> <code>stencil width 3</code></p>

Source: [`grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_tracer_frame.typ`](https://github.com/EarthSciML/earthscidiscretizations/blob/main/grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_tracer_frame.typ)

**Discretization** — one `makearray` tiling the output axes; boundary conditions live in the face regions of this same rule (esm-spec §9.6.8; later regions overwrite earlier, §4.3.2). Imported stencils appear as template applications `name<param=value, ...>`, defined in the Stencils section of this page:

<div class="math-block">$$\left[\begin{array}{l|l} \operatorname{varcoeff\_face\_laplacian\_lev\_interior\_tracer}(K = K, u = u) &amp; [2, \mathit{NLEV} - 1] \times [1, \mathit{NTRACER}] \\ \left[\sum\frac{\operatorname{interface\_flux\_K\_lev\_tracer}(K = K, n = 2, s = s, u = u) - G_{\mathrm{bot}}[s]}{\mathit{dz}_{1}}\right]^{\{u, K\}}_{1, s \in \mathit{tracer}} &amp; [1, 1] \times [1, \mathit{NTRACER}] \\ \left[\sum\frac{G_{\mathrm{top}}[s] - \operatorname{interface\_flux\_K\_lev\_tracer}(K = K, n = \mathit{NLEV}, s = s, u = u)}{\mathit{dz}_{\mathit{NLEV}}}\right]^{\{u, K\}}_{1, s \in \mathit{tracer}} &amp; [\mathit{NLEV}, \mathit{NLEV}] \times [1, \mathit{NTRACER}] \end{array}\right]$$</div>

The three-region layout of the conservative variable-coefficient vertical diffusion d/dz(K du/dz) for a field with a trailing TRACER axis, u shaped [lev, tracer], with the coefficient K [lev_nodes] on the interfaces and the geometry shared by every tracer, whatever closes the two domain faces: a makearray tiling [lev, tracer] whose interior [2, NLEV-1] x [1, NTRACER] is stencils/varcoeff_face_laplacian_lev_interior_tracer, (G[k+1,s] - G[k,s]) / dz[k], and whose two boundary layers are L u[1,s] = ( G[2,s] - G_bot[s] ) / dz[1] and L u[NLEV,s] = ( G_top[s] - G[NLEV,s] ) / dz[NLEV], every interior flux G[n,s] being stencils/interface_flux_K_lev_tracer (the column grid's two-point law grids/column_nonuniform_1d/stencils/interface_flux_K_two_point applied to column s). Column s is exactly the 1-D layout of grids/column_nonuniform_1d applied to u[:, s]. The two domain-face fluxes are PARAMETERS, [tracer]-shaped arrays read element by element as G_bot[s] (the diffusive flux K du/dz of tracer s through the surface face z = ze[1]) and G_top[s] (through the top face z = ze[NLEV+1]), positive in the +lev (upward) direction and bound by the applying match rule to whatever its boundary condition makes them -- the per-tracer free names kdudz_bot / kdudz_top in rules/varcoeff_face_laplacian_lev_flux_bc_tracer. So the interior stencil, the region layout and the boundary-layer flux differences have one spelling on this grid, and the telescoping identity sum_k dz[k] L u[k,s] = G_top[s] - G_bot[s] holds tracer by tracer for every closure. A match-less rule template (it never fires on its own); the bc: tag names what it takes, not a closure it imposes. The two boundary layers are faqs with output_idx [1, s] -- an explicit singleton lev dimension -- because the be1fafe6d Rust binding fills the WHOLE array with NaN when a rank-1 [tracer] value is assigned to a [1, NTRACER] region, although esm-spec 4.3.2 says singleton dimensions are excluded from the shape match. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). Minimum extent NLEV = 2, as on the column grid.

**References**

- LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).

**Conformance fixtures**

- none yet

