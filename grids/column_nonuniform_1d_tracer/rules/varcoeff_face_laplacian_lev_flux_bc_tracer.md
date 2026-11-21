---
esm: "2.3.0"
license: "MIT"
name: "varcoeff_face_laplacian_lev_flux_bc_tracer"
references:
  - citation: "LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids)."
  - citation: "Hong, S.-Y., Noh, Y., and Dudhia, J. (2006). A new vertical diffusion package with an explicit treatment of entrainment processes. Mon. Wea. Rev. 134, 2318-2341. (Surface-flux boundary condition of the YSU scheme.)"
    doi: "10.1175/MWR3199.1"
tags:
  - "esd:rule"
  - "family:finite_volume"
  - "grid:column_nonuniform_1d_tracer"
  - "op:varcoeff_laplacian"
  - "axes:lev,tracer"
  - "order:2"
  - "stencil_width:3"
  - "bc:prescribed_flux"
  - "spec:9.6.8"
---

```esm:description
grids/column_nonuniform_1d/rules/varcoeff_face_laplacian_lev_flux_bc applied to a field with a trailing TRACER axis, u shaped [lev, tracer]: the same compound match D( K*D(u, lev), lev ), the same stencil and the same prescribed-flux closure, column by column, with the coefficient K [lev_nodes] and the geometry shared by every tracer and the two prescribed face fluxes `kdudz_bot` / `kdudz_top` supplied per tracer ([tracer]). The body is this grid's three-region layout rules/varcoeff_face_laplacian_lev_tracer_frame with its face-flux parameters bound to those two free names, so L u[1,s] = ( G[2,s] - kdudz_bot[s] ) / dz[1] and L u[NLEV,s] = ( kdudz_top[s] - G[NLEV,s] ) / dz[NLEV], every interior flux G[n,s] being the column grid's two-point law (grids/column_nonuniform_1d/stencils/interface_flux_K_two_point) applied to column s. Column s of the result is exactly what the 1-D rule gives for u[:, s], so the telescoping identity sum_k dz[k] L u[k,s] = kdudz_top[s] - kdudz_bot[s] holds tracer by tracer. Match-scoped with where {u: {shape: [lev, tracer]}, K: {shape: [lev_nodes]}}, so it is disjoint from the 1-D rule (where u: [lev]) and a model may import both: each fires only on its own field shape (tests/conformance/ast/varcoeff_face_laplacian_lev_flux_bc_tracer imports both and lowers one field of each shape). The tracer axis comes from this grid (grids/column_nonuniform_1d_tracer/grid.typ); a consumer whose axis has another name imports this rule with `rename: {tracer: <name>}` and binds NTRACER. Motivation: one mount of a column tracer-mixing component for every chemical species of a WRF-Chem mechanism (WRF-Chem's vertmx is called once per species with the same exch_h). Second-order in L2 on a smoothly stretched column, as the 1-D rule; MMS convergence with a variable K and unequal per-tracer face fluxes is pinned by problems/heat_column_varcoeff_flux_bc_tracer.
```

```esm:expression_template_imports
- ref: "./varcoeff_face_laplacian_lev_tracer_frame.typ"
```

:::rule{params='["u", "K"]' priority=10 where='{K: {shape: ["lev_nodes"]}, u: {shape: ["lev", "tracer"]}}'}
$$
\text{rule } \operatorname{varcoeff\_face\_laplacian\_lev\_flux\_bc\_tracer}: \; \frac{\partial}{\partial\mathit{lev}} \left(K \frac{\partial u}{\partial\mathit{lev}}\right) \Rightarrow \operatorname{varcoeff\_face\_laplacian\_lev\_tracer\_frame}(G_{\mathrm{bot}} = \mathit{kdudz\_bot}, G_{\mathrm{top}} = \mathit{kdudz\_top}, K = K, u = u)
$$
:::
