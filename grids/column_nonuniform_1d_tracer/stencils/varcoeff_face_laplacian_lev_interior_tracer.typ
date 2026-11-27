#import "@preview/esm:0.1.0": *
#show: esm-document.with(esm: "2.4.0", license: "MIT", name: "varcoeff_face_laplacian_lev_interior_tracer", tags: (
  "esd:stencil",
  "family:finite_volume",
  "grid:column_nonuniform_1d_tracer",
  "op:varcoeff_laplacian",
  "order:2",
  "stencil_width:3",
))

#description[Interior [2, NLEV-1] of the variable-coefficient vertical diffusion for a [lev, tracer] field: (G[k+1,s] - G[k,s]) / dz[k] for every tracer s, G being stencils/interface_flux_K_lev_tracer. Column s is grids/column_nonuniform_1d/stencils/varcoeff_face_laplacian_lev_interior applied to u[:, s].]

#esm-block("operator_imports", (
  (
    ref: "./interface_flux_K_lev_tracer.typ",
  ),
  (
    ref: "../grid.typ",
  ),
))

$ op("varcoeff_face_laplacian_lev_interior_tracer")(u, K) := lr([sum (op("interface_flux_K_lev_tracer")(K = K, n = k + 1, s = s, u = u) - op("interface_flux_K_lev_tracer")(K = K, n = k, s = s, u = u))/italic("dz")_k])^({u, K})_(2 <= k <= italic("NLEV") - 1, s in upright("tracer")) $
