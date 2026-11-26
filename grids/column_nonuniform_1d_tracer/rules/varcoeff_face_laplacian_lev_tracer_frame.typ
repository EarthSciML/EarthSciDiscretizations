#import "@preview/esm:0.1.0": *
#show: esm-document.with(esm: "2.4.0", license: "MIT", name: "varcoeff_face_laplacian_lev_tracer_frame", references: (
  (
    citation: "LeVeque, R. J. (2007). Finite Difference Methods for Ordinary and Partial Differential Equations. SIAM. Section 2.13 (nonuniform grids).",
  ),
), tags: (
  "esd:rule",
  "family:finite_volume",
  "grid:column_nonuniform_1d_tracer",
  "op:varcoeff_laplacian",
  "axes:lev,tracer",
  "order:2",
  "stencil_width:3",
  "bc:supplied_face_flux",
  "spec:9.6.8",
))

#description("The three-region layout of the conservative variable-coefficient vertical diffusion d/dz(K du/dz) for a field with a trailing TRACER axis, u shaped [lev, tracer], with the coefficient K [lev_nodes] on the interfaces and the geometry shared by every tracer, whatever closes the two domain faces: a makearray tiling [lev, tracer] whose interior [2, NLEV-1] x [1, NTRACER] is stencils/varcoeff_face_laplacian_lev_interior_tracer, (G[k+1,s] - G[k,s]) / dz[k], and whose two boundary layers are L u[1,s] = ( G[2,s] - G_bot[s] ) / dz[1] and L u[NLEV,s] = ( G_top[s] - G[NLEV,s] ) / dz[NLEV], every interior flux G[n,s] being stencils/interface_flux_K_lev_tracer (the column grid's two-point law grids/column_nonuniform_1d/stencils/interface_flux_K_two_point applied to column s). Column s is exactly the 1-D layout of grids/column_nonuniform_1d applied to u[:, s]. The two domain-face fluxes are PARAMETERS, [tracer]-shaped arrays read element by element as G_bot[s] (the diffusive flux K du/dz of tracer s through the surface face z = ze[1]) and G_top[s] (through the top face z = ze[NLEV+1]), positive in the +lev (upward) direction and bound by the applying match rule to whatever its boundary condition makes them -- the per-tracer free names kdudz_bot / kdudz_top in rules/varcoeff_face_laplacian_lev_flux_bc_tracer. So the interior stencil, the region layout and the boundary-layer flux differences have one spelling on this grid, and the telescoping identity sum_k dz[k] L u[k,s] = G_top[s] - G_bot[s] holds tracer by tracer for every closure. A match-less rule template (it never fires on its own); the bc: tag names what it takes, not a closure it imposes. The two boundary layers are faqs with output_idx [1, s] -- an explicit singleton lev dimension -- because the be1fafe6d Rust binding fills the WHOLE array with NaN when a rank-1 [tracer] value is assigned to a [1, NTRACER] region, although esm-spec 4.3.2 says singleton dimensions are excluded from the shape match. `zc` and `dz` are the consumer's free-name geometry observeds (grids/column_nonuniform_1d/grid.esm). Minimum extent NLEV = 2, as on the column grid.")

#esm-block("operator_imports", (
  (
    ref: "../stencils/varcoeff_face_laplacian_lev_interior_tracer.typ",
  ),
  (
    ref: "../stencils/interface_flux_K_lev_tracer.typ",
  ),
  (
    ref: "../grid.typ",
  ),
))

$ op("varcoeff_face_laplacian_lev_tracer_frame")(u, K, G_(upright("bot")), G_(upright("top"))) := mat(delim: "[", augment: #1, op("varcoeff_face_laplacian_lev_interior_tracer")(K = K, u = u), [2, italic("NLEV") - 1] times [1, italic("NTRACER")]; lr([sum (op("interface_flux_K_lev_tracer")(K = K, n = 2, s = s, u = u) - G_(upright("bot"))[s])/italic("dz")_1])^({u, K})_(1, s in italic("tracer")), [1, 1] times [1, italic("NTRACER")]; lr([sum (G_(upright("top"))[s] - op("interface_flux_K_lev_tracer")(K = K, n = italic("NLEV"), s = s, u = u))/italic("dz")_(italic("NLEV"))])^({u, K})_(1, s in italic("tracer")), [italic("NLEV"), italic("NLEV")] times [1, italic("NTRACER")]) $
