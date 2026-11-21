#import "@preview/esm:0.1.0": *
#show: esm-document.with(esm: "2.3.0", license: "MIT", name: "interface_flux_K_lev_tracer", tags: (
  "esd:stencil",
  "family:finite_volume",
  "grid:column_nonuniform_1d_tracer",
  "op:interface_flux",
  "order:2",
  "stencil_width:2",
))

#description[The interface flux of grids/column\_nonuniform\_1d/stencils/interface\_flux\_K\_lev for ONE member s of a trailing tracer axis: G\[n,s\] = K\[n\] (u\[n,s\]-u\[n-1,s\])/(zc\[n\]-zc\[n-1\]). The arithmetic is the column grid's one two-point flux law, grids/column\_nonuniform\_1d/stencils/interface\_flux\_K\_two\_point (imported across grids; it reads no grid names), applied with K\_face = K\[n\], u\_lo = u\[n-1,s\], u\_hi = u\[n,s\], z\_lo = zc\[n-1\], z\_hi = zc\[n\]; this template only fixes which values an interior interface of tracer s reads. The coefficient K and the geometry are shared by every tracer; only u carries the second index. Column s is exactly interface\_flux\_K\_lev applied to u\[:,s\]. The interface index n and the tracer index s are template parameters bound at each call site (the enclosing aggregate's loop symbols k or k+1 and s in stencils/varcoeff\_face\_laplacian\_lev\_interior\_tracer, the literal 2 and the metaparameter NLEV in the two boundary layers of rules/varcoeff\_face\_laplacian\_lev\_tracer\_frame). \`zc\` is the consumer's free-name layer-centre observed (grids/column\_nonuniform\_1d/grid.esm).]

#esm-block("expression_template_imports", (
  (
    ref: "../../column_nonuniform_1d/stencils/interface_flux_K_two_point.esm",
  ),
))

$ op("interface_flux_K_lev_tracer")(u, K, n, s) := op("interface_flux_K_two_point")(K_(upright("face")) = K_n, u_(upright("hi")) = u_(n, s), u_(upright("lo")) = u_(n - 1, s), z_(upright("hi")) = italic("zc")_n, z_(upright("lo")) = italic("zc")_(n - 1)) $
