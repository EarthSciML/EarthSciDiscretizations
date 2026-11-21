#import "@preview/esm:0.1.0": *
#show: esm-document.with(esm: "2.3.0", license: "MIT", name: "varcoeff_face_laplacian_lev_flux_bc_tracer_ast_fixture", tags: (
  "conformance",
  "ast",
  "grid:column_nonuniform_1d_tracer",
  "rule:varcoeff_face_laplacian_lev_flux_bc_tracer",
))

#description[AST-expansion conformance fixture for grids/column\_nonuniform\_1d\_tracer/rules/varcoeff\_face\_laplacian\_lev\_flux\_bc\_tracer, and for its disjointness from the 1-D rule grids/column\_nonuniform\_1d/rules/varcoeff\_face\_laplacian\_lev\_flux\_bc: one model imports BOTH rules (the 1-D rule FIRST, so it would win the equal-priority tie on any node both matched) and writes the same compound D(K\*D(., wrt: lev), wrt: lev) for two states sharing the interface coefficient K (shape \[lev\_nodes\]): u shaped \[lev, tracer\] with the per-tracer prescribed fluxes kdudz\_bot / kdudz\_top as \[tracer\] parameters, and w shaped \[lev\]. NLEV=8 and NTRACER=3 are bound at the import edges. The 1-D rule's scalar free names kdudz\_bot / kdudz\_top are rebound at its import edge to w\_kdudz\_bot / w\_kdudz\_top (esm-spec 9.7.7 rebind), which is how a consumer mixing both field shapes keeps the two rules' free names apart. Each compound must fire exactly one rule through the fixpoint (priority 10): du/dt only the tracer rule (where u: \[lev, tracer\]; the 1-D rule's where u: \[lev\] excludes it), dw/dt only the 1-D rule (the tracer rule's where excludes it). The u equation lowers to the tracer rule's body, a reference to the tracer grid's three-region layout rules/varcoeff\_face\_laplacian\_lev\_tracer\_frame over \[lev, tracer\] (interior \[2,7\] x \[1,3\] holding the per-tracer interior stencil, layers 1 and 8 the per-tracer boundary layers as faqs over the tracer axis with an explicit singleton lev output index) with its face-flux parameters bound to kdudz\_bot / kdudz\_top; every interior flux is interface\_flux\_K\_lev\_tracer, an application of the column grid's two-point law interface\_flux\_K\_two\_point. The w equation lowers to the 1-D rule's three-region makearray reading w\_kdudz\_bot / w\_kdudz\_top. Geometry follows the column grid's edge-array contract exactly as in the 1-D fixture: ze from mesh/mms\_column.esm (NLEV=8), zc and dz by column\_cell\_center / column\_cell\_thickness from grids/column\_nonuniform\_1d/grid.esm, shared by both states. The golden pins the byte-identical post-lowering AST across bindings.]

#model("VarcoeffFaceLaplacianLevTracerAst")

#var($K$, description: "Interface diffusivity 1 + ze.", shape: (
  "lev_nodes",
), units: "m^2/s")

#var($italic("dz")$, description: "Layer thicknesses, grid-derived by applying column_cell_thickness to ze.", shape: (
  "lev",
), units: "m")

#param($italic("kdudz_bot")$, default: 0, description: "Per-tracer prescribed K du/dz on the surface interface (the rule's free name).", shape: (
  "tracer",
), units: "K*m/s")

#param($italic("kdudz_top")$, default: 0, description: "Per-tracer prescribed K du/dz on the top interface (the rule's free name).", shape: (
  "tracer",
), units: "K*m/s")

#var($u$, description: "Three tracers over the lev index set.", shape: (
  "lev",
  "tracer",
), units: "K")

#var($w$, description: "A single [lev] field beside the tracers, lowered by the 1-D rule.", shape: (
  "lev",
), units: "K")

#param($w_(upright("kdudz_bot"))$, default: 0, description: "Prescribed K dw/dz on the surface interface: the 1-D rule's free name kdudz_bot, rebound at its import edge so it does not collide with the tracer rule's [tracer] kdudz_bot.", units: "K*m/s")

#param($w_(upright("kdudz_top"))$, default: 0, description: "Prescribed K dw/dz on the top interface: the 1-D rule's free name kdudz_top, rebound at its import edge.", units: "K*m/s")

#var($italic("zc")$, description: "Layer centres, grid-derived by applying column_cell_center to ze.", shape: (
  "lev",
), units: "m")

#var($italic("ze")$, description: "Interface coordinates of the quadratic-stretch column at NLEV=8, ze[kn] = 0.6 ((kn-1)/8) + 0.4 ((kn-1)/8)^2 for kn=1..9: the grid's shared MMS test column mms_column_ze (mesh/mms_column.esm), bound at NLEV=8 at its import edge.", shape: (
  "lev_nodes",
), units: "m")

$ u(t_0) = 0 $

$ w(t_0) = 0 $

$ (dif w)/(dif t) = (partial)/(partial italic("lev")) (K (partial w)/(partial italic("lev"))) $

$ (dif u)/(dif t) = (partial)/(partial italic("lev")) (K (partial u)/(partial italic("lev"))) $

$ K = 1 + italic("ze") $

$ italic("ze") = op("mms_column_ze") $

$ italic("zc") = op("column_cell_center")(italic("ze") = italic("ze")) $

$ italic("dz") = op("column_cell_thickness")(italic("ze") = italic("ze")) $

#esm-block("expression_template_imports", (
  (
    bindings: (
      NLEV: 8,
    ),
    rebind: (
      kdudz_bot: "w_kdudz_bot",
      kdudz_top: "w_kdudz_top",
    ),
    ref: "../../../../grids/column_nonuniform_1d/rules/varcoeff_face_laplacian_lev_flux_bc.esm",
  ),
  (
    bindings: (
      NLEV: 8,
      NTRACER: 3,
    ),
    ref: "../../../../grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer.md",
  ),
  (
    bindings: (
      NLEV: 8,
    ),
    ref: "../../../../grids/column_nonuniform_1d/mesh/mms_column.esm",
  ),
  (
    bindings: (
      NLEV: 8,
    ),
    only: (
      "column_cell_center",
      "column_cell_thickness",
    ),
    ref: "../../../../grids/column_nonuniform_1d/grid.esm",
  ),
))
