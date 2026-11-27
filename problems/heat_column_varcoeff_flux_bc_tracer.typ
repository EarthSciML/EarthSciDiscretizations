#import "@preview/esm:0.1.0": *
#show: esm-document.with(esm: "2.4.0", license: "MIT", name: "heat_column_varcoeff_flux_bc_tracer", references: (
  (
    citation: "Roache, P. J. (2002). Code verification by the method of manufactured solutions. J. Fluids Eng. 124(1), 4-10.",
    doi: "10.1115/1.1436090",
  ),
), tags: (
  "esd:problem",
  "grid:column_nonuniform_1d_tracer",
  "equation:heat",
  "mms",
))

#description[MMS driver for grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer: the 1-D heat equation u_t = d/dz(K du/dz) solved for TWO tracers at once, u shaped [lev, tracer] with NTRACER = 2 bound at the rule's import edge, on the same smoothly stretched unit column as problems/heat_column_varcoeff_flux_bc (the grid's shared MMS column mms_column_ze, quadratic stretch z(xi) = 0.6 xi + 0.4 xi^2) with the same shared interface coefficient K = alpha. Tracer s carries its OWN manufactured solution u_s(z,t) = a_s z + b_s exp(-alpha pi^2 t) cos(pi z), a = (1, -2), b = (1, 0.5) (spelled a_s = 4 - 3s, b_s = 1.5 - s/2), and therefore its OWN prescribed fluxes kdudz_bot[s] = kdudz_top[s] = alpha a_s, observeds tied to alpha. Tracer 1 is exactly problems/heat_column_varcoeff_flux_bc; tracer 2 has the opposite-signed, doubled ramp and half the cosine amplitude. Because the operator is linear, column by column, and annihilates every linear ramp (each discrete interface gradient of a_s z is a_s, cancelled at the faces by the prescribed alpha a_s), tracer 1's error field is exactly the 1-D problem's and tracer 2's is exactly half of it, so the (relative, esm-spec 6.6.5) L2 error over the [lev, tracer] field falls at the 1-D problem's second order. A per-tracer mix-up (reading another tracer's u, or another tracer's face flux) breaks the ramp cancellation and gives an O(1) error. The per-tracer column totals sum_k u[k,s] dz[k] (the grid's tracer_column_total template) are conserved exactly by the discrete operator because each tracer's two face fluxes are equal; the test pins them at t = 0 and t = 0.1 to the same values, tracer 1's as the field max and tracer 2's as the field min (an element assertion on a shaped observed, total[s], is not evaluated by every binding). The ic applies the closed per-layer centre template mms_column_zc (grids/column_nonuniform_1d/mesh/mms_column_centre.esm) inside its aggregate over [lev, tracer]; the 6.6.5 test references spell the layer centre 0.5 (z((i-1)/NLEV) + z(i/NLEV)) inline from literals and NLEV, because the inline-test oracle does not resolve template references. This problem's two face fluxes are equal for each tracer and K is constant, which is what makes the conservation and ramp-cancellation checks exact; problems/heat_column_varcoeff_forced_tracer covers a variable K with unequal, time-dependent per-tracer face fluxes.]

#model("Heat")

#var($K$, description: "Interface diffusivity, the constant alpha on every interface, shared by both tracers.", shape: (
  "lev_nodes",
), units: "m^2/s")

#param($alpha$, default: 0.01, description: "Diffusivity shared by both tracers; the constant interface coefficient and the decay-rate scale.", units: "m^2/s")

#var($italic("dz")$, description: "Layer thicknesses, grid-derived by column_cell_thickness.", shape: (
  "lev",
), units: "m")

#var($italic("kdudz_bot")$, description: "Prescribed K du/dz on the surface interface, per tracer: alpha a_s, the manufactured wall gradient a_s times alpha.", shape: (
  "tracer",
), units: "K*m/s")

#var($italic("kdudz_top")$, description: "Prescribed K du/dz on the top interface, per tracer: alpha a_s.", shape: (
  "tracer",
), units: "K*m/s")

#var($italic("total")$, description: "Per-tracer column total sum_k u[k,s] dz[k], the grid's tracer_column_total.", shape: (
  "tracer",
), units: "K*m")

#var($u$, description: "Two tracers over the lev (layer) axis.", shape: (
  "lev",
  "tracer",
), units: "K")

#var($italic("zc")$, description: "Layer centres, grid-derived by column_cell_center.", shape: (
  "lev",
), units: "m")

#var($italic("ze")$, description: "Interface coordinates of the quadratic-stretch column, the shared mms_column_ze.", shape: (
  "lev_nodes",
), units: "m")

$ u(t_0) = lr([sum ((4 - 3 s) op("mms_column_zc")(i = i) + 1 (1.5 - 0.5 s) cos(3.141592653589793 op("mms_column_zc")(i = i)))])_(i in italic("lev"), s in italic("tracer")) $

$ (dif u)/(dif t) = (partial)/(partial italic("lev")) (K (partial u)/(partial italic("lev"))) $

$ K = lr([sum alpha])_(italic("kn") in italic("lev_nodes")) $

$ italic("kdudz_bot") = lr([sum alpha (4 thin upright("K") thin upright("m")^(-1) - 3 s)])_(s in italic("tracer")) $

$ italic("kdudz_top") = lr([sum alpha (4 thin upright("K") thin upright("m")^(-1) - 3 s)])_(s in italic("tracer")) $

$ italic("total") = op("tracer_column_total")(u = u) $

$ italic("ze") = op("mms_column_ze") $

$ italic("zc") = op("column_cell_center")(italic("ze") = italic("ze")) $

$ italic("dz") = op("column_cell_thickness")(italic("ze") = italic("ze")) $

#esm-block("operator_imports", (
  (
    bindings: (
      NTRACER: 2,
    ),
    ref: "../grids/column_nonuniform_1d_tracer/rules/varcoeff_face_laplacian_lev_flux_bc_tracer.md",
  ),
  (
    bindings: (
      NTRACER: 2,
    ),
    ref: "../grids/column_nonuniform_1d_tracer/grid.typ",
  ),
  (
    ref: "../grids/column_nonuniform_1d/mesh/mms_column.esm",
  ),
  (
    ref: "../grids/column_nonuniform_1d/mesh/mms_column_centre.esm",
  ),
))

#esm-block("tests", (
  (
    assertions: (
      (
        expected: 0,
        reduce: "L2_error",
        reference: $lr([sum ((4 - 3 s) (0.5 ((0.6 (i - 1)/italic("NLEV") + 0.4 (i - 1)/italic("NLEV") (i - 1)/italic("NLEV")) + (0.6 i/italic("NLEV") + 0.4 i/italic("NLEV") i/italic("NLEV")))) + 1 (1.5 - 0.5 s) cos(3.141592653589793 (0.5 ((0.6 (i - 1)/italic("NLEV") + 0.4 (i - 1)/italic("NLEV") (i - 1)/italic("NLEV")) + (0.6 i/italic("NLEV") + 0.4 i/italic("NLEV") i/italic("NLEV"))))))])_(i in italic("lev"), s in italic("tracer"))$,
        time: 0,
        tolerance: (
          abs: 1.0e-9,
        ),
        variable: "u",
      ),
      (
        expected: 0,
        reduce: "L2_error",
        reference: $lr([sum ((4 - 3 s) (0.5 ((0.6 (i - 1)/italic("NLEV") + 0.4 (i - 1)/italic("NLEV") (i - 1)/italic("NLEV")) + (0.6 i/italic("NLEV") + 0.4 i/italic("NLEV") i/italic("NLEV")))) + 0.9901789403074717 (1.5 - 0.5 s) cos(3.141592653589793 (0.5 ((0.6 (i - 1)/italic("NLEV") + 0.4 (i - 1)/italic("NLEV") (i - 1)/italic("NLEV")) + (0.6 i/italic("NLEV") + 0.4 i/italic("NLEV") i/italic("NLEV"))))))])_(i in italic("lev"), s in italic("tracer"))$,
        time: 0.1,
        tolerance: (
          abs: 1.0e-5,
        ),
        variable: "u",
      ),
      (
        expected: 0.49996745019280364,
        reduce: "max",
        time: 0,
        tolerance: (
          abs: 1.0e-12,
          rel: 0,
        ),
        variable: "total",
      ),
      (
        expected: -1.000016274903598,
        reduce: "min",
        time: 0,
        tolerance: (
          abs: 1.0e-12,
          rel: 0,
        ),
        variable: "total",
      ),
      (
        expected: 0.49996745019280364,
        reduce: "max",
        time: 0.1,
        tolerance: (
          abs: 1.0e-9,
          rel: 0,
        ),
        variable: "total",
      ),
      (
        expected: -1.000016274903598,
        reduce: "min",
        time: 0.1,
        tolerance: (
          abs: 1.0e-9,
          rel: 0,
        ),
        variable: "total",
      ),
    ),
    description: "Two manufactured solutions u_s(z,t) = a_s z + b_s exp(-alpha pi^2 t) cos(pi z), a = (1, -2), b = (1, 0.5), alpha = 0.01, on the quadratic-stretch column with per-tracer prescribed K du/dz = alpha a_s at both interfaces. The t=0 assertion pins the ic wiring; the t=0.1 assertion bounds the relative L2 error over the [lev, tracer] field at the default NLEV = 64 (3.08e-6 in the Julia reference). The four `total` assertions (max = tracer 1, min = tracer 2) pin the per-tracer column totals to the same value at t=0 and t=0.1: the discrete operator conserves each tracer's total exactly because its two face fluxes are equal; the values are sum_k (a_s zc_k + b_s cos(pi zc_k)) dz_k on the NLEV=64 column (0.49996745019280364, -1.0000162749035979). Per-resolution errors are pinned in tests/conformance/convergence/heat_column_varcoeff_flux_bc_tracer/golden/errors.json.",
    id: "mms_column_flux_bc_tracer_decay",
    time_span: (
      end: 0.1,
      start: 0,
    ),
  ),
))
