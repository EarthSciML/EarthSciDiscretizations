# EarthSciDiscretizations — Rule Catalog (Phase 0 Inventory)

Authoritative manifest of discretization schemes and grids targeted as `.esm`
content under `discretizations/`. Each row corresponds to one scheme, grid,
staggering convention, or coordinate transform that we want represented as a
validated rule file (or set of files) governed by the EarthSciSerialization
§7 rule schema.

This is a Phase 0 inventory: the catalog drives all subsequent per-scheme
content beads. Correctness of pre-Gas-Town code in `src/` is **not**
guaranteed — the `audit_status` column flags items that need review before
being promoted to authoritative content.

## Header metadata

| Field | Value |
|---|---|
| Inventory date (UTC) | 2026-04-19T22:52Z |
| EarthSciDiscretizations source SHA | `398b5ca88f3abb901ac334d512ca530ed1f4b8b3` (branch tip at inventory time) |
| MOL PR #531 head SHA | `35cc9143dc553ac7d3619738bd77b250c1ed162f` |
| MOL PR #531 base SHA | `ec83208e6f0819f0ab7ff00b77c775b9d7cc18aa` |
| CAM5 SciDoc retrieval | 2026-04-19 (https://ncar.github.io/CAM_SciDoc/doc/build/html/cam5_scientific_guide/dynamics.html) |
| FV3 reference | Harris et al. (2021), GFDL FV3 Technical Memorandum |
| MPAS reference | Skamarock et al. (2012), MWR 140; Thuburn et al. (2009), JCP 228; Ringler et al. (2010), JCP 229 |

## Priority bucket counts

Counts include intentional duplicates between Section F (audit lens) and the
source-specific sections — the audit row for an existing item is its own
deliverable (a review bead) distinct from the content bead implied by the
source section.

| Bucket | Count | Definition |
|---|---:|---|
| P0 — gating | 7 | needed immediately for pilot ESM models |
| P1 — core | 66 | widely-used; first 12 months of content work |
| P2 — common but deferrable | 40 | used in production models, not v1 critical |
| P3 — specialized | 7 | narrow or research-grade |
| Deferred | 8 | out-of-scope for v1 |
| **Total rows** | **128** | |

Per-source counts (a row may map to multiple sources; counts are per section):

| Source | Section | Rows |
|---|---|---:|
| MOL PR #531 | A | 15 |
| CAM5 SciDoc §4.1 (FV core) | B | 14 |
| CAM5 SciDoc §4.2 (HOMME / SEM) | C | 8 |
| FV3 Technical Memorandum | D | 13 |
| MPAS papers + technical notes | E | 10 |
| Pre-Gas-Town `src/` audit | F | 13 |
| Deep-dive (other earth-science) | G.1–G.7 | 49 |
| Deferred (out-of-scope v1) | G.8 | 8 listed |

---

## How to read this catalog

Columns:

| Column | Meaning |
|---|---|
| `name` | Canonical short name used as filename stem (e.g. `centered_2nd_uniform`). |
| `family` | One of `finite_difference` / `finite_volume` / `spectral` / `semi_lagrangian` / `dg` / `grid` / `limiter` / `bc` / `time` / `solver`. |
| `kind` | One of `scheme` (stencil/rule), `grid` (topology), `staggering`, `coord_transform`, `bc`, `limiter`, `time_integrator`, `solver`. |
| `rfc_section` | ESS RFC §7 etc.; `extension` if beyond current RFC. |
| `source_ref` | URL or repo+path:line for reference definition. |
| `priority` | P0 / P1 / P2 / P3 / Deferred. |
| `earth_sci_relevance` | Model classes that use this. |
| `schema_gaps` | Schema features needed beyond ESS main (cite bead IDs where known). |
| `complexity` | S / M / L / XL. |
| `audit_status` | `not_in_repo` / `looks_ok` / `needs_review` / `known_wrong` for pre-Gas-Town code. |
| `target_path` | Intended `.esm` path. |
| `depends_on` | Other rows this builds on. |
| `notes` | Caveats, review flags, justification (mandatory for P3+). |

---

## Section A — MethodOfLines.jl PR #531

Source: https://github.com/SciML/MethodOfLines.jl/pull/531 (head SHA above).
Each row corresponds to a named scheme, rule, or strategy introduced or
modified by this PR. Files are referenced relative to the MOL repo at PR head.

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `mol_array_discretization_strategy` | `finite_difference` | `scheme` | extension | `src/array_discretization.jl`, `src/discretization/array_fd/*` | P1 | atmos.chem, ocean.bgc, generic-PDE | needs ArrayOp template emission flag in §7 | L | not_in_repo | `discretizations/finite_difference/mol_array_discretization.esm` | `mol_scalarized_discretization_strategy` | New strategy in PR #531 — emits a single ArrayOp expression for the centered interior instead of N scalar equations. Fundamental for performance; all subsequent FD rules below should compose with this. |
| `mol_scalarized_discretization_strategy` | `finite_difference` | `scheme` | §7 | `src/MOL_discretization.jl`, `src/interface/disc_strategy_types.jl` | P2 | reference / fallback path | none | M | not_in_repo | `discretizations/finite_difference/mol_scalarized_discretization.esm` | — | Existing per-grid-point scalar emission. Kept as the validation oracle for ArrayDiscretization and as fallback for unsupported patterns. |
| `centered_2nd_uniform` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/centered_difference/centered_difference.jl`; mirrored as `discretizations/finite_difference/centered_2nd_uniform.json` | P0 | atmos.dyn, atmos.chem.transport, ocean.dyn, generic-PDE | none | S | looks_ok | `discretizations/finite_difference/centered_2nd_uniform.esm` | — | Already present as JSON in repo; promote to `.esm` once schema lands. Foundation rule. |
| `centered_arbitrary_order_uniform` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/centered_difference/centered_difference.jl`, `fornberg_calculate_weights.jl` | P1 | atmos.dyn (high-order), spectral-equivalent FD | needs `order` parameter in stencil rule | M | not_in_repo | `discretizations/finite_difference/centered_arbitrary_order_uniform.esm` | `fornberg_weights` | Even-order centered differences via Fornberg recursion (4th, 6th, 8th, …). |
| `centered_arbitrary_order_nonuniform` | `finite_difference` | `scheme` | §7 | `src/discretization/array_fd/precompute.jl`, `fornberg_calculate_weights.jl` | P1 | atmos.dyn (stretched grids), ocean.dyn | needs per-point precomputed weight tensor | L | not_in_repo | `discretizations/finite_difference/centered_arbitrary_order_nonuniform.esm` | `fornberg_weights` | Non-uniform variant — grid-spacing-dependent weights cached on the grid. |
| `upwind_arbitrary_order_uniform` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/upwind_difference/upwind_difference.jl`, `upwind_diff_weights.jl` | P1 | atmos.chem.transport, atmos.dyn (advection) | needs sign-dependent stencil branching | M | partial | `discretizations/finite_difference/upwind_arbitrary_order_uniform.esm` | `centered_arbitrary_order_uniform` | Existing repo has only `upwind_1st.json`. PR #531 supports arbitrary odd order. |
| `upwind_arbitrary_order_nonuniform` | `finite_difference` | `scheme` | §7 | `src/discretization/array_fd/rules_upwind.jl` | P2 | atmos.chem.transport on stretched grids | sign branching + per-point weights | L | not_in_repo | `discretizations/finite_difference/upwind_arbitrary_order_nonuniform.esm` | `upwind_arbitrary_order_uniform`, `centered_arbitrary_order_nonuniform` | Non-uniform upwind. |
| `mixed_2nd_cross_derivative` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/2nd_order_mixed_deriv/2nd_order_mixed_deriv.jl`, `src/discretization/array_fd/rules_mixed.jl` | P1 | atmos.dyn (cross-metric Laplacian), ocean.dyn | needs 2D outer-product stencil rule | M | not_in_repo | `discretizations/finite_difference/mixed_2nd_cross_derivative.esm` | `centered_2nd_uniform` | Dxy(u) — needed for any cross-metric term, including our cubed-sphere Laplacian. |
| `nonlinear_laplacian_div_form` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/nonlinear_laplacian/nonlinear_laplacian.jl`, `src/discretization/array_fd/rules_nonlinlap.jl` | P1 | atmos.chem (nonlinear diffusion), porous-flow, hydrology | needs coefficient-evaluation hook in stencil | L | not_in_repo | `discretizations/finite_difference/nonlinear_laplacian_div_form.esm` | `centered_2nd_uniform` | Dx(f(u)·Dx(u)) — nonlinear diffusion. PR provides both uniform + non-uniform. |
| `spherical_laplacian_radial` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/spherical_laplacian/spherical_laplacian.jl`, `src/discretization/array_fd/rules_spherical.jl` | P2 | atmos (1D radial), planetary, stellar | needs radial coordinate metric | M | not_in_repo | `discretizations/finite_difference/spherical_laplacian_radial.esm` | `nonlinear_laplacian_div_form` | r⁻²·Dr(r²·Dr(u)). 1D radial Laplacian with metric correction. |
| `weno5_jiang_shu_uniform` | `finite_volume` | `scheme` | §7 | `src/discretization/schemes/WENO/WENO.jl`, `src/discretization/array_fd/rules_weno.jl` | P1 | atmos.chem.transport (sharp fronts), ocean.dyn (eddies) | needs WENO smoothness-indicator rule type | XL | not_in_repo | `discretizations/finite_volume/weno5_jiang_shu_uniform.esm` | `centered_2nd_uniform` | Jiang-Shu 5th-order WENO; PR supports uniform grid only. Critical for chemistry transport with sharp gradients. |
| `half_offset_centered_difference` | `finite_difference` | `scheme` | §7 | `src/discretization/schemes/half_offset_centred_difference.jl`, `half_offset_weights.jl` | P2 | atmos.dyn (wave equation, staggered velocities) | needs half-integer index selector | M | not_in_repo | `discretizations/finite_difference/half_offset_centered.esm` | `centered_arbitrary_order_uniform`, `staggered_grid_arakawa_c` | Centered difference between staggered points. PR adds odd-order support on uniform grid. |
| `function_scheme_user_supplied` | `finite_difference` | `scheme` | extension | `src/discretization/schemes/function_scheme/function_scheme.jl` | P3 | model-specific custom stencils | needs callable-rule type in §7 | M | not_in_repo | `discretizations/finite_difference/function_scheme.esm` | — | Justification (P3): user supplies a Julia callable. Powerful escape hatch but resists declarative validation; needed only for research models that don't fit declarative rules. |
| `integral_expansion` | `finite_volume` | `scheme` | extension | `src/discretization/schemes/integral_expansion/integral_expansion.jl` | P3 | atmos.rad (column-integrated quantities), ocean.bgc | needs integral-operator type in §7 | M | not_in_repo | `discretizations/finite_volume/integral_expansion.esm` | — | Justification (P3): turns ∫ operators into trapezoidal/Simpson sums. Needed for column-integrated diagnostics; not gating for transport/dynamics. |
| `fornberg_weights` | `finite_difference` | `scheme` | extension | `src/discretization/schemes/fornberg_calculate_weights.jl`, `extrapolation_weights.jl` | P1 | utility used by all FD rules | none | S | not_in_repo | `discretizations/finite_difference/fornberg_weights.esm` | — | Fornberg (1988) recursion for FD weights at arbitrary order on arbitrary point sets. Underlies all FD rule families. |

---

## Section B — CAM5 dynamical-core schemes (FV core, §4.1)

Source: NCAR CAM5 Scientific Guide §4.1 (Finite-Volume Dynamical Core).
Retrieval date 2026-04-19. The page is not git-versioned; `source_ref` cites
section numbers within the SciDoc.

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `cam5_fv_ffsl_advection` | `semi_lagrangian` | `scheme` | extension | CAM5 SciDoc §4.1.3 (FFSL); Lin & Rood (1996) MWR | P1 | atmos.dyn, atmos.chem.transport | needs flux-form SL rule type | XL | needs_review | `discretizations/semi_lagrangian/cam5_fv_ffsl_advection.esm` | `cam5_fv_ppm_reconstruction`, `cam5_fv_dimensional_split` | "Flux-form semi-Lagrangian" — 2D dimensionally split, conservative. |
| `cam5_fv_ppm_reconstruction` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.3 (PPM); Colella & Woodward (1984) | P0 | atmos.dyn, atmos.chem.transport, ocean.dyn | none beyond MPAS-style polynomial-reconstruction | L | needs_review | `discretizations/finite_volume/ppm_reconstruction.esm` | `centered_2nd_uniform`, `colella_woodward_limiter` | Already partially in repo (`src/operators/reconstruction.jl`). Boundary interfaces use a low-order fallback (cell-mean) — flag in audit. |
| `cam5_fv_advective_form_operators` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.3 | P2 | atmos.dyn (alternative to FFSL inner op) | none | M | not_in_repo | `discretizations/finite_volume/cam5_advective_form.esm` | `cam5_fv_ffsl_advection` | Constant cell-mean inner operator alternative. |
| `cam5_fv_2d_lagrangian_dynamics` | `finite_volume` | `scheme` | extension | CAM5 SciDoc §4.1.4 | P1 | atmos.dyn | needs Lagrangian-vertical coupling rule | XL | not_in_repo | `discretizations/finite_volume/cam5_2d_lagrangian.esm` | `hybrid_eta_vertical_coord`, `cam5_fv_ffsl_advection` | Vertically Lagrangian / horizontally Eulerian shallow-water-like core. |
| `cam5_fv_div_damping_2nd` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.5 | P1 | atmos.dyn (numerical stability) | none | S | not_in_repo | `discretizations/finite_volume/div_damping_2nd.esm` | `fv_divergence_cubed_sphere` | ν·∇²(div) damping of divergent modes. |
| `cam5_fv_div_damping_4th` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.5 | P1 | atmos.dyn (scale-selective damping) | needs ∇⁴ stencil composition rule | M | not_in_repo | `discretizations/finite_volume/div_damping_4th.esm` | `cam5_fv_div_damping_2nd` | Scale-selective fourth-order. |
| `cam5_fv_laplacian_wind_damping` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.5 | P2 | atmos.dyn (polar night jet control) | none | S | not_in_repo | `discretizations/finite_volume/laplacian_wind_damping.esm` | `cam5_fv_div_damping_2nd` | Optional ∇² of wind components. |
| `cam5_fv_mass_energy_conserving_remap` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.6; Lin (2004) MWR | P1 | atmos.dyn (Lagrangian-to-Eulerian remap) | needs conservative-remap rule type | L | needs_review | `discretizations/finite_volume/mass_energy_conserving_remap.esm` | `cam5_fv_ppm_reconstruction`, `vertical_remap_lin_2004` | Already partially in `src/operators/vertical_remap.jl`. |
| `cam5_fv_geopotential_conserving_remap` | `finite_volume` | `scheme` | §7 | CAM5 SciDoc §4.1.7 | P2 | atmos.dyn (alternate vertical remap) | none | M | not_in_repo | `discretizations/finite_volume/geopotential_conserving_remap.esm` | `cam5_fv_mass_energy_conserving_remap` | Alternative remap preserving model-lid geopotential. |
| `cam5_fv_negative_tracer_fixer` | `finite_volume` | `scheme` | extension | CAM5 SciDoc §4.1.9 | P1 | atmos.chem (positivity), ocean.bgc | needs neighbor-borrowing rule type | M | not_in_repo | `discretizations/finite_volume/negative_tracer_fixer.esm` | `cam5_fv_ffsl_advection` | Borrows mass from neighbors when tracer goes negative. Critical for chemistry. |
| `cam5_fv_global_energy_fixer` | `finite_volume` | `scheme` | extension | CAM5 SciDoc §4.1.10 | P2 | atmos.dyn (global energy budget) | needs global-reduction rule type | M | not_in_repo | `discretizations/finite_volume/global_energy_fixer.esm` | `cam5_fv_2d_lagrangian_dynamics` | Uniform DSE adjustment to close energy budget. |
| `cam5_fv_polar_filter_fft` | `finite_volume` | `scheme` | extension | CAM5 SciDoc §4.1 (polar filter) | P2 | atmos.dyn (lat-lon poles only) | needs FFT/spectral-filter rule type | L | not_in_repo | `discretizations/finite_volume/polar_filter_fft.esm` | `lat_lon_regular_grid` | FFT-based filter to stabilize short-wavelength gravity waves at high latitudes. Lat-lon-grid-specific. |
| `cam5_fv_split_explicit_time` | `time` | `time_integrator` | extension | CAM5 SciDoc §4.1 | P1 | atmos.dyn | needs time-integrator rule type | M | not_in_repo | `discretizations/time/cam5_split_explicit.esm` | — | m subcycles for dynamics, 1 large step for tracer transport. |
| `hybrid_eta_vertical_coord` | `grid` | `coord_transform` | extension | CAM5 SciDoc §4.1; Simmons & Burridge (1981) | P1 | atmos.dyn (all production GCMs) | needs vertical-coord rule type | M | not_in_repo | `discretizations/grids/hybrid_eta_vertical.esm` | — | σ-p hybrid with terrain-following Lagrangian surfaces. Universal in operational atmosphere models. |

---

## Section C — CAM5 Spectral Element Core (HOMME, §4.2)

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `homme_continuous_galerkin_sem` | `spectral` | `scheme` | extension | CAM5 SciDoc §4.2.3; Taylor et al. (1997) | P2 | atmos.dyn (CESM, E3SM), ocean.dyn (MPAS-O alt) | needs spectral-element rule type (basis, quadrature, mass matrix) | XL | not_in_repo | `discretizations/spectral/homme_cg_sem.esm` | `cubed_sphere_grid_equiangular`, `equal_angle_gnomonic_projection` | High-order CG SEM on cubed sphere. Used by CAM-SE/E3SM. |
| `equal_angle_gnomonic_projection` | `grid` | `coord_transform` | §7 | CAM5 SciDoc §4.2.3; Sadourny (1972) | P1 | atmos.dyn (cubed-sphere models), ocean.dyn | none beyond cubed-sphere panel definitions | M | not_in_repo (we use **equidistant** in `src/`, not equiangular) | `discretizations/grids/equal_angle_gnomonic.esm` | `cubed_sphere_grid_equiangular` | **Important**: our `CubedSphereGrid` uses gnomonic *equidistant*, not equiangular. Document both variants — they yield different metric tensors. |
| `homme_vector_invariant_momentum` | `spectral` | `scheme` | extension | CAM5 SciDoc §4.2.3; Lauritzen et al. (2014) | P2 | atmos.dyn (SEM cores) | needs Coriolis + KE + vorticity composition | L | not_in_repo | `discretizations/spectral/vector_invariant_momentum.esm` | `homme_continuous_galerkin_sem`, `fv3_kinetic_energy` (analogue) | Alternative to vorticity-divergence form. |
| `homme_surface_pressure_prognostic` | `spectral` | `scheme` | extension | CAM5 SciDoc §4.2.3 | P2 | atmos.dyn (mass conservation) | none | S | not_in_repo | `discretizations/spectral/surface_pressure_prognostic.esm` | `homme_continuous_galerkin_sem` | Advect ps directly (not log ps) for conservation. |
| `homme_vertically_lagrangian_tracer` | `spectral` | `scheme` | extension | CAM5 SciDoc §4.2.3 | P2 | atmos.chem.transport (in SEM cores) | needs Lagrangian-vertical hook | L | not_in_repo | `discretizations/spectral/lagrangian_vertical_tracer.esm` | `homme_continuous_galerkin_sem`, `cam5_fv_mass_energy_conserving_remap` | Same Lagrangian vertical idea as FV core, applied to SEM tracers. |
| `homme_multistage_rk_time` | `time` | `time_integrator` | extension | CAM5 SciDoc §4.2.3 | P2 | atmos.dyn (SEM) | needs RK-stage rule type | M | not_in_repo | `discretizations/time/multistage_rk_dynamics.esm` | `cam5_fv_split_explicit_time` | More stages for dynamics than tracers (gravity-wave stability). |
| `cubed_sphere_grid_equiangular` | `grid` | `grid` | §7 | CAM5 SciDoc §4.2.3; Ronchi et al. (1996) | P1 | atmos.dyn, ocean.dyn | needs cubed-sphere panel topology rule | L | partial (we have **equidistant** variant) | `discretizations/grids/cubed_sphere_equiangular.esm` | `equal_angle_gnomonic_projection` | Equiangular variant; sister to the equidistant version we already have in `src/grids/cubed_sphere.jl`. |
| `cam5_se_quadrature_glln` | `spectral` | `scheme` | extension | CAM5 SciDoc §4.2.3; standard SEM theory | P2 | atmos.dyn (SEM nodes) | needs quadrature-rule type | M | not_in_repo | `discretizations/spectral/glln_quadrature.esm` | `homme_continuous_galerkin_sem` | Gauss-Lobatto-Legendre nodes — defines SEM degrees of freedom. |

---

## Section D — FV3 grid + paired schemes

Source: Harris et al. (2021), "GFDL FV3 Technical Memorandum". The
cubed-sphere FV core used in GFS/SHiELD/UFS.

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `fv3_cubed_sphere_grid` | `grid` | `grid` | §7 | Harris et al. (2021) §3; `src/grids/cubed_sphere.jl` | P1 | atmos.dyn (GFS, UFS, SHiELD), ocean.dyn | needs cubed-sphere topology + corner-singularity flag | L | needs_review | `discretizations/grids/fv3_cubed_sphere.esm` | `equal_angle_gnomonic_projection` | Existing repo grid uses equidistant gnomonic. FV3 uses equiangular. Document both; flag mismatch. |
| `fv3_cubed_sphere_corner_singularities` | `grid` | `grid` | §7 | Harris et al. (2021) §3.4; `src/grids/panel_connectivity.jl` | P1 | atmos.dyn (cubed-sphere models) | none | M | needs_review | `discretizations/grids/cubed_sphere_corners.esm` | `fv3_cubed_sphere_grid` | 8 corners where 3 panels meet — special handling required. Existing `panel_connectivity.jl` covers normal edges; corner code needs audit. |
| `fv3_super_grid_sin_cos` | `grid` | `grid` | §7 | Harris et al. (2021) §3.5; `src/grids/super_grid.jl` | P1 | atmos.dyn (non-orthogonal flux comp) | needs sub-cell sample-point rule | L | looks_ok | `discretizations/grids/fv3_super_grid.esm` | `fv3_cubed_sphere_grid` | sin(α), cos(α) at 9 sub-positions per cell for non-orthogonal cubed sphere. |
| `fv3_d_grid_staggering` | `grid` | `staggering` | §7 | Harris et al. (2021) §3.2 Table 5.1; `src/staggering.jl` | P1 | atmos.dyn (FV3, FV core) | needs C/D-grid distinction in §7 | M | looks_ok | `discretizations/grids/fv3_d_grid_staggering.esm` | `fv3_cubed_sphere_grid` | Normalized covariant winds at face midpoints. |
| `fv3_c_grid_staggering` | `grid` | `staggering` | §7 | Harris et al. (2021) §5; `src/staggering.jl` | P1 | atmos.dyn (FV core) | needs face-normal velocity staggering | S | looks_ok | `discretizations/grids/fv3_c_grid_staggering.esm` | `fv3_cubed_sphere_grid` | Face-normal contravariant winds. |
| `fv3_cd_grid_transformation` | `finite_volume` | `scheme` | §7 | Harris et al. (2021) §5; `src/operators/wind_ops.jl` | P1 | atmos.dyn (FV3 momentum) | none | M | needs_review | `discretizations/finite_volume/fv3_cd_grid_transform.esm` | `fv3_d_grid_staggering`, `fv3_c_grid_staggering` | C↔D grid algorithm — converts D-grid covariant to C-grid contravariant for advection. |
| `fv3_ppm_two_sided_edge` | `finite_volume` | `scheme` | §7 | Harris et al. (2021) Eq. 6.5–6.6; ported to `discretizations/finite_volume/ppm_edge_cubed_sphere.json` (dsc-9yh) | P1 | atmos.dyn (FV3 advection at panel edges) | needs cross-panel ghost-cell access in stencil | M | ported | `discretizations/finite_volume/ppm_edge_cubed_sphere.json` | `cam5_fv_ppm_reconstruction`, `fv3_cubed_sphere_corner_singularities` | Two-sided extrapolation for PPM at cube edges. Ported 2026-04-27 — 4-cell linear stencil (offsets {-1,0,1,2}, coefficients [-1/4, 3/4, 3/4, -1/4]) on the gnomonic equidistant cubed sphere; FV3 eq. 6.6 monotonicity clamp documented separately as a non-linear post-processing step pending an ESS schema extension. Imperative `src/operators/ppm_edge.jl` deleted alongside the port. |
| `fv3_lin_rood_advection` | `finite_volume` | `scheme` | extension | Lin & Rood (1996) MWR; Putman & Lin (2007) JCP; `src/operators/transport_2d.jl` | P1 | atmos.chem.transport, atmos.dyn | needs dimensional-split rule type | XL | needs_review | `discretizations/finite_volume/lin_rood_advection.esm` | `cam5_fv_ppm_reconstruction`, `fv3_ppm_two_sided_edge` | Dimensionally-split 2D PPM advection on cubed sphere. Existing `src/operators/transport_2d.jl` includes a Lin-Rood path alongside Lax-Friedrichs. |
| `fv_divergence_cubed_sphere` | `finite_volume` | `scheme` | §7 | `src/operators/divergence.jl`; cf. MPAS Skamarock §7.3 | P0 | atmos.dyn, ocean.dyn | none | S | looks_ok | `discretizations/finite_volume/fv_divergence_cubed_sphere.esm` | `fv3_cubed_sphere_grid`, `fv3_c_grid_staggering` | Already in repo. Maps directly to ESS RFC §7.3 divergence example. |
| `fv_gradient_cubed_sphere` | `finite_volume` | `scheme` | §7 | `src/operators/gradient.jl` | P0 | atmos.dyn, ocean.dyn | none | S | looks_ok | `discretizations/finite_volume/fv_gradient_cubed_sphere.esm` | `fv3_cubed_sphere_grid` | Edge gradient using physical center-to-center distance. |
| `fv_laplacian_full_covariant_cubed_sphere` | `finite_volume` | `scheme` | §7 | `src/operators/laplacian.jl` | P1 | atmos.dyn (diffusion on sphere) | needs cross-metric stencil composition | L | needs_review | `discretizations/finite_volume/fv_laplacian_full_covariant.esm` | `mixed_2nd_cross_derivative`, `fv3_cubed_sphere_grid` | Existing 9-point covariant Laplacian. Cross-metric expansion is complex; needs careful audit. |
| `fv3_kinetic_energy` | `finite_volume` | `scheme` | extension | Harris et al. (2021) Eq. 6.3; `src/operators/kinetic_energy.jl` | P1 | atmos.dyn (vector-invariant momentum) | needs covariant/contravariant pairing rule | L | needs_review | `discretizations/finite_volume/fv3_kinetic_energy.esm` | `fv3_d_grid_staggering`, `fv3_super_grid_sin_cos` | Upstream-biased KE on non-orthogonal grid. |
| `fv3_vorticity_corner_and_cellmean` | `finite_volume` | `scheme` | §7 | Harris et al. (2021) §3,5,6; `src/operators/vorticity.jl` | P1 | atmos.dyn (FV3 PV flux) | needs corner-staggered output | M | looks_ok | `discretizations/finite_volume/fv3_vorticity.esm` | `fv3_d_grid_staggering`, `fv3_cubed_sphere_grid` | Two variants in code: corner-point (interpolated) and cell-mean (Stokes' theorem). |

---

## Section E — MPAS grid + paired schemes

Source: Skamarock et al. (2012), MWR 140; Thuburn et al. (2009), JCP 228;
Ringler et al. (2010), JCP 229; Skamarock & Gassmann (2011) MWR 139.

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `mpas_voronoi_mesh` | `grid` | `grid` | §7 | Skamarock et al. (2012); Du et al. (1999) SCVT | P1 | atmos.dyn (MPAS-A), ocean.dyn (MPAS-O) | needs unstructured-mesh topology in §7 (cells, edges, vertices, dual) | XL | not_in_repo | `discretizations/grids/mpas_voronoi.esm` | — | Spherical Centroidal Voronoi Tessellation. Variable-resolution capable. |
| `mpas_delaunay_dual_mesh` | `grid` | `grid` | §7 | Ringler et al. (2010) | P1 | atmos.dyn (vorticity points), ocean.dyn | needs primal/dual coupling | M | not_in_repo | `discretizations/grids/mpas_dual_mesh.esm` | `mpas_voronoi_mesh` | Triangular dual mesh — vertices of primal are cell centers of dual. |
| `mpas_c_grid_staggering` | `grid` | `staggering` | §7 | Skamarock et al. (2012); Thuburn et al. (2009) | P1 | atmos.dyn (MPAS-A), ocean.dyn | needs unstructured-C-grid staggering rule | M | not_in_repo | `discretizations/grids/mpas_c_grid_staggering.esm` | `mpas_voronoi_mesh`, `mpas_delaunay_dual_mesh` | Scalars at cell centers, normal velocities at edge midpoints, vorticity at vertices. |
| `mpas_divergence_flux_form` | `finite_volume` | `scheme` | §7.3 | Skamarock et al. (2012) Eq. 24; ESS RFC §7.3 worked example | P0 | atmos.dyn, ocean.dyn (canonical MPAS example) | none — this IS the §7 worked example | S | not_in_repo | `discretizations/finite_volume/mpas_divergence_flux_form.esm` | `mpas_c_grid_staggering` | The RFC §7.3 reference example. Must be the first MPAS rule produced. |
| `mpas_gradient_edge_difference` | `finite_volume` | `scheme` | §7 | Skamarock et al. (2012) Eq. 21 | P1 | atmos.dyn, ocean.dyn | none | S | not_in_repo | `discretizations/finite_volume/mpas_gradient.esm` | `mpas_c_grid_staggering` | Edge-centered gradient = (φ_{c2} − φ_{c1}) / d_e. |
| `mpas_vorticity_circulation` | `finite_volume` | `scheme` | §7 | Skamarock et al. (2012); Ringler et al. (2010) | P1 | atmos.dyn, ocean.dyn | needs vertex-staggered output | M | not_in_repo | `discretizations/finite_volume/mpas_vorticity.esm` | `mpas_c_grid_staggering`, `mpas_delaunay_dual_mesh` | Circulation around dual-cell boundary / dual-cell area. |
| `mpas_nonlinear_coriolis` | `finite_volume` | `scheme` | extension | Thuburn et al. (2009); Ringler et al. (2010) | P1 | atmos.dyn (MPAS-A), ocean.dyn | needs flux-reconstruction-coefficient table per edge | XL | not_in_repo | `discretizations/finite_volume/mpas_nonlinear_coriolis.esm` | `mpas_flux_reconstruction`, `mpas_vorticity_circulation` | TRSK-family nonlinear Coriolis term — energy/enstrophy conserving. |
| `mpas_flux_reconstruction` | `finite_volume` | `scheme` | extension | Skamarock & Gassmann (2011) | P1 | atmos.dyn, ocean.dyn | needs per-edge stencil-coefficient input | L | not_in_repo | `discretizations/finite_volume/mpas_flux_reconstruction.esm` | `mpas_c_grid_staggering` | Tangential velocity reconstruction from neighboring normal components (TRSK weights). |
| `mpas_laplacian_div_grad` | `finite_volume` | `scheme` | §7 | Skamarock et al. (2012) | P1 | atmos.dyn (mixing), ocean.dyn | none | M | not_in_repo | `discretizations/finite_volume/mpas_laplacian.esm` | `mpas_divergence_flux_form`, `mpas_gradient_edge_difference` | ∇·∇φ via composition of MPAS divergence and gradient. |
| `mpas_advection_2nd_3rd_order` | `finite_volume` | `scheme` | §7 | Skamarock & Gassmann (2011) | P1 | atmos.chem.transport (MPAS), ocean.dyn | needs sign-dependent edge stencil | L | not_in_repo | `discretizations/finite_volume/mpas_advection.esm` | `mpas_c_grid_staggering`, `mpas_flux_reconstruction` | 2nd/3rd-order upwind-biased advection on Voronoi cells. |

---

## Section F — Audit of pre-Gas-Town `src/` content

Items that already exist in the repo as Julia code or JSON rule files.
The `audit_status` column is the **single most important** column here:
content is *not* guaranteed correct.

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `centered_2nd_uniform` (existing JSON) | `finite_difference` | `scheme` | §7 | `discretizations/finite_difference/centered_2nd_uniform.json` | P0 | (same as Section A) | none | S | looks_ok | already exists | — | Already in repo. Validate against §7 schema once it lands. |
| `upwind_1st` (existing JSON) | `finite_difference` | `scheme` | §7 | `discretizations/finite_difference/upwind_1st.json` | P0 | atmos.chem.transport | none | S | looks_ok | already exists | — | Already in repo. Donor-cell first-order upwind. |
| `periodic_bc_wrap_x` (existing JSON) | `bc` | `bc` | §7 | `discretizations/finite_difference/periodic_bc.json` | P1 | atmos.chem.transport (channel models), ocean.dyn (zonal channels) | needs `region` and `where`-guard semantics | S | looks_ok | already exists | — | Already in repo. Periodic wrap rule. |
| `cubed_sphere_grid_equidistant` | `grid` | `grid` | §7 | `src/grids/cubed_sphere.jl` | P1 | atmos.dyn (this repo's pilot models), ocean.dyn | needs cubed-sphere topology rule | L | needs_review | `discretizations/grids/cubed_sphere_equidistant.esm` | `equal_angle_gnomonic_projection` (sister) | Existing grid uses equidistant gnomonic. Mismatch with FV3 (equiangular) — must be flagged in catalog and pilot models. |
| `cubed_sphere_panel_connectivity` | `grid` | `grid` | §7 | `src/grids/panel_connectivity.jl` | P1 | (any cubed-sphere model) | needs panel-edge-rotation table | M | looks_ok | `discretizations/grids/cubed_sphere_panel_connectivity.esm` | `cubed_sphere_grid_equidistant` | Edge rotations + 2×2 rotation matrices. |
| `cubed_sphere_metric_tensors` | `grid` | `coord_transform` | §7 | `src/grids/metric_tensors.jl` | P1 | (any cubed-sphere model) | needs metric-tensor field type | M | looks_ok | `discretizations/grids/cubed_sphere_metric.esm` | `cubed_sphere_grid_equidistant` | J, g^{ξξ}, g^{ηη}, g^{ξη} on the gnomonic cube. |
| `cubed_sphere_ghost_cells` | `grid` | `bc` | §7 | `src/ghost_cells.jl` | P1 | (any cubed-sphere model) | needs ghost-region rule type | M | needs_review | `discretizations/grids/cubed_sphere_ghost_cells.esm` | `cubed_sphere_panel_connectivity` | Multi-panel halo exchange. Boundary corner cases historically tricky — recent `audit-arrayop-fixes` PR (#15) suggests recent issues; revisit before promotion. |
| `colella_woodward_limiter` | `limiter` | `limiter` | §7 | `src/operators/reconstruction.jl`; Colella & Woodward (1984) | P1 | atmos.dyn (PPM-based transport), atmos.chem.transport | needs limiter rule type | M | needs_review | `discretizations/limiters/colella_woodward.esm` | `cam5_fv_ppm_reconstruction` | Implementation uses cell-mean as boundary fallback (low-order); flag in catalog. |
| `lax_friedrichs_flux` | `finite_volume` | `scheme` | §7 | `src/operators/flux_1d.jl`, `transport_2d.jl`; mirrored as `discretizations/finite_volume/lax_friedrichs_flux.json` | P2 | atmos.dyn (reference / debugging) | none for the in-panel Cartesian core; the cubed-sphere wrapper depends on schema follow-ups (boundary_policy + face-stagger / time-varying bindings) tracked off dsc-35x | S | looks_ok | `discretizations/finite_volume/lax_friedrichs_flux.json` | — | Now in repo as the 2-point cartesian flux stencil F_{i+1/2} = max(c,0)·q_i + min(c,0)·q_{i+1} with face-Courant `$c` carried as a per-face binding (analogous to limiter `$r`); upwind branch encoded directly via `abs`. Cubed-sphere wrapper (`flux_1d` ArrayOp on CubedSphereGrid) awaits schema infra and remains in `src/operators/flux_1d.jl` for now. PPM-based flux-form transport sibling now mirrored as `flux_1d_ppm` (dsc-r1i); LF-1D cubed-sphere wrapper still pending its own follow-up. |
| `flux_1d_ppm` | `finite_volume` | `scheme` | §7 | `src/operators/flux_1d.jl` (`flux_1d_ppm!`, `flux_1d_ppm_arrayop`, `_build_ppm_face_expr_xi/eta`, `_ppm_limit_cw84_sym`); Colella & Woodward (1984) JCP §4; mirrored as `discretizations/finite_volume/flux_1d_ppm.json` | P2 | atmos.dyn, atmos.chem.transport (PPM transport reference) | needs face-staggered Courant binding contract + ghost-extended input contract + (cubed-sphere wrapper only) panel-boundary distance handling — same dsc-35x follow-ups as `lax_friedrichs_flux` plus a `ghost_width` rule-level field | M | looks_ok | `discretizations/finite_volume/flux_1d_ppm.json` | `lax_friedrichs_flux`, `ppm_reconstruction`, `colella_woodward_limiter` | High-order PPM sibling of `lax_friedrichs_flux` — same `op=flux` shape and per-face `$c`/`$v` binding contract, but composes 4th-order PPM edge interpolation (CW84 eq. 1.6), the CW84 §4 monotonicity limiter (closed-form ifelse AST, identical to `vertical_remap.json` `limiter` block), the Courant-fraction flux integral over the swept upwind volume, and `ifelse(c >= 0, int_left, int_right)` upwind selection. 6-point cartesian stencil at offsets `{-3, -2, -1, 0, 1, 2}` relative to the right cell of the interface. Layer-A canonical fixture under `flux_1d_ppm/fixtures/canonical/` pins single-face flux on a smooth sinusoidal profile (limiter inactive); both walker fixtures declare `applicable: false` until schema gaps land — the same hand-pinned values are exercised inline in `test/test_flux_1d_ppm_rule.jl` against `_ppm_limit_cw84_sym` + the closed-form Courant integral. Imperative path retained in `src/operators/flux_1d.jl` (callers in `transport_2d.jl` + tests) until a future bead lands the schema extensions and migrates callers — same precedent as `vertical_remap_lin_2004`. |
| `transport_2d` | `finite_volume` | `scheme` | §7 | `src/operators/transport_2d.jl` (`transport_2d` ArrayOp); mirrored as `discretizations/finite_volume/transport_2d.json` | P2 | atmos.chem.transport (cubed-sphere reference) | needs cubed_sphere C-grid face-staggered bindings in walker harness | S | looks_ok | already exists (`discretizations/finite_volume/transport_2d.json`) | `cubed_sphere_grid_equidistant`, `lax_friedrichs_flux` | First-order Rusanov / local Lax-Friedrichs flux-form 2D transport on the gnomonic cubed sphere with C-grid Courant bindings. 5-point in-panel stencil (cells (0,0), (±1,0), (0,±1)). Sibling to the 1D `lax_friedrichs_flux` row — this is the in-panel 2D specialization that ports the `transport_2d` ArrayOp from `src/operators/transport_2d.jl`. The Lin-Rood and unsplit-PPM operators in the same source file remain pending separate beads (see catalog rows `fv3_lin_rood_advection` and `cam5_fv_ppm_reconstruction`). |
| `vertical_remap_lin_2004` | `finite_volume` | `scheme` | §7 | `src/operators/vertical_remap.jl`; Lin (2004) MWR | P1 | atmos.dyn (Lagrangian-vertical models) | needs conservative-remap rule type | L | looks_ok | `discretizations/finite_volume/vertical_remap_lin_2004.esm` | `cam5_fv_ppm_reconstruction`, `colella_woodward_limiter` | Conservative PPM-based vertical remap. Already used by FV3 and CAM5-FV. |
| `fv3_d_grid_wind_ops` | `finite_volume` | `scheme` | §7 | `src/operators/wind_ops.jl` | P1 | atmos.dyn (FV3 winds) | needs covariant↔contravariant transform rule | L | needs_review | `discretizations/finite_volume/fv3_d_grid_wind_ops.esm` | `fv3_cd_grid_transformation`, `fv3_super_grid_sin_cos` | Wraps the C-D wind transformation + flux computation with sub-grid sin(α) upwind selection. |
| `cubed_sphere_super_grid` | `grid` | `grid` | §7 | `src/grids/super_grid.jl` | P1 | atmos.dyn (FV3-style models) | needs sub-cell sample-point rule | M | looks_ok | (covered by `fv3_super_grid_sin_cos`) | `cubed_sphere_grid_equidistant` | Same as Section D row but recorded here for audit completeness. |
| `discretization_pipeline_arrayop` | `finite_volume` | `scheme` | extension | `src/discretization.jl`, `src/equation_discretizer.jl` | P1 | (engine glue) | needs ArrayOp template type in §7 (matches MOL PR #531) | L | needs_review | (engine; may not need `.esm`) | `mol_array_discretization_strategy` | Existing pipeline uses ArrayOp directly. Verify it converges with MOL's approach so we share schema. |

---

## Section G — Deep dive: other earth-science discretizations

These are not in the four required sources but are widely used in atmospheric
chemistry / dynamics, ocean dynamics, land-surface, and reactive-transport
models. Including them now lets follow-up beads be filed against a single
source-of-truth manifest.

### G.1 — Limiters and high-resolution schemes

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `weno3_jiang_shu_uniform` | `finite_volume` | `scheme` | §7 | Jiang & Shu (1996) JCP | P2 | atmos.chem.transport (low-cost) | (same as WENO5) | L | not_in_repo | `discretizations/finite_volume/weno3.esm` | `weno5_jiang_shu_uniform` | Cheaper than WENO5; useful when WENO5 is overkill. |
| `weno7_jiang_shu_uniform` | `finite_volume` | `scheme` | §7 | Balsara & Shu (2000) JCP | P3 | atmos.dyn (research, sharp features) | (same as WENO5) | XL | not_in_repo | `discretizations/finite_volume/weno7.esm` | `weno5_jiang_shu_uniform` | Justification (P3): higher-order WENO is research-grade for atmosphere; cost rarely worth it. |
| `eno_3rd_order_uniform` | `finite_volume` | `scheme` | §7 | Harten, Engquist, Osher, Chakravarthy (1987) | P3 | atmos.dyn (research) | needs ENO stencil-selection rule | L | not_in_repo | `discretizations/finite_volume/eno_3rd.esm` | — | Justification (P3): largely superseded by WENO. Keep for completeness. |
| `tvd_minmod_limiter` | `limiter` | `limiter` | §7 | Roe (1986) Ann. Rev. Fluid Mech. | P1 | atmos.chem.transport, ocean.dyn | needs slope-limiter rule type | S | not_in_repo | `discretizations/limiters/minmod.esm` | — | Most-conservative TVD limiter; baseline. |
| `tvd_superbee_limiter` | `limiter` | `limiter` | §7 | Roe (1986) | P1 | atmos.chem.transport (sharp fronts) | (as minmod) | S | not_in_repo | `discretizations/limiters/superbee.esm` | `tvd_minmod_limiter` | Sharper than minmod; can compress smooth profiles — use with care. |
| `tvd_van_leer_limiter` | `limiter` | `limiter` | §7 | van Leer (1974) JCP | P1 | atmos.chem.transport, ocean.dyn | (as minmod) | S | not_in_repo | `discretizations/limiters/van_leer.esm` | `tvd_minmod_limiter` | Smooth limiter — common default. |
| `tvd_mc_limiter` | `limiter` | `limiter` | §7 | van Leer (1977) JCP | P2 | atmos.chem.transport | (as minmod) | S | not_in_repo | `discretizations/limiters/mc.esm` | `tvd_minmod_limiter` | Monotonized central — sharper than van Leer, smoother than superbee. |
| `fct_zalesak` | `limiter` | `limiter` | §7 | Zalesak (1979) JCP | P2 | atmos.chem.transport (positivity) | needs antidiffusive-flux rule type | L | not_in_repo | `discretizations/limiters/fct_zalesak.esm` | `lax_wendroff` | Flux-Corrected Transport — positivity-preserving, monotonic. |
| `donor_cell_upwind_1st` | `finite_volume` | `scheme` | §7 | Courant, Isaacson, Rees (1952) | P1 | atmos.chem.transport (reference, debugging) | none | S | not_in_repo | `discretizations/finite_volume/donor_cell_upwind.esm` | — | First-order upwind in flux form. Universal reference. |
| `lax_wendroff` | `finite_volume` | `scheme` | §7 | Lax & Wendroff (1960) | P2 | atmos.chem.transport (linear test) | none | S | not_in_repo | `discretizations/finite_volume/lax_wendroff.esm` | — | Second-order linear scheme; oscillatory — pair with FCT. |

### G.2 — Grid staggerings

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `arakawa_a_grid` | `grid` | `staggering` | §7 | Arakawa & Lamb (1977) | P2 | atmos.dyn (collocated; spectral cores often A-equivalent) | needs collocated-staggering rule | S | not_in_repo | `discretizations/grids/arakawa_a.esm` | — | All variables co-located. |
| `arakawa_b_grid` | `grid` | `staggering` | §7 | Arakawa & Lamb (1977) | P2 | ocean.dyn (older POM, MOM4) | needs corner-velocity staggering | S | not_in_repo | `discretizations/grids/arakawa_b.esm` | — | Velocities at cell corners. |
| `arakawa_c_grid` | `grid` | `staggering` | §7 | Arakawa & Lamb (1977) | P1 | atmos.dyn, ocean.dyn (MITgcm, MOM6, MPAS) | (already covered for cubed-sphere & MPAS) | S | partial | `discretizations/grids/arakawa_c.esm` | — | Universal C-grid baseline; specializes to FV3-C and MPAS-C. |
| `arakawa_d_grid` | `grid` | `staggering` | §7 | Arakawa & Lamb (1977) | P2 | atmos.dyn (FV3) | (already covered for cubed-sphere) | S | partial | `discretizations/grids/arakawa_d.esm` | `fv3_d_grid_staggering` | Covered for cubed-sphere; needed for general Cartesian D-grid too. |
| `arakawa_e_grid` | `grid` | `staggering` | §7 | Arakawa & Lamb (1977) | P3 | atmos.dyn (NMM, some operational mesoscale) | needs rotated/E-grid staggering | M | not_in_repo | `discretizations/grids/arakawa_e.esm` | `arakawa_b_grid` | Justification (P3): rare today; legacy NMM/Eta models. |

### G.3 — Grids and topologies

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `lat_lon_regular_grid` | `grid` | `grid` | §7 | textbook (e.g. Durran 2010) | P1 | atmos.chem (CTMs like GEOS-Chem), ocean.bgc | none | S | not_in_repo | `discretizations/grids/lat_lon_regular.esm` | — | Universal CTM grid. Most existing chemistry models live here. |
| `lat_lon_reduced_grid` | `grid` | `grid` | §7 | ECMWF IFS notes | P2 | atmos.dyn (IFS reduced Gaussian) | needs variable-longitude-count rule | M | not_in_repo | `discretizations/grids/lat_lon_reduced.esm` | `lat_lon_regular_grid` | Fewer cells near poles to keep grid quasi-uniform. |
| `gaussian_grid` | `grid` | `grid` | §7 | spectral GCM textbooks | P2 | atmos.dyn (spectral models — IFS, GFS legacy) | needs Gauss-quadrature node rule | M | not_in_repo | `discretizations/grids/gaussian.esm` | — | Latitudes at Gaussian quadrature nodes. |
| `icosahedral_geodesic_grid` | `grid` | `grid` | §7 | Williamson (1968); Sadourny et al. (1968); ICON | P2 | atmos.dyn (ICON, Nicam), ocean.dyn | needs hexagonal-cell topology rule | L | not_in_repo | `discretizations/grids/icosahedral_geodesic.esm` | — | Geodesic mesh from icosahedron subdivision. ICON model. |
| `rotated_pole_lat_lon` | `grid` | `coord_transform` | §7 | regional climate models (HIRHAM, COSMO) | P2 | atmos.dyn (regional) | needs pole-rotation transform | S | not_in_repo | `discretizations/grids/rotated_pole.esm` | `lat_lon_regular_grid` | Lat-lon with shifted pole — quasi-uniform over a region. |
| `yin_yang_overset_grid` | `grid` | `grid` | §7 | Kageyama & Sato (2004); JMA NICAM-Y | P3 | atmos.dyn (research overset), planetary | needs overset-mesh communication rule | XL | not_in_repo | `discretizations/grids/yin_yang.esm` | `lat_lon_regular_grid` | Justification (P3): two overlapping low-latitude lat-lon patches; specialized. |
| `sigma_vertical_coord` | `grid` | `coord_transform` | §7 | Phillips (1957) | P1 | atmos.dyn (early GCM legacy), ocean.dyn (POM) | needs vertical-coord rule | S | not_in_repo | `discretizations/grids/sigma_vertical.esm` | — | Pure terrain-following. |
| `z_vertical_coord` | `grid` | `coord_transform` | §7 | textbook | P1 | ocean.dyn (z-coordinate ocean models) | (as sigma) | S | not_in_repo | `discretizations/grids/z_vertical.esm` | — | Geopotential height vertical levels. |
| `z_star_vertical_coord` | `grid` | `coord_transform` | §7 | Adcroft & Campin (2004) | P2 | ocean.dyn (MOM6, NEMO option) | (as sigma) | M | not_in_repo | `discretizations/grids/z_star_vertical.esm` | `z_vertical_coord` | Stretched-z to follow free surface. |
| `isopycnal_vertical_coord` | `grid` | `coord_transform` | §7 | Bleck (2002) HYCOM | P2 | ocean.dyn (HYCOM) | needs Lagrangian-density layering | L | not_in_repo | `discretizations/grids/isopycnal_vertical.esm` | — | Density-following ocean layers. |
| `ale_vertical_coord` | `grid` | `coord_transform` | §7 | Adcroft (MOM6 ALE) | P2 | ocean.dyn (MOM6) | needs general Lagrangian-Eulerian remap | L | not_in_repo | `discretizations/grids/ale_vertical.esm` | `vertical_remap_lin_2004` | Generalized vertical coord with remap. |

### G.4 — High-order and unstructured

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `dg_modal_legendre` | `dg` | `scheme` | extension | Cockburn & Shu (1998) | P2 | atmos.dyn (research, NUMA, DG-cubed-sphere) | needs DG basis + flux + Riemann-solver rule type | XL | not_in_repo | `discretizations/dg/dg_modal_legendre.esm` | — | Discontinuous Galerkin with modal Legendre basis. |
| `dg_nodal_lagrange` | `dg` | `scheme` | extension | Hesthaven & Warburton (2008) | P2 | atmos.dyn (research) | (as dg_modal_legendre) | XL | not_in_repo | `discretizations/dg/dg_nodal_lagrange.esm` | `dg_modal_legendre` | Nodal DG — basis at GLL points. |
| `spectral_element_method_general` | `spectral` | `scheme` | extension | Patera (1984); Karniadakis & Sherwin (2005) | P2 | atmos.dyn (CAM-SE/HOMME), ocean.dyn (MPAS-O alt) | (as homme) | XL | partial | `discretizations/spectral/spectral_element.esm` | `homme_continuous_galerkin_sem` | General SEM rule type — HOMME is one instance. |
| `lagrangian_semi_implicit_gcm` | `semi_lagrangian` | `scheme` | extension | Robert (1981); Staniforth & Cote (1991) | P2 | atmos.dyn (IFS, UM) | needs SI + SL coupling rule | XL | not_in_repo | `discretizations/semi_lagrangian/lagrangian_semi_implicit.esm` | — | Semi-Lagrangian + semi-implicit time stepping; backbone of IFS/UM. |
| `immersed_boundary_method` | `finite_volume` | `bc` | extension | Mittal & Iaccarino (2005) | P3 | coastal-ocean, complex topography | needs cut-cell or forcing-term rule | XL | not_in_repo | `discretizations/finite_volume/immersed_boundary.esm` | — | Justification (P3): mostly coastal/engineering CFD; rare in geophysical climate models, but emerging in coastal-ocean. |
| `amr_block_structured` | `grid` | `grid` | extension | Berger & Colella (1989) | Deferred | atmos.dyn (research), reactive flow | needs hierarchical-mesh rule | XL | not_in_repo | (deferred) | — | Deferred: out of scope for v1; tag for v2 once base content stabilizes. |
| `amr_octree` | `grid` | `grid` | extension | p4est (Burstedde et al. 2011) | Deferred | atmos.dyn (research) | (as block-structured AMR) | XL | not_in_repo | (deferred) | — | Deferred: same reason. |

### G.5 — Time integrators (selected, paired with spatial schemes)

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `forward_euler_time` | `time` | `time_integrator` | extension | textbook | P1 | reference / debugging | needs time-integrator rule type | S | not_in_repo | `discretizations/time/forward_euler.esm` | — | Trivial baseline. |
| `rk4_time` | `time` | `time_integrator` | extension | textbook | P1 | atmos.chem (chemistry boxes), generic | (as forward_euler) | S | not_in_repo | `discretizations/time/rk4.esm` | — | Classic 4-stage RK. |
| `ssp_rk3_time` | `time` | `time_integrator` | extension | Shu & Osher (1988) | P1 | atmos.dyn (with WENO/DG) | (as forward_euler) | S | not_in_repo | `discretizations/time/ssp_rk3.esm` | — | Strong Stability Preserving RK — pairs with WENO/DG. |
| `backward_euler_time` | `time` | `time_integrator` | extension | textbook | P1 | atmos.chem (stiff reactions), reactive transport | needs implicit-solve rule type | M | not_in_repo | `discretizations/time/backward_euler.esm` | — | Stiff baseline; needed for any reaction-diffusion. |
| `dirk_time` | `time` | `time_integrator` | extension | Hairer & Wanner (1996) | P2 | atmos.chem (stiff), atmos.dyn (semi-implicit) | (as backward_euler) | L | not_in_repo | `discretizations/time/dirk.esm` | `backward_euler_time` | Diagonally Implicit RK — higher-order stiff integration. |
| `crank_nicolson_time` | `time` | `time_integrator` | extension | textbook | P1 | atmos.dyn (diffusion), ocean.dyn | (as backward_euler) | S | not_in_repo | `discretizations/time/crank_nicolson.esm` | `backward_euler_time` | 2nd-order trapezoidal; common with parabolic terms. |
| `strang_splitting` | `time` | `time_integrator` | extension | Strang (1968) | P1 | atmos.chem (transport ⊕ chemistry), ocean.bgc | needs operator-splitting rule type | M | not_in_repo | `discretizations/time/strang_splitting.esm` | — | 2nd-order operator splitting. Universal in chemistry-transport. |
| `predictor_corrector_time` | `time` | `time_integrator` | extension | textbook | P2 | atmos.dyn (legacy), ocean.dyn | (as forward_euler) | S | not_in_repo | `discretizations/time/predictor_corrector.esm` | — | Adams-Bashforth-Moulton family. |

### G.6 — Solvers (paired with implicit time integrators)

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `multigrid_geometric` | `solver` | `solver` | extension | Brandt (1977) | P2 | atmos.dyn (elliptic pressure solve), ocean.dyn | needs solver rule type | XL | not_in_repo | `discretizations/solvers/multigrid_geometric.esm` | — | Standard elliptic solver; needed by any semi-implicit core. |
| `gmres_krylov` | `solver` | `solver` | extension | Saad & Schultz (1986) | P2 | atmos.chem (implicit chemistry), generic | (as multigrid) | L | not_in_repo | `discretizations/solvers/gmres.esm` | — | Generic Krylov solver; baseline for nonsymmetric systems. |

### G.7 — Boundary conditions (beyond periodic)

| name | family | kind | rfc_section | source_ref | priority | earth_sci_relevance | schema_gaps | complexity | audit_status | target_path | depends_on | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `dirichlet_bc` | `bc` | `bc` | §7 | textbook | P1 | universal | none | S | not_in_repo | `discretizations/bc/dirichlet.esm` | — | Specified value at boundary. |
| `neumann_bc` | `bc` | `bc` | §7 | textbook | P1 | universal | none | S | not_in_repo | `discretizations/bc/neumann.esm` | — | Specified gradient/flux. |
| `robin_bc` | `bc` | `bc` | §7 | textbook | P1 | atmos.chem (deposition), atmos.dyn (radiation BC) | none | S | not_in_repo | `discretizations/bc/robin.esm` | `dirichlet_bc`, `neumann_bc` | Mixed; common in deposition/exchange parameterizations. |
| `sponge_layer_bc` | `bc` | `bc` | extension | textbook (e.g. Klemp & Lilly 1978) | P2 | atmos.dyn (model lid), regional models (lateral) | needs damping-profile rule type | M | not_in_repo | `discretizations/bc/sponge_layer.esm` | `cam5_fv_div_damping_2nd` | Newtonian damping near upper/lateral boundary; named in CAM5 docs as a complementary scheme. |
| `radiation_open_bc` | `bc` | `bc` | extension | Orlanski (1976) | P2 | atmos.dyn (regional), ocean.dyn (open boundaries) | needs phase-speed rule type | M | not_in_repo | `discretizations/bc/radiation_open.esm` | — | Outflow boundary for waves. |
| `no_flux_bc` | `bc` | `bc` | §7 | textbook | P1 | atmos.chem (sealed top, surface), ocean.bgc | (as neumann) | S | not_in_repo | `discretizations/bc/no_flux.esm` | `neumann_bc` | Specialized Neumann (zero gradient/flux). Universal in chemistry models. |

### G.8 — Deferred (out of scope for v1)

| name | reason |
|---|---|
| `amr_block_structured` | Adaptive mesh refinement is a large content area; tag for v2. |
| `amr_octree` | Same reason. |
| `amr_patch_based` | Same reason; not even listed in main table. |
| `dg_hp_adaptive` | hp-adaptive DG is research-grade; not v1. |
| `dynamic_load_balancing` | Engine concern, not a discretization. |
| `multilayer_neural_pde_correction` | ML/PDE hybrid — explicitly out of scope for v1. |
| `lattice_boltzmann` | LBM exists in geophysics research but is rare; defer. |
| `sph_smoothed_particle_hydrodynamics` | SPH used in some ocean coastal work; defer. |

---

## Authoring policy (AST-first)

All rule coefficients, edge interpolants, limiter ratios, and
reconstruction expressions are authored **directly as ExpressionNode
ASTs** (ESS §7). ESD evaluates them through
`EarthSciDiscretizations.eval_coeff`, a thin passthrough to the
EarthSciSerialization tree-walk evaluator — ESD does not carry a shadow
evaluator. Rule authors MUST NOT introduce per-binding helper functions
that reimplement math a rule should own. See `../AGENTS.md` "Authoring
discretization rules" and ESS `esm-spec.md` §9.2 for the `call`-op
decision tree.

## Cross-cutting notes for downstream beads

1. **Cubed-sphere variant mismatch.** `src/grids/cubed_sphere.jl` uses
   *equidistant* gnomonic; FV3 and CAM-SE use *equiangular*. Document both
   in the catalog (rows present), and the per-scheme follow-up beads must
   resolve which variant a given scheme assumes.
2. **`audit_status: needs_review` items** should not be promoted to
   authoritative `.esm` content without an explicit review bead. Each such
   row should spawn a paired audit bead before content work begins.
3. **MPAS divergence (`mpas_divergence_flux_form`)** is the canonical RFC
   §7.3 worked example. It must be the first MPAS rule produced — it is
   the ground truth that the schema is calibrated against.
4. **MOL ArrayDiscretization (`mol_array_discretization_strategy`)** should
   shape the §7 schema's emission model so our ArrayOp pipeline and MOL's
   share a serialization format. Coordinate with EarthSciSerialization.
5. **P0 rows** (`centered_2nd_uniform`, `upwind_1st`, `cam5_fv_ppm_reconstruction`,
   `fv_divergence_cubed_sphere`, `fv_gradient_cubed_sphere`) gate any pilot
   ESM model build. Promote these first.
6. **Limiter-as-first-class.** The catalog treats limiters as standalone
   rows because PPM, WENO, and TVD schemes all compose with a separately
   namable limiter. The §7 schema should support `limiter`-typed rules as a
   distinct kind from `scheme`-typed rules.
7. **Time integrators and solvers** appear here for completeness but may
   live in a sibling repo (e.g. EarthSciTimeSteppers). Decision is out of
   scope for this Phase 0 inventory; flag for follow-up.

---

## Provenance

- Inventory authored by polecat `quartz` under bead `dsc-dmy` (2026-04-19).
- This file is the Phase 0 deliverable. Each row is the seed for one or
  more follow-up beads (per-scheme content + per-scheme audit where
  flagged). The `target_path` column is the promised landing location for
  the resulting `.esm` file.
