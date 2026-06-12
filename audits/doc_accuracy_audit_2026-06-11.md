# Documentation accuracy audit — 2026-06-11

Scope: every documentation surface in the repo (top-level README/AGENTS, `docs/`,
`discretizations/` READMEs and SELECTOR_KINDS, conformance READMEs under `tests/`,
the python/rust/typescript binding READMEs, and docstrings/comments in `src/`)
compared against the actual code. Only contradictions between existing docs and
current code are listed — not missing documentation. Each finding was verified
against the source before inclusion.

Audit method: five parallel review agents, one per documentation slice, each
instructed to grep broadly before reporting to avoid false positives.

---

## High severity — docs that would fail or materially mislead

### H1. `FVCubedSphere` pipeline advertised but deleted
- `README.md:8-10` and `docs/src/index.md:13` both headline "the `FVCubedSphere`
  discretization pipeline" and an "operator registry".
- Reality: `FVCubedSphere` appears nowhere in `src/` — it was removed in commit
  `953119a` ("chore(cleanup): delete ESD shadow assembly", which deleted
  `src/discretization.jl` and its export). There is no operator registry either;
  only `project_initial_condition` (`src/bc_handler.jl:5`) survives from that
  pipeline.

### H2. `discretizations/README.md:59-61` understates the stencil lowerer
- Claim: "Currently supports `cartesian` selector kind; other kinds raise
  `ArgumentError`".
- Reality: `src/stencil_lowering.jl:29-32,165-193` supports five kinds —
  `cartesian`, `arakawa`, `latlon`, `cubed_sphere`, `vertical`. Only
  `indirect`/`reduction` raise.

### H3. `docs/README.md:37` — plot-render instruction writes to the wrong directory
- Claim: from `docs/`, run `python3 ../tools/render_doc_plots.py` to populate
  `static/plots/` before `hugo --source .`.
- Reality: the script defaults to the cwd-relative path `docs/static/plots`
  (`tools/render_doc_plots.py:1056`), so run from `docs/` it creates
  `docs/docs/static/plots/` and Hugo never sees the plots. Needs
  `--out static/plots` or cwd = repo root.

### H4. `python/README.md:23-25` claims stencil application the package excludes by design
- Claim: the rules runtime includes "structured stencil application".
- Reality: `python/src/earthsci_discretizations/rules/__init__.py:14-18`
  explicitly states stencil application is owned by `earthsci_toolkit` and "ESD
  does not duplicate them here"; the package exports only `eval_coeff`,
  `load_rule`, `Rule`, `StencilEntry`.

---

## Medium severity — inaccurate or inverted claims

### Stale "pending / blocked / in flight" status (the work has landed)

| Doc | Claim | Reality |
|---|---|---|
| `docs/src/tutorial.md:147-178` | Layer B blocked on ESS `esm-4gw`; exemplar fixture shown as `applicable: false` | `discretizations/finite_difference/centered_2nd_uniform/fixtures/convergence/input.esm` is `applicable: true`; the blocker landed |
| `discretizations/README.md:7-9` | Rules "parsed as opaque JSON" until "the ESS rule engine … lands" | Engine landed: `src/rule_eval.jl:38-41`, `src/stencil_lowering.jl:13-18`, `test/walk_esd_tests.jl` all drive it; same README lines 45-47 contradict the claim |
| `discretizations/finite_difference/README.md:21-22` | "non-uniform vertical metrics are deferred" | `centered_2nd_nonuniform_vertical.json` exists with per-cell `dz[k]` coefficients; the cited SELECTOR_KINDS decision #3 says it's available |
| `discretizations/finite_difference/README.md:27-29` | latlon Layer B marks `applicable: false` | The convergence fixture has no `applicable` key (active, `expected_min_order: 1.9`); walker runner `2d_latlon_sphere` is live (`test/walk_esd_tests.jl:617-621,942`) |
| `discretizations/finite_volume/README.md:103-107` | `divergence_arakawa_c` Layer-B convergence is `applicable: false` | Inverted: the convergence fixture is active (`2d_arakawa_periodic` runner, `walk_esd_tests.jl:622-627`); it's the Layer-A **canonical** fixture that is `applicable: false` |
| `tests/conformance/grids/mpas/README.md:21-22,92-100,137-156` | Julia MPAS runtime "still landing" (dsc-7j0); conformance test is `@test_skip`; "three bindings" wired | `src/grids/mpas.jl` is wired (`src/EarthSciDiscretizations.jl:16,190`); `test/test_mpas_conformance.jl` (232 lines) runs the full comparison with no `@test_skip` |
| `tests/conformance/grids/latlon/README.md:3-4,26-27,130-131` | "The Julia `lat_lon` runtime is not yet on `main`" (dsc-1ts) | `src/grids/latlon.jl:49` + `grids.lat_lon` wiring and three unit-test files are on main; only the dedicated conformance runner is missing |
| `AGENTS.md:99` | `scripts/setup_polecat_env.sh` "no longer required for CI" | `.github/workflows/Tests.yml:52-54` still runs it as a mandatory step |

### Incorrect API / symbol references

- **M1.** `docs/GRIDS_API.md:80,126-133` — Python module path documented as
  `earthsci.grids` (`import earthsci`); the actual package is
  `earthsci_discretizations.grids` (`python/pyproject.toml:6`,
  `python/src/earthsci_discretizations/__init__.py`). Rust and TS match the doc;
  Python doesn't.
- **M2.** `discretizations/SELECTOR_KINDS.md:34-44` (decision #9) — lists
  `vertices_on_edge` among MPAS connectivity tables "verbatim" from
  `src/grids/mpas.jl`; no such field exists in any binding or in
  `grids/mpas.schema.json` (real fields: `src/grids/mpas.jl:97-118`). A rule
  referencing it would hit an unbound symbol.
- **M3.** `discretizations/SELECTOR_KINDS.md:184` (decision #12) — claims
  `cos_lat` is bound via `{op:"index", ...}` in
  `centered_2nd_uniform_latlon.json`; the rule binds it as a bare per-row scalar
  symbol (confirmed by `test/walk_esd_tests.jl:617-620`).
- **M4.** `discretizations/grids/mpas/README.md:86` — regeneration workflow says
  to re-run `test/test_mpas_fixtures.jl`, which does not exist (actual:
  `test_mpas_grid.jl` / `test_mpas_trait.jl` / `test_mpas_conformance.jl`).
- **M5.** `tests/conformance/grids/latlon/rules/centered_2nd_uniform_latlon/README.md:3-7`
  — claims the fixture verifies Julia `eval_coeff` plus Python/Rust siblings;
  only the two TypeScript tests consume this fixture dir (the landed Python/Rust
  rule-conformance tests use different harnesses).

### Broken links / stale counts in docs

- **M6.** `docs/src/fv_method.md:42,117` and `docs/src/tutorial.md:5,106,184,191,218`
  — Hugo `{{< ref >}}` shortcodes inside Documenter-built pages
  (`docs/make.jl:20-26`); Documenter doesn't process them and Hugo doesn't read
  these files, so every such link renders broken.
- **M7.** `docs/README.md:54-59` — says 11 stencil diagrams and 3 convergence
  plots with "the other 8 rules" gated; `tools/render_doc_plots.py:44-67`
  defines 13 rules, 6 convergence plots, 7 gated.

### Docstrings contradicting implementations (`src/`)

- **M8.** `src/stencil_lowering.jl:76-78` — latlon doc example says entries
  contribute `index(operand, lon_arg, lat_arg)`; the code sorts axis names
  (`lat` < `lon`) and emits `index(operand, lat_arg, lon_arg)`
  (`stencil_lowering.jl:510-531,575-588`; `test/test_stencil_lowering.jl:754`).
  The example contradicts the sorted-order rule stated one sentence earlier.
- **M9.** `src/stencil_lowering.jl:82-87` — same transposition for cubed-sphere:
  doc says `index(operand, xi_arg, eta_arg)`; code sorts `eta` < `xi` and emits
  `index(operand, eta_arg, xi_arg)` (`stencil_lowering.jl:639-661,704-718`).
- **M10.** `src/ghost_cells.jl:269-271` — `ghost_fill_arrayop` docstring
  recommends preparing input for `fv_laplacian_extended`, which doesn't exist;
  the real `fv_laplacian` (`src/operators/laplacian.jl:25-79`) takes a
  **non**-ghost-extended interior array, so following the docstring gives wrong
  shapes.
- **M11.** `src/grids/cartesian.jl:362-377` — `as_meshes` docstring (and its own
  error message) says loading Meshes.jl enables a `Meshes.CartesianGrid` view;
  the single method unconditionally throws and there is no package extension
  (`Project.toml` has no `[weakdeps]`/`[extensions]`), so installing Meshes.jl
  changes nothing. Related: `docs/GRIDS_API.md:433-439` claims `as_meshes` for
  three grid families plus `as_climacore`, none of which exist.
- **M12.** `src/grids/metric_tensors.jl:172-191` — `compute_coord_jacobian`
  claims "a smooth regularization" / "avoid a hard clamp discontinuity", but the
  regularized formula applies only inside an `abs(det) < 1e-10` branch, creating
  a ~10⁶× jump in scale at the threshold; the "≈ 1/det for |det| >> √ε"
  justification can never hold inside that branch.
- **M13.** `docs/src/index.md:12` — operator list still includes divergence and
  gradient; the `fv_divergence`/`fv_gradient_xi/eta` builders were retired
  (`src/EarthSciDiscretizations.jl:32-35`) and `src/operators/` has no
  divergence/gradient operator.

---

## Low severity — stale naming / cosmetic

- **L1.** `AGENTS.md:61-64` — binding table names the ESS package
  `earthsci_serialization`; all three bindings actually import
  `earthsci_toolkit` / `earthsci-toolkit`
  (`python/src/earthsci_discretizations/rules/evaluator.py:19`,
  `rust/src/rule_eval.rs:16`, `typescript/src/rules/index.ts:15`). The Julia row
  also names `EarthSciSerialization.evaluate`; `src/rule_eval.jl:38-41` actually
  calls `parse_expression` + `evaluate_expr`.
- **L2.** `docs/GRIDS_API.md:108-116` — Julia example calls
  `EarthSciDiscretizations.grids.cubed_sphere(...)`; the `grids` submodule
  (`src/EarthSciDiscretizations.jl:164-263`) has no `cubed_sphere` (Python/Rust/
  TS do). Same nonexistent function cited as fixture provenance in
  `discretizations/grids/cubed_sphere/README.md:45-47` (along with a
  nonexistent `to_esm(::CubedSphereGrid)`).
- **L3.** `docs/REPO_LAYOUT.md:37-40` — `discretizations/` tree omits the
  existing `grids/` subdirectory (grid-family schemas).
- **L4.** `docs/rule-catalog-roadmap.md:129,166-173,218-221` — `source_ref`s
  cite deleted files (`src/operators/divergence.jl`, `gradient.jl`,
  `kinetic_energy.jl`, `vorticity.jl`, `wind_ops.jl`, `vertical_remap.jl`,
  `src/discretization.jl`, `src/equation_discretizer.jl`), some marked "Already
  in repo". Mitigated by the file's own stale-snapshot disclaimer (lines 3-10).
- **L5.** `docs/GRIDS_API.md:675` and `discretizations/SELECTOR_KINDS.md:348` —
  both point to `docs/rule-catalog.md` as the hand-authored inventory/roadmap;
  that file is the gitignored autogenerated manifest; the roadmap is
  `docs/rule-catalog-roadmap.md`.
- **L6.** `discretizations/SELECTOR_KINDS.md:31` (decision #3) — cites the
  missing `1d_vertical_column` runner as the Layer-B blocker; that runner is
  implemented (`test/walk_esd_tests.jl:613-616`); the actual remaining blocker
  is different (the runner binds uniform `h` only).
- **L7.** `discretizations/SELECTOR_KINDS.md:64-65` — says MPAS connectivity
  tables "must exist in the grid's `connectivity` block; see
  `grids/mpas.schema.json`"; that schema defines no `connectivity` block (the
  grid is loader-backed).
- **L8.** `discretizations/finite_difference/README.md:28-29` — latlon
  coefficients quoted without the centered-difference factor 2 (actual:
  `±1/(2·R·cos_lat·dlon)`, `±1/(2·R·dlat)`; SELECTOR_KINDS.md:183 has it right).
- **L9.** `discretizations/finite_difference/README.md:30-36` —
  `covariant_laplacian_cubed_sphere` described with its pre-migration
  `stencil`/`selectors` structure; it is now `kind: "cross_metric"` with 8
  `terms` (commit `2e70560`).
- **L10.** `tests/conformance/grids/cartesian/rules/ppm_reconstruction/README.md:3-4,48-49`
  — Rust sibling marked "in flight";
  `rust/tests/ppm_reconstruction_rule_conformance.rs` landed.
- **L11.** `rust/README.md:25-28` — advertises `--features geo,proj`; the
  features are declared in `Cargo.toml:16-21` but no `#[cfg(feature)]` code
  exists anywhere in the crate, so enabling them adds nothing.
- **L12.** `python/README.md:15-17` — advertises the `[xarray]` extra; the extra
  is declared in `pyproject.toml:25` but neither `xarray` nor `pyproj` is
  imported anywhere in the package.
- **L13.** `src/operators/reconstruction.jl:9-11` — module docstring says the
  helpers are used by `vertical_remap.jl`, which was ported to a JSON rule and
  deleted; only `flux_1d.jl` still calls `_ppm_limit_cw84`.
- **L14.** `src/EarthSciDiscretizations.jl:48-49` — include comment labels
  `bc_handler.jl` "Boundary-condition handler"; the file contains only
  initial-condition projection.

---

## Verified clean (spot-checked, no mismatch)

All `@example` blocks in `docs/src/index.md`/`grid.md`/`operators.md` run
against current signatures; the `centered_2nd_uniform` JSON quoted in
`fv_method.md`/`tutorial.md` matches the rule file; the C-grid staggering table
matches `src/staggering.jl`; panel-connectivity prose matches
`src/grids/panel_connectivity.jl`; all TypeScript README npm scripts and the
`earthsci.grids` namespace exist; every conformance regenerate script and
golden referenced by the vertical/duo/arakawa/cubed-sphere/cartesian READMEs
exists; documented grid-builder signatures and defaults across
cartesian/vertical/latlon/arakawa/duo/mpas match their implementations;
`eval_coeff`'s thin-passthrough claim checks out against the ESS source;
`audits/maintainability_audit_2026-05-01.md` is a properly dated snapshot.

## Cross-cutting observations

1. **The dominant failure mode is stale status, not wrong math.** Roughly half
   the findings are "doc says X is pending/blocked/skipped" where X has since
   landed (Layer-B fixtures, MPAS/latlon runtimes, the stencil lowerer's
   selector coverage, Rust conformance tests). Docs written at bead-creation
   time aren't being revisited when the bead lands.
2. **Deleted-pipeline ghosts.** The `FVCubedSphere`/operator-registry removal
   (commit `953119a`) left dangling references in README.md, docs/src/index.md,
   the roadmap, and two docstrings (`fv_laplacian_extended`,
   `vertical_remap.jl`).
3. **Speculative API surface.** Several docs describe integrations that were
   specified but never implemented (`as_meshes`/`as_climacore`, rust `geo`/
   `proj` features, python `xarray` extra, Julia `grids.cubed_sphere`).
