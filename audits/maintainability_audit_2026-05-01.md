# EarthSciDiscretizations — Maintainability Audit

**Date:** 2026-05-01
**Bead:** dsc-ztvz
**Reviewer:** polecat onyx (EarthSciDiscretizations)
**Scope:** Review-only. No commits, no MRs. Deliverable lives in worktree only.
**Working tree:** branch `polecat/onyx-momzdtrd` against `main`

---

## Executive summary

ESD is in **good structural shape** post-dsc-c4j / dsc-7km / dsc-a8m / dsc-o05.
Single-pathway compliance for production rule evaluation is clean: every
binding's `rule_eval` is a thin passthrough to its language-native ESS
evaluator, no shadow op-dispatch tables, no parallel rule registries.

There are, however, **two live single-pathway anti-patterns** that the
absolute rule from CLAUDE.md flags by name:

1. `tools/render_doc_plots.py` re-implements rule math in numpy to draw
   convergence plots (the "homebrew doc-build simulator" anti-pattern).
2. `test/walk_esd_tests.jl` carries a local AST evaluator (`_eval_limiter_ast`)
   plus a 1-D MUSCL TVD sweep (`_run_tvd_advection`) used by Layer-B'.

Both are documented as "illustrative" / "allowed by dsc-8vu risk section,"
but the absolute rule is unconditional: doc builds and validators MUST go
through an official ESS runner. These should be tracked as in-flight
debt, not retired by self-disclaimer.

The walker harness itself (`test/walk_esd_tests.jl`, `test/test_esd_walker.jl`,
~2K lines) is structurally sound: rule discovery shares the production
`load_rules` path with no enum divergence; layer dispatch is fixture-driven;
12 fixture-declared `applicable: false` skips are correctly gated. Coverage
is uneven — Layer A 5%, Layer B 21%, Layer B' 5%, Layer D 5% — but most
gaps are blocked on ESS harness extensions, not ESD work.

Cross-binding alignment is consistent across Julia / Python / Rust /
TypeScript (the bead's "Python + Rust mirrors" framing is now incomplete:
TypeScript exists). The principal cross-binding gap is **Python rule-test
parity** (Julia has 9 rule test files, Python has 0) plus a **floating ESS
dep** in `python/pyproject.toml` that pins `@main` while Julia pins `0.0.3`.

The doc surface is mostly current; `discretizations/SELECTOR_KINDS.md` and
`docs/GRIDS_API.md` reflect post-dsc-4ze + post-dsc-a8m state.

**Top three concerns**, in priority order:

1. **`render_doc_plots.py` shadow numerics** (P0 by absolute rule, P1
   pragmatically — it gates doc convergence plots).
2. **Walker limiter shadow eval (`_eval_limiter_ast` + `_run_tvd_advection`)**
   (P1 — scoped to Layer-B', but a documented exception that should be
   retired when ESS lands `max`/`min`).
3. **Python ESS pin floats `@main`** (P1 — Python CI breakage decoupled
   from Julia/Rust/TS pins).

Everything else is P2 or below. A 30-line summary lives in the bd notes
on dsc-ztvz.

---

## Findings

### F1 (P0/P1 — single-pathway): `tools/render_doc_plots.py` re-implements rule math

- **File:** `tools/render_doc_plots.py:818–1015`
- **Current state:** Six `convergence_*` functions
  (`convergence_centered_2nd`, `…_vertical`, `…_latlon`,
  `convergence_upwind_1st`, `convergence_ppm_reconstruction`,
  `convergence_weno5_advection_2d`) re-implement the rule formula in
  numpy and produce L∞-error / order-fit plots written to
  `docs/static/plots/rules/<rule>-convergence.png`. The script's docstring
  (lines 17–22) self-disclaims authority and points to the ESS walker as
  the oracle.
- **Why it matters:** CLAUDE.md ABSOLUTE rule lists "Homebrew doc-build
  simulator" as a named anti-pattern (table at the top of the file), with
  the explicit rider "Doc builds may simulate, but only via an official
  ESS runner." A self-disclaimer in a docstring is not the relief valve
  the rule provides for. The risk is real: divergence between this script
  and the rule JSON silently mis-shapes user-facing plots, and the
  WENO5 reconstruction implementation here (lines 960–980 — Jiang–Shu
  weights, β-smoothness indicators, ε regularization) is a fully fledged
  numerical kernel, not "illustrative shape."
- **Severity:** **P0 by the absolute rule; pragmatic P1** (the script
  is plot-generation, not a runtime).
- **Effort:** **M.** Replace the convergence plotters with a thin
  consumer of the ESS walker's per-rule convergence output
  (the walker already emits `(n, l_inf)` pairs through Layer B; expose
  them as a JSON artifact and have the doc-plot script just plot them).
  The PPM/weno5 helpers in the script can be deleted in the same change.
  The grid-visualisation portions of the file
  (`SphericalVoronoi`-driven panel diagrams, etc.) are unrelated and
  should stay.
- **Bead:** file as `audit-from-dsc-ztvz: retire shadow numerics in
  render_doc_plots.py; consume ESS walker output instead.`

### F2 (P1 — single-pathway): walker carries `_eval_limiter_ast` + `_run_tvd_advection`

- **File:** `test/walk_esd_tests.jl:614–836`
- **Current state:** A self-contained AST evaluator
  (`_eval_limiter_ast`, lines 625–655) supporting `+ - * / abs max min`,
  used by `run_monotonicity_check` (lines 668–764) and
  `_run_tvd_advection` (lines 774–836) to drive Layer-B' Sweby +
  TVD checks for the slope-ratio limiters
  (`flux_limiter_minmod`, `flux_limiter_superbee`). The header comment
  (lines 614–623) explicitly cites the dsc-8vu risk section and frames
  this as a stand-in until ESS evaluators gain `max` / `min`.
- **Why it matters:** This is the "test-path evaluator" anti-pattern named
  in CLAUDE.md ("Replace by driving the canonical pipeline + evaluating
  the resulting AST"). The walker is the audit harness for the
  single-pathway invariant, so any local evaluator inside it is an
  internal contradiction even if guarded.
- **Severity:** **P1.** Scoped to Layer-B' (limiters), bounded fix
  (one bead in ESS Julia evaluator), no production exposure.
- **Effort:** **S** once ESS lands `max` / `min`: replace the body of
  `_eval_limiter_ast` calls with `EarthSciDiscretizations.eval_coeff`,
  delete the local evaluator, keep `_run_tvd_advection` as a fixture
  consumer that calls `eval_coeff` per slope-ratio.
- **Bead:** file as `audit-from-dsc-ztvz: retire walker-local AST eval
  once ESS gains max/min`. Cross-link to dsc-8vu (already exists per
  the comment).

### F3 (P1 — cross-binding): Python ESS dependency floats `@main`

- **File:** `python/pyproject.toml:21`
- **Current state:**
  `earthsci-toolkit @ git+https://github.com/EarthSciML/EarthSciSerialization.git@main#subdirectory=packages/earthsci_toolkit`.
  Julia `Project.toml:14` pins `EarthSciSerialization = "0.0.3"`. Rust
  `Cargo.toml:31` and TypeScript `package.json` both use a local
  `path:` dependency on the workspace checkout.
- **Why it matters:** Three different version-resolution regimes across
  bindings means a single ESS commit can break Python CI in isolation,
  and `dsc-m20` (hatchling pin in flight) would not have surfaced if the
  pin had been a tagged version. With dsc-3kz already landing the
  hatchling floor, the pyproject is otherwise clean.
- **Severity:** **P1.**
- **Effort:** **S.** Replace `@main` with the same tag/SHA that
  `Project.toml` resolves to, or land an ESS-native release tag
  (`v0.0.3`-style) and align all three binding pins.
- **Bead:** file as `audit-from-dsc-ztvz: align Python ESS pin with
  Julia compat (avoid @main)`.

### F4 (P2 — cross-binding): Python has no rule-level tests

- **Files:** `python/tests/` (16 files, all grid- or fixture-level).
  Julia `test/` has 9 rule-specific tests
  (`test_rule_catalog.jl`, `test_flux_1d_ppm_rule.jl`,
  `test_flux_limiters_rule.jl`, `test_lax_friedrichs_flux_rule.jl`,
  `test_ppm_reconstruction_rule.jl`, `test_weno5_advection_rule.jl`,
  `test_weno5_advection_2d_rule.jl`, `test_rules.jl`, plus the walker
  itself). Rust has `lat_lon_rule_conformance.rs` and
  `ppm_reconstruction_rule_conformance.rs`. TypeScript has
  `rules.test.ts`, `centered_2nd_uniform_latlon.conformance.test.ts`,
  `ppm_reconstruction.conformance.test.ts`,
  `ppm_reconstruction.convergence.test.ts`.
- **Severity:** **P2.** The Python passthrough is exercised through
  conformance fixtures, but a divergence in Python ESS behaviour for a
  rule shape that no Python test covers will land silently.
- **Effort:** **M.** Add a Python `test_rule_catalog.py` that loads
  every rule in `discretizations/`, calls
  `earthsci_discretizations.eval_coeff` against the rule's canonical
  fixture inputs, and compares against the canonical fixture outputs
  byte-equivalently. This is the same shape the Julia
  `test_rule_catalog.jl` uses.
- **Bead:** file as `audit-from-dsc-ztvz: extend Python rule test suite
  to match Julia coverage`.

### F5 (P2 — coverage): 8 FV3 rules have no fixture coverage at all

- **Files:** `discretizations/finite_volume/fv3_*.json`
  (`fv3_sinsg_flux_xi`, `…_eta`, `fv3_vorticity_cellmean`,
  `fv3_kinetic_energy_cell`, `fv3_d_to_c_xi`, `fv3_d_to_c_eta`,
  `fv3_vorticity_corner`, `fv3_absolute_vorticity_cellmean`).
  No matching fixture directories under
  `discretizations/finite_volume/`.
- **Severity:** **P2.** dsc-247 ports were "reference artifacts." But
  Layer-A canonical fixtures cost minutes per rule and would catch
  a JSON-schema regression cheaply.
- **Effort:** **M** for all 8 (S per rule). Generate
  `fixtures/canonical/{input,expected}.esm` by running them through
  `discretize` once and committing the byte-equivalent output.
- **Bead:** file as `audit-from-dsc-ztvz: add Layer-A canonical
  fixtures to FV3 rule ports`.

### F6 (P2 — coverage): Layer A and Layer D are extremely sparse

- **Files:** Walker counts —
  Layer A 2/38 (`centered_2nd_uniform`,
  `centered_2nd_uniform_vertical` per the per-rule sweep done by the
  walker-coverage subagent),
  Layer D 2/38 (`divergence_arakawa_c`, `flux_limiter_minmod`).
- **Severity:** **P2.** Layer B (8/38 passing, several blocked on ESS
  harness extensions) carries most of the audit weight today, so this
  is opportunity rather than risk. Layer-A canonical fixtures are
  cheap (rerun `discretize` and snapshot); Layer-D conservation
  fixtures need fresh authoring per rule shape.
- **Effort:** **L** if applied broadly; **M** scoped to high-impact
  rules (Arakawa-C divergence, WENO-5 1D + 2D, PPM reconstruction).
- **Bead:** file as `audit-from-dsc-ztvz: lift Layer-A / Layer-D
  fixture coverage on hot-path rules`.

### F7 (P2 — coverage): 12 fixtures declare `applicable: false` awaiting ESS

- **Files:** Counted in the walker-coverage subagent table (rules:
  `covariant_laplacian_cubed_sphere`, `nn_diffusion_mpas`, `periodic_bc`,
  `flux_limiter_minmod`, `flux_limiter_superbee`, `flux_1d_ppm`,
  `lax_friedrichs_flux`, `lax_friedrichs_flux_cubed_sphere_eta`,
  `lax_friedrichs_flux_cubed_sphere_xi`, `ppm_edge_cubed_sphere`,
  `transport_2d`, `vertical_remap`).
- **Severity:** **P2.** Walker correctly skips, no false negatives.
  Genuine blockers split into three groups:
    - ESS multi-axis selectors / cubed_sphere metric bindings
      (covariant_laplacian, transport_2d, ppm_edge_cubed_sphere,
       lax_friedrichs cubed_sphere variants).
    - ESS MPAS unstructured harness
      (nn_diffusion_mpas — tracked under dsc-nc5 per the comment).
    - Architectural (periodic_bc index-rewrite,
      vertical_remap phase-hook) — Layer-B is genuinely N/A.
- **Effort:** Out-of-scope for ESD; tracked as ESS work.
- **Bead:** check if dsc-nc5 (MPAS MMS), esm-4gw (multi-axis walker),
  cubed_sphere harness extension are all filed against ESS; file
  any missing.

### F8 (P3 — operators helpers): `src/operators/reconstruction.jl` retains PPM helpers

- **File:** `src/operators/reconstruction.jl:27–159`
- **Current state:** `_ppm_limit_cw84`, `_ppm_limit_cw84_sym`,
  `ppm_reconstruction_arrayop`, `ppm_flux_integral`. Comment at
  `src/EarthSciDiscretizations.jl:34–36` documents
  `fv_divergence`, `fv_gradient_xi/eta`, `ppm_reconstruction!` /
  `ppm_reconstruction` retirement (dsc-o05) but the helpers above
  remain. Verified call sites: `flux_1d.jl:149,155,180,185,588,593,
  638,642`. Per AGENTS.md this is acceptable: helpers drive the
  schema-gated ArrayOp pipeline, the canonical math is in
  `discretizations/finite_volume/ppm_reconstruction.json`.
- **Severity:** **P3.** Acceptable today, but every Julia helper that
  could plausibly be derived from a JSON rule is an attractive
  shortcut for future contributors.
- **Effort:** Out of scope. Could be revisited if the rule engine grows
  ArrayOp-shaped output (then the helpers become obsolete).

### F9 (P3 — docs): rule-catalog.md is hand-maintained, not generated

- **File:** `docs/rule-catalog.md` (60K, last touch 2026-04-28).
- **Severity:** **P3.** Drift risk vs `discretizations/`; today the
  walker can list rules but doesn't emit the catalog table. A simple
  `julia --project=docs make.jl` could regenerate it.
- **Effort:** **S** to write the generator; **L** to migrate full
  authoring once.

### F10 (P3 — workflow): `Tests.yml` excludes binding paths

- **Files:** `.github/workflows/Tests.yml`, `Python.yml`, `Rust.yml`,
  `TypeScript.yml`. Each binding has its own workflow, gated by paths.
- **Severity:** **P3.** Intentional and correct (per-binding CI).
  Caveat: a refactor that crosses binding boundaries will fail to
  trigger all four workflows; `paths-ignore` in `Tests.yml` is the
  load-bearing config to watch.

### F11 (P3 — naming): scattered binding-name inconsistencies

- **Files:**
  - `src/grids/latlon.jl` (filename) vs
    `EarthSciDiscretizations.grids.lat_lon` (public name) vs
    Python `lat_lon` vs Rust `lat_lon` vs TypeScript `lat_lon`. The
    snake_case public name is consistent across bindings; only the
    Julia filename is `latlon.jl`. Acceptable.
  - Bead description says "Python + Rust mirrors"; actual surface is
    Python + Rust + TypeScript. Documentation is accurate
    (AGENTS.md:49); only the bead phrasing was outdated.

---

## Top 10 priorities

1. **F1 — Retire `tools/render_doc_plots.py` shadow numerics** (P1 pragmatic / P0 by absolute rule).
2. **F2 — Retire walker-local `_eval_limiter_ast` + `_run_tvd_advection`** once ESS lands `max`/`min` (P1).
3. **F3 — Pin Python ESS dep to a tag/SHA, not `@main`** (P1).
4. **F4 — Add Python rule-level test parity** (P2).
5. **F5 — Add Layer-A canonical fixtures to 8 FV3 rules** (P2).
6. **F6 — Lift Layer-A / Layer-D fixture coverage on hot-path rules** (P2).
7. **F7 — Track ESS unblock beads (multi-axis, cubed_sphere metrics, MPAS MMS)** (P2).
8. **F8 — Re-evaluate PPM helpers in `reconstruction.jl` once ArrayOp output is rule-derivable** (P3).
9. **F9 — Generate `docs/rule-catalog.md` from `discretizations/`** (P3).
10. **F11 — Update bead description templates ("Python + Rust" → "Python + Rust + TypeScript")** (P3).

---

## Single-pathway compliance — line-item summary

| Surface | Compliance | Evidence |
|---|---|---|
| Production rule_eval (Julia) | ✅ thin passthrough | `src/rule_eval.jl:38–41` |
| Production rule_eval (Python) | ✅ thin passthrough | `python/src/earthsci_discretizations/rules/evaluator.py:53–61` |
| Production rule_eval (Rust) | ✅ thin passthrough | `rust/src/rule_eval.rs:49–59` |
| Production rule_eval (TypeScript) | ✅ thin passthrough | `typescript/src/rules/index.ts:20–34` |
| Production operators helpers (Julia) | ✅ schema-gated ArrayOp | `src/operators/*.jl` |
| Walker rule discovery | ✅ shares `load_rules` | `test/walk_esd_tests.jl:44–47` |
| Layer A / B / C / D | ✅ ESS-driven | Layer B uses `EarthSciSerialization.verify_mms_convergence` |
| Layer B' (limiters) | ⚠️ local AST eval | `_eval_limiter_ast`, `_run_tvd_advection` (F2) |
| `regenerate_golden.{jl,py}` | ✅ ESS-driven | per per-script docstrings |
| `tools/render_doc_plots.py` | ⚠️ shadow numerics | `convergence_*` functions (F1) |
| Per-rule-shape dispatch outside engine | ✅ none found | grep clean |

---

## Build / dependency health

| Binding | Manifest | Pin style | Notes |
|---|---|---|---|
| Julia | `Project.toml:14` | `EarthSciSerialization = "0.0.3"` | Aligned with workspace; fine |
| Python | `pyproject.toml:21` | `git+https://...@main` | **F3** — floats |
| Rust | `Cargo.toml:31` | `path = "../../EarthSciSerialization/..."` | Workspace checkout; MSRV 1.88 (dsc-hru landed) |
| TypeScript | `package.json` | `file:../../EarthSciSerialization/...` | Workspace checkout |
| Python build | `pyproject.toml:2` | `hatchling>=1.27` | dsc-3kz landed |
| Rust edition | `Cargo.toml:4` | `2021` (uses ESS's `2024`) | dsc-hru bumped MSRV for ESS edition 2024 |

---

## Fixture coverage at a glance

| Family | Rules | Layer-A fixt | Layer-B fixt | Layer-B' | Layer-C | Layer-D |
|---|---|---|---|---|---|---|
| finite_difference | 8 | 2 | 4 (partial) | 0 | 0 | 0 |
| finite_volume | 23 | 0 | several (with applicable:false subset) | 2 | 1 | 2 |
| grids (schemas) | 7 | n/a | n/a | n/a | n/a | n/a |
| spectral | 0 | — | — | — | — | — |

(Counts cross-checked against `find discretizations -name '*.json' -not -path '*/fixtures/*'` and the walker subagent's per-rule sweep.)

---

## What's NOT a problem

- Production rule_eval surfaces in all four bindings — clean
  passthrough, no per-rule dispatch.
- `src/EarthSciDiscretizations.jl` include order — documented and
  correct; no cyclic deps.
- God-module exports in `EarthSciDiscretizations.jl` — intentional
  façade; submodule `grids` cleanly layers atop.
- `discretizations/SELECTOR_KINDS.md` — current.
- `docs/GRIDS_API.md` §4.3 — reflects post-dsc-4ze state, lists all
  four bindings.
- Recent ports (dsc-247 FV3, dsc-9yh PPM edge) are committed to
  JSON; no Julia-side residue.
- Retired ArrayOp builders (dsc-o05) — fully removed; only the
  comment marker survives.

---

## Notes for follow-on bead authors

- F1 + F2 are the only single-pathway findings. Both are bounded.
  Don't widen scope.
- F4 + F5 + F6 are coverage extensions; each can land as an
  independent bead.
- F7 wants beads filed against ESS, not ESD. Confirm dsc-nc5 +
  multi-axis-walker-extension already exist; file what's missing.
- Top-of-funnel beads to file from this audit:
  - `audit-from-dsc-ztvz: retire shadow numerics in render_doc_plots.py`
  - `audit-from-dsc-ztvz: retire walker-local AST eval once ESS gains max/min`
  - `audit-from-dsc-ztvz: pin Python ESS dep to a tag (not @main)`
  - `audit-from-dsc-ztvz: extend Python rule test suite for parity with Julia`
  - `audit-from-dsc-ztvz: add Layer-A canonical fixtures to FV3 rule ports`
  - `audit-from-dsc-ztvz: lift Layer-A / Layer-D fixture coverage on hot-path rules`
  - `audit-from-dsc-ztvz: generate docs/rule-catalog.md from discretizations/`

End of audit.
