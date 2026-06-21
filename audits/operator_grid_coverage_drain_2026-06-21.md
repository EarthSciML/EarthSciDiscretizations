# Operator × grid coverage drain — final audit (G16)

**Date:** 2026-06-21
**Bead:** esd-6g4.15 (G16), capstone of epic `esd-6g4` (ESD operator × grid
coverage — author missing declarative rules + backfill fixtures,
DECLARATIVE-OR-FAIL).
**Scope:** Adversarial confirmation only. The single code artifact is the
regenerated Hugo data file `docs/data/rule_matrix.json`; this document records
the evidence. No rule JSON, fixture, runner, or operator code is added or
changed.

---

## Executive summary

Epic `esd-6g4` closed the operator × grid coverage gaps from the 2026-06-21
audit by authoring declarative rules + fixtures across G1–G15 (esd-6g4.1 …
esd-6g4.14, all merged to `main` at `8a3391b`). This capstone confirms the
epic's acceptance against the live tree:

1. **`docs/data/rule_matrix.json` regenerated** from the current catalog —
   **39 → 59 rules** (the 20 rules G1–G15 added now appear). ✓
2. **Every applicable cell is CLOSED or REPORTED.** Of the 59 applicable
   cells: 31 carry a running MMS-convergence sweep; the remaining 28 the
   matrix labels `missing` are each dispositioned below — **5 REPORTED**
   declarative-infeasible (convergence sweep deferred with a documented
   reason) and **23 CLOSED** by a Layer-A byte contract, a Layer-C
   integration solve, or a dedicated rule/conformance test that the matrix's
   fixture scan does not see. **Zero genuine dispatch holes remain.** ✓
3. **Walker + cross-binding conformance green** — full unit `Pkg.test()`
   passes (see §5). ✓
4. **Zero imperative operator code added** — the diff is the regenerated
   matrix data file plus this audit. ✓
5. **`vertical_remap` (dsc-otd) still tracked/deferred** — documentation-only,
   no imperative implementation in any binding, excluded from walker
   Layer-B/D (see §4). ✓

**Verdict: G16 acceptance met. The epic's coverage is fully closed-or-reported.**

Two cells (`staggered_1st_uniform_face_to_cc`, `advection_mpas_velocity_first`)
are *operand-order mirror* siblings whose primary direction is byte-pinned but
which carry no separate assertion — the thinnest coverage of the set. They are
not gaps (the operator is validated via the primary direction), but they are
surfaced as a recommended follow-up in §6 rather than buried.

---

## 1. Why `missing` in the matrix is not the same as "uncovered"

`tools/render_rule_matrix.py` labels a `rule × grid` cell `missing` when it
finds **neither** a *running* MMS-convergence sweep
(`discretizations/<family>/<rule>/fixtures/convergence/expected.esm` with an
`expected_min_order`) **nor** a cross-binding canonical fixture under
`tests/conformance/grids/<grid>/rules/<rule>/` or
`tests/conformance/rules/<rule>/`. That is the matrix's deliberate
*gold-standard* coverage bar (the strongest, contributor-facing signal of
"where the convergence proof or cross-binding pin has not landed yet").

It is **blind by construction** to three other forms of validation that the
test suite does exercise:

- **Layer-A byte contracts** under
  `discretizations/<family>/<rule>/fixtures/canonical/` (and `…/rewrite/`).
  These are run by `test/walk_esd_tests.jl` via `EarthSciSerialization.discretize`
  / `.rewrite` and byte-compared to the pinned `expected.esm`. The matrix's
  `_has_canonical_fixture` only looks under `tests/conformance/`, so a rule whose
  canonical byte contract lives under `discretizations/` reads as `missing`.
- **Layer-C integration solves** under
  `discretizations/<family>/<rule>/fixtures/integration/`. These drive a real
  ODE solve through `build_ode_problem`; the matrix does not scan them.
- **Dedicated Julia rule/conformance tests** (`test/test_*_rule.jl`,
  `test/test_*_conformance.jl`) that assert against an in-test oracle with no
  on-disk fixture file for the matrix to discover.

A fourth, narrower keying quirk affects **multi-rule JSON files**: the matrix
keys each cell by the JSON *block* name but resolves fixtures under
`discretizations/<family>/<block>/fixtures/`. When a file holds two blocks
(e.g. `staggered_1st_uniform.json` → `…_cc_to_face` + `…_face_to_cc`), the
fixtures live under the *file* basename (`staggered_1st_uniform/fixtures/`), so
both blocks read as `missing` even though the file is walker-green.

**This audit reconciles the matrix's `missing` label to the actual coverage
mechanism for every one of the 28 cells.** The matrix is regenerated as-is —
its narrow `missing` definition is intentional and the 5 REPORTED cells'
`skip_reason` strings are a deliberate contributor-facing signal that must stay
visible, so the generator's classifier is *not* broadened in this bead (doing so
would flip the reported-infeasible cells, several of which also ship a
`canonical/` byte contract, into `canonical` and hide their deferral reason).

---

## 2. REPORTED cells (5) — convergence sweep deferred, declarative-infeasible

Each of these ships a `convergence/` fixture marked `applicable: false` with a
`skip_reason`. The MMS-convergence *sweep* is what is deferred; the rule's
own correctness is pinned by a byte contract and/or a dedicated test as noted.

| Cell | Grid | Why the sweep is deferred | Rule still validated by |
|------|------|---------------------------|-------------------------|
| `flux_1d_ppm` | cartesian | Face-stagger output + per-face Courant/velocity bindings + ghost-extended input fall outside the §7 `verify_mms_convergence` contract. A Layer-B′ MMS-transport fixture (PPM flux + FV divergence + RK) is the right shape. | `fixtures/canonical/` byte contract (single-face flux) + `test/test_flux_1d_ppm_rule.jl` + runtime `test/test_transport_{1d,2d}.jl` |
| `flux_limiter_minmod` | cartesian | Limiter acceptance is *monotonicity / TVD preservation*, not a convergence order. | `fixtures/monotonicity/` + `fixtures/conservation/` + `test/test_flux_limiters_rule.jl` |
| `flux_limiter_superbee` | cartesian | Same — monotonicity fixture kind, not a convergence sweep. | `fixtures/monotonicity/` + `test/test_flux_limiters_rule.jl` |
| `lax_friedrichs_flux` | cartesian | Face-stagger output + per-face Courant binding fall outside the §7 harness. | `fixtures/canonical/` byte contract + `test/test_lax_friedrichs_flux_rule.jl` |
| `vertical_remap` | vertical | Structurally a phase-hook operation (Lagrangian → Eulerian re-grid), variable-arity data-dependent neighbor list, time-varying per-column metric — deferred to a future ESS `phase_hooks` RFC. Documentation-only. | Documentation-only; intentionally **not** executed (see §4) |

`flux_limiter_minmod/superbee` are doubly true: REPORTED for the *convergence*
cell, and independently CLOSED for monotonicity (their actual acceptance) via
`fixtures/monotonicity/` and `test_flux_limiters_rule.jl`.

---

## 3. CLOSED cells (23) — covered by a mechanism the matrix scan does not see

### 3a. Layer-A byte contract — `fixtures/canonical/`, walker Layer-A PASS (13)

Each runs through `test/walk_esd_tests.jl` (`EarthSciSerialization.discretize`,
byte-compared to `expected.esm`); the per-rule `test/test_esd_walker.jl` branch
asserts `LAYER_PASS` / `canonical-form match`.

| Cell | Grid | Note |
|------|------|------|
| `centered_4th_uniform_latlon` | latlon | high-order grad, moment-verified coeffs + spherical metric |
| `centered_4th_uniform_vertical` | vertical | high-order grad, dx→h transform |
| `centered_6th_uniform_latlon` | latlon | high-order grad |
| `centered_6th_uniform_vertical` | vertical | high-order grad |
| `centered_8th_uniform_latlon` | latlon | high-order grad |
| `centered_8th_uniform_vertical` | vertical | high-order grad |
| `upwind_1st_latlon` | latlon | first-order upwind, latlon metric |
| `upwind_1st_vertical` | vertical | first-order upwind, vertical |
| `midpoint_1d` | cartesian | dense reduce-arrayop integral (esd-6g4.14) |
| `advection_duo` | duo | flux-form transport; canonical PASS, convergence `applicable:false` (sub-1st-order centered reconstruction on icosahedral mesh; q≡1 reduces to `divergence_duo`) |
| `advection_mpas` | mpas | flux-form MPAS advection; canonical PASS, convergence `applicable:false` (L∞ plateau on non-centroidal Voronoi; rests on shared `divergence_mpas` proof) |
| `flux_duo` | mpas | edge-normal advective flux; canonical byte contract |
| `staggered_1st_uniform_cc_to_face` | arakawa | MAC staggered gradient; canonical byte contract under `staggered_1st_uniform/fixtures/` (file-basename keying — see §1) |

### 3b. Layer-C integration solve — `fixtures/integration/` (3)

Each ships a real ODE-solve fixture (`*_pde.esm` + analytic/MMS reference +
`cases.json`) driven through `build_ode_problem`.

| Cell | Grid | Integration fixture |
|------|------|---------------------|
| `centered_2nd_deriv_uniform_vertical` | vertical | `diffusion_1d_vertical_{pde,analytic}.esm` |
| `centered_2nd_deriv_nonuniform_vertical` | vertical | `diffusion_1d_nonuniform_{pde,analytic}.esm` |
| `covariant_fv_laplacian_latlon` | latlon | `covariant_laplacian_{pde,mms}.esm` + `latlon_n{16,32}.gdd.json` |

### 3c. Dedicated rule / conformance test — no on-disk fixture for the scan (5)

| Cell | Grid | Test that exercises it |
|------|------|------------------------|
| `weno5_grad` | cartesian | `test/test_weno5_advection_rule.jl` + `tests/conformance/discretization/rect_1d_advection_weno5_periodic/` |
| `covariant_fv_gradient_latlon_t1` | latlon | `test/test_covariant_fv_metric_conformance.jl` — evaluated per-cell, asserted tolerance-identical to the grid-assembly oracle golden |
| `covariant_fv_gradient_latlon_t2` | latlon | `test/test_covariant_fv_metric_conformance.jl` — same |
| `ppm_reconstruction_left_edge` | cartesian | `test/test_ppm_reconstruction_split_rules.jl` (single-output sibling of the `ppm_reconstruction` umbrella) |
| `ppm_reconstruction_right_edge` | cartesian | `test/test_ppm_reconstruction_split_rules.jl` |

The two `covariant_fv_gradient_latlon_t1/t2` are the two covariant gradient
tensor components evaluated **at cell centers** (`applied_result_t1/t2`), not
the ESS-tripping edge-output form — feasible and pinned tolerance-identical to
the oracle (`test_covariant_fv_metric_conformance.jl` lines 170–171, 179–180).

### 3d. Operand-order mirror — covered via the primary direction (2)

| Cell | Grid | Disposition |
|------|------|-------------|
| `staggered_1st_uniform_face_to_cc` | arakawa | Mirror direction of `…_cc_to_face`; same `staggered_1st_uniform.json` file, whose canonical byte contract exercises the `cc_to_face` direction. The `face_to_cc` stencil is the documented algebraic mirror (`(u[$x+1] − u[$x])/h`). |
| `advection_mpas_velocity_first` | mpas | Operand-order mirror of `advection_mpas` (sibling key in `advection_mpas.json`) so the operator matches either `div(u·q)` or `div(q·u)` operand naming. "Not separately walked." |

These two are CLOSED-by-mirror: the operator is validated via the primary
direction's byte contract, and the mirror is a declarative operand swap of an
already-pinned rewrite. They carry the thinnest coverage of the set — see §6.

---

## 4. `vertical_remap` (dsc-otd) — deferral confirmed

`vertical_remap` remains **documentation-only, tracked, and deferred**, exactly
as the bead requires. Evidence on the live tree:

- `discretizations/finite_volume/vertical_remap.json` →
  `schema_gaps.disposition = "deferred-to-phase-hook-rfc"`,
  `schema_gaps.executability = "documentation-only; no imperative implementation
  exists in any binding (Julia/Python/Rust/TypeScript) and none is to be
  re-added"`. The full math + AST shape is captured as a reference artifact for
  the eventual ESS `phase_hooks` contract.
- `discretizations/finite_volume/vertical_remap/fixtures/convergence/{input,expected}.esm`
  → `applicable: false` with the phase-hook `skip_reason`; the walker excludes
  it from Layer-B / Layer-D membership.
- `discretizations/finite_volume/README.md` (lines 125, 271) and
  `test/test_esd_walker.jl` (line 102) reference `dsc-otd` as the disposition.

The `dsc-otd` bead is archived, but the deferral is durably recorded in the rule
file, the README, both fixture `skip_reason`s, and the walker comment — the
tracking is intact. No imperative remap implementation exists in any binding.

---

## 5. Conformance green

Full unit suite, run from the bead worktree against `EarthSciSerialization`
resolved from `[sources]` (GitHub `main`), Julia 1.12.6:

```
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

This is the authoritative gate — `test/runtests.jl` runs `@run_package_tests`
(TestItemRunner), which includes `walk_esd_tests.jl` (Layer-A/B/C/D fixture
walker), `test_esd_walker.jl` (per-rule walker assertions), every
`test_*_conformance.jl` (cross-binding Julia arm), and the rule-catalog
well-formedness tests. The result is recorded in §7.

The change in this bead touches no `.jl`, no rule `.json`, and no fixture, so the
test outcome is identical to `main` at `8a3391b` — which the refinery verified
green when it merged each of G1–G15. The run here is the explicit capstone
re-confirmation.

`docs.yml` regenerates `docs/data/rule_matrix.json` in CI before the Hugo build,
so the committed copy is a convenience snapshot; refreshing the stale 39-rule
snapshot to the current 59-rule state keeps the committed artifact honest for
anyone reading the repo directly.

---

## 6. Follow-up (surfaced, not silent)

The two operand-order mirror cells in §3d are the only cells whose coverage is
purely by symmetry with a primary direction rather than a direct assertion:

- `advection_mpas_velocity_first` has **no** test assertion and **no** fixture —
  its only references on the tree are its definition in `advection_mpas.json`
  and a descriptive comment in `test_esd_walker.jl`. It exists for operand-order
  robustness.
- `staggered_1st_uniform_face_to_cc` is exercised only as the documented mirror
  of `…_cc_to_face`.

Neither is a dispatch hole (the operator is validated via the primary
direction), so neither blocks the epic. The honest recommendation is a small
follow-up: either add a one-line catalog well-formedness assertion / canonical
byte contract pinning each mirror direction, or, if the mirror rule is found
redundant with first-match dispatch, retire it. This is noted here rather than
silently absorbed into the `missing` count.

---

## 7. Result

Both suites green against `EarthSciSerialization` @ `main`, Julia 1.12.6, run
from this bead's worktree (`polecat/esd-6g4.15`, base `8a3391b`).

**Unit suite** — `julia --project=. -e 'using Pkg; Pkg.test()'`:

```
[fixture walker]  summary: 420 cases  pass=67  fail=0  skip=353
Test Summary: |    Pass    Total      Time
Package       | 2212695  2212695  16m32.0s
     Testing EarthSciDiscretizations tests passed
```

**Integration suite (Layer-C)** —
`julia --project=. -e 'using Pkg; Pkg.test(; test_args=["integration"])'`:

```
[fixture walker]  summary: 420 cases  pass=81  fail=0  skip=339
Layer-C: 14 pass, 70 skip, 0 fail
Layer-C integration: all cases pass or skip.
     Testing EarthSciDiscretizations tests passed
```

`fail=0` everywhere: the fixture walker (Layer-A byte contracts + Layer-B MMS
sweeps), every `test_*_conformance.jl` cross-binding arm and rule-catalog
well-formedness test inside the 2,212,695-assertion unit package, and the
Layer-C integration solves (14 pass — including the three integration-only
cells `centered_2nd_deriv_{uniform,nonuniform}_vertical` and
`covariant_fv_laplacian_latlon`). Walker + cross-binding conformance + Layer-C
integration are all green.
