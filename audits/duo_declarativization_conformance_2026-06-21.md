# DUO declarativization — conformance & integration audit (D5)

**Date:** 2026-06-21
**Bead:** esd-heg.10 (D5), capstone of epic `esd-heg` (DUO/unstructured-grid full
declarativization — retire imperative grid construction).
**Scope:** Adversarial verification only. This audit records the evidence; the
only code artifact is a new level-1 cross-language conformance fixture (see §6).

---

## Executive summary

The DUO icosahedral-triangular grid is now constructed **entirely** by declarative
semiring-FAQ rules evaluated through the landed ESS M3 engine. D4b (esd-heg.9)
deleted the imperative builders; this audit confirms the four D5 acceptance
criteria hold against the live code:

1. **Byte-identical grid across Julia / Rust / Python (+ TypeScript) at
   `icos_level0/1/2`** — verified. ✓
2. **`nn_diffusion_duo` solve + convergence (`duo_level2/3/4`)** — verified. ✓
3. **MPAS regression green (shared unstructured pipeline)** — verified. ✓
4. **Zero imperative DUO grid-construction code remains** (grep gate) — verified. ✓

**Verdict: D5 acceptance met.** One coverage gap was found and closed: the
cross-language conformance harness pinned only level 0 and level 2; the bead
names `level0/1/2`. A `medium` (level 1) fixture was added and all four bindings
re-verified against it (§6).

---

## 1. Conformance architecture (what "byte-identical" means here)

DUO correctness is pinned by **two complementary layers**. The distinction
matters for reading the "byte-identical" claim precisely:

### 1a. M3 determinism contract — bit-exact, Julia-reference
The FAQ output is pinned bit-for-bit against canonical goldens produced by the
landed ESS M3 relational engine:

| Contract | Golden | Levels | Assertion |
|----------|--------|--------|-----------|
| Primal topology (value-invention) | `tests/conformance/grids/duo/topology/golden.json` | 0/1/2 | serialized-string `==` (`test_duo_topology_faq.jl`) |
| Voronoi-dual topology | `tests/conformance/grids/duo/voronoi_dual_topology/golden.json` | 1/2/3 | serialized-string `==` |
| Primal geometry | `discretizations/grids/duo/faq/fixtures/canonical/primal_geometry_level0.json` | 0/1/2 | `reinterpret(Int64,·)` 0-ULP (`test_duo_primal_geometry_faq.jl`) |
| Subdivision / level-fold | `…/canonical/mesh_level{0,1,2}.json` | 0/1/2 | 0-ULP |
| Edge + dual geometry | `…/canonical/edge_dual_geometry_level1.json` | 1/2/3 | strict Float64 bit identity incl. ±0.0 |

The **integer topology** (connectivity, edge sets, dedup/rank/inversion-join) is
produced by `EarthSciSerialization.Relational` and is therefore byte-identical
across bindings *by construction* — this is ESS's cross-binding determinism
contract (CONFORMANCE_SPEC §5.5/§5.7). The geometry FAQ reproduces the float ops
exactly (squares as `x*x`, ordered folds, no clipping), so within a binding it is
0-ULP against the reference.

### 1b. Cross-binding accessor conformance — numerically equivalent
Across **Julia / Python / Rust / TypeScript**, geometry cannot be 0-ULP: `acos`,
`atan`, `asin` in the L'Huilier area differ by a few ULP between platform libms.
Per the 2026-04-20 mayor correction on bead `dsc-gbo`, cross-binding conformance
therefore compares **accessor outputs** (`cell_centers`, `neighbors`,
`metric_eval`) at pinned cells under a per-family relative tolerance `1e-12`
(tight enough that real algorithm drift, which lands >1e-9, is caught; loose
enough to absorb libm variance in `[1e-14, 1e-13]`). Integer fields (neighbors,
counts) are compared with exact `==`.

This is the honest reading of the bead's "byte-identical … across
Julia/Rust/Python": **integer topology is byte-identical; geometry is
ULP-tolerant (1e-12) and 0-ULP within the reference binding.**

---

## 2. Grep gate — zero imperative DUO construction

Searched `src/`, `python/src/`, `rust/src/` for definitions of the retired
imperative builders:

- `_subdivide_icosahedron`, `_icosahedron_vertices`, `_midpoint!`,
  `_build_connectivity`, `_spherical_triangle_area`, `_duo_arc` → **no
  definitions** (survive only in explanatory comments/docstrings).
- `_unstructured_const_arrays`, `_rewrite_unstructured_arrayop!` → **gone as
  code** (comment references only, `src/ode_problem.jl`).
- `_extract_connectivity(::DuoGrid)` → **deleted**; replaced by
  `_grid_primitive_arrays(grid::DuoGrid)` which calls the FAQ evaluators
  (`duo_edge_length_faq` / `duo_face_circumcenters_faq` / `duo_dual_edge_length_faq`).
  Only `_extract_connectivity(::MpasGrid)` remains (MPAS mesh-as-data, in scope
  for MPAS not DUO).

The three production constructors (`build_duo_grid` in `src/grids/duo.jl`,
`python/.../grids/duo.py`, `rust/src/grids/duo.rs`) all delegate to the same FAQ
chain: `duo_subdivide_faq → primal_geometry_faq → circumcenter_geo_faq →
primal_topology_faq`. The residual `for`/`push!` loops in `src/subdivide_faq.jl`,
`src/topology_faq.jl`, `src/primal_geometry_faq.jl` are **ESS-front-door
marshalling** — they iterate cells/faces and hand every arithmetic leaf to
`eval_coeff` and every dedup/rank/join to `Relational`; no DUO geometry or
connectivity is computed by hand. This matches the epic's "single arithmetic
pathway, no shadow evaluator" framing (AGENTS.md / GRIDS_API §4.3).

**Gate: PASS.**

---

## 3. Byte-identity / cross-language conformance — results

- **Julia unit suite** (`Pkg.test()`): **2,210,330 / 2,210,330 pass, 0 fail**
  (14m45s). Covers every DUO FAQ byte/0-ULP `@testitem` in §1a plus the
  cross-language accessor `@testitem`.
- **Determinism re-check:** re-running `regenerate_golden.jl` reproduced
  `small.json` and `realistic.json` **byte-for-byte unchanged** (git shows no
  diff) — the committed goldens are current and the build is deterministic.
- **Cross-binding accessor conformance** (`small` L0, `medium` L1, `realistic`
  L2), all green:
  - Julia: `test_duo_conformance.jl` (in the unit suite) + standalone 80-cell
    check of `medium`.
  - Python: `test_duo_conformance.py` — 3 fixtures pass.
  - Rust: `cargo test --test duo_conformance` — `duo_matches_golden` ok.
  - TypeScript: `npm test -- duo.conformance` — 3 tests pass.

---

## 4. Convergence — `nn_diffusion_duo` solve

Integration suite (`Pkg.test(; test_args=["integration"])`): **400 cases,
pass=73, fail=0, skip=327, exit 0** — no `*=FAIL` markers anywhere in the run.

`[finite_difference/nn_diffusion_duo]  A=SKIP  B=PASS  B'=SKIP  C=PASS  D=SKIP`

**Layer C is the real ODE solve.** `_run_unstructured_mms_convergence`
(`test/integration_cases/esd_field_pipeline.jl`) builds the icosahedral ladder
`duo_level{2,3,4}_unit.gdd.json` (R = 1), seeds the manufactured solution
`u = sin(lon)·cos(lat)` (exact `Δu = −2u/R²`, decay `exp(−2t/R²)`), runs
`build_ode_problem` + `solve` per level, and confirms the minimum observed
L-inf convergence order clears the declared `min_order = 1.2`
(`duo_nn_diffusion_convergence.esm`; documented target `expected_min_order = 1.9`).
`C=PASS` ⇒ order ≥ 1.2 on the `duo_level2/3/4` ladder the bead names.

The unit-suite walker independently reports `nn_diffusion_duo B=PASS` (the
structural Layer-B min-order check on the R = 6.371e6 ladder).

---

## 5. MPAS regression — shared unstructured pipeline

`nn_diffusion_duo` and `nn_diffusion_mpas` dispatch through the **same**
`unstructured_ode` runner (`walk_esd_tests.jl`), and `_run_unstructured_mms_convergence`
(`test/integration_cases/esd_field_pipeline.jl`) branches on `family ∈ {mpas,duo}`
in one function. The DUO→MPAS bridge `_duo_voronoi_dual` (`src/grids/mpas.jl`) is
itself now fully declarative (esd-heg.9): `build_duo_grid` +
`voronoi_dual_topology_faq` + `duo_dual_geometry_faq`, no imperative geometry.

- `test_mpas_conformance.jl` (in the unit suite, 0 fail): MPAS accessor
  conformance + the shared-pipeline rule-coefficient check
  (`dv_edge/(dc_edge·area_c)`, `center = -Σ reduction`).
- Integration suite (same run as §4): `[finite_difference/nn_diffusion_mpas]
  C=PASS` — MPAS rides the same Layer-C ODE-solve convergence runner as DUO and
  passes. `[finite_volume/divergence_mpas] A=PASS B=PASS`,
  `[finite_volume/gradient_mpas] A=PASS B=PASS`, `[finite_volume/advection_mpas]
  A=PASS`.

Both `nn_diffusion_duo` and `nn_diffusion_mpas` reaching `C=PASS` through the
single `family ∈ {mpas,duo}` branch is the direct evidence that DUO
declarativization did not regress the shared unstructured pipeline.

**MPAS regression: green.**

---

## 6. Coverage gap found & closed — level-1 cross-language fixture

The bead names `icos_level0/1/2`. The cross-binding accessor harness
(`tests/conformance/grids/duo/`) pinned only `small` (L0, all 20 cells) and
`realistic` (L2, sampled). Level 1 was covered byte-exactly in the Julia-only M3
goldens (§1a) but **not** in the cross-language harness.

Closed by adding fixture **`medium`** (level 1, R = 1.0 unit sphere, all 80
cells) to `fixtures.json`, regenerating `golden/medium.json` from the Julia
reference, and re-verifying all four bindings (§3). R = 1.0 keeps the metric
magnitudes O(1) so first-subdivision determinism is checked at the same scale as
the L0 baseline, while `realistic` continues to confirm Earth-radius scaling at
depth. Counts: Nc = 80, Nv = 42, Ne = 120 (Euler V−E+F = 2). All bindings
iterate `fixtures.json`, so the fixture extends Julia/Python/Rust/TypeScript
uniformly. README updated.

Level-1 correctness across bindings is independently implied by the
already-green level 2 (L2 = subdivide(L1)); the new fixture makes the contract
explicit at the level the bead names.

---

## 7. Verdict

All four D5 acceptance criteria are met. Imperative DUO grid construction is
fully retired; the grid is byte-identical (integer topology) and ULP-equivalent
(geometry, 1e-12) across Julia/Rust/Python/TypeScript at level 0/1/2;
`nn_diffusion_duo` converges; MPAS rides the same declarative pipeline and
regresses green. The epic `esd-heg` goal — "ZERO imperative DUO
grid-construction code remains" — holds.
