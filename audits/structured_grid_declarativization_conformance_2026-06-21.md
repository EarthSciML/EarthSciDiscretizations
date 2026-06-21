# Structured-grid declarativization — conformance & integration audit (S6)

**Date:** 2026-06-21
**Bead:** esd-3we.6 (S6), capstone of epic `esd-3we` (structured-grid full
declarativization — retire imperative cartesian / vertical / lat_lon / arakawa
construction).
**Scope:** Adversarial verification of the structured-grid campaign. This audit
records the evidence; the only code artifact is a new Julia cross-language
conformance `@testitem` for lat_lon (see §6) that closes a missing harness arm.

---

## Executive summary

The structured-grid families are now constructed in the **Julia binding** by
declarative semiring-FAQ rules evaluated through the landed ESS M3 engine
(`cartesian_construction.esm` / `vertical_construction.esm` /
`latlon_construction.esm` / `arakawa_construction.esm`, materialized by
`src/{cartesian,vertical,latlon,arakawa}_faq.jl`, 1045 LoC). S5 (esd-3we.5)
deleted the imperative Julia coordinate/metric builders and routed every bulk
accessor through the FAQ. This audit checks the three S6 acceptance criteria
against the live code:

1. **Byte-identical structured grids across bindings (conformance)** — verified.
   ✓ The Julia FAQ output reproduces the cross-binding goldens
   (cartesian / vertical / arakawa already pinned; **lat_lon Julia arm added
   here**, §6) at the documented tolerance.
2. **PDE solve/verify regression green** — verified. ✓ Layer-C integration
   suite green (latlon diffusion/advection etc.; §5).
3. **Grep gate — no imperative structured construction remains** — **met in
   Julia only.** ✓ Julia bulk construction is fully FAQ-routed; **✗ python/ and
   rust/ still build coordinates/metrics imperatively** (2121 + 4701 LoC, zero
   FAQ routing). This is the explicitly-deferred tail of S5, tracked by the open
   bead **esd-dru**, not a regression (§2).

**Verdict (honest):** The conformance + integration verification is **green**,
and one real coverage gap was found and closed (the lat_lon Julia cross-binding
conformance arm, §6). The campaign goal — "ALL ESD grid construction
declarative" — holds for the **Julia reference binding** (and for the
unstructured families via the sibling epic `esd-heg`), but is **not yet met
cross-binding**: the python/rust imperative source retirement is pending
`esd-dru`. **Epic `esd-3we` cannot be fully closed until `esd-dru` lands**; S6's
own deliverable (the conformance/integration audit) is complete. See §7.

---

## 1. Conformance architecture (what "byte-identical" means here)

Structured-grid correctness is pinned by **two complementary layers**, the same
two-layer contract the DUO audit (D5) documented:

### 1a. M3 determinism contract — byte-exact, Julia-reference
The Julia FAQ output is pinned bit-for-bit against canonical construction
goldens via `reinterpret(UInt64, ·)` 0-ULP equality:

| Family | Construction test | Golden | Tags |
|--------|-------------------|--------|------|
| cartesian | `test_cartesian_construction_faq.jl` | `tests/conformance/grids/cartesian/construction/golden.json` | `[:conformance,:grid,:cartesian,:construction,:faq]` |
| vertical | `test_vertical_construction_faq.jl` | `…/vertical/construction/golden.json` | `…:vertical:construction:faq` |
| lat_lon | `test_latlon_construction_faq.jl` | `…/latlon/construction/golden.json` | `…:lat_lon:construction:faq` |
| arakawa | `test_arakawa_construction_faq.jl` | `…/arakawa/construction/golden.json` | `…:arakawa:construction:faq` |

Each declarative rule lives binding-neutrally under
`discretizations/grids/<family>/rules/<family>_construction.esm`. Post-S5 the
constructor itself routes through the FAQ, so the test's historical
"FAQ matches imperative builder (ULP)" comparison is now **tautological** (the
honest reading S5's commit gave); the substantive contract is the byte
`golden.json` pin, which is unaffected.

### 1b. Cross-binding accessor conformance — numerically equivalent
Across **Julia / Python / Rust / TypeScript**, each family has a harness at
`tests/conformance/grids/<family>/` with ONE `fixtures.json` + committed
`golden/*.json`, iterated by every binding's conformance test. Structured-grid
geometry is closed-form (affine for cartesian/arakawa/vertical, `sin`/`cos` on
cell-centre latitudes for lat_lon — **no `acos`**, unlike DUO), so the tolerance
is tight: **`relative = 1e-14`** for cartesian / lat_lon / arakawa and
**`0.0` (strict bytes)** for vertical. Integer fields (neighbours, counts) are
compared with exact `==`.

The honest reading of "byte-identical … across Julia/Rust/Python" for structured
grids: **integer topology is exact; geometry is ≤1e-14 across bindings and
0-ULP within the Julia reference** (the construction golden of §1a).

---

## 2. Grep gate — imperative structured construction

Searched `src/`, `python/src/`, `rust/src/` for the structured-grid coordinate /
metric / level builders.

### Julia — retired (S5). **PASS.**
All four `src/grids/{cartesian,vertical,latlon,arakawa}.jl` route their bulk
accessors through `_construction_faq(g)` → `<family>_construction_faq`
(`src/{cartesian,vertical,latlon,arakawa}_faq.jl`):

- `cell_centers(g, ::Symbol)`, `cell_volume(g)`, `metric_g/ginv/jacobian/
  dgij_dxk(g)`, `coord_jacobian(g)` all slice the memoized FAQ NamedTuple; the
  imperative affine/trig/midpoint loops are **deleted** (S5, net −139 LoC). The
  axis synthesis (`_uniform_axis`/`_nonuniform_axis`) delegates to `_faq_axis`.
- Residual host-side code is **pure-structural** (neighbour/boundary index maps)
  — sanctioned by the FAQ convention exactly as in the DUO audit — plus the
  scalar single-cell query accessors `cell_area(g,j,i)` / `metric_eval(g,name,
  j,i)` (`src/grids/latlon.jl`), which evaluate one closed-form value on demand.
  These are query-time accessors, not the bulk builder; they are consistent with
  the FAQ math (verified by the §6 conformance test, which exercises exactly
  these scalar accessors against the cross-binding golden) and were retained by
  S5's documented scope.

### Python / Rust — imperative construction REMAINS. **Cross-binding criterion NOT met.**
`python/src/.../grids/` and `rust/src/grids/` contain **only `duo_*_faq`**; there
is **no** `cartesian_faq` / `vertical_faq` / `latlon_faq` / `arakawa_faq` in
either binding. The structured constructors compute geometry by hand:

| Binding | Imperative structured LoC | FAQ / `eval_coeff` refs |
|---------|---------------------------|--------------------------|
| python  | 2121 (`cartesian.py` 519, `vertical.py` 506, `lat_lon.py` 614, `arakawa.py` 482) | **0** |
| rust    | 4701 (`cartesian.rs` 1169, `vertical.rs` 1353, `lat_lon.rs` 1273, `arakawa.rs` 906) | **0** |

e.g. `python/.../grids/lat_lon.py` computes `dlon = 2.0*math.pi/n` and the
spherical-rectangle `cell_area` with `import math`, no `eval_coeff`. This is the
**explicitly-deferred tail of S5**, whose commit reads: *"python/rust deferred …
blocked on porting the FAQ to those bindings … Filed as esd-dru."* It is the
structural analogue of the DUO epic's W1/W2 (esd-ohd/esd-un6) front-door port
that had to precede the DUO cross-binding deletion D4b. **esd-dru is open.**

**Gate: PASS for Julia (the GRIDS_API §4.3 reference binding); FAIL
cross-binding, by the known, tracked `esd-dru` deferral — not a regression.**

---

## 3. Byte-identity / cross-language conformance — results

- **Julia unit suite** (`Pkg.test()`): **2,212,112 pass / 2,212,112 total /
  0 fail** (17m10s, exit 0), **+908 over the pre-edit baseline of 2,211,204** —
  the delta is the new lat_lon conformance `@testitem` (§6). Covers every
  structured construction-FAQ byte/0-ULP `@testitem` (§1a) plus the
  cross-language accessor `@testitem`s for cartesian / vertical / arakawa **and
  the lat_lon arm added by this bead** (§6).
- **Cross-binding accessor conformance** — the four families' Julia
  `@testitem`s compare the FAQ-backed Julia accessors against the **same**
  committed `golden/*.json` that `python/tests/`, `rust/tests/`,
  `typescript/tests/` load:
  - cartesian — `test_cartesian_conformance.jl` (small_1d/small_2d/realistic_3d/
    nonuniform_2d), 1e-14.
  - vertical — `test_vertical_conformance.jl` (sigma n16/n64, z, eta, theta),
    strict bytes.
  - lat_lon — **`test_latlon_conformance.jl` (NEW, §6)**: small (regular) +
    realistic (regular) + **reduced (reduced_gaussian)**, 1e-14.
  - arakawa — `test_arakawa_conformance.jl` (small/realistic, staggers A–E),
    1e-14.
- **Live Python cross-check (this audit).** `python/tests/` run against the
  shared goldens on this box (`PYTHONPATH=src pytest`):
  - `test_lat_lon_conformance.py` — **3 passed** (small/realistic/reduced).
  - `test_vertical_conformance.py` + `test_arakawa_conformance.py` — together
    with lat_lon, **11 passed**.
  - `test_cartesian_conformance.py` — not run live here: its runner uses
    `zip(..., strict=True)` (a Python ≥3.10 idiom) and this box's interpreter is
    3.9.21; the grid accessors are unaffected and CI (Python 3.11) covers it.

  This is a **live** confirmation that the Python imperative accessors ≡ shared
  golden for lat_lon/vertical/arakawa, which — combined with the Julia-FAQ ≡ same
  golden result — makes **Julia ≡ Python** directly verified (not merely
  CI-cited) for three of the four families, and CI-gated for cartesian.
- **rust / typescript** structured-grid source is **byte-unchanged on this
  branch** (the diff is the audit + the Julia `@testitem` + harness docs only),
  so their conformance against the shared goldens is the standing CI-gated green
  (`Rust.yml` / vitest). Full cross-binding chain: **Julia-FAQ ≡ shared golden**
  (verified here, incl. the new lat_lon arm) **· Python ≡ same golden** (live,
  above) **· Rust/TS ≡ same golden** (CI, code unchanged).

---

## 4. Cross-binding parity — the chain in one line

For each structured family the four bindings are pinned to ONE shared
`golden/*.json`. Parity holds transitively:

> **Julia-FAQ ≡ golden** (§3, live, incl. new lat_lon arm) **· Python ≡ golden**
> (§3, live for lat_lon/vertical/arakawa) **· Rust/TS ≡ golden** (CI, source
> byte-unchanged on this branch) ⟹ **all four bindings agree** to ≤1e-14
> (structured geometry is closed-form; integer topology is exact).

The only binding whose accessor output *changed* in this campaign is Julia (S5
rerouted it through the FAQ); that is exactly the arm re-verified live here.

---

## 5. PDE solve/verify — integration regression

Integration suite (`Pkg.test(; test_args=["integration"])`, exit 0):
**400 cases, pass=74, fail=0, skip=326; Layer-C: 13 pass, 67 skip, 0 fail** — no
`*=FAIL` markers anywhere in the run. The S5 reroute touched `ode_problem.jl`
(`_inject_grids_spatial` per-cell widths → `_faq_axis`; `_construct_curvilinear_
grid` returns a FAQ-materialized `LatLonGrid`), so the Layer-C latlon solves are
the direct regression witnesses that the FAQ reroute did not perturb the PDE
pipeline:

```
[finite_difference/centered_2nd_uniform_latlon]  A=PASS  B=PASS  B'=SKIP  C=PASS  D=SKIP
[finite_difference/centered_4th_uniform_latlon]  A=PASS  B=SKIP  ...
[finite_difference/upwind_1st_latlon]            A=PASS  B=SKIP  ...
[finite_volume/weno5_advection]                  A=PASS  B=PASS  B'=SKIP  C=PASS  D=SKIP
[finite_difference/nn_diffusion_duo]             A=SKIP  B=PASS  B'=SKIP  C=PASS  D=SKIP
[finite_difference/nn_diffusion_mpas]            A=SKIP  B=PASS  B'=SKIP  C=PASS  D=SKIP
```

`centered_2nd_uniform_latlon C=PASS` is the live ODE-solve on the FAQ-materialized
LatLonGrid — the latlon diffusion regression the bead names. Higher-order latlon
+ `upwind_1st_latlon` fire their Layer-A/declarative rules (`A=PASS`) on the same
grid; they carry no Layer-C convergence fixture (`C=SKIP`), which is a
pre-existing fixture-coverage state, not a regression. Advection rides
`weno5_advection C=PASS` and the unstructured `advection_mpas A=PASS`. The
unstructured `nn_diffusion_{duo,mpas} C=PASS` confirm the shared pipeline is
unaffected.

**PDE regression: green** (0 fail across 400 cases).

---

## 6. Coverage gap found & closed — lat_lon Julia conformance arm

**The gap.** cartesian, vertical, and arakawa each ship a Julia "cross-language
conformance" `@testitem` that checks the Julia accessor runtime against the
shared cross-binding `golden/*.json`. **lat_lon shipped none.** Its
`golden/{small,realistic,reduced}.json` was generated by the Python binding
(`regenerate_golden.py`) when Python landed first; the Rust and TypeScript
suites load it, but the **Julia** binding — which landed declaratively via the
construction FAQ in S3 (esd-3we.3), after the harness was authored — was never
wired into this layer. `test_latlon_construction_faq.jl` pins the Julia FAQ
output against a *Julia-generated construction* golden (§1a, a Julia-internal
determinism check), but **nothing checked that the Julia lat_lon accessors
reproduce the cross-binding (Python/Rust/TypeScript) golden.** For a bead whose
charge is "byte-identical … across Julia/Python/Rust … lat_lon (regular +
reduced_gaussian)", that arm being absent is a genuine hole.

**The closure.** Added `test/test_latlon_conformance.jl` — a Julia
`@testitem` (tags `[:conformance,:grid,:lat_lon]`) that loads the same
`fixtures.json` + `golden/{small,realistic,reduced}.json`, builds the Julia
`LatLonGrid` from each fixture's opts (regular via nlon/nlat; reduced_gaussian
via `nlon_per_row` + user-supplied `lat_edges` DATA), and compares
`cell_centers`, `neighbors` (W/E/S/N, pole→`nothing`), the seven
`metric_eval` tensor components + Jacobian, and the spherical-rectangle `area`
at every fixture query point, under the harness `1e-14` tolerance. It mirrors
`test_cartesian_conformance.jl` (0-indexed golden → 1-indexed Julia accessors at
the boundary) and the Python arm `python/tests/test_lat_lon_conformance.py`.

All three fixtures pass — small (14 qps), realistic (17 qps), **reduced (20 qps,
the reduced_gaussian case the bead names)** — confirming the Julia FAQ-backed
lat_lon binding is byte-equivalent (≤1e-14) to the Python/Rust/TypeScript golden.
This makes "byte-identical structured grids across Julia/Python/Rust" a
*verified Julia-side contract for all four structured families*, where before it
was unverified for lat_lon.

---

## 7. Verdict

**S6's conformance + integration verification is complete and green**, and a real
lat_lon conformance gap was closed:

- Structured-grid construction in the **Julia reference binding** is fully
  declarative (S1–S5); the FAQ output is byte-pinned (§1a) and cross-binding
  conformant (§1b, §3, §6).
- **PDE solve/verify integration is green** (§5) — the S5 reroute did not
  regress the pipeline.
- The unstructured families (DUO/MPAS) are already declarative via the sibling
  epic `esd-heg` (audit `duo_declarativization_conformance_2026-06-21.md`).

**But the epic's stated goal — declarativize structured construction across all
three bindings (`src/grids/*` + python/ + rust/, the ~9100-LoC charter) — is not
yet met:** python (2121 LoC) and rust (4701 LoC) still build structured grids
imperatively (§2). This is the known, filed `esd-dru` deferral, the structural
analogue of DUO's W1/W2.

**Recommended sequencing:** land `esd-dru` (port the structured FAQ to python +
rust, rewire the constructors, delete the imperative arithmetic — mirror W2
`esd-un6`), then re-run this S6 grep gate to certify "zero imperative structured
construction across all bindings." Only then is epic `esd-3we` closeable.
**S6 itself (this audit + the lat_lon conformance closure) is done.**
