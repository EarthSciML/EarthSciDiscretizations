# Cross-rig periodic-fold regression — 2026-06-22 (esd-2my)

**Status: RESOLVED via Option A (user-approved 2026-06-22). The 7 periodic-grid
Layer-A canonical goldens were regenerated to the new declarative
makearray-region AST; no engine primitives added. Full ESD suite (Tests +
Integration) + Runic green against ESS `main` incl ess-8ne. See "Resolution"
below.** The escalation that produced this doc (STOP-AND-REPORT, per the esd-2my
contract "⚠ if ANY periodic/BC test regresses: STOP and report; do NOT patch
around") is retained for the causal record.

Scope: verify the full ESD suite is green against ESS `main` now that
**ess-8ne** (ESS `d9049910`, "retire imperative periodic fold; periodic wraps
via declarative makearray") is resolved, and add a periodic cross-rig
regression golden.

Resolution confirmed: ESD `[sources]` (`rev = "main"`,
`https://github.com/EarthSciML/EarthSciSerialization.git`) resolves to ESS
package tree `46cea1a8bb19a8de79b134b89e4cf16ef2db62d7`, which maps **exactly**
to commit `d9049910` (ess-8ne). The installed source carries no
`_apply_periodic_folding!` definition (only test comments noting its
retirement); ess-wg0 (`1d98e014`, constant-first reach scan) and ess-hjg
(makearray-region BC lowering) are present.

---

## The regression

`Pkg.test()` (unit suite, Tests.yml path) is **RED**:
`2212898 passed, 14 failed, 0 errored, 0 broken`.

All 14 failures are in `test/test_esd_walker.jl`, testitem *"walker: discovers
seeded rules; per-layer outcomes match the catalog ledger"*, and all are
**Layer-A canonical byte-contract mismatches** — `14 = 7 rules × 2 assertions`
(`r.layer_a.outcome == LAYER_PASS` and `occursin("canonical-form match", …)`).

The 7 affected rules are exactly the rules that discretize on a **periodic
grid**:

| Rule | Family | Grid periodicity | Canonical-form change |
|---|---|---|---|
| `laplacian_2nd_uniform_cartesian` | finite_difference | 2-D cartesian periodic | plain `index(u,i,j)` → `makearray` + `regions` |
| `mixed_deriv_2nd_uniform` | finite_difference | 2-D cartesian periodic | plain stencil → `makearray` + `regions` |
| `centered_2nd_uniform_latlon` | finite_difference | latlon (periodic in lon) | nested ghost → constant-first `index(u, -1+lat, lon)` |
| `centered_4th_uniform_latlon` | finite_difference | latlon (periodic in lon) | nested ghost → constant-first ghost |
| `centered_6th_uniform_latlon` | finite_difference | latlon (periodic in lon) | nested ghost → constant-first ghost |
| `centered_8th_uniform_latlon` | finite_difference | latlon (periodic in lon) | nested ghost → constant-first ghost |
| `upwind_1st_latlon` | finite_difference | latlon (periodic in lon) | nested ghost → constant-first ghost |

The **vertical** high-order rules in the *same* walker set
(`pass_layer_a_canonical_only_g14`: `centered_{4,6,8}th_uniform_vertical`,
`upwind_1st_vertical`) are **non-periodic** and **still pass** — the breakage is
isolated to periodic grids.

Representative diff (from the test log), `laplacian_2nd_uniform_cartesian`:

```
canonical-form mismatch at byte 650 (actual=10575B, expected=3543B)
  actual:   …"rhs":{…"expr":{"args":[{"op":"makearray","regions":[[[1,4],[1,4]],…
  expected: …"rhs":{…"expr":{"args":[{"args":[{"args":["u","i","j"],"op":"index"…
```

`centered_2nd_uniform_latlon`:

```
canonical-form mismatch at byte 914 (actual=2288B, expected=2792B)
  actual:   …{"args":["u",{"args":[-1,"lat"],"op":"+"},"lon"],"op":"index"},…
  expected: …{"args":["u",{"args":[{"args":[{"args":[-1,"lat"],"op":"+"},1],…
```

The committed goldens (`fixtures/canonical/expected.esm`) are **pristine** (no
diff vs `origin/main`) and were authored by `88c0fa0`
("author Layer-A canonical byte contracts for the newly activated rules") — i.e.
against an ESS that **predates** ess-8ne. The laplacian golden contains
`0 makearray / 0 regions / 13 index` nodes: the **pre-fold** stencil form.

---

## Causal chain (why this is ess-8ne, not pre-existing / not env)

1. ESD's last green run (esd-2bo, `7cd1838`, "flip 2-D periodic corner
   @test_broken to @test now ess-wg0 landed") ran against ESS with **ess-wg0
   present and ess-8ne absent** (the fold still in place).
2. The **only** ESS delta between then and now is ess-8ne (`d9049910`).
3. ess-8ne moves periodic wrapping from the imperative `_apply_periodic_folding!`
   (applied *after* AST construction, invisible to the Layer-A canonical AST)
   into the **declarative makearray-region path** (now *visible* in the AST as
   `makearray`/`regions` for 2-D and as the ess-wg0 constant-first ghost index
   for latlon). The Layer-A canonical serialization therefore changes for every
   periodic-grid rule.
4. The change is isolated to periodic grids (vertical unaffected), matching
   ess-8ne's semantics precisely.

So the Layer-A **AST** byte-contracts are no longer byte-identical. ESS verified
**numerical** byte-identity (cartesian decay `exp(-0.4)=0.6703200460356393`,
latlon `λ=-16/π²`), which is a different surface from these AST goldens — which
is why ESS-green hid ESD-red (the ess-e7u → ess-1nm risk the bead names).

---

## Resolution — Option A applied + verified (2026-06-22)

This was **not** a numerical regression — it is the *expected* canonical-AST
consequence of moving periodic from an imperative fold to the declarative
makearray path. The user approved **Option A** (accept + regenerate): the new
declarative-periodic AST is the intended canonical Layer-A form, so the 7
periodic-grid `expected.esm` goldens were re-authored against ess-8ne. This is
a deliberate, independently-verified golden update — **not** a patch-around.

**What was done.** Each of the 7 periodic-grid `fixtures/canonical/expected.esm`
was regenerated through the *same* `EarthSciSerialization.discretize` →
`canonical_doc_json` path the walker test (`test/walk_esd_tests.jl`) uses, so
the bytes are identical to what the test recomputes. The environment resolved
ESS to exactly ess-8ne (`d9049910`; subdir tree `46cea1a8…`, `origin/main`
HEAD). Vertical / non-periodic goldens were **not** touched (`git status`
confirmed only these 7 `expected.esm` changed).

**Form verification (the declarative form, not a degenerate/zero-ghost form).**

| Rule | New bytes | Declarative marker |
|---|---|---|
| `laplacian_2nd_uniform_cartesian` | 10575 | `makearray` + 13-region decomposition (interior `[[1,4],[1,4]]` + 12 boundary regions) |
| `mixed_deriv_2nd_uniform` | 4989 | `makearray` + `regions` |
| `centered_2nd_uniform_latlon` | 2288 | constant-first ghost `index(u, -1+lat, lon)` (±1) |
| `centered_4th_uniform_latlon` | 3452 | constant-first ghost (±1, ±2) |
| `centered_6th_uniform_latlon` | 4616 | constant-first ghost (±1, ±2, ±3) |
| `centered_8th_uniform_latlon` | 5836 | constant-first ghost (±1, ±2, ±3, ±4) |
| `upwind_1st_latlon` | 2154 | constant-first ghost (−1) |

The two byte counts the escalation predicted from the test log (laplacian
10575 B, centered_2nd 2288 B) match exactly, and the centered_2nd golden carries
the exact fragment `{"args":["u",{"args":[-1,"lat"],"op":"+"},"lon"],"op":"index"}`
quoted in the diff above — confirming byte-identity with ess-8ne.

**Suite verification.** Full unit suite (`Pkg.test()`): `2212918 passed,
0 failed` (was `2212898 passed, 14 failed`). The delta reconciles exactly: the
14 walker Layer-A assertions now pass, plus the 6 new assertions from the
regression golden below. Integration suite + Runic format check: green.

**Cross-rig periodic regression golden** (`test/test_periodic_wrap_golden.jl`),
also landed: it builds the 2-D grid-level periodic Laplacian via
`build_ode_problem` and verifies the operator wraps both axes + corners to the
analytic periodic operator (`max_err < 1e-9`); the corner `du[1,1]` is
byte-pinned at `-30.20450580068021` and asserted > 100 away from the `-231`
zero-ghost fallback, so a silent un-wrap regression fails loudly. This closes
the coverage gap (no prior ESD end-to-end test exercised a grid-level periodic
ghost) that let the ess-wg0 2-D reach bug go unnoticed for a refinery cycle.
