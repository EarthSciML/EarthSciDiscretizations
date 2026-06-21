# dimensional_split (§7.5) & flux_form_semi_lagrangian (§7.7) — DECLARATIVE-OR-FAIL verdict

**Bead:** esd-6g4.6 (G6) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Oracle source pin:** EarthSciSerialization GitHub `main` @
`e3c7907f205f6219ad0abd57feeb4926ea629cf9` (the rev ESD `Project.toml`
`[sources]` resolves: `{rev = "main", subdir = "packages/EarthSciSerialization.jl"}`).
All ESS line numbers below are at that commit. ESS spec/RFC line numbers are at
the same commit (`docs/content/rfcs/discretization.md`, `esm-spec.md`).

---

## 0. Verdict (TL;DR)

> **DECLARATIVE — INFEASIBLE over the existing ESS engine.** Neither
> `dimensional_split` (§7.5) nor `flux_form_semi_lagrangian` (§7.7) can be
> authored as a declarative rule that the **existing** ESS engine *evaluates*.
> Per the bead's `DECLARATIVE-OR-FAIL` contract, the correct outcome is this
> precise gap report + escalation — **no rule JSON is committed and no
> imperative executor is hand-coded.**

Both kinds are **spec'd and schema-valid** (they round-trip losslessly in all
bindings), but the ESS rule engine **defers** them — it parses neither into an
evaluable object and has **no operator-sequencing / time-split / sweep
vocabulary** to express what they mean. The ESS RFC itself declares their
expansion/execution **out of scope for the engine** and delegates it to "a
runtime that owns the sweep loop (e.g. ESD)" — i.e. exactly the imperative
executor this bead forbids ESD from hand-coding. That collision is the gap.

This is **not** the same as a missing arrayop footprint (contrast esd-zk9.1's
covariant-FV verdict, which was FEASIBLE because a diagonal corner is just a
two-arg `op:index`). The missing capability here is **temporal**: an ordered,
stateful multi-sweep over a timestep Δt. That is an integration concern, and it
has no spatial-AST encoding.

---

## 1. What the bead asked

> "These scheme kinds are spec'd but have 0 authored rules and the engine defers
> them (scheme_expansion else:continue). Author a real declarative rule using
> each kind (e.g. Lin-Rood dim-split; CAM5 FFSL). If the engine cannot evaluate
> the kind declaratively (needs an executor it lacks), STOP and report … NO
> imperative/hand-coded operator."

The operative test (verbatim): *the rule MUST be declarative JSON
(arrayop/stencil composed over the **EXISTING** ESS engine vocabulary)*. Two
candidate declarative paths exist; **both fail this test**, for independent
reasons (§4, §5). The one thing that *is* declaratively expressible (§6) is the
single-shot **spatial** dimension-by-dimension flux RHS — which is **already in
the catalog** (`weno5_advection.json`) and is *not* what §7.5/§7.7 add.

---

## 2. What §7.5 / §7.7 are, per the spec

Both are **scheme kinds** discriminated by a top-level `"kind"` field, distinct
from the `"stencil"` / `"multi_output_stencil"` kinds the engine handles.

### §7.5 `dimensional_split` (`discretization.md:1482-1562`)
A composite that applies an `inner_rule` as "a **sequence of 1D applications
along orthogonal axes**" (`discretization.md:1486`). Canonical shape
(`discretization.md:1490-1496`):

```json
{ "kind": "dimensional_split", "axes": ["x", "y"],
  "inner_rule": "centered_2nd_uniform", "splitting": "strang",
  "order_of_sweeps": "alternating" }
```

The **splitting semantics are explicitly temporal** (`discretization.md:1515-1526`):

> - `lie`: apply the inner rule along `axes[0]` for time Δt, then `axes[1]` for Δt, …
> - `strang`: apply … along `axes[0]` for Δt/2, then `axes[1]` for Δt/2, … then
>   `axes[N−1]` for Δt, then `axes[N−2]` for Δt/2, … back to `axes[0]` for Δt/2.

`order_of_sweeps` is "per-timestep axis traversal direction" (`:1511`). These are
**operator-splitting time-integration orderings**, not spatial stencils.

### §7.7 `flux_form_semi_lagrangian` (`discretization.md:1695-1773`)
"A sub-cell reconstruction is combined with flux-form remapping so that advection
conserves the volume-integrated tracer mass" (`discretization.md:1699-1700`).
Canonical shape (`discretization.md:1707-1713`):

```json
{ "kind": "flux_form_semi_lagrangian",
  "applies_to": { "op": "adv", "args": ["$u"], "dim": "$x" },
  "grid_family": "cartesian", "reconstruction": { "order": "PPM" },
  "remap": { "semantics": "conservative", "flux_form": "lin_rood_1996" },
  "cfl_policy": "conservative", "dimensions": ["x"] }
```

CAM5/Lin-Rood is assembled by nesting a 1D FFSL rule as the `inner_rule` of a 2D
`dimensional_split` (`discretization.md:1733-1763`).

---

## 3. The engine defers both kinds (the `else: continue`)

`packages/EarthSciSerialization.jl/src/scheme_expansion.jl:333-341` — the entire
scheme-kind dispatch:

```julia
kind_raw = _is_dict_like(v) ? _getkey(v, "kind"; default=nothing) : nothing
kind = kind_raw === nothing ? "stencil" : String(kind_raw)
if kind == "multi_output_stencil"
    out[sk] = parse_multi_output_stencil_scheme(sk, v)
elseif kind == "stencil"
    out[sk] = parse_scheme(sk, v)
else
    # Defer other kinds (cross_metric, dimensional_split, grid_dispatch, …).
    continue          # <-- dimensional_split & flux_form_semi_lagrangian land here
end
```

- The engine handles **exactly two** kinds: `"stencil"` (default) and
  `"multi_output_stencil"`. Everything else hits `else: continue` and is
  **silently dropped from the scheme registry** (`out`). `dimensional_split` is
  named in the deferral comment (`:340`); `flux_form_semi_lagrangian` is not even
  named but falls into the same arm.
- There is **no Julia type, no parse path, and no expand path** for either kind.
  Greps over `src/*.jl` at the pinned commit for
  `DimensionalSplit|FluxFormSemi|FFSL|expand_dimensional|expand_ffsl|parse_dimensional|parse_ffsl`
  return nothing. (`flux_form`, `cfl_policy`, `lin_rood` appear nowhere in `src/`;
  `reconstruction`/`remap` appear only in unrelated docstrings and the separate
  `conservative_regrid.jl` grid-to-grid feature.)
- Consequence for a rule that *uses* such a scheme: because the scheme never
  enters `out`, a downstream `use:` reference fails at `rule_engine.jl`
  `E_SCHEME_MISMATCH` ("references a scheme not declared in `discretizations`").
  So a `kind: dimensional_split` rule is not merely un-expanded — it is
  **un-referenceable**.

**The spec confirms this is by design**, not an oversight
(`discretization.md:1498-1501`):

> "The rule engine is **not required to materialize** dimensional-split schemes
> into a single stencil body at loader time; runtimes that **own the sweep loop
> (e.g. ESD)** consume the declaration directly and orchestrate the 1D
> applications themselves."

and (`discretization.md:1556-1559`):

> "All five bindings … parse and round-trip `kind: "dimensional_split"` schemes
> losslessly. **Structured expansion semantics (inlining the N·Δt/2 sweep into a
> single AST) is out of scope for v0.2.0 and tracked separately in ESD.**"

and for FFSL (`discretization.md:1716-1718`, `:1768-1770`):

> "v0.2.0 treats parse + round-trip + structural rejection as the binding
> contract; concrete runtime dispatch (**execution of the reconstruction / remap
> pair**) lives in ESD." … "Execution of the reconstruction / remap pair is
> **out of scope** for v0.2.0 and is tracked in ESD (`cam5_fv_ffsl_advection`, P1)."

---

## 4. Attempt A — author the rules with their native `kind` (FAILS)

Author `cam5_ffsl_1d` (§7.7) + `cam5_fv_ffsl_advection` / `lin_rood_advection`
(§7.5) exactly as the RFC worked sketches show (`discretization.md:1528-1554`,
`1733-1763`). This produces **schema-valid** JSON (the discriminator + conditional
`required`/`not.required` contracts in `esm-schema.json` accept it).

**Why it fails the bead's test:** "evaluated by the **existing** engine." It is
not. As shown in §3, the engine defers the kind (`else: continue`) — the scheme
is dropped, never expanded, never evaluated. In the Julia walker this rule would
register **0 fixtures it can satisfy**: a Layer-A `canonical` fixture would have
to declare `applicable: false` (the discretize pipeline produces nothing for the
deferred kind), and Layer-B MMS-convergence cannot run (no evaluable operator).
The result is a catalog entry that is **parse-only** — exactly the v0.2.0 binding
contract, and exactly *not* "a real declarative rule … evaluated by the engine"
that the acceptance criteria require.

Committing such a parse-only rule would also contradict the bead's "STOP … do
NOT implement it" once infeasibility is established. It is recorded here as an
*attempt*, not committed as a rule.

---

## 5. Attempt B — inline the whole scheme as one `replacement` arrayop (FAILS)

ESD's house style for multi-axis schemes is to "express their structure inside
the `replacement` itself … rather than via any external scheme-expansion step"
(`discretizations/README.md`), as `weno5_advection.json` does (a 293 KB fully
inlined `div` arrayop). Could a Strang-split / FFSL update be inlined the same
way, bypassing the deferred `kind` entirely? **No** — three independent blockers:

1. **No primitive to bind & stencil-index an intermediate.** A Strang sweep is
   `q^{n+1} = L_x(Δt/2) ∘ L_y(Δt) ∘ L_x(Δt/2) q^n`. Stage 2 must read the
   **array output of stage 1 at neighbor cells** (`±offset`), i.e. stencil-index
   into a *materialized intermediate*. But an `arrayop`'s operands are "the
   **input array operands** that `expr` references" (`esm-spec.md:225`), index
   symbols are node-local with "no cross-node scoping" (`esm-spec.md:214`), and
   the AST has **no `let` / `seq` / temporary-array op** (grep at the pinned
   commit: none). There is no way to name `q_after_x_sweep` as an array and then
   gather `q_after_x_sweep[j±1]` in a later stage. Inlining therefore cannot
   express a multi-stage sweep where a later stage has *stencil* (not pointwise)
   dependence on an earlier stage's output.

2. **Time-update vs. spatial-rewrite mismatch.** Every ESD rule rewrites a
   **spatial** operator — `applies_to: {op: grad|div|adv|laplacian|…}` — into a
   spatial expression (method of lines). §7.5/§7.7 are **full-discrete time
   updates** (`q^{n+1}` from `q^n` over Δt, conservative at CFL > 1). There is no
   `op` for "advance one timestep," and the recent ESD direction is explicitly to
   keep time-stepping *out* of the rule layer (`esd-837`: "make the ODE solver
   test-only — expose the RHS, not a runner"). The Δt/2-then-Δt-then-Δt/2
   ordering, `order_of_sweeps`, and `cfl_policy` are integrator state, not
   stencil coefficients.

3. **The spec forecloses it.** Inlining "the N·Δt/2 sweep into a single AST" is
   named and declared **out of scope** (`discretization.md:1558-1559`). Doing it
   anyway would be inventing an unsanctioned encoding, not "composing over the
   existing vocabulary."

---

## 6. What *is* already covered (so the gap is exactly the new capability)

The **single-shot spatial** dimension-by-dimension flux divergence —
`adv(q,U) = (F^x_E − F^x_W)/dx + (F^y_N − F^y_S)/dy` with per-axis 1D
reconstruction — **is** declaratively expressible and **already shipped**:

- `discretizations/finite_volume/weno5_advection.json` — a 2D
  `dimension_by_dimension` advection inlined as one `div` arrayop.
- Per-axis 1D kernels exist as building blocks: `flux_1d_ppm.json`,
  `ppm_reconstruction.json`, `divergence_arakawa_c.json`.

These cover the *spatial RHS*. What §7.5/§7.7 add on top is **only** the temporal
orchestration (operator-split ordering over Δt; conservative semi-Lagrangian
remap valid for CFL > 1). That orchestration is precisely the infeasible part —
confirming the gap is the genuinely-new capability, not a re-statement of
existing coverage.

---

## 7. The precise gap (for escalation)

Resolving §7.5/§7.7 to *engine-evaluated* rules requires **one** of the
following, none of which a polecat may take unilaterally:

- **(ESS engine front-door)** Add a declarative **operator-sequencing / time-
  split primitive** to the ESS AST + a `scheme_expansion` arm that expands
  `dimensional_split` / `flux_form_semi_lagrangian` (e.g. a `seq`/`sweep` op that
  binds an intermediate array between stages and admits stencil indexing into it,
  plus a Δt-fraction parameter). This is the path that keeps ESD a pure
  passthrough and satisfies "evaluated by the existing engine." **It is new ESS
  engine work** — it does not exist at `e3c7907f`.
- **(Relax the bead constraint)** Accept the RFC's stated division of labor —
  ESD "owns the sweep loop" and executes the split/FFSL imperatively
  (`discretization.md:1500-1501`, `1717-1718`). This **directly contradicts**
  this bead's "do NOT hand-code an executor / NO imperative operator" and ESD's
  single-pathway rule (`AGENTS.md`: ESD carries no evaluator). Requires a
  deliberate architecture decision.
- **(Accept parse-only)** Treat both kinds as catalog entries at the v0.2.0
  binding contract (parse + round-trip + structural rejection only, never
  engine-evaluated). Honest, but does **not** meet this bead's "evaluated by the
  engine" acceptance, so the bead's correct disposition remains *reported*.

**This is an architecture decision spanning ESS and ESD — it belongs to the
mayor** (who owns epic esd-6g4 and coordinates both rigs). The collection point
already exists: `esd-6g4.15` (G16) — "coverage drain — confirm closed/**reported**".

No existing ESS or ESD bead tracks this specific engine-front-door gap (esd-zk9
covers the *different* §7.6 cross_metric kind, and resolved FEASIBLE). If the
mayor chooses the front-door path, file it against ESS as a new engine feature.

---

## 8. Reproduction

```bash
# Authoritative oracle (ESD Project.toml [sources] rev = "main"):
git clone --depth 1 -b main https://github.com/EarthSciML/EarthSciSerialization.git
cd EarthSciSerialization && git rev-parse HEAD     # e3c7907f205f6219ad0abd57feeb4926ea629cf9
S=packages/EarthSciSerialization.jl/src

# 1. The deferral (both kinds → else: continue):
sed -n '333,341p' $S/scheme_expansion.jl

# 2. No expander/parser/type for either kind:
grep -rnE 'DimensionalSplit|FluxFormSemi|FFSL|expand_dimensional|expand_ffsl|parse_dimensional|parse_ffsl|cfl_policy|lin_rood' $S   # (no matches)

# 3. No sequencing/time-split AST op:
grep -rnE '"op": *"(let|seq)"|operator_split|time_split|strang' $S          # (no matches)

# 4. Spec declares expansion/execution out of scope, delegated to ESD runtime:
R=docs/content/rfcs/discretization.md
sed -n '1498,1501p;1556,1559p;1716,1718p;1768,1770p' $R
```

---

## 9. Disposition

- **Verdict:** DECLARATIVE — INFEASIBLE over the existing ESS engine (§0).
- **Action taken:** report only. **No rule JSON, no fixtures, no executor**
  committed — per `DECLARATIVE-OR-FAIL` / "STOP and report … NO imperative
  operator."
- **Escalation:** architecture decision routed to the mayor via the epic's G16
  coverage-drain bead (`esd-6g4.15`); precise front-door change specified in §7.
