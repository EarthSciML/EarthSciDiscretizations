# Arakawa B / D / E stagger operators — DECLARATIVE-OR-FAIL verdict

**Bead:** esd-6g4.13 (G14) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Oracle source pin:** EarthSciSerialization GitHub `main` @
`8ba047ff2590ff621e0cc8dac4f6b049a691d161` (the rev ESD `Project.toml`
`[sources]` resolves: `{rev = "main", subdir = "packages/EarthSciSerialization.jl"}`).
ESD line numbers below are at the esd-6g4.13 worktree HEAD.

---

## 0. Verdict (TL;DR)

> **No new B/D/E operator rule is committed.** The blocker is **not** the
> declarative vocabulary and **not** a schema gate — it is the **validation
> harness**. The arrayop `+ - * / index` vocabulary can *write* every B/D/E
> divergence/gradient AST (corner/edge averaging is just `0.5·(index+index)`),
> and ESS accepts any `stagger` string. But the only Arakawa Layer-B/Layer-D
> runners hardcode the **C-grid face convention** (`face_x`/`face_y`) and reject
> every other stagger, so a B/D/E rule could not be exercised end-to-end — it
> would land as an unverifiable `applicable:false` stub. Per `DECLARATIVE-OR-FAIL`
> (no imperative/unvalidated operator), the correct outcome is this gap report.
>
> Per-grid: **A** — covered by the unstaggered cartesian/latlon operators.
> **C** — done (`divergence_arakawa_c`, `staggered_1st_uniform`). **D** — its
> operators are either byte-identical to C (flux-divergence ⇒ duplication) or
> need cross-face averaging with no validated runner. **B / E** — INFEASIBLE to
> validate (no Corner→CellCenter runner; E additionally has no rotated-axis
> metric).

---

## 1. Grid runtime supports all five staggers; only the operator + runner layer is C-only

`src/grids/arakawa.jl:23` defines `@enum ArakawaStagger ArakawaA ArakawaB
ArakawaC ArakawaD ArakawaE` and `arakawa_variable_locations`
(`src/grids/arakawa.jl:29-34`) gives the (h, u, v) placement:

| stagger | h | u | v | note |
|---|---|---|---|---|
| A | CellCenter | CellCenter | CellCenter | unstaggered |
| B | CellCenter | **Corner** | **Corner** | velocities at vertices |
| C | CellCenter | UEdge (x-face) | VEdge (y-face) | the implemented grid |
| D | CellCenter | **VEdge (y-face)** | **UEdge (x-face)** | C with u/v transposed |
| E | CellCenter | **Corner** | **Corner** | B rotated 45° (`rotated=true`, `arakawa.jl:283`) |

`VarLocation` (`src/staggering.jl:5`) and the Corner accessors are all
first-class, so the *geometry* of B/D/E is fully modeled. The gap is one layer up.

## 2. `stagger` is opaque to ESS — no schema barrier to writing non-C rules

The ESS rule engine's `parse_rule` consumes only `pattern`,
`replacement`/`use`, `where`, `region`, `boundary_policy`, `ghost_width`,
`bindings` — it does **not** read `stagger`, `grid_family`,
`requires_locations`, or `emits_location` (these are ESD-side metadata). So a
rule carrying `stagger:"B"`/`"D"`/`"E"` is parsed exactly like `"C"`; there is
no closed-set rejection of a non-C stagger anywhere in the engine. **Writing**
a B/D/E rule is therefore not blocked. (This mirrors the staggered_1st_uniform
finding in this same bead: the operator AST drives through `discretize`
unchanged.)

## 3. The blocker: the Arakawa runners are C-grid-only

Both arakawa runners hardcode the C-grid face-sampling convention and reject any
other stagger:

- **Layer-B** `_run_layer_b_2d_arakawa_periodic` (`test/walk_esd_tests.jl:1418`)
  samples `Fu[i,j]` at the **face_x** western face and `Fv[i,j]` at the
  **face_y** southern face (`walk_esd_tests.jl:1410-1411`). There is no Corner
  (vertex) sampling and no transposed-face path.
- **Layer-D** `_run_conservation_divergence_2d_periodic`
  (`test/walk_esd_tests.jl:2208`) evaluates a flux read only for
  `stagger == "face_x" && axis == "$x"` or `stagger == "face_y" && axis ==
  "$y"`, and **throws** `"unsupported (stagger=…, axis=…)"` otherwise
  (`walk_esd_tests.jl:2255-2262`).

So any rule whose operands are not C-grid face fluxes has **no runner** to
measure its convergence (Layer-B) or conservation (Layer-D). It could only ship
a Layer-A byte contract plus an `applicable:false` convergence stub — i.e. an
unvalidated operator, which `DECLARATIVE-OR-FAIL` forbids committing.

## 4. Per-grid analysis

### A — covered, no new rule
All variables at CellCenter. The operators are the ordinary unstaggered central
differences already in the catalog (`centered_2nd_uniform`, `…_latlon`, etc.).
Nothing arakawa-specific to author.

### C — done
`discretizations/finite_volume/divergence_arakawa_c.json` (div, face fluxes →
cell center) and `discretizations/finite_difference/staggered_1st_uniform.json`
(grad, cc ↔ face_x; Layer-A activated in esd-6g4.13).

### D — no separate rule warranted
D is "C transposed": the x-velocity `u` is stored on **y-faces** (VEdge) and `v`
on **x-faces** (UEdge) (`arakawa.jl:33`). Two readings, neither yielding a
useful new green rule:

1. **As a flux-divergence** (input = the x-face flux and y-face flux, abstractly
   located): the divergence is the identical 2-point telescoping stencil as C —
   a `divergence_arakawa_d.json` would be **byte-identical** to
   `divergence_arakawa_c.json` (ESS ignores `stagger`), i.e. pure duplication.
2. **As a velocity-divergence** (input = the D-stored `u`,`v`): `∂u/∂x` at the
   cell centre must average `u` across the cell's two **y-faces** before
   differencing in x (and symmetrically for `∂v/∂y`) — a 4-point
   average-and-difference per component. The AST is writable
   (`0.5·(index+index)`), but the Layer-B runner samples on the C-grid faces
   (`walk_esd_tests.jl:1410-1411`), not the transposed faces, so it cannot
   validate this form.

Either way, no separately-validated D operator rule is committed.

### B — INFEASIBLE to validate
`u`,`v` at **Corner** (`arakawa.jl:31`). The standard B-grid divergence averages
the two corners spanning each cell face before differencing
(`div = (ū_E−ū_W)/dx + (v̄_N−v̄_S)/dy`, `ū_E = ½(u[i+1,j]+u[i+1,j+1])`, …). The
fully-expanded form is a fixed linear combination of four corner reads with
constant coefficients — **writable** over the existing vocabulary — but it reads
operands from **Corner** storage and emits at **CellCenter**, a reduction the
runners do not implement: the only stagger values handled are `face_x`/`face_y`
(`walk_esd_tests.jl:2255-2262`). No Corner→CellCenter (vertex) sampling path
exists. Missing primitive to escalate: **a vertex-stagger (Corner→CellCenter)
Layer-B/Layer-D runner**.

### E — INFEASIBLE (strictly stronger than B)
E is B rotated 45° (`rotated=true`, `arakawa.jl:283-284`). It inherits B's
missing Corner→CellCenter runner, **and** its operator lives on the rotated
axes, so its coefficients need a √2-scaled diagonal spacing / rotated-frame
metric that does not exist: `metric_eval` exposes only axis-aligned `:dx`,
`:dy`, `:area` (`arakawa.jl:249-262`), and the rotation is carried solely as a
bare boolean flag (`arakawa.jl:283`). So an E operator cannot even bind the
correct geometric coefficients declaratively today. Missing primitives:
**(a) the Corner→CellCenter runner and (b) a rotated-axis metric / index-set.**

## 5. Bottom line for escalation

To turn B/D/E from "reported" into "green", the engine/harness needs:
1. a **vertex-stagger (Corner→CellCenter) Arakawa runner** (unblocks B, and the
   topological half of E),
2. a **transposed-face (D-grid) sampling path** in the arakawa runner (unblocks
   D's velocity-divergence form), and
3. a **rotated-axis metric / index-set** for the 45° E-grid (unblocks E's
   coefficients).

Until these land, B/D/E operators are **reported infeasible-to-validate**, not
implemented — consistent with the `DIMENSIONAL_SPLIT_FFSL_INFEASIBILITY.md` and
`../latlon/REDUCED_GAUSSIAN_INFEASIBILITY.md` precedents (no rule JSON committed).
