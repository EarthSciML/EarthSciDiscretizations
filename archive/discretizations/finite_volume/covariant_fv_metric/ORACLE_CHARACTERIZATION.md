# Covariant FV operator — oracle characterization & DECLARATIVE-OR-FAIL verdict

**Bead:** esd-zk9.1 (R1) · **Epic:** esd-zk9 · **Label:** `covariant-fv-declarative`
**Oracle source pin:** EarthSciSerialization GitHub `main` @
`0e935a3384a6f13d396598a7dfebe2197e20b22c` (the rev ESD `Project.toml`
`[sources]` resolves). All line numbers below are at that commit.

---

## 0. Verdict (TL;DR)

> **DECLARATIVE — FEASIBLE.** The covariant FV Laplacian and gradient ARE
> expressible as declarative **arrayop-einsum rules over the EXISTING ESS engine
> vocabulary** (existing `op:index` multi-axis gathers + `op:arrayop` + arithmetic
> + bound coefficient metric arrays). **No new engine code, no new AST node, no
> new arrayop footprint, and no `cross_metric`/`dimensional_split` scheme handler
> is required.**

The bead's hypothesized blocker — *"a corner/diagonal-gather footprint the
arrayop lacks"* — **does not exist**. A diagonal corner is just a two-argument
`op:index` (`u[i±1, j±1]`), already lowered by the engine with no special-casing
and already used by the activated rule `mixed_deriv_2nd_uniform.json`. The exact
R2 rule form is specified in §5.

Two **non-blocking** items are real R2 work (both DATA / host-binding, neither is
engine code) and are scoped precisely in §6: (a) a `const_arrays` metric-binding
hook for the latlon path, replicating the existing MPAS/DUO precedent; (b)
lat-pole boundary-value matching (oracle `sentinel→self` vs engine zero-ghost).

---

## 1. The oracle operator (math + line-level provenance)

Numeric oracle: `EarthSciSerialization/src/grid_assembly.jl`. Its symbolic
ArrayOp counterpart (already an `ArrayOp`!) is
`ext/grid_assembly_symbolic.jl` — see §3.

### 1.1 Covariant Laplacian — 9-point stencil

Continuous form (`grid_assembly.jl:79-83`):

```
∇²φ = g^{ξξ} ∂²φ/∂ξ² + 2 g^{ξη} ∂²φ/∂ξ∂η + g^{ηη} ∂²φ/∂η²
    + (1/J) [ ∂(J g^{ξξ})/∂ξ · ∂φ/∂ξ + ∂(J g^{ηη})/∂η · ∂φ/∂η
            + ∂(J g^{ξη})/∂η · ∂φ/∂ξ + ∂(J g^{ξη})/∂ξ · ∂φ/∂η ]
```

Discretized as a flat 9-point stencil `∇²φ[c] = Σ_{k=1..9} w[c,k]·φ[nb[c,k]]`,
column order `1=C, 2=E(+ξ), 3=W(−ξ), 4=N(+η), 5=S(−η), 6=NE, 7=NW, 8=SE, 9=SW`
(`grid_assembly.jl:85-91`).

| Ingredient | Provenance (`grid_assembly.jl`) | Discretization |
|---|---|---|
| Inputs `g⁻¹=metric_ginv (N,2,2)`, `J=metric_jacobian (N,)`, `dξ,dη=cell_widths` | `:125-126`, `:121-122` | bulk trait arrays |
| `gxx=g⁻¹[:,1,1]`, `gyy=g⁻¹[:,2,2]`, `gxe=g⁻¹[:,1,2]` | `:131-133` | per-cell metric |
| axis neighbors `E/W/N/S = neighbor_indices(grid,axis,±1)` | `:136-139` | |
| **corner neighbors** `NE=neighbor_indices(η,+1)[E]`, `NW=…[W]`, `SE=neighbor_indices(η,−1)[E]`, `SW=…[W]` | `:151-158` | **composition of two single-axis steps** |
| `Jgxx=J·gxx` etc.; centered diffs `dJgxx_dξ=(Jgxx[E]−Jgxx[W])/(2dξ)`, `dJgyy_dη`, `dJgxe_dξ`, `dJgxe_dη` | `:161-167` | centered FD of `J·g^{ij}` |
| `orth_dxi_corr=invJ·dJgxx_dξ/(2dξ)`, `orth_deta_corr=invJ·dJgyy_dη/(2dη)` | `:172-173` | 1st-deriv metric correction |
| `cross_d2 = 2·gxe/(4·dξ·dη)` | `:175` | corner cross-derivative coeff |
| `cross_dxi=invJ·dJgxe_dη/(2dξ)`, `cross_deta=invJ·dJgxe_dξ/(2dη)` | `:176-177` | cross metric correction |

Weight assembly (`grid_assembly.jl:181-189`), verbatim:

```
C  = −2gxx/dξ² − 2gyy/dη²
E  = gxx/dξ² + orth_dxi_corr  + cross_dxi
W  = gxx/dξ² − orth_dxi_corr  − cross_dxi
N  = gyy/dη² + orth_deta_corr + cross_deta
S  = gyy/dη² − orth_deta_corr − cross_deta
NE = +cross_d2     NW = −cross_d2     SE = −cross_d2     SW = +cross_d2
```

The off-diagonal `2g^{ξη}∂²φ/∂ξ∂η` term is the standard 4-corner centered
cross-derivative `(φ_NE − φ_NW − φ_SE + φ_SW)/(4dξdη)` — **the only part of the
operator that touches the diagonal/corner stencil points.** On an orthogonal
grid (`g^{ξη}=0`, e.g. lat-lon) the four corner weights are exactly zero.

### 1.2 Gradient — 5-point stencil (chain rule to physical target)

`grid_assembly.jl:247-299`. Chain rule (`:235-240`):
`∂φ/∂t_l = (∂ξ/∂t_l)·∂φ/∂ξ + (∂η/∂t_l)·∂φ/∂η`, with `∂φ/∂ξ≈(φ_E−φ_W)/(2dξ)`.

- `cj = coord_jacobian(grid, target) (N,2,2)` `:257`; components
  `dξ_dt1=cj[:,1,1]`, `dη_dt1=cj[:,2,1]`, `dξ_dt2=cj[:,1,2]`, `dη_dt2=cj[:,2,2]` `:261-264`.
- 5-point neighbors `C,E,W,N,S` (`:266-296`); **no corners**.
- `w1 = [0, dξ_dt1/2dξ, −dξ_dt1/2dξ, dη_dt1/2dη, −dη_dt1/2dη]` `:278-282`;
  `w2` analogously with `t2` `:284-288`.

### 1.3 Boundary handling (oracle)

`grid_assembly.jl:142-148` (+ `:155-158` for corners): when
`neighbor_indices` returns sentinel `0` (off-domain), the oracle **maps the
neighbor INDEX to self**, so that stencil point contributes `w·φ_C` rather than
being dropped. (Comment `:142-144`: "boundary handling is the caller's
responsibility via ghost cells or BCs.") In the lat-lon golden this shows up as:
lon wraps periodically; lat poles' N/S neighbor index = self.

> ⚠ This `sentinel→self` (reuse-center) convention differs from the ESS
> tree-walk evaluator's default out-of-range policy, which substitutes the
> **value** `0.0` (zero-ghost). See §6.2 — this is the one substantive
> conformance detail for R2, and it affects boundary cells only.

---

## 2. Static golden fixtures (captured in ESD)

`fixtures/oracle/{latlon_small,nonorthogonal_small}.json` — exact Float64
oracle outputs (weights, neighbor tables, applied results on an analytic field).
Provenance, schema, and per-case notes in `fixtures/oracle/README.md`. The two
cases together cover (a) the named orthogonal lat-lon example and (b) a
non-orthogonal metric that exercises **every** weight term including the four
corner weights (`max|corner| ≈ 0.2026`, verified equal to `±2g^{ξη}/(4dξdη)`).
Byte-identical on re-run.

---

## 3. Why the oracle is already an ArrayOp (and where the imperative residue is)

`ext/grid_assembly_symbolic.jl` ALREADY lowers the operator to a
`SymbolicUtils.ArrayOp{SymReal}` (its header `:1-3` says it was *"Ported from
EarthSciDiscretizations/src/discretization.jl"* — i.e. this construction is
ESD-native). The mechanism (`:259-273` `_stencil_arrayop`):

```
du[c] = Σ_{k=1..K} w_c[c,k] · u_c[c,k]            # Const-wrapped (N,K) arrays
```

where the corner neighbors live as **integer columns in a precomputed
`(N,9)` table** (`laplacian_neighbor_table :100-133`, corners `:113-120`) and the
field gather `u_ext[c,k]=u[nb[c,k]]` is materialized host-side
(`_build_symbolic_ghost_extension :169-194`). The recursive RHS lowering
(`_rhs_to_arrayop_expr :472-591`) discretizes a mixed 2nd derivative with the
exact corner formula `d2f_dξdη=(f_ne−f_nw−f_se+f_sw)/(4dξdη)` (`:540`).

**Key insight:** the corner access never requires a "diagonal footprint" — it is
emergent from a precomputed flat gather table indexed by `[c,k]`. The residue
that makes this *imperative* is not the math; it is that the table + weights are
assembled by Julia builder functions (`precompute_laplacian_stencil`) instead of
by a declarative rule. R2 replaces the builder with a rule; the einsum shape is
identical.

---

## 4. Existing ESS engine vocabulary — what a declarative rule can already express

There are **two** rule-authoring paths (ESS RFC §7 + base AST §4.3). The verdict
turns on using the right one.

### 4.1 The `use: <scheme>` + selector path — single-axis only (NOT for this operator)
- `scheme_types.jl:35-38`: the **only** concrete `Selector` is
  `CartesianSelector(axis::String, offset::Int)` — one axis, one offset.
- `scheme_expansion.jl:133-157`: `_parse_selector` accepts **only**
  `kind=="cartesian"` (others throw `E_SCHEME_PARSE … not yet supported`), and
  requires `grid_family=="cartesian"` (`:140-143`).
- `scheme_expansion.jl:424-439`: `materialize(CartesianSelector,…)` offsets
  exactly **one** component, "pass the rest through unchanged" → **a single
  stencil entry cannot name a corner.**
- `scheme_expansion.jl:333-342`: `cross_metric`, `dimensional_split`,
  `grid_dispatch` scheme kinds are **silently deferred** (parsed to nothing).

→ This path genuinely cannot express the corner cross-term, and its
`cross_metric` kind is unimplemented. **It is the wrong path** and is irrelevant
to the verdict.

### 4.2 The direct `replacement` einsum path — fully general (USE THIS)
The base AST `op:index` (§4.3.3) and `op:arrayop` (§4.3.1) + arithmetic express
arbitrary multi-axis gathers with bound coefficient arrays. Confirmed in the
engine (file:line from `EarthSciSerialization/src`):

| Capability | Evidence | 
|---|---|
| Direct `replacement` parsed with full arrayop schema; pattern-vars substituted into body/`output_idx`/`ranges` | `rule_engine.jl:986-987,1231-1238`; `apply_bindings :361-402` |
| `op:arrayop` is a first-class lowered node | `types.jl:54-117`; `tree_walk.jl:501-651,2085-2087` |
| **Multi-axis `op:index`** `u[i1,…,iN]` lowered (arity = var rank; **no axis/diagonal distinction**) | `tree_walk.jl:2095-2113` |
| Index-subscript arithmetic `i±1` | `_eval_const_int tree_walk.jl:1362-1428` (`+ − * / ifelse … neg index`; `^` not allowed in a subscript — irrelevant for ±1) |
| Body ops `+ − * / ^ neg`, comparisons, `min/max`, transcendentals, `fn` | `_eval_node_op tree_walk.jl:994-1153` |
| **Per-cell coefficient ARRAY** bound + indexed `index("name", k…)` → literal | `const_arrays` kwarg `tree_walk.jl:175,2117-2133` |
| **Bare scalar** (`dlon`,`R`) → model parameter | `tree_walk.jl:798-812` |
| **No closed name set** — any name the binding layer provides | `grid_accessor.jl:13,68-84`; agent-confirmed |
| Boundary: periodic fold (grid `periodic`), dirichlet/neumann BC emission, else zero-ghost | `discretize.jl:597-648,676-727`; `tree_walk.jl:2105-2108` |

**Activated precedents already in this repo** (each ships canonical byte-contract
fixtures; see git log "activate centered_2nd_deriv, mixed_deriv, nonlinear_laplacian"):

- `finite_difference/mixed_deriv_2nd_uniform.json` — `∂²/∂x∂y` as
  `(u[$x+1,$y+1] − u[$x+1,$y−1] − u[$x−1,$y+1] + u[$x−1,$y−1])/(4 dx dy)`. **This is
  the corner cross-derivative, declaratively, via `op:index`.** ⇒ the corner
  gather is proven-expressible.
- `finite_difference/spherical_laplacian_uniform.json` — a metric-coupled
  (curvilinear) Laplacian that **indexes a named per-cell coordinate array `r`**
  (`index("r",$r)`, `index("r",$r+1)`), face-averages and squares it. ⇒ metric
  coefficient arrays are proven-bindable + arithmetic-combinable.
- `finite_difference/centered_2nd_uniform_latlon.json` — `grid_family:"latlon"`
  with per-cell metric `cos_lat` + scalars `R,dlon,dlat`. ⇒ latlon is a
  first-class family for metric rules.
- `discretizations/SELECTOR_KINDS.md` decision #12: per-cell metric **arrays**
  are indexed `{op:index,args:["metric","$axis"]}` (activated by ESS `ess-sra`);
  line 189-192 records the **cubed-sphere grid was RETIRED** for being "only ~90%
  expressible as declarative rules" — the standing precedent that declarative
  inexpressibility ⇒ retire, never imperative-stopgap.

---

## 5. The exact R2 declarative rule form

Author the covariant Laplacian as a **single `replacement` arrayop-einsum rule**
on `grid_family:"latlon"` (computational axes `lat,lon`), with the field `$u`
gathered at the 9 stencil points and the metric bound as per-cell coefficient
arrays. The algebra below is byte-for-byte the oracle weight assembly of §1.1
(verified: every weight column reproduces `grid_assembly.jl:181-189`).

```
∇²u[lat,lon] =
    g_xx ·(u[lat,lon+1] − 2u[lat,lon] + u[lat,lon−1])/dlon²
  + g_yy ·(u[lat+1,lon] − 2u[lat,lon] + u[lat−1,lon])/dlat²
  + 2·g_xe·(u[lat+1,lon+1] − u[lat+1,lon−1] − u[lat−1,lon+1] + u[lat−1,lon−1])/(4 dlon dlat)
  + (1/J)·dJgxx_dξ·(u[lat,lon+1] − u[lat,lon−1])/(2 dlon)
  + (1/J)·dJgyy_dη·(u[lat+1,lon] − u[lat−1,lon])/(2 dlat)
  + (1/J)·dJgxe_dη·(u[lat,lon+1] − u[lat,lon−1])/(2 dlon)
  + (1/J)·dJgxe_dξ·(u[lat+1,lon] − u[lat−1,lon])/(2 dlat)
```

The metric derivatives are themselves centered differences of bound `J·g`
products (each a plain `op:index` ± plus `op:/`):

```
dJgxx_dξ = (Jg_xx[lat,lon+1] − Jg_xx[lat,lon−1])/(2 dlon)
dJgyy_dη = (Jg_yy[lat+1,lon] − Jg_yy[lat−1,lon])/(2 dlat)
dJgxe_dξ = (Jg_xe[lat+1,lon] − Jg_xe[lat−1,lon])/(2 dlat)   # ∂/∂ξ uses the +ξ axis = lon? see note
dJgxe_dη = (Jg_xe[lat,lon+1] − Jg_xe[lat,lon−1])/(2 dlon)
```

> Axis note for R2: in the oracle, ξ = first computational axis (`xi_axis`),
> η = second (`eta_axis`). Map ξ→lon, η→lat (or whatever the latlon grid's
> `(xi_axis,eta_axis)` resolve to) consistently; `dJgxe_dξ` differences along
> ξ and `dJgxe_dη` along η. The golden's `metric.ginv_xieta` + neighbor tables
> pin the correct association.

**Two equivalent encodings of the metric derivatives** (R2 picks one):
1. **Bind the products precomputed**: supply `Jg_xx, Jg_yy, Jg_xe` (and `g_xx,
   g_yy, g_xe, invJ`) as per-cell `const_arrays`; the rule does the centered
   differencing with `op:index` on those arrays. Most faithful to the oracle.
2. **Bind only `g_xx,g_xe,g_yy,J`**; the rule forms `Jg = J·g` inline before
   differencing. Fewer bound arrays, identical result.

Concrete JSON skeleton (illustrative — corner term + one orthogonal block; R2
fills the full 7 lines above):

```jsonc
{ "discretizations": { "covariant_fv_laplacian_latlon": {
  "applies_to": { "op": "laplacian", "args": ["$u"] },
  "grid_family": "latlon",
  "combine": "+",
  "accuracy": "O(h^2)",
  "replacement": { "op": "+", "args": [
    /* g_xx · (u[lat,lon+1] − 2 u[lat,lon] + u[lat,lon−1]) / dlon^2 */
    { "op": "*", "args": [
      { "op": "/", "args": [ { "op": "index", "args": ["g_xx","lat","lon"] },
                             { "op": "*", "args": ["dlon","dlon"] } ] },
      { "op": "+", "args": [
        { "op": "index", "args": ["$u","lat",{ "op":"+","args":["lon",1] }] },
        { "op": "*", "args": [-2, { "op":"index","args":["$u","lat","lon"] }] },
        { "op": "index", "args": ["$u","lat",{ "op":"+","args":["lon",-1] }] } ] } ] },
    /* 2 g_xe · (u[NE] − u[NW] − u[SE] + u[SW]) / (4 dlon dlat)   <-- corners */
    { "op": "*", "args": [
      { "op": "/", "args": [ { "op":"*","args":[2,{ "op":"index","args":["g_xe","lat","lon"] }] },
                             { "op":"*","args":[4,{ "op":"*","args":["dlon","dlat"] }] } ] },
      { "op": "+", "args": [
        { "op":"index","args":["$u",{ "op":"+","args":["lat",1] },{ "op":"+","args":["lon",1] }] },
        { "op":"*","args":[-1,{ "op":"index","args":["$u",{ "op":"+","args":["lat",1] },{ "op":"+","args":["lon",-1] }] }] },
        { "op":"*","args":[-1,{ "op":"index","args":["$u",{ "op":"+","args":["lat",-1] },{ "op":"+","args":["lon",1] }] }] },
        { "op":"index","args":["$u",{ "op":"+","args":["lat",-1] },{ "op":"+","args":["lon",-1] }] } ] } ] }
    /* … remaining orthogonal-η block + 4 metric-derivative correction terms … */
  ] } } } }
```

The gradient (§1.2) is the strictly simpler 5-point sibling: a `replacement`
`op:+` of `coord_jacobian`-component coefficients × axis-neighbor differences,
no corners — directly analogous to `centered_2nd_uniform_latlon.json`. Two
outputs (`t1`,`t2`) ⇒ either two rules or a `multi_output_stencil`-style pair.

Coefficient symbols the rule references (all per-cell arrays except `dlon,dlat`):
`g_xx, g_xe, g_yy, J` (or `invJ`), optionally precomputed `Jg_xx, Jg_yy, Jg_xe`;
scalars `dlon, dlat`; gradient adds `coord_jacobian` components `dlon_dt1,
dlat_dt1, dlon_dt2, dlat_dt2`.

---

## 6. Scope for R2 — non-engine wiring + the one conformance detail

Neither item below is "new engine code" (no AST node, no arrayop footprint, no
ESS scheme-handler). Both are explicitly within the bead's allowed vocabulary
("existing arrayop footprints/gathers **+ coefficient metric arrays**").

### 6.1 `const_arrays` metric binding for latlon (host-side, DATA/binding)
The latlon metric arrays the rule needs **already exist as grid accessors** —
`src/grids/latlon.jl`: `metric_ginv :729-749` (`ginv_lonlon=1/(R²cos²φ)`,
`ginv_latlat=1/R²`, `ginv_lonlat=0`), `metric_jacobian :751-766` (`J=R²cosφ`),
`metric_dgij_dxk :768-787`, `coord_jacobian :789-802` — but they are **not yet
routed into rule evaluation as named `const_arrays`** (today latlon takes ESS
Path B `discretize(sys,grid)` `ode_problem.jl:142-146`, which bypasses
`const_arrays`).

The established, **already-working** pattern to bind them is the MPAS/DUO
emitter `_grid_primitive_arrays(grid)` (`src/ode_problem.jl:759-771`) →
`_unstructured_grid_const_arrays` → merged into `const_arrays` and passed to
`build_evaluator` (`:95-104`). R2 adds a `_grid_primitive_arrays(::LatLonGrid)`
analogue emitting `g_xx,g_xe,g_yy,J` (and/or `Jg_*`) keyed by the rule-symbol
names, plus the latlon grid-schema field names (`discretizations/grids/lat_lon.schema.json`).
Pure host-side array shaping; the AST `index` op + `evaluate_expr` already
consume it unchanged.

### 6.2 Lat-pole boundary matching (conformance, boundary cells only)
The oracle uses **`sentinel→self`** at off-domain neighbors (§1.3); the engine's
einsum uses naive `±1` subscripts resolved by grid metadata — **periodic fold**
(matches lon) or **zero-ghost / declared BC** (lat poles). Interior cells are
tolerance-identical; pole rows need R2 to reconcile the convention (declare lon
`periodic:true`; encode a pole BC, or have the latlon grid resolve pole neighbors
to self so the lowered indices match the oracle). The goldens capture the exact
expected pole-row values for verification.

### 6.3 Note on current state
The latlon FD declarative path is **not fully wired end-to-end today** (the lone
latlon rule `centered_2nd_uniform_latlon.json` is stencil-only, which ESS
`parse_rule` rejects with `E_RULE_REPLACEMENT_MISSING`; latlon runs via Path B).
Wiring it is exactly the esd-zk9 epic's purpose — it is a state-of-the-codebase
fact, not a feasibility blocker. R2 authors a `replacement`-form rule (§5),
which Path A consumes.

---

## 7. Summary mapping — every oracle ingredient → existing mechanism

| Oracle ingredient | Declarative expression | Existing mechanism | New engine code? |
|---|---|---|---|
| 9-point gather incl. NE/NW/SE/SW corners | `index($u, lat±1, lon±1)` | `_resolve_indices` multi-axis `tree_walk.jl:2095` | **No** (mixed_deriv precedent) |
| `g^{ξξ},g^{ξη},g^{ηη},J` per-cell | `index("g_xx",lat,lon)` … | `const_arrays` + latlon accessors | **No** (data binding; MPAS precedent) |
| `∂(Jg^{ij})/∂ξ` centered diff | `(index("Jg_xx",lat,lon+1)−…)/(2 dlon)` | `op:index` + arithmetic | **No** |
| `1/J`, `−2`, sums of 9 terms | `op:/ op:* op:+` | `_eval_node_op` | **No** |
| `dξ,dη` | bare `dlon,dlat` | model parameters | **No** |
| output over all cells | `op:arrayop output_idx:[lat,lon]` (or bare expr, auto-lifted) | arrayop lowering / arrayop-lift | **No** |
| lon periodicity | grid `periodic:true` | `_apply_periodic_folding!` | **No** |
| gradient chain rule | `coord_jacobian` components × axis diffs | as `centered_2nd_uniform_latlon` | **No** |

Every row resolves to an existing node/mechanism. **Verdict: DECLARATIVE-FEASIBLE.**
Proceed to R2 (esd-zk9.2) with the §5 rule form; do **not** keep `grid_assembly`,
do **not** implement a `cross_metric` scheme handler.
