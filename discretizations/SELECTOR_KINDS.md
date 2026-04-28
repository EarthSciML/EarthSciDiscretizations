# Selector Kinds — Rule Schema Index

This document is the running index of stencil-selector design decisions made
while porting discretization rules across grid families. Each row records
**one** resolved question, the chosen convention, and the bead that pinned it.
New rule authors should read this before introducing a new grid family.

The selector schema lives inside the ESS §7 rule format:

```json
{
  "stencil": [
    {
      "selector": { "kind": "<family>", "axis": "<axis>", "offset": <int> },
      "coeff":    { "...AST..." }
    }
  ]
}
```

`kind` binds the offset to a specific grid-family accessor. `axis` names the
direction along which `offset` is interpreted. `coeff` is an
`ExpressionNode` AST evaluated against per-grid bindings (see ESS §9.2).

## Decisions

| # | Question | Decision | Pinned by |
|---|---|---|---|
| 1 | Selector `kind` naming: per-family (`vertical`, `latlon`, `cubed_sphere`) or structural (`structured_1d`, `structured_2d_metric`)? | **Per-family.** Each kind names exactly one grid-family accessor; structural categories may emerge once 3+ families share a stencil but are premature today. | dsc-mzu |
| 2 | Spacing symbol for uniform 1D structured grids (cartesian, vertical-uniform). | **`h`** is the canonical generic uniform-spacing symbol. Per-axis aliases (`dx`, `dz`, …) are also bound by the MMS harness so existing cartesian rules keep working. New rules SHOULD prefer `h` to make grid-family-agnostic AST reuse possible. | dsc-mzu |
| 3 | Metric for non-uniform vertical grids (`eta`, `theta`, `z` stretched). | **Defer.** The convention will be a per-cell symbolic reference `dz_k` resolved against the grid accessor at evaluation time. Implementing it requires extending the ESS evaluator + grid binding contract; tracked as a follow-up bead. Until then, only uniform vertical (`sigma_uniform`, `z_uniform`) is supported by Layer B. | dsc-mzu (deferred to follow-up) |
| 4 | Boundary handling for non-periodic vertical (top/bottom). | **Compose with a separate BC rule** (parallel to `periodic_bc.json`) rather than inlining BC handling in the stencil. Keeps the centered stencil rule clean and lets the same BC rule cover all 1D structured families. Concrete `vertical_top_bottom_bc.json` is a follow-up. | dsc-mzu (deferred to follow-up) |
| 5 | One rule with `grid_family: ["cartesian", "vertical"]` vs sibling files. | **Sibling files** when the selector `kind` differs (which is always, given decision #1). A `grid_family` array is reserved for cases where the AST and selectors are byte-identical — currently no such case exists. | dsc-mzu |
| 6 | MPAS selector shape: a per-family `kind: "mpas"` with an integer `neighbor` slot, or the upstream ESS structural kinds (`indirect`, `reduction`) with a `table` field naming a connectivity table? | **Structural — adopt ESS verbatim.** MPAS uses `kind: "indirect"` for fixed-arity stencils (one neighbor per stencil row, e.g. edge→cell-pair) and `kind: "reduction"` for variable-valence stencils (one stencil row that lowers to an `arrayop` reduction over a contiguous index range, e.g. cell→edges). The `table` field names a connectivity array declared on the grid (`edges_on_cell`, `cells_on_edge`, `vertices_on_edge`); this is the family-specific accessor of decision #1, just routed through the table identifier instead of the kind name. A bare `kind: "mpas"` with an integer neighbor slot was rejected because it cannot express variable valence (12 pentagons + N hexagons share one mesh) without sentinel padding, and it would diverge from the upstream `ExpressionNode` lowering path defined in ESS RFC `discretization.md` §4 / §7.2. | dsc-6xp |
| 7 | Variable connectivity (cells with 5 vs 6 neighbors): fixed-max with sentinel (`-1`) or true variable-arity? | **True variable-arity via `reduction` selector.** A `reduction` selector carries a `count_expr` (typically `index(n_edges_on_cell, $target)`) and a `k_bound` local index name; expansion (ESS §7.2) lowers exactly one stencil row to an `arrayop` over `k = 0 .. count_expr − 1`. Bindings never see sentinel `-1` entries and never branch on cell valence. Fixed-max-with-sentinel was rejected because it pollutes the AST with conditionals and breaks bit-identity across bindings whose padding conventions differ. | dsc-6xp |
| 8 | Edge-based vs cell-based stencils (MPAS quantities live at cells, edges, or vertices): one selector kind per stagger location, or one shared kind plus a target-location tag on the scheme? | **Shared kind, per-scheme location tag.** `indirect` and `reduction` are agnostic to whether the target is a cell, edge, or vertex; the location is set by the scheme's `emits_location` (`cell_center`, `edge_normal`, `vertex`) and ESS §7.1.1's `$target` chooser binds the scalar `c` / `e` / `v` accordingly. The connectivity table named by the selector resolves the rest (e.g. `cells_on_edge` for edge→cell-pair, `edges_on_cell` for cell→edge-fan). No new selector kinds were introduced for stagger variants. | dsc-6xp |
| 9 | Coefficient symbols available inside MPAS rule `coeff` ASTs. | **The `dsc-7j0` grid-accessor symbol set, snake_case.** Metric arrays: `area_cell`, `dc_edge`, `dv_edge` (and lon/lat/x/y/z per element kind). Connectivity tables: `edges_on_cell`, `cells_on_edge`, `cells_on_cell`, `vertices_on_edge`. Valence array: `n_edges_on_cell`. Symbol names follow the `src/grids/mpas.jl` runtime fields verbatim (and the matching Python/TS bindings); rule authors MUST NOT introduce new bare-name symbols without extending the MPAS grid schema (`grids/mpas.schema.json`). The ESS RFC's camelCase examples (`areaCell`, `nEdgesOnCell`) are illustrative — ESD uses snake_case to match the binding runtime. | dsc-6xp |
| 13 | How to express the diagonal (self) contribution of a Laplace-style operator on MPAS, given that `reduction` selectors fix the operand to `index(operand, table[$target, k])` (RFC §7.2) and cannot point at `$target` itself. | **One `reduction` row for the off-diagonal sum + one `indirect` row for the self-term, with the diagonal weight as an `arrayop` coefficient.** Row 1 sums Σ_k w_k u(n_k) using `kind: "reduction"` over `cells_on_cell`; Row 2 contributes −W(c) u(c) using `kind: "indirect"` with `index_expr: "$target"` and a coeff of `−arrayop(w_k)` reducing over the same `k = 0 .. n_edges_on_cell[$target] − 1` range with `reduce: "+"`. This closes the algebraic identity ∇²u(c) = Σ_k w_k (u(n_k) − u(c)) without inventing a "self" connectivity table or extending the selector schema — `arrayop` is already in the base AST (ESS §4.3.1) and `index_expr` is a free `ExpressionNode` per the RFC `indirect` materialization row. The pattern generalizes: any cell-centered reduction-style operator whose continuous form pairs each neighbor-difference against a self-reference uses this two-row split. | dsc-ion |

## MPAS (unstructured Voronoi) — selector schema

MPAS lacks logical (i, j, k) offsets — neighbors are addressed by integer
**connectivity tables** loaded from the mesh file (`edges_on_cell`,
`cells_on_edge`, `vertices_on_edge`, …). Two structural selector kinds (per
upstream ESS RFC `discretization.md` §4) cover every MPAS stencil pattern:

### `kind: "indirect"` — fixed-arity neighbor

Use when each stencil row references **one** neighbor whose index is a
deterministic function of `$target` and a known integer slot (e.g. the two
cells adjacent to an edge).

```jsonc
{
  "selector": {
    "kind":       "indirect",
    "table":      "cells_on_edge",
    "index_expr": { "op": "index", "args": ["cells_on_edge", "$target", 0] }
  },
  "coeff": <ExpressionNode>
}
```

- `table` — name of a connectivity array declared on the grid (must exist
  in the grid's `connectivity` block; see `grids/mpas.schema.json`).
- `index_expr` — single `ExpressionNode` resolving to the operand's index.
  Composes with `$target` (here: edge index `e`) and integer slots. The
  expansion (ESS §7.2 `indirect` row) substitutes this expression directly
  into the operand reference.

### `kind: "reduction"` — variable-arity reduction

Use when one stencil row stands in for a sum/product/min/max **over a
variable number of neighbors** (e.g. divergence at a cell summing over its
5 or 6 incident edges).

```jsonc
{
  "selector": {
    "kind":       "reduction",
    "table":      "edges_on_cell",
    "count_expr": { "op": "index", "args": ["n_edges_on_cell", "$target"] },
    "k_bound":    "k",
    "combine":    "+"
  },
  "coeff": <ExpressionNode that may reference $target and k>
}
```

- `table` — connectivity table; element `index(table, $target, k)` is the
  k-th neighbor of `$target`.
- `count_expr` — `ExpressionNode` that evaluates to the inclusive valence
  count at `$target`; `index(n_edges_on_cell, $target)` is canonical for
  cell-to-edge fans.
- `k_bound` — local index name introduced by this selector and brought into
  scope inside `coeff` per ESS §7.1.1 (not a `$`-prefixed pattern var).
- `combine` — `"+"`, `"*"`, `"min"`, `"max"`. Almost always matches the
  scheme's outer `combine`; mixed reductions are legal but not yet used.

Expansion (ESS §7.2 `reduction` row) lowers this single stencil row to an
`arrayop` over `k = 0 .. count_expr − 1`. **Variable valence is handled
purely in the AST** — no `-1` sentinel, no per-cell branching at the
binding layer.

### Stagger location and `$target`

The scheme's `emits_location` (`cell_center` / `edge_normal` / `vertex`)
determines whether `$target` resolves to a scalar `c`, `e`, or `v` per ESS
§7.1.1. When the operand sits at a different stagger location from the
emit location (the common case for divergence and gradient), the
`indirect` / `reduction` selector's `table` field bridges the two.

### Worked example — flux divergence at cell centers (MPAS C-grid)

The MPAS divergence operator combines a normal-velocity flux `F` defined at
edge midpoints into a cell-centered tendency:

$$\nabla\!\cdot\!F\;\big|_{c}\;=\;\frac{1}{A_c}\sum_{e\in\mathcal{E}(c)} L_e\, F_e$$

where $A_c$ is `area_cell`, $L_e$ is the dual-edge length `dv_edge`, and
$\mathcal{E}(c)$ is the set of `n_edges_on_cell[c]` edges incident on cell
$c$. Expressed against this document's selector schema:

```jsonc
{
  "discretizations": {
    "mpas_cell_div": {
      "applies_to":   { "op": "div", "args": ["$F"], "dim": "cell" },
      "grid_family":  "unstructured",
      "requires_locations": ["edge_normal"],
      "emits_location":     "cell_center",
      "combine":      "+",
      "accuracy":     "O(dx) on a quasi-uniform Voronoi mesh",
      "stencil": [
        {
          "selector": {
            "kind":       "reduction",
            "table":      "edges_on_cell",
            "count_expr": { "op": "index", "args": ["n_edges_on_cell", "$target"] },
            "k_bound":    "k",
            "combine":    "+"
          },
          "coeff": {
            "op": "/",
            "args": [
              { "op": "index", "args": [
                  "dv_edge",
                  { "op": "index", "args": ["edges_on_cell", "$target", "k"] }
              ]},
              { "op": "index", "args": ["area_cell", "$target"] }
            ]
          }
        }
      ]
    }
  }
}
```

After §7.2 expansion at cell `c` (`$target → c`), this lowers to a single
`arrayop` reduction whose lowered shape is given verbatim in ESS RFC
`discretization.md` §7.3. The lowered form parses against the base-spec
`arrayop` schema (no new AST node).

### Cross-references

- ESS RFC `discretization.md` §4 (selector taxonomy), §6.3 (unstructured
  grid schema), §7.2 (expansion semantics), §7.3 (worked MPAS divergence),
  §7.4 (`staggering_rules` for quantity-location declarations).
- `grids/mpas.schema.json` — MPAS grid family schema (loader-backed).
- `src/grids/mpas.jl` — Julia accessor runtime (canonical symbol names).
- `dsc-7j0` — grid accessor runtime that pins the symbol set.
- ESS `esm-spec.md` — base AST spec; the `index` op (§4.3.3) and `arrayop`
  (§4.3.1) are the only nodes the lowered MPAS stencils require. **No
  amendment to `esm-spec.md` is needed** for MPAS support; the structural
  selectors above lower to existing AST nodes only.

## More structured-grid decisions (latlon)

| # | Question | Decision | Pinned by |
|---|---|---|---|
| 10 | Selector `axis` for `latlon`: pattern variable (`$a`) vs literal (`"lon"`/`"lat"`). | **Literal**, when stencil entries for different axes carry different coefficients (e.g. lon needs `cos_lat`, lat does not). The rule lists 4 stencil entries — two per axis — and the rule engine picks the pair whose `axis` matches the operator's bound `dim`. Pattern-variable axes (`$a`) remain valid for grid families where every axis shares the same coefficient form (e.g. cartesian, vertical). | dsc-8ad |
| 11 | Spacing symbol convention for `latlon`: angular (`dlon`/`dlat` in radians) vs physical (already multiplied by R). | **Angular**, with `R` as a separate binding. `dlon` and `dlat` are per-axis angular increments in radians; `R` is the sphere radius from the grid accessor; `cos_lat` is the per-cell latitude-metric. Coefficients combine them explicitly (e.g. `1/(2 R cos_lat dlon)`). This keeps the AST honest — every dimensional factor appears in the formula rather than being baked into the spacing symbol. | dsc-8ad |
| 12 | Per-cell metric bindings (`cos_lat` for latlon, `dz_k` for stretched vertical). | **Per-family grid-accessor binding**. Per-cell metrics are evaluated against the grid accessor at stencil-application time — distinct from scalar bindings (`R`, `dlon`, `dlat`) which are constant across the sweep. The current ESS MMS harness only supports scalar bindings, so latlon's Layer B fixture marks `applicable: false` until the harness extension lands (follow-up bead dsc-7jc). | dsc-8ad |
| 13 | Selector for `cubed_sphere`: how to encode 2D in-panel offsets, and where does cross-panel ghost handling live? | **Per-axis single-offset selectors composed via a `selectors: [...]` array on each stencil entry** — kind `cubed_sphere`, axes `xi` and `eta`, `offset: <int>` interpreted in-panel. A 9-point stencil entry carries two selectors (one per axis); a 1D-along-axis entry carries one. Cross-panel ghost extension and basis rotation live in the cubed_sphere grid accessor (`src/grids/panel_connectivity.jl`) — selectors do **not** carry a `panel` field. The accessor resolves out-of-panel offsets to the right neighbor cell with the correct rotation matrix. Spacing is uniform isotropic on the gnomonic grid (`dξ = dη = π/(2·Nc)`), so the canonical spacing symbol is `h` (decision #2). | dsc-ap9 |
| 17 | D-grid (and corner) staggering on `cubed_sphere`: how does a stencil entry name *which* sample location it lives on (cell center vs ξ-constant face vs η-constant face vs corner)? | **Categorical `stagger` field on the selector**, mirroring decision #13 for `arakawa`. Values are `cell_center`, `u_edge` (ξ-constant face — V·ê_η at UEdge in FV3 D-grid), `v_edge` (η-constant face — V·ê_ξ at VEdge), `corner` (cell vertex). The full selector form on this family is `{ kind: "cubed_sphere", stagger: <symbol>, axis: <axis>, offset: <int> }`. The `stagger` slot is OPTIONAL and defaults to `cell_center` so existing scalar-only rules (e.g. `covariant_laplacian_cubed_sphere`) need no change. Encoded as a file-local enum per ESS §9.3, named `cubed_sphere_stagger`. The `axis`/`offset` semantics are unchanged — `offset` is the integer step from the output cell along `axis` interpreted at the selector's stagger location. Cross-panel ghost extension (scalar or vector with rotation) still lives in the grid accessor; selectors do not carry a `panel` field. | dsc-247 |
| 18 | Per-cell sub-grid metric bindings available to `cubed_sphere` rules — sin/cos of the angle between coordinate axes (`sin_sg`, `cos_sg`) at the 9 FV3 super-grid sample points, and physical edge lengths/areas (`dx`, `dy`, `area`). | **Per-cell grid-accessor bindings**. `dx[p,i,j]`, `dy[p,i,j]`, `area[p,i,j]` and `lat[p,i,j]` follow the field names on `CubedSphereGrid` verbatim. Sub-grid `sin_sg[p,i,j,pos]`/`cos_sg[p,i,j,pos]` for `pos ∈ 1..9` (1=west mid-edge, 2=south, 3=east, 4=north, 5=cell center, 6–9=corner positions) are exposed as scalar bindings indexed via `{op: "index", args: ["sin_sg", "$target", <pos>]}`. Like `cos_lat` for latlon (decision #12), these are evaluated against the grid accessor at stencil-application time. The Coriolis rotation rate `Omega_rot` is a scalar binding the rule consumer supplies (default 7.292e-5 rad/s for Earth). | dsc-247 |

## More structured-grid decisions (arakawa)

| # | Question | Decision | Pinned by |
|---|---|---|---|
| 13 | Stagger-position selector for the `arakawa` family — how does a stencil entry name *which* C/B/D-grid sample location it lives on? | **Categorical `stagger` field on the selector**, with values `cell_center`, `face_x`, `face_y`, `vertex` (closed set; B/D-grids reuse the same vocabulary). The `axis` and `offset` retain their structured-grid meanings: `offset` is the integer step from the output cell along `axis`. The full selector form is `{ kind: "arakawa", stagger: <symbol>, axis: <axis>, offset: <int> }`. | dsc-cuj |
| 14 | How are the categorical `stagger` values encoded portably across bindings? | **File-local `enums` block per ESS §9.3.** Each Arakawa rule that uses stagger selectors declares an `arakawa_stagger` enum mapping the four symbols to positive integers; the rule's optional `stagger_enum` field names which enum to consult. Bindings lower the selector's string `stagger` to the corresponding integer at load time (mirroring the §4.5 `enum`-op lowering contract), so evaluators only ever see integers. Authors keep human-readable symbols in source. | dsc-cuj (esm-mqc) |
| 15 | Output stagger for an Arakawa rule — is it implicit from the operator (`div` → cell-center) or declared? | **Declared explicitly** via `emits_location: <stagger symbol>` on the scheme, paired with `requires_locations: [<stagger symbol>, …]` for the input fields. This mirrors the unstructured-C-grid worked example in ESS §7.3 (MPAS divergence) and avoids overloading per-operator implicit conventions when adding new ops (curl, grad-on-C, etc.). | dsc-cuj |
| 16 | Mapping a multi-component `applies_to: { op: "div", args: ["$F"] }` to per-stencil-entry components on a structured C-grid. | **Implicit dispatch by the entry's `axis` and `stagger`.** A C-grid's `face_x` location is, by construction, where the x-component of any vector field lives; `face_y` is the y-component. The stencil therefore samples the single pattern variable `$F` and the binding picks the matching axis-component array (e.g. `ux` vs `uy`) from the entry's `(stagger, axis)` pair. No new pattern variables are introduced per component; this matches `arakawa_variable_locations(ArakawaC) = (CellCenter, UEdge, VEdge)` in `src/grids/arakawa.jl`. | dsc-cuj |

## Arakawa (structured staggered) — selector schema

Arakawa A/B/C/D/E grids are structured 2D meshes that stagger different
quantities at four canonical sample locations: cell centers, x-faces,
y-faces, and corners (vertices). The C-grid is the dominant choice in
modern atmosphere/ocean dynamical cores (MITgcm, MOM6, MPAS, FV3, NEMO).

Selectors on this family extend the per-family structured form with a
categorical `stagger` slot:

```jsonc
{
  "selector": {
    "kind":    "arakawa",
    "stagger": "cell_center" | "face_x" | "face_y" | "vertex",
    "axis":    "$x" | "$y",
    "offset":  <int>
  },
  "coeff": <ExpressionNode>
}
```

Categorical `stagger` values are encoded as a file-local `enums` block per
ESS §9.3. Each Arakawa rule SHOULD declare:

```jsonc
{
  "enums": {
    "arakawa_stagger": {
      "cell_center": 1,
      "face_x":      2,
      "face_y":      3,
      "vertex":      4
    }
  },
  "discretizations": { "<rule_name>": { …, "stagger_enum": "arakawa_stagger", … } }
}
```

Bindings lower the selector's string `stagger` to the corresponding integer
at load time (mirroring the §4.5 `enum`-op lowering contract); evaluators
only ever see integers. The `stagger_enum` field on the scheme names the
enum to consult so a rule can carry multiple enums without ambiguity.

Schemes additionally declare quantity locations explicitly:

- `emits_location` — output stagger symbol (e.g. `"cell_center"` for `div`).
- `requires_locations` — list of input stagger symbols, in arg-order.

For multi-component vector arguments (`applies_to: { op: "div", args: ["$F"] }`),
each stencil entry's `axis` and `stagger` together identify which component
array to read: on a C-grid, `face_x` selectors sample the x-component and
`face_y` selectors sample the y-component of `$F`. No per-component pattern
variables are needed.

### Worked example — divergence on Arakawa C (cartesian base)

Two-point centered FV divergence: $\nabla\cdot F = (u_x[i+1,j] - u_x[i,j])/\Delta x +
(u_y[i,j+1] - u_y[i,j])/\Delta y$. See
[`finite_volume/divergence_arakawa_c.json`](finite_volume/divergence_arakawa_c.json)
for the encoded form and `divergence_arakawa_c/fixtures/canonical/` for the
2×2 small-grid sanity fixture.

## More structured-grid decisions (arakawa)

| # | Question | Decision | Pinned by |
|---|---|---|---|
| 13 | Stagger-position selector for the `arakawa` family — how does a stencil entry name *which* C/B/D-grid sample location it lives on? | **Categorical `stagger` field on the selector**, with values `cell_center`, `face_x`, `face_y`, `vertex` (closed set; B/D-grids reuse the same vocabulary). The `axis` and `offset` retain their structured-grid meanings: `offset` is the integer step from the output cell along `axis`. The full selector form is `{ kind: "arakawa", stagger: <symbol>, axis: <axis>, offset: <int> }`. | dsc-cuj |
| 14 | How are the categorical `stagger` values encoded portably across bindings? | **File-local `enums` block per ESS §9.3.** Each Arakawa rule that uses stagger selectors declares an `arakawa_stagger` enum mapping the four symbols to positive integers; the rule's optional `stagger_enum` field names which enum to consult. Bindings lower the selector's string `stagger` to the corresponding integer at load time (mirroring the §4.5 `enum`-op lowering contract), so evaluators only ever see integers. Authors keep human-readable symbols in source. | dsc-cuj (esm-mqc) |
| 15 | Output stagger for an Arakawa rule — is it implicit from the operator (`div` → cell-center) or declared? | **Declared explicitly** via `emits_location: <stagger symbol>` on the scheme, paired with `requires_locations: [<stagger symbol>, …]` for the input fields. This mirrors the unstructured-C-grid worked example in ESS §7.3 (MPAS divergence) and avoids overloading per-operator implicit conventions when adding new ops (curl, grad-on-C, etc.). | dsc-cuj |
| 16 | Mapping a multi-component `applies_to: { op: "div", args: ["$F"] }` to per-stencil-entry components on a structured C-grid. | **Implicit dispatch by the entry's `axis` and `stagger`.** A C-grid's `face_x` location is, by construction, where the x-component of any vector field lives; `face_y` is the y-component. The stencil therefore samples the single pattern variable `$F` and the binding picks the matching axis-component array (e.g. `ux` vs `uy`) from the entry's `(stagger, axis)` pair. No new pattern variables are introduced per component; this matches `arakawa_variable_locations(ArakawaC) = (CellCenter, UEdge, VEdge)` in `src/grids/arakawa.jl`. | dsc-cuj |

## Arakawa (structured staggered) — selector schema

Arakawa A/B/C/D/E grids are structured 2D meshes that stagger different
quantities at four canonical sample locations: cell centers, x-faces,
y-faces, and corners (vertices). The C-grid is the dominant choice in
modern atmosphere/ocean dynamical cores (MITgcm, MOM6, MPAS, FV3, NEMO).

Selectors on this family extend the per-family structured form with a
categorical `stagger` slot:

```jsonc
{
  "selector": {
    "kind":    "arakawa",
    "stagger": "cell_center" | "face_x" | "face_y" | "vertex",
    "axis":    "$x" | "$y",
    "offset":  <int>
  },
  "coeff": <ExpressionNode>
}
```

Categorical `stagger` values are encoded as a file-local `enums` block per
ESS §9.3. Each Arakawa rule SHOULD declare:

```jsonc
{
  "enums": {
    "arakawa_stagger": {
      "cell_center": 1,
      "face_x":      2,
      "face_y":      3,
      "vertex":      4
    }
  },
  "discretizations": { "<rule_name>": { …, "stagger_enum": "arakawa_stagger", … } }
}
```

Bindings lower the selector's string `stagger` to the corresponding integer
at load time (mirroring the §4.5 `enum`-op lowering contract); evaluators
only ever see integers. The `stagger_enum` field on the scheme names the
enum to consult so a rule can carry multiple enums without ambiguity.

Schemes additionally declare quantity locations explicitly:

- `emits_location` — output stagger symbol (e.g. `"cell_center"` for `div`).
- `requires_locations` — list of input stagger symbols, in arg-order.

For multi-component vector arguments (`applies_to: { op: "div", args: ["$F"] }`),
each stencil entry's `axis` and `stagger` together identify which component
array to read: on a C-grid, `face_x` selectors sample the x-component and
`face_y` selectors sample the y-component of `$F`. No per-component pattern
variables are needed.

### Worked example — divergence on Arakawa C (cartesian base)

Two-point centered FV divergence: $\nabla\cdot F = (u_x[i+1,j] - u_x[i,j])/\Delta x +
(u_y[i,j+1] - u_y[i,j])/\Delta y$. See
[`finite_volume/divergence_arakawa_c.json`](finite_volume/divergence_arakawa_c.json)
for the encoded form and `divergence_arakawa_c/fixtures/canonical/` for the
2×2 small-grid sanity fixture.

## When to add a new selector kind

Adding a new family (e.g. `latlon`, `cubed_sphere`):

1. Pick a per-family kind name in lowercase snake_case.
2. Decide the canonical axis names (`(:lon, :lat)`, `(:xi, :eta)`, …).
3. Pick the spacing symbol convention. Reuse `h` for axis-agnostic uniform
   spacing; introduce per-axis symbols (`dlon`, `dlat`) only when metrics
   require it.
4. Add a row above with the rationale and the bead that pinned it.
5. Author the first rule in `discretizations/<family>/` with a paired
   convergence fixture under `<rule>/fixtures/convergence/`.

## See also

- `../docs/rule-catalog.md` — manifest of rules to port across families.
- `README.md` (this directory) — rule-file authoring policy.
- ESS `esm-spec.md` §7 — rule schema; §9.2 — `call`-op decision tree.
