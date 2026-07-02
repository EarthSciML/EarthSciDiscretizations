# reduced_gaussian latlon operators — DECLARATIVE-OR-FAIL verdict

**Bead:** esd-6g4.13 (G14) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Oracle source pin:** EarthSciSerialization GitHub `main` @
`8ba047ff2590ff621e0cc8dac4f6b049a691d161` (the rev ESD `Project.toml`
`[sources]` resolves: `{rev = "main", subdir = "packages/EarthSciSerialization.jl"}`).
ESD line numbers below are at the esd-6g4.13 worktree HEAD.

---

## 0. Verdict (TL;DR)

> **DECLARATIVE — INFEASIBLE over the existing engine.** A reduced-Gaussian
> gradient / divergence / Laplacian rule cannot be authored as a declarative
> arrayop over the existing `+ - * / index` vocabulary. The blocker is the
> **longitude operand index**: on a reduced-Gaussian grid the number of
> longitude cells varies per latitude row (`nlon_per_row[j]`), so a `lon±1`
> neighbour read needs a **row-dependent ragged wrap modulus** — a primitive
> the arrayop `index` op does not have. Per the bead's `DECLARATIVE-OR-FAIL`
> contract, the correct outcome is this gap report — **no rule JSON is
> committed and no imperative operator is hand-coded.**

The **latitude half** of a horizontal operator is expressible at the AST level
(see §3) but is **not exercisable end-to-end**: the only latlon Layer-B runner
is rectangular and variant-blind (§4). So the combo as a whole is reported, not
implemented.

---

## 1. What the bead asked

> "reduced_gaussian latlon operator rules. Declarative rules + fixtures.
> DECLARATIVE-OR-FAIL: the rule MUST be declarative JSON (arrayop/stencil
> composed over the EXISTING ESS engine vocabulary). If it CANNOT be expressed
> declaratively, do NOT implement it — STOP and report the precise gap."

## 2. The grid is a flat *ragged* (jagged) row-major array — not a rectangle

`LatLonGrid` stores `nlon_per_row::Vector{Int}` (`src/grids/latlon.jl:52`) and
has **no** scalar `nlon`; `n_cells(g) = sum(g.nlon_per_row)`
(`src/grids/latlon.jl:65`). Cell `(j, i)` lives at flat index
`row_offset(g, j) + i`, the prefix sum of `nlon_per_row`
(`src/grids/latlon.jl:97`). The `:reduced_gaussian` variant is documented at
`src/grids/latlon.jl:20` with the worked example `nlon_per_row = [4, 8, 12, 12,
8, 4]` (`src/grids/latlon.jl:13`). Per-cell longitude width is genuinely
row-varying: `cell_widths(:lon)` is the ragged array `dlon_j = 2π/nlon_per_row[j]`.

## 3. Latitude derivative — expressible at the AST level, but see §4

The lat term of `centered_2nd_uniform_latlon.json` is
`±1/(2·R·dlat)·index($u, lat±1, lon)`. On reduced_gaussian:

- The metric `dlat` is row-dependent only if `lat_edges` are non-uniform
  (Gaussian). `dlat` is currently a constant scalar binding
  (`discretizations/SELECTOR_KINDS.md:187`: "scalar bindings (`R`, `dlon`,
  `dlat`) which are constant across the sweep"), but it could be lifted to a
  per-cell indexed array `{op:index, args:["dlat","lat"]}` exactly as
  `centered_2nd_nonuniform_vertical.json` does for `dz[k]` — that primitive
  already exists.
- The lat *operand* `index($u, lat±1, lon)` is structurally sound: rows are
  ordered and the cross-row nearest-centre remap already exists in the grid
  runtime.

So the lat-derivative does **not** require new arrayop vocabulary.

## 4. Longitude derivative — INFEASIBLE (the load-bearing gap)

The obstacle is the **operand index** `index($u, lat, lon±1)`, not the
coefficient:

- The arrayop `index` op addresses operands by a fixed logical offset with a
  **single, per-axis-constant wrap modulus**. The latlon Layer-B lowering
  resolves a periodic longitude neighbour with one global `n_lon`; the topology
  classifier returns `"2d_latlon_sphere"` for `grid_family=="latlon"`
  **regardless of `variant`** (`test/walk_esd_tests.jl:976`), and the runner it
  routes to (`_run_layer_b_2d_latlon_sphere`) hardcodes `n_lat = n_lon = n` and
  a single scalar `dlon`. On a ragged grid the wrap modulus is
  `nlon_per_row[lat]` — a function of the *other* index of the very cell being
  read. There is **no `+ - * / index` form for an index whose wrap extent
  depends on another index.** (The only relational/ragged addressing the engine
  offers is the MPAS `kind:"indirect"`/`"reduction"` connectivity-table
  selectors — `discretizations/SELECTOR_KINDS.md` — a different family-specific
  selector that latlon arrayop rules cannot use.)
- Secondary obstacle (physics): adjacent rows have different cell counts, so a
  convergent horizontal lon-gradient generally needs **inter-row interpolation**
  (the nearest-centre remap is O(1)-accurate, fine for topology, not for a
  second-order stencil). Interpolation is likewise not a fixed-offset arrayop.

So the longitude derivative is doubly blocked: (a) no ragged-wrap `index`,
(b) no interpolation primitive.

## 5. No runner, no matrix column

Even the feasible lat-half cannot be validated: the only latlon operator runner
is variant-blind and rectangular (`test/walk_esd_tests.jl:976`, `:732`), so it
cannot ingest `nlon_per_row`. The rule matrix has no `variant` axis — columns are
family-granular only (`tools/render_rule_matrix.py` `GRID_FAMILIES`), so a
reduced_gaussian operator rule would be indistinguishable from a regular one.
The sole `reduced_gaussian` artifact in the repo,
`tests/conformance/grids/latlon/golden/reduced.json`, is a **grid-construction**
conformance golden (cell centres / neighbours / metrics), not an operator.

## 6. Missing primitives to cite for escalation

1. **Ragged longitude wrap** — an `index` form whose periodic modulus is a
   per-row quantity `nlon_per_row[lat]` (engine-level arrayop extension), and/or
2. an **inter-row interpolation** primitive for a convergent lon-gradient across
   rows of differing width, and
3. a **ragged latlon Layer-B runner** that consumes `nlon_per_row` (ESD-side),
   to exercise even the AST-expressible latitude derivative.

Until (1)+(3) land, reduced_gaussian operators are **reported infeasible**, not
implemented — consistent with the `discretizations/DIMENSIONAL_SPLIT_FFSL_INFEASIBILITY.md`
precedent ("no rule JSON is committed").
