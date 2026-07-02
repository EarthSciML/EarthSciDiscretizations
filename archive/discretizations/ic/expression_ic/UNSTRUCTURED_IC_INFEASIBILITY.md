# `expression_ic` on unstructured grids (duo / mpas) — DECLARATIVE-OR-FAIL verdict

**Bead:** esd-6g4.12 (G13) · **Epic:** esd-6g4 · **Label:** `operator-grid-coverage`
**Oracle source pin:** EarthSciSerialization resolved via ESD `Project.toml`
`[sources]` = `{rev = "main", subdir = "packages/EarthSciSerialization.jl"}`,
`Manifest.toml` `git-tree-sha1 = f2b411a6617fdb85539f2abddec7bb4523149d19`
(ESS `main` HEAD `f137bf74`). All ESS line numbers below are at that revision.

---

## 0. Verdict (TL;DR)

> **latlon + arakawa: FEASIBLE — authored.** **duo + mpas: DECLARATIVE-INFEASIBLE
> over the existing ESS engine.**

The bead asks for IC fixtures on `latlon`, `arakawa`, `duo`, and `mpas`. Two of
the four are delivered as a Layer-A canonical byte golden
(`fixtures/canonical/`, one document carrying both a latlon 8×8 model and an
arakawa 4×4 model). The other two **cannot** be exercised through `expression_ic`
without new ESS engine support, because the rule's materialization model is
**structured-grid-only**. Per the bead's `DECLARATIVE-OR-FAIL` contract, the
correct outcome for duo/mpas is this precise gap report + escalation — **no
forced fixture and no hand-coded IC executor is committed.**

This is **not** the same class of gap as a deferred scheme kind (contrast
esd-6g4.6's dimensional_split/FFSL verdict, which is *temporal*). The IC arrayop
*does* materialize; what is missing is a **coordinate model for the unstructured
cell layout axis**.

---

## 1. What `expression_ic` actually does

`expression_ic.json` declares (verbatim):

> "The engine substitutes each spatial-dim symbol with the cell-centre
> coordinate: x→coord_x($i), y→coord_y($j), z→coord_z($k), etc. (uniform:
> (idx - 0.5)*h; nonuniform: via coord_x const_array)."

So the contract is explicitly **per-spatial-dimension**: every loop index is a
grid *dimension*, and each dimension symbol is replaced by a cell-centre
coordinate `index(coord_<dim>, idx)`. The materialization is owned by the ESS
engine (`expression_ic.json`: "ESS `ic-arrayop-engine` materializes the ic
dispatch"), specifically `EarthSciSerialization._try_materialize_ic_arrayop!`
(`src/discretize.jl:1052`). Its essential structure:

```
dim_sizes = get(gmeta, "dim_sizes", nothing)
dim_sizes isa AbstractDict || return false          # discretize.jl:1071
for d in 1:nd                                        # nd = length(shape)
    dim_name = String(shape[d])
    idx = _ARRAYOP_INDEX_NAMES[d]                    # i, j, k, …
    ranges[idx] = [1, dim_sizes[dim_name]]
    coord_subst[dim_name] = index(coord_<dim_name>, idx)   # discretize.jl:1095
end
```

Every index in the materialized arrayop is a **grid dimension** with a known
integer size, and the substituted coordinate is `coord_<dim>` — a cell-centre
position array keyed by *that dimension*. This is exactly right for any grid
whose dimensions **are** spatial axes.

## 2. Why latlon + arakawa are FEASIBLE

Both are structured: their `dimensions` are spatial axes (`lon`/`lat`,
`x`/`y`). The authored fixture's two models discretize cleanly (verified end to
end through the walker's `apply_rule_and_diff`):

| model (grid)     | IC RHS                          | materialized arrayop                                              |
|------------------|---------------------------------|------------------------------------------------------------------|
| Mlatlon (latlon) | `sin(lon)·cos(lat)`             | `output_idx=[i,j]`, `ranges={i:[1,8],j:[1,8]}`, `sin(coord_lon[i])·cos(coord_lat[j])` |
| Marakawa (arakawa) | `sin(2π·x)·sin(2π·y)`          | `output_idx=[i,j]`, `ranges={i:[1,4],j:[1,4]}`, `sin(2π·coord_x[i])·sin(2π·coord_y[j])` |

`coord_lon`/`coord_lat`/`coord_x`/`coord_y` are the cell-centre const_arrays the
host injects at bind time (ESD `_prepare_expression_ics!`,
`src/ode_problem.jl`). This is the same path the existing cartesian + vertical
goldens use (`test/test_bc_ic_goldens.jl`), now extended to two more families.

## 3. Why duo + mpas are INFEASIBLE through `expression_ic`

### 3.1 The layout axis is not a spatial coordinate

DUO and MPAS are unstructured. Both expose a **flat 1-D `cell` layout axis**;
physical coordinates are **separate per-cell arrays**, not dimensions:

- `src/grids/mpas.jl`: *"MPAS uses a flat 1-D cell axis (`:cell`); spatial
  coordinates are exposed via `cell_centers(g, :lon)` / `(g, :lat)`"* and the
  accessor throws *"`:cell` is the layout axis"* if asked to treat `cell` as a
  coordinate (`src/grids/mpas.jl:843`).
- `src/grids/duo.jl`: same — *"duo: cell_centers needs a coordinate axis (:lon
  or :lat); :cell is the layout axis"* (`src/grids/duo.jl:553`); `lon`, `lat`
  are `(Nc,)` cell-centre arrays.

`expression_ic` substitutes the *dimension name*. On these grids the only
"dimension" available to a state field is `cell`, and `coord_cell` is a
layout-index "coordinate" — **not** lon/lat. The substitution the rule performs
therefore has no physical meaning here.

### 3.2 The standard declaration doesn't even materialize

The §6 declarative configs for these families carry **no
`dimensions:[{name,size}]` array** — only scalar counts (`n_cells`, `n_edges`,
…) plus a loader ref (`discretizations/grids/mpas/x1.642.esm`,
`discretizations/grids/duo/icos_level0.esm`). With no `dim_sizes` dict for the
`cell` axis, `_try_materialize_ic_arrayop!` hits `dim_sizes isa AbstractDict ||
return false` (`discretize.jl:1071`) and the IC equation is left **un-lifted**.
Empirically, `discretize` on a duo/mpas model with `initialization_equations:
[{lhs:u, rhs:sin(cell)}]` returns the RHS **unchanged** (`{"op":"sin","args":
["cell"]}`) — no arrayop, no coordinate binding.

### 3.3 Forcing it produces a physically wrong contract

Injecting a synthetic `dimensions:[{name:"cell",size:N}]` *does* make the
materializer fire, but it emits `cell → index(coord_cell, i)` — i.e. it treats
the **cell index** as a spatial coordinate sampled at `(i-0.5)*h`. That is not
the cell's lon/lat; it is meaningless as an initial condition. A genuine
unstructured IC (`u[cell] = f(lon[cell], lat[cell])`) would need the per-cell
`lon`/`lat` const_arrays indexed by the loop variable — a **different**
mechanism than the structured `coord_<dim>` substitution, and one the rule does
not describe.

### 3.4 The repo already handles unstructured ICs a different way

The existing unstructured operator fixtures confirm this split. The DUO heat-
equation PDE fixture leaves its expression IC **empty** and documents that the
IC is injected numerically by the runner, not by `discretize`
(`discretizations/finite_difference/nn_diffusion_duo/fixtures/integration/duo_nn_diffusion_pde.esm`):

> "MMS IC u(lon,lat)=sin(lon)*cos(lat) requires DUO cell-center coordinate
> access. Injected by the GRIDS-F3 helper (run_mms_convergence_case) at bind
> time … placeholder 'default: 0.0' above is overridden by the runner."

So in this codebase, unstructured ICs are a **bind-time numerical sampling**
concern (the runner reads per-cell lon/lat and fills `u0`), categorically
distinct from the `expression_ic` discretize→arrayop path that structured grids
use.

## 4. What closing the gap would require (NOT done here)

A real duo/mpas `expression_ic` path needs **GENERIC ESS engine support**: an
unstructured branch of the ic-arrayop materializer that

1. derives `ranges = [1, n_cells]` from the unstructured cell count (not a
   `dimensions` entry), and
2. binds spatial symbols in the IC RHS to **per-cell coordinate const_arrays**
   (`lon`/`lat` indexed by the cell loop variable), rather than the structured
   `coord_<dim>` cell-centre arrays.

That is an ESS engine change with an ESD const_arrays binding contract (cf. the
`nn_diffusion` unstructured primitive-array path), not a fixture. It is out of
scope for a fixtures-only bead and is exactly the imperative/engine work
`DECLARATIVE-OR-FAIL` says to escalate rather than hand-roll.

## 5. Resolution

- **Authored:** `fixtures/canonical/{input,expected}.esm` (latlon + arakawa) +
  `fixtures/canonical/regenerate_golden.jl`; walker Layer-A registers a PASS for
  `expression_ic` (`test/test_esd_walker.jl`).
- **Escalated:** duo/mpas `expression_ic` support → epic esd-6g4 drain (mayor
  decision: ESS engine work + ESD const_arrays binding, or accept that
  unstructured ICs stay a runner-injected numerical concern).
