"""
    arakawa_faq.jl — arakawa staggered-grid construction via the landed M1
    elementwise FAQ path (`EarthSciSerialization` AST evaluator).

Declarative companion: `discretizations/grids/arakawa/rules/arakawa_construction.esm`
expresses this construction as a semiring-FAQ document (RFC `semiring-faq-unified-ir`
§5.1/§5.2: cartesian interval index-sets + elementwise `arrayop` bodies — the M1
path, NO value-invention, NO recursion, all CONST cadence). It is the staggered-grid
sibling of the structured-grid FAQ template `cartesian_construction.esm`
(`src/cartesian_faq.jl`) and of the M3 value-invention DUO path
(`discretizations/grids/duo/rules/primal_topology.esm`, `src/topology_faq.jl`).

This module is the thin **evaluation bridge**: it shapes an `ArakawaGrid` over a
`CartesianBase` into the per-axis elementwise inputs and routes every piece of grid
ARITHMETIC — the affine spacings `dx`/`dy`, the two staggered 1-D coordinate maps,
and the cell-volume product — through ESS's AST evaluator (`eval_coeff`, the single
pathway in `src/rule_eval.jl`). ESD hosts no shadow evaluator (AGENTS.md single-pathway
invariant; GRIDS_API §4.3), so the determinism contract — and therefore the
cross-binding byte-identity of the output — lives entirely in `EarthSciSerialization`
(`CONFORMANCE_SPEC.md` §5.5). The only ESD-side logic is the *structural* part — the
cartesian-product flattening of the two axis maps into the four staggered locations,
the row-major neighbor linearization, the boundary masks, and the static A–E
variable-location lookup table — exactly as `cartesian_faq.jl` materializes the
per-cell gathers structurally while routing the arithmetic through ESS.

## The arakawa stagger transform is separable

Arakawa staggering over a cartesian base is the cartesian product of **two staggered
1-D affine maps per axis** (`src/grids/arakawa.jl`):

- a **center** map `center(i) = lo + (i - 0.5)*d` (cell centres), and
- a **node**   map `node(i)   = lo + (i - 1)*d`   (cell faces / corners).

The four staggered locations are then pure structural gathers of those two maps —
the same `(node, center)` factorisation the imperative accessors compute pointwise:

| location      | x source   | y source   | imperative accessor    |
|---------------|------------|------------|------------------------|
| `CellCenter`  | `center_x` | `center_y` | `arakawa_cell_center`  |
| `UEdge`       | `node_x`   | `center_y` | `arakawa_x_edge`       |
| `VEdge`       | `center_x` | `node_y`   | `arakawa_y_edge`       |
| `Corner`      | `node_x`   | `node_y`   | `arakawa_corner`       |

Because `node` is exactly the cartesian-template edge map `lo + (i-1)*d`, the ESS
evaluation is bit-identical to the imperative builder (the multiply/add order and the
operands match to the ULP); the `center` map adds only the `i - 0.5` half-cell offset.
"""

# AST nodes mirroring the equation bodies of
# `discretizations/grids/arakawa/rules/arakawa_construction.esm`. Held as module
# consts so the per-element evaluation does not re-allocate the trees; `eval_coeff`
# (src/rule_eval.jl) re-parses each call, the same single pathway the rule catalog uses.
const _ARK_D_NODE = Dict{String, Any}(   # spacing d = (hi - lo)/n
    "op" => "/",
    "args" => Any[Dict{String, Any}("op" => "-", "args" => Any["hi", "lo"]), "n"],
)
const _ARK_NODE_NODE = Dict{String, Any}(   # node/edge map:  lo + (i - 1)*d
    "op" => "+",
    "args" => Any[
        "lo",
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[Dict{String, Any}("op" => "-", "args" => Any["i", 1]), "d"],
        ),
    ],
)
const _ARK_CENTER_NODE = Dict{String, Any}(   # center map:  lo + (i - 0.5)*d
    "op" => "+",
    "args" => Any[
        "lo",
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[Dict{String, Any}("op" => "-", "args" => Any["i", 0.5]), "d"],
        ),
    ],
)
const _ARK_PROD_NODE = Dict{String, Any}("op" => "*", "args" => Any["a", "b"])

# Static variable-location codes mirroring the `.esm` lookup table
# (CellCenter=0, UEdge=1, VEdge=2, Corner=3). `_ARK_LOC_BY_CODE[code+1]` decodes.
const _ARK_LOC_BY_CODE = (CellCenter, UEdge, VEdge, Corner)

"""
    _ark_center_map(T, n, lo, d) -> Vector{T}

The `n` cell-centre coordinates of one axis via the M1 elementwise FAQ:
`center[i] = lo + (i - 0.5)*d` (`arakawa_construction.esm` `centers` equation),
evaluated by ESS (`eval_coeff`) so the result is bit-identical across bindings.
"""
function _ark_center_map(::Type{T}, n::Int, lo::Float64, d::Float64) where {T}
    out = Vector{T}(undef, n)
    @inbounds for i in 1:n
        out[i] = T(eval_coeff(_ARK_CENTER_NODE, Dict("lo" => lo, "i" => Float64(i), "d" => d)))
    end
    return out
end

"""
    _ark_node_map(T, np1, lo, d) -> Vector{T}

The `np1 = n+1` cell-face/corner coordinates of one axis via the M1 elementwise FAQ:
`node[i] = lo + (i - 1)*d` (`arakawa_construction.esm` `nodes` equation) — the same
affine map as the cartesian-template `edges`. Evaluated by ESS (`eval_coeff`).
"""
function _ark_node_map(::Type{T}, np1::Int, lo::Float64, d::Float64) where {T}
    out = Vector{T}(undef, np1)
    @inbounds for i in 1:np1
        out[i] = T(eval_coeff(_ARK_NODE_NODE, Dict("lo" => lo, "i" => Float64(i), "d" => d)))
    end
    return out
end

# Cartesian-product gather of two axis maps into one staggered location's flattened
# coordinate arrays. Flattened ROW-MAJOR — out[(j-1)*ni + i] — matching the imperative
# arakawa Tier-C layout (src/grids/arakawa.jl: slow axis = j). xmap has length ni
# (the location's first-dim size), ymap length nj.
function _ark_location_coords(
        ::Type{T}, xmap::Vector{T}, ymap::Vector{T},
        ni::Int, nj::Int
    ) where {T}
    xs = Vector{T}(undef, ni * nj)
    ys = Vector{T}(undef, ni * nj)
    @inbounds for j in 1:nj, i in 1:ni
        k = (j - 1) * ni + i
        xs[k] = xmap[i]
        ys[k] = ymap[j]
    end
    return xs, ys
end

# Row-major neighbor map over the (nx, ny) cell-centre set for one axis (`d` = 1 → x,
# 2 → y) and `offset` (∓1), with the 0-SENTINEL convention for off-grid — matching
# the imperative `neighbor_indices` (src/grids/arakawa.jl) exactly.
function _ark_axis_neighbor(nx::Int, ny::Int, d::Int, offset::Int)
    out = Vector{Int}(undef, nx * ny)
    @inbounds for j in 1:ny, i in 1:nx
        ii = d == 1 ? i + offset : i
        jj = d == 2 ? j + offset : j
        k = (j - 1) * nx + i
        out[k] = (1 <= ii <= nx && 1 <= jj <= ny) ? ((jj - 1) * nx + ii) : 0
    end
    return out
end

# Row-major boundary mask over the (nx, ny) cell-centre set for one axis and side —
# matching the imperative `boundary_mask` (src/grids/arakawa.jl).
function _ark_axis_boundary(nx::Int, ny::Int, d::Int, lower::Bool)
    out = falses(nx * ny)
    target = lower ? 1 : (d == 1 ? nx : ny)
    @inbounds for j in 1:ny, i in 1:nx
        v = d == 1 ? i : j
        out[(j - 1) * nx + i] = v == target
    end
    return out
end

# Static A–E variable-location table reconstructed from the integer codes (the
# `.esm` ifelse lookup), INDEPENDENT of the imperative `arakawa_variable_locations`
# so the conformance test is a genuine cross-check rather than an identity.
# h is always CellCenter; u/v per stagger: A→(CC,CC) B→(Corner,Corner) C→(UEdge,VEdge)
# D→(VEdge,UEdge) E→(Corner,Corner).
function _ark_variable_locations(s::ArakawaStagger)
    sc = s === ArakawaA ? 1 :
        s === ArakawaB ? 2 :
        s === ArakawaC ? 3 :
        s === ArakawaD ? 4 : 5
    h = 0
    u = sc == 1 ? 0 : sc == 2 ? 3 : sc == 3 ? 1 : sc == 4 ? 2 : 3
    v = sc == 1 ? 0 : sc == 2 ? 3 : sc == 3 ? 2 : sc == 4 ? 1 : 3
    return (_ARK_LOC_BY_CODE[h + 1], _ARK_LOC_BY_CODE[u + 1], _ARK_LOC_BY_CODE[v + 1])
end

# Static location-shape table over a (nx, ny) base — the "shape arithmetic" being
# declarativized out of src/staggering.jl, reconstructed here from the location code.
function _ark_location_shape(loc::VarLocation, nx::Int, ny::Int)
    return loc === CellCenter ? (nx, ny) :
        loc === UEdge ? (nx + 1, ny) :
        loc === VEdge ? (nx, ny + 1) :
        (nx + 1, ny + 1)   # Corner
end

"""
    arakawa_construction_faq(grid::ArakawaGrid) -> NamedTuple

Materialize an arakawa staggered grid's full construction — the four staggered
location coordinate arrays, the cell-centred widths/volume, the row-major neighbor
maps and boundary masks, and the static A–E variable-location/shape table — from the
declarative M1 elementwise FAQ, evaluating all arithmetic through the landed ESS AST
evaluator. Currently supports a `CartesianBase` (the conformance scope; lat-lon-base
composition has not landed in any binding). Returns a NamedTuple whose fields match
the imperative `src/grids/arakawa.jl` accessors so the conformance harness can assert
ULP/byte equality:

- `nx`, `ny` :: `Int`; `dx`, `dy` :: `T` — the base cell counts and affine spacings
  (`metric_eval(g, :dx/:dy)`).
- `center_x`/`center_y` :: `Vector{T}` and `node_x`/`node_y` :: `Vector{T}` — the two
  staggered 1-D axis maps (`length` `nx`/`ny` and `nx+1`/`ny+1`).
- `coords` :: `NamedTuple` `(cell_center, u_edge, v_edge, corner)`, each an `(xs, ys)`
  pair flattened row-major over that location's shape, matching `arakawa_cell_center`
  / `arakawa_x_edge` / `arakawa_y_edge` / `arakawa_corner`.
- `cell_widths` :: `NTuple{2,Vector{T}}` (per axis, `fill(dx)`/`fill(dy)`), `cell_volume`
  :: `Vector{T}` (`fill(dx*dy)`, the product routed through ESS) — over the `nx*ny`
  cell-centre set, matching `cell_widths(g, axis)` / `cell_volume(g)`.
- `neighbor_minus`/`neighbor_plus` :: `NTuple{2,Vector{Int}}` — row-major linear id of
  the ∓1 axis neighbor with the `0` sentinel, matching `neighbor_indices(g, axis, ∓1)`.
- `boundary_lower`/`boundary_upper` :: `NTuple{2,Vector{Bool}}` — matching
  `boundary_mask(g, axis, :lower/:upper)`.
- `stagger` :: `ArakawaStagger`, `rotated` :: `Bool` (E); `variable_locations` ::
  `(h, u, v)` of `VarLocation`, `variable_shapes` :: `(h, u, v)` of `(ni, nj)`,
  `location_shapes` :: `(cell_center, u_edge, v_edge, corner)` of `(ni, nj)` — the
  static A–E table, matching `arakawa_variable_locations` / `variable_shape` /
  `arakawa_shape`.

The coordinate/volume arithmetic rides ESS's determinism contract; the structural
arrays (the cartesian-product flattening into the four staggered locations, the
neighbor linearization, the masks, and the static location table) are pure index
logic, the staggered analogue of the per-cell gathers materialized structurally in
`cartesian_construction_faq`.
"""
function arakawa_construction_faq(grid::ArakawaGrid{T, <:CartesianBase}) where {T}
    b = grid.base
    nx = arakawa_nx(b)
    ny = arakawa_ny(b)
    xlo = Float64(b.xlo)
    xhi = Float64(b.xhi)
    ylo = Float64(b.ylo)
    yhi = Float64(b.yhi)

    # --- affine spacings dx/dy (routed through ESS) ------------------------------
    dx = eval_coeff(_ARK_D_NODE, Dict("hi" => xhi, "lo" => xlo, "n" => Float64(nx)))
    dy = eval_coeff(_ARK_D_NODE, Dict("hi" => yhi, "lo" => ylo, "n" => Float64(ny)))

    # --- the two staggered 1-D axis maps (elementwise FAQ) ----------------------
    center_x = _ark_center_map(T, nx, xlo, dx)
    center_y = _ark_center_map(T, ny, ylo, dy)
    node_x = _ark_node_map(T, nx + 1, xlo, dx)
    node_y = _ark_node_map(T, ny + 1, ylo, dy)

    # --- the four staggered locations (cartesian-product gathers, row-major) -----
    cc_x, cc_y = _ark_location_coords(T, center_x, center_y, nx, ny)
    ue_x, ue_y = _ark_location_coords(T, node_x, center_y, nx + 1, ny)
    ve_x, ve_y = _ark_location_coords(T, center_x, node_y, nx, ny + 1)
    co_x, co_y = _ark_location_coords(T, node_x, node_y, nx + 1, ny + 1)

    # --- cell-centred widths / volume -------------------------------------------
    nc = nx * ny
    cw_x = fill(T(dx), nc)
    cw_y = fill(T(dy), nc)
    cell_volume = fill(T(eval_coeff(_ARK_PROD_NODE, Dict("a" => dx, "b" => dy))), nc)

    # --- neighbor maps (row-major linear ids, 0-sentinel) + boundary masks -------
    neighbor_minus = (_ark_axis_neighbor(nx, ny, 1, -1), _ark_axis_neighbor(nx, ny, 2, -1))
    neighbor_plus = (_ark_axis_neighbor(nx, ny, 1, +1), _ark_axis_neighbor(nx, ny, 2, +1))
    boundary_lower = (_ark_axis_boundary(nx, ny, 1, true), _ark_axis_boundary(nx, ny, 2, true))
    boundary_upper = (_ark_axis_boundary(nx, ny, 1, false), _ark_axis_boundary(nx, ny, 2, false))

    # --- static A–E variable-location / shape table -----------------------------
    s = grid.stagger
    h_loc, u_loc, v_loc = _ark_variable_locations(s)

    return (
        nx = nx,
        ny = ny,
        dx = T(dx),
        dy = T(dy),
        center_x = center_x,
        center_y = center_y,
        node_x = node_x,
        node_y = node_y,
        coords = (
            cell_center = (cc_x, cc_y),
            u_edge = (ue_x, ue_y),
            v_edge = (ve_x, ve_y),
            corner = (co_x, co_y),
        ),
        cell_widths = (cw_x, cw_y),
        cell_volume = cell_volume,
        neighbor_minus = neighbor_minus,
        neighbor_plus = neighbor_plus,
        boundary_lower = boundary_lower,
        boundary_upper = boundary_upper,
        stagger = s,
        rotated = (s === ArakawaE),
        variable_locations = (h = h_loc, u = u_loc, v = v_loc),
        variable_shapes = (
            h = _ark_location_shape(h_loc, nx, ny),
            u = _ark_location_shape(u_loc, nx, ny),
            v = _ark_location_shape(v_loc, nx, ny),
        ),
        location_shapes = (
            cell_center = _ark_location_shape(CellCenter, nx, ny),
            u_edge = _ark_location_shape(UEdge, nx, ny),
            v_edge = _ark_location_shape(VEdge, nx, ny),
            corner = _ark_location_shape(Corner, nx, ny),
        ),
    )
end
