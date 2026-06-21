"""
    latlon_faq.jl — lat-lon structured-grid construction via the landed M1
    elementwise FAQ path (`EarthSciSerialization` AST evaluator).

Declarative companion: `discretizations/grids/latlon/rules/latlon_construction.esm`
expresses this construction as a semiring-FAQ document (RFC `semiring-faq-unified-ir`
§5.1/§5.2: interval index-sets + elementwise `arrayop` bodies — the M1 path, NO
value-invention, all CONST cadence). It is the lat-lon instance of the structured-grid
FAQ template established for cartesian in `discretizations/grids/cartesian/rules/
cartesian_construction.esm` (`src/cartesian_faq.jl`), and the rectilinear-sphere
counterpart to the DUO value-invention path in `src/topology_faq.jl`.

This module is the thin **evaluation bridge**: it shapes a `LatLonGrid` into the
per-row elementwise inputs and routes every piece of grid ARITHMETIC — the per-row
affine longitude map (`dlon = 2π/nlon`, `lon = lon_start + (i-½)·dlon`), the
spherical-rectangle cell area `R²·dlon·(sin φ_n − sin φ_s)`, and the closed-form
curvilinear metric (`g_λλ = R²cos²φ`, `g_φφ = R²`, `ginv`, the Jacobian
`R²|cos φ|`, and the lone non-vanishing derivative `∂g_λλ/∂φ = −2R²cos φ sin φ`) —
through ESS's AST evaluator (`eval_coeff`, the single pathway in `src/rule_eval.jl`,
which dispatches `sin`/`cos`/`abs`/`±`/`*`/`/` to `EarthSciSerialization.evaluate_expr`).
ESD hosts no shadow evaluator (AGENTS.md single-pathway invariant; GRIDS_API §4.3), so
the determinism contract — and therefore cross-binding byte-identity of the output —
lives entirely in `EarthSciSerialization` (`CONFORMANCE_SPEC.md` §5.5).

The only ESD-side logic is the *structural* part: the ragged-row-major flattening, the
periodic-longitude neighbor linearization, the reduced-Gaussian nearest-center N/S row
remap (`_latlon_map_i`, the relational rank/join materialized structurally), the
boundary masks, and the constant identity coordinate Jacobian — exactly as
`topology_faq.jl` materializes `cell_neighbors` structurally while routing the
arithmetic through ESS, and as `cartesian_faq.jl` keeps neighbor linearization /
masks as pure index logic. Supplied `lat_edges` / `lat_centers` are consumed as DATA
(the reduced-Gaussian-as-data branch); only their elementwise difference (the lat cell
width) and the sin/cos of the latitudes ride the ESS arithmetic path.
"""

# AST nodes mirroring the equation bodies of
# `discretizations/grids/latlon/rules/latlon_construction.esm`. Held as module consts
# so the per-row evaluation does not re-allocate the trees; `eval_coeff`
# (src/rule_eval.jl) re-parses each call, the same single pathway the rule catalog
# uses. π is inlined as the binary64 literal 2π = 6.283185307179586, identical to the
# imperative `T(2) * T(pi)` (an exact doubling), so `dlon = 2π/n` is bit-for-bit equal.
const _LL_DLON_NODE = Dict{String, Any}(
    "op" => "/", "args" => Any[6.283185307179586, "n"],
)
const _LL_LON_CENTER_NODE = Dict{String, Any}(
    "op" => "+",
    "args" => Any[
        "lon_start",
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[Dict{String, Any}("op" => "-", "args" => Any["i", 0.5]), "dlon"],
        ),
    ],
)
const _LL_LAT_WIDTH_NODE = Dict{String, Any}("op" => "-", "args" => Any["e_hi", "e_lo"])
# area = R*R*dlon*(sin(lat_n) - sin(lat_s)), left-assoc as ((R*R)*dlon)*(Δsin).
const _LL_AREA_NODE = Dict{String, Any}(
    "op" => "*",
    "args" => Any[
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[Dict{String, Any}("op" => "*", "args" => Any["R", "R"]), "dlon"],
        ),
        Dict{String, Any}(
            "op" => "-",
            "args" => Any[
                Dict{String, Any}("op" => "sin", "args" => Any["lat_n"]),
                Dict{String, Any}("op" => "sin", "args" => Any["lat_s"]),
            ],
        ),
    ],
)
const _LL_COS_NODE = Dict{String, Any}("op" => "cos", "args" => Any["lat"])
# g_λλ = R*R*cos_lat*cos_lat (cos evaluated once, reused — matches imperative).
const _LL_GLL_NODE = Dict{String, Any}(
    "op" => "*",
    "args" => Any[
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[Dict{String, Any}("op" => "*", "args" => Any["R", "R"]), "cos_lat"],
        ),
        "cos_lat",
    ],
)
const _LL_GPP_NODE = Dict{String, Any}("op" => "*", "args" => Any["R", "R"])
const _LL_INV_NODE = Dict{String, Any}("op" => "/", "args" => Any[1.0, "x"])
# jacobian = R*R*abs(cos(lat)), left-assoc as (R*R)*abs(cos).
const _LL_JAC_NODE = Dict{String, Any}(
    "op" => "*",
    "args" => Any[
        Dict{String, Any}("op" => "*", "args" => Any["R", "R"]),
        Dict{String, Any}(
            "op" => "abs",
            "args" => Any[Dict{String, Any}("op" => "cos", "args" => Any["lat"])],
        ),
    ],
)
# ∂g_λλ/∂φ = -2*R*R*cos(φ)*sin(φ), left-assoc as ((((-2)*R)*R)*cos)*sin.
const _LL_DGLL_NODE = Dict{String, Any}(
    "op" => "*",
    "args" => Any[
        Dict{String, Any}(
            "op" => "*",
            "args" => Any[
                Dict{String, Any}(
                    "op" => "*",
                    "args" => Any[
                        Dict{String, Any}("op" => "*", "args" => Any[-2.0, "R"]), "R",
                    ],
                ),
                Dict{String, Any}("op" => "cos", "args" => Any["phi"]),
            ],
        ),
        Dict{String, Any}("op" => "sin", "args" => Any["phi"]),
    ],
)

"""
    latlon_construction_faq(g::LatLonGrid) -> NamedTuple

Materialize a lat-lon grid's full construction — coordinates, the closed-form
curvilinear metric, periodic / nearest-center neighbor maps and boundary masks —
from the declarative M1 elementwise FAQ, evaluating all arithmetic through the landed
ESS AST evaluator. Returns a NamedTuple whose fields match the imperative
`LatLonGrid` trait arrays (`src/grids/latlon.jl`) so the conformance harness can
assert byte/ULP equality. All per-cell arrays use the flat ragged-row-major layout
(cell `(j, i)` at flat index `row_offset(g, j) + i`):

- `cell_centers_lon`, `cell_centers_lat` :: `Vector{T}` — per-cell `(lon, lat)`,
  matching `cell_centers(g, :lon)` / `cell_centers(g, :lat)`. `lon` rides the per-row
  affine map; `lat` is the supplied/derived row-center DATA.
- `cell_widths_lon`, `cell_widths_lat` :: `Vector{T}` — `dlon` and the lat-edge
  difference, matching `cell_widths(g, :lon)` / `cell_widths(g, :lat)`.
- `cell_volume` :: `Vector{T}` — spherical-rectangle area `R²·dlon·(sin φ_n − sin φ_s)`,
  matching `cell_volume(g)`.
- `lat_edges`, `lat_centers` :: `Vector{T}` — the 1-D latitude arrays consumed as data.
- `metric_g`, `metric_ginv`, `coord_jacobian` :: `Array{T,3}` `(nc, 2, 2)`;
  `metric_jacobian` :: `Vector{T}`; `metric_dgij_dxk`, `coord_jacobian_second` ::
  `Array{T,4}` `(nc, 2, 2, 2)` — matching `metric_g(g)` / `metric_ginv(g)` /
  `coord_jacobian(g, :lon_lat)` / `metric_jacobian(g)` / `metric_dgij_dxk(g)` /
  `coord_jacobian_second(g, :lon_lat)`.
- `neighbor_lon_minus`, `neighbor_lon_plus`, `neighbor_lat_minus`, `neighbor_lat_plus`
  :: `Vector{Int}` — flat linear id of the offset ∓1 neighbor with the `0` sentinel
  (longitude wraps periodically; latitude uses the reduced-Gaussian nearest-center
  remap and `0` at the poles), matching `neighbor_indices(g, axis, ∓1)`.
- `boundary_lat_lower`, `boundary_lat_upper`, `boundary_lon_lower`,
  `boundary_lon_upper` :: `Vector{Bool}` — first/last latitude row (longitude wraps,
  so its masks are all-false), matching `boundary_mask(g, axis, side)`.

The coordinate/area/metric arithmetic rides ESS's determinism contract; the structural
arrays (ragged flattening, periodic / nearest-center neighbor linearization, masks,
identity coordinate Jacobian) are pure index logic, the lat-lon analogue of the
inversion joins materialized structurally in `primal_topology_faq`.
"""
function latlon_construction_faq(g::LatLonGrid{T}) where {T}
    nlat = g.nlat
    per_row = g.nlon_per_row
    nc = n_cells(g)
    Rf = Float64(g.R)
    lon_start_f = Float64(g.lon_start)

    # Row offsets (0-based flat start of each row), precomputed once.
    row_off = Vector{Int}(undef, nlat + 1)
    row_off[1] = 0
    @inbounds for j in 1:nlat
        row_off[j + 1] = row_off[j] + per_row[j]
    end

    cc_lon = Vector{T}(undef, nc)
    cc_lat = Vector{T}(undef, nc)
    cw_lon = Vector{T}(undef, nc)
    cw_lat = Vector{T}(undef, nc)
    cell_vol = Vector{T}(undef, nc)
    metric_g_ = zeros(T, nc, 2, 2)
    metric_ginv_ = zeros(T, nc, 2, 2)
    metric_jac = Vector{T}(undef, nc)
    metric_dg = zeros(T, nc, 2, 2, 2)

    # --- per-row elementwise FAQ (affine lon, trig metric, area) ----------------
    @inbounds for j in 1:nlat
        n_i = per_row[j]
        base = row_off[j]

        # Longitude affine spacing dlon = 2π/n (ESS).
        dlon = eval_coeff(_LL_DLON_NODE, Dict("n" => Float64(n_i)))

        # Latitude edges/centers consumed as DATA; only the width difference is ESS.
        lat_c = g.lat_centers[j]
        lat_s = Float64(g.lat_edges[j])
        lat_n = Float64(g.lat_edges[j + 1])
        lat_w = eval_coeff(_LL_LAT_WIDTH_NODE, Dict("e_hi" => lat_n, "e_lo" => lat_s))

        # Spherical-rectangle cell area (ESS; sin routed through the evaluator).
        area = eval_coeff(_LL_AREA_NODE,
            Dict("R" => Rf, "dlon" => dlon, "lat_n" => lat_n, "lat_s" => lat_s))

        # Closed-form curvilinear metric — row-constant (depends only on φ, R).
        lat_cf = Float64(lat_c)
        cos_lat = eval_coeff(_LL_COS_NODE, Dict("lat" => lat_cf))
        g_ll = eval_coeff(_LL_GLL_NODE, Dict("R" => Rf, "cos_lat" => cos_lat))
        g_pp = eval_coeff(_LL_GPP_NODE, Dict("R" => Rf))
        inv_ll = g_ll > 0.0 ? eval_coeff(_LL_INV_NODE, Dict("x" => g_ll)) : Inf
        inv_pp = eval_coeff(_LL_INV_NODE, Dict("x" => g_pp))
        jac = eval_coeff(_LL_JAC_NODE, Dict("R" => Rf, "lat" => lat_cf))
        dgll = eval_coeff(_LL_DGLL_NODE, Dict("R" => Rf, "phi" => lat_cf))

        for i in 1:n_i
            k = base + i
            lon = eval_coeff(_LL_LON_CENTER_NODE,
                Dict("lon_start" => lon_start_f, "i" => Float64(i), "dlon" => dlon))
            cc_lon[k] = T(lon)
            cc_lat[k] = lat_c           # supplied/derived row-center data, gathered as-is
            cw_lon[k] = T(dlon)
            cw_lat[k] = T(lat_w)
            cell_vol[k] = T(area)
            metric_g_[k, 1, 1] = T(g_ll)
            metric_g_[k, 2, 2] = T(g_pp)
            metric_ginv_[k, 1, 1] = T(inv_ll)
            metric_ginv_[k, 2, 2] = T(inv_pp)
            metric_jac[k] = T(jac)
            metric_dg[k, 1, 1, 2] = T(dgll)   # ∂g_λλ/∂φ
        end
    end

    # --- identity coordinate Jacobian (computational axes == (lon, lat)) --------
    coord_jac = zeros(T, nc, 2, 2)
    @inbounds for k in 1:nc
        coord_jac[k, 1, 1] = one(T)
        coord_jac[k, 2, 2] = one(T)
    end
    coord_jac_second = zeros(T, nc, 2, 2, 2)

    # --- neighbor maps (0-sentinel, flat ragged-row-major linear ids) -----------
    function _nbr_lon(offset::Int)
        out = Vector{Int}(undef, nc)
        @inbounds for j in 1:nlat
            n_i = per_row[j]
            base = row_off[j]
            for i in 1:n_i
                ii = mod(i - 1 + offset, n_i) + 1   # periodic longitude wrap
                out[base + i] = base + ii
            end
        end
        return out
    end
    function _nbr_lat(offset::Int)
        out = Vector{Int}(undef, nc)
        @inbounds for j in 1:nlat
            n_i = per_row[j]
            base = row_off[j]
            jj = j + offset
            if 1 <= jj <= nlat
                n_t = per_row[jj]
                base_t = row_off[jj]
                for i in 1:n_i
                    i_t = _latlon_map_i(i, n_i, n_t)   # nearest-center rank/join
                    out[base + i] = base_t + i_t
                end
            else
                for i in 1:n_i
                    out[base + i] = 0                  # pole: pole_policy=:none
                end
            end
        end
        return out
    end

    # --- boundary masks (lon wraps → all-false; lat marks first/last row) --------
    function _bd_lat(target::Int)
        out = Vector{Bool}(undef, nc)
        @inbounds for j in 1:nlat
            n_i = per_row[j]
            base = row_off[j]
            mark = j == target
            for i in 1:n_i
                out[base + i] = mark
            end
        end
        return out
    end

    return (
        cell_centers_lon = cc_lon,
        cell_centers_lat = cc_lat,
        cell_widths_lon = cw_lon,
        cell_widths_lat = cw_lat,
        cell_volume = cell_vol,
        lat_edges = copy(g.lat_edges),
        lat_centers = copy(g.lat_centers),
        metric_g = metric_g_,
        metric_ginv = metric_ginv_,
        metric_jacobian = metric_jac,
        metric_dgij_dxk = metric_dg,
        coord_jacobian = coord_jac,
        coord_jacobian_second = coord_jac_second,
        neighbor_lon_minus = _nbr_lon(-1),
        neighbor_lon_plus = _nbr_lon(+1),
        neighbor_lat_minus = _nbr_lat(-1),
        neighbor_lat_plus = _nbr_lat(+1),
        boundary_lat_lower = _bd_lat(1),
        boundary_lat_upper = _bd_lat(nlat),
        boundary_lon_lower = fill(false, nc),
        boundary_lon_upper = fill(false, nc),
    )
end
