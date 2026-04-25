"""
Cubed-sphere grid: `CubedSphereGrid(Nc; R=1.0, Ng=3)`.
"""

struct CubedSphereGrid{T} <: AbstractCubedSphereGrid
    Nc::Int; Ng::Int; R::T
    ξ_centers::Vector{T}; η_centers::Vector{T}
    ξ_edges::Vector{T}; η_edges::Vector{T}
    lon::Array{T, 3}; lat::Array{T, 3}
    area::Array{T, 3}; dx::Array{T, 3}; dy::Array{T, 3}
    dξ::T; dη::T
    rotation_angles::Dict{Tuple{Int, EdgeDirection}, Float64}
    # 2×2 rotation matrices for vector field ghost cell filling:
    # (M11, M12, M21, M22) transforming neighbor (u_ξ, u_η) to local basis
    rotation_matrices::Dict{Tuple{Int, EdgeDirection}, NTuple{4, Float64}}
    # Metric tensor components at cell centers
    J::Array{T, 3}           # Jacobian
    ginv_ξξ::Array{T, 3}     # Inverse metric g^{ξξ}
    ginv_ηη::Array{T, 3}     # Inverse metric g^{ηη}
    ginv_ξη::Array{T, 3}     # Inverse metric g^{ξη}
    # Physical-to-computational coordinate Jacobian: d(ξ,η)/d(lon,lat)
    dξ_dlon::Array{T, 3}
    dξ_dlat::Array{T, 3}
    dη_dlon::Array{T, 3}
    dη_dlat::Array{T, 3}
    # Second-derivative coordinate Jacobian: d²(ξ,η)/d(lon,lat)²
    d2ξ_dlon2::Array{T, 3}
    d2ξ_dlondlat::Array{T, 3}
    d2ξ_dlat2::Array{T, 3}
    d2η_dlon2::Array{T, 3}
    d2η_dlondlat::Array{T, 3}
    d2η_dlat2::Array{T, 3}
    # Derivatives of J·g^{ab} for the Laplacian metric corrections
    dJgxe_dξ::Array{T, 3}   # ∂(J·g^{ξη})/∂ξ at cell centers
    dJgxe_dη::Array{T, 3}   # ∂(J·g^{ξη})/∂η at cell centers
    dJgxx_dξ::Array{T, 3}   # ∂(J·g^{ξξ})/∂ξ at cell centers
    dJgyy_dη::Array{T, 3}   # ∂(J·g^{ηη})/∂η at cell centers
    # Center-to-center physical distances
    dist_xi::Array{T, 3}     # (6, Nc-1, Nc): distance between cell (i,j) and (i+1,j)
    dist_eta::Array{T, 3}    # (6, Nc, Nc-1): distance between cell (i,j) and (i,j+1)
    # Boundary distances: cross-panel center-to-center distances at panel edges
    dist_xi_bnd::Array{T, 3} # (6, 2, Nc): [p, 1=west/2=east, j] boundary distances in ξ
    dist_eta_bnd::Array{T, 3} # (6, Nc, 2): [p, i, 1=south/2=north] boundary distances in η
    # FV3 super-grid angular values: sin/cos of angle between e_ξ and e_η
    # at 9 sub-positions per cell (see super_grid.jl for position numbering)
    sin_sg::Array{T, 4}      # (6, Nc, Nc, 9)
    cos_sg::Array{T, 4}      # (6, Nc, Nc, 9)
end

function CubedSphereGrid(Nc::Int; R = 1.0, Ng::Int = 3)
    T = typeof(float(R))
    dξ = T(π / 2) / Nc; dη = dξ
    ξ_edges = [T(-π / 4) + (i - 1) * dξ for i in 1:(Nc + 1)]
    η_edges = [T(-π / 4) + (j - 1) * dη for j in 1:(Nc + 1)]
    ξ_centers = [(ξ_edges[i] + ξ_edges[i + 1]) / 2 for i in 1:Nc]
    η_centers = [(η_edges[j] + η_edges[j + 1]) / 2 for j in 1:Nc]

    lon = zeros(T, 6, Nc, Nc); lat = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        lo, la = gnomonic_to_lonlat(ξ_centers[i], η_centers[j], p)
        lon[p, i, j] = lo; lat[p, i, j] = la
    end

    area = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        area[p, i, j] = compute_cell_area((ξ_edges[i], ξ_edges[i + 1]), (η_edges[j], η_edges[j + 1]), R, p)
    end

    dx = zeros(T, 6, Nc + 1, Nc)
    for p in 1:6, i in 1:(Nc + 1), j in 1:Nc
        dx[p, i, j] = compute_edge_length(ξ_edges[i], η_edges[j], ξ_edges[i], η_edges[j + 1], R, p)
    end
    dy = zeros(T, 6, Nc, Nc + 1)
    for p in 1:6, i in 1:Nc, j in 1:(Nc + 1)
        dy[p, i, j] = compute_edge_length(ξ_edges[i], η_edges[j], ξ_edges[i + 1], η_edges[j], R, p)
    end

    rotation_angles = Dict{Tuple{Int, EdgeDirection}, Float64}()
    rotation_matrices = Dict{Tuple{Int, EdgeDirection}, NTuple{4, Float64}}()
    for p in 1:6, dir in (West, East, South, North)
        nb = PANEL_CONNECTIVITY[p][dir]
        rotation_angles[(p, dir)] = compute_rotation_angle(p, dir, nb.neighbor_panel, nb.neighbor_edge)
        rotation_matrices[(p, dir)] = compute_edge_rotation_matrix(p, dir)
    end

    # Precompute metric tensor and its inverse at cell centers
    J_arr = zeros(T, 6, Nc, Nc)
    ginv_ξξ = zeros(T, 6, Nc, Nc)
    ginv_ηη = zeros(T, 6, Nc, Nc)
    ginv_ξη = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        Jv, gxx, gee, gxe = gnomonic_metric(ξ_centers[i], η_centers[j], R)
        J_arr[p, i, j] = Jv
        det_g = gxx * gee - gxe^2
        ginv_ξξ[p, i, j] = gee / det_g
        ginv_ηη[p, i, j] = gxx / det_g
        ginv_ξη[p, i, j] = -gxe / det_g
    end

    # Precompute physical-to-computational coordinate Jacobian
    dξ_dlon_arr = zeros(T, 6, Nc, Nc)
    dξ_dlat_arr = zeros(T, 6, Nc, Nc)
    dη_dlon_arr = zeros(T, 6, Nc, Nc)
    dη_dlat_arr = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        jac = compute_coord_jacobian(ξ_centers[i], η_centers[j], p)
        dξ_dlon_arr[p, i, j] = jac.dξ_dlon
        dξ_dlat_arr[p, i, j] = jac.dξ_dlat
        dη_dlon_arr[p, i, j] = jac.dη_dlon
        dη_dlat_arr[p, i, j] = jac.dη_dlat
    end

    # Precompute second-derivative coordinate Jacobian
    d2ξ_dlon2_arr = zeros(T, 6, Nc, Nc)
    d2ξ_dlondlat_arr = zeros(T, 6, Nc, Nc)
    d2ξ_dlat2_arr = zeros(T, 6, Nc, Nc)
    d2η_dlon2_arr = zeros(T, 6, Nc, Nc)
    d2η_dlondlat_arr = zeros(T, 6, Nc, Nc)
    d2η_dlat2_arr = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        jac2 = compute_second_coord_jacobian(ξ_centers[i], η_centers[j], p)
        d2ξ_dlon2_arr[p, i, j] = jac2.d2ξ_dlon2
        d2ξ_dlondlat_arr[p, i, j] = jac2.d2ξ_dlondlat
        d2ξ_dlat2_arr[p, i, j] = jac2.d2ξ_dlat2
        d2η_dlon2_arr[p, i, j] = jac2.d2η_dlon2
        d2η_dlondlat_arr[p, i, j] = jac2.d2η_dlondlat
        d2η_dlat2_arr[p, i, j] = jac2.d2η_dlat2
    end

    # Precompute ∂(J·g^{ξη})/∂ξ and ∂(J·g^{ξη})/∂η for the Laplacian cross-metric correction.
    # Uses analytical gnomonic_metric at ξ±h and η±h stencil points.
    dJgxe_dξ_arr = zeros(T, 6, Nc, Nc)
    dJgxe_dη_arr = zeros(T, 6, Nc, Nc)
    h_metric = T(1.0e-7)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        ξc = ξ_centers[i]; ηc = η_centers[j]
        # Jg^{ξη} at (ξ±h, η) for ∂/∂ξ
        Jp, gxx_p, gee_p, gxe_p = gnomonic_metric(ξc + h_metric, ηc, R)
        Jm, gxx_m, gee_m, gxe_m = gnomonic_metric(ξc - h_metric, ηc, R)
        det_gp = gxx_p * gee_p - gxe_p^2
        det_gm = gxx_m * gee_m - gxe_m^2
        Jgxe_p = Jp * (-gxe_p / det_gp)
        Jgxe_m = Jm * (-gxe_m / det_gm)
        dJgxe_dξ_arr[p, i, j] = (Jgxe_p - Jgxe_m) / (2 * h_metric)
        # Jg^{ξη} at (ξ, η±h) for ∂/∂η
        Jp2, gxx_p2, gee_p2, gxe_p2 = gnomonic_metric(ξc, ηc + h_metric, R)
        Jm2, gxx_m2, gee_m2, gxe_m2 = gnomonic_metric(ξc, ηc - h_metric, R)
        det_gp2 = gxx_p2 * gee_p2 - gxe_p2^2
        det_gm2 = gxx_m2 * gee_m2 - gxe_m2^2
        Jgxe_p2 = Jp2 * (-gxe_p2 / det_gp2)
        Jgxe_m2 = Jm2 * (-gxe_m2 / det_gm2)
        dJgxe_dη_arr[p, i, j] = (Jgxe_p2 - Jgxe_m2) / (2 * h_metric)
    end

    # Precompute ∂(J·g^{ξξ})/∂ξ and ∂(J·g^{ηη})/∂η for the Laplacian orthogonal correction.
    dJgxx_dξ_arr = zeros(T, 6, Nc, Nc)
    dJgyy_dη_arr = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        ξc = ξ_centers[i]; ηc = η_centers[j]
        # ∂(J·g^{ξξ})/∂ξ
        Jp, gxx_p, gee_p, gxe_p = gnomonic_metric(ξc + h_metric, ηc, R)
        Jm, gxx_m, gee_m, gxe_m = gnomonic_metric(ξc - h_metric, ηc, R)
        det_gp = gxx_p * gee_p - gxe_p^2
        det_gm = gxx_m * gee_m - gxe_m^2
        Jgxx_p = Jp * (gee_p / det_gp)  # J·g^{ξξ} = J·g_{ηη}/det(g)
        Jgxx_m = Jm * (gee_m / det_gm)
        dJgxx_dξ_arr[p, i, j] = (Jgxx_p - Jgxx_m) / (2 * h_metric)
        # ∂(J·g^{ηη})/∂η
        Jp2, gxx_p2, gee_p2, gxe_p2 = gnomonic_metric(ξc, ηc + h_metric, R)
        Jm2, gxx_m2, gee_m2, gxe_m2 = gnomonic_metric(ξc, ηc - h_metric, R)
        det_gp2 = gxx_p2 * gee_p2 - gxe_p2^2
        det_gm2 = gxx_m2 * gee_m2 - gxe_m2^2
        Jgyy_p = Jp2 * (gxx_p2 / det_gp2)  # J·g^{ηη} = J·g_{ξξ}/det(g)
        Jgyy_m = Jm2 * (gxx_m2 / det_gm2)
        dJgyy_dη_arr[p, i, j] = (Jgyy_p - Jgyy_m) / (2 * h_metric)
    end

    # Precompute center-to-center physical distances
    dist_xi = zeros(T, 6, Nc - 1, Nc)
    for p in 1:6, i in 1:(Nc - 1), j in 1:Nc
        v1 = gnomonic_to_cart(ξ_centers[i], η_centers[j], p)
        v2 = gnomonic_to_cart(ξ_centers[i + 1], η_centers[j], p)
        dist_xi[p, i, j] = R * acos(clamp(dot(v1, v2), -1.0, 1.0))
    end
    dist_eta = zeros(T, 6, Nc, Nc - 1)
    for p in 1:6, i in 1:Nc, j in 1:(Nc - 1)
        v1 = gnomonic_to_cart(ξ_centers[i], η_centers[j], p)
        v2 = gnomonic_to_cart(ξ_centers[i], η_centers[j + 1], p)
        dist_eta[p, i, j] = R * acos(clamp(dot(v1, v2), -1.0, 1.0))
    end

    # Precompute cross-panel boundary distances
    dist_xi_bnd = zeros(T, 6, 2, Nc)
    for p in 1:6, j in 1:Nc
        # West boundary: distance from neighbor's cell to this panel's cell 1
        nb_w = PANEL_CONNECTIVITY[p][West]
        i_nb, j_nb = transform_ghost_index(nb_w, 1, j, Nc, Nc, West)
        v_nb = gnomonic_to_cart(ξ_centers[clamp(i_nb, 1, Nc)], η_centers[clamp(j_nb, 1, Nc)], nb_w.neighbor_panel)
        v_loc = gnomonic_to_cart(ξ_centers[1], η_centers[j], p)
        dist_xi_bnd[p, 1, j] = R * acos(clamp(dot(v_nb, v_loc), -1.0, 1.0))
        # East boundary: distance from this panel's cell Nc to neighbor's cell
        nb_e = PANEL_CONNECTIVITY[p][East]
        i_nb, j_nb = transform_ghost_index(nb_e, 1, j, Nc, Nc, East)
        v_loc = gnomonic_to_cart(ξ_centers[Nc], η_centers[j], p)
        v_nb = gnomonic_to_cart(ξ_centers[clamp(i_nb, 1, Nc)], η_centers[clamp(j_nb, 1, Nc)], nb_e.neighbor_panel)
        dist_xi_bnd[p, 2, j] = R * acos(clamp(dot(v_loc, v_nb), -1.0, 1.0))
    end
    dist_eta_bnd = zeros(T, 6, Nc, 2)
    for p in 1:6, i in 1:Nc
        # South boundary: distance from neighbor's cell to this panel's cell 1
        nb_s = PANEL_CONNECTIVITY[p][South]
        i_nb, j_nb = transform_ghost_index(nb_s, 1, i, Nc, Nc, South)
        v_nb = gnomonic_to_cart(ξ_centers[clamp(i_nb, 1, Nc)], η_centers[clamp(j_nb, 1, Nc)], nb_s.neighbor_panel)
        v_loc = gnomonic_to_cart(ξ_centers[i], η_centers[1], p)
        dist_eta_bnd[p, i, 1] = R * acos(clamp(dot(v_nb, v_loc), -1.0, 1.0))
        # North boundary: distance from this panel's cell Nc to neighbor's cell
        nb_n = PANEL_CONNECTIVITY[p][North]
        i_nb, j_nb = transform_ghost_index(nb_n, 1, i, Nc, Nc, North)
        v_loc = gnomonic_to_cart(ξ_centers[i], η_centers[Nc], p)
        v_nb = gnomonic_to_cart(ξ_centers[clamp(i_nb, 1, Nc)], η_centers[clamp(j_nb, 1, Nc)], nb_n.neighbor_panel)
        dist_eta_bnd[p, i, 2] = R * acos(clamp(dot(v_loc, v_nb), -1.0, 1.0))
    end

    grid = CubedSphereGrid{T}(
        Nc, Ng, R, ξ_centers, η_centers, ξ_edges, η_edges,
        lon, lat, area, dx, dy, dξ, dη, rotation_angles, rotation_matrices,
        J_arr, ginv_ξξ, ginv_ηη, ginv_ξη,
        dξ_dlon_arr, dξ_dlat_arr, dη_dlon_arr, dη_dlat_arr,
        d2ξ_dlon2_arr, d2ξ_dlondlat_arr, d2ξ_dlat2_arr,
        d2η_dlon2_arr, d2η_dlondlat_arr, d2η_dlat2_arr,
        dJgxe_dξ_arr, dJgxe_dη_arr, dJgxx_dξ_arr, dJgyy_dη_arr,
        dist_xi, dist_eta, dist_xi_bnd, dist_eta_bnd,
        zeros(T, 6, Nc, Nc, 9), zeros(T, 6, Nc, Nc, 9)
    )

    # Compute FV3 super-grid angular values
    compute_super_grid!(grid.sin_sg, grid.cos_sg, grid)

    return grid
end

total_area(grid::CubedSphereGrid) = sum(grid.area)

# ---------------------------------------------------------------------------
# ESS Grid trait — Tier C + Tier M (curvilinear gnomonic cubed sphere)
#
# Flat indexing convention: cell `(p, i, j)` maps to flat index
# `(p - 1) * Nc^2 + (j - 1) * Nc + i`. Tier-M tensor methods reuse the
# eagerly-precomputed metric arrays already stored on the grid struct.
# ---------------------------------------------------------------------------

n_dims(::CubedSphereGrid) = 2
axis_names(::CubedSphereGrid) = (:xi, :eta)
n_cells(g::CubedSphereGrid) = 6 * g.Nc * g.Nc

function _cs_axis_idx(::CubedSphereGrid, axis::Symbol)
    axis === :xi && return 1
    axis === :eta && return 2
    throw(ArgumentError("cubed_sphere: unknown axis :$axis (expected :xi or :eta)"))
end

@inline _cs_flat(p::Int, i::Int, j::Int, Nc::Int) = (p - 1) * Nc * Nc + (j - 1) * Nc + i

function cell_centers(g::CubedSphereGrid{T}, axis::Symbol) where {T}
    _cs_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_centers, axis)) do
        Nc = g.Nc
        out = Vector{T}(undef, 6 * Nc * Nc)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            out[_cs_flat(p, i, j, Nc)] = axis === :xi ? g.ξ_centers[i] : g.η_centers[j]
        end
        return out
    end
end

function cell_widths(g::CubedSphereGrid{T}, axis::Symbol) where {T}
    _cs_axis_idx(g, axis)
    return _grid_memo!(g, (:cell_widths, axis)) do
        n = 6 * g.Nc * g.Nc
        return fill(axis === :xi ? T(g.dξ) : T(g.dη), n)
    end
end

function cell_volume(g::CubedSphereGrid{T}) where {T}
    return _grid_memo!(g, :cell_volume) do
        Nc = g.Nc
        out = Vector{T}(undef, 6 * Nc * Nc)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            out[_cs_flat(p, i, j, Nc)] = g.area[p, i, j]
        end
        return out
    end
end

function neighbor_indices(g::CubedSphereGrid, axis::Symbol, offset::Int)
    d = _cs_axis_idx(g, axis)
    return _grid_memo!(g, (:neighbor_indices, axis, offset)) do
        Nc = g.Nc
        out = Vector{Int}(undef, 6 * Nc * Nc)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            ii = d == 1 ? i + offset : i
            jj = d == 2 ? j + offset : j
            if 1 <= ii <= Nc && 1 <= jj <= Nc
                out[k] = _cs_flat(p, ii, jj, Nc)
            elseif abs(offset) == 1
                # Single-cell hop crossing a panel edge: resolve via the
                # connectivity table.
                edge = if d == 1
                    offset == 1 ? East : West
                else
                    offset == 1 ? North : South
                end
                nb = PANEL_CONNECTIVITY[p][edge]
                # `transform_ghost_index(nb, depth, along, ni, nj, src_dir)`:
                # depth=1 means "one cell off the panel edge" (the immediate
                # neighbour on the adjacent panel).
                ip, jp = if d == 1
                    transform_ghost_index(nb, 1, j, Nc, Nc, edge)
                else
                    transform_ghost_index(nb, 1, i, Nc, Nc, edge)
                end
                out[k] = (1 <= ip <= Nc && 1 <= jp <= Nc) ?
                    _cs_flat(nb.neighbor_panel, ip, jp, Nc) : 0
            else
                # Multi-cell hops past a panel edge are out of scope for the
                # trait; assemblers stack ±1 hops instead. Sentinel = 0.
                out[k] = 0
            end
        end
        return out
    end
end

function boundary_mask(g::CubedSphereGrid, axis::Symbol, side::Symbol)
    _cs_axis_idx(g, axis)
    side in (:lower, :upper) ||
        throw(ArgumentError("cubed_sphere: side must be :lower or :upper; got :$side"))
    # Every cell on a closed cubed sphere has neighbours via panel
    # connectivity, so there are no "boundary" cells.
    return _grid_memo!(g, (:boundary_mask, axis, side)) do
        falses(6 * g.Nc * g.Nc)
    end
end

# Tier M — pull from the eagerly-materialized metric arrays.

function metric_g(g::CubedSphereGrid{T}) where {T}
    return _grid_memo!(g, :metric_g) do
        Nc = g.Nc
        out = zeros(T, 6 * Nc * Nc, 2, 2)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            ginv_xx = g.ginv_ξξ[p, i, j]
            ginv_yy = g.ginv_ηη[p, i, j]
            ginv_xe = g.ginv_ξη[p, i, j]
            det_inv = ginv_xx * ginv_yy - ginv_xe * ginv_xe
            # g_{ij} = inv(g^{ij}) (pointwise 2×2 inverse).
            out[k, 1, 1] = ginv_yy / det_inv
            out[k, 2, 2] = ginv_xx / det_inv
            out[k, 1, 2] = -ginv_xe / det_inv
            out[k, 2, 1] = -ginv_xe / det_inv
        end
        return out
    end
end

function metric_ginv(g::CubedSphereGrid{T}) where {T}
    return _grid_memo!(g, :metric_ginv) do
        Nc = g.Nc
        out = zeros(T, 6 * Nc * Nc, 2, 2)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            out[k, 1, 1] = g.ginv_ξξ[p, i, j]
            out[k, 2, 2] = g.ginv_ηη[p, i, j]
            out[k, 1, 2] = g.ginv_ξη[p, i, j]
            out[k, 2, 1] = g.ginv_ξη[p, i, j]
        end
        return out
    end
end

function metric_jacobian(g::CubedSphereGrid{T}) where {T}
    return _grid_memo!(g, :metric_jacobian) do
        Nc = g.Nc
        out = Vector{T}(undef, 6 * Nc * Nc)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            out[_cs_flat(p, i, j, Nc)] = g.J[p, i, j]
        end
        return out
    end
end

function metric_dgij_dxk(g::CubedSphereGrid{T}) where {T}
    # Reconstruct ∂g_ij/∂x^k from ∂(J g^{ij})/∂x^k via the product rule:
    # ∂(J g^{ij})/∂x^k = (∂J/∂x^k) g^{ij} + J (∂g^{ij}/∂x^k)
    # Then ∂g_ij/∂x^k = -g_im g_jn (∂g^{mn}/∂x^k).
    # The grid stores `dJgxx_dξ`, `dJgyy_dη`, `dJgxe_dξ`, `dJgxe_dη` which
    # are sufficient for the assemblers' first-derivative correction terms;
    # we compute the full ∂g_ij/∂x^k array by finite-differencing g_ij over
    # the same uniform (ξ, η) grid that the metric arrays live on.
    return _grid_memo!(g, :metric_dgij_dxk) do
        Nc = g.Nc
        N = 6 * Nc * Nc
        out = zeros(T, N, 2, 2, 2)
        # Use precomputed 2D arrays of g_ij from metric_g.
        gij = metric_g(g)
        # Centered differences in computational space (ξ, η). Boundary cells
        # use one-sided differences (no panel-crossing — sufficient for the
        # 9-point stencil's correction terms).
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            # ∂/∂ξ
            ip = i + 1; im = i - 1
            ip_c = clamp(ip, 1, Nc); im_c = clamp(im, 1, Nc)
            inv_dξ = T(1) / (T(ip_c - im_c) * T(g.dξ))
            for a in 1:2, b in 1:2
                kp = _cs_flat(p, ip_c, j, Nc); km = _cs_flat(p, im_c, j, Nc)
                out[k, a, b, 1] = (gij[kp, a, b] - gij[km, a, b]) * inv_dξ
            end
            # ∂/∂η
            jp = j + 1; jm = j - 1
            jp_c = clamp(jp, 1, Nc); jm_c = clamp(jm, 1, Nc)
            inv_dη = T(1) / (T(jp_c - jm_c) * T(g.dη))
            for a in 1:2, b in 1:2
                kp = _cs_flat(p, i, jp_c, Nc); km = _cs_flat(p, i, jm_c, Nc)
                out[k, a, b, 2] = (gij[kp, a, b] - gij[km, a, b]) * inv_dη
            end
        end
        return out
    end
end

function coord_jacobian(g::CubedSphereGrid{T}, target::Symbol) where {T}
    target === :lon_lat ||
        throw(ArgumentError("cubed_sphere: coord_jacobian only supports target=:lon_lat; got :$target"))
    return _grid_memo!(g, (:coord_jacobian, target)) do
        Nc = g.Nc
        out = zeros(T, 6 * Nc * Nc, 2, 2)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            out[k, 1, 1] = g.dξ_dlon[p, i, j]
            out[k, 1, 2] = g.dξ_dlat[p, i, j]
            out[k, 2, 1] = g.dη_dlon[p, i, j]
            out[k, 2, 2] = g.dη_dlat[p, i, j]
        end
        return out
    end
end

function coord_jacobian_second(g::CubedSphereGrid{T}, target::Symbol) where {T}
    target === :lon_lat ||
        throw(ArgumentError("cubed_sphere: coord_jacobian_second only supports target=:lon_lat; got :$target"))
    return _grid_memo!(g, (:coord_jacobian_second, target)) do
        Nc = g.Nc
        out = zeros(T, 6 * Nc * Nc, 2, 2, 2)
        @inbounds for p in 1:6, j in 1:Nc, i in 1:Nc
            k = _cs_flat(p, i, j, Nc)
            out[k, 1, 1, 1] = g.d2ξ_dlon2[p, i, j]
            out[k, 1, 1, 2] = g.d2ξ_dlondlat[p, i, j]
            out[k, 1, 2, 1] = g.d2ξ_dlondlat[p, i, j]
            out[k, 1, 2, 2] = g.d2ξ_dlat2[p, i, j]
            out[k, 2, 1, 1] = g.d2η_dlon2[p, i, j]
            out[k, 2, 1, 2] = g.d2η_dlondlat[p, i, j]
            out[k, 2, 2, 1] = g.d2η_dlondlat[p, i, j]
            out[k, 2, 2, 2] = g.d2η_dlat2[p, i, j]
        end
        return out
    end
end
