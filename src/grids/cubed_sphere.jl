"""
Cubed-sphere grid: `CubedSphereGrid(Nc; R=1.0, Ng=3)`.
"""

struct CubedSphereGrid{T} <: AbstractCubedSphereGrid
    Nc::Int; Ng::Int; R::T
    ξ_centers::Vector{T}; η_centers::Vector{T}
    ξ_edges::Vector{T}; η_edges::Vector{T}
    lon::Array{T,3}; lat::Array{T,3}
    area::Array{T,3}; dx::Array{T,3}; dy::Array{T,3}
    dξ::T; dη::T
    rotation_angles::Dict{Tuple{Int,EdgeDirection},Float64}
    # 2×2 rotation matrices for vector field ghost cell filling:
    # (M11, M12, M21, M22) transforming neighbor (u_ξ, u_η) to local basis
    rotation_matrices::Dict{Tuple{Int,EdgeDirection},NTuple{4,Float64}}
    # Metric tensor components at cell centers
    J::Array{T,3}           # Jacobian
    ginv_ξξ::Array{T,3}     # Inverse metric g^{ξξ}
    ginv_ηη::Array{T,3}     # Inverse metric g^{ηη}
    ginv_ξη::Array{T,3}     # Inverse metric g^{ξη}
    # Physical-to-computational coordinate Jacobian: d(ξ,η)/d(lon,lat)
    dξ_dlon::Array{T,3}
    dξ_dlat::Array{T,3}
    dη_dlon::Array{T,3}
    dη_dlat::Array{T,3}
    # Second-derivative coordinate Jacobian: d²(ξ,η)/d(lon,lat)²
    d2ξ_dlon2::Array{T,3}
    d2ξ_dlondlat::Array{T,3}
    d2ξ_dlat2::Array{T,3}
    d2η_dlon2::Array{T,3}
    d2η_dlondlat::Array{T,3}
    d2η_dlat2::Array{T,3}
    # Derivatives of J·g^{ab} for the Laplacian metric corrections
    dJgxe_dξ::Array{T,3}   # ∂(J·g^{ξη})/∂ξ at cell centers
    dJgxe_dη::Array{T,3}   # ∂(J·g^{ξη})/∂η at cell centers
    dJgxx_dξ::Array{T,3}   # ∂(J·g^{ξξ})/∂ξ at cell centers
    dJgyy_dη::Array{T,3}   # ∂(J·g^{ηη})/∂η at cell centers
    # Center-to-center physical distances
    dist_xi::Array{T,3}     # (6, Nc-1, Nc): distance between cell (i,j) and (i+1,j)
    dist_eta::Array{T,3}    # (6, Nc, Nc-1): distance between cell (i,j) and (i,j+1)
    # Boundary distances: cross-panel center-to-center distances at panel edges
    dist_xi_bnd::Array{T,3} # (6, 2, Nc): [p, 1=west/2=east, j] boundary distances in ξ
    dist_eta_bnd::Array{T,3}# (6, Nc, 2): [p, i, 1=south/2=north] boundary distances in η
end

function CubedSphereGrid(Nc::Int; R = 1.0, Ng::Int = 3)
    T = typeof(float(R))
    dξ = T(π / 2) / Nc; dη = dξ
    ξ_edges = [T(-π / 4) + (i - 1) * dξ for i in 1:Nc+1]
    η_edges = [T(-π / 4) + (j - 1) * dη for j in 1:Nc+1]
    ξ_centers = [(ξ_edges[i] + ξ_edges[i+1]) / 2 for i in 1:Nc]
    η_centers = [(η_edges[j] + η_edges[j+1]) / 2 for j in 1:Nc]

    lon = zeros(T, 6, Nc, Nc); lat = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        lo, la = gnomonic_to_lonlat(ξ_centers[i], η_centers[j], p)
        lon[p, i, j] = lo; lat[p, i, j] = la
    end

    area = zeros(T, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        area[p, i, j] = compute_cell_area((ξ_edges[i], ξ_edges[i+1]), (η_edges[j], η_edges[j+1]), R, p)
    end

    dx = zeros(T, 6, Nc + 1, Nc)
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        dx[p, i, j] = compute_edge_length(ξ_edges[i], η_edges[j], ξ_edges[i], η_edges[j+1], R, p)
    end
    dy = zeros(T, 6, Nc, Nc + 1)
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        dy[p, i, j] = compute_edge_length(ξ_edges[i], η_edges[j], ξ_edges[i+1], η_edges[j], R, p)
    end

    rotation_angles = Dict{Tuple{Int,EdgeDirection},Float64}()
    rotation_matrices = Dict{Tuple{Int,EdgeDirection},NTuple{4,Float64}}()
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
    h_metric = T(1e-7)
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
    for p in 1:6, i in 1:Nc-1, j in 1:Nc
        v1 = gnomonic_to_cart(ξ_centers[i], η_centers[j], p)
        v2 = gnomonic_to_cart(ξ_centers[i+1], η_centers[j], p)
        dist_xi[p, i, j] = R * acos(clamp(dot(v1, v2), -1.0, 1.0))
    end
    dist_eta = zeros(T, 6, Nc, Nc - 1)
    for p in 1:6, i in 1:Nc, j in 1:Nc-1
        v1 = gnomonic_to_cart(ξ_centers[i], η_centers[j], p)
        v2 = gnomonic_to_cart(ξ_centers[i], η_centers[j+1], p)
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

    CubedSphereGrid{T}(Nc, Ng, R, ξ_centers, η_centers, ξ_edges, η_edges,
        lon, lat, area, dx, dy, dξ, dη, rotation_angles, rotation_matrices,
        J_arr, ginv_ξξ, ginv_ηη, ginv_ξη,
        dξ_dlon_arr, dξ_dlat_arr, dη_dlon_arr, dη_dlat_arr,
        d2ξ_dlon2_arr, d2ξ_dlondlat_arr, d2ξ_dlat2_arr,
        d2η_dlon2_arr, d2η_dlondlat_arr, d2η_dlat2_arr,
        dJgxe_dξ_arr, dJgxe_dη_arr, dJgxx_dξ_arr, dJgyy_dη_arr,
        dist_xi, dist_eta, dist_xi_bnd, dist_eta_bnd)
end

total_area(grid::CubedSphereGrid) = sum(grid.area)
