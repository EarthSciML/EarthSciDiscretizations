"""
Precomputed finite-volume stencil weights for efficient ODE function evaluation.

Instead of generating O(6·Nc²) symbolic equations, this module precomputes all
FV stencil weights as numerical arrays at grid construction time, then applies
them directly in the ODE function. This eliminates symbolic overhead entirely.
"""

"""
    FVLaplacianStencil

Precomputed weights and neighbor indices for the full covariant Laplacian
on the cubed sphere, including the cross-metric correction terms.

The stencil uses 9 points: center + 4 cardinal + 4 diagonal.
At each cell (p,i,j), the Laplacian is:
    ∇²φ = Σ_k w[p,i,j,k] * φ[nb_p[p,i,j,k], nb_i[p,i,j,k], nb_j[p,i,j,k]]
"""
struct FVLaplacianStencil
    # Stencil weights: (6, Nc, Nc, 9)
    weights::Array{Float64,4}
    # Neighbor panel indices: (6, Nc, Nc, 9) — Int arrays
    nb_p::Array{Int,4}
    nb_i::Array{Int,4}
    nb_j::Array{Int,4}
end

# Stencil point ordering:
# 1 = center, 2 = east(i+1,j), 3 = west(i-1,j),
# 4 = north(i,j+1), 5 = south(i,j-1),
# 6 = NE(i+1,j+1), 7 = NW(i-1,j+1), 8 = SE(i+1,j-1), 9 = SW(i-1,j-1)

"""
    precompute_laplacian_stencil(grid)

Precompute the FV Laplacian stencil weights and neighbor indices for all
cells on the cubed sphere, including panel boundary cells.

The Laplacian uses:
- Orthogonal part: (1/A) Σ_faces [(φ_nb - φ_center)/dist * edge_length]
- Cross-metric correction: 2g^{ξη}·∂²φ/(∂ξ∂η) + gradient-of-metric terms
"""
function precompute_laplacian_stencil(grid::CubedSphereGrid)
    Nc = grid.Nc
    dξ = grid.dξ; dη = grid.dη
    nstencil = 9

    weights = zeros(6, Nc, Nc, nstencil)
    nb_p = zeros(Int, 6, Nc, Nc, nstencil)
    nb_i = zeros(Int, 6, Nc, Nc, nstencil)
    nb_j = zeros(Int, 6, Nc, Nc, nstencil)

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Resolve all 9 neighbor indices (handles cross-panel)
        center = (p, i, j)
        east   = _neighbor_index(grid, p, i + 1, j)
        west   = _neighbor_index(grid, p, i - 1, j)
        north  = _neighbor_index(grid, p, i, j + 1)
        south  = _neighbor_index(grid, p, i, j - 1)
        # Diagonal neighbors (rotation-aware)
        ne = _diagonal_neighbor(grid, p, i, j, +1, 0, +1, east)
        nw = _diagonal_neighbor(grid, p, i, j, -1, 0, +1, west)
        se = _diagonal_neighbor(grid, p, i, j, +1, 0, -1, east)
        sw = _diagonal_neighbor(grid, p, i, j, -1, 0, -1, west)

        nbs = [center, east, west, north, south, ne, nw, se, sw]
        for k in 1:nstencil
            nb_p[p, i, j, k] = nbs[k][1]
            nb_i[p, i, j, k] = nbs[k][2]
            nb_j[p, i, j, k] = nbs[k][3]
        end

        # --- Compute weights ---
        A = grid.area[p, i, j]
        J_c = grid.J[p, i, j]
        gxe = grid.ginv_ξη[p, i, j]
        dJgxe_dξ = grid.dJgxe_dξ[p, i, j]
        dJgxe_dη = grid.dJgxe_dη[p, i, j]

        # Orthogonal part: flux through each of 4 faces
        # East face (between cell center and east neighbor)
        _dist_east = _physical_distance(grid, center, east)
        _dx_east = _edge_length_xi(grid, p, i + 1, j)
        w_east_orth = _dx_east / (_dist_east * A)

        # West face
        _dist_west = _physical_distance(grid, west, center)
        _dx_west = _edge_length_xi(grid, p, i, j)
        w_west_orth = _dx_west / (_dist_west * A)

        # North face
        _dist_north = _physical_distance(grid, center, north)
        _dy_north = _edge_length_eta(grid, p, i, j + 1)
        w_north_orth = _dy_north / (_dist_north * A)

        # South face
        _dist_south = _physical_distance(grid, south, center)
        _dy_south = _edge_length_eta(grid, p, i, j)
        w_south_orth = _dy_south / (_dist_south * A)

        # Cross-metric correction weights
        cross_d2 = 2 * gxe / (4 * dξ * dη)   # weight for ∂²φ/(∂ξ∂η) cross term
        cross_dxi = (1 / J_c) * dJgxe_dη / (2 * dξ)   # weight for ∂φ/∂ξ correction
        cross_deta = (1 / J_c) * dJgxe_dξ / (2 * dη)  # weight for ∂φ/∂η correction

        # Assemble 9-point stencil weights
        # 1: center
        weights[p, i, j, 1] = -(w_east_orth + w_west_orth + w_north_orth + w_south_orth)
        # 2: east (i+1, j)
        weights[p, i, j, 2] = w_east_orth + cross_dxi
        # 3: west (i-1, j)
        weights[p, i, j, 3] = w_west_orth - cross_dxi
        # 4: north (i, j+1)
        weights[p, i, j, 4] = w_north_orth + cross_deta
        # 5: south (i, j-1)
        weights[p, i, j, 5] = w_south_orth - cross_deta
        # 6: NE (i+1, j+1)
        weights[p, i, j, 6] = +cross_d2
        # 7: NW (i-1, j+1)
        weights[p, i, j, 7] = -cross_d2
        # 8: SE (i+1, j-1)
        weights[p, i, j, 8] = -cross_d2
        # 9: SW (i-1, j-1)
        weights[p, i, j, 9] = +cross_d2
    end

    return FVLaplacianStencil(weights, nb_p, nb_i, nb_j)
end

# Helper: physical (great-circle) distance between two cells
function _physical_distance(grid, pij1, pij2)
    v1 = gnomonic_to_cart(grid.ξ_centers[pij1[2]], grid.η_centers[pij1[3]], pij1[1])
    v2 = gnomonic_to_cart(grid.ξ_centers[pij2[2]], grid.η_centers[pij2[3]], pij2[1])
    return grid.R * acos(clamp(dot(v1, v2), -1.0, 1.0))
end

# Helper: edge length at ξ-edge i (length of cell face perpendicular to ξ)
function _edge_length_xi(grid, p, i, j)
    Nc = grid.Nc
    i_clamped = clamp(i, 1, Nc + 1)
    j_clamped = clamp(j, 1, Nc)
    return grid.dx[p, i_clamped, j_clamped]
end

# Helper: edge length at η-edge j (length of cell face perpendicular to η)
function _edge_length_eta(grid, p, i, j)
    Nc = grid.Nc
    i_clamped = clamp(i, 1, Nc)
    j_clamped = clamp(j, 1, Nc + 1)
    return grid.dy[p, i_clamped, j_clamped]
end

"""
    FVGradientStencil

Precomputed weights for the gradient operator in (lon, lat) or (ξ, η)
coordinates, using the chain-rule transformation.
"""
struct FVGradientStencil
    # Weights for ∂φ/∂x using centered differences: (6, Nc, Nc, 5)
    # Points: center=1, east=2, west=3, north=4, south=5
    weights_lon::Array{Float64,4}
    weights_lat::Array{Float64,4}
    # Same neighbor indices as Laplacian (first 5 only)
    nb_p::Array{Int,4}
    nb_i::Array{Int,4}
    nb_j::Array{Int,4}
end

function precompute_gradient_stencil(grid::CubedSphereGrid)
    Nc = grid.Nc; dξ = grid.dξ; dη = grid.dη
    nstencil = 5

    wlon = zeros(6, Nc, Nc, nstencil)
    wlat = zeros(6, Nc, Nc, nstencil)
    nb_p = zeros(Int, 6, Nc, Nc, nstencil)
    nb_i = zeros(Int, 6, Nc, Nc, nstencil)
    nb_j = zeros(Int, 6, Nc, Nc, nstencil)

    for p in 1:6, i in 1:Nc, j in 1:Nc
        center = (p, i, j)
        east   = _neighbor_index(grid, p, i + 1, j)
        west   = _neighbor_index(grid, p, i - 1, j)
        north  = _neighbor_index(grid, p, i, j + 1)
        south  = _neighbor_index(grid, p, i, j - 1)

        nbs = [center, east, west, north, south]
        for k in 1:nstencil
            nb_p[p, i, j, k] = nbs[k][1]
            nb_i[p, i, j, k] = nbs[k][2]
            nb_j[p, i, j, k] = nbs[k][3]
        end

        # Chain-rule coefficients for lon and lat
        dξ_dlon = grid.dξ_dlon[p, i, j]
        dη_dlon = grid.dη_dlon[p, i, j]
        dξ_dlat = grid.dξ_dlat[p, i, j]
        dη_dlat = grid.dη_dlat[p, i, j]

        # ∂φ/∂lon = (dξ/dlon)·∂φ/∂ξ + (dη/dlon)·∂φ/∂η
        # ∂φ/∂ξ ≈ (φ_east - φ_west) / (2dξ)
        # ∂φ/∂η ≈ (φ_north - φ_south) / (2dη)
        wlon[p, i, j, 1] = 0.0  # center
        wlon[p, i, j, 2] = dξ_dlon / (2 * dξ)   # east: +∂/∂ξ
        wlon[p, i, j, 3] = -dξ_dlon / (2 * dξ)  # west: -∂/∂ξ
        wlon[p, i, j, 4] = dη_dlon / (2 * dη)    # north: +∂/∂η
        wlon[p, i, j, 5] = -dη_dlon / (2 * dη)   # south: -∂/∂η

        wlat[p, i, j, 1] = 0.0
        wlat[p, i, j, 2] = dξ_dlat / (2 * dξ)
        wlat[p, i, j, 3] = -dξ_dlat / (2 * dξ)
        wlat[p, i, j, 4] = dη_dlat / (2 * dη)
        wlat[p, i, j, 5] = -dη_dlat / (2 * dη)
    end

    return FVGradientStencil(wlon, wlat, nb_p, nb_i, nb_j)
end

"""
    apply_laplacian!(du, u, stencil)

Apply the precomputed FV Laplacian stencil to the field `u`, storing
the result in `du`. No ghost cells needed — the stencil uses precomputed
cross-panel neighbor indices.
"""
function apply_laplacian!(du, u, stencil::FVLaplacianStencil)
    @inbounds for p in axes(du, 1), i in axes(du, 2), j in axes(du, 3)
        val = 0.0
        for k in 1:9
            val += stencil.weights[p, i, j, k] *
                   u[stencil.nb_p[p, i, j, k], stencil.nb_i[p, i, j, k], stencil.nb_j[p, i, j, k]]
        end
        du[p, i, j] = val
    end
    return du
end

"""
    apply_gradient!(du_lon, du_lat, u, stencil)

Apply precomputed gradient stencil to compute ∂u/∂lon and ∂u/∂lat.
"""
function apply_gradient!(du_lon, du_lat, u, stencil::FVGradientStencil)
    @inbounds for p in axes(du_lon, 1), i in axes(du_lon, 2), j in axes(du_lon, 3)
        val_lon = 0.0
        val_lat = 0.0
        for k in 1:5
            uk = u[stencil.nb_p[p, i, j, k], stencil.nb_i[p, i, j, k], stencil.nb_j[p, i, j, k]]
            val_lon += stencil.weights_lon[p, i, j, k] * uk
            val_lat += stencil.weights_lat[p, i, j, k] * uk
        end
        du_lon[p, i, j] = val_lon
        du_lat[p, i, j] = val_lat
    end
    return (du_lon, du_lat)
end
