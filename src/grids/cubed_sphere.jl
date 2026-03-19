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
    # Metric tensor components at cell centers
    J::Array{T,3}           # Jacobian
    ginv_ξξ::Array{T,3}     # Inverse metric g^{ξξ}
    ginv_ηη::Array{T,3}     # Inverse metric g^{ηη}
    ginv_ξη::Array{T,3}     # Inverse metric g^{ξη}
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
    for p in 1:6, dir in (West, East, South, North)
        nb = PANEL_CONNECTIVITY[p][dir]
        rotation_angles[(p, dir)] = compute_rotation_angle(p, dir, nb.neighbor_panel, nb.neighbor_edge)
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

    CubedSphereGrid{T}(Nc, Ng, R, ξ_centers, η_centers, ξ_edges, η_edges,
        lon, lat, area, dx, dy, dξ, dη, rotation_angles,
        J_arr, ginv_ξξ, ginv_ηη, ginv_ξη)
end

total_area(grid::CubedSphereGrid) = sum(grid.area)
