@testsnippet GridSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Grid construction" setup = [GridSetup] tags = [:grid] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    @test grid.Nc == Nc
    @test grid.Ng == 3
    @test grid.R == 1.0

    @test length(grid.ξ_centers) == Nc
    @test length(grid.η_centers) == Nc
    @test length(grid.ξ_edges) == Nc + 1
    @test length(grid.η_edges) == Nc + 1

    @test size(grid.lon) == (6, Nc, Nc)
    @test size(grid.lat) == (6, Nc, Nc)
    @test size(grid.area) == (6, Nc, Nc)
    @test size(grid.dx) == (6, Nc + 1, Nc)
    @test size(grid.dy) == (6, Nc, Nc + 1)

    # New arrays: coordinate Jacobian
    @test size(grid.dξ_dlon) == (6, Nc, Nc)
    @test size(grid.dξ_dlat) == (6, Nc, Nc)
    @test size(grid.dη_dlon) == (6, Nc, Nc)
    @test size(grid.dη_dlat) == (6, Nc, Nc)

    # New arrays: center-to-center distances
    @test size(grid.dist_xi) == (6, Nc - 1, Nc)
    @test size(grid.dist_eta) == (6, Nc, Nc - 1)
end

@testitem "Total area equals 4piR^2" setup = [GridSetup] tags = [:grid] begin
    for Nc in [4, 8, 16, 32]
        grid = CubedSphereGrid(Nc; R = 1.0)
        @test isapprox(total_area(grid), 4pi; rtol = 1.0e-10)
    end
end

@testitem "Total area scales with R^2" setup = [GridSetup] tags = [:grid] begin
    R = 6.371e6
    grid = CubedSphereGrid(8; R = R)
    @test isapprox(total_area(grid), 4pi * R^2; rtol = 1.0e-10)
end

@testitem "Grid symmetry" setup = [GridSetup] tags = [:grid] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)

    # All panels should have the same total area
    panel_areas = [sum(grid.area[p, :, :]) for p in 1:6]
    for p in 2:6
        @test isapprox(panel_areas[p], panel_areas[1]; rtol = 1.0e-12)
    end

    # 4-fold symmetry within each panel: area[p,i,j] == area[p,j,i]
    for p in 1:6, i in 1:Nc, j in 1:Nc
        @test isapprox(grid.area[p, i, j], grid.area[p, j, i]; rtol = 1.0e-12)
    end
end

@testitem "Gnomonic projection" setup = [GridSetup] tags = [:grid] begin
    # Panel 1 center should be at lon=0, lat=0
    lon1, lat1 = gnomonic_to_lonlat(0.0, 0.0, 1)
    @test isapprox(lon1, 0.0; atol = 1.0e-14)
    @test isapprox(lat1, 0.0; atol = 1.0e-14)

    # Panel 3 center should be at lat=pi/2 (north pole)
    lon3, lat3 = gnomonic_to_lonlat(0.0, 0.0, 3)
    @test isapprox(lat3, pi / 2; atol = 1.0e-14)

    # Panel 6 center should be at lat=-pi/2 (south pole)
    lon6, lat6 = gnomonic_to_lonlat(0.0, 0.0, 6)
    @test isapprox(lat6, -pi / 2; atol = 1.0e-14)

    # Unit vector check: gnomonic_to_cart should return unit vectors
    using LinearAlgebra: norm
    for p in 1:6
        v = gnomonic_to_cart(0.3, -0.2, p)
        @test isapprox(norm(v), 1.0; atol = 1.0e-14)
    end
end

@testitem "Connectivity consistency" setup = [GridSetup] tags = [:grid] begin
    @test verify_connectivity(8) == true

    # All panels have 4 neighbors
    for p in 1:6
        @test length(PANEL_CONNECTIVITY[p]) == 4
    end

    # No self-adjacency
    for p in 1:6
        for dir in (West, East, South, North)
            nb = PANEL_CONNECTIVITY[p][dir]
            @test nb.neighbor_panel != p
        end
    end
end

@testitem "Connectivity geometric verification" setup = [GridSetup] tags = [:grid] begin
    # Verify that ghost cell mapping produces the correct geometric neighbor
    # by comparing 3D Cartesian positions of boundary cells with their ghost-mapped targets.
    using LinearAlgebra: norm
    Nc = 8
    grid = CubedSphereGrid(Nc)
    Ng = grid.Ng

    # Create fields storing the 3D Cartesian position of each cell center
    q_x = zeros(6, Nc, Nc)
    q_y = zeros(6, Nc, Nc)
    q_z = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        v = gnomonic_to_cart(grid.ξ_centers[i], grid.η_centers[j], p)
        q_x[p, i, j] = v[1]
        q_y[p, i, j] = v[2]
        q_z[p, i, j] = v[3]
    end

    qx_ext = extend_with_ghosts(q_x, grid)
    qy_ext = extend_with_ghosts(q_y, grid)
    qz_ext = extend_with_ghosts(q_z, grid)

    # Helper: get ghost position and expected position for a given boundary cell
    function _check_ghost(p, dir, k, qx_ext, qy_ext, qz_ext, grid, Nc, Ng)
        if dir == West
            ghost_pos = [qx_ext[p, Ng, k + Ng], qy_ext[p, Ng, k + Ng], qz_ext[p, Ng, k + Ng]]
        elseif dir == East
            ghost_pos = [qx_ext[p, Nc + Ng + 1, k + Ng], qy_ext[p, Nc + Ng + 1, k + Ng], qz_ext[p, Nc + Ng + 1, k + Ng]]
        elseif dir == South
            ghost_pos = [qx_ext[p, k + Ng, Ng], qy_ext[p, k + Ng, Ng], qz_ext[p, k + Ng, Ng]]
        else # North
            ghost_pos = [qx_ext[p, k + Ng, Nc + Ng + 1], qy_ext[p, k + Ng, Nc + Ng + 1], qz_ext[p, k + Ng, Nc + Ng + 1]]
        end
        nb = PANEL_CONNECTIVITY[p][dir]
        along_k = (dir == West || dir == East) ? k : k
        i_nb, j_nb = EarthSciDiscretizations.transform_ghost_index(nb, 1, along_k, Nc, Nc, dir)
        i_nb = clamp(i_nb, 1, Nc); j_nb = clamp(j_nb, 1, Nc)
        actual = gnomonic_to_cart(grid.ξ_centers[i_nb], grid.η_centers[j_nb], nb.neighbor_panel)
        return norm(ghost_pos - actual)
    end

    # For each edge ghost cell, verify the ghost-stored 3D position matches
    # the actual position of the source cell on the neighbor panel.
    for p in 1:6, dir in (West, East, South, North), k in 1:Nc
        err = _check_ghost(p, dir, k, qx_ext, qy_ext, qz_ext, grid, Nc, Ng)
        @test err < 1.0e-12
    end
end

@testitem "Positive cell areas" setup = [GridSetup] tags = [:grid] begin
    grid = CubedSphereGrid(8)
    @test all(grid.area .> 0)
end

@testitem "Positive edge lengths" setup = [GridSetup] tags = [:grid] begin
    grid = CubedSphereGrid(8)
    @test all(grid.dx .> 0)
    @test all(grid.dy .> 0)
end

@testitem "Edge length symmetry" setup = [GridSetup] tags = [:grid] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    # dx[p, i, 1] should equal dx[p, i, Nc] by panel symmetry
    for p in 1:6, i in 1:(Nc + 1)
        @test isapprox(grid.dx[p, i, 1], grid.dx[p, i, Nc]; rtol = 1.0e-12)
    end

    # dy[p, 1, j] should equal dy[p, Nc, j] by panel symmetry
    for p in 1:6, j in 1:(Nc + 1)
        @test isapprox(grid.dy[p, 1, j], grid.dy[p, Nc, j]; rtol = 1.0e-12)
    end
end

@testitem "Metric tensor properties" setup = [GridSetup] tags = [:grid] begin
    R = 1.0
    # At the center of a panel (xi=0, eta=0), the metric should be isotropic
    J, g_ξξ, g_ηη, g_ξη = gnomonic_metric(0.0, 0.0, R)
    @test isapprox(g_ξξ, g_ηη; rtol = 1.0e-14)
    @test isapprox(g_ξη, 0.0; atol = 1.0e-14)

    # Jacobian must be positive everywhere
    for xi in range(-pi / 4, pi / 4, length = 10)
        for eta in range(-pi / 4, pi / 4, length = 10)
            J, _, _, _ = gnomonic_metric(xi, eta, R)
            @test J > 0
        end
    end
end

@testitem "Coordinate Jacobian consistency" setup = [GridSetup] tags = [:grid] begin
    # At panel 1 center (ξ=0,η=0), lon≈ξ and lat≈η, so the
    # Jacobian d(ξ,η)/d(lon,lat) should be approximately identity
    jac = compute_coord_jacobian(0.0, 0.0, 1)
    @test isapprox(jac.dξ_dlon, 1.0; rtol = 0.01)
    @test isapprox(jac.dη_dlat, 1.0; rtol = 0.01)
    @test abs(jac.dξ_dlat) < 0.01
    @test abs(jac.dη_dlon) < 0.01

    # Jacobian should be finite everywhere
    grid = CubedSphereGrid(8; R = 1.0)
    for p in 1:6
        @test all(isfinite.(grid.dξ_dlon[p, :, :]))
        @test all(isfinite.(grid.dξ_dlat[p, :, :]))
        @test all(isfinite.(grid.dη_dlon[p, :, :]))
        @test all(isfinite.(grid.dη_dlat[p, :, :]))
    end
end

@testitem "Center-to-center distances are positive" setup = [GridSetup] tags = [:grid] begin
    grid = CubedSphereGrid(8; R = 1.0)
    @test all(grid.dist_xi .> 0)
    @test all(grid.dist_eta .> 0)
end
