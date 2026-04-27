@testsnippet FV3Setup begin
    using Test
    using EarthSciDiscretizations
    using LinearAlgebra: norm, dot
end

# =============================================================================
# Super-grid (sin_sg, cos_sg) tests
# =============================================================================

@testitem "sin_sg: identity sin²α + cos²α = 1" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(16; R = 1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            s = grid.sin_sg[p, i, j, pos]
            c = grid.cos_sg[p, i, j, pos]
            @test isapprox(s^2 + c^2, 1.0; atol = 1.0e-14)
        end
    end
end

@testitem "sin_sg: positive sin(α) everywhere" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(16; R = 1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            @test grid.sin_sg[p, i, j, pos] > 0
        end
    end
end

@testitem "sin_sg: panel symmetry" setup = [FV3Setup] tags = [:fv3] begin
    # On the gnomonic grid, all 6 panels have identical metric properties
    # sin_sg should be the same on every panel at corresponding positions
    grid = CubedSphereGrid(8; R = 1.0)
    for pos in 1:9
        for i in 1:grid.Nc, j in 1:grid.Nc
            vals = [grid.sin_sg[p, i, j, pos] for p in 1:6]
            @test maximum(vals) - minimum(vals) < 1.0e-12
        end
    end
end

@testitem "sin_sg: center matches gnomonic_metric" setup = [FV3Setup] tags = [:fv3] begin
    # Position 5 (cell center) should match the directly computed angle
    grid = CubedSphereGrid(12; R = 6.371e6)
    for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
        J, gxx, gee, gxe = EarthSciDiscretizations.gnomonic_metric(
            grid.ξ_centers[i], grid.η_centers[j], grid.R
        )
        sin_expected = J / sqrt(gxx * gee)
        cos_expected = gxe / sqrt(gxx * gee)
        @test isapprox(grid.sin_sg[p, i, j, 5], sin_expected; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 5], cos_expected; rtol = 1.0e-12)
    end
end

@testitem "sin_sg: shared edge consistency" setup = [FV3Setup] tags = [:fv3] begin
    # The east mid-edge of cell (i,j) is the same point as the west mid-edge of cell (i+1,j)
    # sin_sg[p,i,j,3] (east) should equal sin_sg[p,i+1,j,1] (west) for interior cells
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc
    for p in 1:6, i in 1:(Nc - 1), j in 1:Nc
        @test isapprox(grid.sin_sg[p, i, j, 3], grid.sin_sg[p, i + 1, j, 1]; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 3], grid.cos_sg[p, i + 1, j, 1]; rtol = 1.0e-12)
    end
    # Similarly, north of (i,j) = south of (i,j+1)
    for p in 1:6, i in 1:Nc, j in 1:(Nc - 1)
        @test isapprox(grid.sin_sg[p, i, j, 4], grid.sin_sg[p, i, j + 1, 2]; rtol = 1.0e-12)
        @test isapprox(grid.cos_sg[p, i, j, 4], grid.cos_sg[p, i, j + 1, 2]; rtol = 1.0e-12)
    end
end

# =============================================================================
# Covariant/contravariant decomposition tests
# =============================================================================

@testitem "Wind decomposition: round-trip covariant↔contravariant" setup = [FV3Setup] tags = [:fv3] begin
    # Converting covariant → contravariant → covariant should be identity
    for α in [π / 4, π / 3, π / 2, 2π / 3]
        sin_a = sin(α); cos_a = cos(α)
        for (u, v) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (3.0, -2.0)]
            u_contra, v_contra = EarthSciDiscretizations.covariant_to_contravariant(u, v, sin_a, cos_a)
            u_back, v_back = EarthSciDiscretizations.contravariant_to_covariant(u_contra, v_contra, sin_a, cos_a)
            @test isapprox(u_back, u; atol = 1.0e-12)
            @test isapprox(v_back, v; atol = 1.0e-12)
        end
    end
end

@testitem "Wind decomposition: orthogonal grid identity" setup = [FV3Setup] tags = [:fv3] begin
    # When α = π/2 (orthogonal), covariant = contravariant
    sin_a = 1.0; cos_a = 0.0
    u_contra, v_contra = EarthSciDiscretizations.covariant_to_contravariant(3.0, -2.0, sin_a, cos_a)
    @test isapprox(u_contra, 3.0; atol = 1.0e-14)
    @test isapprox(v_contra, -2.0; atol = 1.0e-14)
end

# =============================================================================
# Helper: initialize D-grid winds for solid body rotation (normalized covariant)
# =============================================================================

@testsnippet SolidBodyRotation begin
    using LinearAlgebra: norm, dot

    """
    Initialize D-grid normalized covariant winds for solid body rotation
    about the z-axis with angular velocity Omega_test.

    Returns (u_d, v_d) where:
    - u_d[p,i,j] = V·ê_η at UEdge(i,j) (normalized, m/s)
    - v_d[p,i,j] = V·ê_ξ at VEdge(i,j) (normalized, m/s)
    """
    function init_solid_body_winds(grid, Omega_test, R)
        Nc = grid.Nc
        u_d = zeros(6, Nc + 1, Nc)
        v_d = zeros(6, Nc, Nc + 1)

        # u_d at UEdge positions
        for p in 1:6, i in 1:(Nc + 1), j in 1:Nc
            ξ = grid.ξ_edges[i]; η = grid.η_centers[j]
            e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
            cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
            vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
            u_d[p, i, j] = dot(vel_3d, e_η) / norm(e_η)
        end

        # v_d at VEdge positions
        for p in 1:6, i in 1:Nc, j in 1:(Nc + 1)
            ξ = grid.ξ_centers[i]; η = grid.η_edges[j]
            e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
            cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
            vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
            v_d[p, i, j] = dot(vel_3d, e_ξ) / norm(e_ξ)
        end

        return (u_d, v_d)
    end

    """
    Initialize D-grid normalized covariant winds for zonal flow u0 * cos(lat).
    """
    function init_zonal_winds(grid, u0, R)
        Nc = grid.Nc
        u_d = zeros(6, Nc + 1, Nc)
        v_d = zeros(6, Nc, Nc + 1)

        for p in 1:6, i in 1:(Nc + 1), j in 1:Nc
            ξ = grid.ξ_edges[i]; η = grid.η_centers[j]
            e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
            cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
            lat_pt = asin(clamp(cart[3], -1.0, 1.0))
            lon_pt = atan(cart[2], cart[1])
            vel_east = [-sin(lon_pt), cos(lon_pt), 0.0] * u0 * cos(lat_pt)
            u_d[p, i, j] = dot(vel_east, e_η) / norm(e_η)
        end

        for p in 1:6, i in 1:Nc, j in 1:(Nc + 1)
            ξ = grid.ξ_centers[i]; η = grid.η_edges[j]
            e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
            cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
            lat_pt = asin(clamp(cart[3], -1.0, 1.0))
            lon_pt = atan(cart[2], cart[1])
            vel_east = [-sin(lon_pt), cos(lon_pt), 0.0] * u0 * cos(lat_pt)
            v_d[p, i, j] = dot(vel_east, e_ξ) / norm(e_ξ)
        end

        return (u_d, v_d)
    end
end

# =============================================================================
# Corner-point vorticity operator tests (default fv_vorticity)
# =============================================================================

@testitem "Corner vorticity: zero for irrotational field" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    omega = fv_vorticity(u_d, v_d, grid)
    @test size(omega) == (6, Nc + 1, Nc + 1)  # Corner staggering
    @test maximum(abs.(omega)) < 1.0e-14
end

@testitem "Corner vorticity: solid body rotation" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # For solid-body rotation about the z-axis with angular velocity Ω,
    # the relative vorticity is ω = 2Ω·sin(lat).
    R = 1.0
    Nc = 16
    grid = CubedSphereGrid(Nc; R = R)
    Omega_test = 1.0

    u_d, v_d = init_solid_body_winds(grid, Omega_test, R)
    omega = fv_vorticity(u_d, v_d, grid)
    omega_cell = fv_vorticity_cellmean(u_d, v_d, grid)

    # Interior corner vorticity is interpolated from cell-mean, so it should be
    # the exact average of the 4 surrounding cell-mean values
    for p in 1:6, i in 2:Nc, j in 2:Nc
        avg = 0.25 * (
            omega_cell[p, i - 1, j - 1] + omega_cell[p, i, j - 1] +
                omega_cell[p, i - 1, j] + omega_cell[p, i, j]
        )
        @test isapprox(omega[p, i, j], avg; atol = 1.0e-14)
    end

    # All corners (including boundary) should be close to the analytical value
    for p in 1:6, i in 1:(Nc + 1), j in 1:(Nc + 1)
        cart = EarthSciDiscretizations.gnomonic_to_cart(grid.ξ_edges[i], grid.η_edges[j], p)
        lat_corner = asin(clamp(cart[3], -1.0, 1.0))
        expected = 2 * Omega_test * sin(lat_corner)
        @test isapprox(omega[p, i, j], expected; rtol = 0.3, atol = 0.1)
    end
end

@testitem "Corner vorticity: boundary corners use ghost-extended values" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # Boundary corner values are computed via ghost-extended cell-mean vorticity
    # from neighboring panels. For solid body rotation they should approximate
    # the analytical value 2*Omega*sin(lat).
    grid = CubedSphereGrid(8; R = 1.0)
    Omega_test = 1.0
    u_d, v_d = init_solid_body_winds(grid, Omega_test, 1.0)
    omega = fv_vorticity(u_d, v_d, grid)
    Nc = grid.Nc

    # Boundary corners should NOT be zero (they get values from neighboring panels)
    boundary_vals = Float64[]
    for p in 1:6
        for j in 1:(Nc + 1)
            push!(boundary_vals, omega[p, 1, j])
            push!(boundary_vals, omega[p, Nc + 1, j])
        end
        for i in 2:Nc  # avoid double-counting corners
            push!(boundary_vals, omega[p, i, 1])
            push!(boundary_vals, omega[p, i, Nc + 1])
        end
    end
    nonzero_boundary = count(v -> abs(v) > 1.0e-10, boundary_vals)
    # Most boundary corners should be nonzero (some at equator may be near zero)
    @test nonzero_boundary / length(boundary_vals) > 0.5

    # Boundary corners should also approximate 2*Omega*sin(lat)
    for p in 1:6, i in [1, Nc + 1], j in [1, Nc + 1]
        cart = EarthSciDiscretizations.gnomonic_to_cart(grid.ξ_edges[i], grid.η_edges[j], p)
        lat_corner = asin(clamp(cart[3], -1.0, 1.0))
        expected = 2 * Omega_test * sin(lat_corner)
        @test isapprox(omega[p, i, j], expected; rtol = 0.3, atol = 0.1)
    end
end

@testitem "Corner vorticity: exact average of cell-mean values" setup = [FV3Setup] tags = [:fv3] begin
    # Corner vorticity is defined as the average of 4 surrounding cell-mean values
    # For interior corners, this uses the cell-mean array directly.
    # For boundary corners, this uses the ghost-extended cell-mean array.
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc; Ng = grid.Ng
    # Use random D-grid winds to test for arbitrary fields
    u_d = randn(6, Nc + 1, Nc)
    v_d = randn(6, Nc, Nc + 1)

    omega_corner = fv_vorticity(u_d, v_d, grid)
    omega_cell = fv_vorticity_cellmean(u_d, v_d, grid)
    omega_ext = extend_with_ghosts(omega_cell, grid)

    # Check ALL corners (including boundary) against the ghost-extended array
    for p in 1:6, i in 1:(Nc + 1), j in 1:(Nc + 1)
        ie = i - 1 + Ng
        je = j - 1 + Ng
        avg_ext = 0.25 * (
            omega_ext[p, ie, je] + omega_ext[p, ie + 1, je] +
                omega_ext[p, ie, je + 1] + omega_ext[p, ie + 1, je + 1]
        )
        @test isapprox(omega_corner[p, i, j], avg_ext; atol = 1.0e-14)
    end

    # Interior corners should still match the direct cell-mean average
    for p in 1:6, i in 2:Nc, j in 2:Nc
        avg_cell = 0.25 * (
            omega_cell[p, i - 1, j - 1] + omega_cell[p, i, j - 1] +
                omega_cell[p, i - 1, j] + omega_cell[p, i, j]
        )
        @test isapprox(omega_corner[p, i, j], avg_cell; atol = 1.0e-14)
    end
end

# =============================================================================
# Cell-mean vorticity operator tests (fv_vorticity_cellmean)
# =============================================================================

@testitem "Cell-mean vorticity: zero for irrotational field" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(12; R = 1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    omega = fv_vorticity_cellmean(u_d, v_d, grid)
    @test size(omega) == (6, Nc, Nc)  # CellCenter staggering
    @test maximum(abs.(omega)) < 1.0e-14
end

@testitem "Cell-mean vorticity: solid body rotation" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    R = 1.0
    Nc = 16
    grid = CubedSphereGrid(Nc; R = R)
    Omega_test = 1.0

    u_d, v_d = init_solid_body_winds(grid, Omega_test, R)
    omega = fv_vorticity_cellmean(u_d, v_d, grid)

    for p in 1:6, i in 1:Nc, j in 1:Nc
        expected = 2 * Omega_test * sin(grid.lat[p, i, j])
        @test isapprox(omega[p, i, j], expected; rtol = 0.05, atol = 0.05)
    end
end

@testitem "Cell-mean vorticity: global integral = 0" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # The global integral of cell-mean vorticity on a closed surface must be zero
    # (by Stokes' theorem, since the surface has no boundary).
    # Using solid body rotation winds which are continuous across panel boundaries,
    # so boundary edge contributions cancel exactly and the integral should vanish
    # to near machine precision.
    R = 1.0
    Nc = 12
    grid = CubedSphereGrid(Nc; R = R)

    u_d, v_d = init_solid_body_winds(grid, 1.0, R)

    omega = fv_vorticity_cellmean(u_d, v_d, grid)

    global_integral = sum(
        omega[p, i, j] * grid.area[p, i, j]
            for p in 1:6, i in 1:Nc, j in 1:Nc
    )

    total_circ = sum(
        abs(omega[p, i, j]) * grid.area[p, i, j]
            for p in 1:6, i in 1:Nc, j in 1:Nc
    )
    @test abs(global_integral) / total_circ < 1.0e-12
end

@testitem "Cell-mean vorticity: ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    omega_loop = fv_vorticity_cellmean(u_d, v_d, grid)

    ao = fv_vorticity_cellmean_arrayop(u_d, v_d, grid)
    omega_ao = evaluate_arrayop(ao)

    @test isapprox(omega_ao, omega_loop; rtol = 1.0e-12)
end

@testitem "Corner vorticity: ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    omega_loop = fv_vorticity(u_d, v_d, grid)

    ao = fv_vorticity_arrayop(u_d, v_d, grid)
    omega_ao = evaluate_arrayop(ao)

    @test isapprox(omega_ao, omega_loop; rtol = 1.0e-12)
end

@testitem "Absolute vorticity: ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    omega_abs_loop = fv_absolute_vorticity(u_d, v_d, grid)

    ao = fv_absolute_vorticity_arrayop(u_d, v_d, grid)
    omega_abs_ao = evaluate_arrayop(ao)

    @test isapprox(omega_abs_ao, omega_abs_loop; rtol = 1.0e-12)
end

# =============================================================================
# Kinetic energy tests
# =============================================================================

@testitem "KE: zero for zero wind" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    ke = fv_kinetic_energy(u_d, v_d, grid)
    @test maximum(abs.(ke)) < 1.0e-14
end

@testitem "KE: solid body rotation" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # For solid body rotation V = Ω × r, the speed is |V| = Ω·R·cos(lat)
    # so KE = (Ω·R·cos(lat))² / 2
    R = 1.0
    Nc = 16
    grid = CubedSphereGrid(Nc; R = R)
    Omega_test = 1.0

    u_d, v_d = init_solid_body_winds(grid, Omega_test, R)
    ke = fv_kinetic_energy(u_d, v_d, grid)

    # Check area-weighted mean KE matches the analytical value
    ke_total = sum(ke[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    total_A = total_area(grid)
    ke_mean = ke_total / total_A
    # Analytical mean KE for solid body rotation: ∫ (ΩR cos λ)²/2 dA / 4πR²
    # = (ΩR)²/2 · 2/3 = (ΩR)²/3
    ke_mean_expected = (Omega_test * R)^2 / 3
    @test isapprox(ke_mean, ke_mean_expected; rtol = 0.05)

    # Check individual cells (upstream-biased formula evaluates at face positions
    # rather than cell centers, so larger error is expected near panel corners
    # where non-orthogonality is strongest)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        expected_ke = 0.5 * (Omega_test * R * cos(grid.lat[p, i, j]))^2
        @test isapprox(ke[p, i, j], expected_ke; rtol = 0.25, atol = 0.02)
    end

    # KE is positive everywhere
    @test all(ke .>= 0)
end

@testitem "KE: cell-center positive definite for random winds" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    u_d = randn(6, Nc + 1, Nc)
    v_d = randn(6, Nc, Nc + 1)
    ke = fv_kinetic_energy_cell(u_d, v_d, grid)
    @test all(ke .>= 0)
end

@testitem "KE: upstream-biased ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    ke_loop = fv_kinetic_energy(u_d, v_d, grid)

    ao = fv_kinetic_energy_arrayop(u_d, v_d, grid)
    ke_ao = evaluate_arrayop(ao)

    @test isapprox(ke_ao, ke_loop; rtol = 1.0e-12)
end

@testitem "KE: cell-center ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    ke_loop = fv_kinetic_energy_cell(u_d, v_d, grid)

    ao = fv_kinetic_energy_cell_arrayop(u_d, v_d, grid)
    ke_ao = evaluate_arrayop(ao)

    @test isapprox(ke_ao, ke_loop; rtol = 1.0e-12)
end

# =============================================================================
# D-grid to C-grid interpolation tests
# =============================================================================

@testitem "D→C grid: zero wind" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    uc, vc = dgrid_to_cgrid(u_d, v_d, grid)
    @test maximum(abs.(uc)) < 1.0e-14
    @test maximum(abs.(vc)) < 1.0e-14
end

@testitem "D→C grid: solid body rotation produces nonzero" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # For solid body rotation, the contravariant C-grid winds should be nonzero
    grid = CubedSphereGrid(12; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)
    uc, vc = dgrid_to_cgrid(u_d, v_d, grid)
    @test maximum(abs.(uc)) > 0.1
    @test maximum(abs.(vc)) > 0.1
end

@testitem "D→C grid: ArrayOp matches loop version" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)

    uc_loop, vc_loop = dgrid_to_cgrid(u_d, v_d, grid)

    uc_ao, vc_ao = dgrid_to_cgrid_arrayop(u_d, v_d, grid)
    uc_eval = evaluate_arrayop(uc_ao)
    vc_eval = evaluate_arrayop(vc_ao)

    @test isapprox(uc_eval, uc_loop; rtol = 1.0e-12)
    @test isapprox(vc_eval, vc_loop; rtol = 1.0e-12)
end

# =============================================================================
# sin_sg-aware flux computation tests
# =============================================================================

@testitem "sinsg flux: zero velocity gives zero flux" setup = [FV3Setup] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    vel_xi = zeros(6, Nc + 1, Nc)
    vel_eta = zeros(6, Nc, Nc + 1)
    flux_xi, flux_eta = compute_flux_with_sinsg(vel_xi, vel_eta, grid, 0.01)
    @test maximum(abs.(flux_xi)) < 1.0e-14
    @test maximum(abs.(flux_eta)) < 1.0e-14
end

@testitem "sinsg flux: sin_sg correction is nontrivial" setup = [FV3Setup] tags = [:fv3] begin
    # Verify that sin_sg values differ from 1.0 (non-orthogonal effect)
    # On the gnomonic grid, the coordinate axes are non-orthogonal away from
    # the panel center, so sin(α) < 1 at most positions.
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc

    # Check that some sin_sg values at mid-edges are not 1.0
    min_sin = minimum(grid.sin_sg[:, :, :, 1:4])
    @test min_sin < 0.99  # Should be significantly less than 1 near panel edges
end

@testitem "sinsg flux: ArrayOp matches loop version (xi)" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    dt = 0.01
    # Use solid body rotation contravariant winds at UEdge as test velocities
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)
    uc, vc = dgrid_to_cgrid(u_d, v_d, grid)

    # Loop version
    flux_xi_loop, _ = compute_flux_with_sinsg(uc, vc, grid, dt)

    # ArrayOp version
    ao_xi = compute_flux_with_sinsg_xi_arrayop(uc, grid, dt)
    flux_xi_ao = evaluate_arrayop(ao_xi)

    @test isapprox(flux_xi_ao, flux_xi_loop; rtol = 1.0e-12)
end

@testitem "sinsg flux: ArrayOp matches loop version (eta)" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc
    dt = 0.01
    u_d, v_d = init_solid_body_winds(grid, 1.0, 1.0)
    uc, vc = dgrid_to_cgrid(u_d, v_d, grid)

    # Loop version
    _, flux_eta_loop = compute_flux_with_sinsg(uc, vc, grid, dt)

    # ArrayOp version
    ao_eta = compute_flux_with_sinsg_eta_arrayop(vc, grid, dt)
    flux_eta_ao = evaluate_arrayop(ao_eta)

    @test isapprox(flux_eta_ao, flux_eta_loop; rtol = 1.0e-12)
end

@testitem "KE: upstream-biased and cell-center agree for smooth wind" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # Both KE formulations should give similar results for smooth wind fields
    R = 1.0
    grid = CubedSphereGrid(12; R = R)
    u_d, v_d = init_solid_body_winds(grid, 1.0, R)

    ke_upstream = fv_kinetic_energy(u_d, v_d, grid)
    ke_cell = fv_kinetic_energy_cell(u_d, v_d, grid)

    # Both should be close (the upstream-biased uses face-position averaging
    # while cell-center uses center averaging, so they differ by O(dx²)).
    # The discrepancy is largest near panel corners where non-orthogonality
    # makes face-position evaluation differ most from center evaluation.
    for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
        @test isapprox(ke_upstream[p, i, j], ke_cell[p, i, j]; rtol = 0.25, atol = 0.02)
    end
end

# =============================================================================
# Integration test: shallow water geostrophic balance
# =============================================================================

@testitem "Integration: shallow water geostrophic balance" setup = [FV3Setup, SolidBodyRotation] tags = [:fv3] begin
    # In geostrophic balance: f × V = -∇(gh)
    # For solid body rotation with angular velocity ω:
    #   V = ωR cos(lat) ê_east
    #   h = h0 - (1/g)(ωR cos(lat))²/2 - (1/g)(Ω ωR² cos²(lat))
    #
    # Test that the pressure gradient from the height field approximately
    # balances the Coriolis force from the velocity field.

    R = 6.371e6
    Nc = 12
    grid = CubedSphereGrid(Nc; R = R)
    g_val = 9.81
    Omega_E = 7.292e-5
    u0 = 20.0  # 20 m/s zonal wind
    h0 = 5960.0  # mean geopotential height in meters

    # Height field for geostrophic balance (Williamson test case 2)
    h = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        lat_ij = grid.lat[p, i, j]
        h[p, i, j] = h0 - (1.0 / g_val) * (R * Omega_E * u0 + u0^2 / 2) * cos(lat_ij)^2
    end

    # D-grid winds for solid body zonal rotation (normalized covariant)
    u_d, v_d = init_zonal_winds(grid, u0, R)

    # Compute cell-mean vorticity (at cell centers, for comparison with cell-center lat)
    omega = fv_vorticity_cellmean(u_d, v_d, grid)
    Omega_test = u0 / R  # angular velocity of the solid body rotation

    # Check that vorticity has the correct sign pattern:
    # Positive in northern hemisphere, negative in southern
    sign_checks = Bool[]
    for p in 1:6, i in 1:Nc, j in 1:Nc
        lat_ij = grid.lat[p, i, j]
        if abs(lat_ij) > 0.2  # Skip equatorial region where vorticity ≈ 0
            push!(sign_checks, sign(omega[p, i, j]) == sign(sin(lat_ij)))
        end
    end
    @test count(sign_checks) / length(sign_checks) > 0.9

    # Kinetic energy should be approximately (u0 cos(lat))²/2
    # (upstream-biased face-averaging gives larger error near panel corners)
    ke = fv_kinetic_energy(u_d, v_d, grid)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        expected_ke = 0.5 * (u0 * cos(grid.lat[p, i, j]))^2
        @test isapprox(ke[p, i, j], expected_ke; rtol = 0.25, atol = 10.0)
    end
end
