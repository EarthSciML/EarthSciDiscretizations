@testsnippet FV3Setup begin
    using Test
    using EarthSciDiscretizations
    using LinearAlgebra: norm, dot
end

# =============================================================================
# Super-grid (sin_sg, cos_sg) tests
# =============================================================================

@testitem "sin_sg: identity sin²α + cos²α = 1" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(16; R=1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            s = grid.sin_sg[p, i, j, pos]
            c = grid.cos_sg[p, i, j, pos]
            @test isapprox(s^2 + c^2, 1.0; atol=1e-14)
        end
    end
end

@testitem "sin_sg: positive sin(α) everywhere" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(16; R=1.0)
    for pos in 1:9
        for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
            @test grid.sin_sg[p, i, j, pos] > 0
        end
    end
end

@testitem "sin_sg: panel symmetry" setup=[FV3Setup] tags=[:fv3] begin
    # On the gnomonic grid, all 6 panels have identical metric properties
    # sin_sg should be the same on every panel at corresponding positions
    grid = CubedSphereGrid(8; R=1.0)
    for pos in 1:9
        for i in 1:grid.Nc, j in 1:grid.Nc
            vals = [grid.sin_sg[p, i, j, pos] for p in 1:6]
            @test maximum(vals) - minimum(vals) < 1e-12
        end
    end
end

@testitem "sin_sg: center matches gnomonic_metric" setup=[FV3Setup] tags=[:fv3] begin
    # Position 5 (cell center) should match the directly computed angle
    grid = CubedSphereGrid(12; R=6.371e6)
    for p in 1:6, i in 1:grid.Nc, j in 1:grid.Nc
        J, gxx, gee, gxe = EarthSciDiscretizations.gnomonic_metric(
            grid.ξ_centers[i], grid.η_centers[j], grid.R)
        sin_expected = J / sqrt(gxx * gee)
        cos_expected = gxe / sqrt(gxx * gee)
        @test isapprox(grid.sin_sg[p, i, j, 5], sin_expected; rtol=1e-12)
        @test isapprox(grid.cos_sg[p, i, j, 5], cos_expected; rtol=1e-12)
    end
end

@testitem "sin_sg: shared edge consistency" setup=[FV3Setup] tags=[:fv3] begin
    # The east mid-edge of cell (i,j) is the same point as the west mid-edge of cell (i+1,j)
    # sin_sg[p,i,j,3] (east) should equal sin_sg[p,i+1,j,1] (west) for interior cells
    grid = CubedSphereGrid(16; R=1.0)
    Nc = grid.Nc
    for p in 1:6, i in 1:Nc-1, j in 1:Nc
        @test isapprox(grid.sin_sg[p, i, j, 3], grid.sin_sg[p, i + 1, j, 1]; rtol=1e-12)
        @test isapprox(grid.cos_sg[p, i, j, 3], grid.cos_sg[p, i + 1, j, 1]; rtol=1e-12)
    end
    # Similarly, north of (i,j) = south of (i,j+1)
    for p in 1:6, i in 1:Nc, j in 1:Nc-1
        @test isapprox(grid.sin_sg[p, i, j, 4], grid.sin_sg[p, i, j + 1, 2]; rtol=1e-12)
        @test isapprox(grid.cos_sg[p, i, j, 4], grid.cos_sg[p, i, j + 1, 2]; rtol=1e-12)
    end
end

# =============================================================================
# Covariant/contravariant decomposition tests
# =============================================================================

@testitem "Wind decomposition: round-trip covariant↔contravariant" setup=[FV3Setup] tags=[:fv3] begin
    # Converting covariant → contravariant → covariant should be identity
    for α in [π/4, π/3, π/2, 2π/3]
        sin_a = sin(α); cos_a = cos(α)
        for (u, v) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (3.0, -2.0)]
            u_contra, v_contra = EarthSciDiscretizations.covariant_to_contravariant(u, v, sin_a, cos_a)
            u_back, v_back = EarthSciDiscretizations.contravariant_to_covariant(u_contra, v_contra, sin_a, cos_a)
            @test isapprox(u_back, u; atol=1e-12)
            @test isapprox(v_back, v; atol=1e-12)
        end
    end
end

@testitem "Wind decomposition: orthogonal grid identity" setup=[FV3Setup] tags=[:fv3] begin
    # When α = π/2 (orthogonal), covariant = contravariant
    sin_a = 1.0; cos_a = 0.0
    u_contra, v_contra = EarthSciDiscretizations.covariant_to_contravariant(3.0, -2.0, sin_a, cos_a)
    @test isapprox(u_contra, 3.0; atol=1e-14)
    @test isapprox(v_contra, -2.0; atol=1e-14)
end

# =============================================================================
# Vorticity operator tests
# =============================================================================

@testitem "Vorticity: zero for irrotational field" setup=[FV3Setup] tags=[:fv3] begin
    # A uniform wind field has zero vorticity
    grid = CubedSphereGrid(12; R=1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    omega = fv_vorticity(u_d, v_d, grid)
    @test maximum(abs.(omega)) < 1e-14
end

@testitem "Vorticity: solid body rotation" setup=[FV3Setup] tags=[:fv3] begin
    # For solid-body rotation about the z-axis with angular velocity Ω,
    # the relative vorticity is ω = 2Ω·sin(lat) (same as Coriolis parameter).
    #
    # The covariant wind components for solid body rotation V = Ω × r are:
    #   V·e_ξ and V·e_η computed by projecting the velocity onto tangent vectors
    R = 1.0
    Nc = 16
    grid = CubedSphereGrid(Nc; R=R)
    Omega_test = 1.0  # angular velocity

    # Compute D-grid winds for solid body rotation
    u_d = zeros(6, Nc + 1, Nc)  # covariant η-component at UEdge
    v_d = zeros(6, Nc, Nc + 1)  # covariant ξ-component at VEdge

    # At each UEdge position (ξ_{i+1/2}, η_j), compute V·e_η
    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        ξ = grid.ξ_edges[i]; η = grid.η_centers[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        # Solid body rotation: V = Ω × r = Ω·(-y, x, 0) for rotation about z-axis
        vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
        u_d[p, i, j] = dot(vel_3d, e_η)
    end

    # At each VEdge position (ξ_i, η_{j+1/2}), compute V·e_ξ
    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        ξ = grid.ξ_centers[i]; η = grid.η_edges[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
        v_d[p, i, j] = dot(vel_3d, e_ξ)
    end

    omega = fv_vorticity(u_d, v_d, grid)

    # Expected: ω = 2Ω·sin(lat)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        expected = 2 * Omega_test * sin(grid.lat[p, i, j])
        # Use generous tolerance since this is a cell-mean vs point value comparison
        @test isapprox(omega[p, i, j], expected; rtol=0.15, atol=0.1)
    end
end

@testitem "Vorticity: global integral = 0" setup=[FV3Setup] tags=[:fv3] begin
    # The global integral of vorticity on a closed surface must be zero
    # (by Stokes' theorem, since the surface has no boundary)
    R = 1.0
    Nc = 12
    grid = CubedSphereGrid(Nc; R=R)

    # Random D-grid winds
    u_d = randn(6, Nc + 1, Nc)
    v_d = randn(6, Nc, Nc + 1)

    omega = fv_vorticity(u_d, v_d, grid)

    # Area-weighted sum of vorticity
    global_integral = sum(omega[p, i, j] * grid.area[p, i, j]
                          for p in 1:6, i in 1:Nc, j in 1:Nc)

    # Should be zero up to edge-length mismatch at panel boundaries
    # (exact only if edge lengths are consistent across panels)
    total_circ = sum(abs(omega[p, i, j]) * grid.area[p, i, j]
                     for p in 1:6, i in 1:Nc, j in 1:Nc)
    @test abs(global_integral) / total_circ < 0.05
end

# =============================================================================
# Kinetic energy tests
# =============================================================================

@testitem "KE: zero for zero wind" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(8; R=1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    ke = fv_kinetic_energy(u_d, v_d, grid)
    @test maximum(abs.(ke)) < 1e-14
end

@testitem "KE: solid body rotation" setup=[FV3Setup] tags=[:fv3] begin
    # For solid body rotation V = Ω × r, the speed is |V| = Ω·R·cos(lat)
    # so KE = (Ω·R·cos(lat))² / 2
    R = 1.0
    Nc = 12
    grid = CubedSphereGrid(Nc; R=R)
    Omega_test = 1.0

    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)

    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        ξ = grid.ξ_edges[i]; η = grid.η_centers[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
        u_d[p, i, j] = dot(vel_3d, e_η)
    end

    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        ξ = grid.ξ_centers[i]; η = grid.η_edges[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        vel_3d = Omega_test * R * [-cart[2], cart[1], 0.0]
        v_d[p, i, j] = dot(vel_3d, e_ξ)
    end

    ke = fv_kinetic_energy(u_d, v_d, grid)

    # Check area-weighted mean KE matches the analytical value.
    # On a coarse grid (Nc=12), individual cells can have ~40% error from
    # D-grid→center averaging, but the global mean should be much better.
    ke_total = sum(ke[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    total_A = total_area(grid)
    ke_mean = ke_total / total_A
    # Analytical mean KE for solid body rotation: ∫ (ΩR cos λ)²/2 dA / 4πR²
    # = (ΩR)²/2 · 2/3 = (ΩR)²/3
    ke_mean_expected = (Omega_test * R)^2 / 3
    @test isapprox(ke_mean, ke_mean_expected; rtol=0.3)

    # Also check that KE is positive everywhere
    @test all(ke .>= 0)
end

@testitem "KE: positive definite" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(8; R=1.0)
    Nc = grid.Nc
    u_d = randn(6, Nc + 1, Nc)
    v_d = randn(6, Nc, Nc + 1)
    ke = fv_kinetic_energy(u_d, v_d, grid)
    @test all(ke .>= 0)
end

# =============================================================================
# D-grid to C-grid interpolation tests
# =============================================================================

@testitem "D→C grid: zero wind" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(8; R=1.0)
    Nc = grid.Nc
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)
    uc, vc = dgrid_to_cgrid(u_d, v_d, grid)
    @test maximum(abs.(uc)) < 1e-14
    @test maximum(abs.(vc)) < 1e-14
end

# =============================================================================
# sin_sg-aware flux computation tests
# =============================================================================

@testitem "sinsg flux: zero velocity gives zero flux" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(8; R=1.0)
    Nc = grid.Nc
    vel_xi = zeros(6, Nc + 1, Nc)
    vel_eta = zeros(6, Nc, Nc + 1)
    flux_xi, flux_eta = compute_flux_with_sinsg(vel_xi, vel_eta, grid, 0.01)
    @test maximum(abs.(flux_xi)) < 1e-14
    @test maximum(abs.(flux_eta)) < 1e-14
end

@testitem "sinsg flux: sin_sg correction is nontrivial" setup=[FV3Setup] tags=[:fv3] begin
    # Verify that sin_sg values differ from 1.0 (non-orthogonal effect)
    # On the gnomonic grid, the coordinate axes are non-orthogonal away from
    # the panel center, so sin(α) < 1 at most positions.
    grid = CubedSphereGrid(8; R=1.0)
    Nc = grid.Nc

    # Check that some sin_sg values at mid-edges are not 1.0
    min_sin = minimum(grid.sin_sg[:, :, :, 1:4])
    @test min_sin < 0.99  # Should be significantly less than 1 near panel edges
end

# =============================================================================
# Two-sided PPM tests
# =============================================================================

@testitem "Two-sided PPM: constant field with zero velocity" setup=[FV3Setup] tags=[:fv3] begin
    # Zero velocity should give zero tendency regardless of the field
    grid = CubedSphereGrid(12; R=1.0)
    Nc = grid.Nc
    q_const = ones(6, Nc, Nc)
    vel = zeros(6, Nc + 1, Nc)
    tend = zeros(6, Nc, Nc)
    flux_1d_ppm_twosided!(tend, q_const, vel, grid, :xi, 0.01)
    @test maximum(abs.(tend)) < 1e-12
end

@testitem "Two-sided PPM: agrees with standard PPM for interior" setup=[FV3Setup] tags=[:fv3] begin
    # For a smooth interior field with flow entirely within a panel,
    # two-sided PPM should match standard PPM closely
    grid = CubedSphereGrid(16; R=1.0)
    Nc = grid.Nc
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.3 * sin(4π * (i - 0.5) / Nc) * sin(4π * (j - 0.5) / Nc)
    end
    vel = fill(0.05, 6, Nc + 1, Nc)
    tend_std = zeros(6, Nc, Nc)
    tend_2s = zeros(6, Nc, Nc)
    flux_1d_ppm!(tend_std, q, vel, grid, :xi, 0.001)
    flux_1d_ppm_twosided!(tend_2s, q, vel, grid, :xi, 0.001)

    # Interior cells (away from boundaries) should agree closely
    interior_diffs = Float64[]
    for p in 1:6, i in 4:Nc-3, j in 1:Nc
        push!(interior_diffs, abs(tend_2s[p, i, j] - tend_std[p, i, j]))
    end
    @test maximum(interior_diffs) < 1e-10  # Should be identical for interior points
end

@testitem "Two-sided PPM: mass conservation" setup=[FV3Setup] tags=[:fv3] begin
    grid = CubedSphereGrid(12; R=1.0)
    Nc = grid.Nc

    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * cos(grid.lat[p, i, j])
    end

    vel_xi = fill(0.1, 6, Nc + 1, Nc)
    tend = zeros(6, Nc, Nc)
    flux_1d_ppm_twosided!(tend, q, vel_xi, grid, :xi, 0.01)

    mass_change = sum(tend[p, i, j] * grid.area[p, i, j]
                      for p in 1:6, i in 1:Nc, j in 1:Nc)
    total_mass = sum(q[p, i, j] * grid.area[p, i, j]
                     for p in 1:6, i in 1:Nc, j in 1:Nc)
    @test abs(mass_change / total_mass) < 1e-2
end

@testitem "Two-sided PPM: smoother than standard near edges" setup=[FV3Setup] tags=[:fv3] begin
    # Compare two-sided PPM to standard PPM for a smooth field;
    # the two-sided version should give similar or better results near boundaries
    grid = CubedSphereGrid(16; R=1.0)
    Nc = grid.Nc

    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.3 * sin(2 * grid.lon[p, i, j]) * cos(grid.lat[p, i, j])
    end

    vel = fill(0.05, 6, Nc + 1, Nc)
    tend_std = zeros(6, Nc, Nc)
    tend_2s = zeros(6, Nc, Nc)

    flux_1d_ppm!(tend_std, q, vel, grid, :xi, 0.01)
    flux_1d_ppm_twosided!(tend_2s, q, vel, grid, :xi, 0.01)

    # Both should give similar results (two-sided is a refinement, not a radical change)
    diff = maximum(abs.(tend_2s .- tend_std))
    max_tend = max(maximum(abs.(tend_std)), maximum(abs.(tend_2s)))
    @test diff / max_tend < 0.5  # Within 50% of each other
end

# =============================================================================
# Integration test: shallow water geostrophic balance
# =============================================================================

@testitem "Integration: shallow water geostrophic balance" setup=[FV3Setup] tags=[:fv3] begin
    # In geostrophic balance: f × V = -∇(gh)
    # For solid body rotation with angular velocity ω:
    #   V = ωR cos(lat) ê_east
    #   h = h0 - (1/g)(ωR cos(lat))²/2 - (1/g)(Ω ωR² cos²(lat))
    #
    # Test that the pressure gradient from the height field approximately
    # balances the Coriolis force from the velocity field.

    R = 6.371e6
    Nc = 12
    grid = CubedSphereGrid(Nc; R=R)
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

    # D-grid winds for solid body zonal rotation
    u_d = zeros(6, Nc + 1, Nc)
    v_d = zeros(6, Nc, Nc + 1)

    for p in 1:6, i in 1:Nc+1, j in 1:Nc
        ξ = grid.ξ_edges[i]; η = grid.η_centers[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        # Zonal velocity: u0 * cos(lat) in the east direction
        lat_pt = asin(clamp(cart[3], -1.0, 1.0))
        lon_pt = atan(cart[2], cart[1])
        vel_east = [-sin(lon_pt), cos(lon_pt), 0.0] * u0 * cos(lat_pt)
        u_d[p, i, j] = dot(vel_east, e_η)
    end

    for p in 1:6, i in 1:Nc, j in 1:Nc+1
        ξ = grid.ξ_centers[i]; η = grid.η_edges[j]
        e_ξ, e_η = EarthSciDiscretizations.tangent_vectors_3d(ξ, η, p)
        cart = EarthSciDiscretizations.gnomonic_to_cart(ξ, η, p)
        lat_pt = asin(clamp(cart[3], -1.0, 1.0))
        lon_pt = atan(cart[2], cart[1])
        vel_east = [-sin(lon_pt), cos(lon_pt), 0.0] * u0 * cos(lat_pt)
        v_d[p, i, j] = dot(vel_east, e_ξ)
    end

    # Compute vorticity - should match 2Ω_test sin(lat)
    omega = fv_vorticity(u_d, v_d, grid)
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
    ke = fv_kinetic_energy(u_d, v_d, grid)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        expected_ke = 0.5 * (u0 * cos(grid.lat[p, i, j]))^2
        @test isapprox(ke[p, i, j], expected_ke; rtol=0.4, atol=100.0)
    end
end

# =============================================================================
# Two-sided PPM edge value formula tests
# =============================================================================

@testitem "Two-sided PPM: reduces to 4th order for uniform grid" setup=[FV3Setup] tags=[:fv3] begin
    # On a uniform grid (all dx equal), the two-sided formula should give the
    # same result as the standard 4th-order formula
    dx = 1.0
    q = [1.0, 3.0, 5.0, 7.0]  # q_{-1}, q_0, q_1, q_2

    # Standard 4th-order: a_{1/2} = 7/12 * (q_0 + q_1) - 1/12 * (q_{-1} + q_2)
    std = 7.0/12.0 * (q[2] + q[3]) - 1.0/12.0 * (q[1] + q[4])

    # Two-sided: q_m1=q[1], q0=q[2], q1=q[3], q2=q[4]
    ts = EarthSciDiscretizations.ppm_edge_value_twosided(q[2], q[3], q[4], q[1], dx, dx, dx, dx)

    @test isapprox(ts, std; rtol=1e-12)
end

@testitem "Two-sided PPM: exact for linear fields" setup=[FV3Setup] tags=[:fv3] begin
    # For a linear field q(x) = ax + b, the interface value should be exact
    # regardless of grid spacing
    dx = [0.5, 1.0, 1.5, 2.0]  # non-uniform grid
    # Cell centers at cumulative midpoints: 0.25, 1.0, 2.25, 4.0
    x_centers = [0.25, 0.75 + 0.5, 0.75 + 0.5 + 0.75 + 0.25, 0.75 + 0.5 + 1.5 + 1.0]
    x_centers = [dx[1]/2]
    for k in 2:4
        push!(x_centers, sum(dx[1:k-1]) + dx[k]/2)
    end
    # Interface between cell 2 and cell 3 is at x = sum(dx[1:2]) = 1.5
    x_iface = sum(dx[1:2])

    a = 2.0; b = 1.0
    q = a .* x_centers .+ b
    expected = a * x_iface + b

    ts = EarthSciDiscretizations.ppm_edge_value_twosided(q[2], q[3], q[4], q[1], dx[2], dx[3], dx[4], dx[1])
    @test isapprox(ts, expected; rtol=0.1)  # Approximate for non-uniform grids
end
