@testsnippet FixesSetup begin
    using Test
    using EarthSciDiscretizations
    using LinearAlgebra: norm, dot
end

# ============================================================
# Task 1: Analytical Jacobian
# ============================================================

@testitem "Analytical Jacobian matches finite differences" setup=[FixesSetup] tags=[:fixes] begin
    # The analytical forward Jacobian should agree with finite-difference
    # computation to high precision (within FD truncation error)
    for panel in 1:6
        for (ξ, η) in [(0.0, 0.0), (0.3, -0.2), (-0.5, 0.4), (0.7, 0.7)]
            fwd = compute_forward_jacobian(ξ, η, panel)

            # Finite-difference reference
            h = 1e-7
            lon_p, lat_p = gnomonic_to_lonlat(ξ + h, η, panel)
            lon_m, lat_m = gnomonic_to_lonlat(ξ - h, η, panel)
            dlon_dξ_fd = (lon_p - lon_m) / (2h)
            dlat_dξ_fd = (lat_p - lat_m) / (2h)

            lon_p, lat_p = gnomonic_to_lonlat(ξ, η + h, panel)
            lon_m, lat_m = gnomonic_to_lonlat(ξ, η - h, panel)
            dlon_dη_fd = (lon_p - lon_m) / (2h)
            dlat_dη_fd = (lat_p - lat_m) / (2h)

            @test isapprox(fwd.dlon_dξ, dlon_dξ_fd, atol=1e-5)
            @test isapprox(fwd.dlon_dη, dlon_dη_fd, atol=1e-5)
            @test isapprox(fwd.dlat_dξ, dlat_dξ_fd, atol=1e-5)
            @test isapprox(fwd.dlat_dη, dlat_dη_fd, atol=1e-5)
        end
    end
end

@testitem "Analytical Jacobian at panel center" setup=[FixesSetup] tags=[:fixes] begin
    # Panel 1 center (0,0): lon = 0, lat = 0
    # dlon/dξ = 1, dlon/dη = 0, dlat/dξ = 0, dlat/dη = 1
    fwd = compute_forward_jacobian(0.0, 0.0, 1)
    @test isapprox(fwd.dlon_dξ, 1.0, atol=1e-12)
    @test isapprox(fwd.dlon_dη, 0.0, atol=1e-12)
    @test isapprox(fwd.dlat_dξ, 0.0, atol=1e-12)
    @test isapprox(fwd.dlat_dη, 1.0, atol=1e-12)
end

@testitem "Coordinate Jacobian inverse is consistent" setup=[FixesSetup] tags=[:fixes] begin
    # Verify J_fwd * J_inv = I (away from poles)
    for panel in [1, 2, 4, 5]  # skip pole panels 3, 6
        for (ξ, η) in [(0.0, 0.0), (0.3, -0.2)]
            fwd = compute_forward_jacobian(ξ, η, panel)
            inv = compute_coord_jacobian(ξ, η, panel)

            # [dlon/dξ  dlon/dη] [dξ/dlon  dξ/dlat]   [1 0]
            # [dlat/dξ  dlat/dη] [dη/dlon  dη/dlat] = [0 1]
            I11 = fwd.dlon_dξ * inv.dξ_dlon + fwd.dlon_dη * inv.dη_dlon
            I12 = fwd.dlon_dξ * inv.dξ_dlat + fwd.dlon_dη * inv.dη_dlat
            I21 = fwd.dlat_dξ * inv.dξ_dlon + fwd.dlat_dη * inv.dη_dlon
            I22 = fwd.dlat_dξ * inv.dξ_dlat + fwd.dlat_dη * inv.dη_dlat

            @test isapprox(I11, 1.0, atol=1e-8)
            @test isapprox(I12, 0.0, atol=1e-8)
            @test isapprox(I21, 0.0, atol=1e-8)
            @test isapprox(I22, 1.0, atol=1e-8)
        end
    end
end

# ============================================================
# Task 3: Vector field rotation at panel boundaries
# ============================================================

@testitem "Rotation matrices are orthogonal" setup=[FixesSetup] tags=[:fixes] begin
    grid = CubedSphereGrid(8)
    for p in 1:6, dir in (West, East, South, North)
        M11, M12, M21, M22 = grid.rotation_matrices[(p, dir)]
        # Orthogonality: M * M^T = I
        @test isapprox(M11^2 + M12^2, 1.0, atol=1e-10)
        @test isapprox(M21^2 + M22^2, 1.0, atol=1e-10)
        @test isapprox(M11*M21 + M12*M22, 0.0, atol=1e-10)
        # det = ±1
        det = M11*M22 - M12*M21
        @test isapprox(abs(det), 1.0, atol=1e-10)
    end
end

@testitem "Vector ghost cells preserve uniform rotation field" setup=[FixesSetup] tags=[:fixes] begin
    # A solid body rotation velocity field should be smooth across panel boundaries.
    # Test: create a uniform angular velocity field and check ghost cell continuity.
    Nc = 8
    grid = CubedSphereGrid(Nc)

    # Create a vector field from solid body rotation (Ω × r in Cartesian → project to ξ,η)
    uξ = zeros(6, Nc, Nc)
    uη = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        ξ = grid.ξ_centers[i]; η = grid.η_centers[j]
        e_ξ, e_η = tangent_vectors_3d(ξ, η, p)
        # Solid body rotation about z-axis: v = Ω × r = (-y, x, 0)
        cart = gnomonic_to_cart(ξ, η, p)
        v_cart = [-cart[2], cart[1], 0.0]
        # Project onto unit tangent vectors
        nξ = norm(e_ξ); nη = norm(e_η)
        uξ[p, i, j] = dot(v_cart, e_ξ) / nξ^2  # contravariant component
        uη[p, i, j] = dot(v_cart, e_η) / nη^2
    end

    uξ_ext, uη_ext = extend_with_ghosts_vector(uξ, uη, grid)
    Ng = grid.Ng

    # Check that ghost cells have reasonable values (not NaN, not wildly different)
    @test all(isfinite, uξ_ext)
    @test all(isfinite, uη_ext)

    # Interior values should be preserved
    for p in 1:6, i in 1:Nc, j in 1:Nc
        @test uξ_ext[p, i+Ng, j+Ng] == uξ[p, i, j]
        @test uη_ext[p, i+Ng, j+Ng] == uη[p, i, j]
    end
end

# ============================================================
# Task 2: Laplacian cross-metric correction
# ============================================================

@testitem "Laplacian of constant is zero (with cross-metric fix)" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    phi = ones(6, Nc, Nc) * 42.0
    lapl = fv_laplacian(phi, grid)
    result = evaluate_arrayop(lapl)
    @test all(abs.(result) .< 1e-10)
end

@testitem "Laplacian convergence improves with resolution" setup=[FixesSetup] tags=[:fixes] begin
    errors = Float64[]
    for Nc in [8, 16]
        grid = CubedSphereGrid(Nc)
        # Test function: cos(2ξ) * cos(2η)
        phi = zeros(6, Nc, Nc)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
        end
        lapl = fv_laplacian(phi, grid)
        result = evaluate_arrayop(lapl)
        push!(errors, maximum(abs.(result)))
    end
    # Error should decrease with resolution
    @test errors[2] < errors[1]
end

# ============================================================
# Task 4: Diagonal neighbor resolution
# ============================================================

@testitem "Diagonal neighbors at panel boundaries are on correct physical location" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    # Test at the east boundary of panel 1 (i=Nc)
    p = 1; i = Nc; j = Nc ÷ 2

    pip1 = EarthSciDiscretizations._neighbor_index(grid, p, i+1, j)
    pim1 = EarthSciDiscretizations._neighbor_index(grid, p, i-1, j)

    # Diagonal neighbor: (i+1, j+1) — should be physically NE of cell (p,i,j)
    pip1jp1 = EarthSciDiscretizations._diagonal_neighbor(grid, p, i, j, +1, 0, +1, pip1)

    # The diagonal neighbor should be a valid cell
    @test 1 <= pip1jp1[1] <= 6
    @test 1 <= pip1jp1[2] <= Nc
    @test 1 <= pip1jp1[3] <= Nc

    # Physical location should be approximately NE
    center_cart = gnomonic_to_cart(grid.ξ_centers[i], grid.η_centers[j], p)
    diag_cart = gnomonic_to_cart(grid.ξ_centers[pip1jp1[2]], grid.η_centers[pip1jp1[3]], pip1jp1[1])

    # The diagonal should be at a reasonable distance (not on the opposite side of the sphere)
    dist = acos(clamp(dot(center_cart, diag_cart), -1.0, 1.0))
    @test dist < π / 2  # Less than 90° away
    @test dist > 0      # Not the same cell
end

# ============================================================
# Task 5: PPM Courant number
# ============================================================

@testitem "PPM CFL warning fires for large Courant" setup=[FixesSetup] tags=[:fixes] begin
    # ppm_flux_integral should warn when |Courant| > 1
    ql = 1.0; qr = 2.0; qi = 1.5
    @test_logs (:warn, r"CFL violation") EarthSciDiscretizations.ppm_flux_integral(ql, qr, qi, 1.5)
    # Should still return a finite value
    @test isfinite(EarthSciDiscretizations.ppm_flux_integral(ql, qr, qi, 1.5))
end

# ============================================================
# Task 7: Precomputed FV stencils
# ============================================================

@testitem "Precomputed Laplacian stencil has correct structure" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    stencil = precompute_laplacian_stencil(grid)

    @test size(stencil.weights) == (6, Nc, Nc, 9)
    @test size(stencil.nb_p) == (6, Nc, Nc, 9)

    # All neighbor indices should be valid
    @test all(1 .<= stencil.nb_p .<= 6)
    @test all(1 .<= stencil.nb_i .<= Nc)
    @test all(1 .<= stencil.nb_j .<= Nc)

    # Center neighbor (index 1) should be self-referencing
    for p in 1:6, i in 1:Nc, j in 1:Nc
        @test stencil.nb_p[p, i, j, 1] == p
        @test stencil.nb_i[p, i, j, 1] == i
        @test stencil.nb_j[p, i, j, 1] == j
    end
end

@testitem "Precomputed Laplacian of constant is zero" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    stencil = precompute_laplacian_stencil(grid)

    u = ones(6, Nc, Nc) * 42.0
    du = zeros(6, Nc, Nc)
    apply_laplacian!(du, u, stencil)

    @test maximum(abs.(du)) < 1e-10
end

@testitem "Precomputed gradient stencil basic test" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    stencil = precompute_gradient_stencil(grid)

    @test size(stencil.weights_lon) == (6, Nc, Nc, 5)
    @test size(stencil.weights_lat) == (6, Nc, Nc, 5)

    # Gradient of constant should be zero
    u = ones(6, Nc, Nc) * 42.0
    du_lon = zeros(6, Nc, Nc)
    du_lat = zeros(6, Nc, Nc)
    apply_gradient!(du_lon, du_lat, u, stencil)
    @test maximum(abs.(du_lon)) < 1e-10
    @test maximum(abs.(du_lat)) < 1e-10
end

# ============================================================
# Task 7: Ghost-extended ArrayOp Laplacian
# ============================================================

@testitem "Extended Laplacian of constant is zero" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    phi = ones(6, Nc, Nc) * 7.0
    phi_ext = extend_with_ghosts(phi, grid)
    lapl = fv_laplacian_extended(phi_ext, grid)
    result = evaluate_arrayop(lapl)
    @test size(result) == (6, Nc, Nc)
    @test maximum(abs.(result)) < 1e-10
end

@testitem "Extended Laplacian covers all cells including boundaries" setup=[FixesSetup] tags=[:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    phi = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
    end
    phi_ext = extend_with_ghosts(phi, grid)
    lapl = fv_laplacian_extended(phi_ext, grid)
    result = evaluate_arrayop(lapl)

    # Should have (6, Nc, Nc) output covering ALL cells
    @test size(result) == (6, Nc, Nc)
    # All values should be finite
    @test all(isfinite, result)
    # Should not be zero everywhere
    @test maximum(abs.(result)) > 1e-3
end

# ============================================================
# Task 8: Double-Laplacian deduplication
# ============================================================

@testitem "∂²u/∂lon² + ∂²u/∂lat² produces single Laplacian" setup=[FixesSetup] begin
    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using Symbolics
    using DomainSets
    using OrdinaryDiffEqDefault
    using SciMLBase

    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    bcs = [u(0, lon, lat) ~ exp(-10 * (lon^2 + lat^2))]
    domains = [t ∈ Interval(0.0, 0.1),
               lon ∈ Interval(-π, π),
               lat ∈ Interval(-π/2, π/2)]

    # Each ∂²u/∂x² maps to (1/n_spatial) of the full covariant Laplacian.
    # So ∂²u/∂lon² + ∂²u/∂lat² (n_spatial=2) gives 1/2 + 1/2 = 1 Laplacian,
    # and a single ∂²u/∂lon² gives 1/2 Laplacian.
    # Therefore, RHS of (both with κ) should be 2× RHS of (one with κ).
    eq_one  = [D(u(t, lon, lat)) ~ 0.1 * Dlon(Dlon(u(t, lon, lat)))]
    eq_both = [D(u(t, lon, lat)) ~ 0.1 * (Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat))))]

    @named sys_one  = PDESystem(eq_one,  bcs, domains, [t, lon, lat], [u(t, lon, lat)])
    @named sys_both = PDESystem(eq_both, bcs, domains, [t, lon, lat], [u(t, lon, lat)])

    Nc = 4
    disc = FVCubedSphere(Nc; R=1.0)

    prob_one  = discretize(sys_one, disc)
    prob_both = discretize(sys_both, disc)

    du_one  = similar(prob_one.u0);  prob_one.f(du_one, prob_one.u0, prob_one.p, 0.0)
    du_both = similar(prob_both.u0); prob_both.f(du_both, prob_both.u0, prob_both.p, 0.0)

    ratio = maximum(abs.(du_both)) / maximum(abs.(du_one))
    @test isapprox(ratio, 2.0, rtol=0.1)

    sol = solve(prob_both)
    @test sol.retcode == SciMLBase.ReturnCode.Success
end
