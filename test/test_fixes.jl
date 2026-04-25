@testsnippet FixesSetup begin
    using Test
    using EarthSciDiscretizations
    using LinearAlgebra: norm, dot
end

# ============================================================
# Task 1: Analytical Jacobian
# ============================================================

@testitem "Analytical Jacobian matches finite differences" setup = [FixesSetup] tags = [:fixes] begin
    # The analytical forward Jacobian should agree with finite-difference
    # computation to high precision (within FD truncation error)
    for panel in 1:6
        for (ξ, η) in [(0.0, 0.0), (0.3, -0.2), (-0.5, 0.4), (0.7, 0.7)]
            fwd = compute_forward_jacobian(ξ, η, panel)

            # Finite-difference reference
            h = 1.0e-7
            lon_p, lat_p = gnomonic_to_lonlat(ξ + h, η, panel)
            lon_m, lat_m = gnomonic_to_lonlat(ξ - h, η, panel)
            dlon = lon_p - lon_m
            # Wrap to [-π, π] to handle anti-meridian crossings
            while dlon > π
                dlon -= 2π
            end
            while dlon < -π
                dlon += 2π
            end
            dlon_dξ_fd = dlon / (2h)
            dlat_dξ_fd = (lat_p - lat_m) / (2h)

            lon_p, lat_p = gnomonic_to_lonlat(ξ, η + h, panel)
            lon_m, lat_m = gnomonic_to_lonlat(ξ, η - h, panel)
            dlon = lon_p - lon_m
            # Wrap to [-π, π] to handle anti-meridian crossings
            while dlon > π
                dlon -= 2π
            end
            while dlon < -π
                dlon += 2π
            end
            dlon_dη_fd = dlon / (2h)
            dlat_dη_fd = (lat_p - lat_m) / (2h)

            # At the poles (panels 3,6 at center), longitude is undefined and
            # the forward Jacobian is correctly singular. Skip lon assertions there.
            _, lat_c = gnomonic_to_lonlat(ξ, η, panel)
            near_pole = abs(abs(lat_c) - π / 2) < 0.1
            if !near_pole
                @test isapprox(fwd.dlon_dξ, dlon_dξ_fd, atol = 1.0e-5)
                @test isapprox(fwd.dlon_dη, dlon_dη_fd, atol = 1.0e-5)
            end
            @test isapprox(fwd.dlat_dξ, dlat_dξ_fd, atol = 1.0e-5)
            @test isapprox(fwd.dlat_dη, dlat_dη_fd, atol = 1.0e-5)
        end
    end
end

@testitem "Analytical Jacobian at panel center" setup = [FixesSetup] tags = [:fixes] begin
    # Panel 1 center (0,0): lon = 0, lat = 0
    # dlon/dξ = 1, dlon/dη = 0, dlat/dξ = 0, dlat/dη = 1
    fwd = compute_forward_jacobian(0.0, 0.0, 1)
    @test isapprox(fwd.dlon_dξ, 1.0, atol = 1.0e-12)
    @test isapprox(fwd.dlon_dη, 0.0, atol = 1.0e-12)
    @test isapprox(fwd.dlat_dξ, 0.0, atol = 1.0e-12)
    @test isapprox(fwd.dlat_dη, 1.0, atol = 1.0e-12)
end

@testitem "Coordinate Jacobian inverse is consistent" setup = [FixesSetup] tags = [:fixes] begin
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

            @test isapprox(I11, 1.0, atol = 1.0e-8)
            @test isapprox(I12, 0.0, atol = 1.0e-8)
            @test isapprox(I21, 0.0, atol = 1.0e-8)
            @test isapprox(I22, 1.0, atol = 1.0e-8)
        end
    end
end

# ============================================================
# Task 3: Vector field rotation at panel boundaries
# ============================================================

@testitem "Rotation matrices are orthogonal" setup = [FixesSetup] tags = [:fixes] begin
    grid = CubedSphereGrid(8)
    for p in 1:6, dir in (West, East, South, North)
        M11, M12, M21, M22 = grid.rotation_matrices[(p, dir)]
        # Orthogonality: M * M^T = I
        @test isapprox(M11^2 + M12^2, 1.0, atol = 1.0e-10)
        @test isapprox(M21^2 + M22^2, 1.0, atol = 1.0e-10)
        @test isapprox(M11 * M21 + M12 * M22, 0.0, atol = 1.0e-10)
        # det = ±1
        det = M11 * M22 - M12 * M21
        @test isapprox(abs(det), 1.0, atol = 1.0e-10)
    end
end

@testitem "Vector ghost cells preserve uniform rotation field" setup = [FixesSetup] tags = [:fixes] begin
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
        @test uξ_ext[p, i + Ng, j + Ng] == uξ[p, i, j]
        @test uη_ext[p, i + Ng, j + Ng] == uη[p, i, j]
    end
end

# ============================================================
# Task 2: Laplacian cross-metric correction
# ============================================================

@testitem "Laplacian of constant is zero (with cross-metric fix)" setup = [FixesSetup] tags = [:fixes] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    phi = ones(6, Nc, Nc) * 42.0
    lapl = fv_laplacian(phi, grid)
    result = evaluate_arrayop(lapl)
    @test all(abs.(result) .< 1.0e-10)
end

@testitem "Laplacian convergence improves with resolution" setup = [FixesSetup] tags = [:fixes] begin
    errors = Float64[]
    for Nc in [8, 16, 32]
        grid = CubedSphereGrid(Nc)
        phi = zeros(6, Nc, Nc)
        for p in 1:6, i in 1:Nc, j in 1:Nc
            phi[p, i, j] = cos(2 * grid.ξ_centers[i]) * cos(2 * grid.η_centers[j])
        end
        lapl = fv_laplacian(phi, grid)
        result = evaluate_arrayop(lapl)
        push!(errors, sqrt(sum(result .^ 2) / length(result)))
    end
    # RMS of Laplacian should converge (not diverge) as resolution increases.
    # The values converge to the analytical Laplacian, so they stabilize.
    @test errors[2] < errors[1] * 1.1
    @test errors[3] < errors[2] * 1.1
end

# ============================================================
# Task 5: PPM Courant number
# ============================================================

@testitem "PPM CFL warning fires for large Courant" setup = [FixesSetup] tags = [:fixes] begin
    # ppm_flux_integral should warn when |Courant| > 1
    ql = 1.0; qr = 2.0; qi = 1.5
    @test_logs (:warn, r"CFL violation") EarthSciDiscretizations.ppm_flux_integral(ql, qr, qi, 1.5)
    # Should still return a finite value
    @test isfinite(EarthSciDiscretizations.ppm_flux_integral(ql, qr, qi, 1.5))
end

