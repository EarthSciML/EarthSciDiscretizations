@testsnippet Transport2DSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
end

@testitem "2D Lax-Friedrichs transport of constant field with zero velocity is zero" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q = fill(7.0, 6, Nc, Nc)
    courant_xi = fill(0.0, 6, Nc, Nc)
    courant_eta = fill(0.0, 6, Nc, Nc)

    ao = transport_2d(q, courant_xi, courant_eta, grid)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc - 2, Nc - 2)
    @test all(abs.(tendency) .< 1.0e-12)
end

@testitem "2D Lax-Friedrichs transport is linear in q" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)

    q1 = randn(6, Nc, Nc)
    q2 = randn(6, Nc, Nc)
    courant_xi = fill(0.1, 6, Nc, Nc)
    courant_eta = fill(0.1, 6, Nc, Nc)

    ao1 = transport_2d(q1, courant_xi, courant_eta, grid)
    ao2 = transport_2d(q2, courant_xi, courant_eta, grid)
    ao_sum = transport_2d(q1 + q2, courant_xi, courant_eta, grid)

    t1 = evaluate_arrayop(ao1)
    t2 = evaluate_arrayop(ao2)
    t_sum = evaluate_arrayop(ao_sum)

    @test isapprox(t_sum, t1 + t2; rtol = 1.0e-10)
end

@testitem "Lin-Rood 2D transport of constant field with zero velocity is zero" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = fill(5.0, 6, Nc, Nc)
    vel_xi = fill(0.0, 6, Nc + 1, Nc)
    vel_eta = fill(0.0, 6, Nc, Nc + 1)
    dt = 0.01

    tendency = zeros(6, Nc, Nc)
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)

    @test all(abs.(tendency) .< 1.0e-14)
end

@testitem "Lin-Rood 2D transport conserves mass" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 16
    grid = CubedSphereGrid(Nc; R = 1.0)

    # Non-trivial smooth initial condition
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * cos(2π * grid.ξ_centers[i] / (π / 2)) *
            cos(2π * grid.η_centers[j] / (π / 2))
    end

    # Non-zero velocity in both directions
    vel_xi = fill(0.5, 6, Nc + 1, Nc)
    vel_eta = fill(0.3, 6, Nc, Nc + 1)
    dt = 0.01

    tendency = zeros(6, Nc, Nc)
    transport_2d_linrood!(tendency, q, vel_xi, vel_eta, grid, dt)

    # Total mass tendency should be zero (conservation on closed sphere)
    # All boundary fluxes are matched, so this should hold to machine precision
    total = sum(tendency[p, i, j] * grid.area[p, i, j] for p in 1:6, i in 1:Nc, j in 1:Nc)
    @test abs(total) < 1.0e-10
end

@testitem "PPM ArrayOp produces ArrayOp type" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    using SymbolicUtils: isarrayop

    Nc = 4
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = zeros(6, Nc, Nc)
    vel_xi = fill(0.1, 6, Nc + 1, Nc)
    dt = 0.01

    q_ext = extend_with_ghosts(q, grid)
    courant = compute_courant_numbers(vel_xi, dt, grid, :xi)

    ao = flux_1d_ppm_arrayop(q_ext, courant, vel_xi, grid, :xi)
    @test isarrayop(ao)
end

@testitem "PPM ArrayOp of constant field with zero velocity is zero" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = fill(3.0, 6, Nc, Nc)
    vel_xi = fill(0.0, 6, Nc + 1, Nc)
    dt = 0.01

    q_ext = extend_with_ghosts(q, grid)
    courant = compute_courant_numbers(vel_xi, dt, grid, :xi)

    ao = flux_1d_ppm_arrayop(q_ext, courant, vel_xi, grid, :xi)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc, Nc)
    @test all(abs.(tendency) .< 1.0e-10)
end

@testitem "PPM ArrayOp matches loop-based PPM" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    # Smooth non-trivial field
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * sin(2π * grid.ξ_centers[i] / (π / 2))
    end

    vel_xi = fill(0.3, 6, Nc + 1, Nc)
    dt = 0.01

    # Loop-based PPM
    tend_loop = zeros(6, Nc, Nc)
    flux_1d_ppm!(tend_loop, q, vel_xi, grid, :xi, dt)

    # ArrayOp-based PPM
    q_ext = extend_with_ghosts(q, grid)
    courant = compute_courant_numbers(vel_xi, dt, grid, :xi)
    ao = flux_1d_ppm_arrayop(q_ext, courant, vel_xi, grid, :xi)
    tend_arrayop = evaluate_arrayop(ao)

    # Interior cells should match (boundary cells may differ due to ghost handling)
    # Compare cells that are well away from panel boundaries
    for p in 1:6, i in 3:(Nc - 2), j in 3:(Nc - 2)
        @test isapprox(tend_arrayop[p, i, j], tend_loop[p, i, j]; rtol = 1.0e-10)
    end
end

@testitem "flux_to_tendency_arrayop matches loop version" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    # Create some flux data
    flux = randn(6, Nc + 1, Nc)

    ao = flux_to_tendency_arrayop(flux, grid, :xi)
    tend_ao = evaluate_arrayop(ao)

    # Compare with loop-based version
    tend_loop = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        tend_loop[p, i, j] = -(
            flux[p, i + 1, j] * grid.dx[p, i + 1, j] -
                flux[p, i, j] * grid.dx[p, i, j]
        ) / grid.area[p, i, j]
    end

    @test isapprox(tend_ao, tend_loop; rtol = 1.0e-12)
end

@testitem "advective_tendency_arrayop matches loop version" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    tend_flux = randn(6, Nc, Nc)
    q = randn(6, Nc, Nc)
    vel = randn(6, Nc + 1, Nc)

    ao = advective_tendency_arrayop(tend_flux, q, vel, grid, :xi)
    tend_ao = evaluate_arrayop(ao)

    # Compare with loop-based version
    tend_loop = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        c_def = (
            vel[p, i + 1, j] * grid.dx[p, i + 1, j] -
                vel[p, i, j] * grid.dx[p, i, j]
        ) / grid.area[p, i, j]
        tend_loop[p, i, j] = tend_flux[p, i, j] + q[p, i, j] * c_def
    end

    @test isapprox(tend_ao, tend_loop; rtol = 1.0e-12)
end

@testitem "compute_courant_numbers_arrayop matches loop version" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)
    dt = 0.01

    # ξ-direction
    vel_xi = randn(6, Nc + 1, Nc) * 0.1
    courant_loop = compute_courant_numbers(vel_xi, dt, grid, :xi)
    ao_xi = compute_courant_numbers_arrayop(vel_xi, dt, grid, :xi)
    courant_ao = evaluate_arrayop(ao_xi)
    @test isapprox(courant_ao, courant_loop; rtol = 1.0e-12)

    # η-direction
    vel_eta = randn(6, Nc, Nc + 1) * 0.1
    courant_loop_eta = compute_courant_numbers(vel_eta, dt, grid, :eta)
    ao_eta = compute_courant_numbers_arrayop(vel_eta, dt, grid, :eta)
    courant_ao_eta = evaluate_arrayop(ao_eta)
    @test isapprox(courant_ao_eta, courant_loop_eta; rtol = 1.0e-12)
end

@testitem "transport_2d_ppm_arrayop of constant field with zero velocity is zero" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = fill(5.0, 6, Nc, Nc)
    vel_xi = fill(0.0, 6, Nc + 1, Nc)
    vel_eta = fill(0.0, 6, Nc, Nc + 1)
    dt = 0.01

    q_ext = extend_with_ghosts(q, grid)
    courant_xi = compute_courant_numbers(vel_xi, dt, grid, :xi)
    courant_eta = compute_courant_numbers(vel_eta, dt, grid, :eta)

    ao = transport_2d_ppm_arrayop(q_ext, courant_xi, courant_eta, vel_xi, vel_eta, grid)
    tendency = evaluate_arrayop(ao)

    @test size(tendency) == (6, Nc, Nc)
    @test all(abs.(tendency) .< 1.0e-10)
end

@testitem "transport_2d_ppm_arrayop matches sum of 1D PPM ArrayOps" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    # Smooth non-trivial field
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 1.0 + 0.5 * sin(2π * grid.ξ_centers[i] / (π / 2)) *
            cos(2π * grid.η_centers[j] / (π / 2))
    end

    vel_xi = fill(0.3, 6, Nc + 1, Nc)
    vel_eta = fill(0.2, 6, Nc, Nc + 1)
    dt = 0.01

    q_ext = extend_with_ghosts(q, grid)
    courant_xi = compute_courant_numbers(vel_xi, dt, grid, :xi)
    courant_eta = compute_courant_numbers(vel_eta, dt, grid, :eta)

    # 2D ArrayOp
    ao_2d = transport_2d_ppm_arrayop(q_ext, courant_xi, courant_eta, vel_xi, vel_eta, grid)
    tend_2d = evaluate_arrayop(ao_2d)

    # Sum of 1D ArrayOps
    ao_xi = flux_1d_ppm_arrayop(q_ext, courant_xi, vel_xi, grid, :xi)
    ao_eta = flux_1d_ppm_arrayop(q_ext, courant_eta, vel_eta, grid, :eta)
    tend_xi = evaluate_arrayop(ao_xi)
    tend_eta = evaluate_arrayop(ao_eta)

    @test isapprox(tend_2d, tend_xi + tend_eta; rtol = 1.0e-10)
end

@testitem "ppm_reconstruction_arrayop produces ArrayOp types" setup = [Transport2DSetup] tags = [:transport, :arrayop] begin
    using SymbolicUtils: isarrayop

    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = randn(6, Nc, Nc)
    q_ext = extend_with_ghosts(q, grid)

    ql_ao, qr_ao = ppm_reconstruction_arrayop(q_ext, grid, :xi)
    @test isarrayop(ql_ao)
    @test isarrayop(qr_ao)

    ql_ao_eta, qr_ao_eta = ppm_reconstruction_arrayop(q_ext, grid, :eta)
    @test isarrayop(ql_ao_eta)
    @test isarrayop(qr_ao_eta)
end

@testitem "ghost_fill_arrayop returns extended array" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = randn(6, Nc, Nc)
    q_ext = ghost_fill_arrayop(q, grid)
    q_ext_ref = extend_with_ghosts(q, grid)

    @test size(q_ext) == size(q_ext_ref)
    @test q_ext ≈ q_ext_ref
end
