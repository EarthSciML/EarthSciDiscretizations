@testsnippet Transport2DSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: evaluate_arrayop
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


@testitem "ghost_fill_arrayop returns extended array" setup = [Transport2DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc; R = 1.0)

    q = randn(6, Nc, Nc)
    q_ext = ghost_fill_arrayop(q, grid)
    q_ext_ref = extend_with_ghosts(q, grid)

    @test size(q_ext) == size(q_ext_ref)
    @test q_ext ≈ q_ext_ref
end
