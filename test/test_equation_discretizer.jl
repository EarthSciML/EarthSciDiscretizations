@testsnippet DiscretizerSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciDiscretizations: _neighbor_index, _identify_lhs_dv, _get_dv_value,
        _is_initial_condition
    using ModelingToolkit
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using Symbolics
    using DomainSets
    using SciMLBase
    using OrdinaryDiffEqDefault
end

@testitem "Ghost cell neighbor index resolution" setup = [DiscretizerSetup] tags = [:discretizer] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc

    # Interior indices should pass through unchanged
    @test _neighbor_index(grid, 1, 4, 4) == (1, 4, 4)
    @test _neighbor_index(grid, 3, 1, 1) == (3, 1, 1)
    @test _neighbor_index(grid, 6, Nc, Nc) == (6, Nc, Nc)

    # Boundary indices should map to neighbor panels
    # West boundary of panel 1 (i=0) -> panel 5 East edge
    p, i, j = _neighbor_index(grid, 1, 0, 4)
    @test p == 5  # Panel 5 is West neighbor of panel 1
    @test 1 <= i <= Nc
    @test 1 <= j <= Nc

    # East boundary of panel 1 (i=Nc+1) -> panel 2 West edge
    p, i, j = _neighbor_index(grid, 1, Nc + 1, 4)
    @test p == 2
    @test 1 <= i <= Nc
    @test 1 <= j <= Nc

    # South boundary of panel 1 (j=0) -> panel 6
    p, i, j = _neighbor_index(grid, 1, 4, 0)
    @test p == 6
    @test 1 <= i <= Nc
    @test 1 <= j <= Nc

    # North boundary of panel 1 (j=Nc+1) -> panel 3
    p, i, j = _neighbor_index(grid, 1, 4, Nc + 1)
    @test p == 3
    @test 1 <= i <= Nc
    @test 1 <= j <= Nc

    # Verify all boundary cells on all panels resolve to valid indices
    for panel in 1:6
        for idx in 1:Nc
            for boundary_val in [0, Nc + 1]
                # ξ-boundaries
                p, i, j = _neighbor_index(grid, panel, boundary_val, idx)
                @test 1 <= p <= 6
                @test 1 <= i <= Nc
                @test 1 <= j <= Nc
                # η-boundaries
                p, i, j = _neighbor_index(grid, panel, idx, boundary_val)
                @test 1 <= p <= 6
                @test 1 <= i <= Nc
                @test 1 <= j <= Nc
            end
        end
    end
end

@testitem "LHS DV identification" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..) v(..)

    dvs = [u(t, lon, lat), v(t, lon, lat)]

    # Identify u from D(u(t, lon, lat))
    lhs_u = D(u(t, lon, lat))
    @test _identify_lhs_dv(lhs_u, dvs) === dvs[1]

    # Identify v from D(v(t, lon, lat))
    lhs_v = D(v(t, lon, lat))
    @test _identify_lhs_dv(lhs_v, dvs) === dvs[2]
end

@testitem "Dimension identification" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @test identify_dimension(:t) == :t
    @test identify_dimension(:lon) == :lon
    @test identify_dimension(:lat) == :lat
    @test identify_dimension(:x) == :xi
    @test identify_dimension(:y) == :eta
    @test identify_dimension(:z) == :vertical
    # Physical coordinates are distinguished from computational ones
    @test identify_dimension(:ξ) == :xi
    @test identify_dimension(:η) == :eta
    @test identify_dimension(:λ) == :lon
    @test identify_dimension(:φ) == :lat
end

@testitem "Multi-variable discretization" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..) v(..)
    Dlon = Differential(lon)

    # Coupled system: du/dt = v, dv/dt = -u
    eq1 = D(u(t, lon, lat)) ~ v(t, lon, lat)
    eq2 = D(v(t, lon, lat)) ~ -u(t, lon, lat)

    bcs = [
        u(0, lon, lat) ~ cos(lat),
        v(0, lon, lat) ~ 0.0,
    ]
    domains = [
        t ∈ Interval(0.0, 1.0),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        [eq1, eq2], bcs, domains,
        [t, lon, lat], [u(t, lon, lat), v(t, lon, lat)]
    )

    disc = FVCubedSphere(4; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    @test prob isa ODEProblem
    # Two variables × 6 panels × 4 × 4 = 192
    @test length(prob.u0) == 2 * 6 * 4 * 4

    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success
end

@testitem "Nonlinear derivative terms" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)

    # PDE with nonlinear term: du/dt = d(u^2)/dlon
    eq = [D(u(t, lon, lat)) ~ Dlon(u(t, lon, lat)^2)]
    bcs = [u(0, lon, lat) ~ 1.0 + 0.1 * cos(lon)]
    domains = [
        t ∈ Interval(0.0, 0.01),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        eq, bcs, domains,
        [t, lon, lat], [u(t, lon, lat)]
    )

    disc = FVCubedSphere(4; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    @test prob isa ODEProblem

    # The RHS should not be identically zero (which would happen if
    # the nonlinearity were discarded)
    du = prob.f(prob.u0, prob.p, 0.0)
    @test !all(iszero, du)
end

@testitem "Coupled variable stencils" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..) v(..)
    Dlon = Differential(lon)

    # du/dt = dv/dlon — the stencil for u's equation must reference v's array
    eq1 = D(u(t, lon, lat)) ~ Dlon(v(t, lon, lat))
    eq2 = D(v(t, lon, lat)) ~ Dlon(u(t, lon, lat))

    bcs = [
        u(0, lon, lat) ~ cos(lon) * cos(lat),
        v(0, lon, lat) ~ sin(lon) * cos(lat),
    ]
    domains = [
        t ∈ Interval(0.0, 0.01),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        [eq1, eq2], bcs, domains,
        [t, lon, lat], [u(t, lon, lat), v(t, lon, lat)]
    )

    disc = FVCubedSphere(4; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    @test prob isa ODEProblem
    @test length(prob.u0) == 2 * 6 * 4 * 4

    # RHS should be nonzero (v varies spatially)
    du = prob.f(prob.u0, prob.p, 0.0)
    @test !all(iszero, du)
end

@testitem "Metric-corrected Laplacian of constant is zero" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    # Laplacian of a constant should be zero even with metric corrections
    eq = [D(u(t, lon, lat)) ~ Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat)))]
    bcs = [u(0, lon, lat) ~ 42.0]  # Constant IC
    domains = [
        t ∈ Interval(0.0, 1.0),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        eq, bcs, domains,
        [t, lon, lat], [u(t, lon, lat)]
    )

    disc = FVCubedSphere(8; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    # Evaluate RHS at t=0 with constant field
    du = prob.f(prob.u0, prob.p, 0.0)
    @test maximum(abs, du) < 1.0e-10
end

@testitem "Mixed derivative discretization" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    # PDE with mixed derivative: du/dt = d²u/(dlon dlat)
    eq = [D(u(t, lon, lat)) ~ Dlon(Dlat(u(t, lon, lat)))]
    bcs = [u(0, lon, lat) ~ sin(lon) * sin(lat)]
    domains = [
        t ∈ Interval(0.0, 0.01),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        eq, bcs, domains,
        [t, lon, lat], [u(t, lon, lat)]
    )

    disc = FVCubedSphere(8; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    @test prob isa ODEProblem

    # The mixed derivative of sin(lon)*sin(lat) = cos(lon)*cos(lat), which is nonzero
    du = prob.f(prob.u0, prob.p, 0.0)
    @test !all(iszero, du)
end

@testitem "IC identification with multiple BCs" setup = [DiscretizerSetup] tags = [:discretizer] begin
    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    eq = [D(u(t, lon, lat)) ~ 0.01 * (Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat))))]

    # Multiple BCs — IC is listed second
    bcs = [
        u(t, 0, lat) ~ 0.0,        # Spatial BC (not an IC)
        u(0, lon, lat) ~ cos(lat),
    ]   # IC

    domains = [
        t ∈ Interval(0.0, 1.0),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        eq, bcs, domains,
        [t, lon, lat], [u(t, lon, lat)]
    )

    disc = FVCubedSphere(4; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    # Should use the cos(lat) IC, not the zero BC
    @test maximum(prob.u0) > 0.5
end

@testitem "Grid has metric tensor arrays" setup = [DiscretizerSetup] tags = [:discretizer] begin
    grid = CubedSphereGrid(8; R = 1.0)
    Nc = grid.Nc

    @test size(grid.J) == (6, Nc, Nc)
    @test size(grid.ginv_ξξ) == (6, Nc, Nc)
    @test size(grid.ginv_ηη) == (6, Nc, Nc)
    @test size(grid.ginv_ξη) == (6, Nc, Nc)

    # Second-derivative Jacobian arrays
    @test size(grid.d2ξ_dlon2) == (6, Nc, Nc)
    @test size(grid.d2ξ_dlondlat) == (6, Nc, Nc)
    @test size(grid.d2ξ_dlat2) == (6, Nc, Nc)
    @test size(grid.d2η_dlon2) == (6, Nc, Nc)
    @test size(grid.d2η_dlondlat) == (6, Nc, Nc)
    @test size(grid.d2η_dlat2) == (6, Nc, Nc)

    # Jacobian should be positive everywhere
    @test all(grid.J .> 0)

    # Inverse metric diagonal terms should be positive
    @test all(grid.ginv_ξξ .> 0)
    @test all(grid.ginv_ηη .> 0)

    # At panel center (ξ=0, η=0), the metric should be isotropic
    # and cross-term should be zero
    center = Nc ÷ 2
    for p in 1:6
        @test isapprox(
            grid.ginv_ξξ[p, center, center],
            grid.ginv_ηη[p, center, center]; rtol = 0.1
        )
        @test abs(grid.ginv_ξη[p, center, center]) < 0.1
    end
end

@testitem "Forward Jacobian is well-conditioned everywhere" setup = [DiscretizerSetup] tags = [:discretizer] begin
    grid = CubedSphereGrid(16; R = 1.0)
    Nc = grid.Nc

    for p in 1:6, i in 1:Nc, j in 1:Nc
        fwd = compute_forward_jacobian(grid.ξ_centers[i], grid.η_centers[j], p)
        # Forward Jacobian should always be finite
        @test isfinite(fwd.dlon_dξ)
        @test isfinite(fwd.dlon_dη)
        @test isfinite(fwd.dlat_dξ)
        @test isfinite(fwd.dlat_dη)
    end

    # Inverse Jacobian should be finite after regularization (even near poles)
    for p in [3, 6], i in 1:Nc, j in 1:Nc
        jac = compute_coord_jacobian(grid.ξ_centers[i], grid.η_centers[j], p)
        @test isfinite(jac.dξ_dlon)
        @test isfinite(jac.dξ_dlat)
        @test isfinite(jac.dη_dlon)
        @test isfinite(jac.dη_dlat)
    end
end

@testitem "Second-derivative chain rule: Laplacian of cos(lat) on equatorial panel" setup = [DiscretizerSetup] tags = [:discretizer] begin
    # Test that ∂²cos(lat)/∂lon² + ∂²cos(lat)/∂lat² gives the correct value.
    # For f = cos(lat):
    #   ∂f/∂lon = 0, so ∂²f/∂lon² = 0
    #   ∂f/∂lat = -sin(lat), ∂²f/∂lat² = -cos(lat)
    @parameters lon lat
    @variables u(..)
    Dlon = Differential(lon)
    Dlat = Differential(lat)

    eq = [D(u(t, lon, lat)) ~ Dlon(Dlon(u(t, lon, lat))) + Dlat(Dlat(u(t, lon, lat)))]
    bcs = [u(0, lon, lat) ~ cos(lat)]
    domains = [
        t ∈ Interval(0.0, 0.01),
        lon ∈ Interval(-π, π),
        lat ∈ Interval(-π / 2, π / 2),
    ]
    @named sys = PDESystem(
        eq, bcs, domains,
        [t, lon, lat], [u(t, lon, lat)]
    )

    disc = FVCubedSphere(8; R = 1.0)
    prob = SciMLBase.discretize(sys, disc)

    du = prob.f(prob.u0, prob.p, 0.0)

    # The ODE state is flattened; check that the tendency is non-trivial.
    # For cos(lat), ∂²f/∂lat² = -cos(lat), so du/dt should have O(1) magnitude.
    # By the divergence theorem, ∫_sphere ∇²f dA = 0, so the global sum should vanish.
    @test !all(iszero, du)
    @test maximum(abs.(du)) > 0.1  # Tendency should be O(1) for a unit-sphere Laplacian
    # Global integral of Laplacian vanishes (divergence theorem)
    @test abs(sum(du)) / (length(du) * maximum(abs.(du))) < 0.1
end
