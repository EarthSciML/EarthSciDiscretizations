@testsnippet ArrayOpSetup begin
    using Test
    using EarthSciDiscretizations
    using SymbolicUtils: SymReal, isarrayop
    using SymbolicUtils
    using Symbolics: unwrap, wrap
end

@testitem "All operators return ArrayOps" setup=[ArrayOpSetup] tags=[:arrayop] begin
    Nc = 4
    grid = CubedSphereGrid(Nc; R=1.0)

    F_xi = zeros(6, Nc + 1, Nc)
    F_eta = zeros(6, Nc, Nc + 1)
    phi = zeros(6, Nc, Nc)
    c_xi = zeros(6, Nc + 1, Nc)
    c_eta = zeros(6, Nc, Nc + 1)

    @test isarrayop(fv_divergence(F_xi, F_eta, grid))
    @test isarrayop(fv_gradient_xi(phi, grid))
    @test isarrayop(fv_gradient_eta(phi, grid))
    @test isarrayop(fv_laplacian(phi, grid))
    @test isarrayop(flux_1d(phi, c_xi, grid, :xi))
    @test isarrayop(flux_1d(phi, c_eta, grid, :eta))
    @test isarrayop(transport_2d(phi, c_xi, c_eta, grid))
    @test isarrayop(ghost_fill_arrayop(phi, grid))
end

@testitem "evaluate_arrayop round-trips Const data" setup=[ArrayOpSetup] tags=[:arrayop] begin
    arr = rand(3, 4)
    c = const_wrap(arr)
    idx = get_idx_vars(2)
    i, j = idx[1], idx[2]
    expr = unwrap(wrap(c[i, j]))
    ao = make_arrayop([i, j], expr, Dict(i => 1:1:3, j => 1:1:4))
    result = evaluate_arrayop(ao)
    @test isapprox(result, arr; rtol=1e-14)
end
