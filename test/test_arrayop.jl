@testsnippet ArrayOpSetup begin
    using Test
    using EarthSciDiscretizations
    using SymbolicUtils: SymReal, isarrayop
    using SymbolicUtils
    using Symbolics: unwrap, wrap
end

@testitem "All operators return ArrayOps" setup = [ArrayOpSetup] tags = [:arrayop] begin
    Nc = 4
    grid = CubedSphereGrid(Nc; R = 1.0)

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

    # ghost_fill_arrayop returns a ghost-extended array (not an ArrayOp)
    q_ext = ghost_fill_arrayop(phi, grid)
    Ng = grid.Ng
    @test size(q_ext) == (6, Nc + 2Ng, Nc + 2Ng)

    # New ArrayOp operators
    courant_xi = compute_courant_numbers(c_xi, 0.01, grid, :xi)
    courant_eta = compute_courant_numbers(c_eta, 0.01, grid, :eta)
    @test isarrayop(compute_courant_numbers_arrayop(c_xi, 0.01, grid, :xi))
    @test isarrayop(compute_courant_numbers_arrayop(c_eta, 0.01, grid, :eta))
    @test isarrayop(transport_2d_ppm_arrayop(q_ext, courant_xi, courant_eta, c_xi, c_eta, grid))
    ql, qr = ppm_reconstruction_arrayop(q_ext, grid, :xi)
    @test isarrayop(ql)
    @test isarrayop(qr)
end

@testitem "evaluate_arrayop round-trips Const data" setup = [ArrayOpSetup] tags = [:arrayop] begin
    arr = rand(3, 4)
    c = const_wrap(arr)
    idx = get_idx_vars(2)
    i, j = idx[1], idx[2]
    expr = unwrap(wrap(c[i, j]))
    ao = make_arrayop([i, j], expr, Dict(i => 1:1:3, j => 1:1:4))
    result = evaluate_arrayop(ao)
    @test isapprox(result, arr; rtol = 1.0e-14)
end
