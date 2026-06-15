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

    phi = zeros(6, Nc, Nc)
    c_xi = zeros(6, Nc + 1, Nc)
    c_eta = zeros(6, Nc, Nc + 1)

    # ghost_fill_arrayop returns a ghost-extended array (not an ArrayOp)
    q_ext = ghost_fill_arrayop(phi, grid)
    Ng = grid.Ng
    @test size(q_ext) == (6, Nc + 2Ng, Nc + 2Ng)

    # ArrayOp transport utilities
    courant_xi = compute_courant_numbers(c_xi, 0.01, grid, :xi)
    courant_eta = compute_courant_numbers(c_eta, 0.01, grid, :eta)
    @test isarrayop(compute_courant_numbers_arrayop(c_xi, 0.01, grid, :xi))
    @test isarrayop(compute_courant_numbers_arrayop(c_eta, 0.01, grid, :eta))
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
