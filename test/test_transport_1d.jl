@testsnippet Transport1DSetup begin
    using Test
    using EarthSciDiscretizations

    # LF xi tendency for interior cells (physical ci=2..Nc-1), avoids cross-panel
    # ghost access. courant[p,k,j] = Courant at interface between cells k-1 and k.
    function _lf_xi_tendency_interior(q, courant, grid)
        Nc = grid.Nc
        T = zeros(6, Nc - 2, Nc)
        for p in 1:6, k in 1:(Nc - 2), j in 1:Nc
            ci = k + 1
            cw = courant[p, ci, j]
            ce = courant[p, ci + 1, j]
            Fw = (cw + abs(cw)) / 2 * q[p, ci - 1, j] + (cw - abs(cw)) / 2 * q[p, ci, j]
            Fe = (ce + abs(ce)) / 2 * q[p, ci, j] + (ce - abs(ce)) / 2 * q[p, ci + 1, j]
            T[p, k, j] = -(Fe * grid.dx[p, ci + 1, j] - Fw * grid.dx[p, ci, j]) / grid.area[p, ci, j]
        end
        return T
    end
end

@testitem "LF cubed-sphere xi: constant field with zero velocity → zero tendency" setup = [Transport1DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    q = fill(3.0, 6, Nc, Nc)
    courant = fill(0.0, 6, Nc, Nc)
    tendency = _lf_xi_tendency_interior(q, courant, grid)
    @test size(tendency) == (6, Nc - 2, Nc)
    @test all(abs.(tendency) .< 1.0e-12)
end

@testitem "LF cubed-sphere xi: sign for uniform positive flow" setup = [Transport1DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = Float64(i)
    end
    courant = fill(0.5, 6, Nc, Nc)
    tendency = _lf_xi_tendency_interior(q, courant, grid)
    @test all(tendency .< 0.0)
end

@testitem "LF cubed-sphere xi: linearity in q" setup = [Transport1DSetup] tags = [:transport] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    q1 = randn(6, Nc, Nc)
    q2 = randn(6, Nc, Nc)
    courant = fill(0.3, 6, Nc, Nc)
    t1 = _lf_xi_tendency_interior(q1, courant, grid)
    t2 = _lf_xi_tendency_interior(q2, courant, grid)
    t_sum = _lf_xi_tendency_interior(q1 + q2, courant, grid)
    @test isapprox(t_sum, t1 + t2; rtol = 1.0e-10)
end
