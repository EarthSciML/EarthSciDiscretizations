@testsnippet Transport2DSetup begin
    using Test
    using EarthSciDiscretizations
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
