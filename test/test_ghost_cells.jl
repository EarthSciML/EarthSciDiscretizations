@testsnippet GhostSetup begin
    using Test
    using EarthSciDiscretizations
end

@testitem "Ghost cell extension" setup=[GhostSetup] tags=[:ghost] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    Ng = grid.Ng

    # Create a field with unique values: q = 100*p + 10*i + j
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = 100.0 * p + 10.0 * i + j
    end

    q_ext = extend_with_ghosts(q, grid)

    # Check extended array size
    @test size(q_ext) == (6, Nc + 2Ng, Nc + 2Ng)

    # Interior should be preserved exactly
    for p in 1:6, i in 1:Nc, j in 1:Nc
        @test q_ext[p, i + Ng, j + Ng] == q[p, i, j]
    end

    # Ghost cells should be non-zero (filled from neighbor panels)
    for p in 1:6
        # West ghosts
        for g in 1:Ng, j in 1:Nc
            @test q_ext[p, g, j + Ng] != 0.0
        end
        # East ghosts
        for g in 1:Ng, j in 1:Nc
            @test q_ext[p, Nc + Ng + g, j + Ng] != 0.0
        end
        # South ghosts
        for i in 1:Nc, g in 1:Ng
            @test q_ext[p, i + Ng, g] != 0.0
        end
        # North ghosts
        for i in 1:Nc, g in 1:Ng
            @test q_ext[p, i + Ng, Nc + Ng + g] != 0.0
        end
    end
end

@testitem "Ghost cells with smooth field" setup=[GhostSetup] tags=[:ghost] begin
    Nc = 16
    grid = CubedSphereGrid(Nc)
    Ng = grid.Ng

    # Create q = cos(lat), which is smooth and bounded in [-1, 1]
    q = zeros(6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        q[p, i, j] = cos(grid.lat[p, i, j])
    end

    q_ext = extend_with_ghosts(q, grid)

    # All ghost values should also be in [-1, 1] since cos(lat) is bounded
    for p in 1:6
        for g in 1:Ng, j in 1:Nc
            @test -1.0 <= q_ext[p, g, j + Ng] <= 1.0
            @test -1.0 <= q_ext[p, Nc + Ng + g, j + Ng] <= 1.0
        end
        for i in 1:Nc, g in 1:Ng
            @test -1.0 <= q_ext[p, i + Ng, g] <= 1.0
            @test -1.0 <= q_ext[p, i + Ng, Nc + Ng + g] <= 1.0
        end
    end
end

@testitem "Stagger-aware grid sizes" setup=[GhostSetup] tags=[:ghost] begin
    Nc = 8; Ng = 3
    @test grid_size(CellCenter, Nc) == (Nc, Nc)
    @test grid_size(UEdge, Nc) == (Nc + 1, Nc)
    @test grid_size(VEdge, Nc) == (Nc, Nc + 1)
    @test grid_size(Corner, Nc) == (Nc + 1, Nc + 1)

    @test full_array_size(CellCenter, Nc) == (6, Nc, Nc)
    @test full_array_size(UEdge, Nc) == (6, Nc + 1, Nc)
    @test full_array_size(VEdge, Nc) == (6, Nc, Nc + 1)

    @test ghost_array_size(CellCenter, Nc, Ng) == (6, Nc + 2Ng, Nc + 2Ng)
    @test ghost_array_size(UEdge, Nc, Ng) == (6, Nc + 1 + 2Ng, Nc + 2Ng)
end

@testitem "Discrete space creation" setup=[GhostSetup] tags=[:ghost] begin
    Nc = 8
    grid = CubedSphereGrid(Nc)
    space = CubedSphereDiscreteSpace(grid, [
        :q => CellCenter,
        :u => UEdge,
        :v => VEdge,
    ])

    @test length(space.var_names) == 3

    q_arr = allocate_variable(space, :q)
    u_arr = allocate_variable(space, :u)
    v_arr = allocate_variable(space, :v)

    @test size(q_arr) == (6, Nc, Nc)
    @test size(u_arr) == (6, Nc + 1, Nc)
    @test size(v_arr) == (6, Nc, Nc + 1)
end
