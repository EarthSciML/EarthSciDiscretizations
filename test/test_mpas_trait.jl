@testsnippet MpasTraitSetup begin
    using Test
    using EarthSciDiscretizations
    using EarthSciSerialization: AbstractGrid, AbstractUnstructuredGrid
    const grids = EarthSciDiscretizations.grids

    # Tiny synthetic mesh: 4 cells arranged as a closed quadrilateral ring
    # (each cell has 2 neighbours via the 2 sides; max_edges=2 keeps the
    # ragged adjacency well-defined for trait testing).
    function _tiny_mpas()
        n_cells = 4
        n_edges = 4
        max_edges = 2
        lon = [0.0, π / 2, π, 3π / 2]
        lat = zeros(n_cells)
        area = [1.0, 1.0, 1.0, 1.0]
        n_eoc = fill(2, n_cells)
        coc = [
            4 1 2 3;     # slot 1: previous neighbour
            2 3 4 1;     # slot 2: next neighbour
        ]
        eoc = [
            4 1 2 3;
            1 2 3 4;
        ]
        coe = [
            1 2 3 4;
            2 3 4 1;
        ]
        dc = fill(π / 2, n_edges)
        dv = fill(π / 2, n_edges)
        return EarthSciDiscretizations.mpas_mesh_data(;
            lon_cell = lon, lat_cell = lat, area_cell = area,
            n_edges_on_cell = n_eoc, cells_on_cell = coc, edges_on_cell = eoc,
            lon_edge = collect(0.0:(π / 2):(3π / 2)), lat_edge = zeros(n_edges),
            cells_on_edge = coe, dc_edge = dc, dv_edge = dv,
            max_edges = max_edges, R = 1.0,
        )
    end
end

@testitem "mpas trait: Tier-C bulk arrays" setup = [MpasTraitSetup] tags = [:grid, :mpas, :trait] begin
    g = grids.mpas(; mesh = _tiny_mpas(), R = 1.0)
    @test g isa AbstractUnstructuredGrid
    @test n_cells(g) == 4
    @test n_dims(g) == 1
    @test axis_names(g) == (:cell,)
    cl = cell_centers(g, :lon)
    cφ = cell_centers(g, :lat)
    @test length(cl) == 4 && length(cφ) == 4
    @test cell_widths(g, :cell) ≈ fill(1.0, 4)   # sqrt(area)
    @test cell_volume(g) == [1.0, 1.0, 1.0, 1.0]
end

@testitem "mpas trait: slot-based neighbour indexing" setup = [MpasTraitSetup] tags = [:grid, :mpas, :trait] begin
    g = grids.mpas(; mesh = _tiny_mpas(), R = 1.0)
    nb1 = neighbor_indices(g, :cell, 1)
    nb2 = neighbor_indices(g, :cell, 2)
    @test nb1 == [4, 1, 2, 3]
    @test nb2 == [2, 3, 4, 1]
    @test_throws ArgumentError neighbor_indices(g, :cell, 0)
    @test_throws ArgumentError neighbor_indices(g, :cell, 3)
    @test_throws ArgumentError neighbor_indices(g, :lon, 1)
end

@testitem "mpas trait: Tier-U adjacency accessors" setup = [MpasTraitSetup] tags = [:grid, :mpas, :trait] begin
    g = grids.mpas(; mesh = _tiny_mpas(), R = 1.0)
    @test cell_neighbor_table(g) == g.mesh.cells_on_cell
    @test cell_valence(g) == g.mesh.n_edges_on_cell
    @test edge_length(g) == g.mesh.dv_edge
    @test cell_distance(g) == g.mesh.dc_edge
    # Closed ring → no boundary cells.
    @test all(.!boundary_mask(g, :cell, :lower))
end

@testitem "mpas trait: cell_centers requires coord axis" setup = [MpasTraitSetup] tags = [:grid, :mpas, :trait] begin
    g = grids.mpas(; mesh = _tiny_mpas(), R = 1.0)
    @test_throws ArgumentError cell_centers(g, :cell)
end
