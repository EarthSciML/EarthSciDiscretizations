@testsnippet CartesianSetup begin
    using Test
    using EarthSciDiscretizations
    const grids = EarthSciDiscretizations.grids
end

@testitem "cartesian: 1D uniform construction" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=10, extent=[(0.0, 1.0)])
    @test g isa CartesianGrid{Float64,1}
    @test ndims(g) == 1
    @test eltype(g) == Float64
    @test g.n == (10,)
    @test g.ghosts == 0
    @test g.uniform == (true,)
    @test length(g.edges[1]) == 11
    @test length(g.centers[1]) == 10
    @test g.edges[1][1] ≈ 0.0
    @test g.edges[1][end] ≈ 1.0
    # Uniform width: 0.1
    @test all(w -> isapprox(w, 0.1; atol=1e-12), g.widths[1])
end

@testitem "cartesian: 2D uniform construction (matrix extent)" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    ext = [0.0 0.0; 2.0 4.0]   # 2x2
    g = grids.cartesian(; nx=4, ny=8, extent=ext)
    @test g isa CartesianGrid{Float64,2}
    @test g.n == (4, 8)
    @test g.extent == ((0.0, 2.0), (0.0, 4.0))
    @test all(isapprox.(g.widths[1], 0.5; atol=1e-12))
    @test all(isapprox.(g.widths[2], 0.5; atol=1e-12))
end

@testitem "cartesian: 3D uniform construction" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=2, ny=3, nz=4,
                          extent=[(0.0, 1.0), (0.0, 1.5), (0.0, 2.0)])
    @test g isa CartesianGrid{Float64,3}
    @test g.n == (2, 3, 4)
    @test g.uniform == (true, true, true)
end

@testitem "cartesian: non-uniform via explicit edges" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    xe = [0.0, 0.1, 0.3, 0.7, 1.0]
    ye = [-1.0, 0.0, 0.5, 1.0]
    g = grids.cartesian(; edges=[xe, ye])
    @test g isa CartesianGrid{Float64,2}
    @test g.n == (4, 3)
    @test g.uniform[1] == false   # widths [0.1, 0.2, 0.4, 0.3]
    @test g.uniform[2] == false   # widths [1.0, 0.5, 0.5]
    @test g.extent == ((0.0, 1.0), (-1.0, 1.0))
    @test g.widths[1] ≈ [0.1, 0.2, 0.4, 0.3]
    @test g.widths[2] ≈ [1.0, 0.5, 0.5]
end

@testitem "cartesian: dtype Float32" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=4, extent=[(0.0, 1.0)], dtype=Float32)
    @test g isa CartesianGrid{Float32,1}
    @test eltype(g.centers[1]) == Float32
end

@testitem "cartesian: ghosts option" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=4, extent=[(0.0, 1.0)], ghosts=3)
    @test g.ghosts == 3
end

# --- Accessor contract: cell_centers / cell_widths / cell_volume / neighbors ---

@testitem "cartesian: cell_centers derive from extent (uniform)" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=4, extent=[(0.0, 1.0)])
    # Centers at 0.125, 0.375, 0.625, 0.875
    @test cell_centers(g, 1)[1] ≈ 0.125
    @test cell_centers(g, 2)[1] ≈ 0.375
    @test cell_centers(g, 3)[1] ≈ 0.625
    @test cell_centers(g, 4)[1] ≈ 0.875
end

@testitem "cartesian: cell_centers in 2D" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=2, ny=2, extent=[(0.0, 2.0), (10.0, 14.0)])
    @test cell_centers(g, 1, 1) == (0.5, 11.0)
    @test cell_centers(g, 2, 2) == (1.5, 13.0)
end

@testitem "cartesian: cell_volume = product of widths" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g3 = grids.cartesian(; nx=2, ny=2, nz=2,
                           extent=[(0.0, 2.0), (0.0, 4.0), (0.0, 8.0)])
    @test cell_volume(g3, 1, 1, 1) ≈ 1.0 * 2.0 * 4.0
    g_nu = grids.cartesian(; edges=[[0.0, 0.5, 1.5]])  # widths 0.5, 1.0
    @test cell_volume(g_nu, 1) ≈ 0.5
    @test cell_volume(g_nu, 2) ≈ 1.0
end

@testitem "cartesian: neighbors omit out-of-bounds sides" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=3, ny=3, extent=[(0.0, 1.0), (0.0, 1.0)])
    # Interior cell (2,2) has 4 neighbors
    nb = neighbors(g, 2, 2)
    @test length(nb) == 4
    @test nb[(1, -1)] == (1, 2)
    @test nb[(1, +1)] == (3, 2)
    @test nb[(2, -1)] == (2, 1)
    @test nb[(2, +1)] == (2, 3)
    # Corner cell (1,1) has 2 neighbors
    nb_corner = neighbors(g, 1, 1)
    @test length(nb_corner) == 2
    @test haskey(nb_corner, (1, +1))
    @test haskey(nb_corner, (2, +1))
    @test !haskey(nb_corner, (1, -1))
    @test !haskey(nb_corner, (2, -1))
end

@testitem "cartesian: metric_eval volume / face_area / dx" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=2, ny=4, extent=[(0.0, 2.0), (0.0, 4.0)])
    @test metric_eval(g, :volume, 1, 1) ≈ 1.0 * 1.0
    @test metric_eval(g, :dx, 1, 1) ≈ 1.0
    @test metric_eval(g, :dy, 1, 1) ≈ 1.0
    @test metric_eval(g, :face_area_x, 1, 1) ≈ 1.0   # = dy
    @test metric_eval(g, :face_area_y, 1, 1) ≈ 1.0   # = dx
    @test metric_eval(g, :jacobian, 1, 1) == 1.0
    gI = metric_eval(g, :g, 1, 1)
    @test gI == ((1.0, 0.0), (0.0, 1.0))
end

@testitem "cartesian: metric_eval rejects unknown axis / name" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g1 = grids.cartesian(; nx=2, extent=[(0.0, 1.0)])
    @test_throws ArgumentError metric_eval(g1, :dy, 1)
    @test_throws ArgumentError metric_eval(g1, :dz, 1)
    @test_throws ArgumentError metric_eval(g1, :bogus, 1)
end

# --- Topology consistency ---

@testitem "cartesian: topology consistency — edges bound centers" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    for nx in (1, 5, 16), ext in ((0.0, 1.0), (-3.0, 7.5))
        g = grids.cartesian(; nx=nx, extent=[ext])
        @test all(g.edges[1][i] < g.centers[1][i] < g.edges[1][i+1] for i in 1:nx)
        @test sum(g.widths[1]) ≈ ext[2] - ext[1]
        @test g.edges[1][1] ≈ ext[1]
        @test g.edges[1][end] ≈ ext[2]
    end
end

@testitem "cartesian: topology consistency — non-uniform partition of unity" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    xe = [0.0, 0.05, 0.2, 0.55, 1.0]
    g = grids.cartesian(; edges=[xe])
    @test sum(g.widths[1]) ≈ xe[end] - xe[1]
    @test g.edges[1] == xe
    # 2D partition of unity: total volume = product of axis spans
    g2 = grids.cartesian(; edges=[[0.0, 0.3, 1.0], [0.0, 0.5, 2.0]])
    total = sum(cell_volume(g2, i, j) for i in 1:g2.n[1], j in 1:g2.n[2])
    @test total ≈ 1.0 * 2.0
end

# --- Metric accuracy: uniform spacing matches analytic dx ---

@testitem "cartesian: metric accuracy — uniform dx exact" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    nx = 7
    extent = (-1.0, 1.0)
    g = grids.cartesian(; nx=nx, extent=[extent])
    expected_dx = (extent[2] - extent[1]) / nx
    for i in 1:nx
        @test metric_eval(g, :dx, i) ≈ expected_dx
    end
end

@testitem "cartesian: metric accuracy — non-uniform widths exact" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    xe = [0.0, 1.0, 1.5, 4.0, 4.25]
    g = grids.cartesian(; edges=[xe])
    for i in 1:length(xe)-1
        @test metric_eval(g, :dx, i) ≈ xe[i+1] - xe[i]
    end
end

# --- Error contract per GRIDS_API.md §9 ---

@testitem "cartesian: missing required option raises ArgumentError" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    @test_throws ArgumentError grids.cartesian()
    @test_throws ArgumentError grids.cartesian(; nx=4)                 # missing extent
    @test_throws ArgumentError grids.cartesian(; nx=4, ny=2)           # 2D missing extent
    @test_throws ArgumentError grids.cartesian(; nx=4, nz=2,
                                                 extent=[(0.0, 1.0), (0.0, 1.0)])  # ny missing for 3D
end

@testitem "cartesian: invalid extent / negative ghosts raise" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    @test_throws DomainError grids.cartesian(; nx=4, extent=[(1.0, 0.0)])   # hi <= lo
    @test_throws ArgumentError grids.cartesian(; nx=0, extent=[(0.0, 1.0)]) # n < 1
    @test_throws ArgumentError grids.cartesian(; nx=4, extent=[(0.0, 1.0)], ghosts=-1)
end

@testitem "cartesian: non-uniform edges must be strictly increasing" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    @test_throws DomainError grids.cartesian(; edges=[[0.0, 0.5, 0.5, 1.0]])   # duplicate
    @test_throws DomainError grids.cartesian(; edges=[[0.0, 0.5, 0.25, 1.0]])  # decreasing
    @test_throws ArgumentError grids.cartesian(; edges=[[0.0]])                # too short
end

@testitem "cartesian: cannot mix `edges` and `nx`/`extent`" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    @test_throws ArgumentError grids.cartesian(; nx=4, edges=[[0.0, 1.0, 2.0]])
    @test_throws ArgumentError grids.cartesian(; extent=[(0.0, 1.0)],
                                                 edges=[[0.0, 1.0, 2.0]])
end

@testitem "cartesian: non-uniform ndim cap (≤ 3)" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    @test_throws ArgumentError grids.cartesian(;
        edges=[[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
end

# --- to_esm: declarative .esm form per the mayor's correction ---

@testitem "cartesian: to_esm has §6 fields (uniform)" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=4, ny=2, extent=[(0.0, 1.0), (0.0, 0.5)], ghosts=2)
    esm = to_esm(g)
    @test esm["family"] == "cartesian"
    @test esm["topology"] == "rectilinear"
    @test esm["dtype"] == "float64"
    @test esm["ndim"] == 2
    @test esm["ghosts"] == 2
    @test esm["n_cells"] == 4 * 2
    @test esm["n"] == [4, 2]
    @test esm["uniform"] == [true, true]
    @test esm["extent"][1] == [0.0, 1.0]
    @test esm["extent"][2] == [0.0, 0.5]
    # Uniform: no explicit edges in the lowering
    @test !haskey(esm, "edges")
    # Provenance present and well-formed
    prov = esm["provenance"]
    @test prov["binding"] == "julia"
    @test haskey(prov, "binding_version")
    @test haskey(prov, "platform")
    @test haskey(prov, "math_lib")
end

@testitem "cartesian: to_esm includes edges for non-uniform" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    xe = [0.0, 0.1, 0.4, 1.0]
    g = grids.cartesian(; edges=[xe])
    esm = to_esm(g)
    @test esm["uniform"] == [false]
    @test haskey(esm, "edges")
    @test esm["edges"][1] == xe
end

@testitem "cartesian: to_esm dtype reflects element type" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g32 = grids.cartesian(; nx=2, extent=[(0.0, 1.0)], dtype=Float32)
    @test to_esm(g32)["dtype"] == "float32"
end

# --- Determinism: same inputs → same .esm payload ---

@testitem "cartesian: determinism — repeated calls produce equal payloads" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    opts = (nx=8, ny=8, extent=[(0.0, 1.0), (0.0, 1.0)], ghosts=1)
    a = to_esm(grids.cartesian(; opts...))
    b = to_esm(grids.cartesian(; opts...))
    # Provenance is fingerprint-deterministic too (same binding/platform/runtime)
    @test a == b
end

@testitem "cartesian: as_meshes errors when Meshes.jl not loaded" setup=[CartesianSetup] tags=[:grid, :cartesian] begin
    g = grids.cartesian(; nx=2, extent=[(0.0, 1.0)])
    @test_throws ArgumentError EarthSciDiscretizations.as_meshes(g)
end
