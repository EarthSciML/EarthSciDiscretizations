@testsnippet CartesianFixturesSetup begin
    using Test
    using JSON
    using EarthSciDiscretizations
    using EarthSciDiscretizations: to_esm,
        cell_centers, cell_widths, cell_volume, neighbors, metric_eval
    const grids = EarthSciDiscretizations.grids

    # Repo root: parent of `src/` two levels above this file's `test/` dir.
    const REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    const FIXTURES_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "cartesian")

    # Strip the binding/platform-specific provenance block so the comparison
    # is binding-neutral. The committed fixtures already have it stripped.
    function _esm_neutral(grid)
        d = to_esm(grid)
        delete!(d, "provenance")
        return d
    end

    # Re-emit the fixture in canonical form (`JSON.print(..., 2)` + trailing
    # newline) for byte-identity comparison with the on-disk file.
    function _canonical_bytes(d::AbstractDict)
        io = IOBuffer()
        JSON.print(io, d, 2)
        print(io, "\n")
        return String(take!(io))
    end

    # Construct a cartesian grid from a fixture's declarative payload.
    function _grid_from_fixture(d::AbstractDict)
        if haskey(d, "edges")
            edges = [Float64.(e) for e in d["edges"] if !isempty(e)]
            return grids.cartesian(; edges = edges)
        end
        n = Int.(d["n"])
        ext = [(Float64(e[1]), Float64(e[2])) for e in d["extent"]]
        ndim = Int(d["ndim"])
        ghosts = Int(get(d, "ghosts", 0))
        if ndim == 1
            return grids.cartesian(; nx = n[1], extent = ext, ghosts = ghosts)
        elseif ndim == 2
            return grids.cartesian(; nx = n[1], ny = n[2], extent = ext, ghosts = ghosts)
        elseif ndim == 3
            return grids.cartesian(;
                nx = n[1], ny = n[2], nz = n[3], extent = ext, ghosts = ghosts
            )
        else
            error("cartesian fixture has unsupported ndim=$ndim")
        end
    end
end

@testitem "cartesian fixtures: layer A byte-diff for every committed .esm" setup = [CartesianFixturesSetup] tags = [:grid, :cartesian, :fixtures] begin
    @assert isdir(FIXTURES_DIR) "missing fixtures dir at $(FIXTURES_DIR)"
    esms = sort(filter(f -> endswith(f, ".esm"), readdir(FIXTURES_DIR)))
    @test !isempty(esms)
    for fname in esms
        path = joinpath(FIXTURES_DIR, fname)
        on_disk = read(path, String)
        d_in = JSON.parsefile(path; dicttype = Dict{String, Any})

        # Round-trip: declaration → grid → declaration, then byte-diff.
        g = _grid_from_fixture(d_in)
        d_out = _esm_neutral(g)
        @test _canonical_bytes(d_out) == on_disk

        # Schema-level invariants every cartesian fixture must satisfy.
        @test d_in["family"] == "cartesian"
        @test d_in["topology"] == "rectilinear"
        @test d_in["dtype"] == "float64"
        @test d_in["ndim"] in (1, 2, 3)
        @test length(d_in["n"]) == d_in["ndim"]
        @test prod(d_in["n"]) == d_in["n_cells"]
    end
end

@testitem "cartesian fixtures: uniform_1d_n16 topology + metrics" setup = [CartesianFixturesSetup] tags = [:grid, :cartesian, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "uniform_1d_n16.esm"); dicttype = Dict{String, Any}
    ))
    @test g isa CartesianGrid{Float64, 1}
    @test g.n == (16,)
    @test g.extent == ((0.0, 1.0),)

    # Cell widths uniform at 1/16.
    @test all(isapprox.(g.widths[1], 1.0 / 16; atol = 1.0e-14))

    # First / last cell centers.
    @test cell_centers(g, 1)[1] ≈ 1 / 32
    @test cell_centers(g, 16)[1] ≈ 1 - 1 / 32

    # Neighbors at an interior cell — both sides present.
    nb = neighbors(g, 8)
    @test nb[(1, -1)] == (7,)
    @test nb[(1, +1)] == (9,)

    # Boundary cell — only one neighbor.
    nb_lo = neighbors(g, 1)
    @test !haskey(nb_lo, (1, -1))
    @test nb_lo[(1, +1)] == (2,)

    # Metric eval — cartesian Jacobian is exactly 1.
    @test metric_eval(g, :jacobian, 8) === 1.0
    @test metric_eval(g, :volume, 8) ≈ 1.0 / 16
    @test metric_eval(g, :dx, 8) ≈ 1.0 / 16
    @test cell_volume(g, 8) ≈ 1.0 / 16
end

@testitem "cartesian fixtures: uniform_2d_n64 topology + metrics" setup = [CartesianFixturesSetup] tags = [:grid, :cartesian, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "uniform_2d_n64.esm"); dicttype = Dict{String, Any}
    ))
    @test g isa CartesianGrid{Float64, 2}
    @test g.n == (64, 64)
    @test all(isapprox.(g.widths[1], 1.0 / 64; atol = 1.0e-14))
    @test all(isapprox.(g.widths[2], 1.0 / 64; atol = 1.0e-14))

    @test cell_centers(g, 1, 1) == (1 / 128, 1 / 128)
    @test cell_widths(g, 32, 32) == (1.0 / 64, 1.0 / 64)
    @test cell_volume(g, 32, 32) ≈ (1.0 / 64)^2
    @test metric_eval(g, :face_area_x, 32, 32) ≈ 1.0 / 64
    @test metric_eval(g, :face_area_y, 32, 32) ≈ 1.0 / 64

    # Identity metric tensor.
    gten = metric_eval(g, :g, 32, 32)
    @test gten == ((1.0, 0.0), (0.0, 1.0))
end

@testitem "cartesian fixtures: uniform_3d_n16 topology + metrics" setup = [CartesianFixturesSetup] tags = [:grid, :cartesian, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "uniform_3d_n16.esm"); dicttype = Dict{String, Any}
    ))
    @test g isa CartesianGrid{Float64, 3}
    @test g.n == (16, 16, 16)
    @test cell_volume(g, 8, 8, 8) ≈ (1.0 / 16)^3
    @test metric_eval(g, :face_area_z, 8, 8, 8) ≈ (1.0 / 16)^2

    # Corner cell — three missing-neighbor sides.
    nb = neighbors(g, 1, 1, 1)
    @test !haskey(nb, (1, -1))
    @test !haskey(nb, (2, -1))
    @test !haskey(nb, (3, -1))
    @test nb[(1, +1)] == (2, 1, 1)
    @test nb[(2, +1)] == (1, 2, 1)
    @test nb[(3, +1)] == (1, 1, 2)
end

@testitem "cartesian fixtures: nonuniform_2d_stretched edges + widths" setup = [CartesianFixturesSetup] tags = [:grid, :cartesian, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "nonuniform_2d_stretched.esm"); dicttype = Dict{String, Any}
    ))
    @test g isa CartesianGrid{Float64, 2}
    @test g.n == (4, 3)
    @test g.uniform == (false, false)

    # Edges round-tripped from the fixture.
    @test g.edges[1] == [0.0, 0.1, 0.3, 0.7, 1.0]
    @test g.edges[2] == [-1.0, 0.0, 0.5, 1.0]

    # Widths derived from edges.
    @test g.widths[1] ≈ [0.1, 0.2, 0.4, 0.3]
    @test g.widths[2] ≈ [1.0, 0.5, 0.5]

    # Cell volume at (i=3, j=1) = 0.4 * 1.0.
    @test cell_volume(g, 3, 1) ≈ 0.4
    # Face area normal to x at (3, 1) is the y-width.
    @test metric_eval(g, :face_area_x, 3, 1) ≈ 1.0
end
