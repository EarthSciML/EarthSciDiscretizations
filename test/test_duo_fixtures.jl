@testsnippet DuoFixturesSetup begin
    using Test
    using JSON
    using EarthSciDiscretizations
    using EarthSciDiscretizations: to_esm,
        cell_centers, neighbors, metric_eval,
        n_cells, n_vertices, n_edges, total_area, family,
        DuoLoader
    const grids = EarthSciDiscretizations.grids

    # Repo root: parent of `src/` two levels above this file's `test/` dir.
    const REPO_ROOT = dirname(dirname(pathof(EarthSciDiscretizations)))
    const FIXTURES_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo")

    # Canonical ladder: level → expected (cells, vertices, edges).
    const EXPECTED = Dict(
        "icos_level0" => (level = 0, n_cells = 20,  n_vertices = 12,  n_edges = 30),
        "icos_level1" => (level = 1, n_cells = 80,  n_vertices = 42,  n_edges = 120),
        "icos_level2" => (level = 2, n_cells = 320, n_vertices = 162, n_edges = 480),
    )

    # Strip the binding-specific provenance block so the comparison is
    # binding-neutral. The committed fixtures already have it stripped.
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

    # Construct a duo grid from a fixture's declarative payload.
    function _grid_from_fixture(d::AbstractDict)
        opts = d["options"]
        ldr = opts["loader"]
        loader = DuoLoader(;
            path = String(ldr["path"]),
            reader = String(ldr["reader"]),
            check = String(ldr["check"]),
        )
        return grids.duo(;
            loader = loader,
            R = Float64(opts["R"]),
            dtype = Float64,
            ghosts = Int(d["ghosts"]),
        )
    end
end

@testitem "duo fixtures: layer A byte-diff for every committed .esm" setup = [DuoFixturesSetup] tags = [:grid, :duo, :fixtures] begin
    @assert isdir(FIXTURES_DIR) "missing fixtures dir at $(FIXTURES_DIR)"
    esms = sort(filter(f -> endswith(f, ".esm"), readdir(FIXTURES_DIR)))
    @test esms == ["icos_level0.esm", "icos_level1.esm", "icos_level2.esm"]

    for fname in esms
        path = joinpath(FIXTURES_DIR, fname)
        on_disk = read(path, String)
        d_in = JSON.parsefile(path; dicttype = Dict{String, Any})

        # Round-trip: declaration → grid → declaration, then byte-diff.
        g = _grid_from_fixture(d_in)
        d_out = _esm_neutral(g)
        @test _canonical_bytes(d_out) == on_disk

        # Schema invariants every duo fixture must satisfy.
        @test d_in["family"] == "duo"
        @test d_in["topology"] == "unstructured"
        @test d_in["dtype"] == "float64"
        @test d_in["ghosts"] == 0
        @test d_in["schema_version"] == "1.0.0"

        opts = d_in["options"]
        @test opts["R"] == 6.371e6
        @test opts["loader"]["reader"] == "builtin_icosahedral"
        @test opts["loader"]["check"] == "strict"
        @test startswith(opts["loader"]["path"], "builtin://icosahedral/")

        # Loader level matches the declared options.level (redundant but
        # guards hand-edits of one without the other).
        declared_level = Int(opts["level"])
        path_level = parse(Int, replace(opts["loader"]["path"], "builtin://icosahedral/" => ""))
        @test declared_level == path_level

        # Dimensions obey the closed-sphere identities for an icosahedral
        # subdivision at level r: Nc=20·4^r, Nv=10·4^r+2, Ne=30·4^r, and
        # the Euler identity V - E + F = 2.
        r = declared_level
        @test d_in["n_cells"]    == 20 * 4^r
        @test d_in["n_vertices"] == 10 * 4^r + 2
        @test d_in["n_edges"]    == 30 * 4^r
        @test d_in["n_vertices"] - d_in["n_edges"] + d_in["n_cells"] == 2
    end
end

@testitem "duo fixtures: no inline geometry arrays" setup = [DuoFixturesSetup] tags = [:grid, :duo, :fixtures] begin
    # Per the 2026-04-20 mayor correction on bead dsc-8qn: fixtures are
    # declarative configs, not serialized geometry blobs.
    for fname in readdir(FIXTURES_DIR)
        endswith(fname, ".esm") || continue
        d = JSON.parsefile(joinpath(FIXTURES_DIR, fname); dicttype = Dict{String, Any})
        for forbidden in ("cells", "edges", "vertices", "faces", "lon", "lat",
                          "cell_cart", "cell_neighbors", "vertex_faces", "area")
            @test !haskey(d, forbidden)
            @test !haskey(d["options"], forbidden)
        end
    end
end

@testitem "duo fixtures: icos_level0 topology + metrics" setup = [DuoFixturesSetup] tags = [:grid, :duo, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "icos_level0.esm"); dicttype = Dict{String, Any}
    ))
    @test g isa DuoGrid{Float64}
    @test family(g) == "duo"
    @test n_cells(g) == 20
    @test n_vertices(g) == 12
    @test n_edges(g) == 30

    # Closed icosahedral mesh: every face has 3 real neighbors (no zeros).
    for c in 1:n_cells(g)
        nb = neighbors(g, c)
        @test all(n > 0 && n <= n_cells(g) for n in nb)
        @test length(unique(nb)) == 3
    end

    # Neighbor reciprocity: c → n implies n → c.
    for c in 1:n_cells(g)
        for n in neighbors(g, c)
            @test c in neighbors(g, n)
        end
    end

    # Total spherical area ≈ 4πR² (L'Huilier quadrature on unit-sphere).
    R = 6.371e6
    @test isapprox(total_area(g), 4π * R^2; rtol = 1.0e-10)

    # Per-cell cartesian coords sit on the sphere of radius R.
    for c in 1:n_cells(g)
        x = metric_eval(g, :x, c)
        y = metric_eval(g, :y, c)
        z = metric_eval(g, :z, c)
        @test isapprox(sqrt(x^2 + y^2 + z^2), R; rtol = 1.0e-12)
        # lon/lat agree with :lon, :lat metric names.
        lon, lat = cell_centers(g, c)
        @test metric_eval(g, :lon, c) == lon
        @test metric_eval(g, :lat, c) == lat
    end

    # Every cell's area is strictly positive.
    for c in 1:n_cells(g)
        @test metric_eval(g, :area, c) > 0
    end

    # Unknown metric raises.
    @test_throws ArgumentError metric_eval(g, :bogus, 1)
end

@testitem "duo fixtures: icos_level2 subdivision identity" setup = [DuoFixturesSetup] tags = [:grid, :duo, :fixtures] begin
    g = _grid_from_fixture(JSON.parsefile(
        joinpath(FIXTURES_DIR, "icos_level2.esm"); dicttype = Dict{String, Any}
    ))
    @test n_cells(g) == 320
    @test n_vertices(g) == 162
    @test n_edges(g) == 480
    @test n_vertices(g) - n_edges(g) + n_cells(g) == 2

    R = 6.371e6
    @test isapprox(total_area(g), 4π * R^2; rtol = 1.0e-10)

    # Latitude range covers the sphere (poles reachable after level-1+).
    lats = cell_centers(g).lat
    @test maximum(lats) >  0.9  # close to +π/2
    @test minimum(lats) < -0.9  # close to -π/2
end
