@testsnippet CubedSphereFixturesSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON

    # Inline topology / metric validation of the committed canonical .esm
    # grid fixtures at discretizations/grids/cubed_sphere/. Consumed by the
    # ESD CI walker (dsc-usq) as Layer A byte-diff-style sanity plus downstream
    # numeric invariants. The committed .esm is the declarative config only;
    # geometry is derived at runtime via the Julia reference binding here.
    const FIXTURES_DIR = joinpath(
        @__DIR__, "..", "discretizations", "grids", "cubed_sphere",
    )
    const FIXTURE_NAMES = ("c24", "c48", "c96")
end

# ---------- Layer A: declarative schema byte-shape ----------

@testitem "cubed_sphere fixtures: schema shape (declarative, no inline geometry)" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    for name in FIXTURE_NAMES
        path = joinpath(FIXTURES_DIR, "$name.esm")
        @test isfile(path)
        doc = JSON.parsefile(path)

        @test doc["family"] == "cubed_sphere"
        @test doc["version"] == "1.0.0"
        @test doc["dtype"] == "float64"
        @test doc["topology"] == "block_structured"
        @test doc["generator"] == "gnomonic_c6"

        params = doc["params"]
        @test haskey(params, "Nc")
        @test haskey(params, "R")
        @test haskey(params, "ghosts")
        @test params["Nc"] isa Integer
        @test params["Nc"] >= 1
        @test params["R"] isa Real && params["R"] > 0
        @test params["ghosts"] isa Integer && params["ghosts"] >= 0

        # Declarative form: no serialized geometry arrays. Per mayor
        # correction 2026-04-20, §6 schema forbids inline cells/edges/vertices.
        for forbidden in ("lon", "lat", "area", "cells", "vertices", "edges")
            @test !haskey(doc, forbidden)
        end

        # Cross-binding byte identity: canonical committed fixture omits
        # provenance (per-build). Binding runtime to_esm adds it back.
        @test !haskey(doc, "provenance")
    end
end

@testitem "cubed_sphere fixtures: resolutions span canonical corpus (Nc = 24/48/96)" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    expected = Dict("c24" => 24, "c48" => 48, "c96" => 96)
    for (name, Nc_expected) in expected
        doc = JSON.parsefile(joinpath(FIXTURES_DIR, "$name.esm"))
        @test Int(doc["params"]["Nc"]) == Nc_expected
        @test Float64(doc["params"]["R"]) == 6371000.0
        @test Int(doc["params"]["ghosts"]) == 0
    end
end

# ---------- Topology invariants derived from the declared config ----------

@testitem "cubed_sphere fixtures: topology invariants (6 panels, Nc×Nc per panel)" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    for name in FIXTURE_NAMES
        doc = JSON.parsefile(joinpath(FIXTURES_DIR, "$name.esm"))
        Nc = Int(doc["params"]["Nc"])
        R = Float64(doc["params"]["R"])
        ghosts = Int(doc["params"]["ghosts"])

        grid = CubedSphereGrid(Nc; R = R, Ng = ghosts)
        @test grid.Nc == Nc
        @test grid.R == R
        @test grid.Ng == ghosts

        # Cell-center arrays are shape (6, Nc, Nc).
        @test size(grid.lon) == (6, Nc, Nc)
        @test size(grid.lat) == (6, Nc, Nc)
        @test size(grid.area) == (6, Nc, Nc)

        # Every panel has Nc×Nc cells; total = 6·Nc².
        @test 6 * Nc * Nc == length(grid.lon)
    end
end

@testitem "cubed_sphere fixtures: panel connectivity involutive at declared Nc" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    for name in FIXTURE_NAMES
        doc = JSON.parsefile(joinpath(FIXTURES_DIR, "$name.esm"))
        Nc = Int(doc["params"]["Nc"])

        # PANEL_CONNECTIVITY is Nc-independent but the test exercises that
        # each declared fixture can actually materialize a grid that obeys
        # the adjacency involutive property: neighbor of neighbor in the
        # opposite direction returns to the origin panel.
        grid = CubedSphereGrid(Nc; R = 1.0, Ng = 0)
        @test grid.Nc == Nc

        for p in 1:6
            for (dir, opp) in ((West, East), (East, West), (South, North), (North, South))
                nb = PANEL_CONNECTIVITY[p][dir]
                nb_back = PANEL_CONNECTIVITY[nb.neighbor_panel][nb.neighbor_edge]
                @test nb_back.neighbor_panel == p
            end
        end
    end
end

# ---------- Metric invariants at representative cells ----------

@testitem "cubed_sphere fixtures: metric positive-definite at cell centers" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    for name in FIXTURE_NAMES
        doc = JSON.parsefile(joinpath(FIXTURES_DIR, "$name.esm"))
        Nc = Int(doc["params"]["Nc"])
        R = Float64(doc["params"]["R"])

        grid = CubedSphereGrid(Nc; R = R, Ng = 0)

        # Sample corner, edge, and interior cells on each panel.
        samples = [
            (1, 1), (Nc, Nc), (1, Nc), (Nc, 1),
            (max(1, Nc ÷ 2), max(1, Nc ÷ 2)),
        ]
        for p in 1:6, (i, j) in samples
            J, gxx, gee, gxe = gnomonic_metric(
                grid.ξ_centers[i], grid.η_centers[j], R,
            )
            @test J > 0
            @test gxx > 0
            @test gee > 0
            # det(g) = gxx·gηη - gξη² > 0 for a Riemannian metric.
            @test gxx * gee - gxe^2 > 0
            # All cell areas are strictly positive.
            @test grid.area[p, i, j] > 0
        end
    end
end

@testitem "cubed_sphere fixtures: total surface area → 4πR² within tolerance" setup =
    [CubedSphereFixturesSetup] tags = [:grid, :cubed_sphere, :fixtures] begin
    # Cubed-sphere tessellation is exact: the sum of all 6·Nc² spherical-quad
    # cell areas must equal the closed-form sphere surface area 4πR². This
    # is the strongest scalar integrity check on the declared config.
    for name in FIXTURE_NAMES
        doc = JSON.parsefile(joinpath(FIXTURES_DIR, "$name.esm"))
        Nc = Int(doc["params"]["Nc"])
        R = Float64(doc["params"]["R"])

        grid = CubedSphereGrid(Nc; R = R, Ng = 0)
        @test isapprox(total_area(grid), 4π * R^2; rtol = 1.0e-10)
    end
end
