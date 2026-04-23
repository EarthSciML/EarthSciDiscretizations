@testsnippet CubedSphereConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON
end

# Mirrors the 0-indexed panel-connectivity logic committed to
# tests/conformance/grids/cubed_sphere/regenerate_golden.jl. Kept inline here
# because Julia's CubedSphereGrid does not yet expose a first-class
# `neighbors(p, i, j)` accessor; the cross-language conformance contract is
# verified by computing neighbors from the same `PANEL_CONNECTIVITY` table.
@testitem "Cubed-sphere cross-language conformance" setup = [CubedSphereConformanceSetup] tags = [:conformance, :grid] begin
    HARNESS_DIR = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "cubed_sphere")
    FIXTURES = JSON.parsefile(joinpath(HARNESS_DIR, "fixtures.json"))
    REL_TOL = Float64(FIXTURES["tolerance"]["relative"])

    DIRECTIONS = (West, East, South, North)
    DIR_KEYS = ("W", "E", "S", "N")

    function neighbor_0idx(p::Int, i::Int, j::Int, dir, Nc::Int)
        if dir == West && i > 0
            return (p, i - 1, j)
        elseif dir == East && i < Nc - 1
            return (p, i + 1, j)
        elseif dir == South && j > 0
            return (p, i, j - 1)
        elseif dir == North && j < Nc - 1
            return (p, i, j + 1)
        end
        nb = PANEL_CONNECTIVITY[p + 1][dir]
        nb_panel = nb.neighbor_panel - 1
        nb_edge = nb.neighbor_edge
        along = (dir == West || dir == East) ? j : i
        if nb.reverse_index
            along = Nc - 1 - along
        end
        if nb_edge == West
            return (nb_panel, 0, along)
        elseif nb_edge == East
            return (nb_panel, Nc - 1, along)
        elseif nb_edge == South
            return (nb_panel, along, 0)
        else
            return (nb_panel, along, Nc - 1)
        end
    end

    function close_rel(a::Real, b::Real, tol::Real)
        scale = max(1.0, abs(a), abs(b))
        return abs(a - b) <= tol * scale
    end

    for fixture in FIXTURES["fixtures"]
        name = fixture["name"]
        opts = fixture["opts"]
        Nc = Int(opts["Nc"])
        R = float(opts["R"])
        ghosts = Int(opts["ghosts"])
        golden = JSON.parsefile(joinpath(HARNESS_DIR, "golden", "$name.json"))

        @testset "fixture=$name" begin
            grid = CubedSphereGrid(Nc; R = R, Ng = ghosts)
            @test 6 * Nc * Nc == golden["n_cells"]

            for (k, qp) in enumerate(fixture["query_points"])
                p, i, j = Int(qp[1]), Int(qp[2]), Int(qp[3])

                lon, lat = gnomonic_to_lonlat(grid.ξ_centers[i + 1], grid.η_centers[j + 1], p + 1)
                gc = golden["cell_centers"][k]
                @test close_rel(lon, gc["lon"], REL_TOL)
                @test close_rel(lat, gc["lat"], REL_TOL)

                for (dkey, dir) in zip(DIR_KEYS, DIRECTIONS)
                    np, ni, nj = neighbor_0idx(p, i, j, dir, Nc)
                    gn = golden["neighbors_$dkey"][k]
                    @test [np, ni, nj] == [Int(gn[1]), Int(gn[2]), Int(gn[3])]
                end

                J, gxx, gyy, gxy = gnomonic_metric(
                    grid.ξ_centers[i + 1], grid.η_centers[j + 1], R,
                )
                det = gxx * gyy - gxy * gxy
                @test close_rel(J, golden["metric_J"][k], REL_TOL)
                @test close_rel(gxx, golden["metric_g_xixi"][k], REL_TOL)
                @test close_rel(gyy, golden["metric_g_etaeta"][k], REL_TOL)
                @test close_rel(gxy, golden["metric_g_xieta"][k], REL_TOL)
                @test close_rel(gyy / det, golden["metric_ginv_xixi"][k], REL_TOL)
                @test close_rel(gxx / det, golden["metric_ginv_etaeta"][k], REL_TOL)
                @test close_rel(-gxy / det, golden["metric_ginv_xieta"][k], REL_TOL)

                area = compute_cell_area(
                    (grid.ξ_edges[i + 1], grid.ξ_edges[i + 2]),
                    (grid.η_edges[j + 1], grid.η_edges[j + 2]),
                    R, p + 1,
                )
                @test close_rel(area, golden["area"][k], REL_TOL)
            end
        end
    end
end
