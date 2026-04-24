@testsnippet CartesianConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON
end

@testitem "Cartesian cross-language conformance" setup = [CartesianConformanceSetup] tags = [:conformance, :grid] begin
    HARNESS_DIR = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "cartesian")
    FIXTURES = JSON.parsefile(joinpath(HARNESS_DIR, "fixtures.json"))
    REL_TOL = Float64(FIXTURES["tolerance"]["relative"])

    # Query points are 0-indexed (binding-neutral); Julia accessors are
    # 1-indexed. Convert at the test boundary so the golden remains shared.
    qp_1idx(qp) = Tuple(Int(x) + 1 for x in qp)

    function close_rel(a::Real, b::Real, tol::Real)
        scale = max(1.0, abs(a), abs(b))
        return abs(a - b) <= tol * scale
    end

    function build_grid(opts)
        kw = Dict{Symbol, Any}()
        if haskey(opts, "edges")
            kw[:edges] = [Vector{Float64}(Float64.(e)) for e in opts["edges"]]
        else
            haskey(opts, "nx") && (kw[:nx] = Int(opts["nx"]))
            haskey(opts, "ny") && (kw[:ny] = Int(opts["ny"]))
            haskey(opts, "nz") && (kw[:nz] = Int(opts["nz"]))
            if haskey(opts, "extent")
                kw[:extent] = [(Float64(e[1]), Float64(e[2])) for e in opts["extent"]]
            end
        end
        kw[:ghosts] = Int(get(opts, "ghosts", 0))
        return EarthSciDiscretizations.grids.cartesian(; kw...)
    end

    for fixture in FIXTURES["fixtures"]
        name = fixture["name"]
        golden = JSON.parsefile(joinpath(HARNESS_DIR, "golden", "$name.json"))

        @testset "fixture=$name" begin
            grid = build_grid(fixture["opts"])
            N = ndims(grid)
            @test N == golden["ndim"]
            @test prod([length(grid.widths[d]) for d in 1:N]) == golden["n_cells"]
            @test length(fixture["query_points"]) == length(golden["cell_centers"])

            for (k, qp) in enumerate(fixture["query_points"])
                idx1 = qp_1idx(qp)

                c = cell_centers(grid, idx1...)
                gc = golden["cell_centers"][k]
                for d in 1:N
                    @test close_rel(Float64(c[d]), Float64(gc[d]), REL_TOL)
                end

                w = cell_widths(grid, idx1...)
                gw = golden["cell_widths"][k]
                for d in 1:N
                    @test close_rel(Float64(w[d]), Float64(gw[d]), REL_TOL)
                end

                v = cell_volume(grid, idx1...)
                @test close_rel(Float64(v), Float64(golden["cell_volume"][k]), REL_TOL)

                # Neighbors: Julia returns a Dict keyed by (axis_1based, side);
                # golden emits a sorted list of {axis, side, index} 0-indexed.
                nbr_dict = neighbors(grid, idx1...)
                nbr_list = Vector{Dict{String, Any}}()
                for d in 1:N
                    for s in (-1, +1)
                        haskey(nbr_dict, (d, s)) || continue
                        nb = nbr_dict[(d, s)]
                        push!(nbr_list, Dict{String, Any}(
                            "axis" => d - 1,
                            "side" => s,
                            "index" => [Int(x) - 1 for x in nb],
                        ))
                    end
                end
                gnbrs = golden["neighbors"][k]
                @test length(nbr_list) == length(gnbrs)
                for (got, exp) in zip(nbr_list, gnbrs)
                    @test got["axis"] == Int(exp["axis"])
                    @test got["side"] == Int(exp["side"])
                    @test got["index"] == [Int(x) for x in exp["index"]]
                end

                # Scalar metrics
                @test close_rel(
                    Float64(metric_eval(grid, :volume, idx1...)),
                    Float64(golden["metric_volume"][k]),
                    REL_TOL,
                )
                @test close_rel(
                    Float64(metric_eval(grid, :jacobian, idx1...)),
                    Float64(golden["metric_jacobian"][k]),
                    REL_TOL,
                )
                width_names = (:dx, :dy, :dz)
                face_names = (:face_area_x, :face_area_y, :face_area_z)
                for d in 1:N
                    @test close_rel(
                        Float64(metric_eval(grid, width_names[d], idx1...)),
                        Float64(golden["metric_$(width_names[d])"][k]),
                        REL_TOL,
                    )
                    @test close_rel(
                        Float64(metric_eval(grid, face_names[d], idx1...)),
                        Float64(golden["metric_$(face_names[d])"][k]),
                        REL_TOL,
                    )
                end

                # Metric tensor (identity for cartesian)
                g = metric_eval(grid, :g, idx1...)
                gg = golden["metric_g"][k]
                for i in 1:N, j in 1:N
                    @test close_rel(Float64(g[i][j]), Float64(gg[i][j]), REL_TOL)
                end
            end
        end
    end
end
