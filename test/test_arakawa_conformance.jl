@testsnippet ArakawaConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON
end

@testitem "Arakawa cross-language conformance" setup = [ArakawaConformanceSetup] tags = [:conformance, :arakawa, :grid] begin
    HARNESS_DIR = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "arakawa")
    FIXTURES = JSON.parsefile(joinpath(HARNESS_DIR, "fixtures.json"))
    REL_TOL = Float64(FIXTURES["tolerance"]["relative"])

    function close_rel(a::Real, b::Real, tol::Real)
        scale = max(1.0, abs(a), abs(b))
        return abs(a - b) <= tol * scale
    end

    stagger_sym(s::AbstractString) = Symbol(s)
    loc_sym(s::AbstractString) = s == "cell_center" ? CellCenter :
        s == "u_edge" ? UEdge :
        s == "v_edge" ? VEdge :
        s == "corner" ? Corner :
        error("unknown location $s")

    function coord_0idx(grid_a::ArakawaGrid, grid_c::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
        if loc === CellCenter
            return cell_centers(grid_a, i + 1, j + 1)
        elseif loc === Corner
            return corners(grid_a, i + 1, j + 1)
        elseif loc === UEdge
            return u_face(grid_c, i + 1, j + 1)
        else
            return v_face(grid_c, i + 1, j + 1)
        end
    end

    function neighbors_wire(grid::ArakawaGrid, loc::VarLocation, i::Int, j::Int)
        w, e, s, n = neighbors(grid, loc, i + 1, j + 1)
        wire(t) = t === nothing ? nothing : [t[2] - 1, t[3] - 1]
        return Dict("W" => wire(w), "E" => wire(e), "S" => wire(s), "N" => wire(n))
    end

    function check_neighbor(got, expected, ctx::AbstractString)
        if expected === nothing
            @test got === nothing
        else
            ex = [Int(expected[1]), Int(expected[2])]
            @test got == ex
        end
    end

    for fixture in FIXTURES["fixtures"]
        name = fixture["name"]
        opts = fixture["opts"]
        base_opts = opts["base"]
        base = CartesianBase(
            xlo = float(base_opts["xlo"]),
            xhi = float(base_opts["xhi"]),
            ylo = float(base_opts["ylo"]),
            yhi = float(base_opts["yhi"]),
            nx = Int(base_opts["nx"]),
            ny = Int(base_opts["ny"]),
        )
        dtype = opts["dtype"] == "float32" ? Float32 : Float64
        ghosts = Int(opts["ghosts"])

        golden = JSON.parsefile(joinpath(HARNESS_DIR, "golden", "$name.json"))

        @testset "fixture=$name" begin
            @test base.nx * base.ny == golden["n_cells"]

            g_a = EarthSciDiscretizations.grids.arakawa(
                base = base, stagger = :A, ghosts = ghosts, dtype = dtype,
            )
            g_c = EarthSciDiscretizations.grids.arakawa(
                base = base, stagger = :C, ghosts = ghosts, dtype = dtype,
            )

            # --- coord + neighbor tables (stagger-independent) ---
            for lname in ("cell_center", "u_edge", "v_edge", "corner")
                loc = loc_sym(lname)
                ctable = golden["coords"][lname]
                ntable = golden["neighbors"][lname]
                for (k, qp) in enumerate(ctable["points"])
                    i, j = Int(qp[1]), Int(qp[2])
                    x, y = coord_0idx(g_a, g_c, loc, i, j)
                    expected_xy = ctable["xy"][k]
                    @test close_rel(x, Float64(expected_xy[1]), REL_TOL)
                    @test close_rel(y, Float64(expected_xy[2]), REL_TOL)

                    got = neighbors_wire(g_a, loc, i, j)
                    for dkey in ("W", "E", "S", "N")
                        check_neighbor(got[dkey], ntable[dkey][k],
                            "loc=$lname dir=$dkey k=$k")
                    end
                end
            end

            # --- metrics (stagger-independent for cartesian) ---
            mtable = golden["metrics"]
            for (k, qp) in enumerate(mtable["points"])
                i, j = Int(qp[1]) + 1, Int(qp[2]) + 1
                @test close_rel(metric_eval(g_a, :dx, i, j),
                    Float64(mtable["dx"][k]), REL_TOL)
                @test close_rel(metric_eval(g_a, :dy, i, j),
                    Float64(mtable["dy"][k]), REL_TOL)
                @test close_rel(metric_eval(g_a, :area, i, j),
                    Float64(mtable["area"][k]), REL_TOL)
            end

            # --- per-stagger variable location + shape tables ---
            for sname in fixture["staggers"]
                s = stagger_sym(sname)
                g_s = EarthSciDiscretizations.grids.arakawa(
                    base = base, stagger = s, ghosts = ghosts, dtype = dtype,
                )
                stab = golden["staggers"][sname]
                h_loc, u_loc, v_loc = arakawa_variable_locations(
                    s === :A ? ArakawaA :
                    s === :B ? ArakawaB :
                    s === :C ? ArakawaC :
                    s === :D ? ArakawaD : ArakawaE,
                )
                @test (g_s.stagger === ArakawaE) == stab["rotated"]
                loc_str(l) = l === CellCenter ? "cell_center" :
                    l === UEdge ? "u_edge" :
                    l === VEdge ? "v_edge" : "corner"
                @test loc_str(h_loc) == stab["variable_locations"]["h"]
                @test loc_str(u_loc) == stab["variable_locations"]["u"]
                @test loc_str(v_loc) == stab["variable_locations"]["v"]

                for var in (:h, :u, :v)
                    sh = variable_shape(g_s, var)
                    ex = stab["variable_shapes"][String(var)]
                    @test collect(sh) == [Int(ex[1]), Int(ex[2])]
                end
                for (lname, loc) in (
                    ("cell_center", CellCenter),
                    ("u_edge", UEdge),
                    ("v_edge", VEdge),
                    ("corner", Corner),
                )
                    sh = arakawa_shape(g_s, loc)
                    ex = stab["location_shapes"][lname]
                    @test collect(sh) == [Int(ex[1]), Int(ex[2])]
                end
            end
        end
    end
end
