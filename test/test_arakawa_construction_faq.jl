@testsnippet ArakawaConstructionFAQSetup begin
    using Test
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON
    import EarthSciSerialization

    const _ARK_STAGGER_SYMS = (:A, :B, :C, :D, :E)

    # Build a CartesianBase from a fixture's opts["base"] dict — identical to
    # tests/conformance/grids/arakawa/construction/regenerate_golden.jl.
    function _ark_build_base(base_opts::AbstractDict)
        base_opts["family"] == "cartesian" || error("only cartesian base supported today")
        return ESD.CartesianBase(
            xlo = Float64(base_opts["xlo"]), xhi = Float64(base_opts["xhi"]),
            ylo = Float64(base_opts["ylo"]), yhi = Float64(base_opts["yhi"]),
            nx = Int(base_opts["nx"]), ny = Int(base_opts["ny"]),
        )
    end

    _ark_dtype(opts) = get(opts, "dtype", "float64") == "float32" ? Float32 : Float64
    _ark_stagger_sym(s::AbstractString) = Symbol(s)

    function _ark_loc_name(loc::VarLocation)
        return loc === CellCenter ? "cell_center" :
            loc === UEdge ? "u_edge" :
            loc === VEdge ? "v_edge" : "corner"
    end

    # Bit-exact equality on floats / float arrays (the FAQ output is bit-identical to
    # the imperative builder, not merely close — this is "match imperative to ULP").
    _ark_biteq(a::Real, b::Real) =
        reinterpret(UInt64, Float64(a)) == reinterpret(UInt64, Float64(b))
    _ark_biteq(a::AbstractArray, b::AbstractArray) =
        size(a) == size(b) && all(_ark_biteq(a[i], b[i]) for i in eachindex(a))

    # Build the imperative reference (xs, ys) row-major over a location's shape from a
    # base-grid coordinate primitive (arakawa_cell_center / x_edge / y_edge / corner).
    function _ark_ref_location(prim, base, ni::Int, nj::Int)
        xs = Vector{Float64}(undef, ni * nj)
        ys = Vector{Float64}(undef, ni * nj)
        for j in 1:nj, i in 1:ni
            x, y = prim(base, i, j)
            xs[(j - 1) * ni + i] = Float64(x)
            ys[(j - 1) * ni + i] = Float64(y)
        end
        return xs, ys
    end

    # --- binding-neutral serialization — MUST match
    # tests/conformance/grids/arakawa/construction/regenerate_golden.jl byte-for-byte.
    # Query-point sampled (the shared fixtures' grids are large): coordinate floats are
    # compact JSON strings; neighbor ids 0-based with -1 sentinel; masks 0/1. ---
    _ark_compact(x) = JSON.json(x)
    # location coords at the 0-based query points -> compact [[x,y],...]
    function _ark_coord_str(xs, ys, ni::Int, points)
        return _ark_compact(
            [
                [
                        Float64(xs[Int(p[2]) * ni + Int(p[1]) + 1]),
                        Float64(ys[Int(p[2]) * ni + Int(p[1]) + 1]),
                    ] for p in points
            ]
        )
    end
    # Tier-C neighbor ids at the 0-based cell-center query points -> compact [..], -1 sentinel
    function _ark_nbr_str(arr, nx::Int, points)
        return _ark_compact(
            [
                (id = arr[Int(p[2]) * nx + Int(p[1]) + 1]; id == 0 ? -1 : Int(id) - 1)
                    for p in points
            ]
        )
    end
    # Tier-C boundary mask at the 0-based cell-center query points -> compact [0/1..]
    function _ark_mask_str(arr, nx::Int, points)
        return _ark_compact([arr[Int(p[2]) * nx + Int(p[1]) + 1] ? 1 : 0 for p in points])
    end
end

@testitem "Arakawa construction FAQ — matches imperative arakawa.jl (ULP), all A–E staggers" setup = [ArakawaConstructionFAQSetup] tags = [:conformance, :grid, :arakawa, :construction, :faq] begin
    base_dir = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "arakawa")
    SPEC = JSON.parsefile(joinpath(base_dir, "fixtures.json"))

    for f in SPEC["fixtures"]
        opts = f["opts"]
        base = _ark_build_base(opts["base"])
        dtype = _ark_dtype(opts)
        ghosts = Int(get(opts, "ghosts", 0))
        nx, ny = base.nx, base.ny

        for sname in f["staggers"]
            g = ESD.grids.arakawa(
                base = base, stagger = _ark_stagger_sym(sname), ghosts = ghosts, dtype = dtype,
            )
            faq = arakawa_construction_faq(g)

            @test faq.nx == nx
            @test faq.ny == ny

            # --- "match imperative arakawa.jl to ULP": the four staggered location
            # coordinate arrays equal the imperative base-grid primitives bit-for-bit. ---
            ref_cc = _ark_ref_location(ESD.arakawa_cell_center, base, nx, ny)
            ref_ue = _ark_ref_location(ESD.arakawa_x_edge, base, nx + 1, ny)
            ref_ve = _ark_ref_location(ESD.arakawa_y_edge, base, nx, ny + 1)
            ref_co = _ark_ref_location(ESD.arakawa_corner, base, nx + 1, ny + 1)
            @test _ark_biteq(faq.coords.cell_center[1], ref_cc[1])
            @test _ark_biteq(faq.coords.cell_center[2], ref_cc[2])
            @test _ark_biteq(faq.coords.u_edge[1], ref_ue[1])
            @test _ark_biteq(faq.coords.u_edge[2], ref_ue[2])
            @test _ark_biteq(faq.coords.v_edge[1], ref_ve[1])
            @test _ark_biteq(faq.coords.v_edge[2], ref_ve[2])
            @test _ark_biteq(faq.coords.corner[1], ref_co[1])
            @test _ark_biteq(faq.coords.corner[2], ref_co[2])

            # --- Tier-C cell-centred arrays equal the imperative trait arrays (ULP). ---
            @test _ark_biteq(faq.coords.cell_center[1], ESD.cell_centers(g, :x))
            @test _ark_biteq(faq.coords.cell_center[2], ESD.cell_centers(g, :y))
            @test _ark_biteq(faq.cell_widths[1], ESD.cell_widths(g, :x))
            @test _ark_biteq(faq.cell_widths[2], ESD.cell_widths(g, :y))
            @test _ark_biteq(faq.cell_volume, ESD.cell_volume(g))
            @test faq.neighbor_minus[1] == ESD.neighbor_indices(g, :x, -1)
            @test faq.neighbor_plus[1] == ESD.neighbor_indices(g, :x, +1)
            @test faq.neighbor_minus[2] == ESD.neighbor_indices(g, :y, -1)
            @test faq.neighbor_plus[2] == ESD.neighbor_indices(g, :y, +1)
            @test faq.boundary_lower[1] == ESD.boundary_mask(g, :x, :lower)
            @test faq.boundary_upper[1] == ESD.boundary_mask(g, :x, :upper)
            @test faq.boundary_lower[2] == ESD.boundary_mask(g, :y, :lower)
            @test faq.boundary_upper[2] == ESD.boundary_mask(g, :y, :upper)

            # --- metric (dx/dy/area) equals the imperative metric_eval. ---
            @test _ark_biteq(faq.dx, ESD.metric_eval(g, :dx, 1, 1))
            @test _ark_biteq(faq.dy, ESD.metric_eval(g, :dy, 1, 1))
            @test _ark_biteq(faq.cell_volume[1], ESD.metric_eval(g, :area, 1, 1))

            # --- static A–E variable-location / shape table matches the imperative
            # arakawa_variable_locations / variable_shape / arakawa_shape. ---
            h_loc, u_loc, v_loc = ESD.arakawa_variable_locations(g.stagger)
            @test faq.variable_locations.h === h_loc
            @test faq.variable_locations.u === u_loc
            @test faq.variable_locations.v === v_loc
            @test faq.rotated == (g.stagger === ArakawaE)
            @test faq.variable_shapes.h == ESD.variable_shape(g, :h)
            @test faq.variable_shapes.u == ESD.variable_shape(g, :u)
            @test faq.variable_shapes.v == ESD.variable_shape(g, :v)
            @test faq.location_shapes.cell_center == ESD.arakawa_shape(g, CellCenter)
            @test faq.location_shapes.u_edge == ESD.arakawa_shape(g, UEdge)
            @test faq.location_shapes.v_edge == ESD.arakawa_shape(g, VEdge)
            @test faq.location_shapes.corner == ESD.arakawa_shape(g, Corner)

            # --- internal consistency of the FAQ construction itself. ---
            @test length(faq.center_x) == nx && length(faq.center_y) == ny
            @test length(faq.node_x) == nx + 1 && length(faq.node_y) == ny + 1
            @test issorted(faq.node_x) && issorted(faq.node_y)             # nodes ascend
            @test all(faq.node_x[i] < faq.center_x[i] < faq.node_x[i + 1] for i in 1:nx)
            @test all(faq.node_y[j] < faq.center_y[j] < faq.node_y[j + 1] for j in 1:ny)
            @test all(>(0), faq.cell_widths[1]) && all(>(0), faq.cell_widths[2])
            @test all(faq.cell_volume .== faq.cell_widths[1] .* faq.cell_widths[2])

            # neighbor 0-sentinel + symmetry over the row-major (nx, ny) cell set.
            for j in 1:ny, i in 1:nx
                k = (j - 1) * nx + i
                @test (faq.neighbor_minus[1][k] == 0) == (i == 1)
                @test (faq.neighbor_plus[1][k] == 0) == (i == nx)
                @test (faq.neighbor_minus[2][k] == 0) == (j == 1)
                @test (faq.neighbor_plus[2][k] == 0) == (j == ny)
                npx = faq.neighbor_plus[1][k]
                npx != 0 && @test faq.neighbor_minus[1][npx] == k
                npy = faq.neighbor_plus[2][k]
                npy != 0 && @test faq.neighbor_minus[2][npy] == k
            end
        end
    end
end

@testitem "Arakawa construction FAQ — byte-identity vs golden" setup = [ArakawaConstructionFAQSetup] tags = [:conformance, :grid, :arakawa, :construction, :faq] begin
    base_dir = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "arakawa")
    SPEC = JSON.parsefile(joinpath(base_dir, "fixtures.json"))
    GOLDEN = JSON.parsefile(joinpath(base_dir, "construction", "golden.json"))
    by_name = Dict(fx["name"] => fx for fx in GOLDEN["fixtures"])

    for f in SPEC["fixtures"]
        opts = f["opts"]
        base = _ark_build_base(opts["base"])
        dtype = _ark_dtype(opts)
        ghosts = Int(get(opts, "ghosts", 0))
        nx, ny = base.nx, base.ny
        qp = f["query_points"]
        gl = by_name[f["name"]]

        # Geometry is stagger-independent; build once at A for coords/neighbors/boundary.
        g0 = ESD.grids.arakawa(base = base, stagger = :A, ghosts = ghosts, dtype = dtype)
        faq0 = arakawa_construction_faq(g0)

        @test gl["nx"] == nx && gl["ny"] == ny && gl["n_cells"] == nx * ny
        @test _ark_biteq(faq0.dx, gl["dx"]) && _ark_biteq(faq0.dy, gl["dy"])
        @test _ark_biteq(faq0.cell_volume[1], gl["cell_volume"])

        # staggered location coords at their query points (compact float strings).
        @test _ark_coord_str(faq0.coords.cell_center..., nx, qp["cell_center"]) ==
            gl["coords"]["cell_center"]
        @test _ark_coord_str(faq0.coords.u_edge..., nx + 1, qp["u_edge"]) ==
            gl["coords"]["u_edge"]
        @test _ark_coord_str(faq0.coords.v_edge..., nx, qp["v_edge"]) ==
            gl["coords"]["v_edge"]
        @test _ark_coord_str(faq0.coords.corner..., nx + 1, qp["corner"]) ==
            gl["coords"]["corner"]

        # Tier-C neighbor ids + boundary masks at the cell-center query points.
        ccp = qp["cell_center"]
        @test _ark_nbr_str(faq0.neighbor_minus[1], nx, ccp) == gl["neighbor_indices"]["x_minus"]
        @test _ark_nbr_str(faq0.neighbor_plus[1], nx, ccp) == gl["neighbor_indices"]["x_plus"]
        @test _ark_nbr_str(faq0.neighbor_minus[2], nx, ccp) == gl["neighbor_indices"]["y_minus"]
        @test _ark_nbr_str(faq0.neighbor_plus[2], nx, ccp) == gl["neighbor_indices"]["y_plus"]
        @test _ark_mask_str(faq0.boundary_lower[1], nx, ccp) == gl["boundary"]["x_lower"]
        @test _ark_mask_str(faq0.boundary_upper[1], nx, ccp) == gl["boundary"]["x_upper"]
        @test _ark_mask_str(faq0.boundary_lower[2], nx, ccp) == gl["boundary"]["y_lower"]
        @test _ark_mask_str(faq0.boundary_upper[2], nx, ccp) == gl["boundary"]["y_upper"]

        # static A–E table per stagger.
        for sname in f["staggers"]
            g = ESD.grids.arakawa(
                base = base, stagger = _ark_stagger_sym(sname), ghosts = ghosts, dtype = dtype,
            )
            faq = arakawa_construction_faq(g)
            gs = gl["staggers"][sname]
            @test _ark_loc_name(faq.variable_locations.h) == gs["variable_locations"]["h"]
            @test _ark_loc_name(faq.variable_locations.u) == gs["variable_locations"]["u"]
            @test _ark_loc_name(faq.variable_locations.v) == gs["variable_locations"]["v"]
            @test collect(faq.variable_shapes.u) == Int.(gs["variable_shapes"]["u"])
            @test collect(faq.variable_shapes.v) == Int.(gs["variable_shapes"]["v"])
            @test collect(faq.location_shapes.corner) == Int.(gs["location_shapes"]["corner"])
        end
    end
end

@testitem "Arakawa construction FAQ document is schema-valid" tags = [:conformance, :grid, :arakawa, :construction, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "arakawa", "rules", "arakawa_construction.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion.
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    model = file.models["ArakawaConstruction"]
    # The interval index sets that seed the construction are present.
    @test haskey(model.index_sets, "cells_x")
    @test haskey(model.index_sets, "nodes_x")
    @test haskey(model.index_sets, "staggers")
    # The two staggered axis maps + the static location table are declared.
    @test haskey(model.variables, "centers")
    @test haskey(model.variables, "nodes")
    @test haskey(model.variables, "u_location")
    @test haskey(model.variables, "v_location")
    @test haskey(model.variables, "neighbor_minus")
    @test haskey(model.variables, "boundary_lower")
end
