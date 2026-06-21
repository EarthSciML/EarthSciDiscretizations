@testsnippet CartesianConstructionFAQSetup begin
    using Test
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON
    import EarthSciSerialization

    const _CART_AXES = (:x, :y, :z)

    # Build a CartesianGrid from a fixture opts dict — the uniform (nx/ny/nz +
    # extent) and non-uniform (edges) paths, identical to
    # tests/conformance/grids/cartesian/construction/regenerate_golden.jl.
    function _cart_build_grid(opts::AbstractDict)
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
        return ESD.grids.cartesian(; kw...)
    end

    # Serialization helpers — MUST match
    # tests/conformance/grids/cartesian/construction/regenerate_golden.jl
    # byte-for-byte so the golden comparison is apples-to-apples. Each construction
    # array family is one compact JSON string of the [axis][...] nesting, the same
    # dense byte form as tests/conformance/grids/duo/topology/golden.json.
    _cart_compact(x) = JSON.json(x)
    _cart_axis_floats(faq, field, N) = _cart_compact([Float64.(getfield(faq, field)[d]) for d in 1:N])
    _cart_axis_nbr(faq, field, N) =
        _cart_compact([[id == 0 ? -1 : Int(id) - 1 for id in getfield(faq, field)[d]] for d in 1:N])
    _cart_axis_mask(faq, field, N) =
        _cart_compact([[b ? 1 : 0 for b in getfield(faq, field)[d]] for d in 1:N])

    # Bit-exact equality on float arrays (the FAQ output is bit-identical to the
    # imperative builder, not merely close — this is "match imperative to ULP").
    _cart_biteq(a::AbstractArray, b::AbstractArray) =
        size(a) == size(b) &&
        all(reinterpret(UInt64, collect(Float64.(a))) .== reinterpret(UInt64, collect(Float64.(b))))
end

@testitem "Cartesian construction FAQ — matches imperative cartesian.jl (ULP) + byte golden" setup = [CartesianConstructionFAQSetup] tags = [:conformance, :grid, :cartesian, :construction, :faq] begin
    base = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "cartesian")
    SPEC = JSON.parsefile(joinpath(base, "fixtures.json"))
    GOLDEN = JSON.parsefile(joinpath(base, "construction", "golden.json"))
    by_name = Dict(f["name"] => f for f in GOLDEN["fixtures"])

    for f in SPEC["fixtures"]
        name = f["name"]
        g = _cart_build_grid(f["opts"])
        N = ndims(g)
        faq = cartesian_construction_faq(g)
        gl = by_name[name]

        # --- counts ---
        @test N == gl["ndim"]
        @test collect(g.n) == Int.(gl["n"])
        @test prod(g.n) == gl["n_cells"]

        # --- "match imperative cartesian.jl to ULP": every FAQ array equals the
        # imperative trait array bit-for-bit (coords, metric, neighbors, boundary). ---
        for d in 1:N
            ax = _CART_AXES[d]
            @test _cart_biteq(faq.edges[d], g.edges[d])
            @test _cart_biteq(faq.centers[d], g.centers[d])
            @test _cart_biteq(faq.widths[d], g.widths[d])
            @test _cart_biteq(faq.cell_centers[d], ESD.cell_centers(g, ax))
            @test _cart_biteq(faq.cell_widths[d], ESD.cell_widths(g, ax))
            @test faq.neighbor_minus[d] == ESD.neighbor_indices(g, ax, -1)
            @test faq.neighbor_plus[d] == ESD.neighbor_indices(g, ax, +1)
            @test faq.boundary_lower[d] == ESD.boundary_mask(g, ax, :lower)
            @test faq.boundary_upper[d] == ESD.boundary_mask(g, ax, :upper)
        end
        @test _cart_biteq(faq.cell_volume, ESD.cell_volume(g))
        @test _cart_biteq(faq.metric_g, ESD.metric_g(g))
        @test _cart_biteq(faq.metric_ginv, ESD.metric_ginv(g))
        @test _cart_biteq(faq.metric_jacobian, ESD.metric_jacobian(g))
        @test _cart_biteq(faq.metric_dgij_dxk, ESD.metric_dgij_dxk(g))
        @test _cart_biteq(faq.coord_jacobian, ESD.coord_jacobian(g, :cartesian))
        @test _cart_biteq(faq.coord_jacobian_second, ESD.coord_jacobian_second(g, :cartesian))

        # --- internal consistency of the FAQ construction itself ---
        for d in 1:N
            e, c, w = faq.edges[d], faq.centers[d], faq.widths[d]
            nd = g.n[d]
            @test length(e) == nd + 1
            @test issorted(e)                                   # affine / supplied edges ascend
            @test all(w[i] > 0 for i in 1:nd)                   # positive widths
            @test all(e[i] < c[i] < e[i + 1] for i in 1:nd)     # center strictly inside its cell
            @test all(w[i] == e[i + 1] - e[i] for i in 1:nd)    # width is the face difference
            # identity metric: g and ginv are the Kronecker delta, derivatives vanish
        end
        nc = prod(g.n)
        @test all(faq.metric_g[k, i, j] == (i == j ? 1.0 : 0.0) for k in 1:nc, i in 1:N, j in 1:N)
        @test all(iszero, faq.metric_dgij_dxk)
        @test all(iszero, faq.coord_jacobian_second)
        @test faq.metric_jacobian == faq.cell_volume            # identity metric ⇒ J = measure

        # --- neighbor 0-sentinel + symmetry (per axis, flattened linear ids) ---
        for d in 1:N
            ci = CartesianIndices(g.n)
            for k in 1:nc
                first_d = ci[k][d] == 1
                last_d = ci[k][d] == g.n[d]
                @test (faq.neighbor_minus[d][k] == 0) == first_d   # 0-sentinel at the lower edge
                @test (faq.neighbor_plus[d][k] == 0) == last_d     # 0-sentinel at the upper edge
                np = faq.neighbor_plus[d][k]
                np != 0 && @test faq.neighbor_minus[d][np] == k    # +/- neighbor symmetry
            end
        end

        # --- byte-identity: the binding-neutral compact serialization equals
        # golden.json. Julia is the reference binding, so the float coordinate
        # strings match exactly (pure-arithmetic, deterministic Ryu repr); the
        # cross-binding contract compares floats at GOLDEN["tolerance"] and the
        # integer neighbor maps / boundary masks exactly. ---
        @test _cart_axis_floats(faq, :edges, N) == gl["edges"]
        @test _cart_axis_floats(faq, :centers, N) == gl["centers"]
        @test _cart_axis_floats(faq, :widths, N) == gl["widths"]
        @test _cart_compact(Float64.(faq.cell_volume)) == gl["cell_volume"]
        @test _cart_compact(Float64.(faq.metric_jacobian)) == gl["metric_jacobian"]
        @test _cart_axis_nbr(faq, :neighbor_minus, N) == gl["neighbor_minus"]
        @test _cart_axis_nbr(faq, :neighbor_plus, N) == gl["neighbor_plus"]
        @test _cart_axis_mask(faq, :boundary_lower, N) == gl["boundary_lower"]
        @test _cart_axis_mask(faq, :boundary_upper, N) == gl["boundary_upper"]
    end
end

@testitem "Cartesian construction FAQ anchors affine coords to ESS evaluation" setup = [CartesianConstructionFAQSetup] tags = [:conformance, :grid, :cartesian, :construction, :faq] begin
    # The affine coordinate map edges[i] = lo + (i-1)*dx is evaluated by the ESS
    # AST evaluator (the single pathway), not a binding-local loop. A uniform axis
    # on [0,1] with n=8 gives dx = 1/8 exactly; every face is a clean dyadic
    # rational, reproduced bit-identically by any binding's ESS evaluator.
    g = ESD.grids.cartesian(nx = 8, extent = [(0.0, 1.0)], ghosts = 0)
    faq = cartesian_construction_faq(g)
    @test faq.edges[1] == [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
    @test faq.centers[1] == [0.0625, 0.1875, 0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375]
    @test all(faq.widths[1] .== 0.125)
    # cell_volume in 1-D is the width; the identity-metric Jacobian equals it.
    @test faq.cell_volume == faq.widths[1]
    @test faq.metric_jacobian == faq.cell_volume
end

@testitem "Cartesian construction FAQ document is schema-valid" tags = [:conformance, :grid, :cartesian, :construction, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "cartesian", "rules", "cartesian_construction.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion.
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    # The two interval index sets that seed the construction are present.
    model = file.models["CartesianConstruction"]
    @test haskey(model.index_sets, "edge_nodes")
    @test haskey(model.index_sets, "cells")
    # The affine generators and the four construction families are declared.
    @test haskey(model.variables, "edges")
    @test haskey(model.variables, "centers")
    @test haskey(model.variables, "neighbor_minus")
    @test haskey(model.variables, "boundary_lower")
end
