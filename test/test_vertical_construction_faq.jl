@testsnippet VerticalConstructionFAQSetup begin
    using Test
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON
    import EarthSciSerialization

    # Build a VerticalGrid from a fixture opts dict — the uniform (coordinate + nz),
    # supplied-levels (coordinate + levels) and hybrid (coordinate + ak/bk/p0) paths,
    # identical to tests/conformance/grids/vertical/construction/regenerate_golden.jl.
    function _vert_build_grid(opts::AbstractDict)
        kw = Dict{Symbol, Any}(:coordinate => Symbol(opts["coordinate"]))
        haskey(opts, "nz") && (kw[:nz] = Int(opts["nz"]))
        haskey(opts, "levels") && (kw[:levels] = Float64.(opts["levels"]))
        haskey(opts, "ak") && (kw[:ak] = Float64.(opts["ak"]))
        haskey(opts, "bk") && (kw[:bk] = Float64.(opts["bk"]))
        haskey(opts, "p0") && (kw[:p0] = Float64(opts["p0"]))
        haskey(opts, "ghosts") && (kw[:ghosts] = Int(opts["ghosts"]))
        return ESD.grids.vertical(; kw...)
    end

    # Serialization helpers — MUST match
    # tests/conformance/grids/vertical/construction/regenerate_golden.jl byte-for-byte
    # so the golden comparison is apples-to-apples.
    _vert_compact(x) = JSON.json(x)
    _vert_floats(v) = _vert_compact(Float64.(v))
    _vert_nbr(v) = _vert_compact([id == 0 ? -1 : Int(id) - 1 for id in v])
    _vert_mask(v) = _vert_compact([b ? 1 : 0 for b in v])

    # Bit-exact equality on float vectors ("match imperative to ULP").
    _vert_biteq(a::AbstractVector, b::AbstractVector) =
        length(a) == length(b) &&
            all(reinterpret(UInt64, Float64(a[i])) == reinterpret(UInt64, Float64(b[i]))
                for i in eachindex(a))
end

@testitem "Vertical construction FAQ — matches imperative vertical.jl (ULP) + byte golden" setup = [VerticalConstructionFAQSetup] tags = [:conformance, :grid, :vertical, :construction, :faq] begin
    base = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "vertical")
    SPEC = JSON.parsefile(joinpath(base, "fixtures.json"))
    GOLDEN = JSON.parsefile(joinpath(base, "construction", "golden.json"))
    by_name = Dict(f["name"] => f for f in GOLDEN["fixtures"])

    for f in SPEC["fixtures"]
        name = f["name"]
        g = _vert_build_grid(f["opts"])
        faq = vertical_construction_faq(g)
        gl = by_name[name]
        nz = ESD.n_cells(g)

        # --- counts ---
        @test nz == gl["n_cells"]
        @test ESD.n_vertices(g) == gl["n_vertices"]
        @test ESD.n_edges(g) == gl["n_edges"]
        @test String(g.coordinate) == gl["coordinate"]

        # --- "match imperative vertical.jl to ULP": every FAQ array equals the
        # imperative trait array / accessor bit-for-bit. ---
        @test _vert_biteq(faq.levels, g.levels)
        @test _vert_biteq(faq.centers, ESD.cell_centers(g))
        @test _vert_biteq(faq.widths, ESD.cell_widths(g))
        @test _vert_biteq(faq.cell_volume, ESD.cell_volume(g))
        @test _vert_biteq(faq.half_levels, ESD.half_levels(g))
        @test _vert_biteq(faq.layer_thickness, ESD.layer_thickness(g))
        @test _vert_biteq(faq.centers, ESD.cell_centers(g, :z))
        @test _vert_biteq(faq.widths, ESD.cell_widths(g, :z))
        @test faq.neighbor_minus == ESD.neighbor_indices(g, :z, -1)
        @test faq.neighbor_plus == ESD.neighbor_indices(g, :z, +1)
        @test faq.boundary_lower == ESD.boundary_mask(g, :z, :lower)
        @test faq.boundary_upper == ESD.boundary_mask(g, :z, :upper)

        # --- named layer metrics: match metric_eval per available name ---
        @test _vert_biteq(faq.metric_dz, [ESD.metric_eval(g, :dz, k) for k in 1:nz])
        @test _vert_biteq(faq.metric_z, [ESD.metric_eval(g, :z, k) for k in 1:nz])
        sigma_like = g.coordinate in (:sigma, :hybrid_sigma_theta, :eta)
        @test (faq.metric_sigma !== nothing) == sigma_like
        if sigma_like
            @test _vert_biteq(faq.metric_sigma, [ESD.metric_eval(g, :sigma, k) for k in 1:nz])
        end
        has_ak = length(g.ak) > 0
        has_bk = length(g.bk) > 0
        @test (faq.metric_pressure !== nothing) == (has_ak && has_bk)
        if has_ak && has_bk
            @test _vert_biteq(faq.metric_pressure, [ESD.metric_eval(g, :pressure, k) for k in 1:nz])
        end
        @test (faq.metric_ak !== nothing) == has_ak
        if has_ak
            @test _vert_biteq(faq.metric_ak, [ESD.metric_eval(g, :ak, k) for k in 1:nz])
        end
        @test (faq.metric_bk !== nothing) == has_bk
        if has_bk
            @test _vert_biteq(faq.metric_bk, [ESD.metric_eval(g, :bk, k) for k in 1:nz])
        end

        # --- pressure_coefficients ---
        pc = ESD.pressure_coefficients(g)
        @test _vert_biteq(faq.ak, pc.ak)
        @test _vert_biteq(faq.bk, pc.bk)
        @test faq.p0 == pc.p0

        # --- internal consistency of the FAQ construction itself ---
        @test length(faq.levels) == nz + 1
        @test all(faq.widths[k] > 0 for k in 1:nz)
        @test all(faq.widths[k] == abs(faq.levels[k + 1] - faq.levels[k]) for k in 1:nz)
        decreasing = faq.levels[2] < faq.levels[1]
        if decreasing
            @test all(faq.levels[k + 1] < faq.levels[k] for k in 1:nz)   # sigma-like descend
        else
            @test all(faq.levels[k + 1] > faq.levels[k] for k in 1:nz)   # z/theta ascend
        end
        @test all(
            min(faq.levels[k], faq.levels[k + 1]) < faq.centers[k] <
                max(faq.levels[k], faq.levels[k + 1]) for k in 1:nz)      # center strictly inside

        # --- neighbor 0-sentinel + symmetry ---
        for k in 1:nz
            @test (faq.neighbor_minus[k] == 0) == (k == 1)               # 0-sentinel at bottom
            @test (faq.neighbor_plus[k] == 0) == (k == nz)               # 0-sentinel at top
            np = faq.neighbor_plus[k]
            np != 0 && @test faq.neighbor_minus[np] == k                 # +/- neighbor symmetry
        end

        # --- byte-identity: the binding-neutral compact serialization equals
        # golden.json. Julia is the reference binding; the vertical family has no
        # transcendental math (pure rationals + a single multiply-add), so the
        # float strings match exactly (tolerance 0.0). ---
        @test _vert_floats(faq.levels) == gl["levels"]
        @test _vert_floats(faq.centers) == gl["centers"]
        @test _vert_floats(faq.widths) == gl["widths"]
        @test _vert_floats(faq.cell_volume) == gl["cell_volume"]
        @test _vert_floats(faq.metric_dz) == gl["metric_dz"]
        @test _vert_floats(faq.metric_z) == gl["metric_z"]
        if haskey(gl, "metric_sigma")
            @test _vert_floats(faq.metric_sigma) == gl["metric_sigma"]
        end
        if haskey(gl, "metric_pressure")
            @test _vert_floats(faq.metric_pressure) == gl["metric_pressure"]
        end
        if haskey(gl, "metric_ak")
            @test _vert_floats(faq.metric_ak) == gl["metric_ak"]
        end
        if haskey(gl, "metric_bk")
            @test _vert_floats(faq.metric_bk) == gl["metric_bk"]
        end
        @test _vert_nbr(faq.neighbor_minus) == gl["neighbor_minus"]
        @test _vert_nbr(faq.neighbor_plus) == gl["neighbor_plus"]
        @test _vert_mask(faq.boundary_lower) == gl["boundary_lower"]
        @test _vert_mask(faq.boundary_upper) == gl["boundary_upper"]
        if haskey(gl, "ak")
            @test _vert_floats(faq.ak) == gl["ak"]
        end
        if haskey(gl, "bk")
            @test _vert_floats(faq.bk) == gl["bk"]
        end
        @test Float64(faq.p0) == gl["p0"]
    end
end

@testitem "Vertical construction FAQ anchors level synthesis to ESS evaluation" setup = [VerticalConstructionFAQSetup] tags = [:conformance, :grid, :vertical, :construction, :faq] begin
    # Uniform sigma synthesis levels[k] = 1 - (k-1)/nz. With nz=16 every interface
    # is a clean dyadic rational (k/16 has ≤4 fractional bits), reproduced
    # bit-identically by any binding's ESS evaluator.
    g = ESD.grids.vertical(coordinate = :sigma, nz = 16)
    faq = vertical_construction_faq(g)
    @test faq.levels == [
        1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625,
        0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625, 0.0,
    ]
    @test all(faq.widths .== 0.0625)
    @test faq.cell_volume == faq.widths            # 1-D measure is the thickness
    @test faq.metric_dz == faq.widths
    @test faq.metric_z == faq.centers
    @test faq.metric_sigma == faq.centers          # sigma-like ⇒ sigma metric = center
    @test faq.metric_pressure === nothing          # no hybrid coefficients
    @test isempty(faq.ak) && isempty(faq.bk)

    # Eta hybrid synthesis levels[k] = ak[k]/p0 + bk[k] (the multiply-add path),
    # bit-identical to the imperative generator and exercising the hybrid metrics.
    ge = ESD.grids.vertical(
        coordinate = :eta,
        ak = [0.0, 219.4067, 489.5209],
        bk = [1.0, 0.9851122, 0.963406],
        p0 = 1.0e5,
    )
    faqe = vertical_construction_faq(ge)
    @test faqe.levels[1] == 1.0                     # ak=0, bk=1 ⇒ surface sigma 1
    @test faqe.levels == ge.levels                  # ESS synthesis == imperative
    @test faqe.metric_pressure !== nothing
    @test faqe.metric_pressure == [ESD.metric_eval(ge, :pressure, k) for k in 1:ESD.n_cells(ge)]
    @test faqe.metric_ak == [ESD.metric_eval(ge, :ak, k) for k in 1:ESD.n_cells(ge)]
    @test faqe.metric_bk == [ESD.metric_eval(ge, :bk, k) for k in 1:ESD.n_cells(ge)]
end

@testitem "Vertical construction FAQ document is schema-valid" tags = [:conformance, :grid, :vertical, :construction, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "vertical", "rules", "vertical_construction.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion.
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    # Sigma model: the two interval index sets + the generative/derived variables.
    sm = file.models["VerticalSigmaConstruction"]
    @test haskey(sm.index_sets, "interfaces")
    @test haskey(sm.index_sets, "layers")
    @test haskey(sm.variables, "levels")
    @test haskey(sm.variables, "centers")
    @test haskey(sm.variables, "neighbor_minus")
    @test haskey(sm.variables, "boundary_lower")
    # Hybrid model: the eta synthesis + the pressure/ak/bk layer-average metrics.
    hm = file.models["VerticalHybridConstruction"]
    @test haskey(hm.variables, "levels")
    @test haskey(hm.variables, "metric_pressure")
    @test haskey(hm.variables, "metric_ak")
    @test haskey(hm.variables, "metric_bk")
end
