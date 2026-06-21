@testsnippet LatLonConstructionFAQSetup begin
    using Test
    using EarthSciDiscretizations
    const ESD = EarthSciDiscretizations
    using JSON
    import EarthSciSerialization

    # Build a LatLonGrid from a fixture opts dict — the regular (nlon/nlat) and
    # reduced_gaussian (nlon_per_row + supplied lat_edges) paths, identical to
    # tests/conformance/grids/latlon/construction/regenerate_golden.jl.
    function _ll_build_grid(opts::AbstractDict)
        variant = Symbol(get(opts, "variant", "regular"))
        dtype = get(opts, "dtype", "float64") == "float32" ? Float32 : Float64
        kw = Dict{Symbol, Any}(:variant => variant, :dtype => dtype)
        kw[:R] = Float64(opts["R"])
        kw[:ghosts] = Int(get(opts, "ghosts", 0))
        haskey(opts, "pole_policy") && (kw[:pole_policy] = Symbol(opts["pole_policy"]))
        if variant === :regular
            kw[:nlon] = Int(opts["nlon"])
            kw[:nlat] = Int(opts["nlat"])
        else
            kw[:nlon_per_row] = Int[Int(x) for x in opts["nlon_per_row"]]
            haskey(opts, "lat_edges") &&
                (kw[:lat_edges] = Float64[Float64(x) for x in opts["lat_edges"]])
        end
        return ESD.grids.lat_lon(; kw...)
    end

    # Serialization helpers — MUST match
    # tests/conformance/grids/latlon/construction/regenerate_golden.jl byte-for-byte.
    _ll_compact(x) = JSON.json(x)
    _ll_floats(v) = _ll_compact([Float64(x) for x in v])
    _ll_comp(a, idx...) = _ll_compact([Float64(a[k, idx...]) for k in 1:size(a, 1)])
    _ll_nbr(v) = _ll_compact([id == 0 ? -1 : Int(id) - 1 for id in v])
    _ll_mask(v) = _ll_compact([b ? 1 : 0 for b in v])

    # Bit-exact equality on float arrays (the FAQ output is bit-identical to the
    # imperative builder, not merely close — this is "match imperative to ULP").
    _ll_biteq(a::AbstractArray, b::AbstractArray) =
        size(a) == size(b) &&
            all(reinterpret(UInt64, collect(Float64.(a))) .== reinterpret(UInt64, collect(Float64.(b))))
end

@testitem "Lat-lon construction FAQ — matches imperative latlon.jl (ULP) + byte golden" setup = [LatLonConstructionFAQSetup] tags = [:conformance, :grid, :latlon, :construction, :faq] begin
    base = joinpath(@__DIR__, "..", "tests", "conformance", "grids", "latlon")
    SPEC = JSON.parsefile(joinpath(base, "fixtures.json"))
    GOLDEN = JSON.parsefile(joinpath(base, "construction", "golden.json"))
    by_name = Dict(f["name"] => f for f in GOLDEN["fixtures"])

    for f in SPEC["fixtures"]
        name = f["name"]
        g = _ll_build_grid(f["opts"])
        faq = latlon_construction_faq(g)
        gl = by_name[name]
        nc = ESD.n_cells(g)

        # --- counts ---
        @test g.nlat == gl["nlat"]
        @test nc == gl["n_cells"]
        @test Int[x for x in g.nlon_per_row] == Int.(gl["nlon_per_row"])

        # --- "match imperative latlon.jl to ULP": every FAQ array equals the
        # imperative trait array bit-for-bit (coords, trig metric, neighbors, boundary). ---
        @test _ll_biteq(faq.cell_centers_lon, ESD.cell_centers(g, :lon))
        @test _ll_biteq(faq.cell_centers_lat, ESD.cell_centers(g, :lat))
        @test _ll_biteq(faq.cell_widths_lon, ESD.cell_widths(g, :lon))
        @test _ll_biteq(faq.cell_widths_lat, ESD.cell_widths(g, :lat))
        @test _ll_biteq(faq.cell_volume, ESD.cell_volume(g))
        @test _ll_biteq(faq.lat_edges, g.lat_edges)
        @test _ll_biteq(faq.lat_centers, g.lat_centers)
        @test _ll_biteq(faq.metric_g, ESD.metric_g(g))
        @test _ll_biteq(faq.metric_ginv, ESD.metric_ginv(g))
        @test _ll_biteq(faq.metric_jacobian, ESD.metric_jacobian(g))
        @test _ll_biteq(faq.metric_dgij_dxk, ESD.metric_dgij_dxk(g))
        @test _ll_biteq(faq.coord_jacobian, ESD.coord_jacobian(g, :lon_lat))
        @test _ll_biteq(faq.coord_jacobian_second, ESD.coord_jacobian_second(g, :lon_lat))
        @test faq.neighbor_lon_minus == ESD.neighbor_indices(g, :lon, -1)
        @test faq.neighbor_lon_plus == ESD.neighbor_indices(g, :lon, +1)
        @test faq.neighbor_lat_minus == ESD.neighbor_indices(g, :lat, -1)
        @test faq.neighbor_lat_plus == ESD.neighbor_indices(g, :lat, +1)
        @test faq.boundary_lat_lower == ESD.boundary_mask(g, :lat, :lower)
        @test faq.boundary_lat_upper == ESD.boundary_mask(g, :lat, :upper)
        @test faq.boundary_lon_lower == ESD.boundary_mask(g, :lon, :lower)
        @test faq.boundary_lon_upper == ESD.boundary_mask(g, :lon, :upper)

        # --- internal consistency of the FAQ construction itself ---
        @test length(faq.cell_volume) == nc
        @test all(faq.cell_volume[k] > 0 for k in 1:nc)             # positive areas
        @test all(faq.cell_widths_lon[k] > 0 for k in 1:nc)         # positive lon width
        @test all(faq.cell_widths_lat[k] > 0 for k in 1:nc)         # positive lat width
        @test length(faq.lat_edges) == g.nlat + 1
        @test issorted(faq.lat_edges)                               # latitude edges ascend
        # orthogonal lat-lon metric: g is diagonal, g_φφ = R², J = R²|cos φ|, ginv inverts g
        R2 = Float64(g.R)^2
        @test all(faq.metric_g[k, 1, 2] == 0.0 && faq.metric_g[k, 2, 1] == 0.0 for k in 1:nc)
        @test all(faq.metric_g[k, 2, 2] == R2 for k in 1:nc)
        @test all(isapprox(faq.metric_g[k, 1, 1] * faq.metric_ginv[k, 1, 1], 1.0; rtol = 1e-12)
                  for k in 1:nc)
        @test all(faq.metric_jacobian[k] ≥ 0 for k in 1:nc)
        @test all(iszero, faq.coord_jacobian_second)
        @test all(faq.coord_jacobian[k, 1, 1] == 1.0 && faq.coord_jacobian[k, 2, 2] == 1.0
                  for k in 1:nc)
        # longitude wraps → no lon boundary cells
        @test !any(faq.boundary_lon_lower)
        @test !any(faq.boundary_lon_upper)

        # --- neighbor periodicity + pole 0-sentinel (flat ragged-row-major ids) ---
        # longitude is periodic: every cell has a + and − neighbor, and they are
        # mutually inverse (no 0-sentinel ever in longitude).
        @test all(faq.neighbor_lon_minus[k] != 0 && faq.neighbor_lon_plus[k] != 0 for k in 1:nc)
        for k in 1:nc
            @test faq.neighbor_lon_minus[faq.neighbor_lon_plus[k]] == k
        end
        # latitude: 0-sentinel exactly on the first row's − side and last row's + side.
        row_off = [sum(@view(g.nlon_per_row[1:(j - 1)])) for j in 1:(g.nlat + 1)]
        for j in 1:g.nlat
            for i in 1:g.nlon_per_row[j]
                k = row_off[j] + i
                @test (faq.neighbor_lat_minus[k] == 0) == (j == 1)
                @test (faq.neighbor_lat_plus[k] == 0) == (j == g.nlat)
            end
        end

        # --- spot-check the per-(j,i) accessors against the bulk FAQ arrays ---
        for j in 1:g.nlat, i in (1, g.nlon_per_row[j])
            k = row_off[j] + i
            lon, lat = ESD.cell_centers(g, j, i)
            @test faq.cell_centers_lon[k] == lon
            @test faq.cell_centers_lat[k] == lat
            @test faq.cell_volume[k] == ESD.cell_area(g, j, i)
            @test faq.metric_jacobian[k] == ESD.metric_eval(g, :J, j, i)
            @test faq.metric_g[k, 1, 1] == ESD.metric_eval(g, :g_lonlon, j, i)
        end

        # --- byte-identity: the binding-neutral serialization equals golden.json. ---
        @test _ll_floats(faq.cell_centers_lon) == gl["cell_centers_lon"]
        @test _ll_floats(faq.cell_centers_lat) == gl["cell_centers_lat"]
        @test _ll_floats(faq.cell_widths_lon) == gl["cell_widths_lon"]
        @test _ll_floats(faq.cell_widths_lat) == gl["cell_widths_lat"]
        @test _ll_floats(faq.cell_volume) == gl["cell_volume"]
        @test _ll_floats(faq.lat_edges) == gl["lat_edges"]
        @test _ll_floats(faq.lat_centers) == gl["lat_centers"]
        @test _ll_comp(faq.metric_g, 1, 1) == gl["metric_g_lonlon"]
        @test _ll_comp(faq.metric_g, 2, 2) == gl["metric_g_latlat"]
        @test _ll_comp(faq.metric_ginv, 1, 1) == gl["metric_ginv_lonlon"]
        @test _ll_floats(faq.metric_jacobian) == gl["metric_jacobian"]
        @test _ll_comp(faq.metric_dgij_dxk, 1, 1, 2) == gl["dg_lonlon_dlat"]
        @test _ll_nbr(faq.neighbor_lon_minus) == gl["neighbor_lon_minus"]
        @test _ll_nbr(faq.neighbor_lon_plus) == gl["neighbor_lon_plus"]
        @test _ll_nbr(faq.neighbor_lat_minus) == gl["neighbor_lat_minus"]
        @test _ll_nbr(faq.neighbor_lat_plus) == gl["neighbor_lat_plus"]
        @test _ll_mask(faq.boundary_lat_lower) == gl["boundary_lat_lower"]
        @test _ll_mask(faq.boundary_lat_upper) == gl["boundary_lat_upper"]
    end
end

@testitem "Lat-lon construction FAQ anchors coords + trig metric to ESS evaluation" setup = [LatLonConstructionFAQSetup] tags = [:conformance, :grid, :latlon, :construction, :faq] begin
    # The per-row affine longitude map and the trig metric are evaluated by the ESS
    # AST evaluator (the single pathway), not a binding-local loop. A regular 4×4 grid
    # on R=2, lon_start=-π gives dlon = 2π/4 = π/2; the first-row centers are clean
    # multiples of π/4, reproduced bit-identically by any binding's ESS evaluator.
    g = ESD.grids.lat_lon(variant = :regular, nlon = 4, nlat = 4, R = 2.0)
    faq = latlon_construction_faq(g)
    dlon = 2 * π / 4
    @test faq.cell_widths_lon[1] == dlon
    @test faq.cell_centers_lon[1] == -π + 0.5 * dlon
    @test faq.cell_centers_lon[2] == -π + 1.5 * dlon
    @test faq.cell_centers_lon[3] == -π + 2.5 * dlon
    @test faq.cell_centers_lon[4] == -π + 3.5 * dlon
    # closed-form trig metric at each row center: g_λλ = R² cos²φ, J = R² |cos φ|.
    for j in 1:g.nlat
        k = (j - 1) * 4 + 1
        φ = g.lat_centers[j]
        @test faq.metric_g[k, 1, 1] == 2.0 * 2.0 * cos(φ) * cos(φ)
        @test faq.metric_g[k, 2, 2] == 2.0 * 2.0
        @test faq.metric_jacobian[k] == 2.0 * 2.0 * abs(cos(φ))
        @test faq.metric_dgij_dxk[k, 1, 1, 2] == -2.0 * 2.0 * 2.0 * cos(φ) * sin(φ)
    end
    # north/south symmetry: equal-angle rows j and nlat+1-j have ±φ ⇒ equal cos² metric
    # (isapprox, not bit-equal: -π/2 + k·dlat accumulates a sub-ULP asymmetry).
    @test faq.metric_g[1, 1, 1] ≈ faq.metric_g[(g.nlat - 1) * 4 + 1, 1, 1] rtol = 1e-12
end

@testitem "Lat-lon construction FAQ reduced-Gaussian remap matches imperative" setup = [LatLonConstructionFAQSetup] tags = [:conformance, :grid, :latlon, :construction, :faq] begin
    # The nearest-center N/S row remap (_latlon_map_i) is the reduced-Gaussian
    # signature: rows of differing nlon. The FAQ neighbor maps must reproduce the
    # imperative neighbor_indices exactly, including the rank/join remap.
    g = ESD.grids.lat_lon(
        variant = :reduced_gaussian,
        nlon_per_row = [4, 8, 12, 16, 16, 12, 8, 4],
        lat_edges = collect(range(-π / 2, π / 2; length = 9)),
        R = 1.0,
    )
    faq = latlon_construction_faq(g)
    @test faq.neighbor_lat_minus == ESD.neighbor_indices(g, :lat, -1)
    @test faq.neighbor_lat_plus == ESD.neighbor_indices(g, :lat, +1)
    @test faq.neighbor_lon_minus == ESD.neighbor_indices(g, :lon, -1)
    @test faq.neighbor_lon_plus == ESD.neighbor_indices(g, :lon, +1)
    # rows differ in nlon ⇒ at least one cross-row remap is NOT the identity column.
    @test g.nlon_per_row[1] != g.nlon_per_row[2]
end

@testitem "Lat-lon construction FAQ document is schema-valid" tags = [:conformance, :grid, :latlon, :construction, :faq, :schema] begin
    import EarthSciSerialization
    path = joinpath(@__DIR__, "..", "discretizations", "grids", "latlon", "rules", "latlon_construction.esm")
    # `load` performs JSON-schema validation and throws on a schema error, so a
    # successful load is itself the schema-validity assertion.
    file = EarthSciSerialization.load(path)
    @test file isa EarthSciSerialization.EsmFile
    res = EarthSciSerialization.validate(file)
    @test isempty(res.schema_errors)
    # The four interval index sets that seed the construction are present.
    model = file.models["LatLonConstruction"]
    @test haskey(model.index_sets, "lat_nodes")
    @test haskey(model.index_sets, "rows")
    @test haskey(model.index_sets, "lon_cells")
    @test haskey(model.index_sets, "cells")
    # The coordinate, metric, neighbor and boundary families are declared.
    @test haskey(model.variables, "lat_centers")
    @test haskey(model.variables, "lon_centers")
    @test haskey(model.variables, "cell_area")
    @test haskey(model.variables, "metric_g_lonlon")
    @test haskey(model.variables, "dg_lonlon_dlat")
    @test haskey(model.variables, "neighbor_lon_plus")
    @test haskey(model.variables, "neighbor_lat_plus")
    @test haskey(model.variables, "boundary_lat_lower")
end
