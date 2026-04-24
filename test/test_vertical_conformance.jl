@testsnippet VerticalConformanceSetup begin
    using Test
    using EarthSciDiscretizations
    using JSON

    const VERTICAL_CONFORMANCE_DIR = joinpath(
        @__DIR__, "..", "tests", "conformance", "grids", "vertical",
    )

    # Build the vertical grid from the fixture's declarative opts. Keys in
    # fixtures.json are strings; the Julia constructor takes the matching
    # kwargs as Symbols with the conventional types.
    function build_vertical_grid(opts::AbstractDict)
        kw = Dict{Symbol, Any}()
        # coordinate is required by fixtures.json for every vertical fixture.
        kw[:coordinate] = Symbol(opts["coordinate"])
        haskey(opts, "nz") && (kw[:nz] = Int(opts["nz"]))
        haskey(opts, "levels") &&
            (kw[:levels] = Float64[Float64(x) for x in opts["levels"]])
        haskey(opts, "ak") &&
            (kw[:ak] = Float64[Float64(x) for x in opts["ak"]])
        haskey(opts, "bk") &&
            (kw[:bk] = Float64[Float64(x) for x in opts["bk"]])
        haskey(opts, "p0") && (kw[:p0] = Float64(opts["p0"]))
        haskey(opts, "transition") && (kw[:transition] = Float64(opts["transition"]))
        haskey(opts, "ghosts") && (kw[:ghosts] = Int(opts["ghosts"]))
        return EarthSciDiscretizations.grids.vertical(; kw...)
    end

    # Strip binding-specific provenance fields declared by fixtures.json so
    # the remaining document is comparable across bindings.
    function strip_provenance(doc::AbstractDict, strip_fields)
        out = Dict{String, Any}()
        for (k, v) in doc
            if k == "provenance" && v isa AbstractDict
                keep = Dict{String, Any}()
                for (pk, pv) in v
                    pk in strip_fields && continue
                    keep[pk] = pv
                end
                out[k] = keep
            else
                out[k] = v
            end
        end
        return out
    end

    # Render a finite number to Python's `repr(float)` form so the canonical
    # output of `to_esm` plus the re-serialized golden agree byte-for-byte on
    # the Julia side. `Base.string(::Float64)` already emits the shortest
    # round-trip representation and appends `.0` to integer-valued floats,
    # matching Python's `json.dumps` rules for the values actually present in
    # the vertical fixture set (no scientific-notation floats appear).
    function _vconf_fmt_number(x::Real)
        if x isa Integer
            return string(x)
        end
        xf = Float64(x)
        isfinite(xf) || error(
            "vertical conformance canonicalization: non-finite number $xf",
        )
        return string(xf)
    end

    # Minimal JSON serializer that sorts object keys and uses 2-space indent,
    # matching Python's `json.dumps(..., sort_keys=True, indent=2)`. Both the
    # emitted `to_esm()` dict and the re-parsed golden pass through this same
    # canonicalizer, so byte-equality reduces to value-equality of the two
    # stripped documents.
    function _vconf_serialize(v, indent::Int)::String
        if v === nothing
            return "null"
        elseif v isa Bool
            return v ? "true" : "false"
        elseif v isa Integer
            return string(v)
        elseif v isa Real
            return _vconf_fmt_number(v)
        elseif v isa AbstractString
            return JSON.json(v)
        elseif v isa AbstractVector
            isempty(v) && return "[]"
            inner = " "^(indent + 2)
            parts = [inner * _vconf_serialize(x, indent + 2) for x in v]
            return "[\n" * join(parts, ",\n") * "\n" * " "^indent * "]"
        elseif v isa AbstractDict
            isempty(v) && return "{}"
            keys_sorted = sort!(String[String(k) for k in keys(v)])
            inner = " "^(indent + 2)
            parts = [
                inner * JSON.json(k) * ": " *
                    _vconf_serialize(v[k], indent + 2)
                    for k in keys_sorted
            ]
            return "{\n" * join(parts, ",\n") * "\n" * " "^indent * "}"
        end
        error(
            "vertical conformance canonicalization: unsupported value type " *
                "$(typeof(v))",
        )
    end

    function canonicalize_vertical(doc::AbstractDict, strip_fields)
        return _vconf_serialize(strip_provenance(doc, strip_fields), 0) * "\n"
    end
end

@testitem "vertical cross-language conformance" setup = [VerticalConformanceSetup] tags = [:grid, :vertical, :conformance] begin
    fixtures = JSON.parsefile(joinpath(VERTICAL_CONFORMANCE_DIR, "fixtures.json"))
    rel_tol = Float64(fixtures["tolerance"]["relative"])
    @test rel_tol == 0.0

    golden_dir = normpath(joinpath(VERTICAL_CONFORMANCE_DIR, fixtures["golden_dir"]))
    strip_fields = Set(String.(fixtures["provenance_strip_fields"]))

    for fixture in fixtures["fixtures"]
        name = fixture["name"]
        @testset "fixture=$name" begin
            grid = build_vertical_grid(fixture["opts"])
            emitted = canonicalize_vertical(to_esm(grid), strip_fields)

            golden_path = joinpath(golden_dir, "$name.esm")
            golden_raw = JSON.parsefile(golden_path)
            expected = canonicalize_vertical(golden_raw, strip_fields)

            @test emitted == expected
        end
    end
end

@testitem "vertical conformance: every committed fixture is in the corpus" setup = [VerticalConformanceSetup] tags = [:grid, :vertical, :conformance] begin
    fixtures = JSON.parsefile(joinpath(VERTICAL_CONFORMANCE_DIR, "fixtures.json"))
    golden_dir = normpath(joinpath(VERTICAL_CONFORMANCE_DIR, fixtures["golden_dir"]))
    committed = sort!(
        String[
            splitext(f)[1] for f in readdir(golden_dir)
                if endswith(f, ".esm")
        ]
    )
    corpus = sort!(String[String(f["name"]) for f in fixtures["fixtures"]])
    @test committed == corpus
end
