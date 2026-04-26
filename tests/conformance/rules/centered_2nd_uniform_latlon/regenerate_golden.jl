#!/usr/bin/env julia
#=
Regenerate golden coefficient + stencil-application values for the
centered_2nd_uniform_latlon cross-binding rule conformance harness.

The Julia binding is the reference evaluator per docs/GRIDS_API.md §4.3:
this script reads `fixtures.json`, evaluates each stencil entry's
coefficient AST under each fixture's bindings via
`EarthSciDiscretizations.eval_coeff` (a thin passthrough to
`EarthSciSerialization.evaluate`), applies the stencil to a closed-form
analytic field for the stencil cases, and writes the results to
`golden/coefficients.json`.

Run from the repo root:

    julia --project=. tests/conformance/rules/centered_2nd_uniform_latlon/regenerate_golden.jl

JSON I/O is intentionally hand-rolled here so the regenerator does not
add a runtime dependency on `JSON.jl`. The fixture / golden surface is
small and structured.
=#

using EarthSciDiscretizations: eval_coeff

const HERE = @__DIR__
const REPO_ROOT = abspath(joinpath(HERE, "..", "..", "..", ".."))
const FIXTURES_PATH = joinpath(HERE, "fixtures.json")
const GOLDEN_PATH = joinpath(HERE, "golden", "coefficients.json")

# ---------------------------------------------------------------------------
# Minimal JSON parser / writer (scalars, arrays, objects). Sufficient for
# the constrained shape of fixtures.json + golden output.
# ---------------------------------------------------------------------------

mutable struct JSONReader
    s::String
    pos::Int
end

_eof(r::JSONReader) = r.pos > lastindex(r.s)
_peek(r::JSONReader) = r.s[r.pos]

function _skip_ws!(r::JSONReader)
    while !_eof(r) && _peek(r) in (' ', '\t', '\r', '\n')
        r.pos = nextind(r.s, r.pos)
    end
end

function _expect!(r::JSONReader, c::Char)
    _skip_ws!(r)
    @assert !_eof(r) && _peek(r) == c "expected $c at pos $(r.pos), got $(_eof(r) ? "EOF" : string(_peek(r)))"
    r.pos = nextind(r.s, r.pos)
end

function _read_string!(r::JSONReader)::String
    _expect!(r, '"')
    io = IOBuffer()
    while !_eof(r)
        c = _peek(r)
        r.pos = nextind(r.s, r.pos)
        if c == '"'
            return String(take!(io))
        elseif c == '\\'
            esc = _peek(r); r.pos = nextind(r.s, r.pos)
            if esc == '"'; print(io, '"')
            elseif esc == '\\'; print(io, '\\')
            elseif esc == '/'; print(io, '/')
            elseif esc == 'b'; print(io, '\b')
            elseif esc == 'f'; print(io, '\f')
            elseif esc == 'n'; print(io, '\n')
            elseif esc == 'r'; print(io, '\r')
            elseif esc == 't'; print(io, '\t')
            else
                error("unsupported escape \\$(esc) at pos $(r.pos)")
            end
        else
            print(io, c)
        end
    end
    error("unterminated string")
end

function _read_number!(r::JSONReader)
    start = r.pos
    while !_eof(r) && (_peek(r) in "+-0123456789.eE")
        r.pos = nextind(r.s, r.pos)
    end
    s = SubString(r.s, start, prevind(r.s, r.pos))
    if any(c -> c in ".eE", s)
        return parse(Float64, s)
    end
    n = tryparse(Int, s)
    return n === nothing ? parse(Float64, s) : n
end

function _read_value!(r::JSONReader)
    _skip_ws!(r)
    c = _peek(r)
    if c == '{'
        r.pos = nextind(r.s, r.pos)
        obj = Pair{String,Any}[]
        _skip_ws!(r)
        if _peek(r) == '}'
            r.pos = nextind(r.s, r.pos)
            return Dict{String,Any}()
        end
        while true
            k = _read_string!(r)
            _expect!(r, ':')
            v = _read_value!(r)
            push!(obj, k => v)
            _skip_ws!(r)
            if _peek(r) == ','
                r.pos = nextind(r.s, r.pos)
            else
                _expect!(r, '}'); break
            end
        end
        return Dict{String,Any}(obj)
    elseif c == '['
        r.pos = nextind(r.s, r.pos)
        arr = Any[]
        _skip_ws!(r)
        if _peek(r) == ']'
            r.pos = nextind(r.s, r.pos)
            return arr
        end
        while true
            push!(arr, _read_value!(r))
            _skip_ws!(r)
            if _peek(r) == ','
                r.pos = nextind(r.s, r.pos)
            else
                _expect!(r, ']'); break
            end
        end
        return arr
    elseif c == '"'
        return _read_string!(r)
    elseif c == 't'
        @assert SubString(r.s, r.pos, r.pos + 3) == "true"
        r.pos += 4; return true
    elseif c == 'f'
        @assert SubString(r.s, r.pos, r.pos + 4) == "false"
        r.pos += 5; return false
    elseif c == 'n'
        @assert SubString(r.s, r.pos, r.pos + 3) == "null"
        r.pos += 4; return nothing
    else
        return _read_number!(r)
    end
end

function _parse_json_file(path::AbstractString)
    r = JSONReader(read(path, String), firstindex(read(path, String)))
    return _read_value!(r)
end

function _write_json(io::IO, v, indent::Int = 0)
    pad(n) = "  "^n
    if v isa AbstractDict
        if isempty(v); print(io, "{}"); return; end
        print(io, "{\n")
        keys_sorted = collect(keys(v))
        for (k, key) in enumerate(keys_sorted)
            print(io, pad(indent + 1), "\"", _escape(key), "\": ")
            _write_json(io, v[key], indent + 1)
            if k < length(keys_sorted); print(io, ","); end
            print(io, "\n")
        end
        print(io, pad(indent), "}")
    elseif v isa AbstractVector
        if isempty(v); print(io, "[]"); return; end
        print(io, "[\n")
        for (k, item) in enumerate(v)
            print(io, pad(indent + 1))
            _write_json(io, item, indent + 1)
            if k < length(v); print(io, ","); end
            print(io, "\n")
        end
        print(io, pad(indent), "]")
    elseif v isa AbstractString
        print(io, "\"", _escape(String(v)), "\"")
    elseif v isa Bool
        print(io, v ? "true" : "false")
    elseif v === nothing
        print(io, "null")
    elseif v isa Integer
        print(io, v)
    elseif v isa AbstractFloat
        # Round-trip preserving repr (matches stdlib JSON.print for finite floats).
        if isnan(v) || !isfinite(v)
            error("cannot serialize non-finite float to JSON: $v")
        end
        # Integer-valued floats: emit ".0" so they round-trip as floats in Python.
        if v == trunc(v) && abs(v) < 1e16
            print(io, Float64(v), "")  # Julia prints e.g. 1.0
        else
            print(io, repr(Float64(v)))
        end
    else
        error("unsupported JSON value type: $(typeof(v))")
    end
end

function _escape(s::AbstractString)
    io = IOBuffer()
    for c in s
        if c == '"'; print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\r'; print(io, "\\r")
        elseif c == '\t'; print(io, "\\t")
        else
            print(io, c)
        end
    end
    return String(take!(io))
end

# ---------------------------------------------------------------------------
# Rule + golden generation.
# ---------------------------------------------------------------------------

function _selector_key(sel::AbstractDict)
    axis = sel["axis"]
    offset = Int(sel["offset"])
    return string(axis, "_", offset >= 0 ? "p" : "m", abs(offset))
end

function _load_rule(path::AbstractString)
    raw = _parse_json_file(path)
    body_dict = raw["discretizations"]
    @assert length(body_dict) == 1
    name = first(keys(body_dict))
    return name, body_dict[name]
end

function _coefficient_block(rule_body, cases)
    out = []
    for case in cases
        bindings = Dict{String,Float64}(string(k) => Float64(v) for (k, v) in case["bindings"])
        coeffs = Dict{String,Float64}()
        for entry in rule_body["stencil"]
            sel = entry["selector"]
            coeffs[_selector_key(sel)] = eval_coeff(entry["coeff"], bindings)
        end
        push!(out, Dict{String,Any}(
            "name" => case["name"],
            "bindings" => case["bindings"],
            "coefficients" => coeffs,
        ))
    end
    return out
end

# Mirrors the per-cell binding convention used by Python's
# `apply_stencil_latlon`: cell-center coordinates from the `regular`
# lat-lon variant, `cos_lat` evaluated at the anchor cell's latitude.
function _stencil_block(rule_body, cases)
    out = []
    for case in cases
        nlat = Int(case["grid"]["nlat"])
        nlon = Int(case["grid"]["nlon"])
        R = Float64(case["grid"]["R"])
        dlon = 2π / nlon
        dlat = π / nlat

        lon_c(i) = (i + 0.5) * dlon
        lat_c(j) = -π / 2 + (j + 0.5) * dlat
        f(lon, lat) = sin(lon) * cos(lat)

        results = []
        for qp in case["interior_query_points"]
            j = Int(qp[1])
            i = Int(qp[2])
            cos_lat = cos(lat_c(j))
            bindings = Dict{String,Float64}(
                "R" => R, "dlon" => dlon, "dlat" => dlat, "cos_lat" => cos_lat,
            )
            total = 0.0
            for entry in rule_body["stencil"]
                sel = entry["selector"]
                axis = sel["axis"]
                offset = Int(sel["offset"])
                if axis == "lon"
                    ip = mod(i + offset, nlon)
                    val = f(lon_c(ip), lat_c(j))
                elseif axis == "lat"
                    jp = j + offset
                    @assert 0 <= jp < nlat "lat-axis offset escaped grid for j=$(j), offset=$(offset)"
                    val = f(lon_c(i), lat_c(jp))
                else
                    error("unexpected axis $(axis) in centered_2nd_uniform_latlon")
                end
                total += eval_coeff(entry["coeff"], bindings) * val
            end
            push!(results, Dict{String,Any}(
                "j" => j, "i" => i,
                "cos_lat" => cos_lat,
                "lon" => lon_c(i), "lat" => lat_c(j),
                "result" => total,
            ))
        end
        push!(out, Dict{String,Any}(
            "name" => case["name"],
            "grid" => case["grid"],
            "results" => results,
        ))
    end
    return out
end

function main()
    fixtures = _parse_json_file(FIXTURES_PATH)
    rule_path = joinpath(REPO_ROOT, fixtures["rule_path"])
    _, rule_body = _load_rule(rule_path)

    out = Dict{String,Any}(
        "version" => fixtures["version"],
        "rule" => fixtures["rule"],
        "reference_binding" => "julia",
        "coefficient_cases" => _coefficient_block(rule_body, fixtures["coefficient_cases"]),
        "stencil_cases" => _stencil_block(rule_body, fixtures["stencil_cases"]),
    )

    open(GOLDEN_PATH, "w") do io
        _write_json(io, out, 0)
        write(io, "\n")
    end
    println("Wrote $(GOLDEN_PATH)")
end

main()
