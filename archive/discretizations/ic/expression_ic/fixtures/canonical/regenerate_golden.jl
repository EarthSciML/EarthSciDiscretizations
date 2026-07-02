#!/usr/bin/env julia
#=
Regenerate the Layer-A canonical golden (`expected.esm`) for the
`expression_ic` rule on the structured grid families (latlon + arakawa).

Per AGENTS.md "single-pathway rule": this script drives the canonical
pipeline — `EarthSciSerialization.discretize` → ArrayOp — and emits the
canonical whole-document JSON form. It carries NO shadow evaluator and no
binding-local IC math: the `ic`-node → arrayop materialization is owned
entirely by the ESS engine ("ic-arrayop-engine"). The walker
(`test/walk_esd_tests.jl`, `apply_rule_and_diff`) consumes the same
`discretize` output and byte-compares it to `expected.esm`, so this
regenerator and the walker agree by construction.

The canonical emitter below mirrors the walker's `canonical_doc_json`
(sorted keys, minified, RFC §5.4.6 float formatting via
`EarthSciSerialization.format_canonical_float`).

Run from the repo root:

    julia --project=. discretizations/ic/expression_ic/fixtures/canonical/regenerate_golden.jl
=#

using EarthSciDiscretizations  # ensure ESS is resolvable via the project env
using EarthSciSerialization
import JSON

const HERE = @__DIR__
const INPUT = joinpath(HERE, "input.esm")
const EXPECTED = joinpath(HERE, "expected.esm")

# --- Canonical whole-document JSON emitter (mirrors walk_esd_tests.jl) ---
function _canon_emit_string(io::IO, s::String)
    print(io, '"')
    for c in s
        cu = UInt32(c)
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif cu == 0x08
            print(io, "\\b")
        elseif cu == 0x09
            print(io, "\\t")
        elseif cu == 0x0A
            print(io, "\\n")
        elseif cu == 0x0C
            print(io, "\\f")
        elseif cu == 0x0D
            print(io, "\\r")
        elseif cu < 0x20
            print(io, "\\u", string(cu; base = 16, pad = 4))
        else
            print(io, c)
        end
    end
    return print(io, '"')
end
_canon_emit(io::IO, s::AbstractString) = _canon_emit_string(io, String(s))
_canon_emit(io::IO, b::Bool) = print(io, b ? "true" : "false")
_canon_emit(io::IO, ::Nothing) = print(io, "null")
_canon_emit(io::IO, n::Integer) = print(io, string(n))
_canon_emit(io::IO, n::AbstractFloat) =
    print(io, EarthSciSerialization.format_canonical_float(Float64(n)))
_canon_emit(io::IO, t::Tuple) = _canon_emit(io, collect(t))
function _canon_emit(io::IO, x::AbstractDict)
    print(io, "{")
    ks = sort!(String[string(k) for k in keys(x)])
    for (i, k) in enumerate(ks)
        i > 1 && print(io, ",")
        _canon_emit_string(io, k)
        print(io, ":")
        _canon_emit(io, x[k])
    end
    return print(io, "}")
end
function _canon_emit(io::IO, xs::AbstractVector)
    print(io, "[")
    for (i, v) in enumerate(xs)
        i > 1 && print(io, ",")
        _canon_emit(io, v)
    end
    return print(io, "]")
end
canonical_doc_json(doc) = (io = IOBuffer(); _canon_emit(io, doc); String(take!(io)))

function main()
    parsed = JSON.parse(read(INPUT, String))
    out = EarthSciSerialization.discretize(parsed)
    actual = canonical_doc_json(out)
    open(EXPECTED, "w") do io
        println(io, actual)
    end
    return println("wrote ", EXPECTED, " (", length(actual), " bytes)")
end

main()
