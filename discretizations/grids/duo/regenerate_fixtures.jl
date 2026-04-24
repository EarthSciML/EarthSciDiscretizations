#!/usr/bin/env julia
# Regenerate canonical DUO `.esm` fixtures from the Julia reference binding.
#
# Per the 2026-04-20 mayor correction on bead `dsc-8qn`, a `.esm` grid entry
# is a small declarative config (family + options + loader ref + dims) —
# not a serialized geometry blob. The per-binding runtime derives geometry
# from the declaration via accessors.
#
# The `provenance` block is stripped so the payload is binding-neutral:
# `provenance.binding` differs between Julia / Python / Rust / TypeScript,
# which would prevent byte-identity across bindings. What remains is the
# declarative §6 schema payload that cross-binding conformance gates on.
#
# Run from the repo root:
#     julia discretizations/grids/duo/regenerate_fixtures.jl
#
# Activates a temporary env so it runs cleanly without touching the main
# Project.toml deps; JSON is pulled in just for IO.

using Pkg
let env_dir = mktempdir()
    Pkg.activate(env_dir; io = devnull)
    Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "..", "..")); io = devnull)
    Pkg.add("JSON"; io = devnull)
end

using EarthSciDiscretizations
using EarthSciDiscretizations: to_esm
using JSON

const REPO_ROOT = dirname(dirname(dirname(dirname(@__FILE__))))
const OUT_DIR = joinpath(REPO_ROOT, "discretizations", "grids", "duo")

# Canonical ladder: three subdivision levels at Earth radius.
# Level r yields 20·4^r triangular cells (20 / 80 / 320).
const FIXTURES = [
    ("icos_level0", 0),
    ("icos_level1", 1),
    ("icos_level2", 2),
]

function canonical_bytes(d::AbstractDict)
    io = IOBuffer()
    JSON.print(io, d, 2)
    print(io, "\n")
    return String(take!(io))
end

function main()
    for (name, level) in FIXTURES
        loader = (
            path = "builtin://icosahedral/$(level)",
            reader = "builtin_icosahedral",
            check = "strict",
        )
        g = EarthSciDiscretizations.grids.duo(;
            loader = loader, R = 6.371e6, dtype = Float64, ghosts = 0,
        )
        d = to_esm(g)
        delete!(d, "provenance")  # strip binding/platform-specific block
        path = joinpath(OUT_DIR, "$(name).esm")
        write(path, canonical_bytes(d))
        println("wrote $(path) (level=$(level), n_cells=$(d["n_cells"]))")
    end
end

main()
