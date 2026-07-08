#!/usr/bin/env julia
# Regenerate a Layer-A canonical `expected.esm` golden from its sibling
# `input.esm`, byte-identical to what the ESD walker's `apply_rule_and_diff`
# compares against: `canonical_doc_json(EarthSciSerialization.discretize(input))`.
#
# Usage (from repo root, env instantiated):
#   julia --project=. tools/gen_canonical_golden.jl <path/to/canonical/input.esm> [more inputs...]
#
# Writes `expected.esm` (no trailing newline beyond the canonical form) next to
# each input. The walker tolerates a single trailing newline, which we add to
# match the committed style of the other goldens.
using EarthSciDiscretizations
import EarthSciSerialization
include(joinpath(@__DIR__, "..", "test", "walk_esd_tests.jl"))
using JSON

for input_path in ARGS
    parsed = JSON.parse(read(input_path, String))
    out = EarthSciSerialization.discretize(parsed)
    canon = WalkESDTests.canonical_doc_json(out)
    expected_path = joinpath(dirname(input_path), "expected.esm")
    open(expected_path, "w") do io
        write(io, canon)
        write(io, "\n")
    end
    println("wrote $(expected_path) ($(length(canon)) bytes)")
end
