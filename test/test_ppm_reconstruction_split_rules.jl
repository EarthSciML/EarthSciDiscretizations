using Test
using TestItems

# Tests for the single-output sibling rules of the multi-output
# `ppm_reconstruction` umbrella (dsc-x1ry). Each sibling encodes one
# named edge of CW84 eq. (1.6) directly in ESS replacement form so it
# can be consumed by `rule_engine.jl::parse_rule` without multi-output
# dispatch (that would require an ESS rule-engine extension).

@testitem "ppm_reconstruction_right_edge rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction_right_edge", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"q_right_edge\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"O(dx^4)\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"index\"", content)
    # Replacement form, no stencil-form fields.
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"selector\"", content)
    @test !occursin("\"offset\"", content)
end

@testitem "ppm_reconstruction_left_edge rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction_left_edge", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"q_left_edge\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"O(dx^4)\"", content)
    @test occursin("\"replacement\"", content)
    @test occursin("\"arrayop\"", content)
    @test occursin("\"output_idx\"", content)
    @test occursin("\"index\"", content)
    @test !occursin("\"stencil\"", content)
    @test !occursin("\"selector\"", content)
    @test !occursin("\"offset\"", content)
end

@testitem "ppm_reconstruction_right_edge JSON parses and pins CW84 eq. (1.6) shape" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction_right_edge", rules)
    @test idx !== nothing
    rule = rules[idx]

    parsed = JSON.parse(read(rule.path, String))
    @test parsed isa AbstractDict
    spec = parsed["discretizations"]["ppm_reconstruction_right_edge"]
    @test spec["applies_to"]["op"] == "q_right_edge"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["applies_to"]["args"] == ["\$q"]
    @test spec["grid_family"] == "cartesian"
    @test spec["accuracy"] == "O(dx^4)"

    # Replacement form — single arrayop with output_idx [$x] and args [$q].
    rep = spec["replacement"]
    @test rep["op"] == "arrayop"
    @test rep["output_idx"] == ["\$x"]
    @test rep["args"] == ["\$q"]

    # Outer expr is /(plus_terms, 12). Pin denominator to 12 (CW84 eq. 1.6).
    outer = rep["expr"]
    @test outer["op"] == "/"
    @test outer["args"][2] == 12

    # Inner +-tree has 4 multiplicative terms, coefficients {-1, 7, 7, -1}.
    plus = outer["args"][1]
    @test plus["op"] == "+"
    @test length(plus["args"]) == 4
    coeffs = [t["args"][1] for t in plus["args"]]
    @test coeffs == [-1, 7, 7, -1]

    # Stencil offsets for q_right_edge: (-1, 0, +1, +2) on `$x`.
    function offset_of(term)
        idx = term["args"][2]   # index() second arg = axis pattern w/ optional ±k
        if idx isa AbstractDict
            op = idx["op"]
            k = idx["args"][2]
            return op == "+" ? Int(k) : -Int(k)
        else
            return 0    # bare "$x"
        end
    end
    offsets = [offset_of(t["args"][2]) for t in plus["args"]]
    @test offsets == [-1, 0, 1, 2]
end

@testitem "ppm_reconstruction_left_edge JSON parses and pins CW84 eq. (1.6) shifted shape" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "ppm_reconstruction_left_edge", rules)
    @test idx !== nothing
    rule = rules[idx]

    parsed = JSON.parse(read(rule.path, String))
    spec = parsed["discretizations"]["ppm_reconstruction_left_edge"]
    @test spec["applies_to"]["op"] == "q_left_edge"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["applies_to"]["args"] == ["\$q"]
    @test spec["grid_family"] == "cartesian"
    @test spec["accuracy"] == "O(dx^4)"

    rep = spec["replacement"]
    @test rep["op"] == "arrayop"
    @test rep["output_idx"] == ["\$x"]
    @test rep["args"] == ["\$q"]

    outer = rep["expr"]
    @test outer["op"] == "/"
    @test outer["args"][2] == 12

    plus = outer["args"][1]
    @test plus["op"] == "+"
    @test length(plus["args"]) == 4
    coeffs = [t["args"][1] for t in plus["args"]]
    @test coeffs == [-1, 7, 7, -1]

    function offset_of(term)
        idx = term["args"][2]
        if idx isa AbstractDict
            op = idx["op"]
            k = idx["args"][2]
            return op == "+" ? Int(k) : -Int(k)
        else
            return 0
        end
    end
    offsets = [offset_of(t["args"][2]) for t in plus["args"]]
    # Left edge q_{i-1/2}: same coeffs shifted one cell left of q_right_edge.
    @test offsets == [-2, -1, 0, 1]
end

@testitem "ppm_reconstruction edge siblings JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    for name in ("ppm_reconstruction_left_edge", "ppm_reconstruction_right_edge")
        idx = findfirst(r -> r.name == name, rules)
        @test idx !== nothing
        rule = rules[idx]
        parsed = JSON.parse(read(rule.path, String))
        reparsed = JSON.parse(JSON.json(parsed))
        @test reparsed == parsed
    end
end
