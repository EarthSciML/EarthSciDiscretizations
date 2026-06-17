@testsnippet DAEConformanceSetup begin
    import EarthSciSerialization
    using JSON
end

@testitem "DAE binding contract: mixed_dae_observed (dae_support=true)" setup = [DAEConformanceSetup] tags = [:conformance, :dae, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization", "dae_missing")
    spec = JSON.parsefile(joinpath(HARNESS, "mixed_dae_observed.json"))

    expect_enabled = only(filter(e -> e["mode"] == "dae_support_enabled", spec["expect"]))
    @assert expect_enabled["kind"] == "output"

    input = spec["input"]
    out = EarthSciSerialization.discretize(input; dae_support = true)

    @test out["metadata"]["system_class"] == expect_enabled["system_class"]
    info = out["metadata"]["dae_info"]
    @test info["algebraic_equation_count"] == expect_enabled["algebraic_equation_count"]
    for (mname, count) in pairs(expect_enabled["per_model"])
        @test info["per_model"][mname] == count
    end
end

@testitem "DAE binding contract: mixed_dae_observed (dae_support=false)" setup = [DAEConformanceSetup] tags = [:conformance, :dae, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization", "dae_missing")
    spec = JSON.parsefile(joinpath(HARNESS, "mixed_dae_observed.json"))

    expect_disabled = only(filter(e -> e["mode"] == "dae_support_disabled", spec["expect"]))
    @assert expect_disabled["kind"] == "error"

    input = spec["input"]
    want_code = expect_disabled["code"]

    err = try
        EarthSciSerialization.discretize(input; dae_support = false)
        nothing
    catch e
        e
    end

    @test err isa EarthSciSerialization.RuleEngineError
    @test err !== nothing && err.code == want_code
    for needle in get(expect_disabled, "message_must_contain", String[])
        @test occursin(needle, err.message)
    end
end

@testitem "DAE binding contract: pure_ode_baseline (dae_support=true)" setup = [DAEConformanceSetup] tags = [:conformance, :dae, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization", "dae_missing")
    spec = JSON.parsefile(joinpath(HARNESS, "pure_ode_baseline.json"))

    expect_enabled = only(filter(e -> e["mode"] == "dae_support_enabled", spec["expect"]))
    @assert expect_enabled["kind"] == "output"

    input = spec["input"]
    out = EarthSciSerialization.discretize(input; dae_support = true)

    @test out["metadata"]["system_class"] == expect_enabled["system_class"]
    info = out["metadata"]["dae_info"]
    @test info["algebraic_equation_count"] == expect_enabled["algebraic_equation_count"]
    for (mname, count) in pairs(expect_enabled["per_model"])
        @test info["per_model"][mname] == count
    end
end

@testitem "DAE binding contract: pure_ode_baseline (dae_support=false)" setup = [DAEConformanceSetup] tags = [:conformance, :dae, :discretization] begin
    HARNESS = joinpath(@__DIR__, "..", "tests", "conformance", "discretization", "dae_missing")
    spec = JSON.parsefile(joinpath(HARNESS, "pure_ode_baseline.json"))

    expect_disabled = only(filter(e -> e["mode"] == "dae_support_disabled", spec["expect"]))
    @assert expect_disabled["kind"] == "output"

    input = spec["input"]
    out = EarthSciSerialization.discretize(input; dae_support = false)

    @test out["metadata"]["system_class"] == expect_disabled["system_class"]
    info = out["metadata"]["dae_info"]
    @test info["algebraic_equation_count"] == expect_disabled["algebraic_equation_count"]
    for (mname, count) in pairs(expect_disabled["per_model"])
        @test info["per_model"][mname] == count
    end
end
