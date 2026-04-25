module WalkESDTests

using EarthSciDiscretizations: load_rules, RuleFile, eval_coeff
using JSON

export walk_esd_tests,
    discover_rules,
    run_layer_a,
    run_layer_b,
    run_layer_c,
    RuleResult,
    LayerOutcome,
    write_junit

@enum LayerOutcome LAYER_PASS LAYER_FAIL LAYER_SKIP

struct LayerResult
    outcome::LayerOutcome
    reason::String
end

struct RuleResult
    family::Symbol
    name::String
    path::String
    layer_a::LayerResult
    layer_b::LayerResult
    layer_c::LayerResult
end

"""
    discover_rules(catalog) -> Vector{RuleFile}

Return the list of rule files under `catalog`. A non-existent or empty catalog
returns an empty vector (no error), so a repo with no rules yields a green
build with zero walker tests.
"""
function discover_rules(catalog::AbstractString)
    isdir(catalog) || return RuleFile[]
    return load_rules(catalog)
end

"""
    rule_fixtures_dir(rule) -> String

Fixtures live at `discretizations/<family>/<name>/fixtures/`, a sibling to
the `<name>.json` rule file. Absent until authored per rule; `isdir` on the
returned path tells the caller whether fixtures exist.
"""
function rule_fixtures_dir(rule::RuleFile)
    family_dir = dirname(rule.path)
    return joinpath(family_dir, rule.name, "fixtures")
end

"""
    run_layer_a(rule) -> LayerResult

Canonical-form conformance: run the rule engine against `fixtures/canonical/input.esm`
and byte-compare to `fixtures/canonical/expected.esm`. Until the ESS rule engine
(gt-b13f) is wired into this walker, layer A skips with a descriptive reason.
"""
function run_layer_a(rule::RuleFile)
    canonical = joinpath(rule_fixtures_dir(rule), "canonical")
    if !isdir(canonical)
        return LayerResult(LAYER_SKIP, "no canonical fixtures at $(relpath_from_repo(canonical))")
    end
    input = joinpath(canonical, "input.esm")
    expected = joinpath(canonical, "expected.esm")
    if !isfile(input) || !isfile(expected)
        return LayerResult(LAYER_FAIL, "canonical/ present but missing input.esm or expected.esm")
    end
    if !rule_engine_available()
        return LayerResult(LAYER_SKIP, "rule engine (ESS gt-b13f) not yet available to walker")
    end
    return apply_rule_and_diff(rule, input, expected)
end

"""
    run_layer_b(rule) -> LayerResult

MMS convergence: instantiate the rule on a sequence of grids, solve a manufactured
problem, measure L∞ error, verify slope matches the declared order. Requires the
tree-walk evaluator (gt-TBD); skips until available.
"""
function run_layer_b(rule::RuleFile)
    convergence = joinpath(rule_fixtures_dir(rule), "convergence")
    if !isdir(convergence)
        return LayerResult(LAYER_SKIP, "no convergence fixtures at $(relpath_from_repo(convergence))")
    end
    if !tree_walk_evaluator_available()
        return LayerResult(LAYER_SKIP, "tree-walk evaluator (gt-TBD) not yet available to walker")
    end
    return run_mms_convergence(rule, convergence)
end

"""
    run_layer_c(rule) -> LayerResult

Integration benchmarks (Williamson 1992, DCMIP, etc.). Skipped by default; enabled
by `ESD_RUN_INTEGRATION=1`.
"""
function run_layer_c(rule::RuleFile)
    if get(ENV, "ESD_RUN_INTEGRATION", "0") != "1"
        return LayerResult(LAYER_SKIP, "integration-only; set ESD_RUN_INTEGRATION=1 to run")
    end
    integration = joinpath(rule_fixtures_dir(rule), "integration")
    if !isdir(integration)
        return LayerResult(LAYER_SKIP, "no integration fixtures declared")
    end
    if !tree_walk_evaluator_available()
        return LayerResult(LAYER_SKIP, "tree-walk evaluator (gt-TBD) not yet available to walker")
    end
    return run_integration_benchmarks(rule, integration)
end

rule_engine_available() = false
# Layer B routes coefficient evaluation through EarthSciSerialization's tree-walk
# evaluator (see src/rule_eval.jl). Layer A still needs the ESS rule engine
# (canonical-form round-trip) which is a larger surface than the evaluator.
tree_walk_evaluator_available() = true

function apply_rule_and_diff(::RuleFile, ::AbstractString, ::AbstractString)
    return LayerResult(LAYER_SKIP, "layer A executor stub — rule engine not yet wired")
end

"""
    run_mms_convergence(rule, convergence_dir) -> LayerResult

Dispatch MMS convergence on a rule using the ESS AST evaluator for all
rule-side coefficient computations. Rule-specific appliers (stencil
application, BC handling, manufactured-solution sampling) live in small
per-rule drivers registered below. Rules without a registered driver skip
with a descriptive reason — the walker keeps green.
"""
function run_mms_convergence(rule::RuleFile, convergence_dir::AbstractString)
    input_path = joinpath(convergence_dir, "input.esm")
    expected_path = joinpath(convergence_dir, "expected.esm")
    if !(isfile(input_path) && isfile(expected_path))
        return LayerResult(LAYER_FAIL, "convergence/ present but missing input.esm or expected.esm")
    end
    input = JSON.parse(read(input_path, String))
    expected = JSON.parse(read(expected_path, String))

    if rule.family === :finite_difference && rule.name == "centered_2nd_uniform"
        return _mms_centered_2nd_uniform(rule, input, expected)
    end
    return LayerResult(
        LAYER_SKIP,
        "no MMS driver registered for $(rule.family)/$(rule.name); fixture present but walker cannot apply it"
    )
end

# Apply the centered-2nd finite-difference stencil to a smooth periodic
# manufactured solution on a sequence of uniform grids. All stencil
# coefficients come from the rule JSON via the ESS evaluator; the walker
# does not reimplement `-1/(2*dx)` etc. in Julia.
function _mms_centered_2nd_uniform(rule::RuleFile, input, expected)
    raw = JSON.parse(read(rule.path, String))
    spec = raw["discretizations"]["centered_2nd_uniform"]
    stencil = spec["stencil"]
    grids = [Int(g["n"]) for g in input["grids"]]
    min_order = Float64(expected["expected_min_order"])

    # Manufactured: u(x) = sin(2π x) on [0, 1] periodic; du/dx = 2π cos(2π x).
    # Cell-center sampling: x_i = (i − 1/2) * dx, i = 1..n.
    errs = Float64[]
    for n in grids
        dx = 1.0 / n
        bindings = Dict("dx" => dx)
        coeff_pairs = [
            (Int(s["selector"]["offset"]), eval_coeff(s["coeff"], bindings))
                for s in stencil
        ]
        u = [sin(2π * ((i - 0.5) * dx)) for i in 1:n]
        du_num = zeros(n)
        for i in 1:n
            acc = 0.0
            for (off, c) in coeff_pairs
                j = mod1(i + off, n)  # periodic wrap
                acc += c * u[j]
            end
            du_num[i] = acc
        end
        du_exact = [2π * cos(2π * ((i - 0.5) * dx)) for i in 1:n]
        push!(errs, maximum(abs.(du_num .- du_exact)))
    end

    if any(!isfinite, errs) || any(e -> e <= 0, errs)
        return LayerResult(LAYER_FAIL, "non-finite or zero error on some grid; errs=$(errs)")
    end
    orders = [log2(errs[i] / errs[i + 1]) for i in 1:(length(errs) - 1)]
    observed = minimum(orders)
    if observed < min_order
        return LayerResult(
            LAYER_FAIL,
            "observed min order $(round(observed; digits = 3)) below expected $(min_order); errs=$(errs)"
        )
    end
    return LayerResult(
        LAYER_PASS,
        "min order $(round(observed; digits = 3)) >= $(min_order) on grids $(grids)"
    )
end

function run_integration_benchmarks(::RuleFile, ::AbstractString)
    return LayerResult(LAYER_SKIP, "layer C executor stub — integration harness not yet wired")
end

function relpath_from_repo(path::AbstractString)
    try
        return relpath(path)
    catch
        return path
    end
end

"""
    walk_esd_tests(; catalog, junit_path=nothing) -> Vector{RuleResult}

Main entry point. Discovers rules under `catalog`, runs layers A/B/C per rule,
prints a report to stdout, and optionally writes JUnit XML.
"""
function walk_esd_tests(;
        catalog::AbstractString,
        junit_path::Union{Nothing, AbstractString} = nothing,
        io::IO = stdout,
    )
    rules = discover_rules(catalog)
    results = Vector{RuleResult}()
    for rule in rules
        a = run_layer_a(rule)
        b = run_layer_b(rule)
        c = run_layer_c(rule)
        push!(results, RuleResult(rule.family, rule.name, rule.path, a, b, c))
    end
    print_report(io, catalog, results)
    if junit_path !== nothing
        mkpath(dirname(junit_path))
        write_junit(junit_path, results)
    end
    return results
end

function outcome_tag(o::LayerOutcome)
    o == LAYER_PASS && return "PASS"
    o == LAYER_FAIL && return "FAIL"
    return "SKIP"
end

function print_report(io::IO, catalog, results::Vector{RuleResult})
    println(io, "=== ESD Walker: $(catalog) ===")
    if isempty(results)
        println(io, "(no rules discovered; empty catalog treated as green)")
        return
    end
    for r in results
        println(io, "[$(r.family)/$(r.name)]  A=$(outcome_tag(r.layer_a.outcome))  B=$(outcome_tag(r.layer_b.outcome))  C=$(outcome_tag(r.layer_c.outcome))")
        for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("C", r.layer_c))
            if lr.outcome != LAYER_PASS && !isempty(lr.reason)
                println(io, "    layer $tag: $(lr.reason)")
            end
        end
    end
    total = length(results) * 3
    npass = sum(count_outcome(r, LAYER_PASS) for r in results; init = 0)
    nfail = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    nskip = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    return println(io, "summary: $total cases  pass=$npass  fail=$nfail  skip=$nskip")
end

function count_outcome(r::RuleResult, target::LayerOutcome)
    return (r.layer_a.outcome == target) + (r.layer_b.outcome == target) + (r.layer_c.outcome == target)
end

function xml_escape(s::AbstractString)
    buf = replace(s, "&" => "&amp;")
    buf = replace(buf, "<" => "&lt;")
    buf = replace(buf, ">" => "&gt;")
    buf = replace(buf, "\"" => "&quot;")
    buf = replace(buf, "'" => "&apos;")
    return buf
end

"""
    write_junit(path, results)

Emit a JUnit XML summary to `path`. One `<testcase>` per (rule, layer). Skipped
layers carry `<skipped>` with the reason; failed layers carry `<failure>`.
"""
function write_junit(path::AbstractString, results::Vector{RuleResult})
    total = length(results) * 3
    failures = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    skipped = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io, "<testsuites>")
        println(io, "  <testsuite name=\"ESD Walker\" tests=\"$total\" failures=\"$failures\" skipped=\"$skipped\">")
        for r in results
            for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("C", r.layer_c))
                classname = xml_escape("$(r.family).$(r.name)")
                name = "layer_$tag"
                print(io, "    <testcase classname=\"", classname, "\" name=\"", name, "\">")
                if lr.outcome == LAYER_SKIP
                    print(io, "<skipped message=\"", xml_escape(lr.reason), "\"/>")
                elseif lr.outcome == LAYER_FAIL
                    print(io, "<failure message=\"", xml_escape(lr.reason), "\"/>")
                end
                println(io, "</testcase>")
            end
        end
        println(io, "  </testsuite>")
        println(io, "</testsuites>")
    end
    return path
end

end # module
