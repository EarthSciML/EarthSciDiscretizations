module WalkESDTests

using EarthSciDiscretizations: load_rules, RuleFile
using EarthSciSerialization: verify_mms_convergence, MMSEvaluatorError
import EarthSciSerialization
using JSON

export walk_esd_tests,
    discover_rules,
    run_layer_a,
    run_layer_b,
    run_layer_c,
    run_layer_limiter,
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
    layer_limiter::LayerResult
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

Rule-engine round-trip against fixture inputs. Two variants are supported:

- `fixtures/canonical/` — whole-document conformance via
  `EarthSciSerialization.discretize`, byte-compared to the canonical-form
  rendering of `expected.esm`. Used by stencil-coefficient rules.
- `fixtures/rewrite/` — expression-level rewrite via
  `EarthSciSerialization.rewrite`, byte-compared to the canonical JSON of
  `expected.esm`. Used by index-rewrite (boundary-condition) rules whose
  acceptance signature is "this input expression rewrites to that output
  expression" rather than a numeric residual.

If both directories exist, both run and outcomes are AND-combined (any FAIL
wins; otherwise PASS unless all variants SKIP). If neither exists, SKIP.
"""
function run_layer_a(rule::RuleFile)
    fdir = rule_fixtures_dir(rule)
    canonical_dir = joinpath(fdir, "canonical")
    rewrite_dir = joinpath(fdir, "rewrite")
    have_canonical = isdir(canonical_dir)
    have_rewrite = isdir(rewrite_dir)
    if !have_canonical && !have_rewrite
        return LayerResult(LAYER_SKIP, "no canonical or rewrite fixtures at $(relpath_from_repo(fdir))")
    end
    results = LayerResult[]
    have_canonical && push!(results, _run_canonical_variant(rule, canonical_dir))
    have_rewrite && push!(results, _run_rewrite_variant(rule, rewrite_dir))
    return _combine_layer_a(results)
end

function _run_canonical_variant(rule::RuleFile, canonical::AbstractString)
    input = joinpath(canonical, "input.esm")
    expected = joinpath(canonical, "expected.esm")
    if !isfile(input) || !isfile(expected)
        return LayerResult(LAYER_FAIL, "canonical/ present but missing input.esm or expected.esm")
    end
    return apply_rule_and_diff(rule, input, expected)
end

function _run_rewrite_variant(rule::RuleFile, rewrite_dir::AbstractString)
    input = joinpath(rewrite_dir, "input.esm")
    expected = joinpath(rewrite_dir, "expected.esm")
    if !isfile(input) || !isfile(expected)
        return LayerResult(LAYER_FAIL, "rewrite/ present but missing input.esm or expected.esm")
    end
    return apply_rewrite_and_diff(rule, input, expected)
end

# AND-combine: any FAIL wins (failures dominate); else PASS if any variant
# passed; else SKIP. Reasons are joined with "; " so the report shows both.
function _combine_layer_a(results::Vector{LayerResult})
    length(results) == 1 && return results[1]
    any(r -> r.outcome == LAYER_FAIL, results) &&
        return LayerResult(LAYER_FAIL, join((r.reason for r in results if r.outcome == LAYER_FAIL), "; "))
    any(r -> r.outcome == LAYER_PASS, results) &&
        return LayerResult(LAYER_PASS, join((r.reason for r in results), "; "))
    return LayerResult(LAYER_SKIP, join((r.reason for r in results), "; "))
end

"""
    run_layer_b(rule) -> LayerResult

MMS convergence: dispatch to ESS's `verify_mms_convergence`, which evaluates
all stencil coefficients via the AST evaluator and runs the manufactured-solution
sweep on a sequence of grids.
"""
function run_layer_b(rule::RuleFile)
    convergence = joinpath(rule_fixtures_dir(rule), "convergence")
    if !isdir(convergence)
        return LayerResult(LAYER_SKIP, "no convergence fixtures at $(relpath_from_repo(convergence))")
    end
    return run_mms_convergence(rule, convergence)
end

"""
    run_layer_limiter(rule) -> LayerResult

Layer B' — monotonicity / TVD acceptance for slope-ratio limiter rules
(Sweby (1984), Roe (1986)). Limiter rules ship a separate fixture kind under
`fixtures/monotonicity/` because their acceptance is *monotonicity preservation*
on the Sweby region plus *strict TVD on smooth+square-wave ICs over an advection
period* — not asymptotic convergence-order on a manufactured solution. See bead
dsc-8vu and `discretizations/finite_volume/flux_limiter_*.json` for the rule
form.

Skips if `monotonicity/` is absent (the common case for non-limiter rules);
fails if either fixture file is missing or any check fails; passes when all
reference (r, phi) pairs match within tolerance, every Sweby property holds
across the sweep, and TV(q_final) <= TV(q_initial) + tvd_tolerance after one
period.
"""
function run_layer_limiter(rule::RuleFile)
    monotonicity = joinpath(rule_fixtures_dir(rule), "monotonicity")
    if !isdir(monotonicity)
        return LayerResult(LAYER_SKIP, "no monotonicity fixtures at $(relpath_from_repo(monotonicity))")
    end
    return run_monotonicity_check(rule, monotonicity)
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
    return run_integration_benchmarks(rule, integration)
end

rule_engine_available() = true
# Layer B routes coefficient evaluation through EarthSciSerialization's tree-walk
# evaluator (see src/rule_eval.jl). Layer A drives the ESS rule engine end-to-end
# via `EarthSciSerialization.discretize` and a canonical whole-document JSON
# emitter that matches the ESS conformance harness contract.
tree_walk_evaluator_available() = true

"""
    apply_rule_and_diff(rule, input_path, expected_path) -> LayerResult

Read `input_path` (an `.esm` JSON document), run the ESS rule engine on it via
`EarthSciSerialization.discretize`, emit the canonical whole-document JSON form
(sorted keys, `format_canonical_float` for floats, minified per the ESS
conformance contract — see `tests/conformance/discretize/README.md` in ESS),
and byte-compare to the (canonical) contents of `expected_path`. Trailing
newlines on the expected file are tolerated. Returns `LAYER_PASS` on byte
equality, `LAYER_FAIL` with a small diff window on mismatch.
"""
function apply_rule_and_diff(
        ::RuleFile, input_path::AbstractString,
        expected_path::AbstractString
    )
    parsed = try
        JSON.parse(read(input_path, String))
    catch err
        return LayerResult(
            LAYER_FAIL,
            "failed to parse $(relpath_from_repo(input_path)): $(err)"
        )
    end
    out = try
        EarthSciSerialization.discretize(parsed)
    catch err
        io = IOBuffer()
        showerror(io, err)
        return LayerResult(
            LAYER_FAIL,
            "discretize threw on $(relpath_from_repo(input_path)): $(String(take!(io)))"
        )
    end
    actual = canonical_doc_json(out)
    expected = rstrip(read(expected_path, String), '\n')
    if actual == expected
        return LayerResult(
            LAYER_PASS,
            "canonical-form match ($(length(actual)) bytes)"
        )
    end
    return LayerResult(LAYER_FAIL, _byte_diff_message(actual, expected))
end

"""
    apply_rewrite_and_diff(rule, input_path, expected_path) -> LayerResult

Read `input_path` (a JSON document containing `context` and `expression`),
parse the rule via `EarthSciSerialization.parse_rules`, build a
`RuleContext` from the fixture's grids/variables tables, run
`EarthSciSerialization.rewrite` on the input expression, and emit the
result as canonical-form JSON via `EarthSciSerialization.canonical_json`.
Byte-compare to the contents of `expected_path` (trailing newlines tolerated).

This is the rule-engine path for index-rewrite rules (e.g. `periodic_bc`)
whose acceptance signature is "this expression rewrites to that one" rather
than a numeric residual on a stencil sweep. See
`discretizations/finite_difference/periodic_bc/fixtures/rewrite/` for the
seed fixture.
"""
function apply_rewrite_and_diff(
        rule::RuleFile, input_path::AbstractString,
        expected_path::AbstractString
    )
    input_doc = try
        JSON.parse(read(input_path, String))
    catch err
        return LayerResult(
            LAYER_FAIL,
            "failed to parse $(relpath_from_repo(input_path)): $(err)"
        )
    end
    expr_json = get(input_doc, "expression", nothing)
    if expr_json === nothing
        return LayerResult(
            LAYER_FAIL,
            "$(relpath_from_repo(input_path)) missing required `expression` field"
        )
    end
    rule_doc = try
        JSON.parse(read(rule.path, String))
    catch err
        return LayerResult(
            LAYER_FAIL,
            "failed to parse $(relpath_from_repo(rule.path)): $(err)"
        )
    end
    rules_obj = get(rule_doc, "rules", nothing)
    if rules_obj === nothing
        return LayerResult(
            LAYER_FAIL,
            "$(relpath_from_repo(rule.path)) missing top-level `rules` table"
        )
    end
    rules = try
        EarthSciSerialization.parse_rules(rules_obj)
    catch err
        return LayerResult(
            LAYER_FAIL,
            "parse_rules threw on $(relpath_from_repo(rule.path)): $(sprint(showerror, err))"
        )
    end
    ctx = _build_rule_context(get(input_doc, "context", Dict{String, Any}()))
    input_expr = try
        _expr_from_json(expr_json)
    catch err
        return LayerResult(
            LAYER_FAIL,
            "failed to lift `expression` from $(relpath_from_repo(input_path)): $(sprint(showerror, err))"
        )
    end
    out_expr = try
        EarthSciSerialization.rewrite(input_expr, rules, ctx)
    catch err
        return LayerResult(
            LAYER_FAIL,
            "rewrite threw on $(relpath_from_repo(input_path)): $(sprint(showerror, err))"
        )
    end
    actual = EarthSciSerialization.canonical_json(out_expr)
    expected = rstrip(read(expected_path, String), '\n')
    if actual == expected
        return LayerResult(
            LAYER_PASS,
            "rewrite canonical-form match ($(length(actual)) bytes)"
        )
    end
    return LayerResult(LAYER_FAIL, _byte_diff_message(actual, expected))
end

# Lift a JSON-decoded subtree (Dict/Number/String) to the corresponding
# `EarthSciSerialization.Expr` node. Mirrors ESS's private `_parse_expr`,
# but only the fields used by index-rewrite fixtures (`op`, `args`, plus
# `wrt`/`dim` for symmetry with the rule engine's pattern shape) — fancy
# `arrayop`/`makearray` fields are out of scope until a fixture needs them.
function _expr_from_json(v)
    if v isa Bool
        # Bool <: Integer in Julia; reject explicitly so a stray `true` in a
        # fixture surfaces as a parse error rather than silently becoming `1`.
        throw(ArgumentError("boolean literal not valid in expression position"))
    elseif v isa Integer
        return EarthSciSerialization.IntExpr(Int64(v))
    elseif v isa AbstractFloat
        return EarthSciSerialization.NumExpr(Float64(v))
    elseif v isa AbstractString
        return EarthSciSerialization.VarExpr(String(v))
    elseif v isa AbstractDict
        op = String(v["op"])
        args_raw = get(v, "args", Any[])
        args = EarthSciSerialization.Expr[_expr_from_json(a) for a in args_raw]
        wrt = haskey(v, "wrt") && v["wrt"] !== nothing ? String(v["wrt"]) : nothing
        dim = haskey(v, "dim") && v["dim"] !== nothing ? String(v["dim"]) : nothing
        return EarthSciSerialization.OpExpr(op, args; wrt = wrt, dim = dim)
    end
    throw(ArgumentError("cannot parse expression node of type $(typeof(v))"))
end

# Build a RuleContext from the `context` block of a rewrite fixture. Both
# subtables default to empty so a fixture that needs no metadata (a pure
# syntactic rewrite, no guards) parses cleanly.
function _build_rule_context(ctx_json)
    grids = Dict{String, Dict{String, Any}}()
    variables = Dict{String, Dict{String, Any}}()
    grids_raw = get(ctx_json, "grids", Dict{String, Any}())
    for (k, v) in grids_raw
        grids[String(k)] = Dict{String, Any}(String(kk) => vv for (kk, vv) in v)
    end
    vars_raw = get(ctx_json, "variables", Dict{String, Any}())
    for (k, v) in vars_raw
        variables[String(k)] = Dict{String, Any}(String(kk) => vv for (kk, vv) in v)
    end
    return EarthSciSerialization.RuleContext(grids, variables)
end

# Show the first divergence with a small surrounding window, so debugging a
# canonical-form drift does not require diffing two minified strings by eye.
function _byte_diff_message(actual::AbstractString, expected::AbstractString)
    n = min(length(actual), length(expected))
    diff_at = n + 1
    for i in 1:n
        if actual[i] != expected[i]
            diff_at = i
            break
        end
    end
    window = 40
    lo = max(1, diff_at - window)
    hi_a = min(length(actual), diff_at + window)
    hi_e = min(length(expected), diff_at + window)
    return string(
        "canonical-form mismatch at byte ", diff_at,
        " (actual=", length(actual), "B, expected=", length(expected), "B)\n",
        "  actual:   …", actual[lo:hi_a], "…\n",
        "  expected: …", expected[lo:hi_e], "…",
    )
end

# ---------------------------------------------------------------------------
# Canonical whole-document JSON emitter — mirrors the ESS conformance harness
# (`tests/conformance/discretize/conformance_discretize_test.jl`): sorted keys,
# minified, RFC §5.4.6 float formatting via `format_canonical_float`. Kept
# in-tree (not exported from ESS) because ESS itself uses it only for tests.
# ---------------------------------------------------------------------------

function canonical_doc_json(doc)::String
    io = IOBuffer()
    _canon_emit(io, doc)
    return String(take!(io))
end

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

_canon_emit(io::IO, t::Tuple) = _canon_emit(io, collect(t))
_canon_emit(io::IO, s::AbstractString) = _canon_emit_string(io, String(s))
_canon_emit(io::IO, b::Bool) = print(io, b ? "true" : "false")
_canon_emit(io::IO, ::Nothing) = print(io, "null")
_canon_emit(io::IO, n::Integer) = print(io, string(n))
_canon_emit(io::IO, n::AbstractFloat) =
    print(io, EarthSciSerialization.format_canonical_float(Float64(n)))

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

"""
    run_mms_convergence(rule, convergence_dir) -> LayerResult

End-to-end Layer B: parse the rule JSON, the input fixture and the expected
fixture, then route through `EarthSciSerialization.verify_mms_convergence`.
ESS owns the manufactured-solution sweep, the AST coefficient evaluator,
and the order-deficit check. Rules whose manufactured solution is not
registered with ESS skip with a descriptive reason.
"""
function run_mms_convergence(rule::RuleFile, convergence_dir::AbstractString)
    input_path = joinpath(convergence_dir, "input.esm")
    expected_path = joinpath(convergence_dir, "expected.esm")
    if !(isfile(input_path) && isfile(expected_path))
        return LayerResult(LAYER_FAIL, "convergence/ present but missing input.esm or expected.esm")
    end
    rule_json = JSON.parse(read(rule.path, String))
    input_json = JSON.parse(read(input_path, String))
    expected_json = JSON.parse(read(expected_path, String))

    # Fixture-declared non-applicability: rules whose acceptance signature
    # isn't a manufactured-solution convergence sweep (index-rewrite rules,
    # TVD limiters, reconstruction-style rules pending ESS harness extension)
    # ship an `applicable: false` + `skip_reason` marker so the walker
    # surfaces a structured SKIP instead of a missing-fixture FAIL.
    if get(input_json, "applicable", true) === false
        reason = get(input_json, "skip_reason", "fixture declares applicable:false (no reason given)")
        return LayerResult(LAYER_SKIP, "fixture-declared not applicable: $(reason)")
    end

    try
        result = verify_mms_convergence(rule_json, input_json, expected_json)
        observed = round(result.observed_min_order; digits = 3)
        threshold = expected_json["expected_min_order"]
        grids = [Int(g["n"]) for g in input_json["grids"]]
        return LayerResult(
            LAYER_PASS,
            "min order $(observed) >= $(threshold) on grids $(grids)"
        )
    catch err
        if err isa MMSEvaluatorError
            return LayerResult(LAYER_FAIL, sprint(showerror, err))
        end
        rethrow()
    end
end

# ---------------------------------------------------------------------------
# Layer B' (limiter) — monotonicity / TVD evaluator scoped to the walker.
#
# ESS 0.0.3's `evaluate` does not yet implement the `max`/`min` ops that
# slope-ratio limiters rely on (a follow-up bead is filed against ESS), so
# the limiter runner uses a small local AST evaluator confined to this file.
# This mirrors the precedent in `test/test_flux_limiters_rule.jl` and is
# allowed by the dsc-8vu risk section. When ESS gains max/min, this can be
# replaced by `EarthSciDiscretizations.eval_coeff` with no fixture changes.
# ---------------------------------------------------------------------------

function _eval_limiter_ast(node, bindings::Dict{String, Float64})::Float64
    if node isa Number
        return Float64(node)
    elseif node isa AbstractString
        name = startswith(node, "\$") ? node[2:end] : node
        haskey(bindings, name) ||
            throw(ArgumentError("unbound variable: $(name)"))
        return bindings[name]
    elseif node isa AbstractDict
        op = node["op"]
        args = [_eval_limiter_ast(a, bindings) for a in node["args"]]
        if op == "max"
            return maximum(args)
        elseif op == "min"
            return minimum(args)
        elseif op == "+"
            return length(args) == 1 ? args[1] : sum(args)
        elseif op == "-"
            return length(args) == 1 ? -args[1] : args[1] - args[2]
        elseif op == "*"
            return length(args) == 1 ? args[1] : prod(args)
        elseif op == "/"
            return args[1] / args[2]
        elseif op == "abs"
            return abs(args[1])
        else
            throw(ArgumentError("unsupported op in limiter walker: $op"))
        end
    end
    throw(ArgumentError("unrecognized AST node type: $(typeof(node))"))
end

"""
    run_monotonicity_check(rule, monotonicity_dir) -> LayerResult

Run the limiter Sweby + TVD checks on a single rule. Reads
`monotonicity/sweby_check.esm` (phi(r) reference values + sweep bounds) and
`monotonicity/tvd_check.esm` (1D periodic advection params), evaluates the
rule formula AST against both, and returns PASS / FAIL with a one-line
summary. A missing file is FAIL (the dir was authored — we expect both
fixtures); a parse error or property violation surfaces as FAIL with the
specific failure included in the message.
"""
function run_monotonicity_check(rule::RuleFile, monotonicity_dir::AbstractString)
    sweby_path = joinpath(monotonicity_dir, "sweby_check.esm")
    tvd_path = joinpath(monotonicity_dir, "tvd_check.esm")
    if !isfile(sweby_path) || !isfile(tvd_path)
        return LayerResult(
            LAYER_FAIL,
            "monotonicity/ present but missing sweby_check.esm or tvd_check.esm",
        )
    end

    rule_json = try
        JSON.parse(read(rule.path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse rule $(relpath_from_repo(rule.path)): $(sprint(showerror, err))")
    end
    spec = get(get(rule_json, "discretizations", Dict()), rule.name, nothing)
    if !(spec isa AbstractDict) || !haskey(spec, "formula")
        return LayerResult(LAYER_FAIL, "rule $(rule.name) has no `formula` AST under discretizations.$(rule.name)")
    end
    formula = spec["formula"]

    # --- Sweby check ------------------------------------------------------
    sweby = try
        JSON.parse(read(sweby_path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse $(relpath_from_repo(sweby_path)): $(sprint(showerror, err))")
    end
    tol = Float64(get(sweby, "tolerance", 1.0e-12))
    refs = get(sweby, "reference_values", Any[])
    n_ref = 0
    for pair in refs
        r = Float64(pair["r"])
        expected = Float64(pair["phi"])
        actual = try
            _eval_limiter_ast(formula, Dict("r" => r))
        catch err
            return LayerResult(LAYER_FAIL, "AST eval failed at r=$(r): $(sprint(showerror, err))")
        end
        if !isapprox(actual, expected; atol = tol, rtol = 0.0)
            return LayerResult(
                LAYER_FAIL,
                "phi($(r)) = $(actual) but reference is $(expected) (tol=$(tol))",
            )
        end
        n_ref += 1
    end
    props = get(sweby, "tvd_properties", Dict())
    rmin = Float64(get(props, "sweep_r_min", -2.0))
    rmax = Float64(get(props, "sweep_r_max", 5.0))
    rstep = Float64(get(props, "sweep_r_step", 0.05))
    n_sweep = 0
    for r in rmin:rstep:rmax
        rf = Float64(r)
        phi = try
            _eval_limiter_ast(formula, Dict("r" => rf))
        catch err
            return LayerResult(LAYER_FAIL, "AST eval failed at r=$(rf): $(sprint(showerror, err))")
        end
        if phi < -tol
            return LayerResult(LAYER_FAIL, "positivity violated at r=$(rf): phi=$(phi)")
        end
        if phi > 2.0 + tol
            return LayerResult(LAYER_FAIL, "Sweby upper bound violated at r=$(rf): phi=$(phi) > 2")
        end
        if rf <= 0 && !isapprox(phi, 0.0; atol = tol)
            return LayerResult(LAYER_FAIL, "monotonicity-preserving violated: phi($(rf))=$(phi), expected 0")
        end
        n_sweep += 1
    end
    phi_one = _eval_limiter_ast(formula, Dict("r" => 1.0))
    if !isapprox(phi_one, 1.0; atol = tol)
        return LayerResult(LAYER_FAIL, "consistency violated: phi(1)=$(phi_one), expected 1")
    end

    # --- TVD check --------------------------------------------------------
    tvd = try
        JSON.parse(read(tvd_path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse $(relpath_from_repo(tvd_path)): $(sprint(showerror, err))")
    end
    tv_summary = try
        _run_tvd_advection(formula, tvd)
    catch err
        return LayerResult(LAYER_FAIL, "TVD check threw: $(sprint(showerror, err))")
    end
    if !tv_summary.passed
        return LayerResult(
            LAYER_FAIL,
            "TVD violated: TV_final=$(tv_summary.tv_final) > TV_initial=$(tv_summary.tv_initial) + tol=$(tv_summary.tol)",
        )
    end

    return LayerResult(
        LAYER_PASS,
        "Sweby OK ($(n_ref) refs, $(n_sweep) sweep pts); TVD OK (TV $(round(tv_summary.tv_initial; digits=6)) -> $(round(tv_summary.tv_final; digits=6)))",
    )
end

# Container for TVD-check outcomes — keeps `run_monotonicity_check` flat.
struct _TVDSummary
    passed::Bool
    tv_initial::Float64
    tv_final::Float64
    tol::Float64
end

function _run_tvd_advection(formula, tvd::AbstractDict)
    n = Int(tvd["grid"]["n"])
    dx = Float64(tvd["grid"]["dx"])
    u = Float64(tvd["advection"]["velocity"])
    cfl = Float64(tvd["advection"]["cfl"])
    periods = Float64(tvd["advection"]["periods"])
    eps_denom = Float64(get(tvd, "eps_denom", 1.0e-12))
    tvd_tol = Float64(get(tvd, "tvd_tolerance", 1.0e-10))

    # Cell averages of the smooth+square IC: sin(2πx)·1_{[0,0.4]} + 1·1_{[0.6,0.85]}.
    F_smooth(x) = -cos(2π * x) / (2π)
    function cell_avg(i)
        a = (i - 1) * dx
        b = i * dx
        sa = max(a, 0.0)
        sb = min(b, 0.4)
        smooth = sb > sa ? (F_smooth(sb) - F_smooth(sa)) / dx : 0.0
        qa = max(a, 0.6)
        qb = min(b, 0.85)
        square = qb > qa ? (qb - qa) / dx : 0.0
        return smooth + square
    end
    q0 = [cell_avg(i) for i in 1:n]

    modn(j) = mod(j - 1, n) + 1
    tv(q) = sum(abs(q[modn(i + 1)] - q[i]) for i in 1:n)
    phi(r) = _eval_limiter_ast(formula, Dict("r" => Float64(r)))

    function step(q, dt)
        Fs = Vector{Float64}(undef, n)
        for i in 1:n
            qm = q[modn(i - 1)]; qi = q[i]; qp = q[modn(i + 1)]
            num = qi - qm
            den = qp - qi
            den_safe = abs(den) > eps_denom ? den : (den >= 0 ? eps_denom : -eps_denom)
            r_ratio = num / den_safe
            Fs[i] = u * (qi + 0.5 * phi(r_ratio) * (qp - qi))
        end
        qnew = similar(q)
        for i in 1:n
            qnew[i] = q[i] - dt / dx * (Fs[i] - Fs[modn(i - 1)])
        end
        return qnew
    end

    dt = cfl * dx / abs(u)
    T = periods / abs(u)
    q = copy(q0)
    tsim = 0.0
    while tsim < T
        dts = min(dt, T - tsim)
        q = step(q, dts)
        tsim += dts
    end
    tv_initial = tv(q0)
    tv_final = tv(q)
    return _TVDSummary(tv_final <= tv_initial + tvd_tol && all(isfinite, q), tv_initial, tv_final, tvd_tol)
end

"""
    run_integration_benchmarks(rule, integration_dir) -> LayerResult

Layer-C dispatcher. Reads `integration_dir/cases.json` (a JSON manifest listing
benchmark cases applicable to this rule), dispatches each case through
`IntegrationCases.run_case`, and aggregates:

- All cases PASS or SKIP, at least one PASS  → `LAYER_PASS` with per-case summary
- All cases SKIP                              → `LAYER_SKIP` with skip reasons
- Any case FAIL                               → `LAYER_FAIL` with failure detail

Cases that the runner cannot yet drive (Williamson 2 geostrophic, DCMIP 1-1)
declare `kind: "stub"` and emit a SKIP carrying the manifest's `skip_reason`.
"""
function run_integration_benchmarks(rule::RuleFile, integration_dir::AbstractString)
    cases_path = joinpath(integration_dir, "cases.json")
    if !isfile(cases_path)
        return LayerResult(LAYER_SKIP, "no cases.json in $(relpath_from_repo(integration_dir))")
    end
    cases_json = try
        JSON.parse(read(cases_path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse $(relpath_from_repo(cases_path)): $(sprint(showerror, err))")
    end
    cases = get(cases_json, "cases", nothing)
    if !(cases isa AbstractVector) || isempty(cases)
        return LayerResult(LAYER_SKIP, "cases.json declares no cases")
    end

    # Deferred load of IntegrationCases so the unit tests that exercise
    # walker plumbing (e.g. layer C honors ESD_RUN_INTEGRATION) don't pay
    # the OrdinaryDiffEq compile cost unless integration is actually requested.
    cases_module_path = joinpath(@__DIR__, "integration_cases", "IntegrationCases.jl")
    if !isfile(cases_module_path)
        return LayerResult(LAYER_FAIL, "IntegrationCases module missing at $(relpath_from_repo(cases_module_path))")
    end
    if !isdefined(@__MODULE__, :IntegrationCases)
        Base.include(@__MODULE__, cases_module_path)
    end
    # Resolve `run_case` via `invokelatest(getfield, ...)` so this caller can
    # see a binding that `Base.include` defined in a newer world. The dispatcher
    # returns a `(outcome::Symbol, message::String)` tuple — using a Symbol tag
    # avoids a second world-age hop on the comparison side.
    cases_mod = Base.invokelatest(getfield, @__MODULE__, :IntegrationCases)
    run_case = Base.invokelatest(getfield, cases_mod, :run_case)

    summaries = String[]
    n_pass = 0
    n_skip = 0
    n_fail = 0
    for case in cases
        case_dict = case isa AbstractDict ? case : Dict{String, Any}()
        outcome, message = Base.invokelatest(run_case, case_dict, integration_dir)
        push!(summaries, message)
        if outcome === :pass
            n_pass += 1
        elseif outcome === :skip
            n_skip += 1
        else
            n_fail += 1
        end
    end

    summary = join(summaries, "; ")
    if n_fail > 0
        return LayerResult(LAYER_FAIL, summary)
    elseif n_pass > 0
        return LayerResult(
            LAYER_PASS, "$(n_pass)/$(length(cases)) cases pass" *
                (n_skip > 0 ? ", $(n_skip) skipped" : "") * ": $summary"
        )
    else
        return LayerResult(LAYER_SKIP, "$(n_skip)/$(length(cases)) cases skipped: $summary")
    end
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
        lim = run_layer_limiter(rule)
        push!(results, RuleResult(rule.family, rule.name, rule.path, a, b, c, lim))
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
        println(io, "[$(r.family)/$(r.name)]  A=$(outcome_tag(r.layer_a.outcome))  B=$(outcome_tag(r.layer_b.outcome))  B'=$(outcome_tag(r.layer_limiter.outcome))  C=$(outcome_tag(r.layer_c.outcome))")
        for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("B'", r.layer_limiter), ("C", r.layer_c))
            if lr.outcome != LAYER_PASS && !isempty(lr.reason)
                println(io, "    layer $tag: $(lr.reason)")
            end
        end
    end
    total = length(results) * 4
    npass = sum(count_outcome(r, LAYER_PASS) for r in results; init = 0)
    nfail = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    nskip = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    return println(io, "summary: $total cases  pass=$npass  fail=$nfail  skip=$nskip")
end

function count_outcome(r::RuleResult, target::LayerOutcome)
    return (r.layer_a.outcome == target) +
        (r.layer_b.outcome == target) +
        (r.layer_limiter.outcome == target) +
        (r.layer_c.outcome == target)
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
    total = length(results) * 4
    failures = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    skipped = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io, "<testsuites>")
        println(io, "  <testsuite name=\"ESD Walker\" tests=\"$total\" failures=\"$failures\" skipped=\"$skipped\">")
        for r in results
            for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("limiter", r.layer_limiter), ("C", r.layer_c))
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
