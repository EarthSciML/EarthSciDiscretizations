module WalkESDTests

using EarthSciDiscretizations: load_rules, RuleFile, eval_coeff
import EarthSciSerialization
using JSON

export walk_esd_tests,
    discover_rules,
    run_layer_a,
    run_layer_b,
    run_layer_c,
    run_layer_limiter,
    run_layer_d,
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
    layer_d::LayerResult
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
    if !isfile(input)
        return LayerResult(LAYER_FAIL, "canonical/ present but missing input.esm")
    end
    skip = _fixture_applicable_skip(input)
    skip === nothing || return skip
    if !isfile(expected)
        return LayerResult(LAYER_FAIL, "canonical/ present but missing expected.esm")
    end
    return apply_rule_and_diff(rule, input, expected)
end

function _run_rewrite_variant(rule::RuleFile, rewrite_dir::AbstractString)
    input = joinpath(rewrite_dir, "input.esm")
    expected = joinpath(rewrite_dir, "expected.esm")
    if !isfile(input)
        return LayerResult(LAYER_FAIL, "rewrite/ present but missing input.esm")
    end
    skip = _fixture_applicable_skip(input)
    skip === nothing || return skip
    if !isfile(expected)
        return LayerResult(LAYER_FAIL, "rewrite/ present but missing expected.esm")
    end
    return apply_rewrite_and_diff(rule, input, expected)
end

# Mirror Layer-B's `applicable:false` honoring (`run_mms_convergence`) for
# Layer-A: a canonical/rewrite fixture whose input declares the rule isn't
# yet drivable end-to-end through ESS (e.g. cubed_sphere selectors pending
# dispatch, face-staggered binding contracts) skips with the fixture's
# declared `skip_reason` rather than blowing up `discretize`.
function _fixture_applicable_skip(input_path::AbstractString)
    parsed = try
        JSON.parse(read(input_path, String))
    catch
        return nothing
    end
    parsed isa AbstractDict || return nothing
    get(parsed, "applicable", true) === false || return nothing
    reason = get(parsed, "skip_reason", "fixture declares applicable:false (no reason given)")
    return LayerResult(LAYER_SKIP, "fixture-declared not applicable: $(reason)")
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

MMS convergence (Layer-B). Defers to `run_mms_convergence`, which honors
`applicable:false` markers and otherwise SKIPs pending the canonical-pipeline
replacement for the retired ESS `verify_mms_convergence` (esm-4t5 — see
`run_mms_convergence` docstring).
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
    run_layer_d(rule) -> LayerResult

Layer D — discrete conservation for finite-volume rules. Verifies that the
rule, applied to a periodic domain, conserves the integrated quantity to
roundoff: the global sum of the rule output (or the per-cell update derived
from it) satisfies `|∑ Δq| < tol`. Mathematically this follows from flux
telescoping on a closed domain — the layer is a regression guard against
sign errors, missing wraps, and stencil drift.

Looks for a `conservation/` fixture directory at
`discretizations/<family>/<rule>/fixtures/conservation/conservation_check.esm`.
Skips when the directory is absent (the common case for non-FV rules); fails
when the fixture file is missing or the conservation property is violated.

Two fixture `kind`s are supported (see fixture files for full schema):

- `conservation_divergence_2d_periodic` — fluxes are sampled on Arakawa-C-style
  staggered faces of a 2D periodic rectangle; the rule's stencil is applied at
  every cell and the global sum of the divergence is checked.
- `conservation_muscl_1d_periodic` — runs the MUSCL upwind scheme
  `F_{i+1/2} = u·(q_i + 0.5·phi(r_i)·(q_{i+1}-q_i))` (same form as the layer-B'
  TVD check) on a 1D periodic domain, then checks that `∑(q_final - q_0) ≈ 0`.
"""
function run_layer_d(rule::RuleFile)
    conservation = joinpath(rule_fixtures_dir(rule), "conservation")
    if !isdir(conservation)
        return LayerResult(LAYER_SKIP, "no conservation fixtures at $(relpath_from_repo(conservation))")
    end
    return run_conservation_check(rule, conservation)
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

Layer B: read the convergence fixture, honor `applicable:false` markers, and
otherwise SKIP pending the per-topology canonical-pipeline replacement.

ESS retired `verify_mms_convergence` + `mms_evaluator` in esm-4t5 (the
2026-04-29 single-evaluation-pathway directive). Layer B previously routed
the rule + fixture through that ESS harness; the harness was a shadow
evaluator outside the canonical `discretize → ArrayOp → simulate` pipeline
and the directive forbade it. The replacement drives Layer B through the
canonical pipeline + an official ESS simulation runner per AGENTS.md
(`discretize` → MTK.System / Symbolics scalarize → evaluate `f!(du, u, p, 0)`
at the manufactured-solution sample → L_inf vs analytic).

The replacement lands per topology family (one follow-up bead per family
under dsc-kswm). Until the corresponding family's runner lands, fixtures
in that family SKIP with the unified `_LAYER_B_PIPELINE_PENDING` reason,
suffixed with the topology key so the report names the missing piece.
"""
const _LAYER_B_PIPELINE_PENDING = "Layer-B awaits canonical-pipeline replacement: ESS retired verify_mms_convergence in esm-4t5 (2026-04-29 single-pathway directive); Layer-B replacement via discretize → ArrayOp → official ESS simulation runner tracked at ESD/dsc-kswm"

# ---------------------------------------------------------------------------
# Generic MMS catalog. Indexed by `mms_kind` declared in the convergence
# fixture (a fixture attribute, NOT a rule identity), each entry maps a
# cell-center coordinate to (initial-condition, analytic-derivative) pairs.
# Layer-B runners look up the MMS by its declared `mms_kind` and never
# branch on the rule's name or shape — extension is by adding entries here
# and the corresponding `mms_kind` field in fixture inputs.
# ---------------------------------------------------------------------------
const _LAYER_B_MMS_CATALOG = Dict{String, NamedTuple{(:ic, :derivative), Tuple{Function, Function}}}(
    # 1D periodic sine: u(x) = sin(2π x); u'(x) = 2π cos(2π x). Used by
    # centered_2nd_uniform / upwind_1st convergence fixtures (the
    # historic Layer-B MMS for 1D Cartesian rules).
    "sin_2pi_x_periodic" => (
        ic = x -> sin(2π * x),
        derivative = x -> 2π * cos(2π * x),
    ),
    # 1D vertical-column sine on the unit column z ∈ [0, 1]: u(z) = sin(2π z);
    # u'(z) = 2π cos(2π z). Naturally periodic on [0, 1] (sin(0) = sin(2π) = 0),
    # so the centered_2nd_uniform_vertical canonical fixture (which declares
    # periodic=true on its single dimension) and any future non-periodic
    # vertical fixture both consume this MMS shape — only the boundary
    # treatment in the topology runner differs.
    "sin_2pi_z_unit_column" => (
        ic = z -> sin(2π * z),
        derivative = z -> 2π * cos(2π * z),
    ),
)

# ---------------------------------------------------------------------------
# Topology dispatch. The Layer-B runner classifies each (rule, fixture)
# pair into a topology key derived from generic attributes:
#   - rule.json `grid_family`
#   - rule.json `applies_to.args` length (1 for scalar-input rules)
#   - canonical fixture's grid `dimensions` shape (1D vs 2D, periodic etc.)
# A canonical fixture is required because Layer-B reuses it as the ESM
# template, parameterizing only the grid `size`. Rules that lack a
# canonical fixture cannot be Layer-B-driven yet (they ship neither the
# ESS-format rule inline nor a parameterizable model template).
#
# This dispatch is the legitimate generic kind: it routes on ESM/rule
# schema attributes, not on rule identity. Per-rule-shape dispatch (e.g.
# `if rule.name == "centered_2nd_uniform" then ...`) is forbidden.
# ---------------------------------------------------------------------------
const _LAYER_B_SUPPORTED_TOPOLOGIES = Set{String}([
    # 1D vertical column (dsc-yz0m) — first per-topology runner. Drives the
    # canonical fixture's discretized RHS through `discretize` →
    # generic per-cell scalarization → `EarthSciSerialization.build_evaluator`
    # (the official ESS tree-walk simulation runner) → `f!(du, u, p, 0)` and
    # compares against the MMS analytic derivative at cell centers.
    "1d_vertical_column",
])

# Map topology_key → follow-up bead tracking the implementation. Used in
# SKIP reasons so the walker report points to the open work.
const _LAYER_B_TOPOLOGY_TRACKING = Dict{String, String}(
    "1d_cartesian_periodic"  => "ESD/dsc-k1li",
    "1d_vertical_column"     => "ESD/dsc-yz0m",
    "2d_cartesian_periodic"  => "ESD/dsc-vst2",
    "2d_arakawa_periodic"    => "ESD/dsc-70zp",
    "2d_latlon_sphere"       => "ESD/dsc-mps8",
    "fv_cell_average_1d"     => "ESD/dsc-a7b2",
    "stencil_form_rule"      => "ESD/dsc-y0jj (prerequisite for ESD/dsc-k1li)",
    "unsupported"            => "no follow-up bead — out of Layer-B scope",
)

function run_mms_convergence(rule::RuleFile, convergence_dir::AbstractString)
    input_path = joinpath(convergence_dir, "input.esm")
    if !isfile(input_path)
        return LayerResult(LAYER_FAIL, "convergence/ present but missing input.esm")
    end
    input_json = JSON.parse(read(input_path, String))

    # Fixture-declared non-applicability: rules whose acceptance signature
    # isn't a manufactured-solution convergence sweep (index-rewrite rules,
    # TVD limiters, reconstruction-style rules pending ESS harness extension)
    # ship an `applicable: false` + `skip_reason` marker so the walker
    # surfaces a structured SKIP. Honored before any further fixture handling
    # so structurally skip-only fixtures (e.g. cubed_sphere η-sibling) need
    # not duplicate the expected payload.
    if get(input_json, "applicable", true) === false
        reason = get(input_json, "skip_reason", "fixture declares applicable:false (no reason given)")
        return LayerResult(LAYER_SKIP, "fixture-declared not applicable: $(reason)")
    end

    # Validate `expected.esm` exists and parses an `expected_min_order`.
    # Per-topology runners use this when promoting from SKIP to PASS/FAIL;
    # validating it here surfaces malformed expected fixtures before any
    # topology runner lands.
    expected_path = joinpath(convergence_dir, "expected.esm")
    if !isfile(expected_path)
        return LayerResult(LAYER_FAIL, "convergence/ present but missing expected.esm")
    end
    expected_json = try
        JSON.parse(read(expected_path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse $(relpath_from_repo(expected_path)): $(sprint(showerror, err))")
    end
    raw_min_order = get(expected_json, "expected_min_order", nothing)
    if raw_min_order === nothing
        return LayerResult(LAYER_FAIL, "$(relpath_from_repo(expected_path)) missing required `expected_min_order`")
    end
    expected_min_order = Float64(raw_min_order)
    if !isfinite(expected_min_order)
        return LayerResult(LAYER_FAIL, "$(relpath_from_repo(expected_path)) `expected_min_order` is not finite: $raw_min_order")
    end

    # Classify into a topology family using only generic ESM/rule schema
    # attributes (see `_LAYER_B_SUPPORTED_TOPOLOGIES` docstring above).
    topology_key = _layer_b_topology_key(rule, input_json)
    if !(topology_key in _LAYER_B_SUPPORTED_TOPOLOGIES)
        bead = get(_LAYER_B_TOPOLOGY_TRACKING, topology_key, "no follow-up bead recorded")
        return LayerResult(LAYER_SKIP, _LAYER_B_PIPELINE_PENDING * " (topology=$topology_key; tracked at $(bead))")
    end

    # Generic MMS lookup. Only fixtures that declare a registered `mms_kind`
    # can be driven; bare descriptive `manufactured_solution` strings SKIP
    # until the fixture is upgraded.
    mms_kind = String(get(input_json, "mms_kind", ""))
    mms = get(_LAYER_B_MMS_CATALOG, mms_kind, nothing)
    if mms === nothing
        return LayerResult(LAYER_SKIP, _LAYER_B_PIPELINE_PENDING * " (mms_kind=$(repr(mms_kind)) not registered in Layer-B MMS catalog)")
    end

    # Per-topology runner. Each family's implementation lands in a
    # follow-up bead under dsc-kswm. The dispatcher is the only place that
    # routes on topology_key — `_run_layer_b_canonical_grid` itself walks
    # the canonical-form AST through MTK / Symbolics generically.
    grids_raw = get(input_json, "grids", Any[])
    if !(grids_raw isa AbstractVector) || isempty(grids_raw)
        return LayerResult(LAYER_FAIL, "convergence/input.esm declares no grids")
    end
    grids = Int[Int(get(g, "n", 0)) for g in grids_raw]
    if any(<=(0), grids)
        return LayerResult(LAYER_FAIL, "convergence/input.esm has non-positive grid sizes: $grids")
    end

    errors = Float64[]
    for n in grids
        err = try
            _run_layer_b_canonical_grid(topology_key, rule, mms, n)
        catch e
            return LayerResult(LAYER_FAIL, "canonical pipeline threw at n=$n: $(sprint(showerror, e))")
        end
        if !(err isa Real) || !isfinite(err)
            return LayerResult(LAYER_FAIL, "canonical pipeline at n=$n returned non-finite/non-numeric: $(typeof(err))")
        end
        push!(errors, Float64(err))
    end

    if any(iszero, errors)
        return LayerResult(LAYER_FAIL, "zero error in convergence sweep (likely degenerate fixture): $errors")
    end
    if length(errors) < 2
        return LayerResult(LAYER_FAIL, "need ≥ 2 grids for convergence-order calc, got $(length(errors))")
    end
    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    min_order = minimum(orders)
    if min_order >= expected_min_order
        return LayerResult(
            LAYER_PASS,
            "min order $(round(min_order; digits = 2)) >= expected $(expected_min_order) over $(length(grids)) grids " *
                "(orders=$(round.(orders; digits = 2)))",
        )
    else
        return LayerResult(
            LAYER_FAIL,
            "min order $(round(min_order; digits = 2)) below expected $(expected_min_order) " *
                "(orders=$(round.(orders; digits = 2)), errors=$(round.(errors; sigdigits = 3)))",
        )
    end
end

"""
    _layer_b_topology_key(rule, input_json) -> String

Classify a (rule, convergence-fixture) pair into a topology family name.
The classifier reads only generic ESM/rule schema attributes — never the
rule's identity — so adding a new rule that fits an existing family
auto-routes to the right runner without code changes here.

Topology keys (closed set):

- `"1d_cartesian_periodic"` — `grid_family=cartesian`, single-arg pattern,
  canonical fixture declares one periodic dimension.
- `"1d_vertical_column"` — single-axis vertical-grid rules (`grid_family=vertical`).
  The seed fixture (`centered_2nd_uniform_vertical`) declares `periodic=true`
  on its single dimension because the MMS sin(2π z) is naturally periodic on
  the unit column; the runner reads the dimension's `periodic` flag and
  selects either modular wrap (periodic) or clamp (non-periodic, e.g. when
  a future fixture pairs vertical staggering with one-sided BCs).
- `"2d_cartesian_periodic"` — two periodic Cartesian axes.
- `"2d_latlon_sphere"` — `grid_family=latlon` (or sphere variant).
- `"fv_cell_average_1d"` — single arg pattern but `sampling=cell_average`
  in the convergence fixture (FV reconstruction-style).
- `"stencil_form_rule"` — rule.json declares `stencil` instead of
  `replacement` (cannot be passed through ESS rule engine until lowered
  to replacement form).
- `"unsupported"` — multi-arg patterns, missing canonical fixture, or any
  shape this classifier hasn't enumerated.

This function returns the key name; only keys present in
`_LAYER_B_SUPPORTED_TOPOLOGIES` get a real runner. Others SKIP with a
reason that names the tracking bead.
"""
function _layer_b_topology_key(rule::RuleFile, input_json::AbstractDict)
    rule_doc = try
        JSON.parse(read(rule.path, String))
    catch
        return "unsupported"
    end
    spec = get(get(rule_doc, "discretizations", Dict{String, Any}()), rule.name, nothing)
    spec isa AbstractDict || return "unsupported"

    grid_family = String(get(spec, "grid_family", ""))

    # Vertical-family rules ride the canonical pipeline via their canonical
    # fixture's embedded rule (replacement-form), which `discretize` already
    # consumes for Layer A. The top-level rule.json may still carry only a
    # `stencil` form (the eventual stencil-form lowering is dsc-y0jj's
    # scope); routing on `grid_family` here keeps the vertical runner
    # reachable in the meantime. Other stencil-only rules still gate on the
    # stencil-form check below.
    grid_family == "vertical" && return "1d_vertical_column"

    # `stencil`-form rules cannot ride the ESS rule engine until their
    # stencil entries are lowered to a `replacement` AST. Tracked
    # separately so a future bead can either author the lowering or
    # ship per-rule canonical fixtures with replacement-form rules.
    if haskey(spec, "stencil") && !haskey(spec, "replacement")
        return "stencil_form_rule"
    end
    applies_to = get(spec, "applies_to", nothing)
    applies_to isa AbstractDict || return "unsupported"
    args = get(applies_to, "args", Any[])
    n_args = length(args)
    sampling = String(get(input_json, "sampling", "cell_center"))

    # Lat-lon sphere is its own topology family, regardless of arg count
    # (Y_l_m manufactured solutions need spherical geometry handling).
    grid_family == "latlon" && return "2d_latlon_sphere"

    # Arakawa C-grid (staggered face-centered fluxes on 2D Cartesian) gets
    # its own family because the canonical-pipeline runner has to handle
    # face-staggered velocity fields, distinct from the cell-centered 2D
    # Cartesian path.
    grid_family == "arakawa" && return "2d_arakawa_periodic"

    # FV cell-average rules (PPM / WENO reconstruction) need their own
    # IC sampling (cell averages, not point samples) and exact-derivative
    # comparison, distinct from FD point-sample rules.
    if sampling == "cell_average"
        return "fv_cell_average_1d"
    end

    # 2D Cartesian rules: diagnose by the convergence fixture's `axis`
    # field (a 2D MMS sweep over a single axis at a time, as used by
    # weno5_advection_2d). The classifier deliberately does NOT pattern-
    # match on the manufactured_solution description string — that field
    # is free-form prose, not machine-readable.
    if grid_family == "cartesian" && haskey(input_json, "axis")
        return "2d_cartesian_periodic"
    end

    # 1D Cartesian: the default for single-arg cartesian rules. Periodic
    # vs vertical (non-periodic column) is decided by the canonical
    # fixture's grid declaration if available; otherwise default to
    # "1d_cartesian_periodic" — the runner gates on canonical-fixture
    # presence and SKIPs there if absent.
    if grid_family == "cartesian" && n_args == 1
        canonical_input = joinpath(dirname(rule.path), rule.name, "fixtures", "canonical", "input.esm")
        if isfile(canonical_input)
            canonical = try
                JSON.parse(read(canonical_input, String))
            catch
                return "1d_cartesian_periodic"
            end
            grids = get(canonical, "grids", Dict{String, Any}())
            if grids isa AbstractDict && !isempty(grids)
                g = first(values(grids))
                if g isa AbstractDict
                    dims = get(g, "dimensions", Any[])
                    if dims isa AbstractVector && length(dims) == 1
                        d = first(dims)
                        if d isa AbstractDict
                            return get(d, "periodic", false) === true ?
                                "1d_cartesian_periodic" : "1d_vertical_column"
                        end
                    elseif dims isa AbstractVector && length(dims) >= 2
                        return "2d_cartesian_periodic"
                    end
                end
            end
        end
        # No canonical fixture available — assume 1D periodic by default;
        # the runner will SKIP with "no canonical fixture" reason.
        return "1d_cartesian_periodic"
    end
    return "unsupported"
end

"""
    _run_layer_b_canonical_grid(topology_key, rule, mms, n) -> Float64

Drive one resolution of the canonical convergence pipeline. Returns the
L_inf error of the rule's discretized output vs the analytic derivative
at cell centers on a grid of size `n`.

Per-topology implementations live in dsc-kswm follow-up beads. The
dispatcher routes on `topology_key` only — the per-topology functions
walk the canonical AST generically (no per-rule-shape dispatch inside).
Until the matching family lands, this dispatcher errors with a reference
to the tracking bead.
"""
function _run_layer_b_canonical_grid(topology_key::AbstractString, rule::RuleFile, mms, n::Int)
    if topology_key == "1d_vertical_column"
        return _run_layer_b_1d_vertical_column(rule, mms, n)
    end
    error("Layer-B canonical-pipeline runner for topology=$(topology_key) not implemented; " *
          "tracked in a dsc-kswm follow-up bead. The dispatcher should not have been reached " *
          "for this topology key.")
end

"""
    _run_layer_b_1d_vertical_column(rule, mms, n) -> Float64

Layer-B runner for 1D vertical-column rules (dsc-yz0m). Drives the rule's
canonical fixture through the canonical pipeline:

1. Read `fixtures/canonical/input.esm`, set the single dimension's `size`
   to `n`, run `EarthSciSerialization.discretize` (rule engine + canonical-
   form rewrite). After this step the equation RHS is the rule's
   replacement AST with `index(u, k+offset)` ops in shape-[k] form.
2. Identify the column variable (the unique state var with `shape=[<dim>]`)
   and scalarize the model into per-cell scalar variables `<u>_1 .. <u>_n`,
   substituting `index(<u>, expr)` → `<u>_<wrap(eval(expr))>` and the
   dimension symbol → cell index `i`. Periodic wrap (or clamp) is selected
   from the canonical dimension's `periodic` flag. The grid spacing
   parameter `h` is bound as a constant `1/n` on the unit column.
3. Hand the scalarized ESM to `EarthSciSerialization.discretize` (no rules
   to apply, but the call is the canonical-pipeline tag) and then to
   `EarthSciSerialization.build_evaluator` — the official ESS tree-walk
   simulation runner — which compiles the per-cell RHS into `f!(du, u, p, t)`.
4. Sample the MMS initial condition at cell centers `z_i = (i - 0.5)/n`,
   call `f!(du, u0, p, 0.0)`, and return the L_inf error vs
   `mms.derivative(z_i)`.

The scalarization in step 2 walks the AST generically — it dispatches on
op name (`index` / dim symbol / arithmetic ops), not on rule identity, so
adding another vertical-column rule with a different RHS shape requires no
changes here.
"""
function _run_layer_b_1d_vertical_column(rule::RuleFile, mms, n::Int)
    canonical = joinpath(dirname(rule.path), rule.name, "fixtures", "canonical", "input.esm")
    isfile(canonical) ||
        error("1d_vertical_column runner requires canonical/input.esm at $(relpath_from_repo(canonical))")
    template = JSON.parse(read(canonical, String))

    # Locate grid + dimension + the single column state variable.
    models_doc = get(template, "models", nothing)
    (models_doc isa AbstractDict && !isempty(models_doc)) ||
        error("canonical fixture has no models")
    model_name, model = first(pairs(models_doc))
    grid_name = String(get(model, "grid", ""))
    grids = get(template, "grids", Dict{String, Any}())
    grid = get(grids, grid_name, nothing)
    grid isa AbstractDict ||
        error("canonical fixture references unknown grid $(repr(grid_name))")
    dims = get(grid, "dimensions", Any[])
    (dims isa AbstractVector && length(dims) == 1) ||
        error("1d_vertical_column expects a single-dimension grid; got $(length(dims))")
    dim = dims[1]
    dim_name = String(get(dim, "name", ""))
    periodic = get(dim, "periodic", false) === true

    # Parameterize grid size to `n`.
    dim["size"] = n

    discretized = EarthSciSerialization.discretize(template)
    disc_model = first(values(get(discretized, "models", Dict{String, Any}())))
    disc_vars = get(disc_model, "variables", Dict{String, Any}())

    col_var = nothing
    for (vname, vspec) in disc_vars
        vspec isa AbstractDict || continue
        get(vspec, "type", "") == "state" || continue
        shape = get(vspec, "shape", nothing)
        if shape isa AbstractVector && length(shape) == 1 && String(shape[1]) == dim_name
            col_var = String(vname)
            break
        end
    end
    col_var === nothing &&
        error("could not identify column state variable (no var with shape=[$(dim_name)])")

    # Cell centers on the unit column z ∈ [0, 1].
    h = 1.0 / n
    z(i) = (i - 0.5) * h

    # Build the scalarized ESM. Vars: u_1 .. u_n (state, IC = mms.ic(z_i)),
    # plus `h` as a parameter so the discretized RHS's free `h` symbol binds.
    scalar_vars = Dict{String, Any}(
        "h" => Dict{String, Any}(
            "type" => "parameter",
            "default" => h,
            "units" => "1",
        ),
    )
    for i in 1:n
        scalar_vars["$(col_var)_$(i)"] = Dict{String, Any}(
            "type" => "state",
            "default" => mms.ic(z(i)),
            "units" => "1",
        )
    end

    scalar_equations = Any[]
    for eq in get(disc_model, "equations", Any[])
        lhs = eq["lhs"]
        rhs = eq["rhs"]
        # The column equation has lhs `D(<col_var>, wrt=t)`; expand into
        # one scalar derivative per cell. Other equations (none expected
        # for the seed fixture) pass through unchanged.
        if !(lhs isa AbstractDict) || String(get(lhs, "op", "")) != "D"
            push!(scalar_equations, eq)
            continue
        end
        lhs_args = get(lhs, "args", Any[])
        (length(lhs_args) >= 1 && String(lhs_args[1]) == col_var) || begin
            push!(scalar_equations, eq)
            continue
        end
        for i in 1:n
            push!(
                scalar_equations,
                Dict{String, Any}(
                    "lhs" => Dict{String, Any}(
                        "op" => "D",
                        "args" => Any["$(col_var)_$(i)"],
                        "wrt" => "t",
                    ),
                    "rhs" => _scalarize_per_cell(rhs, dim_name, i, col_var, n, periodic),
                ),
            )
        end
    end

    scalar_esm = Dict{String, Any}(
        "esm" => "0.2.0",
        "metadata" => Dict{String, Any}(
            "name" => "layer_b_$(rule.name)_n$(n)",
        ),
        "models" => Dict{String, Any}(
            String(model_name) => Dict{String, Any}(
                "variables" => scalar_vars,
                "equations" => scalar_equations,
            ),
        ),
    )

    # Canonical pipeline tag + official ESS tree-walk runner.
    discretized_scalar = EarthSciSerialization.discretize(scalar_esm)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(discretized_scalar)
    du = similar(u0)
    f!(du, u0, p, 0.0)

    err = 0.0
    for i in 1:n
        idx = var_map["$(col_var)_$(i)"]
        analytic = mms.derivative(z(i))
        err = max(err, abs(du[idx] - analytic))
    end
    return err
end

# Generic per-cell scalarization of a discretized AST. Walks the AST and
# substitutes `index(<col_var>, expr)` with the corresponding scalar var
# name and `<dim_name>` (the dimension symbol) with the cell index `i`.
# Other ops recurse without rule-shape branching, so rules whose RHS uses
# different op shapes ride this same path.
function _scalarize_per_cell(node, dim_name::String, i::Int, col_var::String, n::Int, periodic::Bool)
    if node isa AbstractDict
        op = String(get(node, "op", ""))
        args = get(node, "args", Any[])
        if op == "index" && length(args) == 2
            varname = args[1] isa AbstractString ? String(args[1]) : ""
            if varname == col_var
                idx_value = _eval_index_expr(args[2], dim_name, i)
                wrapped = if periodic
                    mod(idx_value - 1, n) + 1
                else
                    clamp(idx_value, 1, n)
                end
                return "$(col_var)_$(wrapped)"
            end
        end
        new_node = Dict{String, Any}(String(k) => v for (k, v) in node)
        new_node["args"] = Any[
            _scalarize_per_cell(a, dim_name, i, col_var, n, periodic) for a in args
        ]
        return new_node
    elseif node isa AbstractString
        return String(node) == dim_name ? i : node
    end
    return node
end

# Tiny integer-arithmetic evaluator for index expressions like `k+1` or
# `k-1`. Only the ops that appear in finite-difference index arithmetic
# need to be supported; anything else is a fixture or rule bug and surfaces
# as an error so the walker FAILs loudly rather than silently mis-indexing.
function _eval_index_expr(node, dim_name::String, i::Int)::Int
    if node isa Bool
        # Reject explicitly so a stray boolean does not coerce to 0/1.
        error("boolean literal in index expression")
    elseif node isa Integer
        return Int(node)
    elseif node isa AbstractString
        s = String(node)
        s == dim_name && return i
        error("unknown symbol $(repr(s)) in index expression (dim=$(repr(dim_name)))")
    elseif node isa AbstractDict
        op = String(get(node, "op", ""))
        args = get(node, "args", Any[])
        if op == "+"
            return sum(_eval_index_expr(a, dim_name, i) for a in args)
        elseif op == "-" && length(args) == 1
            return -_eval_index_expr(args[1], dim_name, i)
        elseif op == "-" && length(args) == 2
            return _eval_index_expr(args[1], dim_name, i) - _eval_index_expr(args[2], dim_name, i)
        elseif op == "*"
            acc = 1
            for a in args
                acc *= _eval_index_expr(a, dim_name, i)
            end
            return acc
        elseif op == "neg" && length(args) == 1
            return -_eval_index_expr(args[1], dim_name, i)
        end
        error("unsupported op $(repr(op)) in index expression")
    end
    error("cannot evaluate index node of type $(typeof(node))")
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
            eval_coeff(formula, Dict("\$r" => r))
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
            eval_coeff(formula, Dict("\$r" => rf))
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
    phi_one = eval_coeff(formula, Dict("\$r" => 1.0))
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
        "Sweby OK ($(n_ref) refs, $(n_sweep) sweep pts); TVD OK (TV $(round(tv_summary.tv_initial; digits = 6)) -> $(round(tv_summary.tv_final; digits = 6)))",
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
    phi(r) = eval_coeff(formula, Dict("\$r" => Float64(r)))

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

# ---------------------------------------------------------------------------
# Layer D (discrete conservation) — generic over fixture `kind`. Two kinds
# are wired up; new kinds can be added without touching the dispatcher.
# ---------------------------------------------------------------------------

function run_conservation_check(rule::RuleFile, conservation_dir::AbstractString)
    fixture_path = joinpath(conservation_dir, "conservation_check.esm")
    if !isfile(fixture_path)
        return LayerResult(
            LAYER_FAIL,
            "conservation/ present but missing conservation_check.esm",
        )
    end
    fixture = try
        JSON.parse(read(fixture_path, String))
    catch err
        return LayerResult(
            LAYER_FAIL,
            "failed to parse $(relpath_from_repo(fixture_path)): $(sprint(showerror, err))",
        )
    end
    kind = get(fixture, "kind", nothing)
    if kind == "conservation_divergence_2d_periodic"
        return _run_conservation_divergence_2d_periodic(rule, fixture)
    elseif kind == "conservation_muscl_1d_periodic"
        return _run_conservation_muscl_1d_periodic(rule, fixture)
    else
        return LayerResult(
            LAYER_FAIL,
            "$(relpath_from_repo(fixture_path)) has unsupported kind $(repr(kind))",
        )
    end
end

# Synthesize a periodic flux field from a deterministic pseudo-random seed.
# A linear congruential generator is used so the fixture is bit-reproducible
# across Julia versions without depending on a stdlib RNG implementation.
function _seeded_flux_field(seed::UInt64, n::Int)
    state = seed == 0 ? UInt64(0x9E3779B97F4A7C15) : seed
    out = Vector{Float64}(undef, n)
    for i in 1:n
        # Numerical Recipes LCG (Knuth) — full-period on UInt64.
        state = (state * 6364136223846793005 + 1442695040888963407) % typemax(UInt64)
        # Map the high 53 bits to [-1, 1).
        u = (state >> 11) / Float64(1 << 53)
        out[i] = 2.0 * u - 1.0
    end
    return out
end

function _run_conservation_divergence_2d_periodic(rule::RuleFile, fixture::AbstractDict)
    rule_json = try
        JSON.parse(read(rule.path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse rule $(relpath_from_repo(rule.path)): $(sprint(showerror, err))")
    end
    spec = get(get(rule_json, "discretizations", Dict()), rule.name, nothing)
    if !(spec isa AbstractDict) || !haskey(spec, "stencil")
        return LayerResult(
            LAYER_FAIL,
            "rule $(rule.name) has no `stencil` under discretizations.$(rule.name)",
        )
    end
    stencil = spec["stencil"]
    if !(stencil isa AbstractVector) || isempty(stencil)
        return LayerResult(LAYER_FAIL, "rule $(rule.name) has empty stencil")
    end

    grid = get(fixture, "grid", nothing)
    grid isa AbstractDict ||
        return LayerResult(LAYER_FAIL, "fixture missing grid object")
    nx = Int(get(grid, "nx", 0))
    ny = Int(get(grid, "ny", 0))
    dx = Float64(get(grid, "dx", 0.0))
    dy = Float64(get(grid, "dy", 0.0))
    (nx > 0 && ny > 0 && dx > 0 && dy > 0) ||
        return LayerResult(LAYER_FAIL, "fixture grid must declare positive nx, ny, dx, dy")

    seed = UInt64(get(fixture, "seed", 42))
    Fx = reshape(_seeded_flux_field(seed, nx * ny), nx, ny)
    Fy = reshape(_seeded_flux_field(seed + 1, nx * ny), nx, ny)
    bindings = Dict{String, Float64}("dx" => dx, "dy" => dy)

    # Evaluate the rule's stencil at every cell; sum the per-cell divergence.
    # Periodic wrap is applied on both axes (modulo nx / ny). The stencil
    # entries are interpreted via their `selector.stagger` / `selector.axis`
    # so we don't hard-code the rule's algebraic shape.
    div_sum = 0.0
    for j in 1:ny, i in 1:nx
        cell_div = 0.0
        for entry in stencil
            sel = entry["selector"]
            kind = String(get(sel, "kind", ""))
            kind == "arakawa" ||
                return LayerResult(LAYER_FAIL, "selector kind $(repr(kind)) not supported by layer-D 2D periodic runner")
            stagger = String(get(sel, "stagger", ""))
            axis = String(get(sel, "axis", ""))
            offset = Int(get(sel, "offset", 0))
            coeff = eval_coeff(entry["coeff"], bindings)
            value = if stagger == "face_x" && axis == "\$x"
                Fx[mod(i - 1 + offset, nx) + 1, j]
            elseif stagger == "face_y" && axis == "\$y"
                Fy[i, mod(j - 1 + offset, ny) + 1]
            else
                return LayerResult(
                    LAYER_FAIL,
                    "unsupported (stagger=$(stagger), axis=$(axis)) for layer-D 2D periodic runner",
                )
            end
            cell_div += coeff * value
        end
        div_sum += cell_div
    end

    tol = Float64(get(fixture, "tolerance", 1.0e-12))
    abs_sum = abs(div_sum)
    if abs_sum > tol
        return LayerResult(
            LAYER_FAIL,
            "global divergence sum |∑ ∇·F| = $(abs_sum) exceeds tol=$(tol) on $(nx)×$(ny) periodic grid",
        )
    end
    return LayerResult(
        LAYER_PASS,
        "conservation OK: |∑ ∇·F| = $(abs_sum) ≤ $(tol) on $(nx)×$(ny) periodic grid",
    )
end

function _run_conservation_muscl_1d_periodic(rule::RuleFile, fixture::AbstractDict)
    rule_json = try
        JSON.parse(read(rule.path, String))
    catch err
        return LayerResult(LAYER_FAIL, "failed to parse rule $(relpath_from_repo(rule.path)): $(sprint(showerror, err))")
    end
    spec = get(get(rule_json, "discretizations", Dict()), rule.name, nothing)
    if !(spec isa AbstractDict) || !haskey(spec, "formula")
        return LayerResult(
            LAYER_FAIL,
            "rule $(rule.name) has no `formula` AST under discretizations.$(rule.name)",
        )
    end
    formula = spec["formula"]

    grid = get(fixture, "grid", nothing)
    grid isa AbstractDict ||
        return LayerResult(LAYER_FAIL, "fixture missing grid object")
    n = Int(get(grid, "n", 0))
    dx = Float64(get(grid, "dx", 0.0))
    (n > 0 && dx > 0) ||
        return LayerResult(LAYER_FAIL, "fixture grid must declare positive n, dx")

    advection = get(fixture, "advection", Dict())
    u = Float64(get(advection, "velocity", 1.0))
    cfl = Float64(get(advection, "cfl", 0.4))
    periods = Float64(get(advection, "periods", 1.0))
    eps_denom = Float64(get(fixture, "eps_denom", 1.0e-12))
    tol = Float64(get(fixture, "tolerance", 1.0e-12))

    # Reuse the layer-B' smooth+square initial condition by default; explicit
    # `initial_condition.values` lets a fixture override with a literal array.
    ic_values = get(fixture, "initial_condition_values", nothing)
    q0 = if ic_values isa AbstractVector
        length(ic_values) == n ||
            return LayerResult(LAYER_FAIL, "initial_condition_values length $(length(ic_values)) != n=$(n)")
        Float64[Float64(v) for v in ic_values]
    else
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
        Float64[cell_avg(i) for i in 1:n]
    end

    modn(j) = mod(j - 1, n) + 1
    phi(r) = eval_coeff(formula, Dict("\$r" => Float64(r)))
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
    nsteps = 0
    while tsim < T
        dts = min(dt, T - tsim)
        q = step(q, dts)
        tsim += dts
        nsteps += 1
    end

    delta_sum = 0.0
    for i in 1:n
        delta_sum += q[i] - q0[i]
    end
    abs_sum = abs(delta_sum)
    if abs_sum > tol
        return LayerResult(
            LAYER_FAIL,
            "MUSCL update violates conservation: |∑ Δq| = $(abs_sum) > tol=$(tol) after $(nsteps) steps",
        )
    end
    return LayerResult(
        LAYER_PASS,
        "conservation OK: |∑ Δq| = $(abs_sum) ≤ $(tol) over $(nsteps) MUSCL steps on n=$(n) periodic grid",
    )
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
        d = run_layer_d(rule)
        push!(results, RuleResult(rule.family, rule.name, rule.path, a, b, c, lim, d))
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
        println(io, "[$(r.family)/$(r.name)]  A=$(outcome_tag(r.layer_a.outcome))  B=$(outcome_tag(r.layer_b.outcome))  B'=$(outcome_tag(r.layer_limiter.outcome))  C=$(outcome_tag(r.layer_c.outcome))  D=$(outcome_tag(r.layer_d.outcome))")
        for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("B'", r.layer_limiter), ("C", r.layer_c), ("D", r.layer_d))
            if lr.outcome != LAYER_PASS && !isempty(lr.reason)
                println(io, "    layer $tag: $(lr.reason)")
            end
        end
    end
    total = length(results) * 5
    npass = sum(count_outcome(r, LAYER_PASS) for r in results; init = 0)
    nfail = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    nskip = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    return println(io, "summary: $total cases  pass=$npass  fail=$nfail  skip=$nskip")
end

function count_outcome(r::RuleResult, target::LayerOutcome)
    return (r.layer_a.outcome == target) +
        (r.layer_b.outcome == target) +
        (r.layer_limiter.outcome == target) +
        (r.layer_c.outcome == target) +
        (r.layer_d.outcome == target)
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
    total = length(results) * 5
    failures = sum(count_outcome(r, LAYER_FAIL) for r in results; init = 0)
    skipped = sum(count_outcome(r, LAYER_SKIP) for r in results; init = 0)
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
        println(io, "<testsuites>")
        println(io, "  <testsuite name=\"ESD Walker\" tests=\"$total\" failures=\"$failures\" skipped=\"$skipped\">")
        for r in results
            for (tag, lr) in (("A", r.layer_a), ("B", r.layer_b), ("limiter", r.layer_limiter), ("C", r.layer_c), ("D", r.layer_d))
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
