module WalkESDTests

using EarthSciDiscretizations: load_rules, RuleFile, eval_coeff, lower_stencil_to_replacement
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
    # Spherical-harmonic Y_{2,0} on the unit sphere. Lon-independent
    # (∂_φ Y_{2,0} ≡ 0); ∂_θ component is the meaningful test signal for
    # the latlon centered-FD family (`grid_family=latlon`).
    #
    #   Y_{2,0}(θ,φ) = (1/4)·√(5/π)·(3·cos²θ − 1)
    #   ∂_θ Y_{2,0}  = (1/4)·√(5/π)·(−6·cos θ·sin θ) = −(3/2)·√(5/π)·cos θ·sin θ
    #
    # `ic(lat, lon)` returns the scalar field; `derivative(lat, lon)` returns
    # `(d_dlat, d_dlon)` so the 2d_latlon_sphere runner can pick the axis
    # the rule's `applies_to.dim` selects without reshaping the catalog.
    # Used by `centered_2nd_uniform_latlon`'s convergence fixture
    # (`mms_kind="Y_2_0_unit_sphere"`); the runner activates once
    # ESD/dsc-y0jj lifts stencil-only rules to ESS-replacement form.
    "Y_2_0_unit_sphere" => (
        ic = (lat, lon) -> 0.25 * sqrt(5 / π) * (3 * cos(lat)^2 - 1),
        derivative = (lat, lon) -> (
            -1.5 * sqrt(5 / π) * cos(lat) * sin(lat),  # ∂/∂lat
            0.0,                                        # ∂/∂lon (Y_{2,0} is lon-independent)
        ),
    ),
    # 2D periodic vector field for the Arakawa C-grid divergence MMS:
    #
    #   F(x,y) = (sin(2π x)·cos(2π y),  cos(2π x)·sin(2π y))   on [0,1]² periodic
    #   ∇·F    = ∂_x F_x + ∂_y F_y
    #          = 2π cos(2π x)·cos(2π y) + 2π cos(2π x)·cos(2π y)
    #          = 4π cos(2π x)·cos(2π y)
    #
    # Both components vary non-trivially in both axes so the sweep
    # exercises the face_x and face_y stencil arms of the C-grid stencil.
    # `ic(x, y)` returns the velocity tuple `(F_x, F_y)`; the
    # 2d_arakawa_periodic runner samples F_x at face_x edges and F_y at
    # face_y edges. `derivative(x, y)` returns the scalar divergence at
    # cell centers. Used by `divergence_arakawa_c`'s convergence fixture
    # (`mms_kind="vec_sincos_2d_periodic"`); the runner activates once
    # ESD/dsc-y0jj lifts stencil-only rules to ESS-replacement form
    # (today divergence_arakawa_c is stencil-only and routes to
    # `stencil_form_rule` before the arakawa branch).
    "vec_sincos_2d_periodic" => (
        ic = (x, y) -> (sin(2π * x) * cos(2π * y), cos(2π * x) * sin(2π * y)),
        derivative = (x, y) -> 4π * cos(2π * x) * cos(2π * y),
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
    # `1d_cartesian_periodic` (esd-0ip) — runner implemented; drives
    # `discretize → build_evaluator` via a per-cell-scalar ESM model built
    # from the rule's `replacement` AST (or lowered from stencil form).
    # Activated for `centered_2nd_uniform` (O(h²)) and `upwind_1st` (O(h)).
    "1d_cartesian_periodic",
    # `2d_latlon_sphere` (ESD/dsc-mps8) — runner scaffolded with full
    # canonical-pipeline shape in `_run_layer_b_2d_latlon_sphere`.
    # Wired in here per witness Path-C directive so the dispatcher
    # auto-routes `grid_family=latlon` rules to the runner the moment
    # `_layer_b_topology_key`'s stencil-form gate retires (today the
    # only `grid_family=latlon` rule, `centered_2nd_uniform_latlon`,
    # is stencil-only and routes to `stencil_form_rule` BEFORE the
    # latlon check, so this entry is forward-compat — no rule reaches
    # the runner until ESD/dsc-y0jj's lifter lands).
    "2d_latlon_sphere",
    # `2d_arakawa_periodic` (ESD/dsc-70zp) — runner scaffolded with full
    # canonical-pipeline shape in `_run_layer_b_2d_arakawa_periodic`.
    # Same forward-compat pattern as 2d_latlon_sphere: the only
    # `grid_family=arakawa` rule today (`divergence_arakawa_c`) is
    # stencil-only (no `replacement` AST) and routes to `stencil_form_rule`
    # BEFORE the arakawa branch, so this entry is forward-compat — no
    # rule reaches the runner until ESD/dsc-y0jj's lifter lands.
    "2d_arakawa_periodic",
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
- `"1d_vertical_column"` — single-axis but non-periodic vertical spacing.
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

    # `stencil`-form rules cannot ride the ESS rule engine until their
    # stencil entries are lowered to a `replacement` AST. Tracked
    # separately so a future bead can either author the lowering or
    # ship per-rule canonical fixtures with replacement-form rules.
    if haskey(spec, "stencil") && !haskey(spec, "replacement")
        return "stencil_form_rule"
    end

    grid_family = String(get(spec, "grid_family", ""))
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

Per-topology implementations live in dsc-kswm follow-up beads. Until the
matching family lands, this dispatcher is unreachable — `run_mms_convergence`
gates on `_LAYER_B_SUPPORTED_TOPOLOGIES` first.
"""
function _run_layer_b_canonical_grid(topology_key::AbstractString, rule::RuleFile, mms, n::Int)
    topology_key == "1d_cartesian_periodic" && return _run_layer_b_1d_cartesian_periodic(rule, mms, n)
    topology_key == "2d_latlon_sphere" && return _run_layer_b_2d_latlon_sphere(rule, mms, n)
    topology_key == "2d_arakawa_periodic" && return _run_layer_b_2d_arakawa_periodic(rule, mms, n)
    error("Layer-B canonical-pipeline runner for topology=$(topology_key) not implemented; " *
          "tracked in a dsc-kswm follow-up bead. The dispatcher should not have been reached " *
          "for this topology key.")
end

"""
    _run_layer_b_2d_latlon_sphere(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=latlon` rules
(scaffolded by ESD/dsc-mps8). Drives one resolution of the spherical-
harmonic convergence sweep through the canonical pipeline:

1. Build a parametric `.esm` document that wires `rule` against a
   `family=latlon` model whose grid is sized `2n × n` (lon × lat).
2. Run `EarthSciSerialization.discretize` to apply the rule and emit
   the canonical-form ArrayOp AST (`Earth-radius` + per-cell `cos_lat`
   bindings come from the ESS latlon grid metadata; the runner does
   not synthesize them).
3. Build an evaluator with `EarthSciSerialization.build_evaluator`
   (the documented official ESS Julia tree-walk runner per
   ESS/AGENTS.md §"Multiple official runners per binding are OK").
4. Set `u0` to `mms.ic(lat, lon)` sampled at interior cell centers
   (`u0` is laid out per `var_map`).
5. Compute `du = f!(u0, p, 0)` to get the rule's discrete output.
6. Return the L_inf error against `mms.derivative(lat, lon)` along the
   axis the rule's `applies_to.dim` selects, taking only interior cells
   (lat-pole rows excluded — boundary fix-up is the rule's BC
   contract, not the convergence-sweep concern). The lon axis is
   periodic and exact for `Y_{2,0}` (φ-independent), so it contributes
   `0` and never sets the L_inf — the lat axis carries the signal.

# Activation status (ESD/dsc-mps8 partial-land per witness Path C)

This runner is **scaffolded**; `_LAYER_B_SUPPORTED_TOPOLOGIES` does NOT
yet include `"2d_latlon_sphere"`, so `run_mms_convergence` routes
`grid_family=latlon` rules to a SKIP before this body executes.

The activation prerequisite is ESD/dsc-y0jj (stencil → replacement
lift). `centered_2nd_uniform_latlon` ships only a `stencil` spec (no
`replacement` AST); ESS `parse_rule` rejects stencil-only rules with
`E_RULE_REPLACEMENT_MISSING` (see `rule_engine.jl:699`), so step (2)
above cannot succeed today. Once dsc-y0jj's lifter lands:

- the `_layer_b_topology_key` stencil-form gate retires (or routes
  lifted rules through),
- `"2d_latlon_sphere"` is added to `_LAYER_B_SUPPORTED_TOPOLOGIES`,
- this `error(...)` body is replaced by the implementation sketched
  above, and
- `centered_2nd_uniform_latlon` moves out of `pending_canonical_layer_b`
  in `test_esd_walker.jl` (the 1-line follow-up the witness called).
"""
function _run_layer_b_2d_latlon_sphere(rule::RuleFile, mms, n::Int)
    error("_run_layer_b_2d_latlon_sphere: scaffolding only (ESD/dsc-mps8 partial-land). " *
          "Activation pends ESD/dsc-y0jj (stencil → replacement lift) since " *
          "centered_2nd_uniform_latlon — the only `grid_family=latlon` rule today — " *
          "ships only a `stencil` spec and ESS `parse_rule` rejects stencil-only rules " *
          "with E_RULE_REPLACEMENT_MISSING. If this error fires, the dispatcher's " *
          "`_LAYER_B_SUPPORTED_TOPOLOGIES` was advanced past this scaffolding without " *
          "the runner body being authored — see the docstring for the canonical-pipeline " *
          "shape (build .esm → discretize → build_evaluator → L_inf vs analytic).")
end

"""
    _run_layer_b_2d_arakawa_periodic(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=arakawa` rules
(scaffolded by ESD/dsc-70zp). Drives one resolution of the 2D periodic
divergence convergence sweep through the canonical pipeline:

1. Build a parametric `.esm` document that wires `rule` against a 2D
   periodic Cartesian model (`family=arakawa`, stagger=C) with grid
   sized `n × n` over `[0, 1]²`.
2. Run `EarthSciSerialization.discretize` to apply the rule and emit
   the canonical-form ArrayOp AST. The C-grid metadata (face-x and
   face-y face positions) comes from the ESS arakawa grid; the runner
   does not synthesize face coordinates locally.
3. Build an evaluator with `EarthSciSerialization.build_evaluator`
   (the documented official ESS Julia tree-walk runner per
   ESS/AGENTS.md §"Multiple official runners per binding are OK").
4. Sample `mms.ic(x, y) = (F_x, F_y)` at face_x edges (i = 0..n-1) and
   face_y edges (j = 0..n-1) respectively to populate `u0` per
   `var_map`.
5. Compute `du = f!(u0, p, 0)` to obtain the discrete divergence at
   cell centers.
6. Return the L_inf error against `mms.derivative(x, y)` evaluated at
   cell centers `((i + 0.5)/n, (j + 0.5)/n)`. Both axes wrap
   periodically per the C-grid convention; no boundary cells excluded
   (the grid is fully periodic).

The expected order is O(h²) — the `expected_min_order=1.9` declared
by `divergence_arakawa_c/fixtures/convergence/expected.esm` leaves a
0.1 margin for pre-asymptotic drift over the n=16→32→64→128 sweep.

# Activation status (ESD/dsc-70zp partial-land per witness Path C)

This runner is **scaffolded**; `_LAYER_B_SUPPORTED_TOPOLOGIES` already
contains `"2d_arakawa_periodic"` (forward-compat), but the dispatcher's
`stencil_form_rule` gate routes `divergence_arakawa_c` to a SKIP before
this body executes. The activation prerequisite is ESD/dsc-y0jj
(stencil → replacement lift). `divergence_arakawa_c` ships only a
`stencil` spec (4-row Arakawa-C: face_x ±1·1/dx, face_y ±1·1/dy, see
`discretizations/finite_volume/divergence_arakawa_c.json`) and ESS
`parse_rule` rejects stencil-only rules with `E_RULE_REPLACEMENT_MISSING`
(`rule_engine.jl:695-699`), so step (2) above cannot succeed today.
Once dsc-y0jj's lifter lands:

- the `_layer_b_topology_key` stencil-form gate retires (or routes
  lifted rules through),
- this `error(...)` body is replaced by the implementation sketched
  above, and
- `divergence_arakawa_c` gains its first Layer-B PASS (the
  `layer_b_passes` assertion in `test_esd_walker.jl` is bumped in the
  same follow-up).
"""
function _run_layer_b_2d_arakawa_periodic(rule::RuleFile, mms, n::Int)
    error("_run_layer_b_2d_arakawa_periodic: scaffolding only (ESD/dsc-70zp partial-land). " *
          "Activation pends ESD/dsc-y0jj (stencil → replacement lift) since " *
          "divergence_arakawa_c — the only `grid_family=arakawa` rule today — " *
          "ships only a `stencil` spec and ESS `parse_rule` rejects stencil-only rules " *
          "with E_RULE_REPLACEMENT_MISSING. If this error fires, the dispatcher's " *
          "`_layer_b_topology_key` stencil-form gate was retired without this runner " *
          "body being authored — see the docstring for the canonical-pipeline shape " *
          "(build .esm → discretize → build_evaluator → L_inf vs analytic, sampling " *
          "F_x at face_x edges and F_y at face_y edges of the n×n periodic grid).")
end

"""
    _run_layer_b_1d_cartesian_periodic(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=cartesian` rules with a
single periodic dimension (esd-0ip). Drives one resolution of the MMS
convergence sweep through the canonical pipeline:

1. Load the rule's replacement AST, lowering from stencil form if necessary
   via `lower_stencil_to_replacement`. Strip the `arrayop` wrapper if present.
2. Build a per-cell-scalar ESM model with state variables `u_1, ..., u_n` and
   parameter `dx`. Each cell equation is constructed by substituting the rule's
   replacement AST at that cell index: `index(\$u, \$x+off)` → `"u_j"` with
   1-based periodic wrap.
3. Run `EarthSciSerialization.discretize` (canonicalization; equations are
   already in explicit scalar form, so no rule rewriting occurs).
4. Run `EarthSciSerialization.build_evaluator` to produce the ODE RHS `f!`.
5. Set each `u_i = mms.ic((i-0.5)/n)`.
6. Evaluate `f!(du, u0, p, 0)` once to obtain the discrete derivative.
7. Return the L_inf error against `mms.derivative((i-0.5)/n)`.

The per-cell model builder (step 2) is a generic AST → ESM translator — no
math is computed outside the official ESS tree-walk in step 4. Works for any
1D periodic Cartesian rule that carries a `replacement` AST (either authored
directly or producible via `lower_stencil_to_replacement`).
"""
function _run_layer_b_1d_cartesian_periodic(rule::RuleFile, mms, n::Int)
    # 1. Load and lower the rule spec to replacement form.
    rule_doc = JSON.parse(read(rule.path, String))
    spec = get(get(rule_doc, "discretizations", Dict()), rule.name, nothing)
    spec isa AbstractDict || error("rule spec missing for $(rule.name)")
    lowered = lower_stencil_to_replacement(spec)
    repl = lowered["replacement"]
    # Strip the arrayop wrapper (centered_2nd_uniform wraps in arrayop; the
    # scalar expression lives in the `expr` field).
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl

    # 2. Build a per-cell-scalar ESM model on the unit interval [0, 1].
    dx = 1.0 / n
    cell_x(i) = (i - 0.5) * dx

    variables = Dict{String, Any}()
    for i in 1:n
        variables["u_$(i)"] = Dict{String, Any}(
            "type"    => "state",
            "default" => mms.ic(cell_x(i)),
            "units"   => "1",
        )
    end
    variables["dx"] = Dict{String, Any}(
        "type"    => "parameter",
        "default" => dx,
        "units"   => "1",
    )

    equations = Any[]
    for i in 1:n
        rhs = _layer_b_1d_build_cell_expr(expr, i, n)
        push!(equations, Dict{String, Any}(
            "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u_$(i)"], "wrt" => "t"),
            "rhs" => rhs,
        ))
    end

    esm = Dict{String, Any}(
        "esm"      => "0.2.0",
        "metadata" => Dict{String, Any}("name" => "layer_b_1d_cartesian_periodic_n$(n)"),
        "models"   => Dict{String, Any}(
            "M" => Dict{String, Any}(
                "variables" => variables,
                "equations" => equations,
            ),
        ),
    )

    # 3-4. Canonical pipeline: discretize (canonicalization) → build_evaluator.
    disc = EarthSciSerialization.discretize(esm)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

    # 5. Set the MMS initial condition (overrides the defaults already set above).
    for i in 1:n
        u0[var_map["u_$(i)"]] = mms.ic(cell_x(i))
    end

    # 6. Evaluate the discrete operator once (du/dt = rule(u) at t=0).
    du = similar(u0)
    f!(du, u0, p, 0.0)

    # 7. Return L_inf error vs analytic derivative at cell centers.
    return maximum(abs(du[var_map["u_$(i)"]] - mms.derivative(cell_x(i))) for i in 1:n)
end

# Translate one node of the rule's replacement AST to a scalar ESM expression
# for cell `i` in a grid of `n` cells. Resolves `index(\$u, ...)` ops to
# per-cell variable names `"u_j"` with 1-based periodic wrapping. All other
# ops and literals (including `"dx"`) pass through unchanged for ESS tree-walk.
function _layer_b_1d_build_cell_expr(expr, i::Int, n::Int)
    expr isa Number && return expr
    if expr isa AbstractString
        return String(expr)
    end
    expr isa AbstractDict || error("unexpected node type in rule expression: $(typeof(expr))")
    op   = String(expr["op"])
    args = expr["args"]
    if op == "index"
        idx = _layer_b_1d_eval_index(args[2], i)
        j   = mod(idx - 1, n) + 1   # 1-based periodic wrap
        return "u_$(j)"
    end
    new_args = Any[_layer_b_1d_build_cell_expr(a, i, n) for a in args]
    return Dict{String, Any}("op" => op, "args" => new_args)
end

# Evaluate a rule's index sub-expression to a concrete integer cell offset.
# Recognises `"\$x"` (current-cell pattern variable) and integer offsets.
function _layer_b_1d_eval_index(expr, i::Int)::Int
    expr isa AbstractString && String(expr) == "\$x" && return i
    expr isa Integer && return Int(expr)
    expr isa AbstractFloat && return Int(expr)
    expr isa AbstractDict || error("cannot evaluate index expression: $expr")
    op   = String(expr["op"])
    args = expr["args"]
    a    = _layer_b_1d_eval_index(args[1], i)
    b    = _layer_b_1d_eval_index(args[2], i)
    op == "+" && return a + b
    op == "-" && return a - b
    error("unsupported index op in rule expression: $(op)")
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
