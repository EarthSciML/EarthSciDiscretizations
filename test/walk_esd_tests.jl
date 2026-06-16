module WalkESDTests

using EarthSciDiscretizations: load_rules, RuleFile, eval_coeff,
    lower_stencil_to_scheme,
    CubedSphereGrid, extend_with_ghosts, gnomonic_metric, metric_eval,
    build_ode_problem, build_duo_grid, build_mpas_grid
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
`applicable:false` markers, drives the per-topology canonical-pipeline
runner where one exists, and otherwise SKIPs naming the missing piece
(see `run_mms_convergence` docstring).
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

Layer B: read the convergence fixture, honor `applicable:false` markers,
and drive the family's runner where one exists; otherwise SKIP pending the
per-topology canonical-pipeline runner.

History: ESS retired `verify_mms_convergence` + `mms_evaluator` in esm-4t5
(the 2026-04-29 single-evaluation-pathway directive) — the old harness was
a shadow evaluator outside the canonical `discretize → ArrayOp → eval`
pipeline. The replacement lands per topology family under dsc-kswm. The
ESS side of the pipeline has since landed for the cartesian foundation
(scheme expansion esm-j1u, arrayop lift + periodic folding, tree-walk
arrayop evaluation): `1d_cartesian_periodic` drives it ArrayOp-natively.
The remaining families' runners still scalarize per cell in this file and
pend either an ESS extension (scheme `grid_family`/selector dispatch beyond
cartesian — cubed-sphere panel ESS/esm-57f, unstructured ESS/esm-bpr) or a
runner port. Fixtures in a family without a runner SKIP with the unified
`_LAYER_B_PIPELINE_PENDING` reason, suffixed with the topology key so the
report names the missing piece.
"""
const _LAYER_B_PIPELINE_PENDING = "Layer-B awaits canonical-pipeline replacement: per-topology runners drive ESS discretize → ArrayOp → eval (tracked at ESD/dsc-kswm); families beyond the implemented set pend ESS scheme dispatch beyond the cartesian foundation (cubed-sphere ESS/esm-57f, unstructured ESS/esm-bpr) or a per-family runner port"

# The ArrayOp-native 1d_cartesian_periodic runner needs the ESS
# `lift_1d_arrayop` discretize kwarg (1D shaped variables otherwise keep the
# bare D(u) LHS, which the tree-walk evaluator cannot expand). Detected at
# load so an older resolved ESS yields a structured SKIP, not a FAIL.
const _ESS_SUPPORTS_LIFT_1D = hasmethod(
    EarthSciSerialization.discretize, Tuple{AbstractDict}, (:lift_1d_arrayop,))

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
    # (`mms_kind="Y_2_0_unit_sphere"`).
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
    # (`mms_kind="vec_sincos_2d_periodic"`).
    "vec_sincos_2d_periodic" => (
        ic = (x, y) -> (sin(2π * x) * cos(2π * y), cos(2π * x) * sin(2π * y)),
        derivative = (x, y) -> 4π * cos(2π * x) * cos(2π * y),
    ),
    # 2D periodic product sine for the cartesian 5-point Laplacian MMS:
    #
    #   u(x,y) = sin(2π x)·sin(2π y)   on [0,1]² periodic
    #   ∇²u    = −2·(2π)²·sin(2π x)·sin(2π y) = −8π²·u
    #
    # `ic(x, y)` returns the scalar field at cell centers; `derivative(x, y)`
    # returns the analytic Laplacian. Used by `laplacian_2nd_uniform_cartesian`'s
    # convergence fixture (`mms_kind="sin2pix_sin2piy_periodic"`).
    "sin2pix_sin2piy_periodic" => (
        ic = (x, y) -> sin(2π * x) * sin(2π * y),
        derivative = (x, y) -> -8 * π^2 * sin(2π * x) * sin(2π * y),
    ),
    # 1D periodic sine in CELL-AVERAGE form, for FV reconstruction rules
    # (`fv_cell_average_1d` topology). Unlike the point-sample kinds above,
    # `ic(xl, xr)` returns the EXACT average of sin(2πx) over the cell
    # [xl, xr]:
    #
    #   (1/(xr−xl)) ∫ sin(2πx) dx = (cos(2π xl) − cos(2π xr)) / (2π (xr−xl))
    #
    # and `derivative(x)` is the point VALUE sin(2πx) at the queried
    # coordinate — reconstruction rules estimate point values (e.g. cell-edge
    # values) from cell averages, so the analytic target is the field itself,
    # evaluated wherever the output's stencil support says it lives. Used by
    # `ppm_reconstruction`'s convergence fixture
    # (`mms_kind="sin_2pi_x_cell_average"`).
    "sin_2pi_x_cell_average" => (
        ic = (xl, xr) -> (cos(2π * xl) - cos(2π * xr)) / (2π * (xr - xl)),
        derivative = x -> sin(2π * x),
    ),
    # 1D periodic advection at unit speed, for flux-form advection rules
    # whose pattern carries an auxiliary velocity field (e.g.
    # `div($U * $q, dim=$x)`): u(x) = sin(2π x + 1) transported by U ≡ 1,
    # so the divergence d(U·u)/dx = 2π cos(2π x + 1). The phase shift
    # matches the fixture's declared `phase_shifted_sine` profile and
    # keeps the field asymmetric w.r.t. the grid so no error component
    # cancels by symmetry. The auxiliary binding lives in
    # `_LAYER_B_MMS_AUX` under the same kind. Used by `weno5_advection`'s
    # convergence fixture (`mms_kind="sin_2pi_x_unit_advection"`).
    "sin_2pi_x_unit_advection" => (
        ic = x -> sin(2π * x + 1.0),
        derivative = x -> 2π * cos(2π * x + 1.0),
    ),
    # 1D periodic sine, SECOND derivative: u(x) = sin(2π x);
    # u''(x) = −(2π)²·sin(2π x). Used by `centered_2nd_deriv_uniform`'s
    # convergence fixture (`mms_kind="sin_2pi_x_second_derivative"`).
    "sin_2pi_x_second_derivative" => (
        ic = x -> sin(2π * x),
        derivative = x -> -4π^2 * sin(2π * x),
    ),
    # 2D periodic product sine, MIXED second derivative:
    # u(x,y) = sin(2π x)·sin(2π y); ∂²u/∂x∂y = 4π²·cos(2π x)·cos(2π y).
    # Used by `mixed_deriv_2nd_uniform`'s convergence fixture
    # (`mms_kind="sin2pix_sin2piy_mixed_deriv"`).
    "sin2pix_sin2piy_mixed_deriv" => (
        ic = (x, y) -> sin(2π * x) * sin(2π * y),
        derivative = (x, y) -> 4π^2 * cos(2π * x) * cos(2π * y),
    ),
    # 1D periodic nonlinear diffusion: u(x) = sin(2π x) with the
    # space-varying coefficient f(x) = 2 + sin(2π x) (strictly positive),
    # so with g = 2πx:
    #
    #   d/dx( f·du/dx ) = f'·u' + f·u''
    #                   = (2π cos g)(2π cos g) + (2 + sin g)(−4π² sin g)
    #                   = 4π²·(cos²g − 2·sin g − sin²g)
    #
    # The coefficient binding lives in `_LAYER_B_MMS_AUX` under the same
    # kind. Used by `nonlinear_laplacian_uniform`'s convergence fixture
    # (`mms_kind="sin_2pi_x_nonlinear_diffusion"`).
    "sin_2pi_x_nonlinear_diffusion" => (
        ic = x -> sin(2π * x),
        derivative = x -> 4π^2 * (cos(2π * x)^2 - 2 * sin(2π * x) - sin(2π * x)^2),
    ),
    # Unstructured-sphere manufactured solution for nn_diffusion_mpas / nn_diffusion_duo.
    #
    #   u(lon, lat) = sin(lon)·cos(lat)   (spherical harmonic Y₁₁, l=1)
    #   ΔΩ u        = −2·sin(lon)·cos(lat) = −2·u   (unit-sphere Laplacian eigenvalue −l(l+1)=−2)
    #
    # The physical Laplacian on a sphere of radius R is ΔΩu / R². The runner
    # divides by R² extracted from the grid, so this catalog entry stays R-independent.
    # ic(lon, lat) → field value; derivative(lon, lat) → unit-sphere Laplacian
    # (the runner applies 1/R² scaling). Used by nn_diffusion_{mpas,duo} convergence
    # fixtures (mms_kind="sin_lon_cos_lat").
    "sin_lon_cos_lat" => (
        ic = (lon, lat) -> sin(lon) * cos(lat),
        derivative = (lon, lat) -> -2.0 * sin(lon) * cos(lat),
    ),
    # Cubed-sphere manufactured solution: φ(ξ,η) = cos(2ξ)·cos(2η) on all 6 panels.
    # The analytic covariant Laplacian ∇²φ = (1/J)[∂_ξ(Jg^{ξξ}∂_ξφ + Jg^{ξη}∂_ηφ)
    #                                              + ∂_η(Jg^{ξη}∂_ξφ + Jg^{ηη}∂_ηφ)]
    # is computed at (ξ,η) via gnomonic_metric + finite-differenced Jg^{ab} derivatives
    # (same approach as test/test_laplacian.jl, which served as the oracle for O(h²)
    # verification at Nc ∈ {8,16,32}).
    #
    # `ic(ξ,η)` returns the scalar field value; `derivative(ξ,η)` returns the
    # analytic covariant Laplacian — panel-independent because the gnomonic metric
    # depends only on (ξ,η), not on the panel index. Used by
    # `covariant_laplacian_cubed_sphere` convergence fixture
    # (`mms_kind="cos2xi_cos2eta_cubed_sphere"`).
    # NOTE: pre-asymptotic behavior for Nc<32 due to gnomonic metric variation
    # near panel corners; fixture uses Nc ∈ {64,128,256} for min order ≥ 1.8.
    "cos2xi_cos2eta_cubed_sphere" => (
        ic = (ξ, η) -> cos(2ξ) * cos(2η),
        derivative = (ξ, η) -> begin
            h = 1.0e-6
            df_dξ = -2 * sin(2ξ) * cos(2η)
            df_dη = -2 * cos(2ξ) * sin(2η)
            d2f_dξ2 = -4 * cos(2ξ) * cos(2η)
            d2f_dη2 = -4 * cos(2ξ) * cos(2η)
            d2f_dξdη = 4 * sin(2ξ) * sin(2η)
            J, g_ξξ, g_ηη, g_ξη = gnomonic_metric(ξ, η, 1.0)
            det_g = g_ξξ * g_ηη - g_ξη^2
            ginv_ξξ = g_ηη / det_g
            ginv_ηη = g_ξξ / det_g
            ginv_ξη = -g_ξη / det_g
            function Jginv_at(ξv, ηv)
                Jv, gxx, gee, gxe = gnomonic_metric(ξv, ηv, 1.0)
                d = gxx * gee - gxe^2
                return (Jv * gee / d, Jv * gxx / d, Jv * (-gxe / d))
            end
            Jgxx_p, _, Jgxe_p = Jginv_at(ξ + h, η)
            Jgxx_m, _, Jgxe_m = Jginv_at(ξ - h, η)
            _, Jgyy_np, Jgxe_np = Jginv_at(ξ, η + h)
            _, Jgyy_nm, Jgxe_nm = Jginv_at(ξ, η - h)
            dJgxx_dξ = (Jgxx_p - Jgxx_m) / (2h)
            dJgyy_dη = (Jgyy_np - Jgyy_nm) / (2h)
            dJgxe_dξ = (Jgxe_p - Jgxe_m) / (2h)
            dJgxe_dη = (Jgxe_np - Jgxe_nm) / (2h)
            ginv_ξξ * d2f_dξ2 + ginv_ηη * d2f_dη2 +
                (1 / J) * (dJgxx_dξ * df_dξ + dJgyy_dη * df_dη) +
                2 * ginv_ξη * d2f_dξdη +
                (1 / J) * (dJgxe_dξ * df_dη + dJgxe_dη * df_dξ)
        end,
    ),
)

# ---------------------------------------------------------------------------
# Auxiliary-field table, keyed by `mms_kind`. Rules whose `applies_to`
# pattern carries pattern variables beyond the transported state (e.g. the
# velocity `$U` in `div($U * $q)`) need concrete document fields to bind
# them to. Each entry maps a STRIPPED pattern-variable name (`"U"` for
# `$U`) to the field's value — a constant, or a function of the cell-center
# coordinate for space-varying fields. Runners declare one shaped state per
# entry the rule's pattern actually uses, with a zero tendency equation,
# and initialize its cells from the value. Pattern variables NOT named here
# bind to the transported state `u`.
# ---------------------------------------------------------------------------
const _LAYER_B_MMS_AUX = Dict{String, Dict{String, Any}}(
    "sin_2pi_x_unit_advection"      => Dict{String, Any}("U" => 1.0),
    "sin_2pi_x_nonlinear_diffusion" => Dict{String, Any}("f" => (x -> 2 + sin(2π * x))),
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
    # `1d_cartesian_periodic` (esd-0ip; ArrayOp-native under dsc-kswm) —
    # runner implemented; builds a PDE ESM document (array state + grad op)
    # and drives `discretize(lift_1d_arrayop=true) → ArrayOp → build_evaluator`
    # with the rule supplied as an ESS scheme + `use:` rule (stencil form,
    # via `lower_stencil_to_scheme`) or a pattern+replacement document rule
    # (replacement form). No per-cell scalarization.
    # Activated for `centered_2nd_uniform` (O(h²)) and `upwind_1st` (O(h)).
    "1d_cartesian_periodic",
    # `1d_vertical_column` (esd-bbp) — runner implemented; reuses the
    # 1D Cartesian per-cell-scalar machinery with parameter "h" (spacing).
    # Activated for `centered_2nd_uniform_vertical` (O(h²)).
    "1d_vertical_column",
    # `2d_latlon_sphere` (esd-bbp) — runner implemented; drives
    # `discretize → build_evaluator` via a per-cell-scalar ESM model on
    # an n×n lat-lon grid, substituting per-row cos_lat parameters.
    # Activated for `centered_2nd_uniform_latlon` (O(h²) on lat axis).
    "2d_latlon_sphere",
    # `2d_arakawa_periodic` (esd-bbp) — runner implemented; evaluates the
    # stencil's coeff ASTs via the canonical pipeline with face values from
    # the MMS IC embedded as literals. Activated for `divergence_arakawa_c`
    # (O(h²) divergence test).
    "2d_arakawa_periodic",
    # `2d_cartesian_periodic` (dsc-vst2) — ArrayOp-native runner; builds a
    # 2D periodic PDE ESM document and drives `discretize → ArrayOp →
    # build_evaluator` with the rule's replacement AST (canonical `$target`
    # component names i/j per RFC §7.1) as a document rule. No 1D-lift
    # kwarg needed — multidimensional equations lift by default. Activated
    # for `laplacian_2nd_uniform_cartesian` (O(h²) 5-point Laplacian).
    "2d_cartesian_periodic",
    # `fv_cell_average_1d` (dsc-a7b2) — ArrayOp-native runner for FV
    # reconstruction rules: the IC is the EXACT cell average of the MMS
    # field, and each named stencil output (multi-output stencil objects
    # lower one output per ESM document via
    # `lower_stencil_to_scheme(...; output=...)`) is compared against the
    # point value at the location its offset-support midpoint implies
    # (±1/2 → cell edge). Activated for `ppm_reconstruction` (CW84
    # unlimited edge interpolation, O(h³)+).
    "fv_cell_average_1d",
    # `cubed_sphere_cross_metric` (esd-ecq) — direct-evaluation runner for
    # cross_metric composite rules on the gnomonic cubed sphere. Evaluates
    # each term of the cross_metric rule by binding per-cell metric arrays
    # (ginv_xi_xi, ginv_eta_eta, ginv_xi_eta, J, dJgxx_dxi, dJgyy_deta,
    # dJgxe_dxi, dJgxe_deta) from CubedSphereGrid.metric_eval and running
    # eval_coeff on each stencil entry. Ghost cells are filled via
    # extend_with_ghosts (panel connectivity aware). Analytic reference from
    # gnomonic_metric + FD metric derivatives. Activated for
    # `covariant_laplacian_cubed_sphere` (O(h²) full covariant Laplacian).
    # The ESS canonical pipeline extension (discretize → build_evaluator with
    # cubed_sphere selector dispatch) is tracked at ESS/ess-cubed-sphere-pipeline.
    "cubed_sphere_cross_metric",
    # `unstructured_ode` (esd-cal) — ODE-pipeline runner for unstructured
    # sphere rules (grid_family=unstructured: nn_diffusion_mpas, nn_diffusion_duo).
    # Loads the MPAS Voronoi or DUO icosahedral mesh via builtin paths,
    # injects per-cell MMS ICs, calls build_ode_problem, evaluates the RHS
    # at t=0, and measures L∞ error vs the analytic Laplacian. Convergence
    # order is computed across the builtin mesh ladder.
    "unstructured_ode",
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
    # dsc-y0jj (the stencil → replacement lowerer) landed; this key fires for
    # stencil rules that lower_stencil_to_replacement rejects for non-unstructured
    # families (not yet classified). Unstructured rules short-circuit to
    # "unstructured_ode" before reaching the lowering gate.
    "stencil_form_rule"           => "ESS/esm-bpr (non-unstructured stencil lowering pending)",
    "cubed_sphere_cross_metric"   => "ESD/esd-ecq",
    "unstructured_ode"            => "ESD/esd-cal",
    "unsupported"                 => "no follow-up bead — out of Layer-B scope",
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
    if topology_key in ("1d_cartesian_periodic", "fv_cell_average_1d") && !_ESS_SUPPORTS_LIFT_1D
        return LayerResult(LAYER_SKIP,
            "$(topology_key) runner requires EarthSciSerialization.discretize(...; " *
            "lift_1d_arrayop=true); the resolved ESS version predates it — update ESS to " *
            "enable the ArrayOp-native Layer-B path")
    end

    # Generic MMS lookup. Only fixtures that declare a registered `mms_kind`
    # can be driven; bare descriptive `manufactured_solution` strings SKIP
    # until the fixture is upgraded.
    mms_kind = String(get(input_json, "mms_kind", ""))
    mms_base = get(_LAYER_B_MMS_CATALOG, mms_kind, nothing)
    if mms_base === nothing
        return LayerResult(LAYER_SKIP, _LAYER_B_PIPELINE_PENDING * " (mms_kind=$(repr(mms_kind)) not registered in Layer-B MMS catalog)")
    end
    # Attach the kind's auxiliary-field bindings (empty for most kinds).
    mms = (; mms_base..., aux = get(_LAYER_B_MMS_AUX, mms_kind, Dict{String, Any}()))

    # Unstructured ODE topology: reads grid_refs (GDD paths), not {n:} dicts.
    # Delegate to a separate sweep function that builds the ODE problem directly.
    if topology_key == "unstructured_ode"
        return _run_layer_b_unstructured_ode_sweep(rule, mms, input_json,
                                                    convergence_dir, expected_min_order)
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

    # Unstructured rules (grid_family=unstructured) use an ODE-pipeline runner
    # that calls build_ode_problem directly and does NOT route through the stencil
    # lowerer. Check grid_family BEFORE the lowering gate so these rules do not
    # fall through to "stencil_form_rule".
    if String(get(spec, "grid_family", "")) == "unstructured"
        return "unstructured_ode"
    end

    # `stencil`-form rules without a `replacement` AST: after esd-t4h all
    # live latlon/vertical rules carry authored `replacement` fields.  The only
    # remaining stencil-without-replacement rules are cubed_sphere spec-only
    # rules (not yet materialised; EINSUM-8 territory).  Multi-output stencil
    # OBJECTS (keyed by output name, e.g. ppm_reconstruction) are drivable one
    # output at a time by the fv_cell_average_1d runner when the fixture
    # declares cell-average sampling; otherwise they remain `stencil_form_rule`.
    if haskey(spec, "stencil") && !haskey(spec, "replacement")
        if spec["stencil"] isa AbstractDict
            sampling_early = String(get(input_json, "sampling", "cell_center"))
            sampling_early == "cell_average" || return "stencil_form_rule"
        else
            # esd-t4h: lower_stencil_to_replacement retired; only cubed_sphere
            # spec-only stencil rules remain here — skip them for now.
            return "stencil_form_rule"
        end
    end

    grid_family = String(get(spec, "grid_family", ""))
    applies_to = get(spec, "applies_to", nothing)
    applies_to isa AbstractDict || return "unsupported"
    args = get(applies_to, "args", Any[])
    n_args = length(args)
    sampling = String(get(input_json, "sampling", "cell_center"))

    # Cubed-sphere cross_metric composite rules (esd-ecq): evaluated via a
    # direct per-cell metric-binding runner (no ESS discretize needed).
    grid_family == "cubed_sphere" && return "cubed_sphere_cross_metric"

    # Lat-lon sphere is its own topology family, regardless of arg count
    # (Y_l_m manufactured solutions need spherical geometry handling).
    grid_family == "latlon" && return "2d_latlon_sphere"

    # Arakawa C-grid (staggered face-centered fluxes on 2D Cartesian) gets
    # its own family because the canonical-pipeline runner has to handle
    # face-staggered velocity fields, distinct from the cell-centered 2D
    # Cartesian path.
    grid_family == "arakawa" && return "2d_arakawa_periodic"

    # Vertical column: single-axis rules on a vertical grid. The runner
    # reuses the 1D Cartesian per-cell-scalar machinery.
    grid_family == "vertical" && return "1d_vertical_column"

    # FV cell-average rules (PPM / WENO reconstruction) need their own
    # IC sampling (cell averages, not point samples) and exact-derivative
    # comparison, distinct from FD point-sample rules.
    if sampling == "cell_average"
        return "fv_cell_average_1d"
    end

    # 2D Cartesian rules: diagnose by the convergence fixture's `axis`
    # field (a 2D MMS sweep over a single axis at a time, as used by
    # weno5_advection_2d), or by the rule's stencil spanning two distinct
    # selector axes (e.g. the 5-point Laplacian's $x/$y entries). The
    # classifier deliberately does NOT pattern-match on the
    # manufactured_solution description string — that field is free-form
    # prose, not machine-readable.
    if grid_family == "cartesian" && haskey(input_json, "axis")
        return "2d_cartesian_periodic"
    end
    if grid_family == "cartesian" && _stencil_unique_axis_count(spec) >= 2
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

# Count the distinct selector axes a rule's stencil spans (generic schema
# attribute — used by the topology classifier to distinguish 1D from 2D
# cartesian rules). Entries without a `selector.axis` (e.g. cubed_sphere's
# plural `selectors`) contribute nothing; rules without a stencil count 0.
function _stencil_unique_axis_count(spec::AbstractDict)::Int
    stencil = get(spec, "stencil", nothing)
    stencil isa AbstractVector || return 0
    axes = Set{String}()
    for entry in stencil
        entry isa AbstractDict || continue
        selector = get(entry, "selector", nothing)
        selector isa AbstractDict || continue
        axis = get(selector, "axis", nothing)
        axis isa AbstractString && push!(axes, String(axis))
    end
    return length(axes)
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
    topology_key == "1d_cartesian_periodic"      && return _run_layer_b_1d_cartesian_periodic(rule, mms, n)
    topology_key == "2d_cartesian_periodic"      && return _run_layer_b_2d_cartesian_periodic(rule, mms, n)
    topology_key == "fv_cell_average_1d"         && return _run_layer_b_fv_cell_average_1d(rule, mms, n)
    topology_key == "1d_vertical_column"         && return _run_layer_b_1d_vertical_column(rule, mms, n)
    topology_key == "2d_latlon_sphere"           && return _run_layer_b_2d_latlon_sphere(rule, mms, n)
    topology_key == "2d_arakawa_periodic"        && return _run_layer_b_2d_arakawa_periodic(rule, mms, n)
    topology_key == "cubed_sphere_cross_metric"  && return _run_layer_b_cubed_sphere_cross_metric(rule, mms, n)
    error("Layer-B canonical-pipeline runner for topology=$(topology_key) not implemented; " *
          "tracked in a dsc-kswm follow-up bead. The dispatcher should not have been reached " *
          "for this topology key.")
end

# ---------------------------------------------------------------------------
# Unstructured ODE runner (esd-cal): nn_diffusion_mpas / nn_diffusion_duo
# ---------------------------------------------------------------------------

function _run_layer_b_unstructured_ode_sweep(rule::RuleFile, mms, input_json::AbstractDict,
                                              convergence_dir::AbstractString,
                                              expected_min_order::Float64)
    esm_rel = get(input_json, "esm", nothing)
    esm_rel === nothing &&
        return LayerResult(LAYER_FAIL, "convergence/input.esm missing 'esm' field for unstructured_ode runner")
    esm_path = joinpath(convergence_dir, String(esm_rel))
    isfile(esm_path) ||
        return LayerResult(LAYER_FAIL, "unstructured_ode runner: esm path not found: $esm_path")

    grid_refs_raw = get(input_json, "grid_refs", Any[])
    isempty(grid_refs_raw) &&
        return LayerResult(LAYER_FAIL, "convergence/input.esm declares no grid_refs for unstructured_ode")
    length(grid_refs_raw) < 2 &&
        return LayerResult(LAYER_FAIL, "unstructured_ode runner needs ≥ 2 grid_refs for convergence order (got $(length(grid_refs_raw)))")

    t_final = Float64(get(input_json, "t_final", 0.0))

    errors = Float64[]
    for gdd_rel in grid_refs_raw
        gdd_path = let p = String(gdd_rel)
            isabspath(p) ? p : joinpath(convergence_dir, p)
        end
        isfile(gdd_path) ||
            return LayerResult(LAYER_FAIL, "unstructured_ode runner: grid_ref not found: $gdd_path")
        err = try
            _run_layer_b_unstructured_ode(mms, esm_path, gdd_path)
        catch e
            return LayerResult(LAYER_FAIL,
                "unstructured ODE pipeline threw at grid=$(basename(gdd_path)): $(sprint(showerror, e))")
        end
        if !(err isa Real) || !isfinite(err)
            return LayerResult(LAYER_FAIL,
                "unstructured ODE runner at grid=$(basename(gdd_path)) returned non-finite: $(typeof(err))")
        end
        push!(errors, Float64(err))
    end

    any(iszero, errors) &&
        return LayerResult(LAYER_FAIL, "zero error in unstructured convergence sweep (degenerate fixture): $errors")

    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    min_order = minimum(orders)
    if min_order >= expected_min_order
        return LayerResult(
            LAYER_PASS,
            "min order $(round(min_order; digits = 2)) >= expected $(expected_min_order) " *
                "over $(length(grid_refs_raw)) grids (orders=$(round.(orders; digits = 2)))",
        )
    else
        return LayerResult(
            LAYER_FAIL,
            "min order $(round(min_order; digits = 2)) below expected $(expected_min_order) " *
                "(orders=$(round.(orders; digits = 2)), errors=$(round.(errors; sigdigits = 3)))",
        )
    end
end

function _run_layer_b_unstructured_ode(mms, esm_path::AbstractString, gdd_path::AbstractString)
    gdd = JSON.parse(read(gdd_path, String))
    domain_spec = get(get(gdd, "grids", Dict{String,Any}()), "domain", Dict{String,Any}())
    family = String(get(domain_spec, "family", ""))

    local lon_c::Vector{Float64}
    local lat_c::Vector{Float64}
    local R_val::Float64

    if family == "mpas"
        loader_spec = get(domain_spec, "loader", nothing)
        loader_spec === nothing && error("MPAS GDD missing loader spec")
        loader_path = String(get(loader_spec, "path", ""))
        startswith(loader_path, "builtin://") || error(
            "unstructured_ode Layer-B runner requires builtin:// MPAS path; got $loader_path")
        grid = build_mpas_grid(
            loader = Dict("path" => loader_path, "reader" => "auto",
                          "check" => String(get(loader_spec, "check", "strict"))),
        )
        lon_c = grid.mesh.lon_cell
        lat_c = grid.mesh.lat_cell
        R_val = Float64(grid.R)
    elseif family == "duo"
        loader_spec = get(domain_spec, "loader", nothing)
        loader_spec === nothing && error("DUO GDD missing loader spec")
        loader_path = String(get(loader_spec, "path", ""))
        grid = build_duo_grid(
            loader = (path = loader_path, reader = String(get(loader_spec, "reader", "auto"))),
        )
        lon_c = grid.cc_lon
        lat_c = grid.cc_lat
        R_val = Float64(grid.R)
    else
        error("unstructured_ode Layer-B runner: unsupported grid family '$family'")
    end

    Nc = length(lon_c)
    extra_ics = Dict{String,Float64}("u[$c]" => mms.ic(lon_c[c], lat_c[c]) for c in 1:Nc)

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path, extra_ics = extra_ics)

    du = similar(prob.u0)
    prob.f(du, prob.u0, prob.p, 0.0)

    inv_R2 = 1.0 / R_val^2
    err = 0.0
    for c in 1:Nc
        idx = get(var_map, "u[$c]", nothing)
        idx === nothing && error("var_map missing key 'u[$c]'")
        analytic = mms.derivative(lon_c[c], lat_c[c]) * inv_R2
        err = max(err, abs(du[idx] - analytic))
    end
    return err
end

"""
    _run_layer_b_cubed_sphere_cross_metric(rule, mms, n) -> Float64

Layer-B runner for `grid_family=cubed_sphere` cross_metric composite rules
(esd-ecq). Directly evaluates the rule on a CubedSphereGrid of resolution `n`
without routing through ESS `discretize → build_evaluator` (that pipeline
extension for cubed_sphere selectors is tracked at ESS/ess-cubed-sphere-pipeline).

Algorithm per interior cell (p, ic, jc) — physical cell index, not output index:
1. For each `term` in `main_spec["terms"]`:
   a. Load the named sub-stencil from `rule.discretizations[term["axis_stencil"]]`.
   b. Evaluate each stencil entry's `coeff` via `eval_coeff` with bindings
      `h = grid.dξ`, `J = grid.J[p,ic,jc]` (Jacobian at the center cell).
   c. Multiply each coeff by the ghost-extended field value at the offset cell.
   d. Sum across entries to get the stencil value.
   e. Weight by `metric_eval(grid, term["metric_component"], p, ic, jc)` and `term["sign"]`.
2. Sum all terms → discrete Laplacian at (p, ic, jc).
3. Compare against `mms.derivative(ξ[ic], η[jc])` for the L_inf error.
   The analytic Laplacian is panel-independent (gnomonic metric depends only on ξ,η).
"""
function _run_layer_b_cubed_sphere_cross_metric(rule::RuleFile, mms, n::Int)
    rule_doc = JSON.parse(read(rule.path, String))
    discretizations = get(rule_doc, "discretizations", Dict{String, Any}())
    main_spec = get(discretizations, rule.name, nothing)
    main_spec isa AbstractDict || error("rule spec missing for $(rule.name)")
    get(main_spec, "kind", "") == "cross_metric" ||
        error("cubed_sphere_cross_metric runner requires kind=\"cross_metric\"; got " * repr(get(main_spec, "kind", missing)))

    Nc = n
    grid = CubedSphereGrid(Nc; R = 1.0)
    h = Float64(grid.dξ)

    # Build test field on all (6, Nc, Nc) cells and extend with ghost cells.
    phi = zeros(Float64, 6, Nc, Nc)
    for p in 1:6, i in 1:Nc, j in 1:Nc
        phi[p, i, j] = mms.ic(grid.ξ_centers[i], grid.η_centers[j])
    end
    phi_ext = extend_with_ghosts(phi, grid)  # (6, Nc+2*Ng, Nc+2*Ng)
    Ng = grid.Ng

    result = zeros(Float64, 6, Nc - 2, Nc - 2)

    terms = get(main_spec, "terms", Any[])
    for term in terms
        axis_stencil = String(term["axis_stencil"])
        metric_name  = String(term["metric_component"])
        sign_val     = Float64(term["sign"])

        stencil_spec = get(discretizations, axis_stencil, nothing)
        stencil_spec isa AbstractDict ||
            error("sub-stencil $(axis_stencil) not found in $(rule.name) rule doc")
        entries = stencil_spec["stencil"]

        for p in 1:6, ii in 1:(Nc - 2), jj in 1:(Nc - 2)
            ic = ii + 1; jc = jj + 1
            mc = metric_eval(grid, metric_name, p, ic, jc)
            J_val = Float64(grid.J[p, ic, jc])
            bindings = Dict{String, Float64}("h" => h, "J" => J_val)

            stencil_sum = 0.0
            for entry in entries
                selectors = entry["selectors"]
                coeff = eval_coeff(entry["coeff"], bindings)
                xi_off = 0; eta_off = 0
                for sel in selectors
                    sel isa AbstractDict || continue
                    ax = String(get(sel, "axis", ""))
                    off = Int(get(sel, "offset", 0))
                    ax == "xi"  && (xi_off  = off)
                    ax == "eta" && (eta_off = off)
                end
                phi_val = phi_ext[p, ic + Ng + xi_off, jc + Ng + eta_off]
                stencil_sum += coeff * phi_val
            end

            result[p, ii, jj] += sign_val * mc * stencil_sum
        end
    end

    err = 0.0
    for p in 1:6, ii in 1:(Nc - 2), jj in 1:(Nc - 2)
        ic = ii + 1; jc = jj + 1
        ref = mms.derivative(grid.ξ_centers[ic], grid.η_centers[jc])
        err = max(err, abs(result[p, ii, jj] - Float64(ref)))
    end
    return err
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
    # 1. Load the rule spec (latlon rules carry authored replacement; esd-t4h).
    rule_doc = JSON.parse(read(rule.path, String))
    spec = get(get(rule_doc, "discretizations", Dict()), rule.name, nothing)
    spec isa AbstractDict || error("rule spec missing for $(rule.name)")
    repl = spec["replacement"]
    expr = (repl isa AbstractDict && get(repl, "op", nothing) == "arrayop") ?
           repl["expr"] : repl

    # 2. Build an n_lat × n_lon grid on the unit sphere.
    #    n_lat = n_lon = n (square angular resolution) — Y_{2,0} is lon-independent
    #    so convergence is driven by the lat axis; uniform n suffices.
    n_lat = n
    n_lon = n
    dlat = Float64(π) / n_lat
    dlon = Float64(2π) / n_lon
    lat_coord(i) = (i - 0.5) * dlat - π / 2
    lon_coord(j) = (j - 0.5) * dlon - π

    # 3. Build per-cell scalar ESM model. cos_lat varies per row, so each row
    #    gets its own parameter cos_lat_i. R, dlat, dlon are global parameters.
    variables = Dict{String, Any}()
    variables["R"]    = Dict{String, Any}("type" => "parameter", "default" => 1.0,  "units" => "1")
    variables["dlat"] = Dict{String, Any}("type" => "parameter", "default" => dlat, "units" => "1")
    variables["dlon"] = Dict{String, Any}("type" => "parameter", "default" => dlon, "units" => "1")
    for i in 1:n_lat
        variables["cos_lat_$(i)"] = Dict{String, Any}(
            "type" => "parameter", "default" => cos(lat_coord(i)), "units" => "1",
        )
        for j in 1:n_lon
            variables["u_$(i)_$(j)"] = Dict{String, Any}(
                "type" => "state", "default" => mms.ic(lat_coord(i), lon_coord(j)), "units" => "1",
            )
        end
    end

    # 4. Build equations: substitute cos_lat → cos_lat_i for each row, then
    #    resolve index($u, lat_arg, lon_arg) to u_i'_j' with periodic wrap.
    equations = Any[]
    for i in 1:n_lat
        row_expr = _layer_b_sub_str(expr, "cos_lat", "cos_lat_$(i)")
        for j in 1:n_lon
            rhs = _layer_b_2d_latlon_build_cell_expr(row_expr, i, j, n_lat, n_lon)
            push!(equations, Dict{String, Any}(
                "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u_$(i)_$(j)"], "wrt" => "t"),
                "rhs" => rhs,
            ))
        end
    end

    esm = Dict{String, Any}(
        "esm"      => "0.2.0",
        "metadata" => Dict{String, Any}("name" => "layer_b_latlon_n$(n)"),
        "models"   => Dict{String, Any}(
            "M" => Dict{String, Any}("variables" => variables, "equations" => equations),
        ),
    )

    # 5-6. Canonical pipeline.
    disc = EarthSciSerialization.discretize(esm)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

    # 7. Set MMS initial condition.
    for i in 1:n_lat, j in 1:n_lon
        u0[var_map["u_$(i)_$(j)"]] = mms.ic(lat_coord(i), lon_coord(j))
    end

    # 8. Evaluate discrete operator.
    du = similar(u0)
    f!(du, u0, p, 0.0)

    # 9. L_inf error over interior lat rows (exclude pole rows 1 and n_lat).
    #    The lat axis carries the convergence signal for Y_{2,0}; the lon terms
    #    vanish (Y_{2,0} is lon-independent). mms.derivative returns (d/dlat, d/dlon).
    err = 0.0
    for i in 2:(n_lat - 1), j in 1:n_lon
        e = abs(du[var_map["u_$(i)_$(j)"]] - mms.derivative(lat_coord(i), lon_coord(j))[1])
        err = max(err, e)
    end
    return err
end

# Recursive AST string-variable substitution: replace every string node equal to
# `old` with `new_name`. Used to bind per-row parameters such as cos_lat → cos_lat_i.
function _layer_b_sub_str(expr, old::String, new_name::String)
    if expr isa AbstractString
        return String(expr) == old ? new_name : String(expr)
    end
    expr isa Number && return expr
    expr isa AbstractDict || return expr
    return Dict{String, Any}(
        "op"   => String(expr["op"]),
        "args" => Any[_layer_b_sub_str(a, old, new_name) for a in expr["args"]],
    )
end

# Build a per-cell scalar ESM expression for the 2D lat-lon rule at cell (i, j).
# Resolves `index($u, lat_arg, lon_arg)` ops to variable names "u_i'_j'" with
# periodic wrap on both axes. All other nodes pass through unchanged.
function _layer_b_2d_latlon_build_cell_expr(expr, i::Int, j::Int, n_lat::Int, n_lon::Int)
    expr isa Number && return expr
    expr isa AbstractString && return String(expr)
    expr isa AbstractDict || error("unexpected node in latlon expr: $(typeof(expr))")
    op   = String(expr["op"])
    args = expr["args"]
    if op == "index"
        # Lowered latlon: args = [variable, lat_arg, lon_arg] (axes sorted "lat" < "lon").
        lat_idx = _layer_b_2d_eval_axis(args[2], "lat", i)
        lon_idx = _layer_b_2d_eval_axis(args[3], "lon", j)
        i_wrap  = mod(lat_idx - 1, n_lat) + 1
        j_wrap  = mod(lon_idx - 1, n_lon) + 1
        return "u_$(i_wrap)_$(j_wrap)"
    end
    return Dict{String, Any}(
        "op"   => op,
        "args" => Any[_layer_b_2d_latlon_build_cell_expr(a, i, j, n_lat, n_lon) for a in args],
    )
end

# Evaluate a 2D axis index sub-expression for axis `axis_name` at current index `cur`.
function _layer_b_2d_eval_axis(expr, axis_name::String, cur::Int)::Int
    if expr isa AbstractString
        s = String(expr)
        s == axis_name && return cur
        error("unexpected axis name '$(s)' in latlon index (expected '$(axis_name)')")
    end
    (expr isa Integer || expr isa AbstractFloat) && return Int(expr)
    op = String(expr["op"])
    a  = _layer_b_2d_eval_axis(expr["args"][1], axis_name, cur)
    b  = _layer_b_2d_eval_axis(expr["args"][2], axis_name, cur)
    op == "+" && return a + b
    op == "-" && return a - b
    error("unsupported axis op '$(op)' in latlon index")
end

"""
    _run_layer_b_2d_arakawa_periodic(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=arakawa` rules (esd-eg5).
Drives one resolution of the 2D periodic divergence convergence sweep through
the canonical `build_ode_problem` pipeline:

1. Locate the rule's integration ESM fixture (`fixtures/integration/<name>_pde.esm`)
   and the pre-built GDD (`discretizations/gdd/arakawa_2d_periodic_n<n>.gdd.json`).
2. Call `build_ode_problem(esm_path; grid_ref=gdd_path)` to construct an
   `ODEProblem` and `var_map` via the standard ESD injection path.
3. Override `u0` with the MMS vector field sampled at face positions:
   `Fu[i,j]` at `x=(i-1)*dx, y=(j-0.5)*dy` (face_x western face of cell i);
   `Fv[i,j]` at `x=(i-0.5)*dx, y=(j-1)*dy` (face_y southern face of cell j).
4. Evaluate `prob.f(du, u0, p, 0)` once to get discrete divergences.
5. Return the L_inf error against `mms.derivative(x, y)` at cell centers.

The expected order is O(h²) — `expected_min_order=1.9` in the convergence
fixture leaves a 0.1 margin over the n=16→32→64→128 sweep.
"""
function _run_layer_b_2d_arakawa_periodic(rule::RuleFile, mms, n::Int)
    disc_dir  = dirname(dirname(rule.path))  # .../discretizations
    esm_path  = joinpath(dirname(rule.path), rule.name, "fixtures", "integration",
                         "$(rule.name)_pde.esm")
    gdd_path  = joinpath(disc_dir, "gdd", "arakawa_2d_periodic_n$(n).gdd.json")

    isfile(esm_path) || error("arakawa runner: integration ESM not found: $esm_path")
    isfile(gdd_path) || error("arakawa runner: GDD not found: $gdd_path")

    prob, var_map = build_ode_problem(esm_path; grid_ref = gdd_path)

    dx = 1.0 / n
    dy = 1.0 / n
    u0 = copy(prob.u0)
    for j in 1:n, i in 1:n
        u0[var_map["Fu[$i,$j]"]] = mms.ic((i - 1) * dx, (j - 0.5) * dy)[1]
        u0[var_map["Fv[$i,$j]"]] = mms.ic((i - 0.5) * dx, (j - 1) * dy)[2]
    end

    du = similar(u0)
    prob.f(du, u0, prob.p, 0.0)

    err = 0.0
    for j in 1:n, i in 1:n
        cx = (i - 0.5) * dx
        cy = (j - 0.5) * dy
        e  = abs(du[var_map["div[$i,$j]"]] - mms.derivative(cx, cy))
        err = max(err, e)
    end
    return err
end

"""
    _run_layer_b_1d_cartesian_periodic(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=cartesian` rules with a
single periodic dimension (esd-0ip; ArrayOp-native under dsc-kswm). Drives
one resolution of the MMS convergence sweep through the canonical pipeline
**without per-cell scalarization**: the ESM document describes the PDE on
the array variable, and ESS does the rest.

1. Load the rule spec from the catalog JSON.
2. Build a PDE ESM document: a 1D periodic cartesian `grids` entry of size
   `n`, an array state `u` with `shape: ["i"]`, parameter `dx = 1/n`, and
   the single equation `D(u, wrt=t) = <applies_to instantiated on u>`
   (e.g. `grad(u, dim="i")`). The rule rides one of two document forms:
   - stencil form → `lower_stencil_to_scheme` emits the
     `discretizations.<name>` scheme + `use:` rule (ESS RFC §7.2.1); or
   - replacement form → the `applies_to` pattern + the authored
     replacement AST as a plain document rule. ESS pattern-matches the
     `dim` bindings in `applies_to` and expands the arrayop accordingly.
3. `EarthSciSerialization.discretize(esm; lift_1d_arrayop=true)` rewrites
   the PDE op via the rule engine / scheme expansion, lifts the equation
   to arrayop form with ranges from the grid, and folds periodic boundary
   accesses.
4. `EarthSciSerialization.build_evaluator` expands the arrayop natively
   (tree-walk) into the ODE RHS `f!`.
5. Set each `u[i] = mms.ic((i-0.5)/n)`, evaluate `f!(du, u0, p, 0)` once,
   and return the L_inf error against `mms.derivative((i-0.5)/n)`.
"""
function _run_layer_b_1d_cartesian_periodic(rule::RuleFile, mms, n::Int)
    # 1. Load the rule spec.
    rule_doc = JSON.parse(read(rule.path, String))
    spec = get(get(rule_doc, "discretizations", Dict()), rule.name, nothing)
    spec isa AbstractDict || error("rule spec missing for $(rule.name)")

    # 2. Build the PDE ESM document on the unit interval [0, 1].
    dx = 1.0 / n
    cell_x(i) = (i - 0.5) * dx

    # Pattern bindings: auxiliary pattern variables (declared in the MMS
    # kind's aux table, e.g. the velocity $U of div($U * $q)) bind to their
    # own shaped fields; everything else binds to the state u.
    aux = hasproperty(mms, :aux) ? mms.aux : Dict{String, Any}()
    bindings = _layer_b_pattern_bindings(spec["applies_to"], aux)
    aux_fields = sort!(unique(String[v for v in values(bindings) if v != "u"]))

    variables = Dict{String, Any}(
        "u" => Dict{String, Any}(
            "type" => "state", "default" => 0.0, "units" => "1",
            "shape" => Any["x"], "location" => "cell_center",
        ),
        "dx" => Dict{String, Any}(
            "type" => "parameter", "default" => dx, "units" => "1",
        ),
    )
    equations = Any[
        Dict{String, Any}(
            "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u"], "wrt" => "t"),
            "rhs" => _layer_b_instantiate_applies_to(spec, bindings, "x"),
        ),
    ]
    for name in aux_fields
        variables[name] = Dict{String, Any}(
            "type" => "state", "default" => 0.0, "units" => "1",
            "shape" => Any["x"], "location" => "cell_center",
        )
        # Frozen auxiliary field: zero tendency, cells set from the aux value.
        push!(equations, Dict{String, Any}(
            "lhs" => Dict{String, Any}("op" => "D", "args" => Any[name], "wrt" => "t"),
            "rhs" => 0.0,
        ))
    end

    esm = Dict{String, Any}(
        "esm"      => "0.4.0",
        "metadata" => Dict{String, Any}("name" => "layer_b_1d_cartesian_periodic_n$(n)"),
        "grids"    => Dict{String, Any}(
            "g" => Dict{String, Any}(
                "family"     => "cartesian",
                "dimensions" => Any[
                    Dict{String, Any}("name" => "x", "size" => n,
                                      "periodic" => true, "spacing" => "uniform"),
                ],
            ),
        ),
        "models" => Dict{String, Any}(
            "M" => Dict{String, Any}(
                "grid"      => "g",
                "variables" => variables,
                "equations" => equations,
            ),
        ),
    )

    if haskey(spec, "stencil")
        # Stencil form → ESS scheme + use: rule (canonical §7.2.1 path).
        scheme, use_rule = lower_stencil_to_scheme(rule.name, spec)
        esm["discretizations"] = Dict{String, Any}(rule.name => scheme)
        esm["rules"] = Any[use_rule]
    else
        # Replacement form → plain pattern + authored replacement AST.
        esm["rules"] = Any[Dict{String, Any}(
            "name"        => rule.name,
            "pattern"     => spec["applies_to"],
            "replacement" => spec["replacement"],
        )]
    end

    # 3-4. Canonical pipeline: discretize (rule rewrite + arrayop lift +
    # periodic folding) → tree-walk build_evaluator.
    disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

    # 5. MMS initial conditions (state + frozen auxiliary fields), single
    # RHS evaluation, L_inf error over the state cells.
    for i in 1:n
        u0[var_map["u[$(i)]"]] = mms.ic(cell_x(i))
    end
    for name in aux_fields
        val = aux[name]
        for i in 1:n
            u0[var_map["$(name)[$(i)]"]] = val isa Function ? Float64(val(cell_x(i))) : Float64(val)
        end
    end
    du = similar(u0)
    f!(du, u0, p, 0.0)
    return maximum(abs(du[var_map["u[$(i)]"]] - mms.derivative(cell_x(i))) for i in 1:n)
end

"""
    _run_layer_b_2d_cartesian_periodic(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=cartesian` rules whose
stencil spans two axes (dsc-vst2; ArrayOp-native under dsc-kswm). Same
document-driven design as `_run_layer_b_1d_cartesian_periodic`, on an
n×n periodic grid:

1. Build a PDE ESM document: 2D periodic cartesian grid (`x`, `y`, both
   size `n`), array state `u` with `shape: ["x", "y"]`, parameters
   `dx = dy = 1/n`, and the equation `D(u, wrt=t) = <applies_to on u>`
   (e.g. `laplacian(u)`).
2. Supply the rule as a pattern + replacement document rule using the
   authored canonical replacement (bare `+` expr with `i`/`j` RFC §7.1
   reserved keywords). The dim-less pattern like `laplacian(\$u)` binds
   no axis variables; canonical component names are resolved positionally
   by the arrayop lift regardless of the document's dimension names.
3. `discretize` rewrites the PDE op, lifts to arrayop (multidimensional
   variables lift by default — no `lift_1d_arrayop` needed), and folds
   periodic boundary accesses; `build_evaluator` expands the arrayop
   natively.
4. Set `u[i,j] = mms.ic(x_i, y_j)`, evaluate `f!` once, and return the
   L_inf error against `mms.derivative(x_i, y_j)`.
"""
function _run_layer_b_2d_cartesian_periodic(rule::RuleFile, mms, n::Int)
    rule_doc = JSON.parse(read(rule.path, String))
    spec = get(get(rule_doc, "discretizations", Dict()), rule.name, nothing)
    spec isa AbstractDict || error("rule spec missing for $(rule.name)")

    expr = spec["replacement"]

    h = 1.0 / n
    cell(c) = (c - 0.5) * h

    esm = Dict{String, Any}(
        "esm"      => "0.4.0",
        "metadata" => Dict{String, Any}("name" => "layer_b_2d_cartesian_periodic_n$(n)"),
        "grids"    => Dict{String, Any}(
            "g" => Dict{String, Any}(
                "family"     => "cartesian",
                "dimensions" => Any[
                    Dict{String, Any}("name" => "x", "size" => n,
                                      "periodic" => true, "spacing" => "uniform"),
                    Dict{String, Any}("name" => "y", "size" => n,
                                      "periodic" => true, "spacing" => "uniform"),
                ],
            ),
        ),
        "rules" => Any[Dict{String, Any}(
            "name"        => rule.name,
            "pattern"     => spec["applies_to"],
            "replacement" => expr,
        )],
        "models" => Dict{String, Any}(
            "M" => Dict{String, Any}(
                "grid" => "g",
                "variables" => Dict{String, Any}(
                    "u" => Dict{String, Any}(
                        "type" => "state", "default" => 0.0, "units" => "1",
                        "shape" => Any["x", "y"], "location" => "cell_center",
                    ),
                    "dx" => Dict{String, Any}(
                        "type" => "parameter", "default" => h, "units" => "1",
                    ),
                    "dy" => Dict{String, Any}(
                        "type" => "parameter", "default" => h, "units" => "1",
                    ),
                ),
                "equations" => Any[
                    Dict{String, Any}(
                        "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u"], "wrt" => "t"),
                        "rhs" => _layer_b_instantiate_applies_to(spec, "u", ["x", "y"]),
                    ),
                ],
            ),
        ),
    )

    disc = EarthSciSerialization.discretize(esm)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

    for j in 1:n, i in 1:n
        u0[var_map["u[$(i),$(j)]"]] = mms.ic(cell(i), cell(j))
    end
    du = similar(u0)
    f!(du, u0, p, 0.0)
    return maximum(
        abs(du[var_map["u[$(i),$(j)]"]] - mms.derivative(cell(i), cell(j)))
        for j in 1:n, i in 1:n)
end

"""
    _run_layer_b_fv_cell_average_1d(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for FV reconstruction rules whose
convergence fixture declares `sampling: "cell_average"` (dsc-a7b2;
ArrayOp-native under dsc-kswm). Reconstruction rules estimate point
values from cell AVERAGES, so this runner differs from the point-sample
families in two ways:

- The IC is the **exact cell average** of the MMS field:
  `u[i] = mms.ic(x_{i-1}, x_i)` over the cell bounds (no quadrature
  error pollutes the convergence order).
- The analytic target for each output is the **point value**
  `mms.derivative(x)` at the location the output's stencil support
  implies: the midpoint of the offset range. CW84's `q_left_edge`
  spans offsets `{-2..1}` → midpoint −1/2 → the cell's left edge
  `x_{i-1/2}`; `q_right_edge` spans `{-1..2}` → +1/2 → `x_{i+1/2}`.
  This is a generic schema attribute (stencil support), not per-rule
  dispatch.

Multi-output stencil objects (the catalog's extension ahead of the RFC —
`stencil` keyed by output name) lower **one output per ESM document**
via `lower_stencil_to_scheme(...; output=...)`; each document drives the
ordinary scheme + `use:` pipeline and the runner returns the worst
L_inf error across outputs. Document-level multi-output emission (one
rewrite producing all named outputs at once, needed for the Layer-A
canonical byte contract) awaits an RFC §7 extension — this runner
covers the numerical-convergence acceptance only.

Outputs documented without stencil rows (e.g. CW84's `parabola.{a_L,
a_R, da, a_6}` block) cannot be driven through the pipeline and are not
exercised here.
"""
function _run_layer_b_fv_cell_average_1d(rule::RuleFile, mms, n::Int)
    rule_doc = JSON.parse(read(rule.path, String))
    spec = get(get(rule_doc, "discretizations", Dict()), rule.name, nothing)
    spec isa AbstractDict || error("rule spec missing for $(rule.name)")
    stencil_field = get(spec, "stencil", nothing)
    stencil_field !== nothing || error(
        "fv_cell_average_1d runner requires a stencil-form rule; $(rule.name) has none")

    outputs = stencil_field isa AbstractDict ?
        Any[String(k) for k in sort!(collect(String.(keys(stencil_field))))] :
        Any[nothing]

    h = 1.0 / n
    worst = 0.0
    for output in outputs
        scheme, use_rule = lower_stencil_to_scheme(rule.name, spec; output = output)

        esm = Dict{String, Any}(
            "esm"      => "0.4.0",
            "metadata" => Dict{String, Any}(
                "name" => "layer_b_fv_cell_average_1d_n$(n)" *
                          (output === nothing ? "" : "_$(output)")),
            "grids"    => Dict{String, Any}(
                "g" => Dict{String, Any}(
                    "family"     => "cartesian",
                    "dimensions" => Any[
                        Dict{String, Any}("name" => "x", "size" => n,
                                          "periodic" => true, "spacing" => "uniform"),
                    ],
                ),
            ),
            "discretizations" => Dict{String, Any}(rule.name => scheme),
            "rules"           => Any[use_rule],
            "models" => Dict{String, Any}(
                "M" => Dict{String, Any}(
                    "grid" => "g",
                    "variables" => Dict{String, Any}(
                        "u" => Dict{String, Any}(
                            "type" => "state", "default" => 0.0, "units" => "1",
                            "shape" => Any["x"], "location" => "cell_center",
                        ),
                        "dx" => Dict{String, Any}(
                            "type" => "parameter", "default" => h, "units" => "1",
                        ),
                    ),
                    "equations" => Any[
                        Dict{String, Any}(
                            "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u"], "wrt" => "t"),
                            "rhs" => _layer_b_instantiate_applies_to(spec, "u", "x"),
                        ),
                    ],
                ),
            ),
        )

        disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)
        f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

        for i in 1:n
            u0[var_map["u[$(i)]"]] = mms.ic((i - 1) * h, i * h)
        end
        du = similar(u0)
        f!(du, u0, p, 0.0)

        # Output location from the stencil support midpoint (in cell-width
        # units relative to the cell center).
        entries = output === nothing ? stencil_field : stencil_field[output]
        offsets = Int[Int(e["selector"]["offset"]) for e in entries]
        mid = (minimum(offsets) + maximum(offsets)) / 2
        err = maximum(
            abs(du[var_map["u[$(i)]"]] - mms.derivative((i - 0.5 + mid) * h))
            for i in 1:n)
        worst = max(worst, err)
    end
    return worst
end

# Map every pattern variable in `applies_to.args` to a document variable
# name: variables whose STRIPPED name appears in the MMS kind's auxiliary
# table keep that name (auxiliary fields, e.g. `$U → U`); every other
# pattern variable binds to the transported state `"u"`. Purely
# declarative — the aux table, not rule identity, decides the split.
function _layer_b_pattern_bindings(applies_to::AbstractDict, aux::AbstractDict)
    bindings = Dict{String, String}()
    walk(node) = if node isa AbstractString && startswith(String(node), "\$")
        s = String(node)
        bindings[s] = haskey(aux, s[2:end]) ? s[2:end] : "u"
    elseif node isa AbstractDict
        foreach(walk, get(node, "args", Any[]))
    end
    foreach(walk, get(applies_to, "args", Any[]))
    return bindings
end

# Instantiate a rule's `applies_to` pattern as a concrete equation RHS:
# each pattern variable in the args tree (including nested operands like
# `div($U * $q)`) is substituted per `bindings`, and every `dim` pattern
# variable — including those on NESTED operators, e.g.
# grad(grad($u, dim=$y), dim=$x) — maps to a document dimension name.
# Dim pattern variables map positionally by sorted name onto `dim_names`
# (sorted $x → dim 1, $y → dim 2), so the equation's dims, the document's
# grid, and the authored replacement's canonical components stay aligned.
# Purely structural — mirrors how ESS pattern matching re-binds the same
# variables during `discretize`.
function _layer_b_instantiate_applies_to(spec::AbstractDict, bindings::Dict{String, String},
                                         dim_names::Vector{String})
    applies_to = get(spec, "applies_to", nothing)
    applies_to isa AbstractDict || error("rule spec has no applies_to object")

    dim_vars = String[]
    _collect_dims(node) = if node isa AbstractDict
        d = get(node, "dim", nothing)
        if d isa AbstractString && startswith(String(d), "\$") && !(String(d) in dim_vars)
            push!(dim_vars, String(d))
        end
        foreach(_collect_dims, get(node, "args", Any[]))
    end
    _collect_dims(applies_to)
    sort!(dim_vars)
    length(dim_vars) <= length(dim_names) || error(
        "applies_to binds $(length(dim_vars)) dim pattern variables but the document " *
        "declares only $(length(dim_names)) dimensions")
    dim_map = Dict{String, String}(v => dim_names[k] for (k, v) in enumerate(dim_vars))

    subst(node) =
        node isa AbstractString ? get(bindings, String(node), String(node)) :
        node isa AbstractDict ? begin
            out = Dict{String, Any}(
                String(k) => (String(k) == "args" ? Any[subst(a) for a in v] : v)
                for (k, v) in node)
            d = get(out, "dim", nothing)
            if d isa AbstractString && startswith(String(d), "\$")
                out["dim"] = dim_map[String(d)]
            end
            out
        end :
        node
    return subst(applies_to)
end

# Single-dimension convenience used by the 1D runners.
function _layer_b_instantiate_applies_to(spec::AbstractDict, bindings::Dict{String, String},
                                         dim_name::String)
    return _layer_b_instantiate_applies_to(spec, bindings, [dim_name])
end

# Convenience for runners without auxiliary fields: bind every pattern
# variable to the single operand.
function _layer_b_instantiate_applies_to(spec::AbstractDict, operand::String,
                                         dim_names::Union{String, Vector{String}})
    applies_to = get(spec, "applies_to", nothing)
    applies_to isa AbstractDict || error("rule spec has no applies_to object")
    bindings = _layer_b_pattern_bindings(applies_to, Dict{String, Any}())
    for k in keys(bindings)
        bindings[k] = operand
    end
    dims = dim_names isa String ? [dim_names] : dim_names
    return _layer_b_instantiate_applies_to(spec, bindings, dims)
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
# Any `$`-prefixed string (the rule's axis pattern variable, e.g. `"$x"` for
# cartesian or `"$k"` for vertical) is treated as the current cell index `i`.
function _layer_b_1d_eval_index(expr, i::Int)::Int
    if expr isa AbstractString
        s = String(expr)
        startswith(s, "\$") && return i
    end
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
    _run_layer_b_1d_vertical_column(rule, mms, n) -> Float64

Layer-B canonical-pipeline runner for `grid_family=vertical` rules (esd-bbp).
Still drives the per-cell-scalar machinery (per-cell ESM → discretize →
build_evaluator) with spacing parameter `"h"` instead of `"dx"`; an
ArrayOp-native port like `_run_layer_b_1d_cartesian_periodic`'s awaits ESS
scheme support for the vertical grid family (the scheme `grid_family` enum
is cartesian/cubed_sphere/unstructured today).

Loads the replacement AST from the rule's canonical Layer-A fixture.  Vertical
rules now carry authored `replacement` ASTs (esd-t4h), but the canonical fixture
is the authoritative cell-center form: `(u[\$x+1] - u[\$x-1]) / (2*h)`.
"""
function _run_layer_b_1d_vertical_column(rule::RuleFile, mms, n::Int)
    # Load the cell-center centered-difference replacement from the canonical fixture.
    canonical_path = joinpath(dirname(rule.path), rule.name, "fixtures", "canonical", "input.esm")
    if !isfile(canonical_path)
        error("_run_layer_b_1d_vertical_column: no canonical fixture at $(relpath_from_repo(canonical_path))")
    end
    canonical_doc = JSON.parse(read(canonical_path, String))
    rules_list = get(canonical_doc, "rules", Any[])
    isempty(rules_list) && error("canonical fixture for $(rule.name) has no rules entries")
    canonical_rule = first(rules_list)
    repl = get(canonical_rule, "replacement", nothing)
    repl isa AbstractDict || error("canonical fixture for $(rule.name) has no replacement object")
    expr = (get(repl, "op", nothing) == "arrayop") ? repl["expr"] : repl

    h = 1.0 / n
    cell_z(k) = (k - 0.5) * h

    variables = Dict{String, Any}()
    for k in 1:n
        variables["u_$(k)"] = Dict{String, Any}(
            "type"    => "state",
            "default" => mms.ic(cell_z(k)),
            "units"   => "1",
        )
    end
    variables["h"] = Dict{String, Any}(
        "type"    => "parameter",
        "default" => h,
        "units"   => "1",
    )

    equations = Any[]
    for k in 1:n
        rhs = _layer_b_1d_build_cell_expr(expr, k, n)
        push!(equations, Dict{String, Any}(
            "lhs" => Dict{String, Any}("op" => "D", "args" => Any["u_$(k)"], "wrt" => "t"),
            "rhs" => rhs,
        ))
    end

    esm = Dict{String, Any}(
        "esm"      => "0.2.0",
        "metadata" => Dict{String, Any}("name" => "layer_b_1d_vertical_n$(n)"),
        "models"   => Dict{String, Any}(
            "M" => Dict{String, Any}(
                "variables" => variables,
                "equations" => equations,
            ),
        ),
    )

    disc = EarthSciSerialization.discretize(esm)
    f!, u0, p, _tspan, var_map = EarthSciSerialization.build_evaluator(disc)

    for k in 1:n
        u0[var_map["u_$(k)"]] = mms.ic(cell_z(k))
    end

    du = similar(u0)
    f!(du, u0, p, 0.0)

    return maximum(abs(du[var_map["u_$(k)"]] - mms.derivative(cell_z(k))) for k in 1:n)
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

# Resolve a canonical index argument (bare "i"/"j" or {"op":"+","args":["i"/"j",off]}).
function _resolve_repl_axis(arg, ci::Int, cj::Int)::Int
    if arg isa AbstractString
        s = String(arg)
        s == "i" && return ci
        s == "j" && return cj
        error("unknown bare axis variable '$s' in canonical replacement")
    end
    arg isa AbstractDict || error("unexpected index arg type $(typeof(arg))")
    op = String(arg["op"])
    base_name = String(arg["args"][1])
    off = Int(arg["args"][2])
    base = (base_name == "i") ? ci : (base_name == "j" ? cj : error("unknown axis var '$base_name'"))
    op == "+" && return base + off
    op == "-" && return base - off
    error("unsupported index offset op '$op'")
end

# Evaluate a canonical bare replacement AST at cell (ci, cj) with periodic wrap.
function _eval_repl_2d(node, fields::Dict{String, Matrix{Float64}}, bindings::Dict{String, Float64}, ci::Int, cj::Int, nx::Int, ny::Int)::Float64
    node isa Number && return Float64(node)
    if node isa AbstractString
        s = String(node)
        s == "i" && return Float64(ci)
        s == "j" && return Float64(cj)
        haskey(bindings, s) && return bindings[s]
        error("unresolved variable '$s' in canonical replacement")
    end
    node isa AbstractDict || error("unexpected node type $(typeof(node))")
    op   = String(node["op"])
    args = node["args"]
    if op == "index"
        length(args) == 3 || error("index expects 3 args, got $(length(args))")
        fname = String(args[1])
        mat   = get(fields, fname, nothing)
        mat === nothing && error("unknown field '$fname' in canonical replacement")
        ii = _resolve_repl_axis(args[2], ci, cj)
        jj = _resolve_repl_axis(args[3], ci, cj)
        return mat[mod(ii - 1, nx) + 1, mod(jj - 1, ny) + 1]
    end
    ev(a) = _eval_repl_2d(a, fields, bindings, ci, cj, nx, ny)
    op == "+" && return sum(ev(a) for a in args)
    op == "*" && return prod(ev(a) for a in args)
    op == "/" && return ev(args[1]) / ev(args[2])
    op == "-" && return length(args) == 1 ? -ev(args[1]) : ev(args[1]) - ev(args[2])
    error("unsupported op '$op' in canonical replacement")
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
    if !(spec isa AbstractDict)
        return LayerResult(
            LAYER_FAIL,
            "rule $(rule.name) has no spec under discretizations.$(rule.name)",
        )
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

    div_sum = 0.0
    if haskey(spec, "stencil")
        stencil = spec["stencil"]
        if !(stencil isa AbstractVector) || isempty(stencil)
            return LayerResult(LAYER_FAIL, "rule $(rule.name) has empty stencil")
        end
        # Evaluate stencil entries at every cell; periodic wrap on both axes.
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
    elseif haskey(spec, "replacement")
        # Canonical bare replacement (esd-eg5): evaluate AST at each cell.
        repl   = spec["replacement"]
        fields = Dict{String, Matrix{Float64}}("\$Fx" => Fx, "\$Fy" => Fy)
        for j in 1:ny, i in 1:nx
            div_sum += _eval_repl_2d(repl, fields, bindings, i, j, nx, ny)
        end
    else
        return LayerResult(
            LAYER_FAIL,
            "rule $(rule.name) has no `stencil` or `replacement` under discretizations.$(rule.name)",
        )
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
