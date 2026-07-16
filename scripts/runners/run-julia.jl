#!/usr/bin/env julia
# run-julia.jl — the official Julia conformance runner for EarthSciDiscretizations.
#
# THIN WRAPPER (AGENTS.md §2, single pathway): every number and byte this script
# emits comes from official EarthSciAST.jl entry points driven over the
# one canonical pipeline (.esm → parse → §9.7 import/metaparameter resolution →
# §9.6.3 rewrite fixpoint → official runner). No evaluator, rule engine, or
# numeric kernel lives here — only argument marshalling and JSON I/O.
#
# CLI (EarthSciAST CONFORMANCE_SPEC.md §4.2):
#   julia scripts/runners/run-julia.jl --output-dir <path> \
#         [--categories ast,simulation,convergence,regridding,reprojection] \
#         [--files <manifest.json>[,<manifest.json>…]] [--verbose]
#
# Outputs (§4.4 layout):
#   <output-dir>/julia_<category>_results.json     one file per category run
#   <output-dir>/julia_summary.json                roll-up
#   <output-dir>/ast/<case>/expanded.golden.json   canonical post-lowering bytes
#                                                  (byte-compared vs the golden)
#
# Category → official ESS pathway:
#   ast          raw §9.7 pipeline (resolve_template_machinery →
#                lower_expression_templates), canonical bytes exactly as
#                $ESS_ROOT/scripts/generate-template-import-goldens.jl emits them
#   simulation   run_pde_tests (§6.6/§6.6.5 inline tests over simulate())
#   convergence  load(problem; metaparameters=…) per manifest resolution →
#                simulate → evaluate_cellwise(reference) → field_reduce norms
#   regridding   run_pde_tests on the fixture (exact-invariant gates) +
#                the regrid_state(1) field + the per-pair A_ij/A_j/W_ij setup
#                arrays via the BuildInspection observability surface
#   reprojection load a runner-built template invocation (manifest fixture:
#                null) → evaluate_expr at every golden point
#
# Environment: the dedicated ESS pde_sim_adapter project (pins the dev'd
# EarthSciAST + OrdinaryDiffEqTsit5 + JSON3) — the same env the ESS
# PDE-simulation conformance adapter runs in.

import Pkg

# ---------------------------------------------------------------------------
# Locate repos + bootstrap the official solver-bearing environment.
# ---------------------------------------------------------------------------
const ESD_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Walk up from the repo root looking for a sibling EarthSciAST checkout. The first
# candidate is ESD_ROOT/../EarthSciAST, the ordinary sibling layout. The walk matters
# inside a git worktree, where ESD_ROOT is <main>/.claude/worktrees/<name> and the
# sibling checkout is several levels up rather than one.
function _locate_ess(start::AbstractString)
    haskey(ENV, "ESS_ROOT") && return normpath(ENV["ESS_ROOT"])
    dir = start
    while true
        cand = normpath(joinpath(dir, "..", "EarthSciAST"))
        isfile(joinpath(cand, "esm-schema.json")) && return cand
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    return normpath(joinpath(start, "..", "EarthSciAST"))
end

const ESS_ROOT = _locate_ess(ESD_ROOT)
isfile(joinpath(ESS_ROOT, "esm-schema.json")) ||
    error("EarthSciAST not found at '$ESS_ROOT'; set ESS_ROOT " *
          "or clone it as a sibling checkout (scripts/ess-locate.sh contract)")

let env = joinpath(ESS_ROOT, "pkg", "EarthSciAST.jl", "scripts", "pde_sim_adapter")
    Pkg.activate(env; io=devnull)
    isfile(joinpath(env, "Manifest.toml")) ||
        Pkg.develop(path=normpath(joinpath(env, "..", "..")); io=devnull)
    Pkg.instantiate(; io=devnull)
end

using EarthSciAST
using EarthSciAST: resolve_template_machinery, lower_expression_templates,
    JSONLikeDict
using JSON3
import OrdinaryDiffEqTsit5
# The spherical geometry kernel (manifold "spherical" on intersect_polygon /
# polygon_intersection_area) lives in the EarthSciASTGeometryOpsExt
# package extension, which loads only when GeometryOps + GeoInterface are in
# the session (both are pinned deps of the pde_sim_adapter environment). The
# planar path never needs them; the spherical regridding conformance case does.
import GeometryOps, GeoInterface
const ODE = OrdinaryDiffEqTsit5
const ESS = EarthSciAST

const ALL_CATEGORIES = ["ast", "simulation", "convergence", "regridding", "reprojection"]
const CONF = joinpath(ESD_ROOT, "tests", "conformance")

# ---------------------------------------------------------------------------
# CLI parsing (§4.2).
# ---------------------------------------------------------------------------
function parse_cli(args)
    output_dir = nothing
    categories = String[]
    files = String[]
    verbose = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--output-dir"
            output_dir = args[i+1]; i += 2
        elseif a == "--categories"
            append!(categories, split(args[i+1], ","; keepempty=false)); i += 2
        elseif a == "--files"
            append!(files, split(args[i+1], ","; keepempty=false)); i += 2
        elseif a == "--verbose"
            verbose = true; i += 1
        else
            error("unknown argument '$a' (usage: --output-dir <path> " *
                  "[--categories <list>] [--files <list>] [--verbose])")
        end
    end
    output_dir === nothing && error("--output-dir is required (CONFORMANCE_SPEC.md §4.2)")
    return (String(output_dir), String.(categories), String.(files), verbose)
end

# ---------------------------------------------------------------------------
# Canonical golden writer — byte-for-byte the scheme the ESS Julia reference
# uses for §9.7 expansion goldens ($ESS_ROOT/scripts/
# generate-template-import-goldens.jl): object keys sorted, arrays in order,
# scalars via JSON3, 2-space indent, trailing newline.
# ---------------------------------------------------------------------------
_norm(x) =
    (x isa AbstractDict || x isa JSON3.Object) ?
        Dict{String,Any}(string(k) => _norm(v) for (k, v) in pairs(x)) :
    (x isa AbstractVector || x isa JSON3.Array) ? Any[_norm(v) for v in x] : x

function _write_sorted(io::IO, x, indent::Int)
    pad = "  "^indent
    pad1 = "  "^(indent + 1)
    if x isa AbstractDict
        isempty(x) && return print(io, "{}")
        print(io, "{\n")
        ks = sort(collect(keys(x)))
        for (i, k) in enumerate(ks)
            print(io, pad1, JSON3.write(string(k)), ": ")
            _write_sorted(io, x[k], indent + 1)
            i < length(ks) && print(io, ",")
            print(io, "\n")
        end
        print(io, pad, "}")
    elseif x isa AbstractVector
        isempty(x) && return print(io, "[]")
        print(io, "[\n")
        for (i, v) in enumerate(x)
            print(io, pad1)
            _write_sorted(io, v, indent + 1)
            i < length(x) && print(io, ",")
            print(io, "\n")
        end
        print(io, pad, "]")
    else
        print(io, JSON3.write(x))
    end
end

function canonical_bytes(doc)::String
    io = IOBuffer()
    _write_sorted(io, _norm(doc), 0)
    print(io, "\n")
    return String(take!(io))
end

# ---------------------------------------------------------------------------
# Manifest discovery.
# ---------------------------------------------------------------------------
function discover_manifests(category::String, files::Vector{String})
    if !isempty(files)
        out = Tuple{String,Any}[]
        for f in files
            p = isabspath(f) ? f : joinpath(ESD_ROOT, f)
            m = JSON3.read(read(p, String))
            String(m["category"]) == category && push!(out, (dirname(abspath(p)), m))
        end
        return out
    end
    dir = joinpath(CONF, category)
    isdir(dir) || return Tuple{String,Any}[]
    out = Tuple{String,Any}[]
    for case in sort(readdir(dir))
        mp = joinpath(dir, case, "manifest.json")
        isfile(mp) || continue
        push!(out, (joinpath(dir, case), JSON3.read(read(mp, String))))
    end
    return out
end

function integrator_opts(manifest)
    haskey(manifest, "integrators") && haskey(manifest["integrators"], "julia") ||
        return (ODE.Tsit5(), 1e-10, 1e-12)
    j = manifest["integrators"]["julia"]
    alg_name = String(get(j, "algorithm", "Tsit5"))
    alg_name == "Tsit5" ||
        error("manifest pins Julia algorithm '$alg_name'; this runner supports Tsit5")
    return (ODE.Tsit5(), Float64(get(j, "reltol", 1e-10)), Float64(get(j, "abstol", 1e-12)))
end

assertion_dicts(results) = [Dict{String,Any}(
        "model" => r.model, "test_id" => r.test_id, "assertion_idx" => r.assertion_idx,
        "variable" => r.variable, "time" => r.time, "reduce" => r.reduce,
        "expected" => r.expected, "actual" => r.actual,
        "rtol" => r.rtol, "atol" => r.atol, "passed" => r.passed,
        "message" => r.message) for r in results]

# ---------------------------------------------------------------------------
# ast — raw §9.7 expansion to canonical bytes + post-lowering gates.
# ---------------------------------------------------------------------------

# Post-lowering gates (the ast manifests' contract): zero
# apply_expression_template nodes, zero unlowered rewrite-target ops (spatial
# D / grad / div / laplacian), zero surviving metaparameter machinery.
function unlowered_violations(node, acc::Vector{String}=String[])
    if node isa AbstractDict
        op = string(get(node, "op", ""))
        op == "apply_expression_template" && push!(acc, "apply_expression_template")
        op in ("grad", "div", "laplacian") && push!(acc, "unlowered $op")
        (op == "D" && string(get(node, "wrt", "t")) != "t") &&
            push!(acc, "unlowered spatial D (wrt=$(node["wrt"]))")
        for (k, v) in pairs(node)
            k in (:metadata, Symbol("metadata")) && continue  # prose may cite op names
            unlowered_violations(v, acc)
        end
    elseif node isa AbstractVector
        for v in node
            unlowered_violations(v, acc)
        end
    end
    return acc
end

function run_ast(output_dir, files, verbose)
    cases = Dict{String,Any}()
    for (case_dir, manifest) in discover_manifests("ast", files)
        case = String(manifest["case"])
        rec = Dict{String,Any}("case" => case, "status" => "ok")
        try
            fixture = joinpath(case_dir, String(manifest["fixture"]))
            raw = JSON3.read(read(fixture, String))
            resolved = resolve_template_machinery(raw, dirname(fixture))
            out = lower_expression_templates(resolved === nothing ? raw : resolved)
            doc = out isa JSONLikeDict ? getfield(out, :data) : out
            doc_n = _norm(doc)
            # The prose-free structural gates.
            violations = String[]
            for (mname, model) in get(doc_n, "models", Dict{String,Any}())
                unlowered_violations(Dict{String,Any}(k => v for (k, v) in model
                                                      if k != "metadata"), violations)
            end
            isempty(violations) ||
                error("post-lowering gate failed: $(join(unique(violations), ", "))")
            bytes = canonical_bytes(doc_n)
            art_dir = joinpath(output_dir, "ast", case)
            mkpath(art_dir)
            art = joinpath(art_dir, "expanded.golden.json")
            write(art, bytes)
            rec["artifact"] = relpath(art, output_dir)
            rec["bytes"] = ncodeunits(bytes)
            rec["golden"] = relpath(joinpath(case_dir, String(manifest["golden"])), ESD_ROOT)
        catch err
            rec["status"] = "error"
            rec["message"] = sprint(showerror, err)
        end
        verbose && println("  ast/$case: $(rec["status"])")
        cases[case] = rec
    end
    return Dict{String,Any}("binding" => "julia", "category" => "ast", "cases" => cases)
end

# ---------------------------------------------------------------------------
# simulation — the problems' inline §6.6/§6.6.5 tests via run_pde_tests.
# ---------------------------------------------------------------------------
function run_simulation(output_dir, files, verbose)
    cases = Dict{String,Any}()
    for (case_dir, manifest) in discover_manifests("simulation", files)
        case = String(manifest["case"])
        rec = Dict{String,Any}("case" => case, "status" => "ok")
        try
            problem = normpath(joinpath(case_dir, String(manifest["problem"])))
            alg, reltol, abstol = integrator_opts(manifest)
            results = ESS.run_pde_tests(problem;
                model_name=String(manifest["model"]),
                alg=alg, reltol=reltol, abstol=abstol)
            wanted = Set(String.(manifest["tests"]))
            results = [r for r in results if r.test_id in wanted]
            isempty(results) && error("no assertions ran for tests $(collect(wanted))")
            rec["assertions"] = assertion_dicts(results)
            rec["passed"] = all(r.passed for r in results)
            rec["status"] = rec["passed"] ? "ok" : "failed"
        catch err
            rec["status"] = "error"
            rec["message"] = sprint(showerror, err)
        end
        verbose && println("  simulation/$case: $(rec["status"])")
        cases[case] = rec
    end
    return Dict{String,Any}("binding" => "julia", "category" => "simulation", "cases" => cases)
end

# ---------------------------------------------------------------------------
# convergence — one load per manifest resolution entry (§9.7.6 binding site 4),
# error norms vs the problem's own §6.6.5 analytic reference at assert_time.
# ---------------------------------------------------------------------------
function run_convergence(output_dir, files, verbose)
    cases = Dict{String,Any}()
    for (case_dir, manifest) in discover_manifests("convergence", files)
        case = String(manifest["case"])
        rec = Dict{String,Any}("case" => case, "status" => "ok")
        try
            problem = normpath(joinpath(case_dir, String(manifest["problem"])))
            model = String(manifest["model"])
            assert_time = Float64(manifest["assert_time"])
            alg, reltol, abstol = integrator_opts(manifest)
            errors = Vector{Dict{String,Any}}()
            for res in manifest["resolutions"]
                n = Int(res["n"])
                bindings = Dict{String,Int}(String(k) => Int(v)
                                            for (k, v) in pairs(res["bindings"]))
                # A resolution entry MAY name its own problem file (meshes are
                # subsystem refs, which §9.7.6 cannot rebind — the MPAS
                # refinement family ships one thin problem file per level).
                res_problem = haskey(res, "problem") ?
                    normpath(joinpath(case_dir, String(res["problem"]))) : problem
                file = ESS.load(res_problem; metaparameters=bindings)
                m = file.models[model]
                # The analytic reference is the problem's own §6.6.5 assertion
                # at assert_time (references are authored over the grid's
                # geometry templates, expanded at load).
                refs = [a for t in m.tests for a in t.assertions
                        if a.time == assert_time && a.reduce == "L2_error"]
                isempty(refs) && error("problem declares no L2_error assertion at " *
                                       "t=$assert_time to take the reference from")
                a = refs[1]
                sim = ESS.simulate(file, (0.0, assert_time);
                                   alg=alg, reltol=reltol, abstol=abstol,
                                   saveat=[assert_time])
                sim.success || error("solver retcode $(sim.retcode) at n=$n")
                cells = ESS._state_cells(sim.var_map, a.variable, model)
                isempty(cells) && error("state '$(a.variable)' has no cells at n=$n")
                state = sim.u[end]
                field = Float64[state[slot] for (_, slot) in cells]
                ref = ESS.evaluate_cellwise(a.reference, [c for (c, _) in cells])
                row = Dict{String,Any}("n" => n)
                for norm in String.(manifest["norms"])
                    row[norm] = ESS.field_reduce(norm, field; reference=ref)
                end
                push!(errors, row)
                verbose && println("  convergence/$case n=$n: ",
                                   join(["$k=$(row[k])" for k in sort(collect(keys(row))) if k != "n"], " "))
            end
            rec["assert_time"] = assert_time
            rec["errors"] = errors
        catch err
            rec["status"] = "error"
            rec["message"] = sprint(showerror, err)
        end
        verbose && println("  convergence/$case: $(rec["status"])")
        cases[case] = rec
    end
    return Dict{String,Any}("binding" => "julia", "category" => "convergence", "cases" => cases)
end

# ---------------------------------------------------------------------------
# regridding — the fixture's exact-invariant inline tests, the regridded field,
# and the per-pair build-time setup arrays (A_ij / A_j / W_ij) read from the
# official BuildInspection surface (`simulate(...; inspect=…)` →
# `build_evaluator` setup-array observability) — the §5.8 per-pair gates.
# ---------------------------------------------------------------------------

# Fetch one named setup array from the build inspection. Flattening prefixes
# each observed with its owning model ("Regrid.A_ij"), so try the qualified
# name, then the bare name, then a unique ".<name>" suffix match.
function _setup_array(insp, model::String, name::String)
    for k in (model * "." * name, name)
        haskey(insp.setup_arrays, k) && return insp.setup_arrays[k]
    end
    hits = [k for k in keys(insp.setup_arrays) if endswith(k, "." * name)]
    length(hits) == 1 && return insp.setup_arrays[hits[1]]
    error("setup array '$name' not exposed by the build (have: " *
          "$(sort(collect(keys(insp.setup_arrays)))))")
end

# Row-major nested lists for JSON emission ([i][j] like the manifest triples).
_rows(m::AbstractMatrix) = [Float64[m[i, j] for j in 1:size(m, 2)] for i in 1:size(m, 1)]

function run_regridding(output_dir, files, verbose)
    cases = Dict{String,Any}()
    for (case_dir, manifest) in discover_manifests("regridding", files)
        case = String(manifest["case"])
        rec = Dict{String,Any}("case" => case, "status" => "ok")
        try
            fixture = joinpath(case_dir, String(manifest["fixture"]))
            model = String(manifest["model"])
            results = ESS.run_pde_tests(fixture; model_name=model,
                                        alg=ODE.Tsit5(), reltol=1e-10, abstol=1e-12)
            rec["assertions"] = assertion_dicts(results)
            rec["passed"] = !isempty(results) && all(r.passed for r in results)
            # regrid_state integrates the constant regridded field from 0 over
            # [0,1], so state(1) IS the regridded field F_tgt.
            file = ESS.load(fixture)
            insp = ESS.BuildInspection()
            sim = ESS.simulate(file, (0.0, 1.0); alg=ODE.Tsit5(),
                               reltol=1e-10, abstol=1e-12, saveat=[1.0],
                               inspect=insp)
            sim.success || error("solver retcode $(sim.retcode)")
            for var in ("regrid_state", "pou_state", "cons_state")
                cells = ESS._state_cells(sim.var_map, var, model)
                rec[var * "_at_1"] = Float64[sim.u[end][slot] for (_, slot) in cells]
            end
            # Per-pair setup arrays (manifest §5.8 gates): the raw overlap-area
            # matrix, its filtered row-sums, and the normalized weights — the
            # build-once geometry the evaluator materializes at setup, emitted
            # verbatim (row-major [i][j], 1-based like the manifest triples).
            rec["A_ij"] = _rows(_setup_array(insp, model, "A_ij"))
            rec["A_j"] = collect(Float64, _setup_array(insp, model, "A_j"))
            rec["W_ij"] = _rows(_setup_array(insp, model, "W_ij"))
            rec["status"] = rec["passed"] ? "ok" : "failed"
        catch err
            rec["status"] = "error"
            rec["message"] = sprint(showerror, err)
        end
        verbose && println("  regridding/$case: $(rec["status"])")
        cases[case] = rec
    end
    return Dict{String,Any}("binding" => "julia", "category" => "regridding", "cases" => cases)
end

# ---------------------------------------------------------------------------
# reprojection — evaluate the library's public templates at every golden point
# through the official pathway: a runner-BUILT invocation document (the
# manifest's `fixture: null` contract) importing the library, loaded so §9.7.3
# body composition inlines the template bodies, then `evaluate_expr` per point.
# ---------------------------------------------------------------------------
function _reproj_wrapper_doc(params)
    crs = Dict{String,Any}(k => Float64(params[k])
                           for k in ("lat_1", "lat_2", "lat_0", "lon_0", "R"))
    mkapply(tpl, args...) = Dict{String,Any}(
        "op" => "apply_expression_template", "args" => Any[], "name" => tpl,
        "bindings" => Base.merge(Dict{String,Any}(a => a for a in args), crs))
    vars = Dict{String,Any}(
        "lon" => Dict("type" => "parameter", "units" => "deg", "default" => 0.0),
        "lat" => Dict("type" => "parameter", "units" => "deg", "default" => 0.0),
        "x" => Dict("type" => "parameter", "units" => "m", "default" => 0.0),
        "y" => Dict("type" => "parameter", "units" => "m", "default" => 0.0),
        "fwd_x" => Dict("type" => "observed", "units" => "m",
            "expression" => mkapply("lambert_conformal_forward_x", "lon", "lat")),
        "fwd_y" => Dict("type" => "observed", "units" => "m",
            "expression" => mkapply("lambert_conformal_forward_y", "lon", "lat")),
        "inv_lon" => Dict("type" => "observed", "units" => "deg",
            "expression" => mkapply("lambert_conformal_inverse_lon", "x", "y")),
        "inv_lat" => Dict("type" => "observed", "units" => "deg",
            "expression" => mkapply("lambert_conformal_inverse_lat", "x", "y")))
    return Dict{String,Any}("esm" => "0.8.0",
        "metadata" => Dict{String,Any}("name" => "lambert_conformal_eval",
            "description" => "Runner-built template invocation (manifest fixture: null)."),
        "models" => Dict{String,Any}("Reproject" => Dict{String,Any}(
            "expression_template_imports" => Any[Dict{String,Any}("ref" => "lambert_conformal.esm")],
            "variables" => vars, "equations" => Any[])))
end

function run_reprojection(output_dir, files, verbose)
    cases = Dict{String,Any}()
    for (case_dir, manifest) in discover_manifests("reprojection", files)
        case = String(manifest["case"])
        rec = Dict{String,Any}("case" => case, "status" => "ok")
        try
            lib = normpath(joinpath(case_dir, String(manifest["library"])))
            gold = JSON3.read(read(joinpath(case_dir, String(manifest["golden"])), String))
            exprs = Dict{String,Dict{String,Any}}()
            for (setname, params) in pairs(gold["parameter_sets"])
                doc = _reproj_wrapper_doc(params)
                file = ESS.load(IOBuffer(JSON3.write(doc)); base_path=dirname(lib))
                m = file.models["Reproject"]
                exprs[String(setname)] = Dict{String,Any}(
                    k => m.variables[k].expression
                    for k in ("fwd_x", "fwd_y", "inv_lon", "inv_lat"))
            end
            fwd = Any[]
            for pt in gold["forward"]
                e = exprs[String(pt["set"])]
                b = Dict("lon" => Float64(pt["lon"]), "lat" => Float64(pt["lat"]))
                push!(fwd, Dict{String,Any}("set" => String(pt["set"]),
                    "lon" => Float64(pt["lon"]), "lat" => Float64(pt["lat"]),
                    "x" => ESS.evaluate_expr(e["fwd_x"], b),
                    "y" => ESS.evaluate_expr(e["fwd_y"], b)))
            end
            inv = Any[]
            for pt in gold["inverse"]
                e = exprs[String(pt["set"])]
                b = Dict("x" => Float64(pt["x"]), "y" => Float64(pt["y"]))
                push!(inv, Dict{String,Any}("set" => String(pt["set"]),
                    "x" => Float64(pt["x"]), "y" => Float64(pt["y"]),
                    "lon" => ESS.evaluate_expr(e["inv_lon"], b),
                    "lat" => ESS.evaluate_expr(e["inv_lat"], b)))
            end
            rt = Any[]
            for pt in gold["roundtrip"]
                e = exprs[String(pt["set"])]
                b = Dict("lon" => Float64(pt["lon"]), "lat" => Float64(pt["lat"]))
                xa = ESS.evaluate_expr(e["fwd_x"], b)
                ya = ESS.evaluate_expr(e["fwd_y"], b)
                b2 = Dict("x" => xa, "y" => ya)
                push!(rt, Dict{String,Any}("set" => String(pt["set"]),
                    "lon" => Float64(pt["lon"]), "lat" => Float64(pt["lat"]),
                    "lon_rt" => ESS.evaluate_expr(e["inv_lon"], b2),
                    "lat_rt" => ESS.evaluate_expr(e["inv_lat"], b2)))
            end
            rec["forward"] = fwd
            rec["inverse"] = inv
            rec["roundtrip"] = rt
        catch err
            rec["status"] = "error"
            rec["message"] = sprint(showerror, err)
        end
        verbose && println("  reprojection/$case: $(rec["status"])")
        cases[case] = rec
    end
    return Dict{String,Any}("binding" => "julia", "category" => "reprojection", "cases" => cases)
end

# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
function main(args)
    output_dir, categories, files, verbose = parse_cli(args)
    isempty(categories) && (categories = copy(ALL_CATEGORIES))
    if !isempty(files)
        # --files overrides category-based discovery: run every category the
        # named manifests belong to (still filtered by --categories if given).
        file_cats = unique(String[String(JSON3.read(read(isabspath(f) ? f :
            joinpath(ESD_ROOT, f), String))["category"]) for f in files])
        categories = [c for c in categories if c in file_cats]
    end
    mkpath(output_dir)
    runners = Dict("ast" => run_ast, "simulation" => run_simulation,
                   "convergence" => run_convergence, "regridding" => run_regridding,
                   "reprojection" => run_reprojection)
    summary = Dict{String,Any}("binding" => "julia",
                               "runner" => "scripts/runners/run-julia.jl",
                               "categories" => Dict{String,Any}())
    exit_code = 0
    for cat in categories
        haskey(runners, cat) || error("unknown category '$cat' (known: $(ALL_CATEGORIES))")
        verbose && println("category: $cat")
        res = runners[cat](output_dir, files, verbose)
        open(joinpath(output_dir, "julia_$(cat)_results.json"), "w") do io
            _write_sorted(io, _norm(res), 0)
            print(io, "\n")
        end
        n_cases = length(res["cases"])
        n_bad = count(r -> r["status"] != "ok", values(res["cases"]))
        n_bad > 0 && (exit_code = 1)
        summary["categories"][cat] = Dict{String,Any}("cases" => n_cases, "failures" => n_bad)
        println("julia/$cat: $n_cases case(s), $n_bad failure(s)")
    end
    open(joinpath(output_dir, "julia_summary.json"), "w") do io
        _write_sorted(io, _norm(summary), 0)
        print(io, "\n")
    end
    return exit_code
end

exit(main(ARGS))
