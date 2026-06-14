import JSON
import EarthSciSerialization
import SciMLBase

"""
    build_ode_problem(esm_path; grid_ref="") -> (ODEProblem, var_map)

Load a PDE component `.esm` file, optionally merge a Grid Discretization
Descriptor (GDD) from `grid_ref`, run the ESS canonical discretization
pipeline, and return a ready-to-solve `SciMLBase.ODEProblem` together
with the state-name → index `var_map::Dict{String,Int}`.

The GDD provides the spatial grid (resolution, extent, boundary conditions)
and the discretization rules to inject into the ESM. No solver is invoked
inside this constructor.

`grid_ref` may be an absolute path or a path relative to `esm_path`'s
directory.

If any model in the ESM declares `initial_conditions.type == "expression"`,
the expression is evaluated at each cell using physical cell-centre
coordinates `(index - 0.5) * spacing` for each spatial dimension.  The
resulting per-cell values are passed to `build_evaluator` via its
`initial_conditions` keyword, so `prob.u0` carries the spatially-varying
field without any post-construction mutation by the caller.
"""
function build_ode_problem(esm_path::AbstractString;
                           grid_ref::AbstractString = "")
    esm = _load_json_mutable(esm_path)

    if !isempty(grid_ref)
        gdd_path = isabspath(grid_ref) ? grid_ref :
                   joinpath(dirname(abspath(esm_path)), grid_ref)
        _merge_gdd!(esm, gdd_path)
    end

    expr_ics = _eval_expression_ics(esm)

    disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)
    f!, u0, p, tspan, var_map = EarthSciSerialization.build_evaluator(disc;
                                    initial_conditions = expr_ics)

    prob = SciMLBase.ODEProblem(f!, u0, tspan, p)
    return prob, var_map
end

# ---------------------------------------------------------------------------
# Expression IC helpers
# ---------------------------------------------------------------------------

# Evaluate an EarthSciSerialization Expr at given variable→Float64 bindings.
# Supports: NumExpr, IntExpr, VarExpr, OpExpr (sin, cos, exp, +, -, *, /, ^).
function _eval_expr(expr::EarthSciSerialization.NumExpr, bindings::Dict{String,Float64})
    return expr.value
end
function _eval_expr(expr::EarthSciSerialization.IntExpr, bindings::Dict{String,Float64})
    return Float64(expr.value)
end
function _eval_expr(expr::EarthSciSerialization.VarExpr, bindings::Dict{String,Float64})
    v = get(bindings, expr.name, nothing)
    v !== nothing && return v
    error("Unbound variable '$(expr.name)' in expression IC (bound: $(sort!(collect(keys(bindings)))))")
end
function _eval_expr(expr::EarthSciSerialization.OpExpr, bindings::Dict{String,Float64})
    op   = expr.op
    args = [_eval_expr(a, bindings) for a in expr.args]
    if     op == "+"   return length(args) == 1 ? args[1] : sum(args)
    elseif op == "*"   return prod(args)
    elseif op == "-"   return length(args) == 1 ? -args[1] : args[1] - args[2]
    elseif op == "/"   return args[1] / args[2]
    elseif op == "^"   return args[1]^args[2]
    elseif op == "sin" return sin(args[1])
    elseif op == "cos" return cos(args[1])
    elseif op == "exp" return exp(args[1])
    else
        error("Unsupported operator '$(op)' in expression IC")
    end
end

# Read expression ICs from all models in a merged ESM dict and evaluate each
# variable's IC at every grid cell using physical cell-centre coordinates.
# Returns a Dict{String,Float64} keyed by "varname[i]" / "varname[i,j]".
function _eval_expression_ics(esm::Dict{String,Any})::Dict{String,Float64}
    result = Dict{String,Float64}()
    grids  = get(esm, "grids", Dict{String,Any}())

    for (_, mspec) in get(esm, "models", Dict{String,Any}())
        ic_spec = get(mspec, "initial_conditions", nothing)
        ic_spec === nothing && continue
        get(ic_spec, "type", "") == "expression" || continue

        grid_name = String(get(mspec, "grid", ""))
        grid = get(grids, grid_name, nothing)
        grid === nothing && continue

        # Dimension sizes from the merged grid.
        dim_size = Dict{String,Int}()
        for d in get(grid, "dimensions", Any[])
            dim_size[String(d["name"])] = Int(d["size"])
        end

        # Spacing from "d<dimname>" parameters (set by _inject_grids! from GDD).
        vars = get(mspec, "variables", Dict{String,Any}())
        dim_spacing = Dict{String,Float64}()
        for (pname, vspec) in vars
            get(vspec, "type", "") == "parameter" || continue
            dflt = get(vspec, "default", nothing)
            dflt === nothing && continue
            pstr = String(pname)
            if length(pstr) >= 2 && pstr[1] == 'd'
                dname = pstr[2:end]
                haskey(dim_size, dname) && (dim_spacing[dname] = Float64(dflt))
            end
        end
        # Fallback: infer spacing from 1/N for any dimension without a "d<dim>" param.
        for (dname, N) in dim_size
            haskey(dim_spacing, dname) || (dim_spacing[dname] = 1.0 / N)
        end

        for (var_name, expr_json) in get(ic_spec, "values", Dict{String,Any}())
            vstr  = String(var_name)
            vspec = get(vars, vstr, nothing)
            vspec === nothing && continue
            shape = get(vspec, "shape", nothing)
            (shape === nothing || isempty(shape)) && continue
            shape_strs = String.(shape)

            ic_expr = EarthSciSerialization.parse_expression(expr_json)

            if length(shape_strs) == 1
                d1 = shape_strs[1]
                N1 = dim_size[d1]; h1 = dim_spacing[d1]
                for i in 1:N1
                    bindings = Dict{String,Float64}(d1 => (i - 0.5) * h1)
                    result["$(vstr)[$i]"] = _eval_expr(ic_expr, bindings)
                end
            elseif length(shape_strs) == 2
                d1 = shape_strs[1]; d2 = shape_strs[2]
                N1 = dim_size[d1]; h1 = dim_spacing[d1]
                N2 = dim_size[d2]; h2 = dim_spacing[d2]
                for i in 1:N1, j in 1:N2
                    bindings = Dict{String,Float64}(
                        d1 => (i - 0.5) * h1,
                        d2 => (j - 0.5) * h2,
                    )
                    result["$(vstr)[$i,$j]"] = _eval_expr(ic_expr, bindings)
                end
            end
        end
    end
    result
end

# ---------------------------------------------------------------------------
# GDD merge helpers
# ---------------------------------------------------------------------------

function _merge_gdd!(esm::Dict{String, Any}, gdd_path::AbstractString)
    gdd = _load_json_mutable(gdd_path)

    gdd_grids = get(gdd, "grids", nothing)
    gdd_grids !== nothing && _inject_grids!(esm, gdd_grids, gdd_path)

    gdd_discs = get(gdd, "discretizations", nothing)
    gdd_discs !== nothing && _inject_rules!(esm, gdd_discs, gdd_path)
end

function _inject_grids!(esm::Dict{String, Any}, gdd_grids, gdd_path::AbstractString)
    esm_grids = get!(esm, "grids", Dict{String, Any}())

    for (domain_key, domain_spec) in gdd_grids
        domain_name = String(domain_key)
        spatial = get(domain_spec, "spatial", Dict())
        bcs     = get(domain_spec, "boundary_conditions", Any[])

        periodic_dims = Set{String}()
        for bc in bcs
            get(bc, "type", "") == "periodic" || continue
            for d in get(bc, "dimensions", Any[])
                push!(periodic_dims, String(d))
            end
        end

        sorted_axes = sort!(String[String(k) for k in keys(spatial)])
        dims = Dict{String, Any}[]
        spacing_vals = Dict{String, Float64}()

        for axis_name in sorted_axes
            axis_spec = spatial[axis_name]
            h  = Float64(axis_spec["grid_spacing"])
            lo = Float64(get(axis_spec, "min", 0.0))
            hi = Float64(get(axis_spec, "max", 1.0))
            N  = round(Int, (hi - lo) / h)
            push!(dims, Dict{String, Any}(
                "name"     => axis_name,
                "size"     => N,
                "periodic" => axis_name in periodic_dims,
                "spacing"  => "uniform",
            ))
            spacing_vals["d$(axis_name)"] = h
        end

        esm_grids[domain_name] = Dict{String, Any}(
            "family"     => "cartesian",
            "dimensions" => dims,
        )

        for (_, mspec) in get(esm, "models", Dict())
            get(mspec, "grid", "") == domain_name || continue
            vars = get(mspec, "variables", Dict())
            for (pname, hval) in spacing_vals
                var_spec = get(vars, pname, nothing)
                if var_spec !== nothing && get(var_spec, "type", "") == "parameter"
                    var_spec["default"] = hval
                end
            end
        end
    end
end

function _inject_rules!(esm::Dict{String, Any}, gdd_discs, gdd_path::AbstractString)
    rules = get!(esm, "rules", Any[])
    discs = get!(esm, "discretizations", Dict{String, Any}())

    for (rule_key, rule_ref) in gdd_discs
        rname = String(rule_key)
        ref   = get(rule_ref, "ref", nothing)

        spec = if ref !== nothing
            ref_path = joinpath(dirname(gdd_path), String(ref))
            rule_doc = _load_json_mutable(ref_path)
            rule_doc["discretizations"][rname]
        else
            rule_ref
        end

        stencil = get(spec, "stencil", nothing)
        if stencil isa AbstractVector && !isempty(stencil)
            first_entry = stencil[1]
            sel   = get(first_entry, "selector", nothing)
            kind  = sel !== nothing ? String(get(sel, "kind", "")) : ""
            if kind == "cartesian"
                scheme, use_rule = lower_stencil_to_scheme(rname, spec)
                discs[rname] = scheme
                push!(rules, use_rule)
            else
                repl = lower_stencil_to_canonical_replacement(spec)
                push!(rules, Dict{String, Any}(
                    "name"        => rname,
                    "pattern"     => spec["applies_to"],
                    "replacement" => repl,
                ))
            end
        elseif haskey(spec, "replacement")
            push!(rules, Dict{String, Any}(
                "name"        => rname,
                "pattern"     => spec["applies_to"],
                "replacement" => spec["replacement"],
            ))
        end
    end
end

function _load_json_mutable(path::AbstractString)::Dict{String, Any}
    return JSON.parse(read(path, String))
end
