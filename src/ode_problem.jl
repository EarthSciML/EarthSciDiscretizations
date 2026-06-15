import JSON
import EarthSciSerialization
import SciMLBase
import ModelingToolkit

# Grid families whose discretize path goes through ESS PDESystem (Path B).
const _CURVILINEAR_FAMILIES = Set{String}(["latlon", "cubed_sphere"])

# Computational axis pair for each curvilinear family (xi_axis, eta_axis).
const _CURVILINEAR_AXES = Dict{String, Tuple{Symbol, Symbol}}(
    "latlon"       => (:lon, :lat),
    "cubed_sphere" => (:xi,  :eta),
)

"""
    build_ode_problem(esm_path; grid_ref="") -> (ODEProblem, var_map)

Load a PDE component `.esm` file, optionally merge a Grid Discretization
Descriptor (GDD) from `grid_ref`, run the ESS canonical discretization
pipeline, and return a ready-to-solve `SciMLBase.ODEProblem` together
with the state-name → index `var_map::Dict{String,Int}`.

When the GDD specifies a curvilinear grid family (`latlon` or `cubed_sphere`),
the function routes through the ESS PDESystem pipeline (Path B): the `.esm`
is loaded via `EarthSciSerialization.load`, flattened to a `PDESystem`, and
discretized by `EarthSciSerialization.discretize(sys, grid)`.  Otherwise the
existing rule-engine pipeline (Path A) is used.  No solver is invoked inside
this constructor.

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
        gdd = _load_json_mutable(gdd_path)
        if _has_curvilinear_domain(gdd)
            return _build_path_b(esm_path, gdd)
        end
        _merge_gdd!(esm, gdd_path)
    end

    # Path A: rule-engine discretize pipeline.
    expr_ics = _eval_expression_ics(esm)

    disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)
    f!, u0, p, tspan, var_map = EarthSciSerialization.build_evaluator(disc;
                                    initial_conditions = expr_ics)

    prob = SciMLBase.ODEProblem(f!, u0, tspan, p)
    return prob, var_map
end

# ---------------------------------------------------------------------------
# Path B: PDESystem → ESS.discretize(sys, grid) curvilinear pipeline (esd-ncg)
# ---------------------------------------------------------------------------

function _has_curvilinear_domain(gdd::Dict{String,Any})::Bool
    for (_, domain_spec) in get(gdd, "grids", Dict{String,Any}())
        family = String(get(domain_spec, "family", "cartesian"))
        family in _CURVILINEAR_FAMILIES && return true
    end
    return false
end

function _build_path_b(esm_path::AbstractString, gdd::Dict{String,Any})
    esmfile = EarthSciSerialization.load(esm_path)
    flat    = EarthSciSerialization.flatten(esmfile)
    sys     = ModelingToolkit.PDESystem(flat)

    # Locate the first curvilinear grid spec in the GDD.
    family      = ""
    domain_spec = Dict{String,Any}()
    for (_, dspec) in get(gdd, "grids", Dict{String,Any}())
        f = String(get(dspec, "family", "cartesian"))
        if f in _CURVILINEAR_FAMILIES
            family      = f
            domain_spec = Dict{String,Any}(String(k) => v for (k, v) in dspec)
            break
        end
    end

    grid              = _construct_curvilinear_grid(family, domain_spec)
    xi_axis, eta_axis = _CURVILINEAR_AXES[family]

    prob = EarthSciSerialization.discretize(sys, grid;
                                            xi_axis  = xi_axis,
                                            eta_axis = eta_axis)

    # Best-effort var_map: one entry per cell of each dependent variable.
    N = EarthSciSerialization.n_cells(grid)
    var_map = Dict{String,Int}()
    for (dv_idx, dv) in enumerate(sys.dvs)
        nm   = String(Symbolics.tosymbol(dv, escape=false))
        base = split(nm, "(")[1]
        for i in 1:N
            var_map["$(base)[$i]"] = (dv_idx - 1) * N + i
        end
    end

    return prob, var_map
end

function _construct_curvilinear_grid(family::String, domain_spec::Dict{String,Any})
    if family == "latlon"
        R       = Float64(get(domain_spec, "R", 1.0))
        spatial = Dict{String,Any}(String(k) => v
                                   for (k, v) in get(domain_spec, "spatial", Dict()))
        lon_sp  = Dict{String,Any}(String(k) => v for (k, v) in spatial["lon"])
        lat_sp  = Dict{String,Any}(String(k) => v for (k, v) in spatial["lat"])
        nlon = round(Int, (Float64(get(lon_sp, "max",  π))    - Float64(get(lon_sp, "min", -π)))    / Float64(lon_sp["grid_spacing"]))
        nlat = round(Int, (Float64(get(lat_sp, "max",  π/2))  - Float64(get(lat_sp, "min", -π/2))) / Float64(lat_sp["grid_spacing"]))
        return _latlon(; nlon = nlon, nlat = nlat, R = R)
    elseif family == "cubed_sphere"
        Nc = Int(get(domain_spec, "Nc", 4))
        R  = Float64(get(domain_spec, "R", 1.0))
        return CubedSphereGrid(Nc; R = R)
    else
        error("_construct_curvilinear_grid: unsupported family '$family'")
    end
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

        # Pick up plain numeric per-cell ICs injected by _inject_grids!
        # (e.g. "dz[k]" entries added directly to the IC dict for non-uniform grids).
        # These must reach build_evaluator so _discover_array_cells finds the array cells.
        for (key, val) in ic_spec
            skey = String(key)
            match(r"^([^\[]+)\[[0-9,]+\]$", skey) === nothing && continue
            val isa Real || continue
            result[skey] = Float64(val)
        end

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

        # For non-uniform grids, reconstruct physical cell centres from the per-cell
        # width ICs injected by _inject_grids! (e.g. "dz[k]" entries in ic_spec).
        # When present, these override the uniform-(k-0.5)*h coordinate for that axis.
        dim_centres = Dict{String, Vector{Float64}}()
        for (dname, N) in dim_size
            pname = "d$(dname)"  # e.g. "dz"
            widths = Float64[]
            for k in 1:N
                val = get(ic_spec, "$(pname)[$k]", nothing)
                val isa Real || (empty!(widths); break)
                push!(widths, Float64(val))
            end
            if length(widths) == N
                cumwidths = cumsum(widths)
                dim_centres[dname] = [cumwidths[k] - widths[k] / 2.0 for k in 1:N]
            end
        end

        for (var_name, expr_json) in get(ic_spec, "values", Dict{String,Any}())
            vstr  = String(var_name)
            vspec = get(vars, vstr, nothing)
            vspec === nothing && continue
            shape = get(vspec, "shape", nothing)
            (shape === nothing || isempty(shape)) && continue
            shape_strs = String.(shape)

            # Variable location determines the physical coordinate of each index.
            # face_x: index i sits at x = (i-1)*h (left face of cell i).
            # face_y: index j sits at y = (j-1)*h (bottom face of cell j).
            # All other locations (cell_center, vertex, …): (i-0.5)*h.
            location = String(get(vspec, "location", "cell_center"))

            ic_expr = EarthSciSerialization.parse_expression(expr_json)

            if length(shape_strs) == 1
                d1 = shape_strs[1]
                N1 = dim_size[d1]; h1 = dim_spacing[d1]
                centres1 = get(dim_centres, d1, nothing)
                c1 = location == "face_x" ? -1.0 : -0.5
                for i in 1:N1
                    z1 = centres1 !== nothing ? centres1[i] : (i + c1) * h1
                    bindings = Dict{String,Float64}(d1 => z1)
                    result["$(vstr)[$i]"] = _eval_expr(ic_expr, bindings)
                end
            elseif length(shape_strs) == 2
                d1 = shape_strs[1]; d2 = shape_strs[2]
                N1 = dim_size[d1]; h1 = dim_spacing[d1]
                N2 = dim_size[d2]; h2 = dim_spacing[d2]
                centres1 = get(dim_centres, d1, nothing)
                centres2 = get(dim_centres, d2, nothing)
                c1 = location == "face_x" ? -1.0 : -0.5
                c2 = location == "face_y" ? -1.0 : -0.5
                for i in 1:N1, j in 1:N2
                    z1 = centres1 !== nothing ? centres1[i] : (i + c1) * h1
                    z2 = centres2 !== nothing ? centres2[j] : (j + c2) * h2
                    bindings = Dict{String,Float64}(
                        d1 => z1,
                        d2 => z2,
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
        # Read grid family from GDD domain_spec; default "cartesian" for back-compat.
        family = String(get(domain_spec, "family", "cartesian"))

        if family == "mpas"
            dims, spacing_vals = _inject_grids_mpas(domain_spec, gdd_path)
            cell_widths_by_axis = Dict{String, Vector{Float64}}()
        elseif family == "duo"
            dims, spacing_vals = _inject_grids_duo(domain_spec, gdd_path)
            cell_widths_by_axis = Dict{String, Vector{Float64}}()
        else
            dims, spacing_vals, cell_widths_by_axis = _inject_grids_spatial(domain_spec)
        end

        esm_grids[domain_name] = Dict{String, Any}(
            "family"     => family,
            "dimensions" => dims,
        )

        for (_, mspec) in get(esm, "models", Dict())
            get(mspec, "grid", "") == domain_name || continue
            vars = get(mspec, "variables", Dict())
            ics  = get(mspec, "initial_conditions", nothing)

            # Inject uniform spacing as parameter defaults.
            for (pname, hval) in spacing_vals
                var_spec = get(vars, pname, nothing)
                if var_spec !== nothing && get(var_spec, "type", "") == "parameter"
                    var_spec["default"] = hval
                end
            end

            # Inject per-cell widths as initial conditions for array state variables
            # named d<axis> (e.g. dz for axis z) declared with shape [axis].
            for (axis_name, widths) in cell_widths_by_axis
                vname = "d$(axis_name)"
                var_spec = get(vars, vname, nothing)
                var_spec === nothing && continue
                get(var_spec, "type", "") == "state" || continue
                shape = get(var_spec, "shape", nothing)
                shape isa AbstractVector && length(shape) == 1 &&
                    String(shape[1]) == axis_name || continue
                # Inject per-cell width ICs: vname[1]=w1, vname[2]=w2, ...
                if ics === nothing
                    mspec["initial_conditions"] = Dict{String, Any}()
                    ics = mspec["initial_conditions"]
                end
                if !isa(ics, AbstractDict)
                    ics = Dict{String, Any}()
                    mspec["initial_conditions"] = ics
                end
                for (k, w) in enumerate(widths)
                    ics["$(vname)[$k]"] = w
                end
            end
        end
    end
end

# Build dims + spacing_vals + cell_widths_by_axis from spatial axis specs
# (cartesian/latlon/vertical/arakawa). Non-uniform axes use "levels" (face positions);
# uniform axes use "grid_spacing" + optional "min"/"max".
function _inject_grids_spatial(domain_spec)
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
    cell_widths_by_axis = Dict{String, Vector{Float64}}()

    for axis_name in sorted_axes
        axis_spec = spatial[axis_name]
        if haskey(axis_spec, "levels")
            # Non-uniform spacing: levels specifies face positions [z0, z1, ..., zN]
            levels = Float64.(axis_spec["levels"])
            N = length(levels) - 1
            N >= 1 || error("_inject_grids_spatial: 'levels' for axis '$axis_name' must have ≥ 2 entries")
            widths = [levels[k+1] - levels[k] for k in 1:N]
            all(w > 0 for w in widths) || error("_inject_grids_spatial: 'levels' must be strictly increasing for axis '$axis_name'")
            push!(dims, Dict{String, Any}(
                "name"     => axis_name,
                "size"     => N,
                "periodic" => axis_name in periodic_dims,
                "spacing"  => "nonuniform",
            ))
            cell_widths_by_axis[axis_name] = widths
            spacing_vals["d$(axis_name)"] = (levels[end] - levels[1]) / N  # average spacing
        else
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
    end

    return dims, spacing_vals, cell_widths_by_axis
end

# Build dims for MPAS unstructured Voronoi grids.
# GDD domain_spec must carry "n_cells" (integer). Optional "n_edges" and
# "max_edges" are stored but not yet consumed by ESS discretize (pends
# ESS/esm-bpr unstructured selector dispatch). "loader" is preserved for
# the future NetCDF mesh-loading hook; no file I/O is performed here.
function _inject_grids_mpas(domain_spec, gdd_path::AbstractString)
    n_cells_raw = get(domain_spec, "n_cells", nothing)
    n_cells_raw === nothing && throw(
        ArgumentError(
            "MPAS GDD domain_spec missing required field 'n_cells'. " *
            "Add \"n_cells\": <integer> to the GDD's grids.<domain> block."
        )
    )
    n_cells = Int(n_cells_raw)
    n_cells > 0 || throw(
        ArgumentError("MPAS GDD 'n_cells' must be a positive integer (got $n_cells)")
    )

    dims = [Dict{String, Any}(
        "name"     => "n_cells",
        "size"     => n_cells,
        "periodic" => false,
        "spacing"  => "unstructured",
    )]

    # n_edges is metadata; included in the grid block for reference but not
    # consumed by the current ESD/ESS pipeline (pends ESS/esm-bpr).
    n_edges_raw = get(domain_spec, "n_edges", nothing)
    if n_edges_raw !== nothing
        push!(dims, Dict{String, Any}(
            "name"     => "n_edges",
            "size"     => Int(n_edges_raw),
            "periodic" => false,
            "spacing"  => "unstructured",
        ))
    end

    # No uniform spacing to inject; MPAS rules read geometry from connectivity
    # tables (area_cell, dc_edge, dv_edge) not from scalar parameters.
    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals
end

# Build dims for DUO icosahedral triangular grids.
# GDD domain_spec must carry "n_cells" (integer = 20 * 4^level). The DUO grid
# has constant valence 3 (all triangular cells), so no n_edges_on_cell array
# is required. "loader" is preserved for the future build_duo_grid integration;
# no grid construction is performed here (rules are gated behind ESS/esm-bpr).
function _inject_grids_duo(domain_spec, gdd_path::AbstractString)
    n_cells_raw = get(domain_spec, "n_cells", nothing)
    n_cells_raw === nothing && throw(
        ArgumentError(
            "DUO GDD domain_spec missing required field 'n_cells'. " *
            "Add \"n_cells\": <integer> (= 20 * 4^level) to the GDD's grids.<domain> block."
        )
    )
    n_cells = Int(n_cells_raw)
    n_cells > 0 || throw(
        ArgumentError("DUO GDD 'n_cells' must be a positive integer (got $n_cells)")
    )

    dims = [Dict{String, Any}(
        "name"     => "n_cells",
        "size"     => n_cells,
        "periodic" => false,
        "spacing"  => "unstructured",
    )]

    # No uniform spacing to inject; DUO rules read geometry from the cell
    # area array (area) and cell_neighbors connectivity table.
    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals
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
            # cubed_sphere entries use plural "selectors" (array); fall back to first selector.
            if sel === nothing
                sels = get(first_entry, "selectors", nothing)
                if sels isa AbstractVector && !isempty(sels) && sels[1] isa AbstractDict
                    sel = sels[1]
                end
            end
            kind  = sel !== nothing ? String(get(sel, "kind", "")) : ""
            # Dispatch table: selector kind → lowering path.
            # Path-A/scheme: cartesian (ESS native scheme expansion via lower_stencil_to_scheme).
            # Path-A/scheme: reduction/indirect (ESS unstructured expansion via lower_stencil_to_scheme; ess-t0z).
            # Path-A/replacement: latlon, vertical (literal axis names; lower_stencil_to_replacement).
            # Path-A/replacement: cubed_sphere (plural selectors; lower_stencil_to_replacement keeps axis names).
            # Path-A/canonical: arakawa and other $-axis kinds (lower_stencil_to_canonical_replacement).
            if kind == "cartesian"
                scheme, use_rule = lower_stencil_to_scheme(rname, spec)
                discs[rname] = scheme
                push!(rules, use_rule)
            elseif kind in ("latlon", "vertical")
                lowered = lower_stencil_to_replacement(spec)
                replacement = lowered["replacement"]
                if kind == "latlon"
                    # Translate literal geographic axis names ("lat", "lon", …) in
                    # the replacement to canonical arrayop loop-variable names
                    # ("i", "j", …) so that build_evaluator's {i→val, j→val}
                    # substitution can resolve all index expressions.  The latlon
                    # lowerer already sorts axis names alphabetically, giving a
                    # stable position→canonical mapping.
                    axis_names_set = Set{String}()
                    for entry in stencil
                        entry isa AbstractDict || continue
                        sel_e = get(entry, "selector", nothing)
                        sel_e isa AbstractDict || continue
                        ax = get(sel_e, "axis", nothing)
                        ax isa AbstractString && !startswith(String(ax), "\$") &&
                            push!(axis_names_set, String(ax))
                    end
                    sorted_axes = sort!(collect(axis_names_set))
                    _CANONICAL = ("i", "j", "k", "l", "m")
                    axis_map = Dict{String,String}(
                        ax => _CANONICAL[d] for (d, ax) in enumerate(sorted_axes)
                    )
                    replacement = _subst_axis_names(replacement, axis_map)
                end
                push!(rules, Dict{String, Any}(
                    "name"        => rname,
                    "pattern"     => spec["applies_to"],
                    "replacement" => replacement,
                ))
            elseif kind in ("reduction", "indirect")
                scheme, use_rule = lower_stencil_to_scheme(rname, spec)
                discs[rname] = scheme
                push!(rules, use_rule)
            elseif kind == "cubed_sphere"
                lowered = lower_stencil_to_replacement(spec)
                push!(rules, Dict{String, Any}(
                    "name"        => rname,
                    "pattern"     => spec["applies_to"],
                    "replacement" => lowered["replacement"],
                ))
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

function _subst_axis_names(expr, mapping::Dict{String,String})
    if expr isa AbstractString
        s = String(expr)
        return get(mapping, s, s)
    elseif expr isa AbstractDict
        out = Dict{String,Any}()
        for (k, v) in expr
            out[String(k)] = _subst_axis_names(v, mapping)
        end
        return out
    elseif expr isa AbstractVector
        return Any[_subst_axis_names(a, mapping) for a in expr]
    else
        return expr
    end
end

function _load_json_mutable(path::AbstractString)::Dict{String, Any}
    return JSON.parse(read(path, String))
end
