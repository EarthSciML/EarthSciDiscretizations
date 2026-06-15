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
                           grid_ref::AbstractString = "",
                           reader_fn = nothing)
    esm = _load_json_mutable(esm_path)

    loaded_grids = Dict{String,Any}()

    if !isempty(grid_ref)
        gdd_path = isabspath(grid_ref) ? grid_ref :
                   joinpath(dirname(abspath(esm_path)), grid_ref)
        gdd = _load_json_mutable(gdd_path)
        if _has_curvilinear_domain(gdd)
            return _build_path_b(esm_path, gdd)
        end
        loaded_grids = _merge_gdd!(esm, gdd_path; reader_fn = reader_fn)
    end

    # Path A: rule-engine discretize pipeline.
    expr_ics = _eval_expression_ics(esm)

    disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)

    # Pre-bind unstructured connectivity into equations (esd-aij Route A).
    isempty(loaded_grids) || _prebind_unstructured!(disc, loaded_grids, expr_ics)

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

function _merge_gdd!(esm::Dict{String, Any}, gdd_path::AbstractString;
                     reader_fn = nothing)
    gdd = _load_json_mutable(gdd_path)

    loaded = Dict{String,Any}()
    gdd_grids = get(gdd, "grids", nothing)
    gdd_grids !== nothing &&
        (loaded = _inject_grids!(esm, gdd_grids, gdd_path; reader_fn = reader_fn))

    gdd_discs = get(gdd, "discretizations", nothing)
    gdd_discs !== nothing && _inject_rules!(esm, gdd_discs, gdd_path)

    return loaded
end

function _inject_grids!(esm::Dict{String, Any}, gdd_grids, gdd_path::AbstractString;
                        reader_fn = nothing)
    esm_grids = get!(esm, "grids", Dict{String, Any}())
    loaded = Dict{String,Any}()

    for (domain_key, domain_spec) in gdd_grids
        domain_name = String(domain_key)
        # Read grid family from GDD domain_spec; default "cartesian" for back-compat.
        family = String(get(domain_spec, "family", "cartesian"))

        if family == "mpas"
            dims, spacing_vals, grid = _inject_grids_mpas(domain_spec, gdd_path;
                                                            reader_fn = reader_fn)
            grid !== nothing && (loaded[domain_name] = grid)
            cell_widths_by_axis = Dict{String, Vector{Float64}}()
        elseif family == "duo"
            dims, spacing_vals, grid = _inject_grids_duo(domain_spec, gdd_path)
            grid !== nothing && (loaded[domain_name] = grid)
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
    return loaded
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
# GDD domain_spec must carry "n_cells" (integer). Optional "n_edges" is
# stored but not yet consumed by ESS discretize (pends ESS/esm-bpr).
# When reader_fn is provided, loads MpasGrid from loader.path and returns
# it as a third value so _prebind_unstructured! can extract connectivity.
function _inject_grids_mpas(domain_spec, gdd_path::AbstractString;
                             reader_fn = nothing)
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

    n_edges_raw = get(domain_spec, "n_edges", nothing)
    if n_edges_raw !== nothing
        push!(dims, Dict{String, Any}(
            "name"     => "n_edges",
            "size"     => Int(n_edges_raw),
            "periodic" => false,
            "spacing"  => "unstructured",
        ))
    end

    # Load MpasGrid when reader_fn is provided; otherwise skip (MPAS pre-bind
    # is deferred until a reader_fn is supplied at build_ode_problem call time).
    grid = nothing
    loader_spec = get(domain_spec, "loader", nothing)
    if loader_spec !== nothing && reader_fn !== nothing
        loader_path_abs = let
            raw = String(get(loader_spec, "path", ""))
            isabspath(raw) ? raw : joinpath(dirname(gdd_path), raw)
        end
        grid = build_mpas_grid(
            loader   = Dict("path"   => loader_path_abs,
                            "reader" => String(get(loader_spec, "reader", "mpas_mesh")),
                            "check"  => String(get(loader_spec, "check",  "strict"))),
            reader_fn = reader_fn,
        )
    end

    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals, grid
end

# Build dims for DUO icosahedral triangular grids.
# GDD domain_spec must carry "n_cells" (integer = 20 * 4^level). When a
# "loader" block is present, loads DuoGrid via build_duo_grid and returns it
# as a third value so _prebind_unstructured! can extract connectivity tables.
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

    # Load the DUO grid from the builtin icosahedral loader if specified.
    grid = nothing
    loader_spec = get(domain_spec, "loader", nothing)
    if loader_spec !== nothing
        loader_path = String(get(loader_spec, "path", ""))
        isempty(loader_path) || (grid = build_duo_grid(loader = loader_spec))
    end

    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals, grid
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
            # Dispatch table: selector kind → lowering path.
            # Path-A/scheme: cartesian (ESS native scheme expansion via lower_stencil_to_scheme).
            # Path-A/scheme: reduction/indirect (ESS unstructured expansion via lower_stencil_to_scheme; ess-t0z).
            # Path-A/replacement: latlon, vertical (literal axis names; lower_stencil_to_replacement).
            # Path-A/canonical: arakawa and other $-axis kinds (lower_stencil_to_canonical_replacement).
            # Path-B seam: cubed_sphere falls to canonical; lower_stencil_to_canonical_replacement
            #   throws for plural selectors, making the gap visible (ESS/esm-57f tracks the extension).
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

# ---------------------------------------------------------------------------
# Unstructured connectivity pre-binding (esd-aij)
# ---------------------------------------------------------------------------

# Return (ctables, scalar_arrs) for a DuoGrid:
# ctables maps name → Matrix{Int} or Vector{Int} (integer connectivity);
# scalar_arrs maps name → Vector{Float64} (per-cell float geometry).
function _extract_connectivity(grid::DuoGrid)
    ctables = Dict{String,Any}(
        "cell_neighbors" => grid.cell_neighbors,  # 3 × n_cells, 1-based, 0-based k
    )
    scalar_arrs = Dict{String,Vector{Float64}}(
        "area" => Float64.(grid.area),
    )
    return ctables, scalar_arrs
end

function _extract_connectivity(grid::MpasGrid)
    mesh = grid.mesh
    ctables = Dict{String,Any}(
        "cells_on_cell"  => mesh.cells_on_cell,   # max_edges × n_cells, 1-based
        "edges_on_cell"  => mesh.edges_on_cell,   # max_edges × n_cells, 1-based
        "n_edges_on_cell" => mesh.n_edges_on_cell, # n_cells,   1-based cell → valence
    )
    scalar_arrs = Dict{String,Vector{Float64}}(
        "area_cell" => mesh.area_cell,
        "dv_edge"   => mesh.dv_edge,
        "dc_edge"   => mesh.dc_edge,
    )
    return ctables, scalar_arrs
end

_grid_n_cells(g::DuoGrid)  = length(g.area)
_grid_n_cells(g::MpasGrid) = g.mesh.n_cells

# Substitute a string variable name with a concrete integer in a dict AST.
# Recurses into "args" lists, "expr" sub-bodies, and range bound expressions
# inside "ranges" dicts (keys are loop-var names and are preserved unchanged).
function _dict_subst(node, var::String, val::Int)
    node isa AbstractString && return String(node) == var ? val : node
    node isa Number          && return node
    node isa AbstractVector  && return Any[_dict_subst(x, var, val) for x in node]
    node isa AbstractDict || return node
    out = Dict{String,Any}()
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            out[ks] = Any[_dict_subst(x, var, val) for x in v]
        elseif ks == "expr" && v !== nothing
            out[ks] = _dict_subst(v, var, val)
        elseif ks == "ranges" && v isa AbstractDict
            # Recurse into range bound expressions; preserve range var name keys.
            new_ranges = Dict{String,Any}()
            for (rv, rbound) in v
                new_ranges[String(rv)] = rbound isa AbstractVector ?
                    Any[_dict_subst(b, var, val) for b in rbound] :
                    _dict_subst(rbound, var, val)
            end
            out[ks] = new_ranges
        else
            out[ks] = v
        end
    end
    return out
end

# Return true if node (a disc expression dict) contains a reduction arrayop.
function _has_reduction(node)::Bool
    node isa AbstractDict || return false
    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "arrayop"
        get(node, "reduce", nothing) !== nothing && return true
    end
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            any(_has_reduction(x) for x in v) && return true
        elseif ks == "expr" && v !== nothing
            _has_reduction(v) && return true
        end
    end
    return false
end

# Evaluate a range bound expression that may reference 1D integer connectivity
# lookups (e.g. index(n_edges_on_cell, c) - 1) to a concrete Int.
# Returns nothing when the expression is not statically evaluable.
function _eval_range_end(node, ctables::Dict{String,Any})::Union{Int,Nothing}
    node isa Integer     && return Int(node)
    node isa AbstractFloat && return Int(round(node))
    node isa AbstractDict || return nothing
    op   = get(node, "op", nothing)
    op isa AbstractString || return nothing
    args = get(node, "args", Any[])
    ops  = String(op)
    if ops == "-" && length(args) == 2
        a = _eval_range_end(args[1], ctables)
        b = _eval_range_end(args[2], ctables)
        (a === nothing || b === nothing) && return nothing
        return a - b
    elseif ops == "+" && !isempty(args)
        total = 0
        for a in args
            v = _eval_range_end(a, ctables); v === nothing && return nothing
            total += v
        end
        return total
    elseif ops == "index" && length(args) >= 2
        tname = args[1]
        tname isa AbstractString || return nothing
        tbl = get(ctables, String(tname), nothing)
        tbl isa Vector{Int} || return nothing
        idx = args[2]
        idx isa Integer || return nothing
        return tbl[Int(idx)]  # 1-based lookup
    end
    return nothing
end

# Replace index(table, int[, int]) lookups with concrete integers in a dict
# expression, where table is in ctables and all supplied index arguments are
# concrete integers (after prior _dict_subst passes).
# Convention: 1D Vector{Int} → 1-based; 2D Matrix{Int}(rows=k+1,cols=c) →
#   first index is 1-based column (the outer cell c), second is 0-based row k.
function _resolve_int_indices(node, ctables::Dict{String,Any})
    node isa Number         && return node
    node isa AbstractString && return node
    node isa AbstractVector && return Any[_resolve_int_indices(x, ctables) for x in node]
    node isa AbstractDict || return node
    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "index"
        args_raw = get(node, "args", Any[])
        if length(args_raw) >= 2
            tname = args_raw[1]
            if tname isa AbstractString && haskey(ctables, String(tname))
                tbl = ctables[String(tname)]
                if length(args_raw) == 2 && args_raw[2] isa Integer && tbl isa Vector{Int}
                    return tbl[Int(args_raw[2])]
                elseif length(args_raw) == 3 && args_raw[2] isa Integer &&
                       args_raw[3] isa Integer && tbl isa Matrix{Int}
                    # index(mat, c, k): c is 1-based outer cell, k is 0-based neighbor slot
                    return tbl[Int(args_raw[3]) + 1, Int(args_raw[2])]
                end
            end
        end
    end
    out = Dict{String,Any}()
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            out[ks] = Any[_resolve_int_indices(x, ctables) for x in v]
        elseif ks == "expr" && v !== nothing
            out[ks] = _resolve_int_indices(v, ctables)
        else
            out[ks] = v
        end
    end
    return out
end

# Find the first genuine reduction arrayop (has "op"="arrayop" AND "reduce") in a
# disc expression dict.  Returns the node or nothing.
function _find_outer_reduction(node)
    node isa AbstractDict || return nothing
    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "arrayop" &&
       get(node, "reduce", nothing) !== nothing
        return node
    end
    for v in get(node, "args", Any[])
        r = _find_outer_reduction(v); r !== nothing && return r
    end
    inner = get(node, "expr", nothing)
    inner !== nothing && (r = _find_outer_reduction(inner); r !== nothing && return r)
    return nothing
end

# Replace `index(operand_name, <any concrete index>)` with 1.0 in a dict AST.
# Used to strip the state-variable operand from per-k body terms so that the
# remaining expression is the pure coefficient.
function _replace_operand_index_with_one(node, operand_name::String)
    node isa Number         && return node
    node isa AbstractString && return node
    node isa AbstractVector && return Any[_replace_operand_index_with_one(x, operand_name) for x in node]
    node isa AbstractDict || return node
    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "index"
        args = get(node, "args", Any[])
        if !isempty(args) && args[1] isa AbstractString && String(args[1]) == operand_name
            return 1.0
        end
    end
    out = Dict{String,Any}()
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            out[ks] = Any[_replace_operand_index_with_one(x, operand_name) for x in v]
        elseif ks == "expr" && v !== nothing
            out[ks] = _replace_operand_index_with_one(v, operand_name)
        else
            out[ks] = v
        end
    end
    return out
end

# Compute sum_k(coeff_k) from the outer reduction's body by substituting the
# state-variable operand with 1.0 in each per-k expansion.  This gives the
# indirect-entry coeff when the stencil uses the same per-k expression for both
# the reduction and the indirect (diagonal) term.
# Returns a dict/number, or nothing if the range cannot be evaluated.
function _compute_indirect_coeff(outer_reduction, ctables::Dict{String,Any},
                                  operand_name::String)
    inner_body = get(outer_reduction, "expr", nothing)
    inner_body === nothing && return nothing
    ranges = get(outer_reduction, "ranges", nothing)
    ranges isa AbstractDict && !isempty(ranges) || return nothing
    k_var  = String(first(keys(ranges)))
    k_range = ranges[k_var]
    k_range isa AbstractVector && length(k_range) == 2 || return nothing
    k_lo = k_range[1] isa Integer ? Int(k_range[1]) : nothing
    k_hi = _eval_range_end(k_range[2], ctables)
    (k_lo === nothing || k_hi === nothing) && return nothing
    reduce_op = String(get(outer_reduction, "reduce", "+"))
    terms = Any[]
    for k_val in k_lo:k_hi
        body_k = _dict_subst(inner_body, k_var, k_val)
        body_k = _resolve_int_indices(body_k, ctables)
        coeff_k = _replace_operand_index_with_one(body_k, operand_name)
        push!(terms, coeff_k)
    end
    isempty(terms) && return 0.0
    length(terms) == 1 && return terms[1]
    return Dict{String,Any}("op" => reduce_op, "args" => Any[terms...])
end

# Replace bare "arrayop" MARKER nodes (op="arrayop" but no "reduce" field) in a
# dict AST with `replacement`.  These markers arise because ESS's _parse_expr only
# reads "op" and "args" from coeff expressions, stripping "reduce"/"expr"/"ranges"
# from any nested reduction arrayop.  The replacement is the concrete computed value
# for the indirect coefficient.
function _replace_marker_arrayops(node, replacement)
    node isa Number         && return node
    node isa AbstractString && return node
    node isa AbstractVector && return Any[_replace_marker_arrayops(x, replacement) for x in node]
    node isa AbstractDict || return node
    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "arrayop" &&
       get(node, "reduce", nothing) === nothing
        return replacement
    end
    out = Dict{String,Any}()
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            out[ks] = Any[_replace_marker_arrayops(x, replacement) for x in v]
        elseif ks == "expr" && v !== nothing
            out[ks] = _replace_marker_arrayops(v, replacement)
        else
            out[ks] = v
        end
    end
    return out
end

# Expand all reduction arrayop nodes in a disc expression dict into explicit
# sums using ctables.  Must be called after _dict_subst has substituted the
# outer cell index, so range bounds can be evaluated with _eval_range_end.
function _expand_reductions(node, ctables::Dict{String,Any})
    node isa Number         && return node
    node isa AbstractString && return node
    node isa AbstractVector && return Any[_expand_reductions(x, ctables) for x in node]
    node isa AbstractDict || return node

    op = get(node, "op", nothing)
    if op isa AbstractString && String(op) == "arrayop"
        reduce_op = get(node, "reduce", nothing)
        if reduce_op !== nothing
            inner_body = get(node, "expr", nothing)
            inner_body === nothing && return node
            ranges = get(node, "ranges", nothing)
            ranges isa AbstractDict && !isempty(ranges) || return node
            k_var  = String(first(keys(ranges)))
            k_range = ranges[k_var]
            (k_range isa AbstractVector && length(k_range) == 2) || return node
            k_lo = k_range[1] isa Integer ? Int(k_range[1]) : nothing
            k_hi = _eval_range_end(k_range[2], ctables)
            (k_lo === nothing || k_hi === nothing) && return node
            terms = Any[]
            for k_val in k_lo:k_hi
                body_k = _dict_subst(inner_body, k_var, k_val)
                body_k = _resolve_int_indices(body_k, ctables)
                body_k = _expand_reductions(body_k, ctables)
                push!(terms, body_k)
            end
            isempty(terms) && return 0.0
            length(terms) == 1 && return terms[1]
            return Dict{String,Any}("op" => String(reduce_op), "args" => Any[terms...])
        end
    end

    out = Dict{String,Any}()
    for (k, v) in node
        ks = String(k)
        if ks == "args" && v isa AbstractVector
            out[ks] = Any[_expand_reductions(x, ctables) for x in v]
        elseif ks == "expr" && v !== nothing
            out[ks] = _expand_reductions(v, ctables)
        else
            out[ks] = v
        end
    end
    return out
end

# Pre-bind unstructured grid connectivity into the disc equations (esd-aij).
#
# For each model whose grid is a loaded DuoGrid or MpasGrid:
#   1. Inject scalar float arrays (area, dv_edge, …) as constant state
#      variables in the disc model and their per-cell values into expr_ics.
#   2. Replace the single outer-arrayop equation per variable with N explicit
#      D(index(v, c)) = concrete_sum equations (one per cell) whose connectivity
#      lookups have been resolved to concrete integers.
#   3. Inject zero ICs for any shape=["n_cells"] state variables not yet in
#      expr_ics so that build_evaluator's _detect_array_vars can find them.
#
# This transforms disc so that build_evaluator's _resolve_indices handles all
# index ops with compile-time-constant indices — no change to ESS required.
function _prebind_unstructured!(disc::Dict{String,Any},
                                 loaded_grids::Dict{String,Any},
                                 expr_ics::Dict{String,Float64})
    models_disc = get(disc, "models", nothing)
    models_disc isa AbstractDict || return

    for (_, mdisc_any) in models_disc
        mdisc_any isa AbstractDict || continue
        mdisc = mdisc_any
        grid_name = String(get(mdisc, "grid", ""))
        grid = get(loaded_grids, grid_name, nothing)
        grid === nothing && continue

        ctables, scalar_arrs = _extract_connectivity(grid)
        n_cells = _grid_n_cells(grid)

        # 1. Inject scalar float arrays as constant state variables + ICs.
        vars_disc = get(mdisc, "variables", nothing)
        if vars_disc === nothing
            mdisc["variables"] = Dict{String,Any}()
            vars_disc = mdisc["variables"]
        end
        for (arr_name, vals) in scalar_arrs
            vars_disc[arr_name] = Dict{String,Any}(
                "type"  => "state",
                "shape" => Any["n_cells"],
            )
            for (c, v) in enumerate(vals)
                expr_ics["$(arr_name)[$c]"] = Float64(v)
            end
        end

        # 2. Pre-unroll outer-arrayop equations.
        eqs_raw = get(mdisc, "equations", nothing)
        eqs_raw isa AbstractVector || continue
        new_eqs = Any[]

        for eq in eqs_raw
            eq isa AbstractDict || (push!(new_eqs, eq); continue)
            lhs = get(eq, "lhs", nothing)
            rhs = get(eq, "rhs", nothing)

            # Only unroll equations whose LHS is an outer arrayop with a D expr.
            lhs isa AbstractDict &&
                String(get(lhs, "op", "")) == "arrayop" || (push!(new_eqs, eq); continue)
            rhs isa AbstractDict                          || (push!(new_eqs, eq); continue)

            out_idxs_raw = get(lhs, "output_idx", nothing)
            out_idxs_raw isa AbstractVector && !isempty(out_idxs_raw) ||
                (push!(new_eqs, eq); continue)
            outer_var    = String(out_idxs_raw[1])
            outer_ranges = get(lhs, "ranges", nothing)
            outer_ranges isa AbstractDict               || (push!(new_eqs, eq); continue)
            haskey(outer_ranges, outer_var)             || (push!(new_eqs, eq); continue)
            o_range = outer_ranges[outer_var]
            o_range isa AbstractVector && length(o_range) == 2 ||
                (push!(new_eqs, eq); continue)
            o_lo = o_range[1] isa Integer ? Int(o_range[1]) : nothing
            o_hi = _eval_range_end(o_range[2], ctables)
            (o_lo === nothing || o_hi === nothing) && (push!(new_eqs, eq); continue)

            # Extract LHS variable name and wrt from the inner D expression.
            lhs_body = get(lhs, "expr", nothing)
            lhs_var  = nothing
            wrt      = "t"
            if lhs_body isa AbstractDict && String(get(lhs_body, "op", "")) == "D"
                wrt = String(get(lhs_body, "wrt", "t"))
                d_args = get(lhs_body, "args", Any[])
                if !isempty(d_args)
                    inner = d_args[1]
                    if inner isa AbstractDict &&
                       String(get(inner, "op", "")) == "index"
                        ia = get(inner, "args", Any[])
                        !isempty(ia) && ia[1] isa AbstractString &&
                            (lhs_var = String(ia[1]))
                    end
                end
            end
            lhs_var === nothing && (push!(new_eqs, eq); continue)

            # Extract the RHS body from the outer lift arrayop.
            rhs_body = get(rhs, "expr", nothing)
            rhs_body === nothing && (push!(new_eqs, eq); continue)
            _has_reduction(rhs_body) || (push!(new_eqs, eq); continue)

            # Expand: one per-cell equation for each outer cell c.
            # ESS scheme_expansion uses "c" as target letter for cell_center schemes;
            # arrayop lifting uses a positional name ("i" for 1-D). Substitute both
            # so no free variable remains in the per-cell equation body.
            for c in o_lo:o_hi
                rhs_c = _dict_subst(rhs_body, outer_var, c)
                outer_var == "c" || (rhs_c = _dict_subst(rhs_c, "c", c))
                # ESS's _parse_expr strips "reduce"/"expr"/"ranges" from any nested
                # reduction arrayop inside an indirect-entry coeff.  The disc therefore
                # contains bare MARKER nodes {"op":"arrayop","args":[var_names...]}
                # (no "reduce") that _expand_reductions cannot expand.  Replace each
                # such MARKER with the concrete indirect coeff: sum_k(coeff_k), computed
                # from the outer genuine-reduction body with the state operand set to 1.
                outer_red = _find_outer_reduction(rhs_c)
                if outer_red !== nothing
                    indirect_coeff = _compute_indirect_coeff(outer_red, ctables, lhs_var)
                    if indirect_coeff !== nothing
                        rhs_c = _replace_marker_arrayops(rhs_c, indirect_coeff)
                    end
                end
                rhs_expanded = _expand_reductions(rhs_c, ctables)
                push!(new_eqs, Dict{String,Any}(
                    "lhs" => Dict{String,Any}(
                        "op"   => "D",
                        "args" => Any[Dict{String,Any}(
                            "op"   => "index",
                            "args" => Any[lhs_var, c],
                        )],
                        "wrt"  => wrt,
                    ),
                    "rhs" => rhs_expanded,
                ))
            end
        end

        mdisc["equations"] = new_eqs

        # 3. Inject zero ICs for shape=["n_cells"] state variables with no ICs
        #    so build_evaluator's _detect_array_vars can discover them as array vars.
        for (vname_any, vmeta_any) in vars_disc
            vname  = String(vname_any)
            vmeta_any isa AbstractDict || continue
            String(get(vmeta_any, "type", "state")) == "state" || continue
            shape = get(vmeta_any, "shape", nothing)
            shape isa AbstractVector && length(shape) == 1 &&
                String(shape[1]) == "n_cells" || continue
            any(startswith(String(k), "$(vname)[") for k in keys(expr_ics)) && continue
            for c in 1:n_cells
                expr_ics["$(vname)[$c]"] = 0.0
            end
        end
    end
end
