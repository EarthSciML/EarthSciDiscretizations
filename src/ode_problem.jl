import JSON
import EarthSciSerialization
import SciMLBase
import ModelingToolkit

# Grid families whose discretize path goes through ESS PDESystem (Path B).
const _CURVILINEAR_FAMILIES = Set{String}(["latlon"])

# Computational axis pair for each curvilinear family (xi_axis, eta_axis).
const _CURVILINEAR_AXES = Dict{String, Tuple{Symbol, Symbol}}(
    "latlon" => (:lon, :lat),
)

"""
    build_ode_problem(esm_path; grid_ref="", reader_fn=nothing,
                      extra_ics=Dict{String,Float64}()) -> (ODEProblem, var_map)

Load a PDE component `.esm` file, optionally merge a Grid Discretization
Descriptor (GDD) from `grid_ref`, run the ESS canonical discretization
pipeline, and return a ready-to-solve `SciMLBase.ODEProblem` together
with the state-name → index `var_map::Dict{String,Int}`.

When the GDD specifies a curvilinear grid family (`latlon`), the function
routes through the ESS PDESystem pipeline (Path B): the `.esm`
is loaded via `EarthSciSerialization.load`, flattened to a `PDESystem`, and
discretized by `EarthSciSerialization.discretize(sys, grid)`.  Otherwise the
existing rule-engine pipeline (Path A) is used.  No solver is invoked inside
this constructor.

# Arguments
- `grid_ref`: absolute path, or a path relative to `esm_path`'s directory, to
  the GDD to merge before discretizing.
- `reader_fn`: optional mesh reader callback forwarded to the unstructured grid
  loaders (`build_mpas_grid`) for non-builtin MPAS meshes.
- `extra_ics`: additional per-cell numeric initial conditions (keyed
  `"var[i]"` / `"var[i,j]"`) merged into the ICs after expression ICs are
  lowered — used e.g. to inject manufactured-solution fields per cell.

# Expression initial conditions
If a model declares `initial_conditions.type == "expression"`,
`_prepare_expression_ics!` rewrites each entry under `values` into an
`initialization_equation` (`lhs = var`, `rhs = the IC expression`) on the ESM,
and ESS's IC-arrayop engine materializes it: `discretize` turns the `ic` node
into an `arrayop` that evaluates the RHS at every cell, substituting each
spatial-dim symbol with an `index(coord_<dim>, <loop-idx>)` lookup, and
`build_evaluator` evaluates those arrayops into `prob.u0` (explicit numeric ICs
still win).  The `coord_<dim>` arrays are injected here as `const_arrays`.

Sampling honours each variable's staggered `location`: a `cell_center`
variable is sampled at cell centres `(index - 0.5) * spacing`, while a
face-staggered component (`face_x`, `face_y`, `vertex`; e.g. Arakawa C-grid
velocities) has its faced axis sampled at the cell face `(index - 1) * spacing`
via a per-variable `coord_<var>_<dim>` array (nonuniform axes use the per-cell
widths analogously).  Sampling a staggered field at the cell centre would
degrade an otherwise O(h²) operator to O(h).
"""
function build_ode_problem(
        esm_path::AbstractString;
        grid_ref::AbstractString = "",
        reader_fn = nothing,
        extra_ics::Dict{String, Float64} = Dict{String, Float64}()
    )
    esm = _load_json_mutable(esm_path)

    loaded_grids = Dict{String, Any}()

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
    # _prepare_expression_ics! mutates esm: converts type:expression ICs to
    # initialization_equations for the ESS IC-arrayop engine, and returns the
    # numeric per-cell ICs (e.g. dz[k]) and coord const_arrays separately.
    numeric_ics, coord_arrays = _prepare_expression_ics!(esm)
    merge!(numeric_ics, extra_ics)

    disc = EarthSciSerialization.discretize(esm; lift_1d_arrayop = true)

    # Bind unstructured (MPAS/DUO) grid geometry declaratively: the authored
    # nn_diffusion rule references the PRIMITIVE mesh arrays directly, and its
    # FVM-coefficient sub-expression num/(den*area) is a CONST-cadence static
    # partition (RFC semiring-faq-unified-ir §6.1, §9) — every input is a
    # build-time-constant grid array, so ESS folds the coefficient at build time
    # and only the state difference (u[nbr] − u[i]) reaches the per-step hot
    # tree. The host supplies the primitive arrays as `const_arrays`; it does NO
    # coefficient const-fold and NO arrayop rewrite (the retired Route-B
    # _unstructured_const_arrays / _rewrite_unstructured_arrayop!).
    const_arrays = Dict{String, AbstractArray{Float64}}(coord_arrays)
    if !isempty(loaded_grids)
        merge!(const_arrays, _unstructured_grid_const_arrays(loaded_grids))
    end

    f!, u0, p, tspan, var_map = EarthSciSerialization.build_evaluator(
        disc;
        initial_conditions = numeric_ics,
        const_arrays = const_arrays
    )

    prob = SciMLBase.ODEProblem(f!, u0, tspan, p)
    return prob, var_map
end

# ---------------------------------------------------------------------------
# Path B: PDESystem → ESS.discretize(sys, grid) curvilinear pipeline (esd-ncg)
# ---------------------------------------------------------------------------

function _has_curvilinear_domain(gdd::Dict{String, Any})::Bool
    for (_, domain_spec) in get(gdd, "grids", Dict{String, Any}())
        family = String(get(domain_spec, "family", "cartesian"))
        family in _CURVILINEAR_FAMILIES && return true
    end
    return false
end

function _build_path_b(esm_path::AbstractString, gdd::Dict{String, Any})
    esmfile = EarthSciSerialization.load(esm_path)
    flat = EarthSciSerialization.flatten(esmfile)
    sys = ModelingToolkit.PDESystem(flat)

    # Locate the first curvilinear grid spec in the GDD.
    family = ""
    domain_spec = Dict{String, Any}()
    for (_, dspec) in get(gdd, "grids", Dict{String, Any}())
        f = String(get(dspec, "family", "cartesian"))
        if f in _CURVILINEAR_FAMILIES
            family = f
            domain_spec = Dict{String, Any}(String(k) => v for (k, v) in dspec)
            break
        end
    end

    grid = _construct_curvilinear_grid(family, domain_spec)
    xi_axis, eta_axis = _CURVILINEAR_AXES[family]

    prob = EarthSciSerialization.discretize(
        sys, grid;
        xi_axis = xi_axis,
        eta_axis = eta_axis
    )

    # Best-effort var_map: one entry per cell of each dependent variable.
    N = EarthSciSerialization.n_cells(grid)
    var_map = Dict{String, Int}()
    for (dv_idx, dv) in enumerate(sys.dvs)
        nm = String(Symbolics.tosymbol(dv, escape = false))
        base = split(nm, "(")[1]
        for i in 1:N
            var_map["$(base)[$i]"] = (dv_idx - 1) * N + i
        end
    end

    return prob, var_map
end

function _construct_curvilinear_grid(family::String, domain_spec::Dict{String, Any})
    if family == "latlon"
        R = Float64(get(domain_spec, "R", 1.0))
        spatial = Dict{String, Any}(
            String(k) => v
                for (k, v) in get(domain_spec, "spatial", Dict())
        )
        lon_sp = Dict{String, Any}(String(k) => v for (k, v) in spatial["lon"])
        lat_sp = Dict{String, Any}(String(k) => v for (k, v) in spatial["lat"])
        nlon = round(Int, (Float64(get(lon_sp, "max", π)) - Float64(get(lon_sp, "min", -π))) / Float64(lon_sp["grid_spacing"]))
        nlat = round(Int, (Float64(get(lat_sp, "max", π / 2)) - Float64(get(lat_sp, "min", -π / 2))) / Float64(lat_sp["grid_spacing"]))
        return _latlon(; nlon = nlon, nlat = nlat, R = R)
    else
        error("_construct_curvilinear_grid: unsupported family '$family'")
    end
end

# ---------------------------------------------------------------------------
# Expression IC helpers
# ---------------------------------------------------------------------------

# Mutates esm: converts InitialConditions type:expression to initialization_equations
# for the ESS IC-arrayop engine (ess-tt6). Returns:
#   numeric_ics   — Dict{String,Float64} of per-cell numeric ICs (e.g. "dz[k]")
#   coord_arrays  — Dict{String,Vector{Float64}} of coordinate const_arrays
#                   passed to build_evaluator: shared cell-centre "coord_<dim>"
#                   plus per-variable "coord_<var>_<dim>" for face-staggered
#                   components (sampled at their own location, not cell centres).
function _prepare_expression_ics!(esm::Dict{String, Any})
    numeric_ics = Dict{String, Float64}()
    coord_arrays = Dict{String, Vector{Float64}}()
    grids = get(esm, "grids", Dict{String, Any}())

    for (_, mspec) in get(esm, "models", Dict{String, Any}())
        ic_spec = get(mspec, "initial_conditions", nothing)
        ic_spec === nothing && continue

        # Collect numeric per-cell ICs injected by _inject_grids! (e.g. "dz[k]").
        for (key, val) in ic_spec
            skey = String(key)
            match(r"^([^\[]+)\[[0-9,]+\]$", skey) === nothing && continue
            val isa Real || continue
            numeric_ics[skey] = Float64(val)
        end

        get(ic_spec, "type", "") == "expression" || continue

        grid_name = String(get(mspec, "grid", ""))
        grid = get(grids, grid_name, nothing)
        grid === nothing && continue

        dim_size = Dict{String, Int}()
        for d in get(grid, "dimensions", Any[])
            dim_size[String(d["name"])] = Int(d["size"])
        end

        vars = get(mspec, "variables", Dict{String, Any}())
        dim_spacing = Dict{String, Float64}()
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
        for (dname, N) in dim_size
            haskey(dim_spacing, dname) || (dim_spacing[dname] = 1.0 / N)
        end

        # Per-dim cell widths (nonuniform) or `nothing` (uniform, use dim_spacing).
        dim_widths = Dict{String, Union{Vector{Float64}, Nothing}}()
        for (dname, N) in dim_size
            pname = "d$(dname)"
            widths = Float64[]
            for k in 1:N
                val = get(ic_spec, "$(pname)[$k]", nothing)
                val isa Real || (empty!(widths); break)
                push!(widths, Float64(val))
            end
            dim_widths[dname] = length(widths) == N ? widths : nothing
        end

        # Build shared coord_<dim> cell-centre const_arrays. These back the
        # cell-centred (and non-staggered) variables, whose IC RHS keeps its
        # bare spatial symbols for ESS to substitute as index(coord_<dim>, idx).
        for (dname, N) in dim_size
            coord_name = "coord_$dname"
            haskey(coord_arrays, coord_name) && continue
            coord_arrays[coord_name] =
                _coord_axis_array(dim_widths[dname], dim_spacing[dname], N, false)
        end

        # Move expression ICs to initialization_equations for ESS IC-arrayop engine.
        init_eqs = get!(mspec, "initialization_equations", Any[])
        if !(init_eqs isa AbstractVector)
            init_eqs = Any[]
            mspec["initialization_equations"] = init_eqs
        end
        for (var_name, expr_json) in get(ic_spec, "values", Dict{String, Any}())
            vstr = String(var_name)
            vspec = get(vars, vstr, nothing)
            shape = vspec === nothing ? nothing : get(vspec, "shape", nothing)
            loc = vspec === nothing ? "" : String(get(vspec, "location", ""))
            faced = _faced_dims(loc)

            rhs = expr_json
            if !isempty(faced) && shape isa AbstractVector
                # Staggered variable: sample each spatial-dim symbol at this
                # variable's own physical location (face vs cell centre per axis),
                # not the shared cell-centre coord_<dim>. ESS materializes the IC
                # arrayop with loop index `_ARRAYOP_INDEX_NAMES[d]` = (i,j,k,…) by
                # shape position, so the per-variable coord array we inject and the
                # index node we splice must use that same letter.
                idx_letters = ("i", "j", "k", "l", "m", "n")
                coord_subst = Dict{String, Any}()
                for (d, dim_any) in enumerate(shape)
                    d <= length(idx_letters) || break
                    dname = String(dim_any)
                    N = get(dim_size, dname, nothing)
                    N === nothing && continue
                    is_face = dname in faced
                    coord_name = is_face ? "coord_$(vstr)_$dname" : "coord_$dname"
                    if is_face
                        coord_arrays[coord_name] = _coord_axis_array(
                            dim_widths[dname], get(dim_spacing, dname, 1.0 / N), N, true
                        )
                    end
                    coord_subst[dname] = Dict{String, Any}(
                        "op" => "index", "args" => Any[coord_name, idx_letters[d]]
                    )
                end
                rhs = _subst_ic_coords(expr_json, coord_subst)
            end

            push!(init_eqs, Dict{String, Any}("lhs" => vstr, "rhs" => rhs))
        end
        delete!(ic_spec, "type")
        delete!(ic_spec, "values")
    end

    return numeric_ics, coord_arrays
end

# Coordinate samples for one spatial axis of length `N`.
# `widths` non-nothing ⇒ nonuniform (per-cell widths); else uniform spacing `h`.
# `face=true` returns the lower (west/south) cell-face positions; `face=false`
# returns cell centres. Mirrors the (idx-0.5)*h / (idx-1)*h Arakawa convention
# in src/grids/arakawa.jl (no domain-min offset, matching the shared coord_<dim>).
function _coord_axis_array(
        widths::Union{Vector{Float64}, Nothing}, h::Float64,
        N::Int, face::Bool
    )
    if widths !== nothing
        edges = cumsum(widths)            # edges[k] = upper face of cell k
        return face ? [edges[k] - widths[k] for k in 1:N] :
            [edges[k] - widths[k] / 2.0 for k in 1:N]
    else
        return face ? [(k - 1) * h for k in 1:N] : [(k - 0.5) * h for k in 1:N]
    end
end

# Spatial dims that a variable's `location` places on a cell face rather than at
# the cell centre. Arakawa C/D-grid locations: "face_<dim>" stages axis <dim>;
# "vertex" stages every horizontal axis. "cell_center" (or unset) ⇒ none.
function _faced_dims(location::AbstractString)::Set{String}
    loc = String(location)
    startswith(loc, "face_") && return Set{String}([loc[(length("face_") + 1):end]])
    loc == "vertex" && return Set{String}(["x", "y"])
    return Set{String}()
end

# Replace bare spatial-dim leaf strings in `args` positions with their per-cell
# coord-index nodes (the same rewrite ESS's _substitute_coord_syms performs, but
# done here so a staggered variable resolves its own coord array). op/wrt/key
# names are untouched; once rewritten, ESS finds no bare symbol left to swap.
function _subst_ic_coords(node, subst::Dict{String, Any})
    node isa AbstractString && return get(subst, String(node), node)
    node isa AbstractDict || return node
    out = Dict{String, Any}()
    for (k, v) in node
        key = String(k)
        if key == "args" && v isa AbstractVector
            out[key] = Any[_subst_ic_coords(a, subst) for a in v]
        elseif v isa AbstractDict
            out[key] = _subst_ic_coords(v, subst)
        elseif v isa AbstractVector
            out[key] = Any[a isa AbstractDict ? _subst_ic_coords(a, subst) : a for a in v]
        else
            out[key] = v
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# GDD merge helpers
# ---------------------------------------------------------------------------

function _merge_gdd!(
        esm::Dict{String, Any}, gdd_path::AbstractString;
        reader_fn = nothing
    )
    gdd = _load_json_mutable(gdd_path)

    loaded = Dict{String, Any}()
    gdd_grids = get(gdd, "grids", nothing)
    gdd_grids !== nothing &&
        (loaded = _inject_grids!(esm, gdd_grids, gdd_path; reader_fn = reader_fn))

    gdd_discs = get(gdd, "discretizations", nothing)
    gdd_discs !== nothing && _inject_rules!(esm, gdd_discs, gdd_path)

    return loaded
end

function _inject_grids!(
        esm::Dict{String, Any}, gdd_grids, gdd_path::AbstractString;
        reader_fn = nothing
    )
    esm_grids = get!(esm, "grids", Dict{String, Any}())
    loaded = Dict{String, Any}()

    for (domain_key, domain_spec) in gdd_grids
        domain_name = String(domain_key)
        # Read grid family from GDD domain_spec; default "cartesian" for back-compat.
        family = String(get(domain_spec, "family", "cartesian"))

        if family == "mpas"
            dims, spacing_vals, grid = _inject_grids_mpas(
                domain_spec, gdd_path;
                reader_fn = reader_fn
            )
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
            "family" => family,
            "dimensions" => dims,
        )

        for (_, mspec) in get(esm, "models", Dict())
            get(mspec, "grid", "") == domain_name || continue
            vars = get(mspec, "variables", Dict())
            ics = get(mspec, "initial_conditions", nothing)

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
    bcs = get(domain_spec, "boundary_conditions", Any[])

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
            widths = [levels[k + 1] - levels[k] for k in 1:N]
            all(w > 0 for w in widths) || error("_inject_grids_spatial: 'levels' must be strictly increasing for axis '$axis_name'")
            push!(
                dims, Dict{String, Any}(
                    "name" => axis_name,
                    "size" => N,
                    "periodic" => axis_name in periodic_dims,
                    "spacing" => "nonuniform",
                )
            )
            cell_widths_by_axis[axis_name] = widths
            spacing_vals["d$(axis_name)"] = (levels[end] - levels[1]) / N  # average spacing
        else
            h = Float64(axis_spec["grid_spacing"])
            lo = Float64(get(axis_spec, "min", 0.0))
            hi = Float64(get(axis_spec, "max", 1.0))
            N = round(Int, (hi - lo) / h)
            push!(
                dims, Dict{String, Any}(
                    "name" => axis_name,
                    "size" => N,
                    "periodic" => axis_name in periodic_dims,
                    "spacing" => "uniform",
                )
            )
            spacing_vals["d$(axis_name)"] = h
        end
    end

    return dims, spacing_vals, cell_widths_by_axis
end

# Build dims for MPAS unstructured Voronoi grids.
# GDD domain_spec must carry "n_cells" (integer). Optional "n_edges" is
# stored but not yet consumed by ESS discretize (pends ESS/esm-bpr).
# When reader_fn is provided, loads MpasGrid from loader.path and returns
# it as a third value so _unstructured_grid_const_arrays can emit its const_arrays.
function _inject_grids_mpas(
        domain_spec, gdd_path::AbstractString;
        reader_fn = nothing
    )
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

    dims = [
        Dict{String, Any}(
            "name" => "n_cells",
            "size" => n_cells,
            "periodic" => false,
            "spacing" => "unstructured",
        ),
    ]

    n_edges_raw = get(domain_spec, "n_edges", nothing)
    if n_edges_raw !== nothing
        push!(
            dims, Dict{String, Any}(
                "name" => "n_edges",
                "size" => Int(n_edges_raw),
                "periodic" => false,
                "spacing" => "unstructured",
            )
        )
    end

    # Load MpasGrid: builtin Voronoi paths need no reader_fn; all others do.
    grid = nothing
    loader_spec = get(domain_spec, "loader", nothing)
    if loader_spec !== nothing
        loader_path_raw = String(get(loader_spec, "path", ""))
        loader_path_abs = startswith(loader_path_raw, "builtin://") ? loader_path_raw :
            (
                isabspath(loader_path_raw) ? loader_path_raw :
                joinpath(dirname(gdd_path), loader_path_raw)
            )
        R_sphere = Float64(get(loader_spec, "sphere_radius", 6.371e6))
        builtin_level = _parse_builtin_voronoi_level(loader_path_abs)
        if builtin_level !== nothing
            grid = build_mpas_grid(
                loader = Dict(
                    "path" => loader_path_abs,
                    "reader" => "auto",
                    "check" => String(get(loader_spec, "check", "strict"))
                ),
                R = R_sphere,
            )
        elseif reader_fn !== nothing
            grid = build_mpas_grid(
                loader = Dict(
                    "path" => loader_path_abs,
                    "reader" => String(get(loader_spec, "reader", "mpas_mesh")),
                    "check" => String(get(loader_spec, "check", "strict"))
                ),
                reader_fn = reader_fn,
                R = R_sphere,
            )
        end
    end

    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals, grid
end

# Build dims for DUO icosahedral triangular grids.
# GDD domain_spec must carry "n_cells" (integer = 20 * 4^level). When a
# "loader" block is present, loads DuoGrid via build_duo_grid and returns it
# as a third value so _unstructured_grid_const_arrays can emit its const_arrays.
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

    dims = [
        Dict{String, Any}(
            "name" => "n_cells",
            "size" => n_cells,
            "periodic" => false,
            "spacing" => "unstructured",
        ),
    ]

    # Load the DUO grid from the builtin icosahedral loader if specified.
    grid = nothing
    loader_spec = get(domain_spec, "loader", nothing)
    if loader_spec !== nothing
        loader_path = String(get(loader_spec, "path", ""))
        if !isempty(loader_path)
            R_sphere = Float64(get(loader_spec, "sphere_radius", 6.371e6))
            grid = build_duo_grid(loader = loader_spec, R = R_sphere)
        end
    end

    spacing_vals = Dict{String, Float64}()
    return dims, spacing_vals, grid
end

function _inject_rules!(esm::Dict{String, Any}, gdd_discs, gdd_path::AbstractString)
    rules = get!(esm, "rules", Any[])
    discs = get!(esm, "discretizations", Dict{String, Any}())

    for (rule_key, rule_ref) in gdd_discs
        rname = String(rule_key)
        ref = get(rule_ref, "ref", nothing)

        spec = if ref !== nothing
            ref_path = joinpath(dirname(gdd_path), String(ref))
            rule_doc = _load_json_mutable(ref_path)
            rule_doc["discretizations"][rname]
        else
            rule_ref
        end

        if haskey(spec, "replacement")
            replacement = spec["replacement"]
            # Authored latlon replacements use geographic axis names ("lat", "lon").
            # Translate to canonical arrayop loop-variable names ("i", "j") so that
            # build_evaluator's {i→val, j→val} substitution resolves all index expressions.
            # Axes sort alphabetically: "lat" → "i", "lon" → "j" (esd-t4h).
            if String(get(spec, "grid_family", "")) == "latlon"
                replacement = _subst_axis_names(
                    replacement,
                    Dict{String, String}("lat" => "i", "lon" => "j")
                )
            end
            push!(
                rules, Dict{String, Any}(
                    "name" => rname,
                    "pattern" => spec["applies_to"],
                    "replacement" => replacement,
                )
            )
        end
    end
    return
end

function _subst_axis_names(expr, mapping::Dict{String, String})
    if expr isa AbstractString
        s = String(expr)
        return get(mapping, s, s)
    elseif expr isa AbstractDict
        out = Dict{String, Any}()
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
# Unstructured connectivity extraction (geometry source for Route B binding)
# ---------------------------------------------------------------------------

# Return (ctables, scalar_arrs) for a DuoGrid:
# ctables maps name → Matrix{Int} or Vector{Int} (integer connectivity);
# scalar_arrs maps name → Vector{Float64} (per-cell float geometry).
function _extract_connectivity(grid::DuoGrid)
    # Build MPAS Voronoi dual to obtain circumcenter-to-circumcenter (dv_edge)
    # and primal-edge-length (dc_edge) geometry needed for the FVM weights.
    mesh = _duo_voronoi_dual(grid.level; R = grid.R)

    # Invert cells_on_edge: (DUO vertex pair) → MPAS edge index.
    Ne = mesh.n_edges
    edge_idx_duo = Dict{Tuple{Int, Int}, Int}()
    sizehint!(edge_idx_duo, Ne)
    for e in 1:Ne
        v1 = mesh.cells_on_edge[1, e]
        v2 = mesh.cells_on_edge[2, e]
        key = v1 < v2 ? (v1, v2) : (v2, v1)
        edge_idx_duo[key] = e
    end

    # Build edges_on_face (3 × Nc, 0-based k): local edge k of DUO triangle c
    # → MPAS edge index. Edge opposite vertex k+1 (1-based) of each triangle:
    #   k=0 (opposite v1): edge (v2, v3)
    #   k=1 (opposite v2): edge (v3, v1)
    #   k=2 (opposite v3): edge (v1, v2)
    Nc = size(grid.faces, 2)
    edges_on_face = Matrix{Int}(undef, 3, Nc)
    for c in 1:Nc
        v1 = grid.faces[1, c]; v2 = grid.faces[2, c]; v3 = grid.faces[3, c]
        edges_on_face[1, c] = edge_idx_duo[v2 < v3 ? (v2, v3) : (v3, v2)]
        edges_on_face[2, c] = edge_idx_duo[v3 < v1 ? (v3, v1) : (v1, v3)]
        edges_on_face[3, c] = edge_idx_duo[v1 < v2 ? (v1, v2) : (v2, v1)]
    end

    # area_eff[c] = (1/4) Σ_k dc_edge[k] * dv_edge[k] is the correct
    # normalization area for the dc/dv FVM formula on triangular DUO cells.
    # It equals the triangle area for equilateral cells and gives second-order
    # accuracy on quasi-uniform meshes where grid.area (circumscribed-circle
    # area) would not cancel the O(h) balance residual from distorted cells.
    area_eff = Vector{Float64}(undef, Nc)
    for c in 1:Nc
        s = 0.0
        for k in 1:3
            e = edges_on_face[k, c]
            s += mesh.dc_edge[e] * mesh.dv_edge[e]
        end
        area_eff[c] = 0.25 * s
    end

    ctables = Dict{String, Any}(
        "cell_neighbors" => grid.cell_neighbors,  # 3 × n_cells, 1-based, 0-based k
        "edges_on_face" => edges_on_face,         # 3 × n_cells, MPAS edge idx, 0-based k
    )
    scalar_arrs = Dict{String, Vector{Float64}}(
        "area" => area_eff,
        "dv_edge" => mesh.dv_edge,
        "dc_edge" => mesh.dc_edge,
    )
    return ctables, scalar_arrs
end

function _extract_connectivity(grid::MpasGrid)
    mesh = grid.mesh
    ctables = Dict{String, Any}(
        "cells_on_cell" => mesh.cells_on_cell,   # max_edges × n_cells, 1-based
        "edges_on_cell" => mesh.edges_on_cell,   # max_edges × n_cells, 1-based
        "n_edges_on_cell" => mesh.n_edges_on_cell, # n_cells,   1-based cell → valence
    )
    scalar_arrs = Dict{String, Vector{Float64}}(
        "area_cell" => mesh.area_cell,
        "dv_edge" => mesh.dv_edge,
        "dc_edge" => mesh.dc_edge,
    )
    return ctables, scalar_arrs
end

# ---------------------------------------------------------------------------
# Declarative unstructured binding: the authored nn_diffusion rule (DUO + MPAS)
# references the PRIMITIVE mesh arrays directly — the incident-edge metrics
# dc_edge/dv_edge, the per-cell normalisation area, the incident-edge table, and
# the neighbour table. Its FVM-coefficient sub-expression num/(den*area) is a
# CONST-cadence static partition (RFC semiring-faq-unified-ir §6.1, §9): every
# input is a build-time-constant grid array, so ESS folds the coefficient at
# build time and only the state difference (u[nbr] − u[i]) reaches the per-step
# hot tree. The host performs NO coefficient const-fold and NO arrayop rewrite —
# the DEFERRED FAQ-IR cleanup, now whole-grid (retiring the prior Route-B
# _unstructured_const_arrays + _rewrite_unstructured_arrayop!).
#
# Two host-side shaping steps remain — pure index logic, the structured-grid FAQ
# convention of `cartesian_faq.jl` / `arakawa_faq.jl` (geometry rides ESS's
# determinism contract; only the array layout is local):
#
#   1. Connectivity tables are transposed from the grid's (slot, cell) storage
#      to the (cell, slot) layout ESS gathers as `index(table, i, k)` (i = cell,
#      k = slot), and stored as Float64 (ESS resolves an integer gather by
#      rounding the looked-up value).
#   2. Each array is keyed by the symbol its family's authored rule references.
#
# No padding/compaction is needed: the authored reduction's per-cell upper bound
# is the cell's true valence — the constant 3 of the closed triangular DUO mesh,
# or `index(n_edges_on_cell, i)` for MPAS — so the gather iterates exactly the
# real neighbours, which occupy the leading slots `1..valence`, and never reads a
# trailing 0-sentinel padding column. Holds for both closed (MPAS-global, DUO
# icosahedral) and open meshes.

# Emit, for every loaded MPAS/DUO grid, the primitive mesh arrays the authored
# nn_diffusion rule references, keyed by their rule-symbol names. The result is
# merged into the `const_arrays` passed to `build_evaluator`.
function _unstructured_grid_const_arrays(loaded_grids::Dict{String, Any})
    out = Dict{String, AbstractArray{Float64}}()
    for (_, grid) in loaded_grids
        (grid isa DuoGrid || grid isa MpasGrid) || continue
        merge!(out, _grid_primitive_arrays(grid))
    end
    return out
end

# (slot, cell) integer connectivity → (cell, slot) Float64, the layout ESS
# gathers as `index(table, i = cell, k = slot)`.
_cell_slot_table(m::AbstractMatrix) = Float64.(permutedims(m, (2, 1)))

# DUO: the authored rule's coefficient is dc_edge[e] / (dv_edge[e] * area_eff[c])
# and its reduction bound is the constant 3 — the closed icosahedral mesh is
# valence-3, so every slot of cell_neighbors/edges_on_face is a real neighbour.
function _grid_primitive_arrays(grid::DuoGrid)
    ctables, scalar_arrs = _extract_connectivity(grid)
    return Dict{String, AbstractArray{Float64}}(
        "cell_neighbors" => _cell_slot_table(ctables["cell_neighbors"]::Matrix{Int}),
        "edges_on_face" => _cell_slot_table(ctables["edges_on_face"]::Matrix{Int}),
        "area" => scalar_arrs["area"],
        "dv_edge" => scalar_arrs["dv_edge"],
        "dc_edge" => scalar_arrs["dc_edge"],
    )
end

# MPAS: the authored rule's coefficient is dv_edge[e] / (dc_edge[e] *
# area_cell[c]); valence varies (5/6), so its reduction bound is
# `index(n_edges_on_cell, i)` and the gather stops at the cell's real-neighbour
# count, never reading the trailing padding columns of cells_on_cell/edges_on_cell.
function _grid_primitive_arrays(grid::MpasGrid)
    ctables, scalar_arrs = _extract_connectivity(grid)
    return Dict{String, AbstractArray{Float64}}(
        "cells_on_cell" => _cell_slot_table(ctables["cells_on_cell"]::Matrix{Int}),
        "edges_on_cell" => _cell_slot_table(ctables["edges_on_cell"]::Matrix{Int}),
        "n_edges_on_cell" => Float64.(ctables["n_edges_on_cell"]::Vector{Int}),
        "area_cell" => scalar_arrs["area_cell"],
        "dv_edge" => scalar_arrs["dv_edge"],
        "dc_edge" => scalar_arrs["dc_edge"],
    )
end

# LatLon: the covariant FV Laplacian/gradient rules
# (discretizations/finite_volume/covariant_fv_{laplacian,gradient}_latlon.json)
# reference the per-cell curvilinear metric directly, exactly as the MPAS/DUO
# rules reference primitive mesh arrays. The inverse metric g^{ij}, the Jacobian
# J = √det(g) = R²|cos φ|, the J·g^{ij} products (whose centered differences are
# the connection-term corrections), and the coordinate-Jacobian components
# ∂(comp)/∂(target) are all build-time-constant grid data, so the host emits them
# as `const_arrays` keyed by the rule-symbol names. ESS gathers them via
# `index(name, lat, lon)` and const-folds the coefficient; the per-step hot tree
# sees only the state stencil. This is a pure host-side binding hook — no AST
# node, no arrayop footprint, no scheme handler (esd-zk9.2 §6.1). The flat
# per-cell layout (lon-fastest within each lat row) is exactly the order the
# rule's `index(name, lat, lon)` → cell resolution expects, and matches the
# oracle goldens' `c = i + (j-1)·Nx` indexing.
function _grid_primitive_arrays(grid::LatLonGrid)
    ginv = metric_ginv(grid)                # (nc,2,2): [:,1,1]=g^{lonlon}, [:,2,2]=g^{latlat}, [:,1,2]=g^{lonlat}
    jac = metric_jacobian(grid)             # (nc,):    R²|cos φ|
    cj = coord_jacobian(grid, :lon_lat)     # (nc,2,2): ∂(computational)/∂(target)
    nc = n_cells(grid)
    g_xx = Float64[ginv[k, 1, 1] for k in 1:nc]
    g_yy = Float64[ginv[k, 2, 2] for k in 1:nc]
    g_xe = Float64[ginv[k, 1, 2] for k in 1:nc]
    return Dict{String, AbstractArray{Float64}}(
        "g_xx" => g_xx,
        "g_yy" => g_yy,
        "g_xe" => g_xe,
        "invJ" => Float64[1.0 / jac[k] for k in 1:nc],
        "Jg_xx" => Float64[jac[k] * g_xx[k] for k in 1:nc],
        "Jg_yy" => Float64[jac[k] * g_yy[k] for k in 1:nc],
        "Jg_xe" => Float64[jac[k] * g_xe[k] for k in 1:nc],
        "dxi_dt1" => Float64[cj[k, 1, 1] for k in 1:nc],
        "deta_dt1" => Float64[cj[k, 2, 1] for k in 1:nc],
        "dxi_dt2" => Float64[cj[k, 1, 2] for k in 1:nc],
        "deta_dt2" => Float64[cj[k, 2, 2] for k in 1:nc],
    )
end
