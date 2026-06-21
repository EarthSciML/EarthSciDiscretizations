"""
MPAS-SCVT mesh generator host driver (epic `esd-e5m`, bead `esd-e5m.4` / D4).

The external fixed-point **loop** of the declarative MPAS-SCVT mesh generator —
the host analogue of `build_ode_problem`. Where D2 (`lloyd_step.esm`) is the
per-iteration declarative FAQ **step** and D3 (`scvt_voronoi_connectivity`) is the
one-time post-convergence topology **leaf**, this driver supplies the RHS-only
glue between them:

  generators₀, ρ, D1 background  ──┐
    repeat (host loop, NO loop in the IR):                       ┌─ esd-e5m.2 (D2)
      vi  = materialize_value_invention(lloyd_step.esm, …)  ◀────┘  declarative STEP
      cₛ  = vi.centroid                       (density-weighted Lloyd centroid)
      gₙ  = R·cₛ / |cₛ|                        (HOST sphere re-projection, √)
      until ‖gₙ − g‖∞ < tol                    (HOST convergence test)
    conn = scvt_voronoi_connectivity(g; R)     (esd-e5m.3 / D3 LEAF, ONCE)  ◀── topology
    MpasMeshData(conn + dual GEOMETRY)                                          ──▶ output

Per the epic's DECLARATIVE-OR-FAIL contract the per-iteration *step* is the
declarative semiring-FAQ (driven here through `materialize_value_invention`, the
value-invention front-door); the *loop* — the convergence test and the sphere
re-projection (which needs a `sqrt`, outside the build-time value-invention op
set, `LLOYD_STEP_CONTRACT.md` §1) — is host RHS-only code; the *topology* is the
allowed irreducible leaf. No recurrence is ever lowered into the IR: the engine
materialises ONE Lloyd iteration, the host iterates it.

The driver emits the same `MpasMeshData` (`src/grids/mpas.jl`) the imperative
`_duo_voronoi_dual` produces, reusing the identical dual TOPOLOGY
(`voronoi_dual_topology_faq`, wrapped by the D3 leaf) and dual GEOMETRY
(`duo_dual_geometry_faq`) FAQ stack — the ONLY difference from the DUO dual is the
source of the generators (Lloyd convergence under a density `ρ`, not the fixed
icosahedral subdivision) and of the triangulation (the spherical-Delaunay leaf,
not the DUO primal faces). This is the only path supporting variable resolution
via a density function.
"""

# Centroid / mass integrand ASTs (the D1 `background_quadrature.esm` pathway,
# `discretizations/grids/mpas/scvt/regenerate_fixtures.jl`). Every background
# quadrature value is produced THROUGH `eval_coeff` (the single-pathway evaluator),
# byte-identical to the D1 golden — never a binding-local shortcut.
_scvt_mk(op, args...) = Dict{String, Any}("op" => op, "args" => collect(Any, args))
# Unit-sphere primal-cell centroid over the three DUO corner unit vectors a,b,c.
const _SCVT_MX = _scvt_mk("+", "a1", "b1", "c1")
const _SCVT_MY = _scvt_mk("+", "a2", "b2", "c2")
const _SCVT_MZ = _scvt_mk("+", "a3", "b3", "c3")
const _SCVT_NRM = _scvt_mk("sqrt",
    _scvt_mk("+", _scvt_mk("*", _SCVT_MX, _SCVT_MX),
        _scvt_mk("*", _SCVT_MY, _SCVT_MY), _scvt_mk("*", _SCVT_MZ, _SCVT_MZ)))
const _SCVT_UCOMP = (_scvt_mk("/", _SCVT_MX, _SCVT_NRM),
    _scvt_mk("/", _SCVT_MY, _SCVT_NRM), _scvt_mk("/", _SCVT_MZ, _SCVT_NRM))
# FAQ q-1 integrand: bg_mass = ρ · area (background_quadrature.esm).
const _SCVT_MASS_AST = _scvt_mk("*", "rho", "area")

_scvt_corner_binding(V, i, j, k) = Dict{String, Float64}(
    "a1" => V[1, i], "a2" => V[2, i], "a3" => V[3, i],
    "b1" => V[1, j], "b2" => V[2, j], "b3" => V[3, j],
    "c1" => V[1, k], "c2" => V[2, k], "c3" => V[3, k])

"""
    scvt_background_quadrature(level; density=nothing, R=6.371e6)
        -> (bg_coord::Matrix{Float64}, bg_mass::Vector{Float64})

The FIXED background integration point set (bead `esd-e5m.1` / D1) the SCVT step
integrates over, materialised through the single-pathway `eval_coeff`:

- `bg_coord` `(Nc, 3)` — the unit-sphere DUO primal-cell centroids at subdivision
  `level` (the `background_quadrature.esm` quadrature points), `Nc = 10·4^level + 2`.
  Row-major `(cells, space)` to match the `lloyd_step.esm` `bg_coord` parameter.
- `bg_mass` `(Nc,)` — the density-weighted quadrature measure `ρ·area`, where
  `area` is the canonical DUO primal-cell spherical area (`build_duo_grid`) and `ρ`
  is sampled at the quadrature point.

`density === nothing` is the uniform measure `ρ ≡ 1` (SCVT → CVT); otherwise
`density(x, y, z) -> Real` is sampled at each unit `bg_coord` (the
variable-resolution mechanism). The background should be FINER than the generator
set so every generator is attended (`LLOYD_STEP_CONTRACT.md` §4); a background as
coarse as the generators starves some groups (empty-group `den = 0`, handled by
the loop as a held generator).
"""
function scvt_background_quadrature(level::Integer; density = nothing, R::Real = 6.371e6)
    level >= 0 || throw(DomainError(level, "scvt_background_quadrature: level must be ≥ 0"))
    (R > 0 && isfinite(R)) || throw(DomainError(R,
        "scvt_background_quadrature: R must be a positive finite number"))
    g = build_duo_grid(;
        loader = (path = "builtin://icosahedral/$(level)", reader = "builtin_icosahedral"), R = R)
    V, F = duo_subdivide_faq(Float64, Int(level))
    Nc = n_cells(g)
    bg_coord = Matrix{Float64}(undef, Nc, 3)
    bg_mass = Vector{Float64}(undef, Nc)
    @inbounds for c in 1:Nc
        b = _scvt_corner_binding(V, F[1, c], F[2, c], F[3, c])
        x = eval_coeff(_SCVT_UCOMP[1], b)
        y = eval_coeff(_SCVT_UCOMP[2], b)
        z = eval_coeff(_SCVT_UCOMP[3], b)
        bg_coord[c, 1] = x
        bg_coord[c, 2] = y
        bg_coord[c, 3] = z
        rho = density === nothing ? 1.0 : Float64(density(x, y, z))
        bg_mass[c] = eval_coeff(_SCVT_MASS_AST, Dict{String, Float64}("rho" => rho, "area" => g.area[c]))
    end
    return bg_coord, bg_mass
end

# Load + select the declarative SCVT step model (`lloyd_step.esm`, D2) and size its
# index sets to the actual background / generator counts. The model is a pure
# function of its parameters, so it is loaded once and reused across iterations.
function _scvt_load_step_model(n_bg::Int, n_gen::Int)
    REPO_ROOT = dirname(dirname(pathof(@__MODULE__)))
    esm_path = joinpath(REPO_ROOT, "discretizations", "grids", "mpas", "scvt", "lloyd_step.esm")
    doc = JSON.parsefile(esm_path; dicttype = Dict{String, Any})
    mj = EarthSciSerialization._select_model_json(doc, "ScvtLloydStep")
    mj["index_sets"]["cells"]["size"] = n_bg
    mj["index_sets"]["generators"]["size"] = n_gen
    return mj
end

"""
    scvt_lloyd_solve(generators, bg_coord, bg_mass; R=6.371e6, tol=1e-10,
                     max_iters=1000, verbose=false) -> NamedTuple

The external Lloyd / SCVT fixed-point loop (RHS-only — NO loop in the IR). Each
iteration evaluates ONE declarative Lloyd step (`lloyd_step.esm`, D2) through the
value-invention front-door, then applies the host sphere re-projection and the
host convergence test:

1. `vi = materialize_value_invention(lloyd_step.esm, {bg_coord, bg_mass, gen})`
   — the declarative ASSIGN → grouped CENTROID → elementwise MOVE step.
2. `gₙ[g] = R · centroid[g] / |centroid[g]|` — re-project the density-weighted
   centroid onto the sphere of radius `R` (the `sqrt` is host, outside the
   value-invention op set). An UNATTENDED generator (`den[g] = 0`, NaN centroid)
   keeps its previous position.
3. Converged when the per-generator displacement `max_g ‖gₙ[g] − g[g]‖ < tol`.

`generators` is `(3, Nc)` Cartesian seed positions; `bg_coord` `(Nc_bg, 3)` /
`bg_mass` `(Nc_bg,)` are the fixed D1 background (`scvt_background_quadrature`).
Returns `(; generators, iterations, displacement, converged)` with the converged
generators as `(3, Nc)` on the sphere of radius `R`.
"""
function scvt_lloyd_solve(generators::AbstractMatrix{<:Real}, bg_coord::AbstractMatrix{<:Real},
        bg_mass::AbstractVector{<:Real}; R::Real = 6.371e6, tol::Real = 1e-10,
        max_iters::Integer = 1000, verbose::Bool = false)
    size(generators, 1) == 3 || throw(ArgumentError(
        "scvt_lloyd_solve: generators must be (3, Nc); got $(size(generators))"))
    size(bg_coord, 2) == 3 || throw(ArgumentError(
        "scvt_lloyd_solve: bg_coord must be (Nc_bg, 3); got $(size(bg_coord))"))
    size(bg_coord, 1) == length(bg_mass) || throw(ArgumentError(
        "scvt_lloyd_solve: bg_coord rows $(size(bg_coord, 1)) must match bg_mass length $(length(bg_mass))"))
    (R > 0 && isfinite(R)) || throw(DomainError(R,
        "scvt_lloyd_solve: R must be a positive finite number"))
    (tol > 0) || throw(DomainError(tol, "scvt_lloyd_solve: tol must be positive"))
    max_iters >= 1 || throw(DomainError(max_iters, "scvt_lloyd_solve: max_iters must be ≥ 1"))
    Nc = size(generators, 2)

    # The step works on UNIT directions (bg_coord is the unit-sphere quadrature
    # set), so the argmin metric and the centroid live in the same unit space; the
    # radius R is restored only at output and for the downstream geometry. The
    # `lloyd_step.esm` `gen` parameter is row-major (generators, space).
    gen = Matrix{Float64}(undef, Nc, 3)
    @inbounds for i in 1:Nc
        x = Float64(generators[1, i]); y = Float64(generators[2, i]); z = Float64(generators[3, i])
        nrm = sqrt(x * x + y * y + z * z)
        nrm > 0 || throw(ArgumentError("scvt_lloyd_solve: generator $i is the zero vector"))
        gen[i, 1] = x / nrm; gen[i, 2] = y / nrm; gen[i, 3] = z / nrm
    end

    bgc = Matrix{Float64}(bg_coord)
    bgm = Vector{Float64}(bg_mass)
    mj = _scvt_load_step_model(size(bgc, 1), Nc)

    gen_next = Matrix{Float64}(undef, Nc, 3)
    displacement = Inf
    converged = false
    iters = 0
    for it in 1:Int(max_iters)
        iters = it
        vi = EarthSciSerialization.materialize_value_invention(
            mj, Dict("bg_coord" => bgc, "bg_mass" => bgm, "gen" => gen),
            Dict{String, Float64}())
        cx = vi.groups["centroid_x"]; cy = vi.groups["centroid_y"]; cz = vi.groups["centroid_z"]
        maxd = 0.0
        @inbounds for g in 1:Nc
            if isnan(cx[g]) || isnan(cy[g]) || isnan(cz[g])
                # Unattended generator (empty group): hold the previous seed.
                gen_next[g, 1] = gen[g, 1]; gen_next[g, 2] = gen[g, 2]; gen_next[g, 3] = gen[g, 3]
                continue
            end
            nrm = sqrt(cx[g]^2 + cy[g]^2 + cz[g]^2)
            nx = cx[g] / nrm; ny = cy[g] / nrm; nz = cz[g] / nrm
            dx = nx - gen[g, 1]; dy = ny - gen[g, 2]; dz = nz - gen[g, 3]
            d = sqrt(dx * dx + dy * dy + dz * dz)
            d > maxd && (maxd = d)
            gen_next[g, 1] = nx; gen_next[g, 2] = ny; gen_next[g, 3] = nz
        end
        gen, gen_next = gen_next, gen
        displacement = maxd
        verbose && @info "scvt_lloyd_solve" iteration = it displacement = maxd
        if maxd < tol
            converged = true
            break
        end
    end

    R_f = Float64(R)
    out = Matrix{Float64}(undef, 3, Nc)
    @inbounds for i in 1:Nc
        out[1, i] = R_f * gen[i, 1]; out[2, i] = R_f * gen[i, 2]; out[3, i] = R_f * gen[i, 3]
    end
    return (generators = out, iterations = iters, displacement = displacement, converged = converged)
end

# Map each canonical edge (column of `edges`, sorted `(a < b)` generator pair) to
# the ≤ 2 incident Delaunay triangles, edge-aligned with `edges`. Mirrors
# `_duo_edge_cells`, but parameterised by the leaf triangulation + the canonical
# edge numbering (`primal_topology_faq(...).edges`) so the result is index-aligned
# with the `edges_on_cell` the D3 leaf emits.
function _scvt_edge_cells(faces::AbstractMatrix{<:Integer}, edges::AbstractMatrix{<:Integer})
    cells_by_key = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    @inbounds for c in 1:size(faces, 2)
        v1 = Int(faces[1, c]); v2 = Int(faces[2, c]); v3 = Int(faces[3, c])
        for (a, b) in ((v2, v3), (v3, v1), (v1, v2))
            key = a < b ? (a, b) : (b, a)
            prev = get(cells_by_key, key, (0, 0))
            cells_by_key[key] = prev[1] == 0 ? (c, 0) : (prev[1], c)
        end
    end
    Ne = size(edges, 2)
    out = Matrix{Int}(undef, 2, Ne)
    @inbounds for e in 1:Ne
        key = (Int(edges[1, e]), Int(edges[2, e]))
        cc = get(cells_by_key, key, (0, 0))
        out[1, e] = cc[1]; out[2, e] = cc[2]
    end
    return out
end

"""
    build_scvt_mesh(; generators, density=nothing, background_level,
                      R=6.371e6, tol=1e-10, max_iters=1000, verbose=false)
        -> MpasMeshData

The MPAS-SCVT mesh generator host driver (bead `esd-e5m.4` / D4): runs the
external Lloyd fixed-point loop (`scvt_lloyd_solve`) over the declarative D2 step,
calls the D3 spherical-topology leaf (`scvt_voronoi_connectivity`) ONCE at
convergence, and assembles a validated `MpasMeshData`.

Arguments:
- `generators` `(3, Nc)` — Cartesian seed cell centres (`Nc ≥ 4`); the count is
  the number of mesh cells. Icosahedral-dual vertices
  (`duo_subdivide_faq(Float64, level)[1]`) are a natural seed.
- `density` — `nothing` (uniform `ρ ≡ 1` → CVT) or `density(x, y, z) -> Real`
  sampled on the unit sphere (variable resolution; must be positive).
- `background_level` — the D1 background quadrature subdivision level
  (`scvt_background_quadrature`). Choose it FINER than the generator set so every
  generator is attended (`LLOYD_STEP_CONTRACT.md` §4).
- `R`, `tol`, `max_iters`, `verbose` — passed to `scvt_lloyd_solve`.

The returned mesh shares the dual TOPOLOGY/GEOMETRY FAQ stack with the imperative
DUO Voronoi dual (`_duo_voronoi_dual`): cells = generators, dual vertices = the
spherical-Delaunay circumcentres, edges = the Delaunay edges. For `Nc` generators
forming a closed mesh, `n_edges = 3·Nc − 6` and `n_vertices = 2·Nc − 4` (Euler).
"""
function build_scvt_mesh(; generators::AbstractMatrix{<:Real}, density = nothing,
        background_level::Integer, R::Real = 6.371e6, tol::Real = 1e-10,
        max_iters::Integer = 1000, verbose::Bool = false)
    size(generators, 1) == 3 || throw(ArgumentError(
        "build_scvt_mesh: generators must be (3, Nc); got $(size(generators))"))
    Nc = size(generators, 2)
    Nc >= 4 || throw(ArgumentError(
        "build_scvt_mesh: need ≥ 4 generators for a closed spherical mesh; got $Nc"))
    R_f = Float64(R)

    # --- D1: the fixed background quadrature (the integration measure) ----------
    bg_coord, bg_mass = scvt_background_quadrature(background_level; density = density, R = R_f)
    size(bg_coord, 1) >= Nc || @warn(
        "build_scvt_mesh: background ($(size(bg_coord, 1)) points) is not finer than the " *
        "generators ($Nc); some generators may be unattended (held at their seed).")

    # --- D2 step driven by the host fixed-point LOOP (RHS-only) -----------------
    sol = scvt_lloyd_solve(generators, bg_coord, bg_mass;
        R = R_f, tol = tol, max_iters = max_iters, verbose = verbose)
    sol.converged || @warn(
        "build_scvt_mesh: Lloyd loop did not converge in $(sol.iterations) iterations " *
        "(displacement $(sol.displacement) ≥ tol $tol); emitting the last iterate.")
    gen = sol.generators                       # (3, Nc), on the sphere of radius R

    # --- D3 topology LEAF, ONCE at convergence ----------------------------------
    conn = scvt_voronoi_connectivity(gen; R = R_f)

    # --- dual GEOMETRY (shared FAQ stack with the DUO Voronoi dual) -------------
    # Canonical edge numbering = `primal_topology_faq(faces, Nc).edges`, the SAME
    # numbering the D3 leaf's `edges_on_cell` references, so every per-edge array
    # below is index-aligned with the leaf connectivity.
    prim = primal_topology_faq(conn.faces, Nc)
    edge_cells = _scvt_edge_cells(conn.faces, prim.edges)
    # The dual-cell area fan is rotation-invariant, so the leaf's CCW ring
    # (`vertices_on_cell`, valence-trimmed) serves as `sorted_vertex_faces`.
    sorted_vertex_faces = Vector{Vector{Int}}(undef, Nc)
    @inbounds for c in 1:Nc
        k = conn.n_edges_on_cell[c]
        sorted_vertex_faces[c] = Int[conn.vertices_on_cell[i, c] for i in 1:k]
    end
    geo = duo_dual_geometry_faq(Float64, gen, conn.faces, conn.circumcenters,
        prim.edges, edge_cells, sorted_vertex_faces, R_f)

    # MPAS edge ↔ Delaunay edge: the two cells are the edge's generator pair.
    Ne = size(prim.edges, 2)
    cells_on_edge = Matrix{Int}(undef, 2, Ne)
    @inbounds for e in 1:Ne
        cells_on_edge[1, e] = prim.edges[1, e]
        cells_on_edge[2, e] = prim.edges[2, e]
    end

    return mpas_mesh_data(;
        lon_cell = geo.lon_cell,
        lat_cell = geo.lat_cell,
        x_cell = gen[1, :], y_cell = gen[2, :], z_cell = gen[3, :],
        area_cell = geo.area_cell,
        n_edges_on_cell = conn.n_edges_on_cell,
        cells_on_cell = conn.cells_on_cell,
        edges_on_cell = conn.edges_on_cell,
        lon_edge = geo.lon_edge,
        lat_edge = geo.lat_edge,
        cells_on_edge = cells_on_edge,
        dc_edge = geo.dc_edge,
        dv_edge = geo.dv_edge,
        max_edges = conn.max_edges,
        n_vertices = conn.n_triangles,
        R = R_f,
    )
end
