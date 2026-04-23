"""
DUO icosahedral triangular grid family (Heikes et al. 2023).

Loader-backed per GRIDS_API.md §10. `.esm` carries a small declarative
config (family, level, reader ref); geometry is derived at runtime by
the accessors. Cross-binding conformance is compared at pinned query
points, not serialized bytes.

Subdivision level `r` yields 20·4^r triangular cells, 10·4^r + 2 vertices,
30·4^r edges. All vertices lie on a sphere of radius R.
"""

struct DuoLoader
    path::String
    reader::String
    check::String
end

DuoLoader(;
    path::AbstractString, reader::AbstractString = "auto",
    check::AbstractString = "strict"
) =
    DuoLoader(String(path), String(reader), String(check))

"""
    DuoGrid{T} <: AbstractGrid

Icosahedral triangular mesh, fully materialized.

Fields:
- `level`: subdivision level r (r=0 → 20 cells, bare icosahedron)
- `R`, `dtype`, `ghosts`: per GRIDS_API §2.7
- `vertices`: (3, Nv) cartesian unit-sphere coords scaled by R
- `faces`: (3, Nc) vertex indices per triangular cell (1-based)
- `lon`, `lat`: (Nc,) cell-center geographic coords (centroid of 3 vertex
  unit vectors, renormalized)
- `cell_cart`: (3, Nc) cell-center cartesian coords on the sphere
- `area`: (Nc,) spherical-triangle area per cell (R² scaled)
- `edges`: (2, Ne) sorted vertex-pair per edge
- `cell_neighbors`: (3, Nc) neighboring cell indices across each of the 3
  edges of cell c; 0 for boundary (bare icosahedron has none, so no zeros
  in a closed mesh)
- `vertex_faces`: ragged; `vertex_faces[v]` is the sorted list of face
  indices incident on vertex v
- `provenance`: §6.4 fingerprint
- `loader`: DuoLoader the grid was generated from
"""
struct DuoGrid{T} <: AbstractGrid
    level::Int
    R::T
    dtype::String
    ghosts::Int
    vertices::Matrix{T}
    faces::Matrix{Int}
    lon::Vector{T}
    lat::Vector{T}
    cell_cart::Matrix{T}
    area::Vector{T}
    edges::Matrix{Int}
    cell_neighbors::Matrix{Int}
    vertex_faces::Vector{Vector{Int}}
    provenance::Dict{String, Any}
    loader::DuoLoader
end

const _DUO_FAMILY_VERSION = "1.0.0"

# --- Base icosahedron ---------------------------------------------------

"""
Return the 12 canonical icosahedron vertices on the unit sphere, as a
(3, 12) matrix of element type `T`. Ordering is fixed and deterministic
so subdivision level-0 output is byte-stable across runs and platforms.
"""
function _icosahedron_vertices(::Type{T}) where {T}
    φ = (T(1) + sqrt(T(5))) / T(2)
    raw = T[
        0  1  φ;
        0 -1  φ;
        0  1 -φ;
        0 -1 -φ;
        1  φ  0;
        -1  φ  0;
        1 -φ  0;
        -1 -φ  0;
        φ  0  1;
        φ  0 -1;
        -φ  0  1;
        -φ  0 -1;
    ]
    V = Matrix{T}(undef, 3, 12)
    for i in 1:12
        v = @view raw[i, :]
        n = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        V[1, i] = v[1] / n
        V[2, i] = v[2] / n
        V[3, i] = v[3] / n
    end
    return V
end

"""
Return the 20 triangular faces of the base icosahedron as (3, 20) Int
matrix of 1-based vertex indices into `_icosahedron_vertices`. Winding
is outward (right-hand rule → outward normal for positive spherical area).
"""
function _icosahedron_faces()
    F = [
        1  9  2;
        1  2 11;
        1 11  6;
        1  6  5;
        1  5  9;
        9  5 10;
        9 10  7;
        9  7  2;
        2  7  8;
        2  8 11;
        11  8 12;
        11 12  6;
        6 12  3;
        6  3  5;
        5  3 10;
        10  3  4;
        10  4  7;
        7  4  8;
        8  4 12;
        12  4  3;
    ]
    return Matrix{Int}(F')
end

# --- Recursive subdivision ---------------------------------------------

# Edge-midpoint cache keyed on sorted vertex-pair → new vertex index.
function _midpoint!(
        V::Vector{NTuple{3, T}}, cache::Dict{Tuple{Int, Int}, Int},
        a::Int, b::Int
    ) where {T}
    key = a < b ? (a, b) : (b, a)
    idx = get(cache, key, 0)
    if idx != 0
        return idx
    end
    va = V[a]; vb = V[b]
    mx = va[1] + vb[1]; my = va[2] + vb[2]; mz = va[3] + vb[3]
    n = sqrt(mx^2 + my^2 + mz^2)
    push!(V, (mx / n, my / n, mz / n))
    idx = length(V)
    cache[key] = idx
    return idx
end

"""
Subdivide the base icosahedron `level` times on the unit sphere. Returns
(vertices 3×Nv, faces 3×Nc) at level `level`.
"""
function _subdivide_icosahedron(::Type{T}, level::Int) where {T}
    V0 = _icosahedron_vertices(T)
    F0 = _icosahedron_faces()
    verts = NTuple{3, T}[(V0[1, i], V0[2, i], V0[3, i]) for i in 1:size(V0, 2)]
    faces = Tuple{Int, Int, Int}[(F0[1, c], F0[2, c], F0[3, c]) for c in 1:size(F0, 2)]

    for _ in 1:level
        cache = Dict{Tuple{Int, Int}, Int}()
        new_faces = Tuple{Int, Int, Int}[]
        sizehint!(new_faces, 4 * length(faces))
        for (a, b, c) in faces
            ab = _midpoint!(verts, cache, a, b)
            bc = _midpoint!(verts, cache, b, c)
            ca = _midpoint!(verts, cache, c, a)
            push!(new_faces, (a, ab, ca))
            push!(new_faces, (b, bc, ab))
            push!(new_faces, (c, ca, bc))
            push!(new_faces, (ab, bc, ca))
        end
        faces = new_faces
    end

    Nv = length(verts); Nc = length(faces)
    V = Matrix{T}(undef, 3, Nv)
    for i in 1:Nv
        v = verts[i]
        V[1, i] = v[1]; V[2, i] = v[2]; V[3, i] = v[3]
    end
    F = Matrix{Int}(undef, 3, Nc)
    for c in 1:Nc
        f = faces[c]
        F[1, c] = f[1]; F[2, c] = f[2]; F[3, c] = f[3]
    end
    return V, F
end

# --- Geometry helpers --------------------------------------------------

_cart_to_lonlat(x::T, y::T, z::T) where {T} = (atan(y, x), asin(clamp(z, -one(T), one(T))))

"""
Spherical excess area of the unit-sphere triangle with vertices a,b,c
(given as 3-tuples). Uses L'Huilier's theorem for numerical stability
on near-degenerate triangles.
"""
function _spherical_triangle_area(a, b, c)
    # Side lengths (great-circle arcs, on unit sphere).
    da = acos(clamp(b[1] * c[1] + b[2] * c[2] + b[3] * c[3], -1.0, 1.0))
    db = acos(clamp(c[1] * a[1] + c[2] * a[2] + c[3] * a[3], -1.0, 1.0))
    dc = acos(clamp(a[1] * b[1] + a[2] * b[2] + a[3] * b[3], -1.0, 1.0))
    s = 0.5 * (da + db + dc)
    # L'Huilier: tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
    t = tan(s / 2) * tan((s - da) / 2) * tan((s - db) / 2) * tan((s - dc) / 2)
    t = max(t, 0.0)
    return 4 * atan(sqrt(t))
end

# --- Connectivity -------------------------------------------------------

"""
Return (edges, cell_neighbors) from a faces matrix (3, Nc).
- `edges`: (2, Ne) sorted vertex pairs (a < b).
- `cell_neighbors`: (3, Nc) where `cell_neighbors[k, c]` is the cell
  sharing the edge opposite vertex `faces[k, c]` in cell c. That edge
  is (faces[k%3+1, c], faces[(k+1)%3+1, c]). 0 if no neighbor.
"""
function _build_connectivity(faces::Matrix{Int})
    Nc = size(faces, 2)
    # edge_key (min,max) → Vector of (cell, local_edge_k)
    edge_map = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for c in 1:Nc
        v1 = faces[1, c]; v2 = faces[2, c]; v3 = faces[3, c]
        # Edge opposite vertex k=1 is (v2, v3); k=2 → (v3, v1); k=3 → (v1, v2)
        for (k, (a, b)) in enumerate(((v2, v3), (v3, v1), (v1, v2)))
            key = a < b ? (a, b) : (b, a)
            list = get!(edge_map, key, Tuple{Int, Int}[])
            push!(list, (c, k))
        end
    end
    neighbors = zeros(Int, 3, Nc)
    edges = Matrix{Int}(undef, 2, length(edge_map))
    e = 0
    for (key, list) in edge_map
        e += 1
        edges[1, e] = key[1]; edges[2, e] = key[2]
        if length(list) == 2
            (c1, k1), (c2, k2) = list[1], list[2]
            neighbors[k1, c1] = c2
            neighbors[k2, c2] = c1
        elseif length(list) > 2
            throw(AssertionError("duo: non-manifold edge $(key) shared by $(length(list)) cells"))
        end
        # length==1 → boundary; closed icosahedral mesh has none.
    end
    return edges, neighbors
end

function _vertex_faces(faces::Matrix{Int}, Nv::Int)
    vf = [Int[] for _ in 1:Nv]
    for c in 1:size(faces, 2), k in 1:3
        push!(vf[faces[k, c]], c)
    end
    for v in 1:Nv
        sort!(vf[v])
    end
    return vf
end

# --- Loader --------------------------------------------------------------

const _BUILTIN_READERS = (
    "builtin_icosahedral",
    "duo_mesh",  # file-backed; honored via builtin path scheme below
    "auto",
)

"""
Parse `builtin://icosahedral/<level>` loader paths. Returns the
subdivision level as an Int, or `nothing` if the path is not a builtin
icosahedral spec.
"""
function _parse_builtin_level(path::AbstractString)
    prefix = "builtin://icosahedral/"
    if !startswith(path, prefix)
        return nothing
    end
    tail = path[(length(prefix) + 1):end]
    lvl = tryparse(Int, tail)
    lvl === nothing && throw(ArgumentError("duo: cannot parse level from loader path $path"))
    lvl < 0 && throw(DomainError(lvl, "duo: subdivision level must be ≥ 0"))
    return lvl
end

"""
Resolve a loader to a subdivision level. Falls back to ArgumentError
when the reader/path combination is not yet implemented — a future
`.duo` mesh reader will land with the ESS file-format spec.
"""
function _resolve_loader_level(loader::DuoLoader)
    lvl = _parse_builtin_level(loader.path)
    if lvl !== nothing
        return lvl
    end
    if loader.reader == "duo_mesh" || loader.reader == "auto"
        throw(
            ArgumentError(
                "duo: .duo mesh-file reader not yet implemented — pending " *
                    "EarthSciSerialization file-format spec. Use " *
                    "builtin://icosahedral/<level> in the meantime."
            )
        )
    end
    throw(ArgumentError("duo: unrecognized loader path $(loader.path) with reader=$(loader.reader)"))
end

# --- Generator ----------------------------------------------------------

"""
    build_duo_grid(; loader, R=6.371e6, dtype=Float64, ghosts=0)

Primary constructor for the DUO family. `loader` is a `DuoLoader` or
NamedTuple/Dict with keys `(path, reader, check)` per GRIDS_API §10.

Only `builtin://icosahedral/<level>` loader paths are accepted today;
`.duo` file readers land with the ESS file-format spec (future bead).
"""
function build_duo_grid(;
        loader,
        R::Real = 6.371e6,
        dtype = Float64,
        ghosts::Int = 0
    )
    ghosts >= 0 || throw(DomainError(ghosts, "duo: ghosts must be ≥ 0"))
    if !(dtype === Float64 || dtype === Float32)
        throw(ArgumentError("duo: dtype must be Float64 or Float32, got $dtype"))
    end

    ldr = _coerce_loader(loader)
    level = _resolve_loader_level(ldr)

    T = dtype
    R_T = T(R)

    V, F = _subdivide_icosahedron(T, level)
    Nv = size(V, 2)
    Nc = size(F, 2)

    # Cell centers: normalized centroid of the three unit vertex vectors.
    cell_cart = Matrix{T}(undef, 3, Nc)
    lon = Vector{T}(undef, Nc)
    lat = Vector{T}(undef, Nc)
    area = Vector{T}(undef, Nc)
    R2 = R_T * R_T
    for c in 1:Nc
        a_i = F[1, c]; b_i = F[2, c]; c_i = F[3, c]
        ax = V[1, a_i]; ay = V[2, a_i]; az = V[3, a_i]
        bx = V[1, b_i]; by = V[2, b_i]; bz = V[3, b_i]
        cx = V[1, c_i]; cy = V[2, c_i]; cz = V[3, c_i]
        mx = ax + bx + cx; my = ay + by + cy; mz = az + bz + cz
        n = sqrt(mx^2 + my^2 + mz^2)
        ux = mx / n; uy = my / n; uz = mz / n
        cell_cart[1, c] = R_T * ux
        cell_cart[2, c] = R_T * uy
        cell_cart[3, c] = R_T * uz
        lon[c], lat[c] = _cart_to_lonlat(ux, uy, uz)
        # Area via L'Huilier on unit sphere, then scale by R².
        a_unit = _spherical_triangle_area(
            (Float64(ax), Float64(ay), Float64(az)),
            (Float64(bx), Float64(by), Float64(bz)),
            (Float64(cx), Float64(cy), Float64(cz)),
        )
        area[c] = T(a_unit * Float64(R2))
    end

    # Scale vertex array to radius R for downstream consumers.
    V_scaled = Matrix{T}(undef, 3, Nv)
    @inbounds for i in 1:Nv
        V_scaled[1, i] = R_T * V[1, i]
        V_scaled[2, i] = R_T * V[2, i]
        V_scaled[3, i] = R_T * V[3, i]
    end

    edges, cell_neighbors = _build_connectivity(F)
    vf = _vertex_faces(F, Nv)

    provenance = Dict{String, Any}(
        "binding" => "julia",
        "family" => "duo",
        "version" => _DUO_FAMILY_VERSION,
        "level" => level,
        "reader" => ldr.reader,
        "path" => ldr.path,
        "check" => ldr.check,
        "dtype" => dtype === Float64 ? "float64" : "float32",
    )

    dtype_str = dtype === Float64 ? "float64" : "float32"
    return DuoGrid{T}(
        level, R_T, dtype_str, ghosts,
        V_scaled, F, lon, lat, cell_cart, area,
        edges, cell_neighbors, vf, provenance, ldr
    )
end

# Accept NamedTuple / Dict / DuoLoader as loader input.
_coerce_loader(ldr::DuoLoader) = ldr
function _coerce_loader(ldr::NamedTuple)
    path = get(ldr, :path, nothing)
    path === nothing && throw(ArgumentError("duo: loader.path is required"))
    return DuoLoader(;
        path = path,
        reader = String(get(ldr, :reader, "auto")),
        check = String(get(ldr, :check, "strict")),
    )
end
function _coerce_loader(ldr::AbstractDict)
    path = get(ldr, "path", get(ldr, :path, nothing))
    path === nothing && throw(ArgumentError("duo: loader.path is required"))
    return DuoLoader(;
        path = String(path),
        reader = String(get(ldr, "reader", get(ldr, :reader, "auto"))),
        check = String(get(ldr, "check", get(ldr, :check, "strict"))),
    )
end

# --- Accessors (per mayor's correction: runtime-derived) ---------------

"""
    cell_centers(grid::DuoGrid) -> NamedTuple{(:lon,:lat)}
    cell_centers(grid::DuoGrid, c::Integer) -> (lon, lat)
"""
cell_centers(g::DuoGrid) = (lon = g.lon, lat = g.lat)
cell_centers(g::DuoGrid, c::Integer) = (g.lon[c], g.lat[c])

"""
    neighbors(grid::DuoGrid, c::Integer) -> NTuple{3,Int}

Returns the three cell indices sharing an edge with cell c (0 = boundary,
not reachable on a closed icosahedral mesh).
"""
neighbors(g::DuoGrid, c::Integer) =
    (g.cell_neighbors[1, c], g.cell_neighbors[2, c], g.cell_neighbors[3, c])

"""
    metric_eval(grid::DuoGrid, name::Symbol, c::Integer)

Accessor for per-cell scalar metrics. Supported names:
- `:area` — spherical-triangle area (m² if R is in m)
- `:lon`, `:lat` — cell-center geographic coords (rad)
- `:x`, `:y`, `:z` — cell-center cartesian coords (same units as R)
"""
function metric_eval(g::DuoGrid, name::Symbol, c::Integer)
    if name === :area
        return g.area[c]
    elseif name === :lon
        return g.lon[c]
    elseif name === :lat
        return g.lat[c]
    elseif name === :x
        return g.cell_cart[1, c]
    elseif name === :y
        return g.cell_cart[2, c]
    elseif name === :z
        return g.cell_cart[3, c]
    else
        throw(ArgumentError("duo: metric_eval: unknown metric $name"))
    end
end

metric_eval(g::DuoGrid, name::AbstractString, c::Integer) =
    metric_eval(g, Symbol(name), c)

n_cells(g::DuoGrid) = size(g.faces, 2)
n_vertices(g::DuoGrid) = size(g.vertices, 2)
n_edges(g::DuoGrid) = size(g.edges, 2)
total_area(g::DuoGrid) = sum(g.area)

# --- .esm lowering (small declarative config per mayor's correction) ---

"""
    to_esm(grid::DuoGrid) -> Dict{String,Any}

Returns a §6-schema-valid declarative config: family + level +
loader ref + tolerances + provenance. No inline geometry arrays.
"""
function to_esm(g::DuoGrid)
    return Dict{String, Any}(
        "family" => "duo",
        "topology" => "unstructured",
        "dtype" => g.dtype,
        "ghosts" => g.ghosts,
        "n_cells" => n_cells(g),
        "n_vertices" => n_vertices(g),
        "n_edges" => n_edges(g),
        "options" => Dict{String, Any}(
            "R" => g.R,
            "level" => g.level,
            "loader" => Dict{String, Any}(
                "path" => g.loader.path,
                "reader" => g.loader.reader,
                "check" => g.loader.check,
            ),
        ),
        "provenance" => g.provenance,
        "schema_version" => _DUO_FAMILY_VERSION,
    )
end

family(::DuoGrid) = "duo"
