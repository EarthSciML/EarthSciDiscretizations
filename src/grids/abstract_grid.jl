"""
Abstract grid types for EarthSciDiscretizations.

The Grid trait is owned by EarthSciSerialization (RFC `docs/rfcs/grid-trait.md`,
ESS `src/abstract_grid.jl`). ESD re-exports the trait root + tier subtypes
and registers concrete grid families against them. Tier-C / Tier-M trait
methods (and the `GridTraitError` exception) come from ESS as well; concrete
families override them on the bulk-array contract.

The legacy supertype `AbstractCubedSphereGrid` is kept here so existing ESD
code that dispatches on it (panel connectivity, ghost-cell utilities) keeps
working — it now sits beneath the ESS curvilinear root.
"""

using EarthSciSerialization:
    AbstractGrid,
    AbstractCurvilinearGrid,
    AbstractStaggeredGrid,
    AbstractVerticalGrid,
    AbstractUnstructuredGrid,
    GridTraitError
import EarthSciSerialization:
    cell_centers,
    cell_volume,
    cell_widths,
    neighbor_indices,
    boundary_mask,
    n_cells,
    n_dims,
    axis_names,
    metric_g,
    metric_ginv,
    metric_jacobian,
    metric_dgij_dxk,
    coord_jacobian,
    coord_jacobian_second

abstract type AbstractCubedSphereGrid <: AbstractCurvilinearGrid end

# ---------------------------------------------------------------------------
# Per-grid lazy materialization cache (RFC §1: don't recompute bulk arrays
# on every call). Indexed by object identity so the cache survives `==`
# overrides without hashing the whole grid.
# ---------------------------------------------------------------------------

const _GRID_TRAIT_CACHE = Base.IdDict{Any, Dict{Any, Any}}()

"""
    _grid_memo!(builder, grid, key)

Materialize and cache a bulk-array trait value for `grid` keyed by `key`.
`builder()` is called the first time the key is seen; subsequent calls
return the cached result.
"""
function _grid_memo!(builder, grid, key)
    cache = get!(_GRID_TRAIT_CACHE, grid) do
        Dict{Any, Any}()
    end
    return get!(builder, cache, key)
end
