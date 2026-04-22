"""
`EarthSciDiscretizations.grids` — per GRIDS_API.md §2.3, one function per
grid family. Family generators are re-exported into this submodule; they
continue to live at the top level for backwards compatibility with
existing ESD consumers.
"""
module grids

import ..EarthSciDiscretizations:
    build_duo_grid,
    DuoGrid,
    DuoLoader

"""
    grids.duo(; loader, R=6.371e6, dtype=Float64, ghosts=0) -> DuoGrid

Loader-backed icosahedral triangular grid (Heikes et al. 2023). Only
`builtin://icosahedral/<level>` loader paths are honored today; `.duo`
mesh files land with the ESS file-format spec.
"""
duo(; kwargs...) = build_duo_grid(; kwargs...)

export duo

end # module grids
