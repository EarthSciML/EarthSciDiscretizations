"""
Discrete variable space for cubed-sphere grids.
"""

struct CubedSphereDiscreteSpace
    grid::CubedSphereGrid
    var_names::Vector{Symbol}
    locations::Dict{Symbol,VarLocation}
end

function CubedSphereDiscreteSpace(grid::CubedSphereGrid, var_specs::Vector{Pair{Symbol,VarLocation}})
    var_names = Symbol[]
    locations = Dict{Symbol,VarLocation}()
    for (name, loc) in var_specs
        push!(var_names, name)
        locations[name] = loc
    end
    CubedSphereDiscreteSpace(grid, var_names, locations)
end

function allocate_variable(space::CubedSphereDiscreteSpace, name::Symbol)
    loc = space.locations[name]
    dims = full_array_size(loc, space.grid.Nc)
    return zeros(dims...)
end
