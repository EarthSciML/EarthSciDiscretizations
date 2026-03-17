"""
Initial condition projection onto the cubed-sphere grid.
"""

function project_initial_condition(f, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    Nc = grid.Nc; ni, nj = grid_size(loc, Nc)
    result = zeros(6, ni, nj)
    for p in 1:6, i in 1:ni, j in 1:nj
        if loc == CellCenter
            result[p, i, j] = f(grid.lon[p, i, j], grid.lat[p, i, j])
        else
            ic = clamp(i, 1, Nc); jc = clamp(j, 1, Nc)
            result[p, i, j] = f(grid.lon[p, ic, jc], grid.lat[p, ic, jc])
        end
    end
    return result
end
