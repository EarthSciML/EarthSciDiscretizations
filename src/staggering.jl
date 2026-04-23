"""
C-grid staggering definitions.
"""

@enum VarLocation CellCenter UEdge VEdge Corner

function grid_size(loc::VarLocation, Nc::Int)
    return loc == CellCenter ? (Nc, Nc) :
        loc == UEdge ? (Nc + 1, Nc) :
        loc == VEdge ? (Nc, Nc + 1) :
        (Nc + 1, Nc + 1)
end

function full_array_size(loc::VarLocation, Nc::Int)
    ni, nj = grid_size(loc, Nc)
    return (6, ni, nj)
end

function ghost_array_size(loc::VarLocation, Nc::Int, Ng::Int)
    ni, nj = grid_size(loc, Nc)
    return (6, ni + 2Ng, nj + 2Ng)
end
