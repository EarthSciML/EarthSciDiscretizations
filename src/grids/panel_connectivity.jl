"""
Panel connectivity for the cubed-sphere grid.

Standard cubed-sphere panel numbering (1-6) with connectivity table
defining neighbor panels and index transformations for each edge.

Panel 1: front, Panel 2: right, Panel 3: top (north pole),
Panel 4: back, Panel 5: left, Panel 6: bottom (south pole).
"""

@enum EdgeDirection West East South North

struct PanelNeighbor
    neighbor_panel::Int
    neighbor_edge::EdgeDirection
    reverse_index::Bool
end

const PANEL_CONNECTIVITY = Dict{Int, Dict{EdgeDirection, PanelNeighbor}}(
    1 => Dict(
        West => PanelNeighbor(5, East, false), East => PanelNeighbor(2, West, false),
        South => PanelNeighbor(6, North, false), North => PanelNeighbor(3, South, false)
    ),
    2 => Dict(
        West => PanelNeighbor(1, East, false), East => PanelNeighbor(4, West, false),
        South => PanelNeighbor(6, East, true), North => PanelNeighbor(3, East, false)
    ),
    3 => Dict(
        West => PanelNeighbor(5, North, true), East => PanelNeighbor(2, North, false),
        South => PanelNeighbor(1, North, false), North => PanelNeighbor(4, North, true)
    ),
    4 => Dict(
        West => PanelNeighbor(2, East, false), East => PanelNeighbor(5, West, false),
        South => PanelNeighbor(6, South, true), North => PanelNeighbor(3, North, true)
    ),
    5 => Dict(
        West => PanelNeighbor(4, East, false), East => PanelNeighbor(1, West, false),
        South => PanelNeighbor(6, West, false), North => PanelNeighbor(3, West, true)
    ),
    6 => Dict(
        West => PanelNeighbor(5, South, false), East => PanelNeighbor(2, South, true),
        South => PanelNeighbor(4, South, true), North => PanelNeighbor(1, South, false)
    ),
)

function transform_indices(neighbor::PanelNeighbor, i, j, Nc)
    nb_edge = neighbor.neighbor_edge
    rev = neighbor.reverse_index
    j_along = j
    i_perp = i
    if rev
        j_along = Nc + 1 - j_along
    end
    if nb_edge == West
        return (Nc + 1 - i_perp, j_along)
    elseif nb_edge == East
        return (i_perp, j_along)
    elseif nb_edge == South
        return (j_along, Nc + 1 - i_perp)
    elseif nb_edge == North
        return (j_along, i_perp)
    end
end

function verify_connectivity(Nc)
    for p in 1:6
        for dir in (West, East, South, North)
            nb = PANEL_CONNECTIVITY[p][dir]
            nb2 = PANEL_CONNECTIVITY[nb.neighbor_panel][nb.neighbor_edge]
            nb2.neighbor_panel != p && error("Connectivity inconsistency at panel $p $dir")
        end
    end
    return true
end
