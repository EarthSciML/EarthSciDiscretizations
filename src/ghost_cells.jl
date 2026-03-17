"""
Ghost cell management for cubed-sphere panel boundaries.
"""

function fill_ghost_cells!(q_ext, q, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    Nc = grid.Nc; Ng = grid.Ng
    ni, nj = grid_size(loc, Nc)
    for p in 1:6, i in 1:ni, j in 1:nj
        q_ext[p, i + Ng, j + Ng] = q[p, i, j]
    end
    for p in 1:6, dir in (West, East, South, North)
        nb = PANEL_CONNECTIVITY[p][dir]
        if dir == West
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, j, ni, nj, dir)
                q_ext[p, g, j + Ng] = q[nb.neighbor_panel, clamp(i_nb, 1, ni), clamp(j_nb, 1, nj)]
            end
        elseif dir == East
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, g, j, ni, nj, dir)
                q_ext[p, ni + Ng + g, j + Ng] = q[nb.neighbor_panel, clamp(i_nb, 1, ni), clamp(j_nb, 1, nj)]
            end
        elseif dir == South
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, i, ni, nj, dir)
                q_ext[p, i + Ng, g] = q[nb.neighbor_panel, clamp(i_nb, 1, ni), clamp(j_nb, 1, nj)]
            end
        elseif dir == North
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, g, i, ni, nj, dir)
                q_ext[p, i + Ng, nj + Ng + g] = q[nb.neighbor_panel, clamp(i_nb, 1, ni), clamp(j_nb, 1, nj)]
            end
        end
    end
    return q_ext
end

function extend_with_ghosts(q, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    Nc = grid.Nc; Ng = grid.Ng
    ni, nj = grid_size(loc, Nc)
    q_ext = zeros(eltype(q), 6, ni + 2Ng, nj + 2Ng)
    fill_ghost_cells!(q_ext, q, grid, loc)
    return q_ext
end

function transform_ghost_index(neighbor::PanelNeighbor, i_perp, j_along, ni, nj, src_dir)
    nb_edge = neighbor.neighbor_edge
    j_max = (src_dir == West || src_dir == East) ? nj : ni
    j_mapped = neighbor.reverse_index ? (j_max + 1 - j_along) : j_along
    nb_edge == West && return (i_perp, j_mapped)
    nb_edge == East && return (ni + 1 - i_perp, j_mapped)
    nb_edge == South && return (j_mapped, i_perp)
    return (j_mapped, nj + 1 - i_perp)  # North
end

function ghost_fill_indices(grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    Nc = grid.Nc; Ng = grid.Ng
    ni, nj = grid_size(loc, Nc)
    ne_i = ni + 2Ng; ne_j = nj + 2Ng
    src_panel = zeros(Int, 6, ne_i, ne_j)
    src_i = zeros(Int, 6, ne_i, ne_j); src_j = zeros(Int, 6, ne_i, ne_j)
    for p in 1:6, i in 1:ni, j in 1:nj
        src_panel[p, i + Ng, j + Ng] = p
        src_i[p, i + Ng, j + Ng] = i; src_j[p, i + Ng, j + Ng] = j
    end
    for p in 1:6, dir in (West, East, South, North)
        nb = PANEL_CONNECTIVITY[p][dir]
        if dir == West
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, j, ni, nj, dir)
                src_panel[p, g, j + Ng] = nb.neighbor_panel
                src_i[p, g, j + Ng] = clamp(i_nb, 1, ni); src_j[p, g, j + Ng] = clamp(j_nb, 1, nj)
            end
        elseif dir == East
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, g, j, ni, nj, dir)
                src_panel[p, ni + Ng + g, j + Ng] = nb.neighbor_panel
                src_i[p, ni + Ng + g, j + Ng] = clamp(i_nb, 1, ni); src_j[p, ni + Ng + g, j + Ng] = clamp(j_nb, 1, nj)
            end
        elseif dir == South
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, i, ni, nj, dir)
                src_panel[p, i + Ng, g] = nb.neighbor_panel
                src_i[p, i + Ng, g] = clamp(i_nb, 1, ni); src_j[p, i + Ng, g] = clamp(j_nb, 1, nj)
            end
        elseif dir == North
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, g, i, ni, nj, dir)
                src_panel[p, i + Ng, nj + Ng + g] = nb.neighbor_panel
                src_i[p, i + Ng, nj + Ng + g] = clamp(i_nb, 1, ni); src_j[p, i + Ng, nj + Ng + g] = clamp(j_nb, 1, nj)
            end
        end
    end
    return (src_panel, src_i, src_j)
end

function ghost_fill_arrayop(u_interior, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    ni, nj = grid_size(loc, grid.Nc)
    idx = get_idx_vars(3); p, i, j = idx[1], idx[2], idx[3]
    u_c = const_wrap(unwrap(u_interior))
    expr = wrap(u_c[p, i, j])
    ranges = Dict(p => 1:1:6, i => 1:1:ni, j => 1:1:nj)
    return make_arrayop(idx, unwrap(expr), ranges)
end
