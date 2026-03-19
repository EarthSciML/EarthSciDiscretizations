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

    # Fill corner ghost cells by resolving through two successive edge crossings.
    # At cube vertices (where 3 panels meet), we pick values from the nearest
    # edge ghost cell, which was already filled above.
    for p in 1:6
        # SW corner: extended indices (1:Ng, 1:Ng)
        for gi in 1:Ng, gj in 1:Ng
            q_ext[p, gi, gj] = q_ext[p, gi, Ng + 1]
        end
        # SE corner: extended indices (ni+Ng+1:ni+2Ng, 1:Ng)
        for gi in 1:Ng, gj in 1:Ng
            q_ext[p, ni + Ng + gi, gj] = q_ext[p, ni + Ng + gi, Ng + 1]
        end
        # NW corner: extended indices (1:Ng, nj+Ng+1:nj+2Ng)
        for gi in 1:Ng, gj in 1:Ng
            q_ext[p, gi, nj + Ng + gj] = q_ext[p, gi, nj + Ng]
        end
        # NE corner: extended indices (ni+Ng+1:ni+2Ng, nj+Ng+1:nj+2Ng)
        for gi in 1:Ng, gj in 1:Ng
            q_ext[p, ni + Ng + gi, nj + Ng + gj] = q_ext[p, ni + Ng + gi, nj + Ng]
        end
    end

    return q_ext
end

"""
    fill_ghost_cells_vector!(uξ_ext, uη_ext, uξ, uη, grid, loc)

Fill ghost cells for a contravariant vector field (uξ, uη), applying the
proper rotation when crossing panel boundaries so that components are expressed
in the local panel's (ξ, η) basis.

The rotation matrix for each edge is precomputed in `grid.rotation_matrices`.
"""
function fill_ghost_cells_vector!(uξ_ext, uη_ext, uξ, uη, grid::CubedSphereGrid,
                                   loc::VarLocation = CellCenter)
    Nc = grid.Nc; Ng = grid.Ng
    ni, nj = grid_size(loc, Nc)

    # Copy interior
    for p in 1:6, i in 1:ni, j in 1:nj
        uξ_ext[p, i + Ng, j + Ng] = uξ[p, i, j]
        uη_ext[p, i + Ng, j + Ng] = uη[p, i, j]
    end

    # Fill edge ghost cells with rotation
    for p in 1:6, dir in (West, East, South, North)
        nb = PANEL_CONNECTIVITY[p][dir]
        M11, M12, M21, M22 = grid.rotation_matrices[(p, dir)]

        if dir == West
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, j, ni, nj, dir)
                i_nb = clamp(i_nb, 1, ni); j_nb = clamp(j_nb, 1, nj)
                vξ = uξ[nb.neighbor_panel, i_nb, j_nb]
                vη = uη[nb.neighbor_panel, i_nb, j_nb]
                uξ_ext[p, g, j + Ng] = M11 * vξ + M12 * vη
                uη_ext[p, g, j + Ng] = M21 * vξ + M22 * vη
            end
        elseif dir == East
            for g in 1:Ng, j in 1:nj
                i_nb, j_nb = transform_ghost_index(nb, g, j, ni, nj, dir)
                i_nb = clamp(i_nb, 1, ni); j_nb = clamp(j_nb, 1, nj)
                vξ = uξ[nb.neighbor_panel, i_nb, j_nb]
                vη = uη[nb.neighbor_panel, i_nb, j_nb]
                uξ_ext[p, ni + Ng + g, j + Ng] = M11 * vξ + M12 * vη
                uη_ext[p, ni + Ng + g, j + Ng] = M21 * vξ + M22 * vη
            end
        elseif dir == South
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, Ng + 1 - g, i, ni, nj, dir)
                i_nb = clamp(i_nb, 1, ni); j_nb = clamp(j_nb, 1, nj)
                vξ = uξ[nb.neighbor_panel, i_nb, j_nb]
                vη = uη[nb.neighbor_panel, i_nb, j_nb]
                uξ_ext[p, i + Ng, g] = M11 * vξ + M12 * vη
                uη_ext[p, i + Ng, g] = M21 * vξ + M22 * vη
            end
        else  # North
            for i in 1:ni, g in 1:Ng
                i_nb, j_nb = transform_ghost_index(nb, g, i, ni, nj, dir)
                i_nb = clamp(i_nb, 1, ni); j_nb = clamp(j_nb, 1, nj)
                vξ = uξ[nb.neighbor_panel, i_nb, j_nb]
                vη = uη[nb.neighbor_panel, i_nb, j_nb]
                uξ_ext[p, i + Ng, nj + Ng + g] = M11 * vξ + M12 * vη
                uη_ext[p, i + Ng, nj + Ng + g] = M21 * vξ + M22 * vη
            end
        end
    end

    # Corner ghost cells: copy from nearest edge ghost (already rotated)
    for p in 1:6
        for gi in 1:Ng, gj in 1:Ng
            uξ_ext[p, gi, gj] = uξ_ext[p, gi, Ng + 1]
            uη_ext[p, gi, gj] = uη_ext[p, gi, Ng + 1]
            uξ_ext[p, ni + Ng + gi, gj] = uξ_ext[p, ni + Ng + gi, Ng + 1]
            uη_ext[p, ni + Ng + gi, gj] = uη_ext[p, ni + Ng + gi, Ng + 1]
            uξ_ext[p, gi, nj + Ng + gj] = uξ_ext[p, gi, nj + Ng]
            uη_ext[p, gi, nj + Ng + gj] = uη_ext[p, gi, nj + Ng]
            uξ_ext[p, ni + Ng + gi, nj + Ng + gj] = uξ_ext[p, ni + Ng + gi, nj + Ng]
            uη_ext[p, ni + Ng + gi, nj + Ng + gj] = uη_ext[p, ni + Ng + gi, nj + Ng]
        end
    end

    return (uξ_ext, uη_ext)
end

function extend_with_ghosts_vector(uξ, uη, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    Nc = grid.Nc; Ng = grid.Ng
    ni, nj = grid_size(loc, Nc)
    T = promote_type(eltype(uξ), eltype(uη))
    uξ_ext = zeros(T, 6, ni + 2Ng, nj + 2Ng)
    uη_ext = zeros(T, 6, ni + 2Ng, nj + 2Ng)
    fill_ghost_cells_vector!(uξ_ext, uη_ext, uξ, uη, grid, loc)
    return (uξ_ext, uη_ext)
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
