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
    # At cube vertices (where 3 panels meet), we trace through two panel
    # boundaries to find the correct source cell on the diagonal neighbor.
    for p in 1:6
        # SW corner: extended indices (1:Ng, 1:Ng)
        for gi in 1:Ng, gj in 1:Ng
            i_virt = gi - Ng  # virtual interior index (negative)
            j_virt = gj - Ng
            q_ext[p, gi, gj] = _resolve_corner_value(q, p, i_virt, j_virt, ni, nj)
        end
        # SE corner: extended indices (ni+Ng+1:ni+2Ng, 1:Ng)
        for gi in 1:Ng, gj in 1:Ng
            i_virt = ni + gi  # past the east boundary
            j_virt = gj - Ng
            q_ext[p, ni + Ng + gi, gj] = _resolve_corner_value(q, p, i_virt, j_virt, ni, nj)
        end
        # NW corner: extended indices (1:Ng, nj+Ng+1:nj+2Ng)
        for gi in 1:Ng, gj in 1:Ng
            i_virt = gi - Ng
            j_virt = nj + gj  # past the north boundary
            q_ext[p, gi, nj + Ng + gj] = _resolve_corner_value(q, p, i_virt, j_virt, ni, nj)
        end
        # NE corner: extended indices (ni+Ng+1:ni+2Ng, nj+Ng+1:nj+2Ng)
        for gi in 1:Ng, gj in 1:Ng
            i_virt = ni + gi
            j_virt = nj + gj
            q_ext[p, ni + Ng + gi, nj + Ng + gj] = _resolve_corner_value(q, p, i_virt, j_virt, ni, nj)
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

"""
    _resolve_corner_value(q, p, i_virt, j_virt, ni, nj, max_depth=2)

Resolve a ghost cell position that may require crossing two panel boundaries
(corner ghost cells). Uses recursive edge crossings with a depth limit.

`i_virt` and `j_virt` are virtual interior indices that may be outside [1:ni]
or [1:nj]. Each recursive call crosses one panel boundary, mapping the
out-of-range index to a position on the neighbor panel.
"""
function _resolve_corner_value(q, p, i_virt, j_virt, ni, nj, max_depth=2)
    if 1 <= i_virt <= ni && 1 <= j_virt <= nj
        return q[p, i_virt, j_virt]
    end
    if max_depth <= 0
        return q[p, clamp(i_virt, 1, ni), clamp(j_virt, 1, nj)]
    end

    # Resolve one out-of-range dimension at a time via edge crossing
    if i_virt < 1
        nb = PANEL_CONNECTIVITY[p][West]
        depth = 1 - i_virt
        i_nb, j_nb = transform_ghost_index(nb, depth, j_virt, ni, nj, West)
        return _resolve_corner_value(q, nb.neighbor_panel, i_nb, j_nb, ni, nj, max_depth - 1)
    elseif i_virt > ni
        nb = PANEL_CONNECTIVITY[p][East]
        depth = i_virt - ni
        i_nb, j_nb = transform_ghost_index(nb, depth, j_virt, ni, nj, East)
        return _resolve_corner_value(q, nb.neighbor_panel, i_nb, j_nb, ni, nj, max_depth - 1)
    elseif j_virt < 1
        nb = PANEL_CONNECTIVITY[p][South]
        depth = 1 - j_virt
        i_nb, j_nb = transform_ghost_index(nb, depth, i_virt, ni, nj, South)
        return _resolve_corner_value(q, nb.neighbor_panel, i_nb, j_nb, ni, nj, max_depth - 1)
    else  # j_virt > nj
        nb = PANEL_CONNECTIVITY[p][North]
        depth = j_virt - nj
        i_nb, j_nb = transform_ghost_index(nb, depth, i_virt, ni, nj, North)
        return _resolve_corner_value(q, nb.neighbor_panel, i_nb, j_nb, ni, nj, max_depth - 1)
    end
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

"""
    ghost_fill_arrayop(u_interior, grid, loc=CellCenter)

Extend `u_interior` with ghost cells from neighbor panels, returning the
ghost-extended array suitable for use with ArrayOp-based operators.

This is the recommended entry point for preparing data for ArrayOp operators
like `flux_1d_ppm_arrayop`, `fv_laplacian_extended`, and
`transport_2d_ppm_arrayop`, which expect ghost-extended input arrays.

Returns the extended array of size (6, ni+2Ng, nj+2Ng).
"""
function ghost_fill_arrayop(u_interior, grid::CubedSphereGrid, loc::VarLocation = CellCenter)
    return extend_with_ghosts(u_interior, grid, loc)
end
