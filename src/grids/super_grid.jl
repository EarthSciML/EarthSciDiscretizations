"""
FV3 super-grid angular values (sin_sg, cos_sg) for non-orthogonal flux computation.

On the gnomonic cubed-sphere, the coordinate vectors e_ξ and e_η are generally
non-orthogonal. The angle α between them varies across each cell. FV3 precomputes
sin(α) and cos(α) at 9 sub-positions within each cell:

```
  9---4---8
  |       |
  1   5   3
  |       |
  6---2---7
```

Positions 1-4: mid-edge (west, south, east, north)
Position 5:   cell center
Positions 6-9: cell corners (SW, SE, NE, NW)

These are critical for computing the correct normal flux component through cell
faces on the non-orthogonal grid, especially at cube edges where α can change
discontinuously.

Reference: Harris et al. (2021), GFDL FV3 Technical Memorandum, Section 3.
"""

"""
    compute_angle_at_point(ξ, η, R)

Compute sin(α) and cos(α) at a point (ξ, η) on the gnomonic grid, where α is
the angle between the coordinate vectors e_ξ and e_η.

From the metric tensor:
    cos(α) = g_ξη / √(g_ξξ · g_ηη)
    sin(α) = J / √(g_ξξ · g_ηη)

where J = √(det(g)) = √(g_ξξ·g_ηη - g_ξη²).
"""
function compute_angle_at_point(ξ, η, R)
    J, g_ξξ, g_ηη, g_ξη = gnomonic_metric(ξ, η, R)
    denom = sqrt(g_ξξ * g_ηη)
    sin_α = J / denom
    cos_α = g_ξη / denom
    return (sin_α, cos_α)
end

"""
    compute_super_grid!(sin_sg, cos_sg, grid)

Compute sin_sg and cos_sg arrays at all 9 sub-positions for every cell.
Arrays have shape (6, Nc, Nc, 9).

The 9 positions within cell (i, j) are:
- 1: west mid-edge  (ξ_{i-1/2}, η_j)
- 2: south mid-edge (ξ_i, η_{j-1/2})
- 3: east mid-edge  (ξ_{i+1/2}, η_j)
- 4: north mid-edge (ξ_i, η_{j+1/2})
- 5: cell center    (ξ_i, η_j)
- 6: SW corner      (ξ_{i-1/2}, η_{j-1/2})
- 7: SE corner      (ξ_{i+1/2}, η_{j-1/2})
- 8: NE corner      (ξ_{i+1/2}, η_{j+1/2})
- 9: NW corner      (ξ_{i-1/2}, η_{j+1/2})
"""
function compute_super_grid!(sin_sg, cos_sg, grid)
    Nc = grid.Nc; R = grid.R
    ξc = grid.ξ_centers; ηc = grid.η_centers
    ξe = grid.ξ_edges; ηe = grid.η_edges

    for p in 1:6, i in 1:Nc, j in 1:Nc
        # Position 1: west mid-edge (ξ_{i-1/2}, η_j)
        s, c = compute_angle_at_point(ξe[i], ηc[j], R)
        sin_sg[p, i, j, 1] = s; cos_sg[p, i, j, 1] = c

        # Position 2: south mid-edge (ξ_i, η_{j-1/2})
        s, c = compute_angle_at_point(ξc[i], ηe[j], R)
        sin_sg[p, i, j, 2] = s; cos_sg[p, i, j, 2] = c

        # Position 3: east mid-edge (ξ_{i+1/2}, η_j)
        s, c = compute_angle_at_point(ξe[i+1], ηc[j], R)
        sin_sg[p, i, j, 3] = s; cos_sg[p, i, j, 3] = c

        # Position 4: north mid-edge (ξ_i, η_{j+1/2})
        s, c = compute_angle_at_point(ξc[i], ηe[j+1], R)
        sin_sg[p, i, j, 4] = s; cos_sg[p, i, j, 4] = c

        # Position 5: cell center (ξ_i, η_j)
        s, c = compute_angle_at_point(ξc[i], ηc[j], R)
        sin_sg[p, i, j, 5] = s; cos_sg[p, i, j, 5] = c

        # Position 6: SW corner (ξ_{i-1/2}, η_{j-1/2})
        s, c = compute_angle_at_point(ξe[i], ηe[j], R)
        sin_sg[p, i, j, 6] = s; cos_sg[p, i, j, 6] = c

        # Position 7: SE corner (ξ_{i+1/2}, η_{j-1/2})
        s, c = compute_angle_at_point(ξe[i+1], ηe[j], R)
        sin_sg[p, i, j, 7] = s; cos_sg[p, i, j, 7] = c

        # Position 8: NE corner (ξ_{i+1/2}, η_{j+1/2})
        s, c = compute_angle_at_point(ξe[i+1], ηe[j+1], R)
        sin_sg[p, i, j, 8] = s; cos_sg[p, i, j, 8] = c

        # Position 9: NW corner (ξ_{i-1/2}, η_{j+1/2})
        s, c = compute_angle_at_point(ξe[i], ηe[j+1], R)
        sin_sg[p, i, j, 9] = s; cos_sg[p, i, j, 9] = c
    end
    return (sin_sg, cos_sg)
end

"""
    compute_super_grid(grid)

Allocate and compute the sin_sg and cos_sg arrays for the grid.
Returns (sin_sg, cos_sg), each of shape (6, Nc, Nc, 9).
"""
function compute_super_grid(grid)
    Nc = grid.Nc
    T = eltype(grid.area)
    sin_sg = zeros(T, 6, Nc, Nc, 9)
    cos_sg = zeros(T, 6, Nc, Nc, 9)
    compute_super_grid!(sin_sg, cos_sg, grid)
    return (sin_sg, cos_sg)
end
