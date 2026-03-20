"""
Metric tensors and geometric quantities for the gnomonic cubed-sphere projection.
"""

function gnomonic_to_lonlat(ξ, η, panel)
    X = tan(ξ); Y = tan(η); R = sqrt(1 + X^2 + Y^2)
    if panel == 1;     x3d = 1/R;  y3d = X/R;  z3d = Y/R
    elseif panel == 2; x3d = -X/R; y3d = 1/R;  z3d = Y/R
    elseif panel == 3; x3d = -Y/R; y3d = X/R;  z3d = 1/R
    elseif panel == 4; x3d = -1/R; y3d = -X/R; z3d = Y/R
    elseif panel == 5; x3d = X/R;  y3d = -1/R; z3d = Y/R
    elseif panel == 6; x3d = Y/R;  y3d = X/R;  z3d = -1/R
    else error("Invalid panel: $panel")
    end
    return (atan(y3d, x3d), asin(clamp(z3d, -1.0, 1.0)))
end

function gnomonic_metric(ξ, η, R)
    X = tan(ξ); Y = tan(η); D2 = 1 + X^2 + Y^2
    secξ2 = 1 + X^2; secη2 = 1 + Y^2
    J = R^2 * secξ2 * secη2 / D2^(3 / 2)
    g_ξξ = R^2 * secξ2^2 * secη2 / D2^2
    g_ηη = R^2 * secη2^2 * secξ2 / D2^2
    g_ξη = -R^2 * X * Y * secξ2 * secη2 / D2^2
    return (J, g_ξξ, g_ηη, g_ξη)
end

function gnomonic_to_cart(ξ, η, panel)
    X = tan(ξ); Y = tan(η); R = sqrt(1 + X^2 + Y^2)
    if panel == 1;     return [1/R, X/R, Y/R]
    elseif panel == 2; return [-X/R, 1/R, Y/R]
    elseif panel == 3; return [-Y/R, X/R, 1/R]
    elseif panel == 4; return [-1/R, -X/R, Y/R]
    elseif panel == 5; return [X/R, -1/R, Y/R]
    elseif panel == 6; return [Y/R, X/R, -1/R]
    end
end

function compute_cell_area(ξ_edges, η_edges, R, panel)
    ξw, ξe = ξ_edges; ηs, ηn = η_edges
    corners = [gnomonic_to_cart(ξw, ηs, panel), gnomonic_to_cart(ξe, ηs, panel),
               gnomonic_to_cart(ξe, ηn, panel), gnomonic_to_cart(ξw, ηn, panel)]
    n = 4
    angles = Float64[]
    for k in 1:n
        v_prev = corners[mod1(k - 1, n)]; v_curr = corners[k]; v_next = corners[mod1(k + 1, n)]
        t1 = cross(v_curr, cross(v_prev - v_curr, v_curr))
        t2 = cross(v_curr, cross(v_next - v_curr, v_curr))
        n1 = norm(t1); n2 = norm(t2)
        push!(angles, (n1 < 1e-15 || n2 < 1e-15) ? π / 2 : acos(clamp(dot(t1, t2) / (n1 * n2), -1.0, 1.0)))
    end
    return R^2 * (sum(angles) - (n - 2) * π)
end

function compute_edge_length(ξ1, η1, ξ2, η2, R, panel)
    v1 = gnomonic_to_cart(ξ1, η1, panel); v2 = gnomonic_to_cart(ξ2, η2, panel)
    return R * acos(clamp(dot(v1, v2), -1.0, 1.0))
end

function compute_rotation_angle(panel_src, edge_src, panel_dst, edge_dst)
    same_type = (edge_src in (West, East)) == (edge_dst in (West, East))
    if same_type
        return 0.0
    else
        if (edge_src == North && edge_dst == East) || (edge_src == East && edge_dst == North) ||
           (edge_src == South && edge_dst == West) || (edge_src == West && edge_dst == South)
            return π / 2
        else
            return -π / 2
        end
    end
end

"""
    _panel_abc(X, Y, panel)

Return the unnormalized Cartesian coordinates (a, b, c) and their derivatives
w.r.t. X = tan(ξ) and Y = tan(η) for the gnomonic projection on `panel`.

The 3D point on the unit sphere is (x, y, z) = (a, b, c) / √(1 + X² + Y²).
Returns (a, b, c, da_dX, da_dY, db_dX, db_dY, dc_dX, dc_dY).
"""
function _panel_abc(X, Y, panel)
    if panel == 1     # (a,b,c) = (1, X, Y)
        return ( 1.0,  X,    Y,    0.0, 0.0,  1.0, 0.0,  0.0, 1.0)
    elseif panel == 2 # (a,b,c) = (-X, 1, Y)
        return (-X,    1.0,  Y,   -1.0, 0.0,  0.0, 0.0,  0.0, 1.0)
    elseif panel == 3 # (a,b,c) = (-Y, X, 1)
        return (-Y,    X,    1.0,  0.0,-1.0,  1.0, 0.0,  0.0, 0.0)
    elseif panel == 4 # (a,b,c) = (-1, -X, Y)
        return (-1.0, -X,    Y,    0.0, 0.0, -1.0, 0.0,  0.0, 1.0)
    elseif panel == 5 # (a,b,c) = (X, -1, Y)
        return ( X,   -1.0,  Y,    1.0, 0.0,  0.0, 0.0,  0.0, 1.0)
    elseif panel == 6 # (a,b,c) = (Y, X, -1)
        return ( Y,    X,   -1.0,  0.0, 1.0,  1.0, 0.0,  0.0, 0.0)
    else
        error("Invalid panel: $panel")
    end
end

"""
    compute_forward_jacobian(ξ, η, panel)

Compute the forward Jacobian d(lon,lat)/d(ξ,η) analytically.
Returns (dlon_dξ, dlon_dη, dlat_dξ, dlat_dη).

Uses the exact derivatives of the gnomonic projection composed with the
spherical coordinate formulas:
    lon = atan2(y, x)   →  dlon = (x dy - y dx) / (x² + y²)
    lat = asin(z)       →  dlat = dz / √(1 - z²)

At the poles (x² + y² = 0, z = ±1), the forward Jacobian is singular because
longitude is undefined. This is an inherent property of (lon, lat) coordinates,
not a deficiency of the projection.
"""
function compute_forward_jacobian(ξ, η, panel)
    X = tan(ξ); Y = tan(η)
    sx = 1 + X^2   # sec²(ξ)
    sy = 1 + Y^2   # sec²(η)
    δ = 1 + X^2 + Y^2

    a, b, c, da_dX, da_dY, db_dX, db_dY, dc_dX, dc_dY = _panel_abc(X, Y, panel)

    ab2 = a^2 + b^2   # = δ - c²; equals (x² + y²) · δ on the unit sphere

    # At the poles ab2 → 0: lon is undefined, return zero Jacobian.
    # The inverse Jacobian (compute_coord_jacobian) handles this gracefully.
    if ab2 < 1e-30
        return (dlon_dξ=0.0, dlon_dη=0.0, dlat_dξ=0.0, dlat_dη=0.0)
    end

    sqrt_ab2 = sqrt(ab2)

    # dlon/dξ = sec²(ξ) · (a · db/dX - b · da/dX) / (a² + b²)
    dlon_dξ = sx * (a * db_dX - b * da_dX) / ab2
    dlon_dη = sy * (a * db_dY - b * da_dY) / ab2

    # dlat/dξ = sec²(ξ) · (dc/dX · δ - c · X) / (δ · √(a² + b²))
    dlat_dξ = sx * (dc_dX * δ - c * X) / (δ * sqrt_ab2)
    dlat_dη = sy * (dc_dY * δ - c * Y) / (δ * sqrt_ab2)

    return (dlon_dξ=dlon_dξ, dlon_dη=dlon_dη, dlat_dξ=dlat_dξ, dlat_dη=dlat_dη)
end

"""
    compute_coord_jacobian(ξ, η, panel)

Compute the Jacobian of the inverse coordinate mapping (lon, lat) → (ξ, η)
at a given point. Returns (dξ_dlon, dξ_dlat, dη_dlon, dη_dlat).

This is needed to transform PDE derivatives written in (lon, lat) coordinates
to the computational (ξ, η) coordinate system via the chain rule:
    ∂/∂lon = (∂ξ/∂lon)·∂/∂ξ + (∂η/∂lon)·∂/∂η
    ∂/∂lat = (∂ξ/∂lat)·∂/∂ξ + (∂η/∂lat)·∂/∂η

Computed by analytically inverting the 2×2 forward Jacobian d(lon,lat)/d(ξ,η).
Near the poles the determinant vanishes because (lon,lat) coordinates are
singular there. A smooth regularization ensures the inverse remains bounded;
for typical even-Nc grids no cell center falls on a pole.
"""
function compute_coord_jacobian(ξ, η, panel)
    fwd = compute_forward_jacobian(ξ, η, panel)
    dlon_dξ = fwd.dlon_dξ; dlon_dη = fwd.dlon_dη
    dlat_dξ = fwd.dlat_dξ; dlat_dη = fwd.dlat_dη

    det = dlon_dξ * dlat_dη - dlon_dη * dlat_dξ

    # Smooth regularization near poles: use det² / (det² + ε²) scaling to
    # avoid a hard clamp discontinuity. For |det| >> ε this is ~1; for
    # |det| → 0 the result smoothly approaches zero.
    ε_reg = 1e-14
    if abs(det) < 1e-10
        scale = det / (det^2 + ε_reg)  # = 1/det · det²/(det² + ε) ≈ 1/det for |det| >> √ε
    else
        scale = 1.0 / det
    end

    dξ_dlon =  dlat_dη * scale
    dξ_dlat = -dlon_dη * scale
    dη_dlon = -dlat_dξ * scale
    dη_dlat =  dlon_dξ * scale

    return (dξ_dlon=dξ_dlon, dξ_dlat=dξ_dlat, dη_dlon=dη_dlon, dη_dlat=dη_dlat)
end

"""
    compute_second_coord_jacobian(ξ, η, panel)

Compute the second derivatives of the coordinate mapping (lon, lat) → (ξ, η):
    d²ξ/dlon², d²ξ/dlon·dlat, d²ξ/dlat², d²η/dlon², d²η/dlon·dlat, d²η/dlat².

Uses centered finite differences of the *analytical* inverse Jacobian
(single level of numerical differentiation, not triple-nested), combined
via the chain rule:
    d²ξ/dlon² = (dξ/dlon)·∂(dξ/dlon)/∂ξ + (dη/dlon)·∂(dξ/dlon)/∂η
"""
function compute_second_coord_jacobian(ξ, η, panel)
    h = 1e-6

    # Evaluate the analytical inverse Jacobian at stencil points
    jac_xp = compute_coord_jacobian(ξ + h, η, panel)
    jac_xm = compute_coord_jacobian(ξ - h, η, panel)
    jac_yp = compute_coord_jacobian(ξ, η + h, panel)
    jac_ym = compute_coord_jacobian(ξ, η - h, panel)

    # Derivatives of inverse Jacobian entries w.r.t. computational coordinates
    d_dξdlon_dξ = (jac_xp.dξ_dlon - jac_xm.dξ_dlon) / (2h)
    d_dξdlon_dη = (jac_yp.dξ_dlon - jac_ym.dξ_dlon) / (2h)
    d_dξdlat_dξ = (jac_xp.dξ_dlat - jac_xm.dξ_dlat) / (2h)
    d_dξdlat_dη = (jac_yp.dξ_dlat - jac_ym.dξ_dlat) / (2h)
    d_dηdlon_dξ = (jac_xp.dη_dlon - jac_xm.dη_dlon) / (2h)
    d_dηdlon_dη = (jac_yp.dη_dlon - jac_ym.dη_dlon) / (2h)
    d_dηdlat_dξ = (jac_xp.dη_dlat - jac_xm.dη_dlat) / (2h)
    d_dηdlat_dη = (jac_yp.dη_dlat - jac_ym.dη_dlat) / (2h)

    # Get the analytical inverse Jacobian at the center point
    jac0 = compute_coord_jacobian(ξ, η, panel)

    # Second derivatives via chain rule
    d2ξ_dlon2    = jac0.dξ_dlon * d_dξdlon_dξ  + jac0.dη_dlon * d_dξdlon_dη
    d2ξ_dlat2    = jac0.dξ_dlat * d_dξdlat_dξ  + jac0.dη_dlat * d_dξdlat_dη
    d2ξ_dlondlat = jac0.dξ_dlon * d_dξdlat_dξ  + jac0.dη_dlon * d_dξdlat_dη
    d2η_dlon2    = jac0.dξ_dlon * d_dηdlon_dξ  + jac0.dη_dlon * d_dηdlon_dη
    d2η_dlat2    = jac0.dξ_dlat * d_dηdlat_dξ  + jac0.dη_dlat * d_dηdlat_dη
    d2η_dlondlat = jac0.dξ_dlon * d_dηdlat_dξ  + jac0.dη_dlon * d_dηdlat_dη

    return (d2ξ_dlon2=d2ξ_dlon2, d2ξ_dlondlat=d2ξ_dlondlat, d2ξ_dlat2=d2ξ_dlat2,
            d2η_dlon2=d2η_dlon2, d2η_dlondlat=d2η_dlondlat, d2η_dlat2=d2η_dlat2)
end

"""
    tangent_vectors_3d(ξ, η, panel)

Compute the 3D tangent vectors ∂r/∂ξ and ∂r/∂η at (ξ, η) on `panel`.
Returns two 3-element vectors (e_ξ, e_η) in Cartesian (x, y, z) coordinates.
"""
function tangent_vectors_3d(ξ, η, panel)
    X = tan(ξ); Y = tan(η)
    sx = 1 + X^2; sy = 1 + Y^2
    δ = 1 + X^2 + Y^2
    denom = δ^(3 / 2)
    a, b, c, da_dX, da_dY, db_dX, db_dY, dc_dX, dc_dY = _panel_abc(X, Y, panel)
    e_ξ = [sx * (da_dX * δ - a * X) / denom,
           sx * (db_dX * δ - b * X) / denom,
           sx * (dc_dX * δ - c * X) / denom]
    e_η = [sy * (da_dY * δ - a * Y) / denom,
           sy * (db_dY * δ - b * Y) / denom,
           sy * (dc_dY * δ - c * Y) / denom]
    return (e_ξ, e_η)
end

"""
    compute_edge_rotation_matrix(panel, dir)

Compute the 2×2 rotation matrix that transforms contravariant vector
components (u_ξ, u_η) from a neighbor panel's basis into the local panel's
basis when crossing through edge `dir`.

Returns (M11, M12, M21, M22) such that:
    u_ξ_local = M11 * u_ξ_nb + M12 * u_η_nb
    u_η_local = M21 * u_ξ_nb + M22 * u_η_nb

The matrix is computed by comparing unit tangent vectors at the midpoint of
the shared edge. On the cubed sphere, these entries are exactly 0, ±1.
"""
function compute_edge_rotation_matrix(panel, dir)
    nb = PANEL_CONNECTIVITY[panel][dir]
    nb_panel = nb.neighbor_panel
    nb_edge = nb.neighbor_edge

    # Representative midpoint on the shared edge (along-edge coord = 0)
    _edge_point(d) = d == West  ? (-π / 4, 0.0) :
                     d == East  ? ( π / 4, 0.0) :
                     d == South ? (0.0, -π / 4) :
                                  (0.0,  π / 4)

    ξ_loc, η_loc = _edge_point(dir)
    ξ_nb,  η_nb  = _edge_point(nb_edge)

    e_ξ_l, e_η_l = tangent_vectors_3d(ξ_loc, η_loc, panel)
    e_ξ_n, e_η_n = tangent_vectors_3d(ξ_nb, η_nb, nb_panel)

    nξl = norm(e_ξ_l); nηl = norm(e_η_l)
    nξn = norm(e_ξ_n); nηn = norm(e_η_n)

    M11 = dot(e_ξ_l, e_ξ_n) / (nξl * nξn)
    M12 = dot(e_ξ_l, e_η_n) / (nξl * nηn)
    M21 = dot(e_η_l, e_ξ_n) / (nηl * nξn)
    M22 = dot(e_η_l, e_η_n) / (nηl * nηn)

    # Round to nearest integer (should be exactly 0, ±1 on cubed sphere)
    _snap(x) = abs(x) < 0.1 ? 0.0 : sign(x) * 1.0
    return (_snap(M11), _snap(M12), _snap(M21), _snap(M22))
end

"""
Wrap angle difference to [-π, π].
"""
function _wrap_angle(dθ)
    while dθ > π; dθ -= 2π; end
    while dθ < -π; dθ += 2π; end
    return dθ
end
