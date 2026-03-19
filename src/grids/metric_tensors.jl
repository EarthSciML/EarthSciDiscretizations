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
    compute_coord_jacobian(ξ, η, panel)

Compute the Jacobian of the inverse coordinate mapping (lon, lat) → (ξ, η)
at a given point. Returns (dξ_dlon, dξ_dlat, dη_dlon, dη_dlat).

This is needed to transform PDE derivatives written in (lon, lat) coordinates
to the computational (ξ, η) coordinate system via the chain rule:
    ∂/∂lon = (∂ξ/∂lon)·∂/∂ξ + (∂η/∂lon)·∂/∂η
    ∂/∂lat = (∂ξ/∂lat)·∂/∂ξ + (∂η/∂lat)·∂/∂η
"""
function compute_coord_jacobian(ξ, η, panel)
    eps = 1e-7

    # Forward Jacobian d(lon,lat)/d(ξ,η) via central differences
    lon_p, lat_p = gnomonic_to_lonlat(ξ + eps, η, panel)
    lon_m, lat_m = gnomonic_to_lonlat(ξ - eps, η, panel)
    # Handle longitude wrapping
    dlon_dξ = _wrap_angle(lon_p - lon_m) / (2 * eps)
    dlat_dξ = (lat_p - lat_m) / (2 * eps)

    lon_p, lat_p = gnomonic_to_lonlat(ξ, η + eps, panel)
    lon_m, lat_m = gnomonic_to_lonlat(ξ, η - eps, panel)
    dlon_dη = _wrap_angle(lon_p - lon_m) / (2 * eps)
    dlat_dη = (lat_p - lat_m) / (2 * eps)

    # Invert the 2×2 forward Jacobian to get d(ξ,η)/d(lon,lat)
    det = dlon_dξ * dlat_dη - dlon_dη * dlat_dξ
    dξ_dlon = dlat_dη / det
    dξ_dlat = -dlon_dη / det
    dη_dlon = -dlat_dξ / det
    dη_dlat = dlon_dξ / det

    return (dξ_dlon=dξ_dlon, dξ_dlat=dξ_dlat, dη_dlon=dη_dlon, dη_dlat=dη_dlat)
end

"""
Wrap angle difference to [-π, π].
"""
function _wrap_angle(dθ)
    while dθ > π; dθ -= 2π; end
    while dθ < -π; dθ += 2π; end
    return dθ
end
