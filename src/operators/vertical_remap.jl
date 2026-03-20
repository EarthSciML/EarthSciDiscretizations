"""
Conservative vertical remapping using PPM (Lin 2004).

Remaps a quantity `q` from an old set of layer pressure thicknesses to a new set,
conserving the column-integrated mass of `q * dp`.
"""

"""
    vertical_remap(q, dp_old, dp_new, Nk)

Remap `q` (6, Nc, Nc, Nk) from layers defined by `dp_old` (6, Nc, Nc, Nk) to
layers defined by `dp_new` (6, Nc, Nc, Nk), returning the remapped field.

Uses PPM with Colella-Woodward limiting for the sub-grid reconstruction
within each old-grid layer, then integrates the parabolic profile
across the new-grid layers for exact conservation.
"""
function vertical_remap(q, dp_old, dp_new, Nk::Int)
    sz = size(q)
    length(sz) >= 4 || error("q must have at least 4 dimensions (6, Nc, Nc, Nk)")
    sz[4] == Nk || error("4th dimension of q must equal Nk=$Nk, got $(sz[4])")

    q_new = similar(q)

    for p in axes(q, 1), i in axes(q, 2), j in axes(q, 3)
        q_col = view(q, p, i, j, :)
        dp_old_col = view(dp_old, p, i, j, :)
        dp_new_col = view(dp_new, p, i, j, :)
        q_new_col = view(q_new, p, i, j, :)

        _remap_column!(q_new_col, q_col, dp_old_col, dp_new_col, Nk)
    end

    return q_new
end

"""
Remap a single column using PPM conservative remapping.
"""
function _remap_column!(q_new, q_old, dp_old, dp_new, Nk)
    # Step 1: Build cumulative pressure coordinates
    p_old = zeros(Nk + 1)  # pressure at old layer interfaces
    p_new = zeros(Nk + 1)  # pressure at new layer interfaces
    for k in 1:Nk
        p_old[k + 1] = p_old[k] + dp_old[k]
        p_new[k + 1] = p_new[k] + dp_new[k]
    end

    # Step 2: PPM reconstruction in each old layer
    # Compute interface values using 4th-order formula
    q_left = zeros(Nk)
    q_right = zeros(Nk)
    for k in 1:Nk
        if k == 1 || k == Nk
            q_left[k] = q_old[k]
            q_right[k] = q_old[k]
        elseif k == 2 || k == Nk - 1
            q_left[k] = 0.5 * (q_old[k - 1] + q_old[k])
            q_right[k] = 0.5 * (q_old[k] + q_old[k + 1])
            q_left[k], q_right[k] = _ppm_limit_cw84(q_left[k], q_right[k], q_old[k])
        else
            ql = (7.0 / 12.0) * (q_old[k - 1] + q_old[k]) - (1.0 / 12.0) * (q_old[k - 2] + q_old[k + 1])
            qr = (7.0 / 12.0) * (q_old[k] + q_old[k + 1]) - (1.0 / 12.0) * (q_old[k - 1] + q_old[k + 2])
            q_left[k], q_right[k] = _ppm_limit_cw84(ql, qr, q_old[k])
        end
    end

    # Step 3: Sweep through new layers, integrating PPM profiles from old layers
    k_old = 1  # current old layer
    for k_new in 1:Nk
        p_bot = p_new[k_new]
        p_top = p_new[k_new + 1]
        dp = p_top - p_bot

        if dp < 1e-30
            q_new[k_new] = q_old[min(k_old, Nk)]
            continue
        end

        integral = 0.0
        p_lo = p_bot

        while p_lo < p_top - 1e-30
            # Find which old layer contains p_lo
            while k_old < Nk && p_old[k_old + 1] <= p_lo + 1e-30
                k_old += 1
            end

            p_hi = min(p_top, p_old[k_old + 1])
            dp_old_k = dp_old[k_old]

            if dp_old_k < 1e-30
                integral += q_old[k_old] * (p_hi - p_lo)
            else
                # Integrate PPM parabola over [p_lo, p_hi] within old layer k_old
                # Map to local coordinate x ∈ [0,1] within old layer
                x_lo = (p_lo - p_old[k_old]) / dp_old_k
                x_hi = (p_hi - p_old[k_old]) / dp_old_k

                ql = q_left[k_old]; qr = q_right[k_old]; qi = q_old[k_old]
                dq = qr - ql
                q6 = 6.0 * (qi - 0.5 * (ql + qr))

                # ∫_x_lo^x_hi [ql + x*(dq + q6*(1-x))] dx
                F(x) = ql * x + 0.5 * x^2 * (dq + q6) - q6 * x^3 / 3.0
                integral += (F(x_hi) - F(x_lo)) * dp_old_k
            end

            p_lo = p_hi
        end

        q_new[k_new] = integral / dp
    end

    return q_new
end

"""
    vertical_remap_tendency(q, dp_old, dp_new, Nk)

Deprecated: use `vertical_remap` instead.
"""
function vertical_remap_tendency(q, dp_old, dp_new, Nk::Int)
    Base.depwarn("`vertical_remap_tendency` is deprecated, use `vertical_remap` instead.", :vertical_remap_tendency)
    return vertical_remap(q, dp_old, dp_new, Nk)
end
