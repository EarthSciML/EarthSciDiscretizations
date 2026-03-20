"""
Piecewise Parabolic Method (PPM) sub-grid reconstruction (Colella & Woodward 1984).

Uses the 4th-order interface interpolation formula and the full Colella-Woodward
monotonicity limiter that preserves cell averages and ensures the parabolic
profile is monotone within each cell.
"""

"""
    ppm_reconstruction!(q_left, q_right, q, grid, dim)

In-place PPM reconstruction. Fills `q_left` and `q_right` with the
limited left and right interface values for each interior cell.
"""
function ppm_reconstruction!(q_left, q_right, q, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        for p in 1:6, j in 1:Nc
            # Step 1: Compute 4th-order interface values at all half-points
            # Interface i+1/2 is between cell i and cell i+1
            ifaces = Vector{Float64}(undef, Nc + 1)
            # Use 4th-order formula for interior interfaces (need i-1,i,i+1,i+2)
            for i in 2:Nc-1
                ifaces[i] = (7.0 / 12.0) * (q[p, i, j] + q[p, i - 1, j]) -
                            (1.0 / 12.0) * (q[p, max(1, i - 2), j] + q[p, min(Nc, i + 1), j])
            end
            # Boundary interfaces: use cell center value as edge estimate
            ifaces[1] = q[p, 1, j]
            ifaces[Nc + 1] = q[p, Nc, j]

            # Step 2: Assign left/right edges and apply CW84 limiter
            for i in 3:Nc-2
                ql = ifaces[i]      # left edge of cell i (interface i-1/2)
                qr = ifaces[i + 1]  # right edge of cell i (interface i+1/2)
                qi = q[p, i, j]

                ql, qr = _ppm_limit_cw84(ql, qr, qi)

                q_left[p, i - 2, j] = ql
                q_right[p, i - 2, j] = qr
            end
        end
    elseif dim == :eta
        for p in 1:6, i in 1:Nc
            ifaces = Vector{Float64}(undef, Nc + 1)
            for jj in 2:Nc-1
                ifaces[jj] = (7.0 / 12.0) * (q[p, i, jj] + q[p, i, jj - 1]) -
                             (1.0 / 12.0) * (q[p, i, max(1, jj - 2)] + q[p, i, min(Nc, jj + 1)])
            end
            if Nc >= 2
                ifaces[1] = q[p, i, 1]
                ifaces[Nc + 1] = q[p, i, Nc]
            end

            for jj in 3:Nc-2
                ql = ifaces[jj]
                qr = ifaces[jj + 1]
                qi = q[p, i, jj]

                ql, qr = _ppm_limit_cw84(ql, qr, qi)

                q_left[p, i, jj - 2] = ql
                q_right[p, i, jj - 2] = qr
            end
        end
    end
end

"""
Colella-Woodward (1984) PPM limiter.

Given initial left and right edge values (ql, qr) and the cell average qi,
apply monotonicity constraints:
1. If cell is a local extremum, flatten to qi.
2. If the parabola's extremum falls within the cell, adjust the offending edge.

Returns the limited (ql, qr) pair that ensures:
- The parabola is monotone between ql and qr (no internal extremum)
- The cell average is preserved: (ql + qr)/2 + q6/6 = qi
  where q6 = 6(qi - (ql+qr)/2)
"""
function _ppm_limit_cw84(ql, qr, qi)
    # Step 1: Local extremum detection
    if (qr - qi) * (qi - ql) <= 0
        ql = qi
        qr = qi
    else
        # Step 2: Overshoot correction
        dq = qr - ql
        q6 = 6.0 * (qi - 0.5 * (ql + qr))
        if dq * q6 > dq^2
            ql = 3.0 * qi - 2.0 * qr
        end
        if -dq^2 > dq * q6
            qr = 3.0 * qi - 2.0 * ql
        end
    end
    return (ql, qr)
end

function ppm_reconstruction(q, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc; T = eltype(q)
    if dim == :xi
        q_left = zeros(T, 6, Nc - 4, Nc); q_right = zeros(T, 6, Nc - 4, Nc)
    else
        q_left = zeros(T, 6, Nc, Nc - 4); q_right = zeros(T, 6, Nc, Nc - 4)
    end
    ppm_reconstruction!(q_left, q_right, q, grid, dim)
    return (q_left, q_right)
end

"""
    ppm_flux_integral(ql, qr, qi, courant)

Compute the PPM flux integral: the area-averaged value swept through an
interface over one time step. Uses the parabolic reconstruction to
integrate over the Courant-number-scaled fraction of the cell.

For positive courant (flow from left to right), integrates the parabola
in the upwind cell from the right edge inward. For negative courant,
integrates from the left edge inward.
"""
function ppm_flux_integral(ql, qr, qi, courant)
    # Parabola: q(x) = ql + x*(dq + q6*(1-x))  for x ∈ [0,1]
    # where dq = qr - ql, q6 = 6*(qi - (ql+qr)/2)
    dq = qr - ql
    q6 = 6.0 * (qi - 0.5 * (ql + qr))
    c = abs(courant)
    if c < 1e-15
        return 0.0
    end
    if c > 1.0
        @warn "CFL violation: |Courant| = $c > 1. PPM flux integral requires |Courant| < 1." maxlog=1
        c = min(c, 1.0)  # Clamp to prevent extrapolation beyond cell
    end
    # Integrate from the upwind edge
    if courant >= 0
        # Integrate from x = 1-c to x = 1 (right portion of upwind cell)
        # ∫_{1-c}^{1} q(x) dx / c
        return qr - 0.5 * c * (dq - q6 * (1.0 - 2.0 / 3.0 * c))
    else
        # Integrate from x = 0 to x = c (left portion of downwind cell)
        # ∫_{0}^{c} q(x) dx / c
        return ql + 0.5 * c * (dq + q6 * (1.0 - 2.0 / 3.0 * c))
    end
end
