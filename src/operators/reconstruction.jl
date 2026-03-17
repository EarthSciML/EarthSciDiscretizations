"""
Piecewise Parabolic Method (PPM) sub-grid reconstruction (Colella-Woodward 1984).
"""

function ppm_reconstruction!(q_left, q_right, q, grid::CubedSphereGrid, dim::Symbol)
    Nc = grid.Nc
    if dim == :xi
        for p in 1:6, i in 3:Nc-2, j in 1:Nc
            qim1 = q[p, i - 1, j]; qi = q[p, i, j]; qip1 = q[p, i + 1, j]; qip2 = q[p, i + 2, j]; qim2 = q[p, i - 2, j]
            qr = (7.0 / 12.0) * (qi + qip1) - (1.0 / 12.0) * (qim1 + qip2)
            ql = (7.0 / 12.0) * (qim1 + qi) - (1.0 / 12.0) * (qim2 + qip1)
            q_right[p, i - 2, j] = _ppm_limit(ql, qr, qim1, qi, qip1, true)
            q_left[p, i - 2, j] = _ppm_limit(ql, qr, qim1, qi, qip1, false)
        end
    elseif dim == :eta
        for p in 1:6, i in 1:Nc, j in 3:Nc-2
            qim1 = q[p, i, j - 1]; qi = q[p, i, j]; qip1 = q[p, i, j + 1]; qip2 = q[p, i, j + 2]; qim2 = q[p, i, j - 2]
            qr = (7.0 / 12.0) * (qi + qip1) - (1.0 / 12.0) * (qim1 + qip2)
            ql = (7.0 / 12.0) * (qim1 + qi) - (1.0 / 12.0) * (qim2 + qip1)
            q_right[p, i, j - 2] = _ppm_limit(ql, qr, qim1, qi, qip1, true)
            q_left[p, i, j - 2] = _ppm_limit(ql, qr, qim1, qi, qip1, false)
        end
    end
end

function _ppm_limit(ql, qr, qim1, qi, qip1, is_right)
    val, lo, hi = is_right ? (qr, min(qi, qip1), max(qi, qip1)) : (ql, min(qim1, qi), max(qim1, qi))
    return clamp(val, lo, hi)
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
