# Oracle golden generation for the covariant FV Laplacian/gradient (esd-zk9.1, R1).
#
# Drives the ESS numeric oracle `precompute_laplacian_stencil` /
# `precompute_gradient_stencil` (ESS src/grid_assembly.jl @ 0e935a3) on two
# curvilinear metrics and serializes the EXACT Float64 stencil weights,
# neighbor tables, and applied results as STATIC JSON committed in ESD.
# These survive the later ESS grid_assembly deletion and are the conformance
# target for the R2 declarative arrayop-einsum rule.
#
# Mock grid: implements the AbstractCurvilinearGrid trait with caller-supplied
# metric arrays (mirroring ESS test/grid_assembly_symbolic_test.jl's
# SymTestCartesianGrid), so the metric — including a NON-ORTHOGONAL g^{xi eta}
# that a real orthogonal lat-lon grid would zero out — is controlled exactly.

using EarthSciSerialization
using EarthSciSerialization: AbstractCurvilinearGrid, n_cells, n_dims, axis_names,
    cell_centers, cell_widths, cell_volume, neighbor_indices, boundary_mask,
    metric_g, metric_ginv, metric_jacobian, metric_dgij_dxk,
    coord_jacobian, coord_jacobian_second,
    precompute_laplacian_stencil, precompute_gradient_stencil
using JSON

const ESS_COMMIT = "0e935a3384a6f13d396598a7dfebe2197e20b22c"

# ---------------------------------------------------------------------------
# Mock curvilinear grid with caller-supplied metric. Flat index c = i + (j-1)*Nx.
# Axes are (:xi, :eta). Periodicity per axis controls boundary sentinels (0),
# which the oracle maps to self (zero-displacement) per its boundary policy.
# ---------------------------------------------------------------------------
struct MockGrid <: AbstractCurvilinearGrid
    Nx::Int
    Ny::Int
    periodic_x::Bool
    periodic_y::Bool
    xis::Vector{Float64}        # cell-center xi per flat cell (N,)
    etas::Vector{Float64}       # cell-center eta per flat cell (N,)
    dxi::Float64
    deta::Float64
    ginv::Array{Float64, 3}      # (N,2,2) inverse metric g^{ij}
    Jac::Vector{Float64}        # (N,) jacobian sqrt(g)
    cj::Array{Float64, 3}        # (N,2,2) coord_jacobian d(xi_k)/d(t_l)
end

_flat(g::MockGrid, i, j) = i + (j - 1) * g.Nx
_wrap(idx, n) = mod(idx - 1, n) + 1

EarthSciSerialization.n_cells(g::MockGrid) = g.Nx * g.Ny
EarthSciSerialization.n_dims(g::MockGrid) = 2
EarthSciSerialization.axis_names(g::MockGrid) = (:xi, :eta)

EarthSciSerialization.cell_centers(g::MockGrid, axis::Symbol) =
    axis === :xi ? copy(g.xis) : axis === :eta ? copy(g.etas) :
    throw(ArgumentError("unknown axis $axis"))

EarthSciSerialization.cell_widths(g::MockGrid, axis::Symbol) =
    axis === :xi ? fill(g.dxi, n_cells(g)) : axis === :eta ? fill(g.deta, n_cells(g)) :
    throw(ArgumentError("unknown axis $axis"))

EarthSciSerialization.cell_volume(g::MockGrid) = g.Jac .* (g.dxi * g.deta)

function EarthSciSerialization.neighbor_indices(g::MockGrid, axis::Symbol, offset::Int)
    N = n_cells(g); out = Vector{Int}(undef, N)
    for j in 1:g.Ny, i in 1:g.Nx
        c = _flat(g, i, j)
        if axis === :xi
            ni = i + offset
            if g.periodic_x
                out[c] = _flat(g, _wrap(ni, g.Nx), j)
            else
                out[c] = (1 <= ni <= g.Nx) ? _flat(g, ni, j) : 0   # boundary sentinel
            end
        elseif axis === :eta
            nj = j + offset
            if g.periodic_y
                out[c] = _flat(g, i, _wrap(nj, g.Ny))
            else
                out[c] = (1 <= nj <= g.Ny) ? _flat(g, i, nj) : 0
            end
        else
            throw(ArgumentError("unknown axis $axis"))
        end
    end
    return out
end

EarthSciSerialization.boundary_mask(g::MockGrid, ::Symbol, ::Symbol) = falses(n_cells(g))
EarthSciSerialization.metric_ginv(g::MockGrid) = g.ginv
EarthSciSerialization.metric_jacobian(g::MockGrid) = g.Jac
function EarthSciSerialization.metric_g(g::MockGrid)
    N = n_cells(g); out = zeros(N, 2, 2)
    @inbounds for c in 1:N
        a = g.ginv[c, 1, 1]; b = g.ginv[c, 1, 2]; d = g.ginv[c, 2, 2]
        det = a * d - b * b
        out[c, 1, 1] = d / det; out[c, 2, 2] = a / det
        out[c, 1, 2] = -b / det; out[c, 2, 1] = -b / det
    end
    return out
end
EarthSciSerialization.metric_dgij_dxk(g::MockGrid) = zeros(n_cells(g), 2, 2, 2)
EarthSciSerialization.coord_jacobian(g::MockGrid, ::Symbol) = g.cj
EarthSciSerialization.coord_jacobian_second(g::MockGrid, ::Symbol) = zeros(n_cells(g), 2, 2, 2)

# ---------------------------------------------------------------------------
# Numeric stencil apply (reference; the public apply path was retired upstream).
# ---------------------------------------------------------------------------
function apply_stencil(weights, neighbors, u)
    N, K = size(neighbors); du = zeros(N)
    @inbounds for c in 1:N, k in 1:K
        du[c] += weights[c, k] * u[neighbors[c, k]]
    end
    return du
end

rows(M) = [collect(@view M[c, :]) for c in 1:size(M, 1)]                 # (N,K) -> Vector of length-K rows
ginv_cols(G) = (
    gxx = [G[c, 1, 1] for c in 1:size(G, 1)],
    gxe = [G[c, 1, 2] for c in 1:size(G, 1)],
    gyy = [G[c, 2, 2] for c in 1:size(G, 1)],
)
cj_cols(C) = (
    dxi_dt1 = [C[c, 1, 1] for c in 1:size(C, 1)], deta_dt1 = [C[c, 2, 1] for c in 1:size(C, 1)],
    dxi_dt2 = [C[c, 1, 2] for c in 1:size(C, 1)], deta_dt2 = [C[c, 2, 2] for c in 1:size(C, 1)],
)

function run_case(g::MockGrid, phi::Vector{Float64}, target::Symbol)
    lap = precompute_laplacian_stencil(g; xi_axis = :xi, eta_axis = :eta)
    grad = precompute_gradient_stencil(g, target; xi_axis = :xi, eta_axis = :eta)
    du_lap = apply_stencil(lap.weights, lap.neighbors, phi)
    gt1 = apply_stencil(grad.weights_t1, grad.neighbors, phi)
    gt2 = apply_stencil(grad.weights_t2, grad.neighbors, phi)
    gi = ginv_cols(g.ginv); cjc = cj_cols(g.cj)
    return Dict(
        "grid" => Dict(
            "Nx" => g.Nx, "Ny" => g.Ny, "periodic_xi" => g.periodic_x,
            "periodic_eta" => g.periodic_y, "dxi" => g.dxi, "deta" => g.deta,
            "flat_index" => "c = i + (j-1)*Nx, i in 1:Nx, j in 1:Ny",
            "xi_centers" => g.xis, "eta_centers" => g.etas
        ),
        "metric" => Dict(
            "ginv_xixi" => gi.gxx, "ginv_xieta" => gi.gxe, "ginv_etaeta" => gi.gyy,
            "jacobian_J" => g.Jac
        ),
        "gradient_target" => String(target),
        "coord_jacobian" => Dict(
            "dxi_dt1" => cjc.dxi_dt1, "deta_dt1" => cjc.deta_dt1,
            "dxi_dt2" => cjc.dxi_dt2, "deta_dt2" => cjc.deta_dt2
        ),
        "field_phi" => phi,
        "laplacian" => Dict(
            "stencil_columns" => ["C", "E(+xi)", "W(-xi)", "N(+eta)", "S(-eta)", "NE", "NW", "SE", "SW"],
            "weights" => rows(lap.weights), "neighbors" => rows(lap.neighbors),
            "applied_result" => du_lap
        ),
        "gradient" => Dict(
            "stencil_columns" => ["C", "E(+xi)", "W(-xi)", "N(+eta)", "S(-eta)"],
            "weights_t1" => rows(grad.weights_t1), "weights_t2" => rows(grad.weights_t2),
            "neighbors" => rows(grad.neighbors),
            "applied_result_t1" => gt1, "applied_result_t2" => gt2
        ),
        "checksums" => Dict(
            "sum_du_lap" => sum(du_lap), "sumabs_du_lap" => sum(abs, du_lap), "maxabs_du_lap" => maximum(abs, du_lap),
            "sumabs_grad_t1" => sum(abs, gt1), "sumabs_grad_t2" => sum(abs, gt2)
        ),
    )
end

# ---------------------------------------------------------------------------
# Case A — latlon (orthogonal spherical metric). xi=lon in [0,2pi) periodic;
# eta=lat in (-pi/2, pi/2) NON-periodic (poles -> boundary sentinel -> self).
# g^{lonlon}=1/(R^2 cos^2 phi), g^{latlat}=1/R^2, g^{lonlat}=0, J=R^2 cos phi.
# coord_jacobian(:lon_lat) = identity (matches ESD LatLonGrid).
# ---------------------------------------------------------------------------
function build_latlon(; Nx = 8, Ny = 6, R = 1.0)
    dxi = 2pi / Nx; deta = pi / Ny
    N = Nx * Ny
    xis = zeros(N); etas = zeros(N)
    ginv = zeros(N, 2, 2); Jac = zeros(N); cj = zeros(N, 2, 2)
    for j in 1:Ny, i in 1:Nx
        c = i + (j - 1) * Nx
        lon = (i - 0.5) * dxi                      # [0,2pi)
        lat = -pi / 2 + (j - 0.5) * deta             # interior, avoids exact poles
        xis[c] = lon; etas[c] = lat
        cphi = cos(lat)
        ginv[c, 1, 1] = 1 / (R^2 * cphi^2); ginv[c, 2, 2] = 1 / R^2; ginv[c, 1, 2] = 0.0; ginv[c, 2, 1] = 0.0
        Jac[c] = R^2 * cphi
        cj[c, 1, 1] = 1.0; cj[c, 2, 2] = 1.0           # identity
    end
    g = MockGrid(Nx, Ny, true, false, xis, etas, dxi, deta, ginv, Jac, cj)
    phi = [sin(xis[c]) * cos(etas[c]) for c in 1:N]
    return g, phi
end

# ---------------------------------------------------------------------------
# Case B — non-orthogonal curvilinear (synthetic smooth SPD metric, g^{xi eta}!=0).
# Exercises EVERY weight term incl. the 4 corner cross-derivative weights.
# Periodic in BOTH axes (clean, no boundary special-casing). xi,eta in [0,2pi).
# g^{xixi}=1+0.2cos(xi), g^{etaeta}=1+0.2sin(eta), g^{xieta}=0.25 sin(xi+eta),
# J = 1/sqrt(det(g^inv))  (physically consistent: J=sqrt(det g)).
# coord_jacobian: a smooth non-identity map so the gradient chain rule is exercised.
# ---------------------------------------------------------------------------
function build_nonorthogonal(; Nx = 8, Ny = 8)
    dxi = 2pi / Nx; deta = 2pi / Ny
    N = Nx * Ny
    xis = zeros(N); etas = zeros(N)
    ginv = zeros(N, 2, 2); Jac = zeros(N); cj = zeros(N, 2, 2)
    for j in 1:Ny, i in 1:Nx
        c = i + (j - 1) * Nx
        xi = (i - 0.5) * dxi; eta = (j - 0.5) * deta
        xis[c] = xi; etas[c] = eta
        gxx = 1.0 + 0.2 * cos(xi)
        gyy = 1.0 + 0.2 * sin(eta)
        gxe = 0.25 * sin(xi + eta)
        ginv[c, 1, 1] = gxx; ginv[c, 2, 2] = gyy; ginv[c, 1, 2] = gxe; ginv[c, 2, 1] = gxe
        det = gxx * gyy - gxe^2                  # > 0 (checked: >= ~0.57)
        Jac[c] = 1 / sqrt(det)
        # smooth non-identity coord_jacobian d(xi_k)/d(t_l)
        cj[c, 1, 1] = 1.0 + 0.1 * cos(eta); cj[c, 2, 1] = 0.15 * sin(xi)
        cj[c, 1, 2] = 0.15 * sin(eta);    cj[c, 2, 2] = 1.0 + 0.1 * cos(xi)
    end
    g = MockGrid(Nx, Ny, true, true, xis, etas, dxi, deta, ginv, Jac, cj)
    phi = [sin(xis[c]) * sin(etas[c]) for c in 1:N]
    return g, phi
end

# ---------------------------------------------------------------------------
outdir = ARGS[1]
mkpath(outdir)

g_ll, phi_ll = build_latlon()
caseA = run_case(g_ll, phi_ll, :lon_lat)
caseA["meta"] = Dict(
    "name" => "covariant_fv_oracle_latlon_small",
    "ess_commit" => ESS_COMMIT, "oracle" => "precompute_laplacian_stencil / precompute_gradient_stencil (ESS src/grid_assembly.jl)",
    "grid_kind" => "latlon (orthogonal spherical, R=1)", "R" => 1.0,
    "note" => "Orthogonal: g^{xieta}=0 so the 4 corner weights are exactly 0 (documents that lat-lon needs no corner terms; the eta metric-derivative correction IS exercised). Boundary: lon periodic, lat poles -> sentinel -> self."
)

g_no, phi_no = build_nonorthogonal()
caseB = run_case(g_no, phi_no, :physical)
caseB["meta"] = Dict(
    "name" => "covariant_fv_oracle_nonorthogonal_small",
    "ess_commit" => ESS_COMMIT, "oracle" => "precompute_laplacian_stencil / precompute_gradient_stencil (ESS src/grid_assembly.jl)",
    "grid_kind" => "non-orthogonal curvilinear (synthetic smooth SPD metric, g^{xieta}!=0)",
    "note" => "Exercises the full 9-point covariant stencil incl. the 4 NE/NW/SE/SW corner cross-derivative weights and the cross/orthogonal metric-derivative corrections. Periodic in both axes."
)

open(joinpath(outdir, "latlon_small.json"), "w") do io
    JSON.print(io, caseA, 2)
end
open(joinpath(outdir, "nonorthogonal_small.json"), "w") do io
    JSON.print(io, caseB, 2)
end

# Console summary for the verdict writeup.
function summary(tag, c)
    lapw = c["laplacian"]["weights"]
    corners = [abs(lapw[r][k]) for r in 1:length(lapw), k in 6:9]
    return println(
        "[$tag] N=$(c["grid"]["Nx"] * c["grid"]["Ny"]) maxabs|corner weight|=",
        maximum(corners), "  sum_du_lap=", c["checksums"]["sum_du_lap"],
        "  maxabs_du_lap=", c["checksums"]["maxabs_du_lap"]
    )
end
summary("latlon", caseA)
summary("nonorthogonal", caseB)
println("WROTE: ", joinpath(outdir, "latlon_small.json"), " and nonorthogonal_small.json")
