using EarthSciSerialization: apply_weno5_reconstruction_periodic_1d
using JSON

"""
    run_weno5_2d_smooth_advection_convergence(name, manifest) -> Tuple{Symbol,String}

Layer-C operational substitute for the 2D WENO5 layer-B sweep until ESS gains a
2D dispatch (see `../convergence/input.esm` skip_reason). The rule applies 1D
Jiang-Shu (1996) WENO5 reconstruction dimension-by-dimension (Shu 1998 §2.2 /
LeVeque 2002 §20). This runner exercises that semantics directly:

1. Resolve the rule JSON from the catalog (sibling of `fixtures/`).
2. For each grid `n` in `manifest["grids"]`, build cell averages of the
   manufactured solution `u(x,y) = sin(2π x + φ_x)·sin(2π y + φ_y)` on the
   [0,1]² periodic square. Cell averages are analytic and separable.
3. For each row `j`, apply `apply_weno5_reconstruction_periodic_1d` with the
   rule's `axes.x.reconstruction_<side>` block to recover face values
   `q_{i+1/2,j}^{side}`. Compare to the cross-section average truth at that
   face (the y-cell-average evaluated at `x_{i+1/2}`). Repeat per column
   along y using `axes.y.reconstruction_<side>`. The L∞ error is the max
   over both axes.
4. Compute refinement orders `log2(e_k / e_{k+1})` and take the minimum.
5. Pass iff `observed_min_order >= manifest["tolerance"]["min_order"]`.

This case fails (NOT skips) if the observed order falls below the threshold;
that is the trigger for escalation per dsc-5od acceptance #5.
"""
function run_weno5_2d_smooth_advection_convergence(
        name::AbstractString, manifest::AbstractDict,
    )
    rule_path = _resolve_weno5_2d_rule_path(manifest)
    rule_path === nothing && return (
        :fail,
        "$name: cannot locate weno5_advection_2d.json in catalog (looked under " *
        "discretizations/finite_volume/)",
    )
    rule_json = JSON.parse(read(rule_path, String))
    spec = rule_json["discretizations"]["weno5_advection_2d"]
    haskey(spec, "axes") || return (
        :fail,
        "$name: rule has no `axes` block (dimension-by-dimension WENO5 expected)",
    )
    axes = spec["axes"]
    haskey(axes, "x") && haskey(axes, "y") || return (
        :fail, "$name: rule.axes must carry both `x` and `y` blocks",
    )

    side_str = String(get(manifest, "reconstruction", "left_biased"))
    side = side_str == "left_biased" ? :left_biased :
           side_str == "right_biased" ? :right_biased :
           return (:fail, "$name: reconstruction must be left_biased or right_biased, got $(repr(side_str))")
    eps = Float64(get(manifest, "weno_epsilon", 1.0e-6))

    ms = manifest["manufactured_solution"]
    phi_x = Float64(get(ms, "phi_x", 1.0))
    phi_y = Float64(get(ms, "phi_y", 0.5))

    grids_raw = manifest["grids"]
    grids = Int[Int(g["n"]) for g in grids_raw]
    length(grids) >= 2 || return (
        :fail, "$name: convergence requires at least two grids; got $(grids)",
    )

    tol = manifest["tolerance"]
    min_order = Float64(tol["min_order"])

    # Cell average of sin(2π·t + φ) over [a, b] is
    # (cos(2π·a + φ) − cos(2π·b + φ)) / (2π · (b − a)).
    sin_avg(a, b, phi) = (cos(2π * a + phi) - cos(2π * b + phi)) / (2π * (b - a))

    errors = Float64[]
    n_per_grid = Int[]
    for n in grids
        h = 1.0 / n
        # Centers and edges along [0, 1] with periodic wrap.
        x_avg = [sin_avg((i - 1) * h, i * h, phi_x) for i in 1:n]
        y_avg = [sin_avg((j - 1) * h, j * h, phi_y) for j in 1:n]

        # qbar[i, j] = (cell avg of sin in x at row j) ·
        #              (cell avg of sin in y at col i).
        # By separability of the manufactured solution, the 2D cell average is
        # the product of the per-axis 1D cell averages.
        qbar = [x_avg[i] * y_avg[j] for i in 1:n, j in 1:n]

        # x-axis sweep: per row j, reconstruct face values at x_{i+1/2}; truth
        # = sin_face_x(i+1/2) · y_avg[j], where sin_face_x at x = i·h is
        # sin(2π·i·h + φ_x).
        err_x = 0.0
        x_face_truth = [sin(2π * (i * h) + phi_x) for i in 1:n]
        for j in 1:n
            row = qbar[:, j]
            qhat = apply_weno5_reconstruction_periodic_1d(
                axes["x"], row, side; eps = eps,
            )
            for i in 1:n
                truth_face = side === :left_biased ?
                    x_face_truth[i] * y_avg[j] :
                    sin(2π * ((i - 1) * h) + phi_x) * y_avg[j]
                err_x = max(err_x, abs(qhat[i] - truth_face))
            end
        end

        # y-axis sweep: per column i, reconstruct face values at y_{j+1/2};
        # truth = x_avg[i] · sin_face_y(j+1/2).
        err_y = 0.0
        y_face_truth = [sin(2π * (j * h) + phi_y) for j in 1:n]
        for i in 1:n
            col = qbar[i, :]
            qhat = apply_weno5_reconstruction_periodic_1d(
                axes["y"], col, side; eps = eps,
            )
            for j in 1:n
                truth_face = side === :left_biased ?
                    x_avg[i] * y_face_truth[j] :
                    x_avg[i] * sin(2π * ((j - 1) * h) + phi_y)
                err_y = max(err_y, abs(qhat[j] - truth_face))
            end
        end

        push!(errors, max(err_x, err_y))
        push!(n_per_grid, n)
    end

    if any(!isfinite, errors) || any(e -> e <= 0, errors)
        return (
            :fail,
            "$name: non-finite or zero error on some grid; errors=$(errors)",
        )
    end

    orders = [log2(errors[i] / errors[i + 1]) for i in 1:(length(errors) - 1)]
    observed_min = minimum(orders)
    rounded_orders = [round(o; digits = 3) for o in orders]
    rounded_min = round(observed_min; digits = 3)

    if observed_min < min_order
        return (
            :fail,
            "$name: observed_min_order $(rounded_min) < threshold $(min_order) " *
            "on grids $(n_per_grid) (orders=$(rounded_orders), errors=$(errors))",
        )
    end
    return (
        :pass,
        "$name: observed_min_order $(rounded_min) >= $(min_order) " *
        "on grids $(n_per_grid) (orders=$(rounded_orders))",
    )
end

# Locate the rule JSON so the runner is independent of the manifest path.
function _resolve_weno5_2d_rule_path(::AbstractDict)
    # IntegrationCases is loaded from `<repo>/test/integration_cases/`.
    repo_root = dirname(dirname(@__DIR__))
    candidate = joinpath(
        repo_root,
        "discretizations", "finite_volume", "weno5_advection_2d.json",
    )
    return isfile(candidate) ? candidate : nothing
end
