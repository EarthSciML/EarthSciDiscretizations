# Staggered B-spline interpolating regridder — Julia evaluator conformance.
#
# ESD owns this declarative regridder PROGRAM (bspline_regrid.esm, co-located
# here); the EarthSciSerialization engine (an ESD dependency) EVALUATES it.
# Sibling to conservative_regrid_overlap_join_conformance_test.jl — the OTHER
# per-staggering regridding kernel (RFC pure-io-data-loaders §5.2: cell-centered
# fields → conservative; cell-edge/face-staggered fields → this B-spline kernel,
# the EarthSciData InterpolatingRegridder).
#
# The kernel is a `sum_product` interpolation FAQ over a fixed-offset B-spline
# stencil whose interpolation-weight factor is the STATIC-UNROLL degree-d Lagrange
# cardinal polynomial in the per-target local coordinate `s` (closed form, no
# de-Boor prefilter solve; no per-binding kernel — the same polynomial-
# reconstruction-as-arrayop idiom as discretizations/finite_volume/flux_1d_ppm.json).
# Degree-1 is the linear/bilinear special case. Each target value is
# F_tgt[j] = Σ_k w_k(s[j])·F_src[base[j]+k]: an index-gather (F_src at the host-
# derived base[j]+k offsets, folded at build time) and a polynomial weighted sum.
#
# Each model is a zero-IC constant-RHS D-equation, so du = f!(u0) at the zero IC
# IS the regridded field (the conservative-regridder assembly precedent). F_src /
# base / s are host-supplied const_arrays (referenced, not declared — the
# integral_probe Check-4 idiom for a non-geometry / non-value-invention model).
#
# RFC pure-io-data-loaders §8 (Regridding invariants): "B-spline: reproduction of
# the spline's polynomial degree exactly on staggered locations." This suite drives
# that invariant through build_evaluator, plus partition-of-unity, the exact
# 4th-order staggered face weights {-1/16, 9/16, 9/16, -1/16}, and the
# bilinear = degree-1-tensor-product identity.

@testsnippet BSplineRegridSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach it
    # through the dependency so the raw document is the exact type build_evaluator
    # expects (the value-invention / const_array front-door).
    const JSON3 = ESS.JSON3

    const _BS_FIXTURE = joinpath(@__DIR__, "bspline_regrid.esm")
    _bs_raw() = JSON3.read(read(_BS_FIXTURE, String))

    # Drive a model through build_evaluator and return the regridded F_tgt vector.
    # Every F_tgt[j] is a constant-RHS D-equation from a zero IC, so du = f!(u0) IS
    # the assembled value. `base` is supplied as Float64 (the index-fold rounds it).
    function _bs_eval(model_name, ntgt, const_arrays)
        ics = Dict("F_tgt[$j]" => 0.0 for j in 1:ntgt)
        f!, u0, p, _, vmap = ESS.build_evaluator(
            _bs_raw(); model_name = model_name,
            initial_conditions = ics, const_arrays = const_arrays
        )
        du = similar(u0); f!(du, u0, p, 0.0)
        return [du[vmap["F_tgt[$j]"]] for j in 1:ntgt]
    end

    # 1-D apply: F_tgt[j] = Σ_k w_k(s[j])·F_src[base[j]+k].
    _bs_eval_1d(model_name, Fsrc, base, s) = _bs_eval(
        model_name, length(base),
        Dict("F_src" => collect(Float64, Fsrc),
             "base" => collect(Float64, base), "s" => collect(Float64, s))
    )

    # 2-D tensor-product apply (bilinear).
    _bs_eval_2d(model_name, Fsrc2, bx, by, sx, sy) = _bs_eval(
        model_name, length(bx),
        Dict("F_src" => Matrix{Float64}(Fsrc2),
             "base_x" => collect(Float64, bx), "base_y" => collect(Float64, by),
             "s_x" => collect(Float64, sx), "s_y" => collect(Float64, sy))
    )

    # The 1-D models bake tgt_locs size = 4 (the equation unrolls j=1..4), so every
    # eval must supply 4 targets' worth of base/s. The reproduction tests pass 4
    # distinct targets; the weight probes replicate one probe target 4× and read [1].
    const _BS_NTGT = 4

    # Recover the (width) interpolation weights at a single target (base0, s0) by
    # feeding unit-impulse source fields — w_k = F_tgt when F_src = e_{base0+k}.
    # This reads the kernel's weights straight out of the evaluator.
    function _bs_weights_1d(model_name, nsrc, base0, s0, width)
        w = zeros(Float64, width)
        for k in 0:(width - 1)
            Fsrc = zeros(Float64, nsrc); Fsrc[base0 + k] = 1.0
            got = _bs_eval_1d(model_name, Fsrc, fill(base0, _BS_NTGT), fill(s0, _BS_NTGT))
            w[k + 1] = got[1]
        end
        return w
    end

    # Host stencil geometry on a UNIFORM source grid (node i at x = (i-1)·h, h=1).
    # LINEAR (width 2): bracketing interval [base, base+1]; base = floor(x)+1.
    _bs_linear_base(x) = floor(Int, x) + 1
    # CUBIC (width 4): centered stencil base..base+3, target in [base+1, base+2];
    # base = floor(x) (so base+1 = floor(x)+1 is the central-left node).
    _bs_cubic_base(x) = floor(Int, x)
    _bs_frac(x) = x - floor(x)
end

@testitem "B-spline interpolating regridder — polynomial-degree reproduction (esd-47z.3)" setup = [BSplineRegridSetup] tags = [:regridding, :bspline, :conformance] begin

    @testset "fixture loads (schema + structural)" begin
        @test isfile(_BS_FIXTURE)
        @test (ESS.load(_BS_FIXTURE); true)
        raw = _bs_raw()
        @test Set(string.(keys(raw["models"]))) ==
            Set(["BSplineRegridLinear1D", "BSplineRegridCubic1D", "BSplineRegridBilinear2D"])
    end

    # ---- DEGREE-1 (LINEAR) — the special case `bilinear` is the tensor product of.
    @testset "linear (degree 1): reproduces linear fields EXACTLY on staggered targets" begin
        P(x) = 2.0 + 3.0x                      # any linear field
        xsrc = Float64.(0:4)                    # 5 uniform nodes, h = 1
        Fsrc = P.(xsrc)
        # targets incl. a cell FACE at x=1.5 (s=0.5, the Arakawa-C staggering point)
        xt = [0.5, 1.5, 2.3, 3.0]
        base = _bs_linear_base.(xt); s = _bs_frac.(xt)
        got = _bs_eval_1d("BSplineRegridLinear1D", Fsrc, base, s)
        @test got ≈ P.(xt) rtol = 0 atol = 1.0e-13
        # the s=0.5 face value is the two-node average (weights {1/2,1/2})
        @test got[2] ≈ (Fsrc[2] + Fsrc[3]) / 2 atol = 1.0e-13
    end

    @testset "linear weights at a face are {1/2, 1/2}; partition unity for any s" begin
        @test _bs_weights_1d("BSplineRegridLinear1D", 5, 2, 0.5, 2) ≈ [0.5, 0.5] atol = 1.0e-13
        for s0 in (0.0, 0.27, 0.5, 0.83)
            @test sum(_bs_weights_1d("BSplineRegridLinear1D", 5, 2, s0, 2)) ≈ 1.0 atol = 1.0e-13
        end
    end

    # the degree is EXACTLY 1: a linear kernel does NOT reproduce a quadratic.
    @testset "linear does NOT reproduce a quadratic (the degree is meaningful)" begin
        Q(x) = x^2
        xsrc = Float64.(0:4); Fsrc = Q.(xsrc)
        xt = fill(1.5, _BS_NTGT)              # model bakes 4 targets; probe is [1]
        got = _bs_eval_1d("BSplineRegridLinear1D", Fsrc, _bs_linear_base.(xt), _bs_frac.(xt))
        @test !isapprox(got[1], Q(1.5); atol = 1.0e-6)     # 2.5 (chord) ≠ 2.25 (exact)
    end

    # ---- DEGREE-3 (CUBIC) — the representative higher-order B-spline stencil.
    @testset "cubic (degree 3): reproduces cubic fields EXACTLY on staggered targets" begin
        P(x) = 2.0 - 0.7x + 0.3x^2 - 0.1x^3
        xsrc = Float64.(0:6); Fsrc = P.(xsrc)       # 7 uniform nodes
        xt = [1.5, 2.5, 3.3, 4.0]                    # interior so base..base+3 ∈ 1..7
        base = _bs_cubic_base.(xt); s = _bs_frac.(xt)
        got = _bs_eval_1d("BSplineRegridCubic1D", Fsrc, base, s)
        @test got ≈ P.(xt) rtol = 0 atol = 1.0e-12
        # also reproduces every lower degree exactly (constant / linear / quadratic)
        for Plow in (x -> 5.0, x -> 1.0 + 2.0x, x -> 0.5x^2 - x + 3.0)
            Fl = Plow.(xsrc)
            @test _bs_eval_1d("BSplineRegridCubic1D", Fl, base, s) ≈ Plow.(xt) atol = 1.0e-12
        end
    end

    @testset "cubic face weights are the 4th-order {-1/16, 9/16, 9/16, -1/16}" begin
        w = _bs_weights_1d("BSplineRegridCubic1D", 7, 3, 0.5, 4)   # base=3 → nodes 3..6, central [4,5]
        @test w ≈ [-1/16, 9/16, 9/16, -1/16] atol = 1.0e-13
        @test sum(w) ≈ 1.0 atol = 1.0e-13                          # partition of unity
        @test w[1] ≈ w[4] && w[2] ≈ w[3]                           # symmetric about the face
        # partition of unity holds for any s (not just the face)
        for s0 in (0.0, 0.31, 0.5, 0.92)
            @test sum(_bs_weights_1d("BSplineRegridCubic1D", 7, 3, s0, 4)) ≈ 1.0 atol = 1.0e-13
        end
    end

    # ---- DEGREE-1 TENSOR PRODUCT (BILINEAR) — the staggered Arakawa-C case.
    @testset "bilinear (degree-1 tensor product): reproduces bilinear fields EXACTLY" begin
        P(x, y) = 1.0 + 2.0x - 0.5y + 0.3x * y       # a + bx + cy + d·xy
        xs = Float64.(0:3); ys = Float64.(0:3)
        Fsrc2 = [P(x, y) for x in xs, y in ys]        # 4×4 uniform grid
        xt = [0.5, 1.5, 2.0, 1.2]; yt = [0.5, 2.0, 1.3, 0.7]
        bx = _bs_linear_base.(xt); by = _bs_linear_base.(yt)
        sx = _bs_frac.(xt);        sy = _bs_frac.(yt)
        got = _bs_eval_2d("BSplineRegridBilinear2D", Fsrc2, bx, by, sx, sy)
        @test got ≈ [P(xt[j], yt[j]) for j in 1:4] rtol = 0 atol = 1.0e-12
    end

    @testset "bilinear IS the degree-1 special case: center = 4-pt avg = linear⊗linear" begin
        # F_src = the four corner basis impulses → recover the four bilinear weights
        # at the cell center (s_x = s_y = 0.5). They must equal {1/4,1/4,1/4,1/4} and
        # be the OUTER PRODUCT of the degree-1 linear weights {1/2,1/2}.
        wbil = Float64[]
        for (ax, ay) in ((0, 0), (1, 0), (0, 1), (1, 1))
            F = zeros(Float64, 2, 2); F[1 + ax, 1 + ay] = 1.0
            # BSplineRegridBilinear2D bakes tgt_locs size = 4; replicate the center
            # probe 4× (base_x=base_y=1, s_x=s_y=0.5) and read the first target.
            got = _bs_eval_2d("BSplineRegridBilinear2D", F,
                fill(1, _BS_NTGT), fill(1, _BS_NTGT), fill(0.5, _BS_NTGT), fill(0.5, _BS_NTGT))
            push!(wbil, got[1])
        end
        @test wbil ≈ [0.25, 0.25, 0.25, 0.25] atol = 1.0e-13
        @test sum(wbil) ≈ 1.0 atol = 1.0e-13                        # partition of unity
        wlin = [0.5, 0.5]                                            # degree-1 linear weights
        wouter = [wlin[1 + ax] * wlin[1 + ay] for (ax, ay) in ((0, 0), (1, 0), (0, 1), (1, 1))]
        @test wbil ≈ wouter atol = 1.0e-13                          # bilinear = linear ⊗ linear
    end
end
