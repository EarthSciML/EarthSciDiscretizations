# Reprojection LAMBERT CONFORMAL CONIC — Julia evaluator conformance.
#
# ESD owns this declarative reprojection model (lambert_conformal.esm, co-located
# here); the EarthSciSerialization engine (an ESD dependency) EVALUATES it. This
# is the second member of the reprojection/ family (RFC pure-io-data-loaders §5.1
# + Appendix B, bead esd-47z.2) and the first with a NON-trivial forward/inverse
# pair — the spherical Lambert Conformal Conic projection, authored as a plain
# declarative .esm over the EXISTING operator set (sin/cos/tan/atan/atan2/log/^/
# sqrt/sign + the built-in `pi` op) with NO new rule-schema, scaffold, validation,
# resolver, or per-binding code. It mirrors reprojection/longlat.esm exactly: the
# projection is PARAMETERIZED by {lat_1, lat_2, lat_0, lon_0, R} (scalar
# `parameter`s driven via parameter_overrides), and the FORWARD (lon,lat → x,y)
# and INVERSE (x,y → lon,lat) maps are SCALAR `observed`s — a coordinate transform
# is a point-wise scalar formula (a binding applies the kernel over a coordinate
# array at the runtime layer; the tree-walk evaluator admits array-shaped
# observeds only under an active geometry/value-invention kernel, which a pure
# transform has not).
#
#   FORWARD (Appendix B), inputs in degrees, x/y in metres:
#     n  = log(cos φ₁/cos φ₂) / log(tan(π/4+φ₂/2)/tan(π/4+φ₁/2))
#     F  = cos φ₁·tan(π/4+φ₁/2)^n / n
#     ρ  = R·F/tan(π/4+φ/2)^n      ρ₀ = R·F/tan(π/4+φ₀/2)^n
#     x  = ρ·sin(n(λ−λ₀))          y  = ρ₀ − ρ·cos(n(λ−λ₀))
#   INVERSE (Appendix B), closed form:
#     ρ = sign(n)·sqrt(x²+(ρ₀−y)²)   θ = atan2(x, ρ₀−y)
#     λ = λ₀ + θ/n                   φ = 2·atan((R·F/ρ)^(1/n)) − π/2
#
# The inverse observeds compose on the forward observeds (rho0_my/rho_inv/
# theta_inv read proj_x/proj_y), so inverse∘forward is the identity BY
# CONSTRUCTION. Each map is read off through build_evaluator via the regridder's
# zero-IC state idiom: a state whose D-equation is d(state)/dt = <observed> from a
# zero IC, so du = f!(u0) IS the map value (the narrow_phase_area / longlat
# precedent). forward_x/forward_y surface the forward map; roundtrip_lon/
# roundtrip_lat surface the round-trip.
#
# ACCEPTANCE (esd-47z.2): the rule is the spherical LCC through the existing
# engine with ZERO new schema/validation/resolver code, and —
#   • The forward map matches an INDEPENDENT proj4/Snyder reference (computed
#     from the published spherical-LCC closed form, identical to PROJ
#     `+proj=lcc +R=…` for the sphere) at representative CONUS points, for BOTH
#     the WRF and NEI2016 parameter sets — two real datasets exercise the ONE
#     parameterized rule.
#   • The projection origin (lon_0, lat_0) maps to the false origin (0, 0).
#   • inverse∘forward = identity to tolerance over the domain (both param sets).
#   • The parameters are LOAD-BEARING: WRF and NEI2016 produce DIFFERENT outputs
#     for the same coordinate (the declarative-or-fail proof the .esm binds a
#     genuinely parameterized transform, not a hard-coded constant).
# Per CONFORMANCE_SPEC §5.8 the other bindings schema-/structurally-validate the
# document while Julia evaluates it here.

@testsnippet LambertConformalSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach it
    # through the dependency rather than adding a redundant top-level dep — and so
    # the raw document is the exact type build_evaluator expects.
    const JSON3 = ESS.JSON3

    const _LCC_FIXTURE = joinpath(@__DIR__, "lambert_conformal.esm")

    # The two real spherical-LCC parameter sets from Appendix A
    # ({lat_1, lat_2, lat_0, lon_0, R}); both flow through the ONE rule.
    const _WRF     = (lat_1 = 30.0, lat_2 = 60.0, lat_0 = 38.999996, lon_0 = -97.0, R = 6370000.0)
    const _NEI2016 = (lat_1 = 33.0, lat_2 = 45.0, lat_0 = 40.0,      lon_0 = -97.0, R = 6370997.0)

    # INDEPENDENT forward reference: (lon, lat, x, y) tuples computed from the
    # published spherical-LCC closed form (Snyder 1987 eqs. 15-1..15-4 = PROJ
    # `+proj=lcc +R=…` for the sphere) in NumPy, OUTSIDE the .esm AST, so they
    # cross-check the declarative document against an independent implementation.
    # Representative CONUS points, well away from the singular far pole.
    const _GOLD_WRF = [
        (-97.0, 39.0, 0.0, 0.43226828519254923),
        (-120.0, 35.0, -2028208.5469169612, -140947.2885132283),
        (-75.0, 40.0, 1795192.9893785133, 356152.97854254674),
        (-90.0, 45.0, 530758.0790195653, 668954.7596800169),
        (-100.0, 25.0, -309853.27839640283, -1541532.7554675443),
        (-80.0, 48.0, 1213070.0930733178, 1097145.5185974706),
        (-110.0, 31.0, -1228247.9439207606, -773888.551179463),
        (-104.0, 42.5, -554207.2432213794, 401411.8990695188),
    ]
    const _GOLD_NEI2016 = [
        (-97.0, 39.0, 0.0, -110589.55965320487),
        (-120.0, 35.0, -2066463.0536515424, -290403.0815928243),
        (-75.0, 40.0, 1845776.353572973, 224515.92418459523),
        (-90.0, 45.0, 549842.4483977047, 575299.4792623082),
        (-100.0, 25.0, -309377.0930418335, -1668877.3176074345),
        (-80.0, 48.0, 1266623.2685378059, 1007635.9929630421),
        (-110.0, 31.0, -1239963.1270495378, -909332.8255445026),
        (-104.0, 42.5, -571190.8151147809, 298694.874228091),
    ]

    # Drive the LCC model through build_evaluator at a single coordinate
    # (lon, lat) and parameter set, returning (forward_x, forward_y,
    # roundtrip_lon, roundtrip_lat). Each surfacing state is a constant-RHS
    # D-equation from a zero IC, so du = f!(u0) IS the surfaced map value.
    function _lcc_eval(lon::Float64, lat::Float64, p)
        raw = JSON3.read(read(_LCC_FIXTURE, String))
        ics = Dict(
            "forward_x" => 0.0, "forward_y" => 0.0,
            "roundtrip_lon" => 0.0, "roundtrip_lat" => 0.0,
        )
        f!, u0, pp, _, vmap = ESS.build_evaluator(
            raw; model_name = "LambertConformal",
            initial_conditions = ics,
            parameter_overrides = Dict(
                "lon" => lon, "lat" => lat,
                "lat_1" => p.lat_1, "lat_2" => p.lat_2, "lat_0" => p.lat_0,
                "lon_0" => p.lon_0, "R" => p.R,
            ),
        )
        du = similar(u0)
        f!(du, u0, pp, 0.0)
        return (
            du[vmap["forward_x"]], du[vmap["forward_y"]],
            du[vmap["roundtrip_lon"]], du[vmap["roundtrip_lat"]],
        )
    end
end

@testitem "reprojection lambert_conformal forward/inverse (esd-47z.2)" setup = [LambertConformalSetup] tags = [:reprojection, :lambert_conformal, :conformance] begin

    @testset "fixture loads (schema + structural)" begin
        @test isfile(_LCC_FIXTURE)
        @test (ESS.load(_LCC_FIXTURE); true)
    end

    # The FORWARD map matches the independent proj4/Snyder reference for BOTH
    # parameter sets. The .esm AST and the NumPy reference both evaluate the
    # SAME closed form, differing only by cross-libm rounding over a ~10-op
    # chain, so they agree to ≪ mm over the ~2000 km CONUS domain; a real
    # transcription error in the AST would be metres-to-km off, far outside the
    # tolerance. This is the two-real-datasets-one-rule acceptance.
    @testset "forward matches proj4 reference (WRF + NEI2016)" begin
        for (p, gold) in ((_WRF, _GOLD_WRF), (_NEI2016, _GOLD_NEI2016))
            for (lon, lat, gx, gy) in gold
                fx, fy, _, _ = _lcc_eval(lon, lat, p)
                @test isapprox(fx, gx; rtol = 1.0e-7, atol = 1.0e-4)
                @test isapprox(fy, gy; rtol = 1.0e-7, atol = 1.0e-4)
            end
        end
    end

    # The projection origin maps to the false origin (0, 0): at (lon_0, lat_0),
    # θ = n(λ−λ₀) = 0 ⇒ x = ρ·sin 0 = 0, and ρ = ρ₀ ⇒ y = ρ₀ − ρ = 0. A
    # structural invariant of a correctly-wired LCC, independent of the
    # reference points.
    @testset "origin (lon_0, lat_0) → (0, 0)" begin
        for p in (_WRF, _NEI2016)
            fx, fy, _, _ = _lcc_eval(p.lon_0, p.lat_0, p)
            @test isapprox(fx, 0.0; atol = 1.0e-6)
            @test isapprox(fy, 0.0; atol = 1.0e-6)
        end
    end

    # ACCEPTANCE — inverse∘forward = identity to tolerance over the domain, for
    # BOTH parameter sets. The inverse observeds compose on the forward
    # observeds, so the round-trip restores the input by construction; this
    # gate confirms it numerically across the CONUS domain.
    @testset "inverse∘forward round-trip = identity (WRF + NEI2016)" begin
        for (p, gold) in ((_WRF, _GOLD_WRF), (_NEI2016, _GOLD_NEI2016))
            for (lon, lat, _gx, _gy) in gold
                _, _, rlon, rlat = _lcc_eval(lon, lat, p)
                @test isapprox(rlon, lon; rtol = 0, atol = 1.0e-9)
                @test isapprox(rlat, lat; rtol = 0, atol = 1.0e-9)
            end
        end
    end

    # The PARAMETERS are LOAD-BEARING (the declarative-or-fail proof that the
    # .esm binds a PARAMETERIZED transform, not a constant): the WRF and NEI2016
    # parameter sets — which share lon_0 = −97 but differ in standard parallels,
    # origin latitude, and radius — produce DIFFERENT projected coordinates for
    # the same geographic point. The cone constant alone differs (n ≈ 0.7156 vs
    # 0.6305), so no output can coincide by accident.
    @testset "parameters are load-bearing (WRF ≠ NEI2016)" begin
        for (lon, lat, _gx, _gy) in _GOLD_WRF
            (lon == -97.0) && continue   # x ≡ 0 on the shared central meridian
            fxW, fyW, _, _ = _lcc_eval(lon, lat, _WRF)
            fxN, fyN, _, _ = _lcc_eval(lon, lat, _NEI2016)
            @test fxW != fxN
            @test fyW != fyN
        end
    end
end
