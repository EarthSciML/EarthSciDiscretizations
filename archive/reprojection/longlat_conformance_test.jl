# Reprojection LONGLAT — Julia evaluator conformance (round-trip identity).
#
# ESD owns this declarative reprojection model (longlat.esm, co-located here);
# the EarthSciSerialization engine (an ESD dependency) EVALUATES it. This is the
# convention-establishing member of the reprojection/ family (RFC
# pure-io-data-loaders §5.1, bead esd-47z.1): the geographic WGS84 long/lat
# "identity" projection, authored as a plain declarative .esm over the EXISTING
# operator set — NO new rule-schema, scaffold, validation, resolver, or
# per-binding code. It mirrors regridding/conservative_regrid_overlap_join.esm:
# the projection's configuration (`lon_0`, the central meridian) is a scalar
# `parameter` driven via parameter_overrides, and the FORWARD
# (geographic→projected) and INVERSE (projected→geographic) maps are SCALAR
# `observed` expressions — a coordinate transform is a point-wise scalar formula
# (the regridder's `narrow_phase_area`/`A_ij_rep` scalar pair, not its
# array-shaped `A_ij`; the tree-walk evaluator only admits array-shaped
# parameters/observeds under an active geometry/value-invention kernel).
#
#   FORWARD:  x = lon − lon_0,   y = lat
#   INVERSE:  lon = x + lon_0,   lat = y
#
# The inverse observeds compose on the forward observeds (geo_lon reads proj_x),
# so inverse∘forward = (lon − lon_0) + lon_0 = lon is the identity BY
# CONSTRUCTION. Each map is read off through build_evaluator via the regridder's
# zero-IC state idiom: a state whose D-equation is d(state)/dt = <observed> from
# a zero IC, so du = f!(u0) at the zero IC IS the map's value (the
# narrow_phase_area precedent). forward_x/forward_y surface the forward map;
# roundtrip_lon/roundtrip_lat surface the round-trip.
#
# ACCEPTANCE (esd-47z.1): the round-trip is the identity through the existing
# engine with ZERO new schema/validation/resolver code —
#   • lon_0 = 0 (WGS84 default): forward IS the geographic identity and the
#     round-trip equals the input EXACTLY (the bead's "forward = inverse =
#     identity").
#   • lon_0 ≠ 0: the parameter is LOAD-BEARING (forward x = lon − lon_0 ≠ lon),
#     proving the format binds a PARAMETERIZED transform (declarative-or-fail),
#     while the inverse∘forward round-trip stays the identity.
# Per CONFORMANCE_SPEC §5.8 the other bindings schema-/structurally-validate the
# document while Julia evaluates it here.

@testsnippet LongLatSetup begin
    using Test
    import EarthSciSerialization as ESS
    # JSON3 is the engine's serialization library (a hard dep of ESS); reach it
    # through the dependency rather than adding a redundant top-level dep — and so
    # the raw document is the exact type build_evaluator expects.
    const JSON3 = ESS.JSON3

    const _LL_FIXTURE = joinpath(@__DIR__, "longlat.esm")

    # Representative geographic coordinates (degrees), chosen exactly
    # representable so the shifted forward map and the round-trip are byte-exact
    # at integer lon_0. The kernel is a scalar point map; a binding applies it
    # over a coordinate array, so the gate sweeps these points one at a time.
    const _COORDS = [
        (-120.0, -45.0), (0.0, 0.0), (75.0, 60.0), (180.0, -90.0), (-179.0, 89.0),
    ]

    # Drive the longlat model through build_evaluator at a single coordinate
    # (lon, lat) and central meridian lon_0, returning
    # (forward_x, forward_y, roundtrip_lon, roundtrip_lat). Each surfacing state
    # is a constant-RHS D-equation from a zero IC, so du = f!(u0) IS the surfaced
    # map value.
    function _ll_eval(lon::Float64, lat::Float64, lon_0::Float64)
        raw = JSON3.read(read(_LL_FIXTURE, String))
        ics = Dict(
            "forward_x" => 0.0, "forward_y" => 0.0,
            "roundtrip_lon" => 0.0, "roundtrip_lat" => 0.0,
        )
        f!, u0, p, _, vmap = ESS.build_evaluator(
            raw; model_name = "LongLat",
            initial_conditions = ics,
            parameter_overrides = Dict("lon" => lon, "lat" => lat, "lon_0" => lon_0),
        )
        du = similar(u0)
        f!(du, u0, p, 0.0)
        return (
            du[vmap["forward_x"]], du[vmap["forward_y"]],
            du[vmap["roundtrip_lon"]], du[vmap["roundtrip_lat"]],
        )
    end
end

@testitem "reprojection longlat round-trip identity (esd-47z.1)" setup = [LongLatSetup] tags = [:reprojection, :longlat, :conformance] begin

    @testset "fixture loads (schema + structural)" begin
        @test isfile(_LL_FIXTURE)
        @test (ESS.load(_LL_FIXTURE); true)
    end

    # lon_0 = 0 (WGS84 default): the forward map IS the geographic identity and
    # the inverse∘forward round-trip equals the input EXACTLY. The bead's
    # "forward = inverse = identity" — exact (subtracting/adding 0.0 is a
    # floating-point no-op) at every representative coordinate.
    @testset "lon_0 = 0: forward AND round-trip are the geographic identity (exact)" begin
        for (lon, lat) in _COORDS
            fx, fy, rlon, rlat = _ll_eval(lon, lat, 0.0)
            @test fx == lon          # forward x is the identity at lon_0 = 0
            @test fy == lat          # forward y is always the identity in latitude
            @test rlon == lon        # ACCEPTANCE: round-trip longitude == input
            @test rlat == lat        # ACCEPTANCE: round-trip latitude  == input
        end
    end

    # The PARAMETER is LOAD-BEARING (the declarative-or-fail proof that the .esm
    # format binds a PARAMETERIZED transform, not a constant): a non-zero central
    # meridian shifts the forward longitude, x = lon − lon_0 ≠ lon, while latitude
    # is untouched. Integer-valued coords and lon_0 ⇒ exact.
    @testset "lon_0 ≠ 0: parameter shifts the forward map (load-bearing)" begin
        lon_0 = 10.0
        for (lon, lat) in _COORDS
            fx, fy, _, _ = _ll_eval(lon, lat, lon_0)
            @test fx == lon - lon_0      # the shift is real — not folded to identity
            @test fx != lon              # lon_0 genuinely moves the forward output
            @test fy == lat              # latitude is the identity regardless of lon_0
        end
    end

    # ACCEPTANCE — the round-trip stays the identity for a non-zero parameter:
    # inverse∘forward = (lon − lon_0) + lon_0 = lon, restoring the input. This is
    # the round-trip a real (non-identity) projection in the family must satisfy.
    @testset "lon_0 ≠ 0: round-trip recovers the input (identity preserved)" begin
        for lon_0 in (10.0, -73.5, 180.0)
            for (lon, lat) in _COORDS
                _, _, rlon, rlat = _ll_eval(lon, lat, lon_0)
                @test isapprox(rlon, lon; rtol = 0, atol = 1.0e-12)
                @test isapprox(rlat, lat; rtol = 0, atol = 1.0e-12)
            end
        end
    end
end
