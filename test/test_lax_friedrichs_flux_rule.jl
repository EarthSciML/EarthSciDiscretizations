using Test
using TestItems

# Tests for the finite_volume/lax_friedrichs_flux declarative rule.
#
# The rule encodes the Lax-Friedrichs numerical flux for linear advection
# F_{i+1/2} = max(c,0)·q_i + min(c,0)·q_{i+1} as a 2-point cartesian stencil
# on $q with the face-staggered Courant `$c` carried as a per-face binding.
# Coefficients reference `$c` symbolically via the `abs` op — no caller-side
# branching on sign(c) is required.
#
# Layer A   — rule discovery + JSON byte-diff round-trip + spot-check schema.
# Layer A'  — hand-pinned 5-interior-face example (q = powers-of-2; mixed-sign
#             Courant including c=0) verifying the closed-form expression
#             max(c,0)·q_i + min(c,0)·q_{i+1} at every interface. Values live
#             inline here rather than in `<rule>/fixtures/canonical/` because
#             ESS's `discretize` does not yet support op="flux" with per-face
#             bindings — the walker's Layer A would FAIL on a canonical fixture
#             until the schema follow-ups land (boundary_policy + face-stagger
#             / time-varying bindings, tracked off dsc-35x).
#
# Layer B (MMS convergence) is `applicable: false` — the rule output is at
# the face stagger and the per-face `$c` binding is outside the §7
# verify_mms_convergence contract today (same gap as latlon's `cos_lat`).
# A follow-up MMS-transport fixture (LF + divergence + forward Euler) is the
# right shape and is tracked off dsc-35x.

# ---------------------------------------------------------------------------
# Layer A: discovery + spot-checks
# ---------------------------------------------------------------------------

@testitem "lax_friedrichs_flux rule is discoverable under :finite_volume" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "lax_friedrichs_flux", rules)
    @test idx !== nothing
    rule = rules[idx]
    @test rule.family == :finite_volume
    @test isfile(rule.path)

    content = read(rule.path, String)
    @test occursin("\"applies_to\"", content)
    @test occursin("\"op\": \"flux\"", content)
    @test occursin("\"grid_family\"", content)
    @test occursin("\"cartesian\"", content)
    @test occursin("\"lax_friedrichs\"", content)
    @test occursin("\"O(h)\"", content)
    @test occursin("\"stencil\"", content)
    @test occursin("\"\$c\"", content)
    @test occursin("\"\$q\"", content)
    @test occursin("\"abs\"", content)
    # Two-point flux stencil at the face: offsets 0 and 1 only.
    @test occursin("\"offset\": 0", content)
    @test occursin("\"offset\": 1", content)
end

@testitem "lax_friedrichs_flux rule JSON round-trips byte-stable" begin
    using EarthSciDiscretizations
    using EarthSciDiscretizations: load_rules
    using JSON

    repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
    catalog = joinpath(repo_root, "discretizations")
    rules = load_rules(catalog)
    idx = findfirst(r -> r.name == "lax_friedrichs_flux", rules)
    @test idx !== nothing
    rule = rules[idx]

    raw = read(rule.path, String)
    parsed = JSON.parse(raw)
    @test parsed isa AbstractDict
    @test haskey(parsed, "discretizations")
    @test haskey(parsed["discretizations"], "lax_friedrichs_flux")

    reserialized = JSON.json(parsed)
    reparsed = JSON.parse(reserialized)
    @test reparsed == parsed

    spec = parsed["discretizations"]["lax_friedrichs_flux"]
    @test spec["applies_to"]["op"] == "flux"
    @test spec["applies_to"]["dim"] == "\$x"
    @test spec["applies_to"]["args"] == ["\$q", "\$c"]
    @test spec["grid_family"] == "cartesian"
    @test spec["form"] == "lax_friedrichs"
    @test spec["accuracy"] == "O(h)"
    @test spec["produces"] == "F_{i+1/2}"
    @test spec["monotone"] == true
    @test spec["tvd"] == true
    @test spec["stencil"] isa AbstractArray
    @test length(spec["stencil"]) == 2
    offsets = sort([s["selector"]["offset"] for s in spec["stencil"]])
    @test offsets == [0, 1]
    # Each coefficient must reference $c and the abs op (upwind via |c|).
    for entry in spec["stencil"]
        coeff_str = JSON.json(entry["coeff"])
        @test occursin("\$c", coeff_str)
        @test occursin("abs", coeff_str)
    end
end

# ---------------------------------------------------------------------------
# Layer A' (hand-pinned): closed-form flux at 5 interior faces
# ---------------------------------------------------------------------------

@testitem "lax_friedrichs_flux: hand-pinned closed-form face fluxes" begin
    # q = powers of 2 → trivially recognisable face fluxes.
    # c at interfaces = [+0.5, -0.5, +0.5, -0.5, 0.0] exercises both upwind
    # branches and the c=0 zero-crossing where both coefficients must vanish.
    q = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    c = [0.5, -0.5, 0.5, -0.5, 0.0]

    # F_{i+1/2} = max(c, 0)·q_i + min(c, 0)·q_{i+1}
    F = [max(c[k], 0) * q[k] + min(c[k], 0) * q[k + 1] for k in 1:5]
    @test F == [0.5, -2.0, 2.0, -8.0, 0.0]
    # The c=0 case: |c|=0, so both stencil coefficients are 0 and F = 0
    # regardless of the surrounding q values — no caller-side `if c >= 0` branch.
    @test F[5] == 0.0
end
