# Steady-state / elliptic operator — DECLARATIVE-OR-FAIL verdict (esd-6g4.5, G5).
#
# Pins the verdict in `discretizations/STEADY_STATE_ELLIPTIC_VERDICT.md`:
# the discrete elliptic OPERATOR/RESIDUAL is expressible over the EXISTING engine
# with ZERO new rules. The existing 2nd-order Laplacian rule is solve-mode-
# agnostic — placed on the RHS of an ALGEBRAIC equation (`f = ∇²u`, i.e. the
# Poisson residual `∇²u − f = 0`), `discretize` rewrites it to the discrete
# stencil and classifies the system as a DAE; the *solve* is left to an external
# DAE/nonlinear assembler (RFC §2 non-goal — never bundled here).
#
# These cases drive the canonical pipeline (`EarthSciSerialization.discretize`)
# on the Laplacian rule's own canonical fixture, repurposed as a steady Poisson
# model. No new rule, fixture, or imperative operator is involved.

@testsnippet EllipticResidualSetup begin
    import EarthSciSerialization
    using JSON

    const LAPLACIAN_CANONICAL_INPUT = joinpath(
        @__DIR__, "..", "discretizations", "finite_difference",
        "laplacian_2nd_uniform_cartesian", "fixtures", "canonical", "input.esm"
    )

    # Steady Poisson model from the Laplacian canonical template: add a forcing
    # field `f` and replace `∂u/∂t = ∇²u` with the algebraic `f = ∇²u`
    # (operator on the RHS — the form `discretize` rewrites; see verdict §4).
    function build_poisson_steady()
        input = JSON.parsefile(LAPLACIAN_CANONICAL_INPUT)
        input["metadata"]["name"] = "elliptic_residual_steady_poisson"
        input["models"]["M"]["variables"]["f"] = Dict(
            "default" => 0.0, "type" => "parameter", "units" => "1",
            "shape" => ["x", "y"], "location" => "cell_center"
        )
        input["models"]["M"]["equations"] = [
            Dict(
                "lhs" => "f",
                "rhs" => Dict("op" => "laplacian", "args" => ["u"])
            ),
        ]
        return input
    end

    # The matched operator subtree is gone from a rewritten equation expression.
    eq_json(eq) = JSON.json(eq)
end

@testitem "elliptic residual: steady Poisson discretizes to a DAE via the existing Laplacian rule" setup = [EllipticResidualSetup] tags = [:elliptic, :dae, :discretization] begin
    out = EarthSciSerialization.discretize(build_poisson_steady(); dae_support = true)

    # Steady/elliptic ⇒ algebraic ⇒ DAE classification.
    @test out["metadata"]["system_class"] == "dae"
    @test out["metadata"]["dae_info"]["algebraic_equation_count"] == 1
    @test out["metadata"]["dae_info"]["per_model"]["M"] == 1

    eq = out["models"]["M"]["equations"][1]
    s = eq_json(eq)
    # The Laplacian operator was rewritten to the discrete stencil (residual form
    # f = ∇²ₕu): no `laplacian` op remains, `index` gathers are materialized.
    @test !occursin("laplacian", s)
    @test occursin("\"index\"", s)
    # Residual shape preserved: forcing on the LHS, stencil sum on the RHS.
    @test eq["lhs"] == "f"
    @test eq["rhs"] isa AbstractDict && eq["rhs"]["op"] == "+"
end

@testitem "elliptic residual: steady Poisson is rejected when DAE support is off" setup = [EllipticResidualSetup] tags = [:elliptic, :dae, :discretization] begin
    # A binding that hands off only to an ODE integrator cannot accept the
    # elliptic (algebraic) system — RFC §12 contract.
    err = try
        EarthSciSerialization.discretize(build_poisson_steady(); dae_support = false)
        nothing
    catch e
        e
    end
    @test err isa EarthSciSerialization.RuleEngineError
    @test err !== nothing && err.code == "E_NO_DAE_SUPPORT"
end

@testitem "elliptic residual: the Laplacian rule is solve-mode-agnostic (ODE form stays ODE)" setup = [EllipticResidualSetup] tags = [:elliptic, :dae, :discretization] begin
    # The unmodified template (∂u/∂t = ∇²u) discretizes through the SAME rule but
    # classifies as an ODE — the rule encodes no solve mode; the equation does.
    out = EarthSciSerialization.discretize(
        JSON.parsefile(LAPLACIAN_CANONICAL_INPUT); dae_support = true
    )
    @test out["metadata"]["system_class"] == "ode"
    @test out["metadata"]["dae_info"]["algebraic_equation_count"] == 0
    s = eq_json(out["models"]["M"]["equations"][1])
    @test !occursin("laplacian", s)   # still rewritten to the discrete stencil
    @test occursin("\"index\"", s)
end
