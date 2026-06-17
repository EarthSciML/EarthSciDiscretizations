# Layer-C integration driver — invoked by the Integration CI job.
#
# Runs walk_esd_tests with ESD_RUN_INTEGRATION=1 set in the environment so
# every rule with integration fixtures is exercised end-to-end: ESM+GDD →
# build_ode_problem → ODE solve → tolerance check.  Exits non-zero if any
# Layer-C case fails (LAYER_FAIL).  LAYER_SKIP is acceptable (stub cases or
# rules without integration fixtures declared).

include(joinpath(@__DIR__, "walk_esd_tests.jl"))
using .WalkESDTests
using EarthSciDiscretizations

if get(ENV, "ESD_RUN_INTEGRATION", "0") != "1"
    println(stderr, "ERROR: ESD_RUN_INTEGRATION must be set to \"1\" to run this script.")
    println(stderr, "Usage: ESD_RUN_INTEGRATION=1 julia --project=. test/run_integration_tests.jl")
    exit(1)
end

repo_root = dirname(dirname(pathof(EarthSciDiscretizations)))
catalog = joinpath(repo_root, "discretizations")

println("Running Layer-C integration suite (catalog: $(catalog))")
results = WalkESDTests.walk_esd_tests(; catalog = catalog)

n_fail = count(r -> r.layer_c.outcome == WalkESDTests.LAYER_FAIL, results)
n_pass = count(r -> r.layer_c.outcome == WalkESDTests.LAYER_PASS, results)
n_skip = count(r -> r.layer_c.outcome == WalkESDTests.LAYER_SKIP, results)

println("\nLayer-C: $n_pass pass, $n_skip skip, $n_fail fail")
for r in filter(r -> r.layer_c.outcome == WalkESDTests.LAYER_FAIL, results)
    println("  FAIL [$(r.family)/$(r.name)]: $(r.layer_c.reason)")
end

if n_fail > 0
    println("\nERROR: $n_fail Layer-C case(s) failed — see output above.")
    exit(1)
end
println("Layer-C integration: all cases pass or skip.")
