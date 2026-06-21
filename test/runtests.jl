# Test entrypoint, run under the package's test target (`Pkg.test`), where the
# [extras]+[targets].test deps — including the test-only ODE solver
# OrdinaryDiffEqDefault — are on the load path.
#
#   Pkg.test()                              → TestItemRunner unit suite (default)
#   Pkg.test(; test_args=["integration"])   → Layer-C integration suite
#
# The integration suite is kept a separate Pkg.test invocation (its own CI job)
# because the ODE solve carries a heavy precompile cost the unit suite shouldn't
# pay. Routing it through the test target is what makes OrdinaryDiffEqDefault —
# a test-only dependency, not a runtime dep of the library — importable by the
# integration cases.
#
# `using TestItemRunner` MUST stay a top-level statement (not folded into the
# `if` below): a macro call is expanded together with its enclosing top-level
# expression, so `@run_package_tests` inside the `if` is resolved before any
# `using` in that same block has run — importing the package here, as a separate
# preceding top-level statement, is what makes the macro defined at that point.
using TestItemRunner

if "integration" in ARGS
    ENV["ESD_RUN_INTEGRATION"] = "1"
    include(joinpath(@__DIR__, "run_integration_tests.jl"))
else
    @run_package_tests
end
