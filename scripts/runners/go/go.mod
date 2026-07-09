// Module for the Go conformance runner (scripts/runners/run-go.sh copies this
// directory to a scratch build dir and rewrites the replace directive to the
// resolved $ESS_ROOT before `go run`, so the sibling-checkout default below is
// only the local-convenience path — scripts/ess-locate.sh contract).
module github.com/EarthSciML/earthscidiscretizations/scripts/runners/go-runner

go 1.21

require github.com/EarthSciML/EarthSciAST/pkg/earthsci-ast-go v0.0.0

replace github.com/EarthSciML/EarthSciAST/pkg/earthsci-ast-go => ../../../../EarthSciAST/pkg/earthsci-ast-go
