// Module for the Go conformance runner (scripts/runners/run-go.sh copies this
// directory to a scratch build dir and rewrites the replace directive to the
// resolved $ESS_ROOT before `go run`, so the sibling-checkout default below is
// only the local-convenience path — scripts/ess-locate.sh contract).
module github.com/EarthSciML/earthscidiscretizations/scripts/runners/go-runner

go 1.21

require github.com/ctessum/EarthSciSerialization/packages/esm-format-go v0.0.0

replace github.com/ctessum/EarthSciSerialization/packages/esm-format-go => ../../../../EarthSciSerialization/packages/esm-format-go
