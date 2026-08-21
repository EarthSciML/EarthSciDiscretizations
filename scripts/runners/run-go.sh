#!/usr/bin/env bash
# run-go.sh — the official Go conformance runner for EarthSciDiscretizations.
#
# Thin shell over scripts/runners/go/main.go (the §4.2 CLI, ast category),
# which drives the official earthsci-ast-go raw §9.7 pipeline
# (esm.ResolveAndLowerReferencePreserving). This wrapper only resolves the ESS checkout and
# builds the runner in a scratch dir so the committed tree stays clean:
# the runner module's `replace` directive is rewritten to the resolved
# $ESS_ROOT (scripts/ess-locate.sh contract), then `go run` forwards the
# §4.2 arguments unchanged.
#
# Usage: scripts/runners/run-go.sh --output-dir <path> [--categories ast]
#                                  [--files <manifest.json>[,…]] [--verbose]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

BUILD="$(mktemp -d "${TMPDIR:-/tmp}/esd-go-runner.XXXXXX")"
trap 'rm -rf "$BUILD"' EXIT

cp "$REPO/scripts/runners/go/main.go" "$REPO/scripts/runners/go/go.mod" "$BUILD/"
(
  cd "$BUILD"
  go mod edit -replace \
    "github.com/EarthSciML/EarthSciAST/pkg/earthsci-ast-go=$ESS_ROOT/pkg/earthsci-ast-go"
  GOFLAGS=-mod=mod ESD_ROOT="$REPO" go run . "$@"
)
