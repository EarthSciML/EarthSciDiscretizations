#!/usr/bin/env bash
# test-conformance.sh — cross-binding conformance orchestrator.
#
# Mirrors the shape of EarthSciSerialization's scripts/test-conformance.sh:
# probe each binding's availability, run the available runners for the
# requested categories into conformance-results/<binding>/, then gate every
# output against the committed goldens (and against the reference binding)
# with scripts/compare-outputs.py.
#
# Only the Julia reference runner exists today; the per-binding probe/run
# tables below are where run-python.py, run-rust, run-go, run-typescript drop
# in as the ESS §9.7 ports land (AGENTS.md "Phasing").
#
# Usage:
#   scripts/test-conformance.sh [--categories ast,simulation,…]
#                               [--bindings julia[,…]]
#                               [--output-dir conformance-results] [--verbose]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

CATEGORIES="ast,simulation,convergence,regridding,reprojection"
BINDINGS="julia"
OUTPUT_DIR="$REPO/conformance-results"
VERBOSE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --categories) CATEGORIES="$2"; shift 2 ;;
    --bindings) BINDINGS="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --verbose) VERBOSE="--verbose"; shift ;;
    *) echo "error: unknown argument '$1'" >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# Per-binding availability probes + runner invocations (§4.2 CLI for each).
# ---------------------------------------------------------------------------
probe_julia() {
  command -v julia >/dev/null 2>&1 &&
    [[ -f "$ESS_ROOT/packages/EarthSciSerialization.jl/Project.toml" ]]
}
run_julia() { # $1 = output dir
  julia "$REPO/scripts/runners/run-julia.jl" \
    --output-dir "$1" --categories "$CATEGORIES" $VERBOSE
}

# Future bindings land here as the ESS §9.7 ports arrive:
#   probe_python / run_python  -> scripts/runners/run-python.py  (earthsci_toolkit)
#   probe_rust   / run_rust    -> scripts/runners/run-rust       (earthsci-toolkit-rs)
#   probe_go     / run_go      -> scripts/runners/run-go         (esm-format-go)
#   probe_typescript / run_typescript -> scripts/runners/run-typescript.js

ran=()
IFS=',' read -ra requested <<<"$BINDINGS"
for binding in "${requested[@]}"; do
  if ! declare -F "probe_$binding" >/dev/null; then
    echo "-- $binding: no runner registered yet (skipping)"
    continue
  fi
  if ! "probe_$binding"; then
    echo "-- $binding: toolchain unavailable (skipping)"
    continue
  fi
  echo "== running $binding runner (categories: $CATEGORIES)"
  mkdir -p "$OUTPUT_DIR/$binding"
  "run_$binding" "$OUTPUT_DIR/$binding"
  ran+=("$binding")
done

if [[ ${#ran[@]} -eq 0 ]]; then
  echo "error: no requested binding is available" >&2
  exit 1
fi

echo "== comparing outputs"
RAN_CSV="$(IFS=,; echo "${ran[*]}")"
python3 "$REPO/scripts/compare-outputs.py" \
  --results-root "$OUTPUT_DIR" --bindings "$RAN_CSV" --categories "$CATEGORIES" $VERBOSE
