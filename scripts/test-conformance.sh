#!/usr/bin/env bash
# test-conformance.sh — cross-binding conformance orchestrator.
#
# Mirrors the shape of EarthSciSerialization's scripts/test-conformance.sh:
# probe each binding's availability, run the available runners for the
# requested categories into conformance-results/<binding>/, then gate every
# output against the committed goldens (and against the reference binding)
# with scripts/compare-outputs.py.
#
# All five §9.7 binding runners are registered below. Each runner skips the
# categories outside its binding's scope (CONFORMANCE_SPEC.md §5.9; the
# manifests' scope_excluded / blocked-upstream notes carry the reasons), so
# requesting every category against every binding is always safe.
#
# Usage:
#   scripts/test-conformance.sh [--categories ast,simulation,…]
#                               [--bindings julia,python,rust,typescript,go]
#                               [--output-dir conformance-results] [--verbose]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

CATEGORIES="ast,simulation,convergence,regridding,reprojection"
BINDINGS="julia,python,rust,typescript,go"
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

# Python (earthsci_toolkit): prefer the binding's own venv interpreter (the
# layout its pytest suite runs in); fall back to whatever python3 is on PATH
# (e.g. a CI image provisioning numpy/scipy globally).
_python_bin() {
  local venv="$ESS_ROOT/packages/earthsci_toolkit/.venv/bin/python3"
  [[ -x "$venv" ]] && { echo "$venv"; return; }
  command -v python3
}
probe_python() {
  [[ -d "$ESS_ROOT/packages/earthsci_toolkit/src/earthsci_toolkit" ]] &&
    [[ -n "$(_python_bin)" ]]
}
run_python() { # $1 = output dir
  "$(_python_bin)" "$REPO/scripts/runners/run-python.py" \
    --output-dir "$1" --categories "$CATEGORIES" $VERBOSE
}

probe_rust() {
  command -v cargo >/dev/null 2>&1 &&
    [[ -f "$ESS_ROOT/packages/earthsci-toolkit-rs/Cargo.toml" ]]
}
run_rust() { # $1 = output dir
  "$REPO/scripts/runners/run-rust.sh" \
    --output-dir "$1" --categories "$CATEGORIES" $VERBOSE
}

probe_typescript() {
  command -v node >/dev/null 2>&1 &&
    [[ -f "$ESS_ROOT/packages/earthsci-toolkit/dist/cjs/index.js" ]]
}
run_typescript() { # $1 = output dir
  node "$REPO/scripts/runners/run-typescript.js" \
    --output-dir "$1" --categories "$CATEGORIES" $VERBOSE
}

probe_go() {
  command -v go >/dev/null 2>&1 &&
    [[ -d "$ESS_ROOT/packages/esm-format-go/pkg/esm" ]]
}
run_go() { # $1 = output dir
  "$REPO/scripts/runners/run-go.sh" \
    --output-dir "$1" --categories "$CATEGORIES" $VERBOSE
}

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
