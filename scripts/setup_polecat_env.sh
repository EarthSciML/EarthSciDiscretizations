#!/usr/bin/env bash
# Resolve EarthSciSerialization.jl into the active Julia environment so this
# project can `using EarthSciSerialization` and call the AST evaluator.
#
# Resolution order:
#   1. $EARTHSCI_SERIALIZATION_PATH (if set and a directory exists there)
#   2. The first matching Gas Town workspace checkout, in priority order:
#        ../../../../EarthSciSerialization/refinery/rig/packages/EarthSciSerialization.jl
#        ../../../../EarthSciSerialization/mayor/rig/packages/EarthSciSerialization.jl
#   3. Fallback: install from https://github.com/EarthSciML/EarthSciSerialization.git@main
#
# Usage:
#   scripts/setup_polecat_env.sh                # resolves into Project.toml's env
#   EARTHSCI_SERIALIZATION_REV=abc123 scripts/setup_polecat_env.sh   # pin URL fallback
#
# After this runs, `julia --project=. -e 'using EarthSciSerialization'` works.
# Re-running is idempotent: Pkg.develop on the same path is a no-op.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

ess_path=""

if [[ -n "${EARTHSCI_SERIALIZATION_PATH:-}" ]]; then
  if [[ -d "$EARTHSCI_SERIALIZATION_PATH" ]]; then
    ess_path="$EARTHSCI_SERIALIZATION_PATH"
  else
    echo "warn: EARTHSCI_SERIALIZATION_PATH=$EARTHSCI_SERIALIZATION_PATH is not a directory; trying workspace defaults" >&2
  fi
fi

if [[ -z "$ess_path" ]]; then
  for candidate in \
    "$repo_root/../../../../EarthSciSerialization/refinery/rig/packages/EarthSciSerialization.jl" \
    "$repo_root/../../../../EarthSciSerialization/mayor/rig/packages/EarthSciSerialization.jl"
  do
    if [[ -f "$candidate/Project.toml" ]]; then
      ess_path="$(cd "$candidate" && pwd)"
      break
    fi
  done
fi

if [[ -n "$ess_path" ]]; then
  echo "Pkg.develop EarthSciSerialization from: $ess_path"
  julia --project=. -e "using Pkg; Pkg.develop(path=\"$ess_path\"); Pkg.instantiate()"
else
  rev="${EARTHSCI_SERIALIZATION_REV:-main}"
  echo "No local EarthSciSerialization.jl checkout found; Pkg.add from GitHub (rev=$rev)"
  julia --project=. -e "using Pkg; Pkg.add(url=\"https://github.com/EarthSciML/EarthSciSerialization.git\", rev=\"$rev\", subdir=\"packages/EarthSciSerialization.jl\"); Pkg.instantiate()"
fi
