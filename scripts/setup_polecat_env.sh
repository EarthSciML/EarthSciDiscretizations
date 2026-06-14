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
#   scripts/setup_polecat_env.sh                # resolves into ./Project.toml's env
#   JULIA_PROJECT=docs scripts/setup_polecat_env.sh   # resolves into docs/Project.toml's env
#   EARTHSCI_SERIALIZATION_REV=abc123 scripts/setup_polecat_env.sh   # pin URL fallback
#
# After this runs, `julia --project=$JULIA_PROJECT -e 'using EarthSciSerialization'` works.
# Re-running is idempotent: Pkg.develop on the same path is a no-op.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

julia_project="${JULIA_PROJECT:-.}"

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

# Snapshot the [sources] block before any Pkg operation. Pkg.develop and
# Pkg.add(url=...) rewrite Project.toml and silently drop [sources], breaking
# fresh-clone Pkg.instantiate() for unregistered packages like ESS.
_project_file="${julia_project%/}/Project.toml"
_sources_block=""
if grep -q '^\[sources\]' "$_project_file" 2>/dev/null; then
    _sources_block=$(awk '
        /^\[sources\]/ { in_sources = 1 }
        in_sources && /^\[/ && !/^\[sources\]/ { in_sources = 0 }
        in_sources { print }
    ' "$_project_file")
fi

if [[ -n "$ess_path" ]]; then
  echo "Pkg.develop EarthSciSerialization from: $ess_path (project=$julia_project)"
  julia --project="$julia_project" -e "using Pkg; Pkg.develop(path=\"$ess_path\"); Pkg.instantiate()"
else
  rev="${EARTHSCI_SERIALIZATION_REV:-main}"
  echo "No local EarthSciSerialization.jl checkout found; Pkg.add from GitHub (rev=$rev, project=$julia_project)"
  julia --project="$julia_project" -e "using Pkg; Pkg.add(url=\"https://github.com/EarthSciML/EarthSciSerialization.git\", rev=\"$rev\", subdir=\"packages/EarthSciSerialization.jl\"); Pkg.instantiate()"
fi

# Restore [sources] if Pkg dropped it. The block was snapshotted above and
# is appended back to keep Project.toml committable without the entry missing.
if [[ -n "$_sources_block" ]] && ! grep -q '^\[sources\]' "$_project_file" 2>/dev/null; then
    printf '\n%s\n' "$_sources_block" >> "$_project_file"
    echo "Re-asserted [sources] block in $_project_file (Pkg dropped it during env setup)." >&2
fi
