#!/usr/bin/env bash
# ess-locate.sh — resolve the EarthSciAST checkout this repo tests against.
#
# Every runner and validation script sources this file (or invokes it) instead of
# hard-coding a path. Resolution order:
#   1. $ESS_ROOT, if set.
#   2. The nearest ../EarthSciAST found by walking up from this repo's root. The
#      first candidate is the ordinary sibling checkout used locally and in CI.
#      The walk matters inside a git worktree, where the repo root is
#      <main>/.claude/worktrees/<name> and the sibling checkout is several levels
#      up rather than one -- resolving only one level up yields a path that does
#      not exist, and the resulting failure is easy to misread as a test failure.
#
# Usage:
#   ESS_ROOT="$(scripts/ess-locate.sh)"          # print the resolved path
#   source scripts/ess-locate.sh                  # export ESS_ROOT into the caller

set -euo pipefail

_esd_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"

if [[ -z "${ESS_ROOT:-}" ]]; then
  _dir="${_esd_root}"
  while true; do
    _cand="$(cd "${_dir}/.." && pwd)/EarthSciAST"
    if [[ -f "${_cand}/esm-schema.json" ]]; then
      ESS_ROOT="${_cand}"
      break
    fi
    _parent="$(dirname "${_dir}")"
    [[ "${_parent}" == "${_dir}" ]] && break
    _dir="${_parent}"
  done
  # Fall back to the plain sibling guess so the error below names the expected path.
  ESS_ROOT="${ESS_ROOT:-$(cd "${_esd_root}/.." && pwd)/EarthSciAST}"
fi

if [[ ! -f "${ESS_ROOT}/esm-schema.json" ]]; then
  echo "error: EarthSciAST not found at '${ESS_ROOT}'" >&2
  echo "       set ESS_ROOT or clone EarthSciML/EarthSciAST as a sibling checkout" >&2
  exit 1
fi

export ESS_ROOT

# When executed (not sourced), print the resolved path.
if [[ "${BASH_SOURCE[0]:-}" == "$0" ]]; then
  echo "${ESS_ROOT}"
fi
