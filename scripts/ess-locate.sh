#!/usr/bin/env bash
# ess-locate.sh — resolve the EarthSciSerialization checkout this repo tests against.
#
# Every runner and validation script sources this file (or invokes it) instead of
# hard-coding a path. Resolution order:
#   1. $ESS_ROOT, if set.
#   2. ../EarthSciSerialization relative to this repo's root (the sibling-checkout
#      convention used locally and in CI).
#
# Usage:
#   ESS_ROOT="$(scripts/ess-locate.sh)"          # print the resolved path
#   source scripts/ess-locate.sh                  # export ESS_ROOT into the caller

set -euo pipefail

_esd_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"

if [[ -z "${ESS_ROOT:-}" ]]; then
  ESS_ROOT="$(cd "${_esd_root}/.." && pwd)/EarthSciSerialization"
fi

if [[ ! -f "${ESS_ROOT}/esm-schema.json" ]]; then
  echo "error: EarthSciSerialization not found at '${ESS_ROOT}'" >&2
  echo "       set ESS_ROOT or clone EarthSciML/EarthSciSerialization as a sibling checkout" >&2
  exit 1
fi

export ESS_ROOT

# When executed (not sourced), print the resolved path.
if [[ "${BASH_SOURCE[0]:-}" == "$0" ]]; then
  echo "${ESS_ROOT}"
fi
