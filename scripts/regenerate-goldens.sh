#!/usr/bin/env bash
# regenerate-goldens.sh [category …] — Julia-reference-only golden regeneration.
#
# Drives the canonical ESS pipeline through the official Julia runner
# (scripts/runners/run-julia.jl) and installs its outputs as the committed
# goldens (AGENTS.md §5: regeneration is Julia-reference-only; a PR that
# changes goldens must say why):
#
#   ast          tests/conformance/ast/<case>/expanded.golden.json
#                (the runner's canonical post-lowering bytes, installed
#                verbatim — byte-compared by every other binding). A golden
#                >= 64 KiB (AST_GOLDEN_MAX_BYTES) is instead installed as
#                expanded.golden.sha256 (a digest of the bytes; the full .json
#                is removed) so multi-MB blobs stay out of git — compare-outputs.py
#                then gates every binding's output on the hash. Regenerate to inspect.
#   convergence  tests/conformance/convergence/<case>/golden/errors.json
#                (format per scripts/check_convergence_order.py docstring)
#   regridding   tests/conformance/regridding/<case>/golden/weights.json
#                (the per-pair A_ij/A_j/W_ij setup arrays the runner reads
#                from the official EarthSciSerialization.jl BuildInspection
#                surface; gated by compare-outputs.py per manifest §5.8)
#
# Manifests whose golden this run fulfills have "status": "pending-golden"
# (or the regridding interim "active-invariants") updated to "active".
# Idempotent: a second run is byte-identical.
#
# Usage: scripts/regenerate-goldens.sh [ast] [convergence] [regridding]
#                                      [--files <manifest.json>[,<manifest.json>…]]
#        (no arguments = all three categories; --files restricts regeneration to
#        the named case manifests — passed through to the Julia runner — leaving
#        every other case's golden untouched)

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

CATEGORIES=()
FILES_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --files) FILES_CSV="$2"; shift 2 ;;
    *) CATEGORIES+=("$1"); shift ;;
  esac
done
[[ ${#CATEGORIES[@]} -eq 0 ]] && CATEGORIES=(ast convergence regridding)
for c in "${CATEGORIES[@]}"; do
  case "$c" in
    ast|convergence|regridding) ;;
    *) echo "error: regenerate-goldens.sh handles categories 'ast', 'convergence', and 'regridding', got '$c'" >&2
       exit 2 ;;
  esac
done

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/esd-goldens.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT

CATS_CSV="$(IFS=,; echo "${CATEGORIES[*]}")"
echo "== running the Julia reference runner (categories: $CATS_CSV)"
if [[ -n "$FILES_CSV" ]]; then
  julia "$REPO/scripts/runners/run-julia.jl" --output-dir "$STAGE" \
    --categories "$CATS_CSV" --files "$FILES_CSV"
else
  julia "$REPO/scripts/runners/run-julia.jl" --output-dir "$STAGE" --categories "$CATS_CSV"
fi

# Install artifacts + refresh manifests. Stdlib-python JSON handling only
# (structural file I/O, no numerics — the numbers all came from the runner).
python3 - "$REPO" "$STAGE" "${CATEGORIES[@]}" <<'PY'
import json
import pathlib
import sys

repo = pathlib.Path(sys.argv[1])
stage = pathlib.Path(sys.argv[2])
categories = sys.argv[3:]


def activate(manifest_path: pathlib.Path) -> None:
    """Flip "status": "pending-golden" (or the regridding interim
    "active-invariants") to "active" (text-level, format-preserving)."""
    text = manifest_path.read_text()
    updated = text.replace('"status": "pending-golden"', '"status": "active"')
    updated = updated.replace('"status": "active-invariants"', '"status": "active"')
    if updated != text:
        manifest_path.write_text(updated)
        print(f"   manifest: {manifest_path.relative_to(repo)} -> status: active")


AST_GOLDEN_MAX_BYTES = 65536  # AST goldens >= this are checked in as a sha256 DIGEST
# (expanded.golden.sha256) rather than full bytes — the canonical bytes are pinned by hash,
# not stored, so the repo never carries multi-MB derived blobs (e.g. WENO/HJ-WENO). Every
# binding still reproduces the exact bytes; compare-outputs.py hashes each and checks the digest.

if "ast" in categories and (stage / "julia_ast_results.json").is_file():
    import hashlib
    results = json.loads((stage / "julia_ast_results.json").read_text())
    for case, rec in sorted(results["cases"].items()):
        if rec["status"] != "ok":
            sys.exit(f"ast/{case}: runner status {rec['status']}: {rec.get('message', '')}")
        case_dir = repo / "tests" / "conformance" / "ast" / case
        data = (stage / rec["artifact"]).read_bytes()
        full = case_dir / json.loads((case_dir / "manifest.json").read_text())["golden"]
        digest = case_dir / "expanded.golden.sha256"
        if len(data) >= AST_GOLDEN_MAX_BYTES:
            digest.write_text(json.dumps({
                "algorithm": "sha256",
                "hex": hashlib.sha256(data).hexdigest(),
                "bytes": len(data),
                "generated_by": "scripts/regenerate-goldens.sh",
                "note": (f"Digest of the canonical post-lowering AST ({len(data)} B). The full "
                         f"golden exceeds the {AST_GOLDEN_MAX_BYTES} B checked-in threshold, so it "
                         "is pinned by hash rather than stored. Every binding must reproduce these "
                         "exact bytes; regenerate locally with scripts/regenerate-goldens.sh to inspect."),
            }, indent=2, ensure_ascii=False) + "\n")
            if full.suffix == ".json" and full.is_file():
                full.unlink()
            print(f"   wrote {digest.relative_to(repo)} (sha256 of {len(data)} B; full golden not checked in)")
        else:
            full.write_bytes(data)
            if digest.is_file():
                digest.unlink()
            print(f"   wrote {full.relative_to(repo)} ({len(data)} bytes)")
        activate(case_dir / "manifest.json")

if "convergence" in categories and (stage / "julia_convergence_results.json").is_file():
    results = json.loads((stage / "julia_convergence_results.json").read_text())
    for case, rec in sorted(results["cases"].items()):
        if rec["status"] != "ok":
            sys.exit(f"convergence/{case}: runner status {rec['status']}: {rec.get('message', '')}")
        case_dir = repo / "tests" / "conformance" / "convergence" / case
        manifest = json.loads((case_dir / "manifest.json").read_text())
        golden = case_dir / manifest.get("golden", "golden/errors.json")
        golden.parent.mkdir(parents=True, exist_ok=True)
        doc = {
            "case": case,
            "binding": "julia",
            "generated_by": "scripts/regenerate-goldens.sh",
            "assert_time": rec["assert_time"],
            "errors": sorted(rec["errors"], key=lambda r: r["n"]),
        }
        golden.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
        print(f"   wrote {golden.relative_to(repo)}")
        activate(case_dir / "manifest.json")

if "regridding" in categories and (stage / "julia_regridding_results.json").is_file():
    results = json.loads((stage / "julia_regridding_results.json").read_text())
    for case, rec in sorted(results["cases"].items()):
        if rec["status"] != "ok":
            sys.exit(f"regridding/{case}: runner status {rec['status']}: {rec.get('message', '')}")
        missing = [k for k in ("A_ij", "A_j", "W_ij") if k not in rec]
        if missing:
            sys.exit(f"regridding/{case}: runner emitted no {missing} "
                     "(BuildInspection setup arrays)")
        case_dir = repo / "tests" / "conformance" / "regridding" / case
        manifest = json.loads((case_dir / "manifest.json").read_text())
        golden_rel = manifest.get("golden", "golden/weights.json")
        if not golden_rel.endswith(".json"):
            golden_rel = "golden/weights.json"
        golden = case_dir / golden_rel
        golden.parent.mkdir(parents=True, exist_ok=True)
        doc = {
            "case": case,
            "binding": "julia",
            "generated_by": "scripts/regenerate-goldens.sh",
            "A_ij": rec["A_ij"],
            "A_j": rec["A_j"],
            "W_ij": rec["W_ij"],
        }
        golden.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
        print(f"   wrote {golden.relative_to(repo)}")
        activate(case_dir / "manifest.json")
PY

echo "== done; verify with scripts/check_convergence_order.py and scripts/test-conformance.sh"
