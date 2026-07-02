#!/usr/bin/env bash
# run-rust.sh — the official Rust conformance runner for EarthSciDiscretizations.
#
# Thin shell over the earthsci-toolkit-rs `canonical_expand` example (the
# official raw §9.7 pipeline resolve_template_machinery →
# lower_expression_templates with the reference canonical-byte writer,
# upstream in $ESS_ROOT/packages/earthsci-toolkit-rs/examples/). This wrapper
# implements the §4.2 CLI: manifest discovery, per-case invocation, artifact
# and result-JSON assembly (stdlib-python structural I/O only — every byte in
# every artifact comes from the cargo example; AGENTS.md §2).
#
# Usage: scripts/runners/run-rust.sh --output-dir <path> [--categories ast]
#                                    [--files <manifest.json>[,…]] [--verbose]
#
# Scope: ast. The Rust binding's simulator (simulate_array) does not yet
# execute §6.6.5 reduction assertions or coordinate-expression ic seeding —
# the simulation/convergence manifests carry the blocked-upstream notes.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

MANIFEST_PATH="$ESS_ROOT/packages/earthsci-toolkit-rs/Cargo.toml"
[[ -f "$MANIFEST_PATH" ]] || {
  echo "error: earthsci-toolkit-rs not found at '$MANIFEST_PATH'" >&2
  exit 1
}

# Build once up front so per-case `cargo run --quiet` is silent and fast.
cargo build --quiet --manifest-path "$MANIFEST_PATH" --example canonical_expand

ESD_ROOT="$REPO" ESS_ROOT="$ESS_ROOT" RUST_MANIFEST_PATH="$MANIFEST_PATH" \
  python3 - "$@" <<'PY'
import json
import os
import subprocess
import sys
from pathlib import Path

ESD_ROOT = Path(os.environ["ESD_ROOT"])
MANIFEST_PATH = os.environ["RUST_MANIFEST_PATH"]
CONF = ESD_ROOT / "tests" / "conformance"

args = sys.argv[1:]
output_dir = None
categories: list[str] = []
files: list[str] = []
verbose = False
i = 0
while i < len(args):
    a = args[i]
    if a == "--output-dir":
        output_dir = args[i + 1]; i += 2
    elif a == "--categories":
        categories += [c for c in args[i + 1].split(",") if c]; i += 2
    elif a == "--files":
        files += [f for f in args[i + 1].split(",") if f]; i += 2
    elif a == "--verbose":
        verbose = True; i += 1
    else:
        sys.exit(f"error: unknown argument '{a}'")
if output_dir is None:
    sys.exit("error: --output-dir is required (CONFORMANCE_SPEC.md §4.2)")
if not categories:
    categories = ["ast"]
if files:
    file_cats = {json.loads((Path(f) if os.path.isabs(f) else ESD_ROOT / f)
                            .read_text())["category"] for f in files}
    categories = [c for c in categories if c in file_cats]

output = Path(output_dir)
output.mkdir(parents=True, exist_ok=True)


def discover(category):
    if files:
        out = []
        for f in files:
            p = Path(f) if os.path.isabs(f) else ESD_ROOT / f
            m = json.loads(p.read_text())
            if m["category"] == category:
                out.append((p.resolve().parent, m))
        return out
    root = CONF / category
    if not root.is_dir():
        return []
    return [(d, json.loads((d / "manifest.json").read_text()))
            for d in sorted(root.iterdir()) if (d / "manifest.json").is_file()]


def unlowered_violations(node, acc):
    """The prose-free structural gates (mirrors run-julia.jl)."""
    if isinstance(node, dict):
        op = str(node.get("op", ""))
        if op == "apply_expression_template":
            acc.append("apply_expression_template")
        if op in ("grad", "div", "laplacian"):
            acc.append(f"unlowered {op}")
        if op == "D" and str(node.get("wrt", "t")) != "t":
            acc.append(f"unlowered spatial D (wrt={node['wrt']})")
        for k, v in node.items():
            if k != "metadata":  # prose may cite op names
                unlowered_violations(v, acc)
    elif isinstance(node, list):
        for v in node:
            unlowered_violations(v, acc)
    return acc


def run_ast():
    cases = {}
    for case_dir, manifest in discover("ast"):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        try:
            fixture = case_dir / manifest["fixture"]
            proc = subprocess.run(
                ["cargo", "run", "--quiet", "--manifest-path", MANIFEST_PATH,
                 "--example", "canonical_expand", "--", str(fixture)],
                capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr.decode("utf-8", "replace").strip())
            data = proc.stdout
            doc = json.loads(data)
            violations = []
            for model in doc.get("models", {}).values():
                unlowered_violations(
                    {k: v for k, v in model.items() if k != "metadata"}, violations)
            if violations:
                raise RuntimeError("post-lowering gate failed: "
                                   + ", ".join(sorted(set(violations))))
            art_dir = output / "ast" / case
            art_dir.mkdir(parents=True, exist_ok=True)
            art = art_dir / "expanded.golden.json"
            art.write_bytes(data)
            rec["artifact"] = str(art.relative_to(output))
            rec["bytes"] = len(data)
            rec["golden"] = str((case_dir / manifest["golden"]).resolve()
                                .relative_to(ESD_ROOT))
        except Exception as err:  # noqa: BLE001 — recorded per-case
            rec["status"] = "error"
            rec["message"] = str(err)
        if verbose:
            print(f"  ast/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "rust", "category": "ast", "cases": cases}


summary = {"binding": "rust", "runner": "scripts/runners/run-rust.sh",
           "categories": {}}
exit_code = 0
for cat in categories:
    if cat != "ast":
        print(f"rust/{cat}: unsupported category (skipping; see the manifests' "
              "blocked-upstream notes for the Rust simulation pathway)")
        continue
    if verbose:
        print(f"category: {cat}")
    res = run_ast()
    (output / f"rust_{cat}_results.json").write_text(
        json.dumps(res, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
    n_bad = sum(1 for r in res["cases"].values() if r["status"] != "ok")
    if n_bad:
        exit_code = 1
    summary["categories"][cat] = {"cases": len(res["cases"]), "failures": n_bad}
    print(f"rust/{cat}: {len(res['cases'])} case(s), {n_bad} failure(s)")
(output / "rust_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
sys.exit(exit_code)
PY
