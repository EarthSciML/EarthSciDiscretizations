#!/usr/bin/env bash
# run-rust.sh — the official Rust conformance runner for EarthSciDiscretizations.
#
# Thin shell over two earthsci-toolkit-rs examples (upstream in
# $ESS_ROOT/packages/earthsci-toolkit-rs/examples/):
#
#   canonical_expand   the official raw §9.7 pipeline
#                      (resolve_template_machinery → lower_expression_templates)
#                      with the reference canonical-byte writer        [ast]
#   pde_conformance    the official simulation / evaluation pathways:
#                      run_pde_tests (§6.6/§6.6.5 inline tests over
#                      simulate, solver pinned by the manifest)        [simulation]
#                      per-resolution loader-API metaparameter load →
#                      simulate → evaluate_cellwise → field_reduce     [convergence]
#                      runner-built wrapper doc importing the library,
#                      §9.7.3 body composition, evaluate per point     [reprojection]
#                      run_pde_tests exact-invariant gates + the
#                      regrid/pou/cons fields at t=1 + the per-pair
#                      A_ij/A_j/W_ij setup arrays via the official
#                      BuildInspection surface
#                      (simulate_with_inspection)                      [regridding]
#
# This wrapper implements the §4.2 CLI: manifest discovery, per-case
# invocation, artifact and result-JSON assembly (stdlib-python structural I/O
# only — every number and byte in every artifact comes from the cargo
# examples; AGENTS.md §2).
#
# Usage: scripts/runners/run-rust.sh --output-dir <path>
#            [--categories ast,simulation,convergence,regridding,reprojection]
#            [--files <manifest.json>[,…]] [--verbose]

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$REPO/scripts/ess-locate.sh"

MANIFEST_PATH="$ESS_ROOT/packages/earthsci-toolkit-rs/Cargo.toml"
[[ -f "$MANIFEST_PATH" ]] || {
  echo "error: earthsci-toolkit-rs not found at '$MANIFEST_PATH'" >&2
  exit 1
}

# Build once up front so per-case `cargo run --quiet` is silent and fast.
# Dev profile, matching the upstream crate's own test/CI builds (the release
# profile would rebuild the vendored s2 geometry kernel from scratch, which
# needs a full OpenSSL dev environment; the numbers are profile-independent —
# Rust f64 arithmetic is IEEE-strict at every opt level).
cargo build --quiet --manifest-path "$MANIFEST_PATH" \
  --example canonical_expand --example pde_conformance

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
CATEGORIES_SUPPORTED = ["ast", "simulation", "convergence", "regridding",
                        "reprojection"]

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
    categories = list(CATEGORIES_SUPPORTED)
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


def rust_out_of_scope(manifest):
    """Manifest-declared reasons this binding must not run the case (the
    same markers compare-outputs.py keys its SKIP rows on)."""
    for key in ("scope_excluded", "blocked_upstream_bindings"):
        if "rust" in manifest.get(key, {}):
            return manifest[key]["rust"]
    return None


def cargo_example(example, *argv):
    """Run one upstream cargo example; return parsed-JSON stdout."""
    proc = subprocess.run(
        ["cargo", "run", "--quiet", "--manifest-path",
         MANIFEST_PATH, "--example", example, "--", *argv],
        capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode("utf-8", "replace").strip())
    return proc.stdout


def integrator_opts(manifest):
    j = manifest.get("integrators", {}).get("rust")
    if j is None:
        return "Erk", 1e-10, 1e-12
    return (str(j.get("solver", "Erk")), float(j.get("reltol", 1e-10)),
            float(j.get("abstol", 1e-12)))


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
            data = cargo_example("canonical_expand", str(fixture))
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


def run_simulation():
    cases = {}
    for case_dir, manifest in discover("simulation"):
        case = manifest["case"]
        reason = rust_out_of_scope(manifest)
        if reason is not None:
            if verbose:
                print(f"  simulation/{case}: skipped ({reason})")
            continue
        rec = {"case": case, "status": "ok"}
        try:
            problem = (case_dir / manifest["problem"]).resolve()
            solver, reltol, abstol = integrator_opts(manifest)
            out = json.loads(cargo_example(
                "pde_conformance", "pde-tests", str(problem),
                "--model", manifest["model"], "--solver", solver,
                "--reltol", repr(reltol), "--abstol", repr(abstol)))
            wanted = set(manifest["tests"])
            assertions = [a for a in out["assertions"] if a["test_id"] in wanted]
            if not assertions:
                raise RuntimeError(f"no assertions ran for tests {sorted(wanted)}")
            rec["assertions"] = assertions
            rec["passed"] = all(a["passed"] for a in assertions)
            rec["status"] = "ok" if rec["passed"] else "failed"
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = str(err)
        if verbose:
            print(f"  simulation/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "rust", "category": "simulation", "cases": cases}


def run_convergence():
    cases = {}
    for case_dir, manifest in discover("convergence"):
        case = manifest["case"]
        reason = rust_out_of_scope(manifest)
        if reason is not None:
            if verbose:
                print(f"  convergence/{case}: skipped ({reason})")
            continue
        rec = {"case": case, "status": "ok"}
        try:
            problem = (case_dir / manifest["problem"]).resolve()
            solver, reltol, abstol = integrator_opts(manifest)
            out = json.loads(cargo_example(
                "pde_conformance", "convergence", str(problem),
                "--model", manifest["model"],
                "--assert-time", repr(float(manifest["assert_time"])),
                "--solver", solver,
                "--reltol", repr(reltol), "--abstol", repr(abstol),
                "--norms", ",".join(manifest["norms"]),
                "--resolutions", json.dumps(manifest["resolutions"])))
            rec["assert_time"] = out["assert_time"]
            rec["errors"] = out["errors"]
            if verbose:
                for row in out["errors"]:
                    print(f"  convergence/{case} n={row['n']}: " + " ".join(
                        f"{k}={row[k]}" for k in sorted(row) if k != "n"))
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = str(err)
        if verbose:
            print(f"  convergence/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "rust", "category": "convergence", "cases": cases}


def run_regridding():
    # Mirrors run-julia.jl's run_regridding record field-for-field: the
    # fixture's exact-invariant inline tests (run_pde_tests), the recorded
    # regrid/pou/cons state fields at t=1, and the per-pair A_ij/A_j/W_ij
    # setup arrays from the upstream BuildInspection surface — all from ONE
    # `pde_conformance regrid` invocation (solver pinning identical to the
    # Julia runner's hardcoded Tsit5/1e-10/1e-12 regridding setup).
    cases = {}
    for case_dir, manifest in discover("regridding"):
        case = manifest["case"]
        reason = rust_out_of_scope(manifest)
        if reason is not None:
            if verbose:
                print(f"  regridding/{case}: skipped ({reason})")
            continue
        rec = {"case": case, "status": "ok"}
        try:
            fixture = (case_dir / manifest["fixture"]).resolve()
            out = json.loads(cargo_example(
                "pde_conformance", "regrid", str(fixture),
                "--model", manifest["model"], "--solver", "Erk",
                "--reltol", "1e-10", "--abstol", "1e-12"))
            for key in ("assertions", "passed", "regrid_state_at_1",
                        "pou_state_at_1", "cons_state_at_1",
                        "A_ij", "A_j", "W_ij"):
                rec[key] = out[key]
            rec["status"] = "ok" if rec["passed"] else "failed"
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = str(err)
        if verbose:
            print(f"  regridding/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "rust", "category": "regridding", "cases": cases}


def run_reprojection():
    cases = {}
    for case_dir, manifest in discover("reprojection"):
        case = manifest["case"]
        reason = rust_out_of_scope(manifest)
        if reason is not None:
            if verbose:
                print(f"  reprojection/{case}: skipped ({reason})")
            continue
        rec = {"case": case, "status": "ok"}
        try:
            library = (case_dir / manifest["library"]).resolve()
            points = (case_dir / manifest["golden"]).resolve()
            out = json.loads(cargo_example(
                "pde_conformance", "reproject", str(library), str(points)))
            rec["forward"] = out["forward"]
            rec["inverse"] = out["inverse"]
            rec["roundtrip"] = out["roundtrip"]
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = str(err)
        if verbose:
            print(f"  reprojection/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "rust", "category": "reprojection", "cases": cases}


RUNNERS = {"ast": run_ast, "simulation": run_simulation,
           "convergence": run_convergence, "regridding": run_regridding,
           "reprojection": run_reprojection}

summary = {"binding": "rust", "runner": "scripts/runners/run-rust.sh",
           "categories": {}}
exit_code = 0
for cat in categories:
    if cat not in RUNNERS:
        print(f"rust/{cat}: unsupported category (skipping; this runner "
              f"implements {CATEGORIES_SUPPORTED}; see the regridding "
              "manifest's blocked-upstream notes)")
        continue
    if verbose:
        print(f"category: {cat}")
    res = RUNNERS[cat]()
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
