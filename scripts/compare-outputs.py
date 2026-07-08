#!/usr/bin/env python3
"""compare-outputs.py — gate runner outputs against goldens and against each other.

Consumes the per-binding result trees the conformance runners write
(conformance-results/<binding>/<binding>_<category>_results.json, per
CONFORMANCE_SPEC.md §4.2/§4.4) and applies each category's comparison contract:

  ast          byte-identity: the runner's canonical post-lowering artifact vs
               the committed golden, and vs the reference binding's artifact
               (CONFORMANCE_SPEC.md §5.5.3 canonical bytes).
  simulation   every inline assertion passed (the tolerances live in the
               problem files, §6.6.4/§6.6.5); actual values cross-checked
               against the reference binding at rtol 1e-9 / atol 1e-12 (pure
               re-evaluation of the same pinned integrator setup).
  convergence  error norms vs the committed golden/errors.json at the
               manifest's error_vs_golden_rtol/atol.
  regridding   every exact-invariant assertion passed; recorded fields carried
               through for cross-binding comparison at the manifest's
               weight_rtol/weight_atol; per-pair A_ij/A_j/W_ij setup arrays vs
               the committed golden/weights.json (areas at area_rtol + the
               5.8.2 sliver-floor atol, weights at weight_rtol/weight_atol) for
               every runner that emits them (SKIP, never silent pass, for a
               binding without an upstream inspection surface).
  reprojection forward/inverse/roundtrip values vs the independent analytic
               reference (golden/points.json) at the manifest tolerances
               (combined |v - ref| <= atol + rtol*|ref| checks).

Pure stdlib; performs no model evaluation — only comparisons of runner-emitted
numbers/bytes (single-pathway doctrine, AGENTS.md §2).

Usage:
  python3 scripts/compare-outputs.py --results-root conformance-results \
      [--bindings julia[,…]] [--categories ast,…] [--verbose]

Writes <results-root>/comparison/summary.json and prints a human summary.
Exit 0 iff every gate passed (pending-golden cases are reported as SKIP).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONF = REPO / "tests" / "conformance"
ALL_CATEGORIES = ["ast", "simulation", "convergence", "regridding", "reprojection"]

# Cross-binding tolerance for simulation actuals vs the reference binding:
# both bindings integrate the same pinned setup, so this is a tight
# re-evaluation band (cf. ESS pde_simulation traj_golden_* rationale).
SIM_CROSS_RTOL = 1e-9
SIM_CROSS_ATOL = 1e-12


def close(actual: float, ref: float, rtol: float, atol: float) -> bool:
    return abs(actual - ref) <= atol + rtol * abs(ref)


class Report:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def add(self, category: str, case: str, binding: str, gate: str,
            status: str, detail: str = "") -> None:
        self.rows.append({"category": category, "case": case, "binding": binding,
                          "gate": gate, "status": status, "detail": detail})
        mark = {"pass": "ok  ", "skip": "SKIP", "fail": "FAIL"}[status]
        line = f"{mark} {category}/{case} [{binding}] {gate}"
        if detail:
            line += f": {detail}"
        print(line)

    @property
    def failed(self) -> bool:
        return any(r["status"] == "fail" for r in self.rows)


def manifests(category: str):
    root = CONF / category
    if not root.is_dir():
        return
    for case_dir in sorted(root.iterdir()):
        mp = case_dir / "manifest.json"
        if mp.is_file():
            yield case_dir, json.loads(mp.read_text())


def binding_results(results_root: Path, binding: str, category: str):
    p = results_root / binding / f"{binding}_{category}_results.json"
    if not p.is_file():
        return None
    return json.loads(p.read_text())


def case_rec(results, case: str):
    if results is None:
        return None
    return results.get("cases", {}).get(case)


def compare_ast(results_root, bindings, rep):
    import hashlib
    for case_dir, manifest in manifests("ast"):
        case = manifest["case"]
        # The committed golden is either the full canonical bytes (small cases) or, for
        # cases whose lowered AST exceeds the regenerate threshold, a sha256 DIGEST
        # (`expanded.golden.sha256`) that pins the bytes without storing them. Either way
        # every binding must reproduce the exact same canonical bytes; the digest just
        # records their hash (regression pin + cross-binding gate, minus the multi-MB blob).
        digest_path = case_dir / "expanded.golden.sha256"
        full_path = case_dir / manifest["golden"]
        golden = None       # full bytes, if stored
        golden_hex = None   # committed sha256 hex, if digest-pinned
        golden_n = None
        if digest_path.is_file():
            dinfo = json.loads(digest_path.read_text())
            golden_hex, golden_n = dinfo["hex"], dinfo.get("bytes")
        elif full_path.is_file():
            golden = full_path.read_bytes()
        reference = manifest.get("reference_binding", "julia")
        ref_artifact = None
        for binding in bindings:
            rec = case_rec(binding_results(results_root, binding, "ast"), case)
            if rec is None:
                rep.add("ast", case, binding, "runner-output", "skip", "no result recorded")
                continue
            if rec["status"] != "ok":
                rep.add("ast", case, binding, "runner-output", "fail",
                        rec.get("message", rec["status"]))
                continue
            artifact = (results_root / binding / rec["artifact"]).read_bytes()
            if binding == reference:
                ref_artifact = artifact
            if golden_hex is not None:
                h = hashlib.sha256(artifact).hexdigest()
                if h == golden_hex:
                    rep.add("ast", case, binding, "golden-digest", "pass",
                            f"sha256 match ({len(artifact)} bytes; golden pinned by digest)")
                else:
                    rep.add("ast", case, binding, "golden-digest", "fail",
                            f"sha256 {h[:12]}… != committed {golden_hex[:12]}… "
                            f"({len(artifact)} B vs {golden_n} B golden)")
            elif golden is not None:
                if artifact == golden:
                    rep.add("ast", case, binding, "golden-bytes", "pass",
                            f"{len(artifact)} bytes identical")
                else:
                    rep.add("ast", case, binding, "golden-bytes", "fail",
                            f"artifact ({len(artifact)} B) != golden ({len(golden)} B)")
            else:
                status = "skip" if manifest.get("status") == "pending-golden" else "fail"
                rep.add("ast", case, binding, "golden-bytes", status,
                        "golden not yet generated" if status == "skip"
                        else "golden missing and not marked pending-golden")
            if binding != reference and ref_artifact is not None:
                rep.add("ast", case, binding, "cross-binding-bytes",
                        "pass" if artifact == ref_artifact else "fail")


def compare_simulation(results_root, bindings, rep):
    for _case_dir, manifest in manifests("simulation"):
        case = manifest["case"]
        reference = manifest.get("reference_binding", "julia")
        ref_actuals = {}
        for binding in bindings:
            if binding in manifest.get("scope_excluded", {}):
                rep.add("simulation", case, binding, "scope", "skip",
                        manifest["scope_excluded"][binding])
                continue
            if binding in manifest.get("blocked_upstream_bindings", {}):
                rep.add("simulation", case, binding, "blocked-upstream", "skip",
                        manifest["blocked_upstream_bindings"][binding])
                continue
            rec = case_rec(binding_results(results_root, binding, "simulation"), case)
            if rec is None:
                rep.add("simulation", case, binding, "runner-output", "skip",
                        "no result recorded")
                continue
            if rec["status"] not in ("ok", "failed"):
                rep.add("simulation", case, binding, "runner-output", "fail",
                        rec.get("message", rec["status"]))
                continue
            bad = [a for a in rec["assertions"] if not a["passed"]]
            if bad:
                first = bad[0]
                rep.add("simulation", case, binding, "inline-assertions", "fail",
                        f"{len(bad)} failing; first: {first['test_id']}"
                        f"[{first['assertion_idx']}] {first.get('message', '')}")
            else:
                detail = "; ".join(
                    f"{a['test_id']}[{a['assertion_idx']}] {a['reduce'] or 'point'}"
                    f"={a['actual']:.6e}" for a in rec["assertions"])
                rep.add("simulation", case, binding, "inline-assertions", "pass", detail)
            key = lambda a: (a["test_id"], a["assertion_idx"])
            if binding == reference:
                ref_actuals = {key(a): a["actual"] for a in rec["assertions"]}
            elif ref_actuals:
                mismatches = [
                    a for a in rec["assertions"]
                    if a["actual"] is not None and key(a) in ref_actuals
                    and not close(a["actual"], ref_actuals[key(a)],
                                  SIM_CROSS_RTOL, SIM_CROSS_ATOL)]
                rep.add("simulation", case, binding, "cross-binding-actuals",
                        "fail" if mismatches else "pass",
                        f"{len(mismatches)} outside rtol={SIM_CROSS_RTOL}" if mismatches else "")


def compare_convergence(results_root, bindings, rep):
    for case_dir, manifest in manifests("convergence"):
        case = manifest["case"]
        golden_path = case_dir / manifest.get("golden", "golden/errors.json")
        golden = json.loads(golden_path.read_text()) if golden_path.is_file() else None
        rtol = manifest["tolerances"]["error_vs_golden_rtol"]
        atol = manifest["tolerances"]["error_vs_golden_atol"]
        for binding in bindings:
            if binding in manifest.get("scope_excluded", {}):
                rep.add("convergence", case, binding, "scope", "skip",
                        manifest["scope_excluded"][binding])
                continue
            if binding in manifest.get("blocked_upstream_bindings", {}):
                rep.add("convergence", case, binding, "blocked-upstream", "skip",
                        manifest["blocked_upstream_bindings"][binding])
                continue
            rec = case_rec(binding_results(results_root, binding, "convergence"), case)
            if rec is None:
                rep.add("convergence", case, binding, "runner-output", "skip",
                        "no result recorded")
                continue
            if rec["status"] != "ok":
                rep.add("convergence", case, binding, "runner-output", "fail",
                        rec.get("message", rec["status"]))
                continue
            if golden is None:
                status = "skip" if manifest.get("status") == "pending-golden" else "fail"
                rep.add("convergence", case, binding, "errors-vs-golden", status,
                        "golden not yet generated" if status == "skip"
                        else "golden missing and not marked pending-golden")
                continue
            grows = {r["n"]: r for r in golden["errors"]}
            bad = []
            for row in rec["errors"]:
                g = grows.get(row["n"])
                if g is None:
                    bad.append(f"n={row['n']} not in golden")
                    continue
                for norm in manifest["norms"]:
                    if not close(row[norm], g[norm], rtol, atol):
                        bad.append(f"n={row['n']} {norm}: {row[norm]} vs {g[norm]}")
            if sorted(grows) != sorted(r["n"] for r in rec["errors"]):
                bad.append(f"resolution sets differ: {sorted(r['n'] for r in rec['errors'])} "
                           f"vs golden {sorted(grows)}")
            rep.add("convergence", case, binding, "errors-vs-golden",
                    "fail" if bad else "pass",
                    "; ".join(bad) if bad else
                    f"{len(rec['errors'])} resolutions within rtol={rtol}")


def _flat_floats(x):
    """Flatten a (possibly nested) list of numbers to a flat float list."""
    if isinstance(x, list):
        out = []
        for v in x:
            out.extend(_flat_floats(v))
        return out
    return [float(x)]


def compare_regridding(results_root, bindings, rep):
    for case_dir, manifest in manifests("regridding"):
        case = manifest["case"]
        reference = manifest.get("reference_binding", "julia")
        tol = manifest.get("tolerances", {})
        wr = tol.get("weight_rtol", 1e-9)
        wa = tol.get("weight_atol", 1e-12)
        # Areas (A_ij / A_j) are gated with the §5.8.2 area band: the calibrated
        # rel tolerance plus the sliver floor atol = factor * L^2.
        ar = tol.get("area_rtol", 1e-9)
        aa = (tol.get("area_atol_factor", 1e-15)
              * tol.get("characteristic_length", 1.0) ** 2)
        golden_rel = manifest.get("golden", "")
        golden_path = case_dir / golden_rel if golden_rel else None
        gold = (json.loads(golden_path.read_text())
                if golden_path is not None and golden_path.is_file() else None)
        ref_fields = {}
        for binding in bindings:
            if binding in manifest.get("scope_excluded", {}):
                rep.add("regridding", case, binding, "scope", "skip",
                        manifest["scope_excluded"][binding])
                continue
            if binding in manifest.get("blocked_upstream_bindings", {}):
                rep.add("regridding", case, binding, "blocked-upstream", "skip",
                        manifest["blocked_upstream_bindings"][binding])
                continue
            rec = case_rec(binding_results(results_root, binding, "regridding"), case)
            if rec is None:
                rep.add("regridding", case, binding, "runner-output", "skip",
                        "no result recorded")
                continue
            if rec["status"] not in ("ok", "failed"):
                rep.add("regridding", case, binding, "runner-output", "fail",
                        rec.get("message", rec["status"]))
                continue
            bad = [a for a in rec["assertions"] if not a["passed"]]
            rep.add("regridding", case, binding, "exact-invariants",
                    "fail" if bad else "pass",
                    f"{len(bad)} failing assertions" if bad else
                    f"{len(rec['assertions'])} assertions (partition-of-unity, conservation)")
            fields = {k: v for k, v in rec.items() if k.endswith("_at_1")}
            if binding == reference:
                ref_fields = fields
            elif ref_fields:
                mism = [k for k in fields
                        if k in ref_fields and any(
                            not close(a, b, wr, wa)
                            for a, b in zip(fields[k], ref_fields[k]))]
                rep.add("regridding", case, binding, "cross-binding-fields",
                        "fail" if mism else "pass", ", ".join(mism))
            # Per-pair setup arrays (§5.8.1-5.8.2): gate the runner-emitted
            # A_ij / A_j / W_ij against golden/weights.json. A runner that does
            # not (yet) emit them — its binding lacks an upstream inspection
            # surface — is reported as SKIP, never silently passed.
            per_pair = {k: rec.get(k) for k in ("A_ij", "A_j", "W_ij")}
            if all(v is not None for v in per_pair.values()):
                if gold is None:
                    status = ("skip" if manifest.get("status")
                              in ("pending-golden", "active-invariants") else "fail")
                    rep.add("regridding", case, binding, "per-pair-weights", status,
                            "golden not yet generated" if status == "skip"
                            else f"golden '{golden_rel}' missing and manifest not "
                                 "marked pending-golden")
                else:
                    bad = []
                    for label, rt, at in (("A_ij", ar, aa), ("A_j", ar, aa),
                                          ("W_ij", wr, wa)):
                        got = _flat_floats(per_pair[label])
                        ref = _flat_floats(gold[label])
                        if len(got) != len(ref):
                            bad.append(f"{label}: {len(got)} entries vs golden "
                                       f"{len(ref)}")
                            continue
                        for k, (a, b) in enumerate(zip(got, ref)):
                            if not close(a, b, rt, at):
                                bad.append(f"{label}[{k}]: {a} vs {b}")
                    rep.add("regridding", case, binding, "per-pair-weights",
                            "fail" if bad else "pass",
                            "; ".join(bad[:3]) if bad else
                            f"A_ij/A_j/W_ij within rtol={wr} (areas rtol={ar})")
            elif manifest.get("blocked_upstream"):
                rep.add("regridding", case, binding, "per-pair-weights", "skip",
                        "blocked-upstream: " + manifest["blocked_upstream"])
            else:
                rep.add("regridding", case, binding, "per-pair-weights", "skip",
                        "runner does not emit per-pair arrays (no inspection "
                        "surface in this binding yet)")


def compare_reprojection(results_root, bindings, rep):
    for case_dir, manifest in manifests("reprojection"):
        case = manifest["case"]
        gold = json.loads((case_dir / manifest["golden"]).read_text())
        tol = manifest["tolerances"]
        for binding in bindings:
            if binding in manifest.get("scope_excluded", {}):
                rep.add("reprojection", case, binding, "scope", "skip",
                        manifest["scope_excluded"][binding])
                continue
            rec = case_rec(binding_results(results_root, binding, "reprojection"), case)
            if rec is None:
                rep.add("reprojection", case, binding, "runner-output", "skip",
                        "no result recorded")
                continue
            if rec["status"] != "ok":
                rep.add("reprojection", case, binding, "runner-output", "fail",
                        rec.get("message", rec["status"]))
                continue
            bad = []
            for got, ref in zip(rec["forward"], gold["forward"]):
                for axis in ("x", "y"):
                    if not close(got[axis], ref[axis],
                                 tol["forward_rtol"], tol["forward_atol_m"]):
                        bad.append(f"fwd {ref['set']}({ref['lon']},{ref['lat']}).{axis}: "
                                   f"{got[axis]} vs {ref[axis]}")
            for got, ref in zip(rec["inverse"], gold["inverse"]):
                for axis in ("lon", "lat"):
                    if not close(got[axis], ref[axis],
                                 tol["inverse_rtol"], tol["inverse_atol_deg"]):
                        bad.append(f"inv {ref['set']}({ref['x']},{ref['y']}).{axis}: "
                                   f"{got[axis]} vs {ref[axis]}")
            for got, ref in zip(rec["roundtrip"], gold["roundtrip"]):
                for axis in ("lon", "lat"):
                    if abs(got[axis + "_rt"] - ref[axis]) > tol["roundtrip_atol_deg"]:
                        bad.append(f"rt {ref['set']}({ref['lon']},{ref['lat']}).{axis}: "
                                   f"{got[axis + '_rt']} vs {ref[axis]}")
            n = len(rec["forward"]) * 2 + len(rec["inverse"]) * 2 + len(rec["roundtrip"]) * 2
            rep.add("reprojection", case, binding, "points-vs-analytic",
                    "fail" if bad else "pass",
                    "; ".join(bad[:3]) if bad else f"{n} coordinate checks within gates")


COMPARATORS = {
    "ast": compare_ast,
    "simulation": compare_simulation,
    "convergence": compare_convergence,
    "regridding": compare_regridding,
    "reprojection": compare_reprojection,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", default=str(REPO / "conformance-results"))
    ap.add_argument("--bindings", default="julia")
    ap.add_argument("--categories", default=",".join(ALL_CATEGORIES))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    results_root = Path(args.results_root)
    bindings = [b for b in args.bindings.split(",") if b]
    categories = [c for c in args.categories.split(",") if c]

    rep = Report()
    for category in categories:
        if category not in COMPARATORS:
            print(f"error: unknown category '{category}'", file=sys.stderr)
            return 2
        COMPARATORS[category](results_root, bindings, rep)

    counts = {s: sum(1 for r in rep.rows if r["status"] == s)
              for s in ("pass", "fail", "skip")}
    print(f"\n{counts['pass']} passed, {counts['fail']} failed, {counts['skip']} skipped")

    out = results_root / "comparison"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps({
        "bindings": bindings, "categories": categories,
        "counts": counts, "gates": rep.rows}, indent=2) + "\n")
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
