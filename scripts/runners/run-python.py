#!/usr/bin/env python3
"""run-python.py — the official Python conformance runner for EarthSciDiscretizations.

THIN WRAPPER (AGENTS.md §2, single pathway): every number and byte this script
emits comes from official earthsci_ast entry points driven over the one
canonical pipeline (.esm → parse → §9.7 import/metaparameter resolution →
§9.6.3 rewrite fixpoint → official runner). No evaluator, rule engine, or
numeric kernel lives here — only argument marshalling and JSON I/O.

CLI (EarthSciAST CONFORMANCE_SPEC.md §4.2):
  run-python.py --output-dir <path> \
      [--categories ast,simulation,convergence] \
      [--files <manifest.json>[,<manifest.json>…]] [--verbose]

Outputs (§4.4 layout, mirroring scripts/runners/run-julia.jl):
  <output-dir>/python_<category>_results.json    one file per category run
  <output-dir>/python_summary.json               roll-up
  <output-dir>/ast/<case>/expanded.golden.json   canonical post-lowering bytes

Category → official earthsci_ast pathway:
  ast          raw §9.7 pipeline (template_imports.resolve_template_machinery →
               lower_expression_templates), canonical bytes byte-identical to
               the Julia reference writer (sorted keys, 2-space indent, JSON3
               numeric normalization: integral floats → integers, Ryu-shortest
               float tokens)
  simulation   pde_inline_tests.run_pde_tests (§6.6/§6.6.5 inline tests over
               simulation.simulate, scipy method pinned by the manifest)
  convergence  parse.load(problem, metaparameters=…) per manifest resolution →
               simulate → evaluate_cellwise(reference) → field_reduce norms

Environment: the earthsci_ast venv ($ESS_ROOT/pkg/earthsci-ast-py/
.venv) or any interpreter with numpy+scipy; the package is imported from the
ESS checkout's src tree (the same layout its own pytest suite runs in).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

ESD_ROOT = Path(__file__).resolve().parent.parent.parent
ESS_ROOT = Path(os.environ.get("ESS_ROOT", ESD_ROOT.parent / "EarthSciAST"))
if not (ESS_ROOT / "esm-schema.json").is_file():
    sys.exit(f"error: EarthSciAST not found at '{ESS_ROOT}'; set ESS_ROOT "
             "or clone it as a sibling checkout (scripts/ess-locate.sh contract)")
sys.path.insert(0, str(ESS_ROOT / "packages" / "earthsci_ast" / "src"))

from earthsci_ast.lower_expression_templates import lower_expression_templates
from earthsci_ast.template_imports import resolve_template_machinery

CONF = ESD_ROOT / "tests" / "conformance"
ALL_CATEGORIES = ["ast", "simulation", "convergence", "regridding", "reprojection"]
INT64_MIN, INT64_MAX = -(2**63), 2**63 - 1

# ---------------------------------------------------------------------------
# Canonical golden writer — byte-for-byte the scheme the Julia reference runner
# (scripts/runners/run-julia.jl) emits: object keys sorted, arrays in order,
# 2-space indent, trailing newline, scalars exactly as JSON3 renders them
# after JSON3.read's numeric normalization (integral floats parse to Int64,
# hence 0.0 → 0; other floats print via Ryu shortest with Julia's positional/
# scientific switch at 1e-5 / 1e6).
# ---------------------------------------------------------------------------


def _norm(x):
    if isinstance(x, dict):
        return {str(k): _norm(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_norm(v) for v in x]
    if isinstance(x, float) and x.is_integer() and INT64_MIN <= x <= INT64_MAX:
        return int(x)  # JSON3.read parses integral floats as Int64 (0.0 → 0)
    return x


def julia_float_repr(v: float) -> str:
    """Render a non-integral finite float exactly as Julia's JSON3.write
    (Base.Ryu shortest digits; scientific iff decimal exponent < -4 or >= 6;
    scientific mantissa always carries a fraction digit; exponent unpadded)."""
    if v != v or v in (float("inf"), float("-inf")):
        raise ValueError(f"non-finite float {v!r} is not valid JSON")
    s = repr(float(v))
    neg = s.startswith("-")
    s = s.lstrip("-")
    if "e" in s or "E" in s:
        mant, exp = re.split(r"[eE]", s)
        e = int(exp)
    else:
        mant, e = s, 0
    ip, _, fp = mant.partition(".")
    digits = (ip + fp).lstrip("0")
    if ip != "0" and ip != "":
        e10 = e + len(ip) - 1
    else:
        lead = len(fp) - len(fp.lstrip("0"))
        e10 = e - lead - 1
    digits = digits.rstrip("0") or "0"
    sign = "-" if neg else ""
    if -4 <= e10 <= 5:  # positional
        if e10 < 0:
            return sign + "0." + "0" * (-e10 - 1) + digits
        if len(digits) <= e10 + 1:
            return sign + digits + "0" * (e10 + 1 - len(digits)) + ".0"
        return sign + digits[: e10 + 1] + "." + digits[e10 + 1:]
    frac = digits[1:] or "0"
    return sign + digits[0] + "." + frac + "e" + str(e10)


def _scalar(x) -> str:
    if x is True:
        return "true"
    if x is False:
        return "false"
    if x is None:
        return "null"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return julia_float_repr(x)
    if isinstance(x, str):
        return json.dumps(x, ensure_ascii=False)
    raise TypeError(f"unsupported scalar {type(x)}")


def _write_sorted(out: list[str], x, indent: int) -> None:
    pad, pad1 = "  " * indent, "  " * (indent + 1)
    if isinstance(x, dict):
        if not x:
            out.append("{}")
            return
        out.append("{\n")
        ks = sorted(x.keys())
        for i, k in enumerate(ks):
            out.append(pad1 + json.dumps(k, ensure_ascii=False) + ": ")
            _write_sorted(out, x[k], indent + 1)
            out.append(",\n" if i < len(ks) - 1 else "\n")
        out.append(pad + "}")
    elif isinstance(x, list):
        if not x:
            out.append("[]")
            return
        out.append("[\n")
        for i, v in enumerate(x):
            out.append(pad1)
            _write_sorted(out, v, indent + 1)
            out.append(",\n" if i < len(x) - 1 else "\n")
        out.append(pad + "]")
    else:
        out.append(_scalar(x))


def canonical_bytes(doc) -> bytes:
    out: list[str] = []
    _write_sorted(out, _norm(doc), 0)
    out.append("\n")
    return "".join(out).encode("utf-8")


# ---------------------------------------------------------------------------
# Manifest discovery (mirrors run-julia.jl).
# ---------------------------------------------------------------------------


def discover_manifests(category: str, files: list[str]):
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
    out = []
    for case_dir in sorted(root.iterdir()):
        mp = case_dir / "manifest.json"
        if mp.is_file():
            out.append((case_dir, json.loads(mp.read_text())))
    return out


def integrator_opts(manifest):
    j = manifest.get("integrators", {}).get("python")
    if j is None:
        return "RK45", 1e-10, 1e-12
    return (str(j.get("method", "RK45")), float(j.get("rtol", 1e-10)),
            float(j.get("atol", 1e-12)))


def python_excluded(manifest):
    """The manifest's own scope/blocked verdict for THIS binding, or None.

    A case whose manifest lists python under scope_excluded /
    blocked_upstream_bindings is recorded as status "skipped" (with the
    manifest's reason) instead of being run — compare-outputs.py reports the
    same SKIP from the manifest before ever reading the runner record, and the
    runner's exit code must not fail the orchestrator over a case the
    conformance contract says this binding cannot run yet."""
    for key in ("scope_excluded", "blocked_upstream_bindings"):
        reason = manifest.get(key, {}).get("python")
        if reason:
            return f"{key}: {reason}"
    return None


def assertion_dicts(results):
    return [{
        "model": r.model, "test_id": r.test_id, "assertion_idx": r.assertion_idx,
        "variable": r.variable, "time": r.time, "reduce": r.reduce,
        "expected": r.expected, "actual": r.actual,
        "rtol": r.rtol, "atol": r.atol, "passed": r.passed,
        "message": r.message} for r in results]


# ---------------------------------------------------------------------------
# ast — raw §9.7 expansion to canonical bytes + post-lowering gates.
# ---------------------------------------------------------------------------


def unlowered_violations(node, acc):
    if isinstance(node, dict):
        op = str(node.get("op", ""))
        if op == "apply_expression_template":
            acc.append("apply_expression_template")
        if op in ("grad", "div", "laplacian"):
            acc.append(f"unlowered {op}")
        if op == "D" and str(node.get("wrt", "t")) != "t":
            acc.append(f"unlowered spatial D (wrt={node['wrt']})")
        for k, v in node.items():
            if k == "metadata":  # prose may cite op names
                continue
            unlowered_violations(v, acc)
    elif isinstance(node, list):
        for v in node:
            unlowered_violations(v, acc)
    return acc


def run_ast(output_dir: Path, files, verbose):
    cases = {}
    for case_dir, manifest in discover_manifests("ast", files):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        excluded = python_excluded(manifest)
        if excluded is not None:
            rec["status"] = "skipped"
            rec["message"] = excluded
            if verbose:
                print(f"  {manifest['category']}/{case}: skipped ({excluded})")
            cases[case] = rec
            continue
        try:
            fixture = case_dir / manifest["fixture"]
            raw = json.loads(fixture.read_text())
            resolved = resolve_template_machinery(raw, str(fixture.parent))
            doc = lower_expression_templates(resolved if resolved is not None else raw)
            doc_n = _norm(doc)
            violations = []
            for model in doc_n.get("models", {}).values():
                unlowered_violations(
                    {k: v for k, v in model.items() if k != "metadata"}, violations)
            if violations:
                raise RuntimeError("post-lowering gate failed: "
                                   + ", ".join(sorted(set(violations))))
            data = canonical_bytes(doc_n)
            art_dir = output_dir / "ast" / case
            art_dir.mkdir(parents=True, exist_ok=True)
            art = art_dir / "expanded.golden.json"
            art.write_bytes(data)
            rec["artifact"] = str(art.relative_to(output_dir))
            rec["bytes"] = len(data)
            rec["golden"] = str((case_dir / manifest["golden"]).resolve()
                                .relative_to(ESD_ROOT))
        except Exception as err:  # noqa: BLE001 — recorded per-case
            rec["status"] = "error"
            rec["message"] = f"{type(err).__name__}: {err}"
        if verbose:
            print(f"  ast/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "python", "category": "ast", "cases": cases}


# ---------------------------------------------------------------------------
# simulation — the problems' inline §6.6/§6.6.5 tests via run_pde_tests.
# ---------------------------------------------------------------------------


def run_simulation(output_dir: Path, files, verbose):
    from earthsci_ast.pde_inline_tests import run_pde_tests

    cases = {}
    for case_dir, manifest in discover_manifests("simulation", files):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        excluded = python_excluded(manifest)
        if excluded is not None:
            rec["status"] = "skipped"
            rec["message"] = excluded
            if verbose:
                print(f"  {manifest['category']}/{case}: skipped ({excluded})")
            cases[case] = rec
            continue
        try:
            problem = (case_dir / manifest["problem"]).resolve()
            method, rtol, atol = integrator_opts(manifest)
            results = run_pde_tests(str(problem), model_name=manifest["model"],
                                    method=method, rtol=rtol, atol=atol)
            wanted = set(manifest["tests"])
            results = [r for r in results if r.test_id in wanted]
            if not results:
                raise RuntimeError(f"no assertions ran for tests {sorted(wanted)}")
            rec["assertions"] = assertion_dicts(results)
            rec["passed"] = all(r.passed for r in results)
            rec["status"] = "ok" if rec["passed"] else "failed"
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = f"{type(err).__name__}: {err}"
        if verbose:
            print(f"  simulation/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "python", "category": "simulation", "cases": cases}


# ---------------------------------------------------------------------------
# convergence — one load per manifest resolution entry (§9.7.6 binding site 4),
# error norms vs the problem's own §6.6.5 analytic reference at assert_time.
# ---------------------------------------------------------------------------


def run_convergence(output_dir: Path, files, verbose):
    from earthsci_ast.parse import load
    from earthsci_ast.pde_inline_tests import (
        evaluate_cellwise, field_reduce, simulate_states, state_cells)

    cases = {}
    for case_dir, manifest in discover_manifests("convergence", files):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        excluded = python_excluded(manifest)
        if excluded is not None:
            rec["status"] = "skipped"
            rec["message"] = excluded
            if verbose:
                print(f"  {manifest['category']}/{case}: skipped ({excluded})")
            cases[case] = rec
            continue
        try:
            problem = (case_dir / manifest["problem"]).resolve()
            model = manifest["model"]
            assert_time = float(manifest["assert_time"])
            method, rtol, atol = integrator_opts(manifest)
            errors = []
            for res in manifest["resolutions"]:
                n = int(res["n"])
                bindings = {str(k): int(v) for k, v in res["bindings"].items()}
                # A resolution entry MAY name its own problem file (meshes are
                # subsystem refs, which §9.7.6 cannot rebind — the MPAS
                # refinement family ships one thin problem file per level);
                # mirrors run-julia.jl's per-resolution override.
                res_problem = ((case_dir / res["problem"]).resolve()
                               if "problem" in res else problem)
                file = load(str(res_problem), metaparameters=bindings)
                m = file.models[model]
                refs = [a for t in m.tests for a in t.assertions
                        if a.time == assert_time and a.reduce == "L2_error"]
                if not refs:
                    raise RuntimeError("problem declares no L2_error assertion at "
                                       f"t={assert_time} to take the reference from")
                a = refs[0]
                sim = simulate_states(file, (0.0, assert_time),
                                      method=method, rtol=rtol, atol=atol,
                                      saveat=[assert_time])
                cells = state_cells(sim.var_map, a.variable, model)
                if not cells:
                    raise RuntimeError(f"state '{a.variable}' has no cells at n={n}")
                state = sim.states[-1]
                field = [float(state[slot]) for _, slot in cells]
                ref = evaluate_cellwise(a.reference, [c for c, _ in cells],
                                        index_sets=file.index_sets)
                row = {"n": n}
                for norm in manifest["norms"]:
                    row[norm] = field_reduce(norm, field, reference=ref)
                errors.append(row)
                if verbose:
                    print(f"  convergence/{case} n={n}: " + " ".join(
                        f"{k}={row[k]}" for k in sorted(row) if k != "n"))
            rec["assert_time"] = assert_time
            rec["errors"] = errors
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = f"{type(err).__name__}: {err}"
        if verbose:
            print(f"  convergence/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "python", "category": "convergence", "cases": cases}


# ---------------------------------------------------------------------------
# regridding — the fixture's exact-invariant inline tests, the regridded field,
# and the per-pair build-time setup arrays (A_ij / A_j / W_ij) read from the
# official BuildInspection surface (`simulate_states(...; inspect=…)` →
# `simulate` setup-array observability) — the §5.8 per-pair gates (mirrors
# run-julia.jl's run_regridding).
# ---------------------------------------------------------------------------


def _setup_array(insp, model: str, name: str):
    """Fetch one named setup array from the build inspection. Flattening
    prefixes each observed with its owning model ("Regrid.A_ij"), so try the
    qualified name, then the bare name, then a unique ".<name>" suffix match
    (identical to run-julia.jl's _setup_array)."""
    for key in (f"{model}.{name}", name):
        if key in insp.setup_arrays:
            return insp.setup_arrays[key]
    hits = [k for k in insp.setup_arrays if k.endswith("." + name)]
    if len(hits) == 1:
        return insp.setup_arrays[hits[0]]
    raise RuntimeError(f"setup array '{name}' not exposed by the build "
                       f"(have: {sorted(insp.setup_arrays)})")


def _rows(m):
    """Row-major nested float lists for JSON emission ([i][j] like the
    manifest triples)."""
    return [[float(v) for v in row] for row in m]


def run_regridding(output_dir: Path, files, verbose):
    from earthsci_ast.parse import load
    from earthsci_ast.pde_inline_tests import (
        run_pde_tests, simulate_states, state_cells)
    from earthsci_ast.simulation import BuildInspection

    cases = {}
    for case_dir, manifest in discover_manifests("regridding", files):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        excluded = python_excluded(manifest)
        if excluded is not None:
            rec["status"] = "skipped"
            rec["message"] = excluded
            if verbose:
                print(f"  {manifest['category']}/{case}: skipped ({excluded})")
            cases[case] = rec
            continue
        try:
            fixture = (case_dir / manifest["fixture"]).resolve()
            model = manifest["model"]
            results = run_pde_tests(str(fixture), model_name=model,
                                    method="LSODA", rtol=1e-12, atol=1e-14)
            rec["assertions"] = assertion_dicts(results)
            rec["passed"] = bool(results) and all(r.passed for r in results)
            # regrid_state integrates the constant regridded field from 0 over
            # [0,1], so state(1) IS the regridded field F_tgt.
            file = load(str(fixture))
            insp = BuildInspection()
            sim = simulate_states(file, (0.0, 1.0), method="LSODA",
                                  rtol=1e-12, atol=1e-14, saveat=[1.0],
                                  inspect=insp)
            for var in ("regrid_state", "pou_state", "cons_state"):
                cells = state_cells(sim.var_map, var, model)
                rec[var + "_at_1"] = [float(sim.states[-1][slot])
                                      for _, slot in cells]
            # Per-pair setup arrays (manifest §5.8 gates): the raw overlap-area
            # matrix, its filtered row-sums, and the normalized weights — the
            # build-once geometry the evaluator materializes at setup, emitted
            # verbatim (row-major [i][j], 1-based like the manifest triples).
            rec["A_ij"] = _rows(_setup_array(insp, model, "A_ij"))
            rec["A_j"] = [float(v) for v in _setup_array(insp, model, "A_j")]
            rec["W_ij"] = _rows(_setup_array(insp, model, "W_ij"))
            rec["status"] = "ok" if rec["passed"] else "failed"
            rec["blocked_upstream"] = str(manifest.get("blocked_upstream", ""))
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = f"{type(err).__name__}: {err}"
        if verbose:
            print(f"  regridding/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "python", "category": "regridding", "cases": cases}


# ---------------------------------------------------------------------------
# reprojection — evaluate the library's public templates at every golden point
# through the official pathway: a runner-BUILT invocation document (the
# manifest's `fixture: null` contract) importing the library, loaded so §9.7.3
# body composition inlines the template bodies, then `evaluate` per point
# (mirrors run-julia.jl's _reproj_wrapper_doc scheme).
# ---------------------------------------------------------------------------


def _reproj_wrapper_doc(params):
    crs = {k: float(params[k]) for k in ("lat_1", "lat_2", "lat_0", "lon_0", "R")}

    def mkapply(tpl, *args):
        return {"op": "apply_expression_template", "args": [], "name": tpl,
                "bindings": {**{a: a for a in args}, **crs}}

    # esm 1.0.0 declares exactly two variable types. The four projection
    # quantities are UNKNOWNS whose defining equations (bare-variable LHS) carry
    # the template invocations; there is no `expression` field on a variable and
    # no `observed` type any more (esm-spec §6.3.1). `observed_definitions` in
    # run_reprojection reads them back — the binding's classification API, not a
    # local re-derivation.
    observed = (
        ("fwd_x", "m", mkapply("lambert_conformal_forward_x", "lon", "lat")),
        ("fwd_y", "m", mkapply("lambert_conformal_forward_y", "lon", "lat")),
        ("inv_lon", "deg", mkapply("lambert_conformal_inverse_lon", "x", "y")),
        ("inv_lat", "deg", mkapply("lambert_conformal_inverse_lat", "x", "y")),
    )
    variables = {
        "lon": {"type": "parameter", "units": "deg", "default": 0.0},
        "lat": {"type": "parameter", "units": "deg", "default": 0.0},
        "x": {"type": "parameter", "units": "m", "default": 0.0},
        "y": {"type": "parameter", "units": "m", "default": 0.0},
    }
    for name, units, _expr in observed:
        variables[name] = {"type": "unknown", "units": units}
    equations = [{"lhs": name, "rhs": expr} for name, _units, expr in observed]
    return {"esm": "1.0.0",
            "metadata": {"name": "lambert_conformal_eval",
                         "description": "Runner-built template invocation "
                                        "(manifest fixture: null)."},
            "models": {"Reproject": {
                "expression_template_imports": [{"ref": "lambert_conformal.esm"}],
                "variables": variables, "equations": equations}}}


def run_reprojection(output_dir: Path, files, verbose):
    from earthsci_ast.classification import observed_definitions
    from earthsci_ast.numpy_interpreter import evaluate
    from earthsci_ast.parse import load

    cases = {}
    for case_dir, manifest in discover_manifests("reprojection", files):
        case = manifest["case"]
        rec = {"case": case, "status": "ok"}
        excluded = python_excluded(manifest)
        if excluded is not None:
            rec["status"] = "skipped"
            rec["message"] = excluded
            if verbose:
                print(f"  {manifest['category']}/{case}: skipped ({excluded})")
            cases[case] = rec
            continue
        try:
            lib = (case_dir / manifest["library"]).resolve()
            gold = json.loads((case_dir / manifest["golden"]).read_text())
            exprs = {}
            for setname, params in gold["parameter_sets"].items():
                doc = _reproj_wrapper_doc(params)
                file = load(json.dumps(doc), base_path=str(lib.parent))
                m = file.models["Reproject"]
                defs = observed_definitions(m)
                exprs[setname] = {k: defs[k]
                                  for k in ("fwd_x", "fwd_y", "inv_lon", "inv_lat")}
            fwd = []
            for pt in gold["forward"]:
                e = exprs[pt["set"]]
                b = {"lon": float(pt["lon"]), "lat": float(pt["lat"])}
                fwd.append({"set": pt["set"], "lon": b["lon"], "lat": b["lat"],
                            "x": evaluate(e["fwd_x"], b),
                            "y": evaluate(e["fwd_y"], b)})
            inv = []
            for pt in gold["inverse"]:
                e = exprs[pt["set"]]
                b = {"x": float(pt["x"]), "y": float(pt["y"])}
                inv.append({"set": pt["set"], "x": b["x"], "y": b["y"],
                            "lon": evaluate(e["inv_lon"], b),
                            "lat": evaluate(e["inv_lat"], b)})
            rt = []
            for pt in gold["roundtrip"]:
                e = exprs[pt["set"]]
                b = {"lon": float(pt["lon"]), "lat": float(pt["lat"])}
                xa = evaluate(e["fwd_x"], b)
                ya = evaluate(e["fwd_y"], b)
                b2 = {"x": xa, "y": ya}
                rt.append({"set": pt["set"], "lon": b["lon"], "lat": b["lat"],
                           "lon_rt": evaluate(e["inv_lon"], b2),
                           "lat_rt": evaluate(e["inv_lat"], b2)})
            rec["forward"] = fwd
            rec["inverse"] = inv
            rec["roundtrip"] = rt
        except Exception as err:  # noqa: BLE001
            rec["status"] = "error"
            rec["message"] = f"{type(err).__name__}: {err}"
        if verbose:
            print(f"  reprojection/{case}: {rec['status']}")
        cases[case] = rec
    return {"binding": "python", "category": "reprojection", "cases": cases}


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--categories", default="")
    ap.add_argument("--files", default="")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    categories = [c for c in args.categories.split(",") if c]
    files = [f for f in args.files.split(",") if f]
    if not categories:
        categories = list(ALL_CATEGORIES)
    if files:
        file_cats = {json.loads((Path(f) if os.path.isabs(f) else ESD_ROOT / f)
                                .read_text())["category"] for f in files}
        categories = [c for c in categories if c in file_cats]

    runners = {"ast": run_ast, "simulation": run_simulation,
               "convergence": run_convergence, "regridding": run_regridding,
               "reprojection": run_reprojection}
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {"binding": "python", "runner": "scripts/runners/run-python.py",
               "categories": {}}
    exit_code = 0
    for cat in categories:
        if cat not in runners:
            print(f"python/{cat}: unsupported category (skipping; this runner "
                  f"implements {ALL_CATEGORIES})")
            continue
        if args.verbose:
            print(f"category: {cat}")
        res = runners[cat](output_dir, files, args.verbose)
        (output_dir / f"python_{cat}_results.json").write_bytes(canonical_bytes(res))
        n_cases = len(res["cases"])
        n_bad = sum(1 for r in res["cases"].values()
                    if r["status"] not in ("ok", "skipped"))
        if n_bad:
            exit_code = 1
        summary["categories"][cat] = {"cases": n_cases, "failures": n_bad}
        print(f"python/{cat}: {n_cases} case(s), {n_bad} failure(s)")
    (output_dir / "python_summary.json").write_bytes(canonical_bytes(summary))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
