#!/usr/bin/env python3
"""check_convergence_order.py — assert observed order of accuracy from committed goldens.

Reads every tests/conformance/convergence/<case>/manifest.json and its
golden/errors.json (produced by scripts/regenerate-goldens.sh through the official
Julia pathway — see AGENTS.md "Goldens"). For each norm, computes the observed order

    p_k = log(e_{k-1} / e_k) / log(n_k / n_{k-1})

between consecutive resolutions and asserts

    | median(p) - expected_order | <= order_tolerance.

Runs stdlib-only (no bindings needed) so it belongs to the fast CI validate job.
Cases whose golden is missing and whose manifest says "status": "pending-golden"
are reported as SKIP; a missing golden without that marker is a failure.

golden/errors.json format:
{
  "case": "<case>", "binding": "julia", "generated_by": "scripts/regenerate-goldens.sh",
  "assert_time": 0.1,
  "errors": [ { "n": 16, "L2_error": 1.2e-4, "Linf_error": 1.9e-4 }, ... ]
}

Exit status: 0 if every case passes or is a marked SKIP; 1 otherwise.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONVERGENCE = REPO / "tests" / "conformance" / "convergence"


def median(values: list[float]) -> float:
    s = sorted(values)
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else 0.5 * (s[mid - 1] + s[mid])


def check_case(case_dir: Path) -> str:
    manifest = json.load(open(case_dir / "manifest.json"))
    case = manifest["case"]
    golden_path = case_dir / manifest.get("golden", "golden/errors.json")

    if not golden_path.exists():
        if manifest.get("status") == "pending-golden":
            print(f"SKIP {case}: golden not yet generated (status: pending-golden)")
            return "skip"
        print(f"FAIL {case}: golden {golden_path.name} missing and not marked pending-golden")
        return "fail"

    golden = json.load(open(golden_path))
    rows = sorted(golden["errors"], key=lambda r: r["n"])
    ns = [r["n"] for r in rows]
    expected_ns = [r["n"] for r in manifest["resolutions"]]
    if ns != sorted(expected_ns):
        print(f"FAIL {case}: golden resolutions {ns} != manifest resolutions {sorted(expected_ns)}")
        return "fail"

    expected = manifest["expected_order"]
    tol = manifest["order_tolerance"]
    ok = True
    for norm in manifest["norms"]:
        errors = [r[norm] for r in rows]
        if any(e <= 0 for e in errors):
            print(f"FAIL {case} [{norm}]: non-positive error in golden: {errors}")
            ok = False
            continue
        orders = [
            math.log(errors[k - 1] / errors[k]) / math.log(ns[k] / ns[k - 1])
            for k in range(1, len(errors))
        ]
        observed = median(orders)
        verdict = abs(observed - expected) <= tol
        detail = ", ".join(f"{n}:{e:.3e}" for n, e in zip(ns, errors))
        pairwise = ", ".join(f"{p:.2f}" for p in orders)
        status = "ok  " if verdict else "FAIL"
        print(
            f"{status} {case} [{norm}]: observed order {observed:.2f} "
            f"(expected {expected} ± {tol}; pairwise [{pairwise}]; errors {detail})"
        )
        ok = ok and verdict
    return "ok" if ok else "fail"


def main() -> int:
    if not CONVERGENCE.is_dir():
        print("no convergence cases found")
        return 0
    results = [check_case(d) for d in sorted(CONVERGENCE.iterdir()) if (d / "manifest.json").exists()]
    print(
        f"{results.count('ok')} passed, {results.count('fail')} failed, "
        f"{results.count('skip')} pending-golden"
    )
    return 1 if "fail" in results else 0


if __name__ == "__main__":
    sys.exit(main())
