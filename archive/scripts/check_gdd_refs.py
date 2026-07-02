#!/usr/bin/env python3
"""
check_gdd_refs.py — resolve {ref} targets in GridDiscretization Descriptor files.

For each GDD file (kind: grid_discretization_descriptor), resolves every {ref}
in the discretizations map to a path relative to the GDD file, verifies the
target exists and is valid JSON, and checks that the loaded file contains the
expected discretization-scheme structure (a top-level "discretizations" map
whose entries each have "applies_to" and "grid_family").

Usage:
    python3 scripts/check_gdd_refs.py [path ...]

    path: one or more GDD .json files, or directories to search recursively.
          Defaults to "discretizations/gdd/" when no arguments are given.

Exit codes:
    0  all refs resolved without error
    1  one or more errors found (printed to stderr)
"""

import json
import os
import sys


def check_rule_file(rule_path: str) -> list[str]:
    """Return a list of error strings for a rule file, empty on success."""
    errors = []
    if not os.path.isfile(rule_path):
        return [f"ref target does not exist: {rule_path}"]
    try:
        with open(rule_path) as f:
            doc = json.load(f)
    except json.JSONDecodeError as exc:
        return [f"ref target is invalid JSON ({rule_path}): {exc}"]

    # Accept §7 discretization-scheme files: {"discretizations": {"name": {...}}}
    if "discretizations" in doc:
        for name, scheme in doc["discretizations"].items():
            missing = [k for k in ("applies_to", "grid_family") if k not in scheme]
            if missing:
                errors.append(
                    f"rule '{name}' in {rule_path} is missing required §7 fields: "
                    + ", ".join(missing)
                )
        return errors

    # Also accept a bare rule object at top level (future standalone format).
    if "applies_to" in doc and "grid_family" in doc:
        return []

    errors.append(
        f"ref target {rule_path!r} has neither a 'discretizations' map nor "
        "top-level 'applies_to'/'grid_family' — not a valid §7 rule doc"
    )
    return errors


def check_gdd(gdd_path: str) -> list[str]:
    """Return error strings for a GDD file, empty on success."""
    errors = []
    try:
        with open(gdd_path) as f:
            gdd = json.load(f)
    except json.JSONDecodeError as exc:
        return [f"{gdd_path}: invalid JSON: {exc}"]

    if gdd.get("kind") != "grid_discretization_descriptor":
        return [f"{gdd_path}: not a GDD (kind != 'grid_discretization_descriptor')"]

    grids = gdd.get("grids")
    if not isinstance(grids, dict) or not grids:
        errors.append(f"{gdd_path}: missing or empty 'grids' block")

    discs = gdd.get("discretizations", {})
    if not isinstance(discs, dict):
        errors.append(f"{gdd_path}: 'discretizations' must be a dict")
        return errors

    gdd_dir = os.path.dirname(os.path.abspath(gdd_path))
    for key, value in discs.items():
        if not isinstance(value, dict):
            errors.append(f"{gdd_path}: discretizations[{key!r}] must be a dict")
            continue
        ref = value.get("ref")
        if ref is None:
            # Inline rule object — check for required fields.
            missing = [k for k in ("applies_to", "grid_family") if k not in value]
            if missing:
                errors.append(
                    f"{gdd_path}: inline rule '{key}' missing §7 fields: "
                    + ", ".join(missing)
                )
            continue

        # Resolve the ref relative to the GDD file.
        if ref.startswith(("http://", "https://")):
            # URL refs are not checked locally.
            continue
        rule_path = os.path.normpath(os.path.join(gdd_dir, ref))
        for err in check_rule_file(rule_path):
            errors.append(f"{gdd_path}: discretizations[{key!r}]: {err}")

    return errors


def collect_gdd_files(paths: list[str]) -> list[str]:
    result = []
    for p in paths:
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                for f in files:
                    if f.endswith(".gdd.json") or f.endswith(".gdd"):
                        result.append(os.path.join(root, f))
        elif os.path.isfile(p):
            result.append(p)
        else:
            print(f"warning: {p!r} does not exist, skipping", file=sys.stderr)
    return result


def main(argv: list[str]) -> int:
    search_paths = argv if argv else ["discretizations/gdd"]
    gdd_files = collect_gdd_files(search_paths)

    if not gdd_files:
        print("no GDD files found", file=sys.stderr)
        return 1

    all_errors = []
    for gdd_path in sorted(gdd_files):
        errs = check_gdd(gdd_path)
        all_errors.extend(errs)
        status = "OK" if not errs else f"FAIL ({len(errs)} error(s))"
        print(f"{gdd_path}: {status}")
        for e in errs:
            print(f"  ERROR: {e}", file=sys.stderr)

    if all_errors:
        print(f"\n{len(all_errors)} error(s) total", file=sys.stderr)
        return 1
    print(f"\n{len(gdd_files)} GDD file(s) checked, all refs resolved OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
