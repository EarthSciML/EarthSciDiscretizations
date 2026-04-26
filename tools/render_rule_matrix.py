#!/usr/bin/env python3
"""
render_rule_matrix — emit the rule × grid applicability matrix as a Hugo data
file (``docs/data/rule_matrix.json``).

Inputs:
  * ``discretizations/<family>/*.json`` — rule metadata (declared
    ``grid_family``, ``applies_to.op`` per rule).
  * ``discretizations/<family>/<rule>/fixtures/convergence/expected.esm`` —
    convergence fixture; ``expected_min_order`` is the order shown in the
    cell. ``applicable: false`` marks a structurally-skipped fixture.
  * ``tests/conformance/grids/<grid>/rules/<rule>/fixtures.json`` and
    ``tests/conformance/rules/<rule>/fixtures.json`` — canonical
    cross-binding fixtures.

Output: ``docs/data/rule_matrix.json``. The Hugo ``matrix`` section reads it
via ``site.Data.rule_matrix`` and renders the table with the layout under
``docs/layouts/matrix/list.html``.

Run from the repo root::

    python3 tools/render_rule_matrix.py

CI invokes this immediately before ``hugo --source docs --minify``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]

# Column order — matches dsc-599 acceptance and ALL_GRID_FAMILIES in
# tools/render_doc_plots.py.
GRID_FAMILIES: tuple[str, ...] = (
    "cartesian",
    "latlon",
    "vertical",
    "cubed_sphere",
    "mpas",
    "arakawa",
    "duo",
)

# Synonyms used inside rule metadata that should map to a column above.
GRID_FAMILY_ALIASES: dict[str, str] = {
    "unstructured": "mpas",  # nn_diffusion_mpas declares grid_family=unstructured
    "lat_lon": "latlon",
}


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)


def _normalize_family(name: str) -> str:
    return GRID_FAMILY_ALIASES.get(name, name)


def _walk_rules() -> list[dict[str, Any]]:
    """Walk ``discretizations/{finite_difference,finite_volume,spectral}/*.json``
    and return one record per declared rule."""
    rules: list[dict[str, Any]] = []
    for family_dir in sorted((REPO_ROOT / "discretizations").iterdir()):
        if not family_dir.is_dir() or family_dir.name == "grids":
            continue
        for rule_path in sorted(family_dir.glob("*.json")):
            doc = _load_json(rule_path)
            # Two shapes appear in the catalog: the canonical "discretizations"
            # block (one rule per file) and the "rules" block used by
            # periodic_bc.json. Both are surfaced.
            blocks = doc.get("discretizations") or doc.get("rules") or {}
            for name, body in blocks.items():
                grid_family = _normalize_family(body.get("grid_family", ""))
                if not grid_family:
                    # Boundary-condition rules (e.g. periodic_wrap_x) live in
                    # the "rules" block without a declared grid_family. They
                    # are not part of rule × grid dispatch, so omit them from
                    # the matrix.
                    continue
                applies_to = body.get("applies_to") or {}
                rules.append(
                    {
                        "name": name,
                        "family": family_dir.name,
                        "rule_path": str(rule_path.relative_to(REPO_ROOT)),
                        "grid_family": grid_family,
                        "applies_to_op": applies_to.get("op"),
                        "accuracy": body.get("accuracy"),
                    }
                )
    rules.sort(key=lambda r: (r["family"], r["name"]))
    return rules


def _convergence_status(family: str, rule: str) -> dict[str, Any]:
    """Return the convergence fixture summary for a rule, if any."""
    expected = (
        REPO_ROOT
        / "discretizations"
        / family
        / rule
        / "fixtures"
        / "convergence"
        / "expected.esm"
    )
    if not expected.is_file():
        return {"present": False}
    body = _load_json(expected)
    if body.get("applicable") is False:
        return {
            "present": True,
            "running": False,
            "skip_reason": body.get("skip_reason"),
        }
    order = body.get("expected_min_order")
    if order is None:
        return {"present": True, "running": False}
    return {
        "present": True,
        "running": True,
        "expected_min_order": order,
        "metric": body.get("metric"),
    }


def _has_canonical_fixture(rule: str, grid: str) -> bool:
    """True if a cross-binding canonical fixture exists for ``rule × grid``."""
    candidates = [
        REPO_ROOT / "tests" / "conformance" / "grids" / grid / "rules" / rule / "fixtures.json",
        REPO_ROOT / "tests" / "conformance" / "rules" / rule / "fixtures.json",
    ]
    return any(p.is_file() for p in candidates)


def _build_cell(rule: dict[str, Any], grid: str, conv: dict[str, Any]) -> dict[str, Any]:
    """Return the matrix cell for one ``rule × grid`` combination."""
    applicable = rule["grid_family"] == grid
    if not applicable:
        return {
            "applicable": False,
            "fixture": "n/a",
            "order": None,
        }

    canonical = _has_canonical_fixture(rule["name"], grid)
    if conv.get("running"):
        return {
            "applicable": True,
            "fixture": "convergence",
            "order": conv.get("expected_min_order"),
            "canonical_fixture": canonical,
        }
    if canonical:
        return {
            "applicable": True,
            "fixture": "canonical",
            "order": None,
        }
    if conv.get("present"):
        # Convergence fixture exists but is structurally skipped (pending
        # harness extension). Still a hole in dispatch coverage.
        return {
            "applicable": True,
            "fixture": "missing",
            "order": None,
            "skip_reason": conv.get("skip_reason"),
        }
    return {
        "applicable": True,
        "fixture": "missing",
        "order": None,
    }


def build_matrix() -> dict[str, Any]:
    rules_meta = _walk_rules()
    out_rules: list[dict[str, Any]] = []
    coverage = {
        "rules_total": 0,
        "applicable_cells": 0,
        "convergence_cells": 0,
        "canonical_cells": 0,
        "missing_cells": 0,
    }
    for rule in rules_meta:
        conv = _convergence_status(rule["family"], rule["name"])
        cells: dict[str, dict[str, Any]] = {}
        for grid in GRID_FAMILIES:
            cell = _build_cell(rule, grid, conv)
            cells[grid] = cell
            if cell["applicable"]:
                coverage["applicable_cells"] += 1
                if cell["fixture"] == "convergence":
                    coverage["convergence_cells"] += 1
                elif cell["fixture"] == "canonical":
                    coverage["canonical_cells"] += 1
                elif cell["fixture"] == "missing":
                    coverage["missing_cells"] += 1
        out_rules.append({**rule, "cells": cells})
        coverage["rules_total"] += 1
    return {
        "schema_version": 1,
        "generated_at": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "grid_families": list(GRID_FAMILIES),
        "rules": out_rules,
        "coverage": coverage,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "docs" / "data" / "rule_matrix.json",
        help="output path for the Hugo data file",
    )
    args = parser.parse_args(argv)

    matrix = build_matrix()
    rendered = json.dumps(matrix, indent=2, sort_keys=False) + "\n"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(rendered)
    print(f"wrote {args.out.relative_to(REPO_ROOT)} "
          f"({len(matrix['rules'])} rules × {len(matrix['grid_families'])} grids)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
