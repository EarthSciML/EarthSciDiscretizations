#!/usr/bin/env python3
"""Regenerate golden accessor outputs for the latlon conformance harness.

The Julia latlon binding is not yet on main (tracked by bead dsc-1ts); per
`docs/GRIDS_API.md` §4.3 Julia is the nominal reference binding for ULP ties,
but for the initial corpus the Python binding is used. Once the Julia binding
lands, this script should be retired (or kept as a cross-check) and the golden
regenerated from Julia.

The golden also carries reference rule-evaluator outputs under ``rule_evals``
so cross-binding rule runtimes (Rust, …) can verify ESS-AST coefficient
agreement against the same canonical values used for accessor outputs. The
reference evaluator is the ESS Python ``earthsci_toolkit.expression.evaluate``
walker (``earthsci_toolkit.parse._parse_expression`` lifts the JSON-decoded
AST into an ``ExprNode`` first).

The script activates an editable install of the Python binding so that it runs
against whatever is on the current branch. Run from the repo root:

    python tests/conformance/grids/latlon/regenerate_golden.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[3]
PYTHON_SRC = REPO_ROOT / "python" / "src"
sys.path.insert(0, str(PYTHON_SRC))

from earthsci_toolkit import grids  # noqa: E402


def _eval_ast(node, bindings: dict[str, float]) -> float:
    """Walk a JSON-decoded ESS expression-AST scalar coefficient node.

    Mirrors the closed-scalar registry implemented by the reference Julia
    evaluator (``EarthSciSerialization.evaluate(::OpExpr, ::Dict{String,Float64})``)
    and the cross-binding Rust evaluator (``earthsci_grids::eval_coeff``):
    arithmetic ``+ - * / ^``, elementary functions ``sin cos tan exp log
    sqrt abs``, the ``π``/``pi``/``e`` constants, and the ``const``
    literal carrier. ``+`` and ``*`` fold left-to-right (matching Julia's
    ``sum``/``prod`` and Rust's ``Iterator::sum``/``product``) so all
    three bindings agree at the bit-exact level on this rule's coefficient
    expressions.

    This evaluator lives inline in the regen script because the ESS Python
    distribution and the ESD Python binding share the ``earthsci_toolkit``
    package name, which makes ``from earthsci_toolkit.expression import
    evaluate`` brittle once ESD's editable tree is on ``sys.path``. The
    Rust conformance test that consumes this golden uses the production
    ``earthsci_grids::eval_coeff`` walker, so any divergence between this
    inline reference and the cross-binding evaluator surfaces immediately.
    """
    if isinstance(node, bool):  # exclude bool from numeric branch
        raise TypeError(f"boolean nodes are not supported in scalar AST eval: {node}")
    if isinstance(node, (int, float)):
        return float(node)
    if isinstance(node, str):
        if node not in bindings:
            raise ValueError(f"unbound variable: {node}")
        return float(bindings[node])
    if not isinstance(node, dict):
        raise TypeError(f"unsupported AST node: {node!r}")
    op = node["op"]
    if op == "const":
        v = node["value"]
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise TypeError(f"scalar `const` expected, got {v!r}")
        return float(v)
    if op in ("enum", "fn"):
        raise NotImplementedError(f"`{op}` op not supported by the scalar regen evaluator")
    args = [_eval_ast(a, bindings) for a in node["args"]]
    if op == "+":
        return sum(args) if args else 0.0
    if op == "-":
        if len(args) == 1:
            return -args[0]
        if len(args) == 2:
            return args[0] - args[1]
        raise ValueError(f"`-` expects 1 or 2 args, got {len(args)}")
    if op == "*":
        if not args:
            return 1.0
        result = 1.0
        for v in args:
            result *= v
        return result
    if op == "/":
        if len(args) != 2:
            raise ValueError(f"`/` expects 2 args, got {len(args)}")
        return args[0] / args[1]
    if op == "^":
        if len(args) != 2:
            raise ValueError(f"`^` expects 2 args, got {len(args)}")
        return args[0] ** args[1]
    if op in ("sin", "cos", "tan", "exp", "abs"):
        if len(args) != 1:
            raise ValueError(f"`{op}` expects 1 arg, got {len(args)}")
        return getattr(math, op)(args[0]) if op != "abs" else abs(args[0])
    if op == "log":
        if len(args) != 1:
            raise ValueError(f"`log` expects 1 arg, got {len(args)}")
        if args[0] <= 0.0:
            raise ValueError(f"log domain error: {args[0]}")
        return math.log(args[0])
    if op == "sqrt":
        if len(args) != 1:
            raise ValueError(f"`sqrt` expects 1 arg, got {len(args)}")
        if args[0] < 0.0:
            raise ValueError(f"sqrt domain error: {args[0]}")
        return math.sqrt(args[0])
    if op in ("π", "pi"):
        if args:
            raise ValueError(f"`{op}` constant takes 0 args, got {len(args)}")
        return math.pi
    if op == "e":
        if args:
            raise ValueError(f"`e` constant takes 0 args, got {len(args)}")
        return math.e
    raise ValueError(f"unsupported op: {op}")

FIXTURES_PATH = HERE / "fixtures.json"
GOLDEN_DIR = HERE / "golden"

DIR_KEYS = ("W", "E", "S", "N")
METRIC_NAMES = (
    "J",
    "g_lonlon",
    "g_latlat",
    "g_lonlat",
    "ginv_lonlon",
    "ginv_latlat",
    "ginv_lonlat",
)

# Rules to evaluate at every fixture query point. Paths are repo-relative so
# bindings can resolve them the same way regardless of working directory.
RULE_PATHS = (
    "discretizations/finite_difference/centered_2nd_uniform_latlon.json",
)


def _load_rule(rel_path: str) -> tuple[str, dict]:
    path = REPO_ROOT / rel_path
    payload = json.loads(path.read_text())
    discs = payload["discretizations"]
    if len(discs) != 1:
        raise ValueError(
            f"{rel_path}: expected exactly one entry under `discretizations`, got {sorted(discs)}"
        )
    (name,) = discs
    return name, discs[name]


def _bindings_for(grid, j: int, i: int) -> dict[str, float]:
    """Bindings consumed by the centered_2nd_uniform_latlon stencil.

    Per the rule's `discretizations` block the coefficient nodes reference
    ``R``, ``cos_lat``, ``dlon`` (per-row), and ``dlat`` (per-row).
    """
    _, lat = grid.cell_centers(j=j, i=i)
    nlon = int(grid.nlon_per_row[j])
    dlon = 2.0 * math.pi / nlon
    lat_s = float(grid.lat_edges[j])
    lat_n = float(grid.lat_edges[j + 1])
    dlat = lat_n - lat_s
    return {
        "R": float(grid.R),
        "cos_lat": math.cos(float(lat)),
        "dlon": dlon,
        "dlat": dlat,
    }


def _eval_rule(rule_def: dict, bindings: dict[str, float]) -> list[float]:
    return [float(_eval_ast(entry["coeff"], bindings)) for entry in rule_def["stencil"]]


def _rule_eval_block(grid, qps: list) -> list[dict]:
    blocks: list[dict] = []
    for rel_path in RULE_PATHS:
        name, rule_def = _load_rule(rel_path)
        per_qp_coeffs: list[list[float]] = []
        per_qp_bindings: list[dict[str, float]] = []
        for qp in qps:
            j, i = int(qp[0]), int(qp[1])
            bindings = _bindings_for(grid, j, i)
            per_qp_bindings.append(bindings)
            per_qp_coeffs.append(_eval_rule(rule_def, bindings))
        blocks.append({
            "rule": name,
            "rule_path": rel_path,
            "stencil_selectors": [entry["selector"] for entry in rule_def["stencil"]],
            "bindings_per_qp": per_qp_bindings,
            "stencil_coeffs": per_qp_coeffs,
        })
    return blocks


def _build_grid(opts: dict):
    variant = opts["variant"]
    kwargs = {
        "variant": variant,
        "R": float(opts["R"]),
        "ghosts": int(opts["ghosts"]),
        "dtype": opts["dtype"],
        "pole_policy": opts["pole_policy"],
    }
    if variant == "regular":
        kwargs["nlon"] = int(opts["nlon"])
        kwargs["nlat"] = int(opts["nlat"])
    else:
        kwargs["nlon_per_row"] = [int(x) for x in opts["nlon_per_row"]]
        if "lat_edges" in opts:
            kwargs["lat_edges"] = [float(x) for x in opts["lat_edges"]]
    return grids.lat_lon(**kwargs)


def _neighbor_entry(nbr):
    if nbr is None:
        return None
    j, i = nbr
    return [int(j), int(i)]


def build_output(fixture: dict) -> dict:
    opts = fixture["opts"]
    grid = _build_grid(opts)
    qps = fixture["query_points"]

    centers = []
    neighbors = {k: [] for k in DIR_KEYS}
    metrics = {name: [] for name in METRIC_NAMES}
    areas = []

    for qp in qps:
        j, i = int(qp[0]), int(qp[1])
        lon, lat = grid.cell_centers(j=j, i=i)
        centers.append({"lon": float(lon), "lat": float(lat)})

        nbrs = grid.neighbors(j, i)
        for dkey in DIR_KEYS:
            neighbors[dkey].append(_neighbor_entry(nbrs[dkey]))

        for name in METRIC_NAMES:
            metrics[name].append(float(grid.metric_eval(name, j, i)))
        areas.append(float(grid.metric_eval("area", j, i)))

    return {
        "fixture": fixture["name"],
        "family": "lat_lon",
        "variant": opts["variant"],
        "opts": opts,
        "n_cells": int(grid.n_cells),
        "query_points": qps,
        "cell_centers": centers,
        "neighbors_W": neighbors["W"],
        "neighbors_E": neighbors["E"],
        "neighbors_S": neighbors["S"],
        "neighbors_N": neighbors["N"],
        "metric_J": metrics["J"],
        "metric_g_lonlon": metrics["g_lonlon"],
        "metric_g_latlat": metrics["g_latlat"],
        "metric_g_lonlat": metrics["g_lonlat"],
        "metric_ginv_lonlon": metrics["ginv_lonlon"],
        "metric_ginv_latlat": metrics["ginv_latlat"],
        "metric_ginv_lonlat": metrics["ginv_lonlat"],
        "area": areas,
        "rule_evals": _rule_eval_block(grid, qps),
    }


def main() -> None:
    spec = json.loads(FIXTURES_PATH.read_text())
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    for fx in spec["fixtures"]:
        out = build_output(fx)
        path = GOLDEN_DIR / f"{fx['name']}.json"
        with path.open("w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
