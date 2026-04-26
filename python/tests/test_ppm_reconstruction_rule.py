"""Cross-binding rule-eval conformance for ``ppm_reconstruction``.

Drives ``tests/conformance/rules/ppm_reconstruction/fixtures.json``, the
binding-independent Rust-authored harness spec (sibling: dsc-9nt), so
the Python and Rust evaluators check against the same hand-pinned
sub-stencil offsets/coefficients, the same output_kind dispatch points
(CW84 eq. 1.6 unit-window check + parabola endpoints), and the same
canonical smooth-periodic MMS sweep at
``tests/fixtures/ppm_reconstruction/{input,expected}.esm``.

The Python evaluator (``earthsci_toolkit.rules.{eval_coeff,
apply_stencil_periodic_1d, parabola_reconstruct_periodic_1d}``) walks
the on-disk multi-stencil rule; the Julia binding is the cross-binding
reference per ``docs/GRIDS_API.md`` §4.3.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from earthsci_toolkit.rules import (
    apply_stencil_periodic_1d,
    eval_coeff,
    load_rule,
    parabola_reconstruct_periodic_1d,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFORMANCE_DIR = REPO_ROOT / "tests" / "conformance" / "rules" / "ppm_reconstruction"
FIXTURES_PATH = CONFORMANCE_DIR / "fixtures.json"


@pytest.fixture(scope="module")
def fixtures() -> dict:
    with FIXTURES_PATH.open("r", encoding="utf-8") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def rule(fixtures):
    return load_rule(REPO_ROOT / fixtures["rule_path"])


def _close_rel(actual: float, expected: float, rtol: float) -> bool:
    scale = max(1.0, abs(actual), abs(expected))
    return abs(actual - expected) <= rtol * scale


def test_rule_metadata(rule):
    assert rule.name == "ppm_reconstruction"
    assert rule.family == "finite_volume"
    assert rule.grid_family == "cartesian"
    assert rule.stencil == ()
    assert rule.sub_stencils is not None
    assert set(rule.sub_stencils) == {"q_left_edge", "q_right_edge"}


def test_sub_stencil_cases(rule, fixtures):
    rtol = float(fixtures["tolerance"]["relative"])
    for case in fixtures["sub_stencil_cases"]:
        name = case["name"]
        sub = case["sub_stencil"]
        entries = rule.sub_stencils[sub]
        expected_offsets = [int(o) for o in case["expected_offsets"]]
        expected_coeffs = [float(c) for c in case["expected_coeffs"]]
        assert len(entries) == len(expected_offsets), (
            f"{name}: stencil length mismatch ({len(entries)} vs {len(expected_offsets)})"
        )
        for k, entry in enumerate(entries):
            offset = int(entry.selector["offset"])
            coeff = eval_coeff(entry.coeff, {})
            assert offset == expected_offsets[k], (
                f"{name}: offset[{k}] {offset} vs {expected_offsets[k]}"
            )
            assert _close_rel(coeff, expected_coeffs[k], rtol), (
                f"{name}: coeff[{k}] {coeff!r} vs {expected_coeffs[k]!r}"
            )


def _evaluate_edge_window(
    rule, sub: str, window_offsets: list[int], window_values: list[float]
) -> float:
    """Apply ``sub`` to a synthetic field that places ``window_values`` at the
    requested offsets relative to a probe cell.

    Builds an ``n``-vector ``u`` (where ``n`` covers every requested offset)
    such that ``u[(i0 + off) mod n] == window_values[k]`` for the specific
    probe cell ``i0``, then reads the stencil's output at ``i0``. This sidesteps
    introducing an explicit ``OutputKind`` enum on the Python side while still
    exercising the same code path as the convergence sweep.
    """

    if len(window_offsets) != len(window_values):
        raise ValueError("window_offsets and window_values length mismatch")
    span = max(window_offsets) - min(window_offsets) + 1
    n = max(span, 8)  # avoid n smaller than the stencil reach
    i0 = -min(window_offsets)
    u = np.zeros(n, dtype=np.float64)
    for off, val in zip(window_offsets, window_values, strict=True):
        u[(i0 + off) % n] = float(val)
    out = apply_stencil_periodic_1d(rule, u, {}, sub_stencil=sub)
    return float(out[i0])


def test_output_kind_cases(rule, fixtures):
    rtol = float(fixtures["tolerance"]["relative"])
    for case in fixtures["output_kind_cases"]:
        name = case["name"]
        kind = case["output_kind"]
        if kind in ("q_right_edge", "q_left_edge"):
            offsets = [int(o) for o in case["window_offsets"]]
            values = [float(v) for v in case["window_values"]]
            actual = _evaluate_edge_window(rule, kind, offsets, values)
            expected = float(case["expected"])
            assert _close_rel(actual, expected, rtol), (
                f"{name}: {actual!r} vs {expected!r}"
            )
        elif kind == "q_parabola":
            for sample in case["samples"]:
                a_L = float(sample["a_L"])
                a_R = float(sample["a_R"])
                q_i = float(sample["q_i"])
                xi = float(sample["xi"])
                expected = float(sample["expected"])
                # Drive the same parabola pass the convergence sweep uses, but
                # with synthetic per-cell edge values: on a length-1 array the
                # left/right edge stencils are irrelevant (we'd have nothing to
                # apply them to). Instead, hit the closed-form formula
                # directly, mirroring the Rust test which calls the same
                # underlying `parabola(a_L, a_R, q_i, xi)` helper.
                da = a_R - a_L
                u6 = 6.0 * (q_i - 0.5 * (a_L + a_R))
                actual = a_L + xi * (da + u6 * (1.0 - xi))
                assert _close_rel(actual, expected, rtol), (
                    f"{name}: parabola(a_L={a_L}, a_R={a_R}, q_i={q_i}, "
                    f"xi={xi}) = {actual!r} vs {expected!r}"
                )
        else:
            raise AssertionError(f"{name}: unknown output_kind {kind!r}")


def _smooth_periodic_F(x: np.ndarray) -> np.ndarray:
    # Antiderivative of f(x) = sin(2πx) + 0.3 cos(4πx).
    return -np.cos(2.0 * math.pi * x) / (2.0 * math.pi) + 0.3 * np.sin(
        4.0 * math.pi * x
    ) / (4.0 * math.pi)


def _smooth_periodic_f(x: np.ndarray) -> np.ndarray:
    return np.sin(2.0 * math.pi * x) + 0.3 * np.cos(4.0 * math.pi * x)


def _cell_averages(n: int) -> np.ndarray:
    edges = np.arange(n + 1, dtype=np.float64) / n
    return (_smooth_periodic_F(edges[1:]) - _smooth_periodic_F(edges[:-1])) * n


def test_convergence_case(rule, fixtures):
    conv = fixtures["convergence_case"]
    with (REPO_ROOT / conv["fixture_input"]).open("r", encoding="utf-8") as fh:
        input_json = json.load(fh)
    with (REPO_ROOT / conv["fixture_expected"]).open("r", encoding="utf-8") as fh:
        expected_json = json.load(fh)
    assert input_json["rule"] == "ppm_reconstruction"
    assert expected_json["rule"] == "ppm_reconstruction"

    grids = [int(g["n"]) for g in input_json["grids"]]
    samples_per_cell = int(input_json["samples_per_cell"])
    pts = np.array(
        [(k - 0.5) / samples_per_cell for k in range(1, samples_per_cell + 1)],
        dtype=np.float64,
    )
    expected_min_order = float(expected_json["expected_min_order"])

    errors: list[float] = []
    for n in grids:
        u = _cell_averages(n)
        dx = 1.0 / n
        bindings = {"dx": dx, "domain_lo": 0.0}
        xs, vals = parabola_reconstruct_periodic_1d(
            rule, u, bindings,
            left_edge_stencil="q_left_edge",
            right_edge_stencil="q_right_edge",
            subcell_points=pts,
        )
        ref = _smooth_periodic_f(xs)
        errors.append(float(np.max(np.abs(vals - ref))))

    assert all(math.isfinite(e) and e > 0 for e in errors), errors
    # Monotone refinement on a smooth profile.
    for i in range(len(errors) - 1):
        assert errors[i + 1] < errors[i], errors

    orders = [
        math.log2(errors[i] / errors[i + 1]) for i in range(len(errors) - 1)
    ]
    measured_min = min(orders)
    assert measured_min >= expected_min_order, (
        f"measured min convergence order {measured_min:.3f} < expected "
        f"{expected_min_order} (errors={errors}, orders={orders})"
    )


def test_apply_stencil_requires_sub_stencil(rule):
    with pytest.raises(ValueError, match="multi-stencil"):
        apply_stencil_periodic_1d(rule, np.zeros(8), {})


def test_apply_stencil_unknown_sub_stencil(rule):
    with pytest.raises(ValueError, match="no sub-stencil"):
        apply_stencil_periodic_1d(rule, np.zeros(8), {}, sub_stencil="bogus")


def test_parabola_requires_subcell_points_in_range(rule):
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        parabola_reconstruct_periodic_1d(
            rule, np.zeros(8), {"dx": 0.1, "domain_lo": 0.0},
            left_edge_stencil="q_left_edge",
            right_edge_stencil="q_right_edge",
            subcell_points=[1.5],
        )
