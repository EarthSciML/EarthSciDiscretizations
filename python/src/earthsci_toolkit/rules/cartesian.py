"""Stencil application + PPM parabolic reconstruction on a 1-D periodic grid.

Mirrors the cartesian half of ``EarthSciSerialization.mms_evaluator``:

* :func:`apply_stencil_periodic_1d` — apply a single sub-stencil (or a
  bare list-form stencil) to a length-``n`` cell-averaged sample vector
  with periodic neighbour wrap. Each entry's coefficient AST is
  evaluated once against ``bindings``.
* :func:`parabola_reconstruct_periodic_1d` — Colella–Woodward (1984) PPM
  sub-cell reconstruction (eqs. 1.5/1.7/1.10): combine the named left-
  and right-edge sub-stencils into per-cell parabolas and sample them at
  the supplied ``ξ ∈ [0, 1]`` sub-cell points.
* :data:`OUTPUT_KINDS` + :func:`reference_samples` — the per-cell
  reference-sample selector used by ESS's ``output_kind`` dispatch
  (``derivative_at_cell_center`` etc.). Provided for parity so that a
  Python harness can be written against the same fixture surface as the
  Julia ``mms_convergence`` driver.

The selector kind is **not** required to be ``"cartesian"`` — the
runtime only consults ``selector.offset``. This matches the Julia path
(``apply_stencil_periodic_1d`` reads ``offset`` directly) and avoids
coupling to the rule's ``$x``/``$y`` axis placeholder syntax.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np

from .evaluator import eval_coeff
from .loader import Rule, StencilEntry

__all__ = [
    "OUTPUT_KINDS",
    "apply_stencil_periodic_1d",
    "parabola_reconstruct_periodic_1d",
    "reference_samples",
    "resolve_sub_stencil",
]


OUTPUT_KINDS: tuple[str, ...] = (
    "derivative_at_cell_center",
    "value_at_cell_center",
    "value_at_edge_left",
    "value_at_edge_right",
)


def resolve_sub_stencil(
    rule: Rule, sub_stencil: str | None
) -> tuple[StencilEntry, ...]:
    """Pick the entry list a 1-D stencil application should walk.

    Mirrors the Julia ``_resolve_substencil`` helper: a list-form rule
    rejects any non-``None`` ``sub_stencil``; a multi-stencil mapping
    requires ``sub_stencil`` and rejects unknown names.
    """

    if rule.sub_stencils is not None:
        if sub_stencil is None:
            available = sorted(rule.sub_stencils)
            raise ValueError(
                f"rule {rule.name!r} has multi-stencil mapping; caller must "
                f"select one via sub_stencil (available: {available})"
            )
        try:
            return rule.sub_stencils[sub_stencil]
        except KeyError as exc:
            available = sorted(rule.sub_stencils)
            raise ValueError(
                f"rule {rule.name!r} has no sub-stencil {sub_stencil!r} "
                f"(available: {available})"
            ) from exc
    if sub_stencil is not None:
        raise ValueError(
            f"sub_stencil={sub_stencil!r} requested but rule {rule.name!r} "
            f"carries a single stencil list"
        )
    return rule.stencil


def apply_stencil_periodic_1d(
    rule: Rule,
    u: np.ndarray | Sequence[float],
    bindings: Mapping[str, float],
    *,
    sub_stencil: str | None = None,
) -> np.ndarray:
    """Apply a 1-D stencil to ``u`` with periodic neighbour wrap.

    Each stencil entry contributes ``coeff * u[(i + offset) mod n]`` at
    every cell ``i``. Coefficients are evaluated once per call against
    ``bindings`` (no per-cell coefficient variation at this layer — that
    is the ``apply_stencil_latlon`` regime).

    Parameters
    ----------
    rule:
        Rule loaded via :func:`load_rule`. Either a list-form rule or a
        multi-stencil rule selected via ``sub_stencil``.
    u:
        Length-``n`` cell-averaged sample vector. Coerced to ``float64``.
    bindings:
        Mapping of variable name to numeric value (e.g. ``{"dx": 0.1}``).
    sub_stencil:
        Required for multi-stencil rules; rejected otherwise.

    Returns
    -------
    numpy.ndarray
        Length-``n`` ``float64`` result.
    """

    entries = resolve_sub_stencil(rule, sub_stencil)
    arr = np.asarray(u, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"u must be 1-D; got shape {arr.shape}")
    n = arr.shape[0]
    if n == 0:
        return arr.copy()

    coeff_pairs: list[tuple[int, float]] = []
    for entry in entries:
        offset_raw = entry.selector.get("offset")
        if offset_raw is None:
            raise ValueError(
                f"stencil entry selector missing 'offset': {entry.selector!r}"
            )
        offset = int(offset_raw)
        coeff_pairs.append((offset, eval_coeff(entry.coeff, bindings)))

    out = np.zeros(n, dtype=np.float64)
    for offset, coeff in coeff_pairs:
        out += coeff * np.roll(arr, -offset)
    return out


def parabola_reconstruct_periodic_1d(
    rule: Rule,
    u_bar: np.ndarray | Sequence[float],
    bindings: Mapping[str, float],
    *,
    left_edge_stencil: str,
    right_edge_stencil: str,
    subcell_points: Iterable[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Sample the CW84 PPM parabola at sub-cell points in every cell.

    Inside cell ``i`` with edge values ``u_L``, ``u_R`` and average
    ``ū``, the reconstruction is ::

        u(ξ) = u_L + ξ · (Δu + u₆ · (1 − ξ)),
        Δu = u_R − u_L,
        u₆ = 6 (ū − ½ (u_L + u_R)).

    The named sub-stencils (typically ``"q_left_edge"`` and
    ``"q_right_edge"``) compute ``u_L``/``u_R`` from ``u_bar`` via
    :func:`apply_stencil_periodic_1d`. ``subcell_points`` are
    normalised positions in ``[0, 1]``; one absolute position and one
    reconstructed value are produced per (cell, ξ).

    ``bindings`` must carry ``"dx"`` and ``"domain_lo"`` so the absolute
    positions are well-defined.

    Returns
    -------
    (xs, vals):
        Two length-``n × len(subcell_points)`` arrays. ``xs`` are the
        absolute sample positions; ``vals`` are the reconstructed
        values.
    """

    if rule.sub_stencils is None:
        raise ValueError(
            f"parabola pass requires a multi-stencil rule; {rule.name!r} "
            f"carries a single stencil list"
        )
    pts = np.asarray(list(subcell_points), dtype=np.float64)
    if pts.size == 0:
        raise ValueError("subcell_points must be non-empty")
    if not np.all((pts >= 0.0) & (pts <= 1.0)):
        raise ValueError(f"subcell_points must lie in [0, 1]; got {pts.tolist()!r}")
    if "dx" not in bindings:
        raise ValueError("parabola pass requires `dx` in bindings")
    if "domain_lo" not in bindings:
        raise ValueError("parabola pass requires `domain_lo` in bindings")

    u_arr = np.asarray(u_bar, dtype=np.float64)
    if u_arr.ndim != 1:
        raise ValueError(f"u_bar must be 1-D; got shape {u_arr.shape}")
    n = u_arr.shape[0]

    u_L = apply_stencil_periodic_1d(
        rule, u_arr, bindings, sub_stencil=left_edge_stencil
    )
    u_R = apply_stencil_periodic_1d(
        rule, u_arr, bindings, sub_stencil=right_edge_stencil
    )

    dx = float(bindings["dx"])
    lo = float(bindings["domain_lo"])

    da = u_R - u_L
    u6 = 6.0 * (u_arr - 0.5 * (u_L + u_R))

    cell_lo = lo + np.arange(n, dtype=np.float64) * dx
    xs = (cell_lo[:, None] + pts[None, :] * dx).reshape(-1)
    vals = (
        u_L[:, None]
        + pts[None, :] * (da[:, None] + u6[:, None] * (1.0 - pts[None, :]))
    ).reshape(-1)
    return xs, vals


def reference_samples(
    sample: Any,
    output_kind: str,
    domain_lo: float,
    dx: float,
    n: int,
    *,
    derivative: Any = None,
) -> np.ndarray:
    """Sample an analytic reference at the per-cell point implied by ``output_kind``.

    Mirrors ESS's ``_reference_samples``. ``sample`` and ``derivative``
    are callables returning floats given an ``x`` coordinate; only the
    one needed by ``output_kind`` is required.
    """

    if output_kind not in OUTPUT_KINDS:
        raise ValueError(
            f"unknown output_kind {output_kind!r}; expected one of {OUTPUT_KINDS}"
        )
    idx = np.arange(1, n + 1, dtype=np.float64)
    if output_kind == "derivative_at_cell_center":
        if derivative is None:
            raise ValueError(
                "output_kind='derivative_at_cell_center' requires a derivative callable"
            )
        xs = domain_lo + (idx - 0.5) * dx
        return np.array([derivative(x) for x in xs], dtype=np.float64)
    if output_kind == "value_at_cell_center":
        xs = domain_lo + (idx - 0.5) * dx
    elif output_kind == "value_at_edge_left":
        xs = domain_lo + (idx - 1.0) * dx
    else:  # value_at_edge_right
        xs = domain_lo + idx * dx
    return np.array([sample(x) for x in xs], dtype=np.float64)
