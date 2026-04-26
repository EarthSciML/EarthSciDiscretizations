"""Stencil application for the ``latlon`` selector kind.

Restricted to the offset-only structured form used by
``centered_2nd_uniform_latlon``: each stencil entry has a selector
``{"kind": "latlon", "axis": "lon"|"lat", "offset": int}`` and a coefficient
AST. Applying the stencil at cell ``(j, i)`` evaluates each coefficient
under per-cell bindings and accumulates ``coeff * field[neighbor]`` under
the rule's ``combine`` operator. Longitude is treated as periodic; lat
indices outside ``[0, nlat)`` raise ``IndexError`` (the
``centered_2nd_uniform_latlon`` rule's claim is restricted to the
non-pole interior — pole handling is owned by the grid accessor's
``pole_policy``, not by the rule).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from .evaluator import eval_coeff
from .loader import Rule, StencilEntry

__all__ = ["apply_stencil_latlon"]


def _resolve_neighbor(
    selector: Mapping[str, Any], j: int, i: int, nlat: int, nlon: int
) -> tuple[int, int]:
    if selector.get("kind") != "latlon":
        raise ValueError(
            f"apply_stencil_latlon: selector kind must be 'latlon', got {selector.get('kind')!r}"
        )
    axis = selector.get("axis")
    offset = selector.get("offset", 0)
    if axis == "lon":
        return j, (i + int(offset)) % nlon
    if axis == "lat":
        jp = j + int(offset)
        if jp < 0 or jp >= nlat:
            raise IndexError(
                f"latlon lat-axis offset {offset} from j={j} leaves grid (nlat={nlat}); "
                "pole handling belongs to the grid accessor pole_policy, not the rule."
            )
        return jp, i
    raise ValueError(
        f"latlon selector axis must be 'lon' or 'lat', got {axis!r}"
    )


def _combine(values: list[float], op: str) -> float:
    if op == "+":
        return float(sum(values)) if values else 0.0
    raise ValueError(f"Unsupported combine op for latlon stencil: {op!r}")


def apply_stencil_latlon(
    rule: Rule,
    field: np.ndarray,
    j: int,
    i: int,
    bindings: Mapping[str, float],
) -> float:
    """Apply ``rule`` at cell ``(j, i)`` of ``field`` under ``bindings``.

    Parameters
    ----------
    rule:
        Rule loaded via :func:`load_rule`. Must declare ``grid_family``
        ``"latlon"`` and use only ``latlon`` selectors.
    field:
        2-D array shaped ``(nlat, nlon)``. Lon is periodic; lat is bounded.
    j, i:
        Cell index. ``i`` is wrapped modulo ``nlon`` when neighbour offsets
        cross the meridian; ``j`` plus a lat offset must stay in
        ``[0, nlat)``.
    bindings:
        Per-cell bindings consumed by the coefficient AST (e.g. ``R``,
        ``dlon``, ``dlat``, ``cos_lat``).
    """

    if rule.grid_family is not None and rule.grid_family != "latlon":
        raise ValueError(
            f"apply_stencil_latlon: rule {rule.name!r} declares grid_family"
            f" {rule.grid_family!r} (expected 'latlon')"
        )
    arr = np.asarray(field)
    if arr.ndim != 2:
        raise ValueError(f"field must be 2-D (nlat, nlon); got shape {arr.shape}")
    nlat, nlon = arr.shape

    contributions: list[float] = []
    for entry in rule.stencil:
        contributions.append(_apply_entry(entry, arr, j, i, nlat, nlon, bindings))
    return _combine(contributions, rule.combine)


def _apply_entry(
    entry: StencilEntry,
    field: np.ndarray,
    j: int,
    i: int,
    nlat: int,
    nlon: int,
    bindings: Mapping[str, float],
) -> float:
    jp, ip = _resolve_neighbor(entry.selector, j, i, nlat, nlon)
    coeff = eval_coeff(entry.coeff, bindings)
    return coeff * float(field[jp, ip])
