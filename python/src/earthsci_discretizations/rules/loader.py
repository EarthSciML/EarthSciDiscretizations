"""Rule-file loader.

Returns a structured view of a rule JSON file (ESS spec §7) without
materialising the per-cell stencil — that lives in :mod:`stencil`. The
loader keeps the original AST nodes intact so the evaluator walks the
exact same JSON the rule author wrote.

Two stencil layouts are supported (cf. ESS §7.5):

* **Single-stencil** rules (e.g. ``centered_2nd_uniform_latlon``) carry
  ``stencil`` as a list of entries. The :class:`Rule` exposes them via
  :attr:`Rule.stencil` and :attr:`Rule.sub_stencils` is ``None``.
* **Multi-stencil** rules (e.g. ``ppm_reconstruction``) carry ``stencil``
  as an object mapping sub-stencil names (``"q_left_edge"``,
  ``"q_right_edge"``, …) to entry lists. The loader stores the whole
  mapping in :attr:`Rule.sub_stencils`; :attr:`Rule.stencil` is the empty
  tuple.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["Rule", "StencilEntry", "load_rule"]


@dataclass(frozen=True)
class StencilEntry:
    """One entry from a rule's ``stencil`` array.

    ``selector`` is the raw selector dict (e.g.
    ``{"kind": "latlon", "axis": "lon", "offset": 1}`` for
    ``centered_2nd_uniform_latlon``, or
    ``{"kind": "cartesian", "axis": "$x", "offset": -2}`` for
    ``ppm_reconstruction``). ``coeff`` is the original coefficient AST —
    pass it to :func:`eval_coeff` with appropriate per-cell bindings.
    """

    selector: Mapping[str, Any]
    coeff: Any


@dataclass(frozen=True)
class Rule:
    """Parsed rule JSON.

    Attributes mirror the ESS §7 schema fields. ``raw`` retains the full
    decoded JSON so callers can read fields not yet promoted to first-class
    attributes (e.g. ``selectors`` for 2D in-panel rules) without a second
    parse.
    """

    name: str
    family: str | None
    grid_family: str | None
    applies_to: Mapping[str, Any]
    combine: str
    accuracy: str | None
    stencil: tuple[StencilEntry, ...]
    sub_stencils: Mapping[str, tuple[StencilEntry, ...]] | None
    raw: Mapping[str, Any] = field(repr=False)


def _parse_entries(raw_entries: Any, source: Path, where: str) -> tuple[StencilEntry, ...]:
    if not isinstance(raw_entries, list):
        raise ValueError(
            f"Stencil {where} must be a JSON array in {source}; "
            f"got {type(raw_entries).__name__}"
        )
    out: list[StencilEntry] = []
    for raw in raw_entries:
        if not isinstance(raw, Mapping):
            raise ValueError(
                f"Stencil entry in {where} must be an object, got {type(raw).__name__}"
            )
        selector = raw.get("selector")
        coeff = raw.get("coeff")
        if not isinstance(selector, Mapping):
            raise ValueError(f"Stencil entry missing 'selector' object in {source}")
        if coeff is None:
            raise ValueError(f"Stencil entry missing 'coeff' in {source}")
        out.append(StencilEntry(selector=dict(selector), coeff=coeff))
    return tuple(out)


def load_rule(path: str | Path) -> Rule:
    """Load a rule JSON file and return a :class:`Rule`.

    The JSON layout matches ``discretizations/<family>/<name>.json``
    where the top level is ``{"discretizations": {<name>: <body>}}``.

    Single-stencil rules expose their entries via :attr:`Rule.stencil`.
    Multi-stencil rules (PPM-style; ESS §7.5) expose every named entry
    list via :attr:`Rule.sub_stencils` and leave :attr:`Rule.stencil`
    empty — the cartesian runtime selects one entry list by name.
    """

    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    discretizations = data.get("discretizations")
    if not isinstance(discretizations, Mapping) or len(discretizations) != 1:
        raise ValueError(
            f"Rule file {p} must contain exactly one entry under 'discretizations'"
        )
    [(name, body)] = discretizations.items()
    if not isinstance(body, Mapping):
        raise ValueError(f"Rule body for {name!r} must be a JSON object")

    raw_stencil = body.get("stencil")
    stencil: tuple[StencilEntry, ...] = ()
    sub_stencils: dict[str, tuple[StencilEntry, ...]] | None = None
    if raw_stencil is None:
        stencil = ()
    elif isinstance(raw_stencil, Mapping):
        sub_stencils = {
            str(key): _parse_entries(value, p, f"sub-stencil {key!r}")
            for key, value in raw_stencil.items()
        }
    else:
        stencil = _parse_entries(raw_stencil, p, "list")

    family_dir = p.parent.name if p.parent.name else None
    return Rule(
        name=str(name),
        family=family_dir,
        grid_family=body.get("grid_family"),
        applies_to=dict(body.get("applies_to") or {}),
        combine=str(body.get("combine") or "+"),
        accuracy=body.get("accuracy"),
        stencil=stencil,
        sub_stencils=sub_stencils,
        raw=dict(body),
    )
