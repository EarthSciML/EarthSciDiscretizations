"""Rule-file loader.

Returns a structured view of a rule JSON file (ESS spec §7) without
materialising the per-cell stencil — that lives in :mod:`stencil`. The
loader keeps the original AST nodes intact so the evaluator walks the
exact same JSON the rule author wrote.
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

    ``selector`` is the raw selector dict (``{"kind": ..., "axis": ..., "offset": ...}``
    for the simple structured kinds used by ``centered_2nd_uniform_latlon``).
    ``coeff`` is the original coefficient AST — pass it to :func:`eval_coeff`
    with appropriate per-cell bindings.
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
    raw: Mapping[str, Any] = field(repr=False)


def load_rule(path: str | Path) -> Rule:
    """Load a rule JSON file and return a :class:`Rule`.

    The JSON layout matches ``discretizations/<family>/<name>.json``
    where the top level is ``{"discretizations": {<name>: <body>}}``.
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

    raw_stencil = body.get("stencil") or []
    entries: list[StencilEntry] = []
    for raw in raw_stencil:
        if not isinstance(raw, Mapping):
            raise ValueError(f"Stencil entry must be an object, got {type(raw).__name__}")
        selector = raw.get("selector")
        coeff = raw.get("coeff")
        if not isinstance(selector, Mapping):
            raise ValueError(f"Stencil entry missing 'selector' object in {p}")
        if coeff is None:
            raise ValueError(f"Stencil entry missing 'coeff' in {p}")
        entries.append(StencilEntry(selector=dict(selector), coeff=coeff))

    family_dir = p.parent.name if p.parent.name else None
    return Rule(
        name=str(name),
        family=family_dir,
        grid_family=body.get("grid_family"),
        applies_to=dict(body.get("applies_to") or {}),
        combine=str(body.get("combine") or "+"),
        accuracy=body.get("accuracy"),
        stencil=tuple(entries),
        raw=dict(body),
    )
