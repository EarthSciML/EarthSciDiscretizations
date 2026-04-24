"""Cross-language conformance for the vertical grid family.

Every binding (python, typescript, julia, rust) must emit byte-identical
canonical .esm content for the same vertical-grid options — per
``docs/GRIDS_API.md`` §3.5 and §4.1. The vertical family has no
transcendental math (pure rationals + verbatim user-provided levels + at
most one IEEE-754 multiply-add for eta), so the §4.2 ULP fallback is not
expected to fire; ``tolerance.relative`` in fixtures.json is 0.0.

The committed ``discretizations/grids/vertical/*.esm`` files are both the
ESD fixture corpus and the conformance golden — duplicating them under
``golden/`` here would leave two sources of truth for the same declarative
config.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from earthsci_toolkit import grids

HARNESS_DIR = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "conformance"
    / "grids"
    / "vertical"
)
FIXTURES_PATH = HARNESS_DIR / "fixtures.json"
_SPEC = json.loads(FIXTURES_PATH.read_text())
GOLDEN_DIR = (HARNESS_DIR / _SPEC["golden_dir"]).resolve()
REL_TOL = float(_SPEC["tolerance"]["relative"])
_STRIP_FIELDS = tuple(_SPEC["provenance_strip_fields"])
_FIXTURE_NAMES = [f["name"] for f in _SPEC["fixtures"]]


def _fixture(name: str) -> dict:
    for f in _SPEC["fixtures"]:
        if f["name"] == name:
            return f
    raise KeyError(name)


def _canonicalize(doc: dict) -> str:
    stripped = dict(doc)
    if "provenance" in stripped:
        prov = {
            k: v for k, v in stripped["provenance"].items() if k not in _STRIP_FIELDS
        }
        stripped["provenance"] = prov
    # Match discretizations/grids/vertical/regenerate_fixtures.py:
    # json.dumps(..., sort_keys=True, indent=2) + "\n".
    return json.dumps(stripped, sort_keys=True, indent=2) + "\n"


@pytest.mark.parametrize("fixture_name", _FIXTURE_NAMES)
def test_conformance_against_golden(fixture_name: str) -> None:
    assert REL_TOL == 0.0, (
        "Vertical family is declarative-only; ULP fallback not expected. "
        f"If {REL_TOL} is intentional, update this assertion and the Rust/TS "
        "runners together (docs/GRIDS_API.md §4.2 requires review across all "
        "bindings)."
    )

    opts = _fixture(fixture_name)["opts"]
    grid = grids.vertical(**opts)
    emitted = _canonicalize(grid.to_esm())

    golden_path = GOLDEN_DIR / f"{fixture_name}.esm"
    golden_raw = json.loads(golden_path.read_text())
    expected = _canonicalize(golden_raw)

    assert emitted == expected, (
        f"python vertical binding disagrees with committed golden "
        f"({fixture_name}.esm) after stripping provenance.{_STRIP_FIELDS}"
    )


def test_all_committed_fixtures_are_covered() -> None:
    """Every committed .esm under discretizations/grids/vertical/ is in the corpus.

    If a new fixture is committed there without being added to
    ``fixtures.json``, a binding could silently diverge on it. This check
    is the tripwire.
    """
    committed = sorted(p.stem for p in GOLDEN_DIR.glob("*.esm"))
    corpus = sorted(_FIXTURE_NAMES)
    assert committed == corpus, (
        "Committed vertical fixtures diverged from the conformance corpus. "
        f"committed={committed}, corpus={corpus}. "
        "Add the new fixture to tests/conformance/grids/vertical/fixtures.json "
        "or remove the orphan .esm file."
    )
