"""Repo-level tests for the discretizations/ rule catalog.

Mirrors `test/test_rule_catalog.jl` (Julia) so a divergence in Python ESS
behaviour for a rule shape that no Python test covers cannot land silently
(provenance: dsc-eqdu, audit dsc-ztvz F4).

These validate that the canonical rule files are discoverable in the
``discretizations/`` catalog and carry the expected schema-level markers
(§7 for schemes, §5.2 for rules; ESS discretization RFC). The end-to-end
rule-engine exercise (.esm → ESS.parse → ESS.rewrite → ESS.evaluate)
lives elsewhere; this module only checks discoverability (the JSON exists
under the expected family directory) and JSON schema markers (text-level
``occursin``-style assertions, matching the Julia version line-for-line).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG = REPO_ROOT / "discretizations"
RULE_FAMILIES = ("finite_difference", "finite_volume")


@dataclass(frozen=True)
class CatalogRule:
    """Discovered rule entry — name, family directory, source file path."""

    name: str
    family: str
    path: Path


def _load_catalog() -> dict[str, CatalogRule]:
    """Discover every top-level rule JSON under the family directories.

    Mirrors the Julia ``load_rules`` (``src/rules.jl``): name = filename
    without the ``.json`` suffix; family = parent directory name. The
    JSON contents are only inspected by the per-rule schema-marker
    assertions below — discoverability does not depend on file structure
    (``periodic_bc.json`` carries a ``"rules"`` key while scheme files
    carry ``"discretizations"``; both are valid catalog members).
    """

    rules: dict[str, CatalogRule] = {}
    for family_dir in RULE_FAMILIES:
        for path in sorted((CATALOG / family_dir).glob("*.json")):
            name = path.stem
            rules[name] = CatalogRule(name=name, family=family_dir, path=path)
    return rules


@pytest.fixture(scope="module")
def catalog() -> dict[str, CatalogRule]:
    return _load_catalog()


def _read(rule: CatalogRule) -> str:
    return rule.path.read_text(encoding="utf-8")


def test_centered_2nd_uniform_well_formed(catalog):
    rule = catalog["centered_2nd_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    # dsc-rar: canonical linear-scheme exemplar — the rule lowers grad(u, dim=x)
    # via a closed arrayop replacement in §4.2 ops (no scheme-specific
    # `stencil`/`selector`/`offset` blobs). The replacement is
    #   (u[$x+1] − u[$x−1]) / (2·dx)
    # wrapped in an arrayop over output index $x; BC handling is delegated to
    # downstream BC rules (e.g. periodic_bc) at lowering time.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cartesian"' in content
    assert '"op": "grad"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"output_idx"' in content
    assert '"index"' in content
    # Closed-form expression built from §4.2 ops only — no stencil
    # coefficient blobs and no off-spec selector/offset fields.
    assert '"stencil"' not in content
    assert '"selector"' not in content
    assert '"offset"' not in content


def test_upwind_1st_well_formed(catalog):
    rule = catalog["upwind_1st"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cartesian"' in content
    assert '"stencil"' in content
    assert '"op": "grad"' in content
    # 1st-order upwind uses a one-sided stencil: offsets -1 and 0.
    assert '"offset": -1' in content
    assert '"offset": 0' in content


def test_periodic_bc_well_formed(catalog):
    rule = catalog["periodic_bc"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    # Periodic BC is a rewrite rule (§5.2), not a scheme (§7): it carries
    # `pattern`/`where`/`replacement` rather than `applies_to`/`stencil`.
    assert '"pattern"' in content
    assert '"where"' in content
    assert '"replacement"' in content
    assert '"dim_is_periodic"' in content
    assert '"mod"' in content


def test_centered_2nd_uniform_vertical_well_formed(catalog):
    rule = catalog["centered_2nd_uniform_vertical"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"vertical"' in content
    assert '"stencil"' in content
    assert '"op": "grad"' in content
    # Vertical centered FD uses per-family selector kind and axis k.
    assert '"kind": "vertical"' in content
    assert '"axis": "$k"' in content
    # Face-staggered MMS dispatch (esm-bhv) — selectors carry the per-face
    # stagger plus an integer offset; the two-point centered stencil reads
    # the cell's own bottom and top faces (offset 0).
    assert '"stagger": "face_bottom"' in content
    assert '"stagger": "face_top"' in content
    assert '"offset": 0' in content


def test_centered_2nd_uniform_latlon_well_formed(catalog):
    rule = catalog["centered_2nd_uniform_latlon"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"latlon"' in content
    assert '"stencil"' in content
    assert '"op": "grad"' in content
    # Latlon centered FD uses literal axis values "lon" / "lat" (per
    # SELECTOR_KINDS.md decision #10) so different stencil entries can carry
    # different metrics. Both axes have offsets -1 and 1.
    assert '"kind": "latlon"' in content
    assert '"axis": "lon"' in content
    assert '"axis": "lat"' in content
    assert '"offset": -1' in content
    assert '"offset": 1' in content
    # Coefficient symbols: angular spacings dlon/dlat, sphere radius R, and
    # the latitude metric cos_lat (lon-axis only) per decisions #11 and #12.
    assert '"dlon"' in content
    assert '"dlat"' in content
    assert '"R"' in content
    assert '"cos_lat"' in content


def test_nn_diffusion_mpas_well_formed(catalog):
    rule = catalog["nn_diffusion_mpas"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"unstructured"' in content
    assert '"stencil"' in content
    # Operator class: scalar Laplacian acting at cell centers.
    assert '"op": "laplacian"' in content
    assert '"emits_location": "cell_center"' in content
    # Two-row formulation closes ∇²u(c) = Σ_k w_k (u(n_k) − u(c)) — Row 1 is a
    # variable-valence neighbor reduction (decision #6/#7 in SELECTOR_KINDS.md);
    # Row 2 is a self-targeting indirect row whose coeff is the diagonal-weight
    # arrayop sum −Σ_k w_k (decision #6 / RFC §7.2 indirect materialization).
    assert '"kind": "reduction"' in content
    assert '"kind": "indirect"' in content
    assert '"table": "cells_on_cell"' in content
    assert '"k_bound": "k"' in content
    assert '"index_expr": "$target"' in content
    assert '"arrayop"' in content
    # Coefficient symbols come from the dsc-7j0 MPAS accessor runtime
    # (SELECTOR_KINDS.md decision #9, snake_case): area_cell, dc_edge,
    # dv_edge, edges_on_cell, n_edges_on_cell.
    assert '"dv_edge"' in content
    assert '"dc_edge"' in content
    assert '"area_cell"' in content
    assert '"edges_on_cell"' in content
    assert '"n_edges_on_cell"' in content


def test_catalog_exposes_seeded_finite_difference_rules(catalog):
    # Superset assertion: the catalog has grown with grid schemas and other
    # families; only require that the canonical FD rules remain discoverable
    # under "finite_difference".
    seeded_fd = (
        "centered_2nd_uniform",
        "centered_2nd_uniform_vertical",
        "centered_2nd_uniform_latlon",
        "nn_diffusion_mpas",
        "periodic_bc",
        "upwind_1st",
    )
    for name in seeded_fd:
        assert name in catalog, f"missing rule {name!r}"
        assert catalog[name].family == "finite_difference"

    # finite_volume/ppm_reconstruction (CW84 §1) is the first FV rule.
    assert "ppm_reconstruction" in catalog
    assert catalog["ppm_reconstruction"].family == "finite_volume"

    # finite_volume/weno5_advection (Jiang-Shu 1996) — 5th-order WENO flux
    # reconstruction, upwind-biased, 1D uniform Cartesian.
    assert "weno5_advection" in catalog
    assert catalog["weno5_advection"].family == "finite_volume"


# ---------------------------------------------------------------------------
# Fornberg-generated uniform FD rules (esd-ec4)
# ---------------------------------------------------------------------------

def test_centered_4th_uniform_generated_header(catalog):
    """centered_4th_uniform carries the fornberg_gen.py generated-file header."""
    rule = catalog["centered_4th_uniform"]
    content = _read(rule)
    assert "_generated_by" in content
    assert "fornberg_gen.py" in content


def test_centered_6th_uniform_generated_header(catalog):
    """centered_6th_uniform carries the fornberg_gen.py generated-file header."""
    rule = catalog["centered_6th_uniform"]
    content = _read(rule)
    assert "_generated_by" in content
    assert "fornberg_gen.py" in content


def test_centered_2nd_deriv_uniform_generated_header(catalog):
    """centered_2nd_deriv_uniform carries the fornberg_gen.py generated-file header."""
    rule = catalog["centered_2nd_deriv_uniform"]
    content = _read(rule)
    assert "_generated_by" in content
    assert "fornberg_gen.py" in content


def test_centered_8th_uniform_well_formed(catalog):
    """centered_8th_uniform: O(h^8) grad stencil, 8 non-zero offsets."""
    rule = catalog["centered_8th_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cartesian"' in content
    assert '"op": "grad"' in content
    assert '"order": 8' in content
    assert '"stencil"' in content
    assert '"kind": "cartesian"' in content
    # Eight non-zero offsets (−4 to −1 and +1 to +4); center weight is zero.
    for off in [-4, -3, -2, -1, 1, 2, 3, 4]:
        assert f'"offset": {off}' in content
    # Denominator 840*dx for the Fornberg [3,-32,168,-672,...] weights.
    assert "840" in content
    assert '"dx"' in content
    assert "_generated_by" in content
    assert "fornberg_gen.py" in content
    assert '"stencil"' in content
    assert '"replacement"' not in content


def test_centered_4th_deriv_uniform_well_formed(catalog):
    """centered_4th_deriv_uniform: O(h^4) d2 stencil, 5 offsets."""
    rule = catalog["centered_4th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"stencil"' in content
    for off in [-2, -1, 0, 1, 2]:
        assert f'"offset": {off}' in content
    assert "12" in content
    assert '"dx"' in content
    assert "_generated_by" in content
    assert '"replacement"' not in content


def test_centered_6th_deriv_uniform_well_formed(catalog):
    """centered_6th_deriv_uniform: O(h^6) d2 stencil, 7 offsets."""
    rule = catalog["centered_6th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"stencil"' in content
    for off in [-3, -2, -1, 0, 1, 2, 3]:
        assert f'"offset": {off}' in content
    assert "180" in content
    assert '"dx"' in content
    assert "_generated_by" in content
    assert '"replacement"' not in content


def test_centered_8th_deriv_uniform_well_formed(catalog):
    """centered_8th_deriv_uniform: O(h^8) d2 stencil, 9 offsets."""
    rule = catalog["centered_8th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"stencil"' in content
    for off in [-4, -3, -2, -1, 0, 1, 2, 3, 4]:
        assert f'"offset": {off}' in content
    assert "5040" in content
    assert '"dx"' in content
    assert "_generated_by" in content
    assert '"replacement"' not in content


def test_catalog_exposes_fornberg_generated_rules(catalog):
    """All Fornberg-generated FD rules are in the finite_difference family."""
    generated_rules = (
        "centered_4th_uniform",
        "centered_6th_uniform",
        "centered_8th_uniform",
        "centered_2nd_deriv_uniform",
        "centered_4th_deriv_uniform",
        "centered_6th_deriv_uniform",
        "centered_8th_deriv_uniform",
    )
    for name in generated_rules:
        assert name in catalog, f"missing generated rule {name!r}"
        assert catalog[name].family == "finite_difference"
