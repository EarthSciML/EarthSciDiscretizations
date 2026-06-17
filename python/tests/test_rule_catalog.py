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

    # esd-t4h: migrated off the legacy `stencil`/`selector`/`offset` blob to a
    # closed arrayop replacement in §4.2 ops. 1st-order upwind reads the cell's
    # own value (`$x`) and its upwind neighbour (`$x-1`) via `index` ops; the
    # one-sidedness now lives inside the replacement AST, not a selector list.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cartesian"' in content
    assert '"op": "grad"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"output_idx"' in content
    assert '"index"' in content
    # Upwind offset -1 appears inside an `index` arg as `$x + -1`.
    assert '"$x"' in content
    assert "-1" in content
    # No legacy stencil coefficient blobs / off-spec selector fields.
    assert '"stencil"' not in content
    assert '"selector"' not in content
    assert '"offset"' not in content


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

    # esd-t4h: migrated off the legacy `stencil`/`selector`/`stagger` blob to a
    # closed replacement AST. The two-point centered vertical stencil lowers to
    #   (u[$k+1] − u[$k−1]) / (2·h)
    # built from §4.2 ops; the per-axis index `$k` and its ±1 offsets live
    # inside the `index` op args rather than a selector list.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"vertical"' in content
    assert '"op": "grad"' in content
    assert '"replacement"' in content
    assert '"index"' in content
    assert '"$k"' in content
    assert '"h"' in content
    # No legacy stencil coefficient blobs / off-spec selector/stagger fields.
    assert '"stencil"' not in content
    assert '"selector"' not in content
    assert '"stagger"' not in content
    assert '"offset"' not in content


def test_centered_2nd_uniform_latlon_well_formed(catalog):
    rule = catalog["centered_2nd_uniform_latlon"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    # esd-t4h: migrated off the legacy `stencil`/`selector`/`offset` blob to a
    # closed replacement AST that sums the four directional terms. Latlon FD
    # reads neighbours via `index` ops over the literal axes "lat" / "lon" (per
    # SELECTOR_KINDS.md decision #10) with ±1 offsets baked into the index args,
    # so each term can carry its own metric coefficient.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"latlon"' in content
    assert '"op": "grad"' in content
    assert '"replacement"' in content
    assert '"index"' in content
    assert '"lon"' in content
    assert '"lat"' in content
    # Coefficient symbols: angular spacings dlon/dlat, sphere radius R, and
    # the latitude metric cos_lat (lon-axis only) per decisions #11 and #12.
    assert '"dlon"' in content
    assert '"dlat"' in content
    assert '"R"' in content
    assert '"cos_lat"' in content
    # No legacy stencil coefficient blobs / off-spec selector/offset fields.
    assert '"stencil"' not in content
    assert '"selector"' not in content
    assert '"offset"' not in content


def test_nn_diffusion_mpas_well_formed(catalog):
    rule = catalog["nn_diffusion_mpas"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    # esd-t4h: migrated off the legacy two-row `stencil`/`selector` blob to a
    # single closed arrayop replacement that reduces over the cell's neighbour
    # ring: ∇²u(c) = Σ_k w_k (u(n_k) − u(c)). The variable-valence reduction is
    # now an `arrayop` with `reduce: "+"` and a `ranges` bound on k; neighbours
    # and metrics are fetched via `index` ops over the connectivity tables.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"unstructured"' in content
    # Operator class: scalar Laplacian acting at cell centers.
    assert '"op": "laplacian"' in content
    assert '"emits_location": "cell_center"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"reduce": "+"' in content
    assert '"ranges"' in content
    assert '"index"' in content
    # Connectivity tables and coefficient symbols come from the dsc-7j0 MPAS
    # accessor runtime (SELECTOR_KINDS.md decision #9, snake_case): the
    # neighbour/edge ring tables plus the per-edge/cell metrics.
    assert '"cells_on_cell"' in content
    assert '"dv_edge"' in content
    assert '"dc_edge"' in content
    assert '"area_cell"' in content
    assert '"edges_on_cell"' in content
    assert '"n_edges_on_cell"' in content
    # No legacy stencil coefficient blobs / off-spec selector fields.
    assert '"stencil"' not in content
    assert '"selector"' not in content


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

# esd-t4h migrated the Fornberg-seeded uniform FD rules off the legacy
# `stencil` blob (with its `_generated_by`/`fornberg_gen.py` provenance header)
# to closed arrayop replacements in §4.2 ops. The high-order weights now live
# inside the replacement `expr` (coefficients × `index` reads over a denominator
# `N·dx`); the generated-header marker no longer exists. These tests therefore
# assert the migrated arrayop shape and accuracy marker.


def test_centered_4th_uniform_well_formed_arrayop(catalog):
    """centered_4th_uniform: O(dx^4) grad lowered to a closed arrayop replacement."""
    rule = catalog["centered_4th_uniform"]
    content = _read(rule)
    assert '"op": "grad"' in content
    assert '"accuracy": "O(dx^4)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    assert '"dx"' in content
    assert '"stencil"' not in content


def test_centered_6th_uniform_well_formed_arrayop(catalog):
    """centered_6th_uniform: O(dx^6) grad lowered to a closed arrayop replacement."""
    rule = catalog["centered_6th_uniform"]
    content = _read(rule)
    assert '"op": "grad"' in content
    assert '"accuracy": "O(dx^6)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    assert '"dx"' in content
    assert '"stencil"' not in content


def test_centered_2nd_deriv_uniform_well_formed_arrayop(catalog):
    """centered_2nd_deriv_uniform: O(dx^2) d2 lowered to a closed arrayop replacement."""
    rule = catalog["centered_2nd_deriv_uniform"]
    content = _read(rule)
    assert '"op": "d2"' in content
    assert '"accuracy": "O(dx^2)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    assert '"dx"' in content
    assert '"stencil"' not in content


def test_centered_8th_uniform_well_formed(catalog):
    """centered_8th_uniform: O(dx^8) grad lowered to a closed arrayop replacement."""
    rule = catalog["centered_8th_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    # esd-t4h: the 8th-order grad weights now live inside the replacement `expr`
    # rather than a `stencil` list. Each Fornberg pair reads u at ±1..±4 via
    # `index` ops scaled by [3, 32, 168, 672] over a denominator of 840·dx.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cartesian"' in content
    assert '"op": "grad"' in content
    assert '"accuracy": "O(dx^8)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"output_idx"' in content
    assert '"index"' in content
    # The Fornberg [3, -32, 168, -672] numerators and the 840*dx denominator are
    # the load-bearing coefficients; the ±1..±4 reads are encoded as index args.
    for coeff in ("3", "32", "168", "672"):
        assert coeff in content
    assert "840" in content
    assert '"dx"' in content
    # No legacy stencil/order blob; the offsets live in `index` op args.
    assert '"stencil"' not in content
    assert '"order"' not in content
    assert '"selector"' not in content


def test_centered_4th_deriv_uniform_well_formed(catalog):
    """centered_4th_deriv_uniform: O(dx^4) d2 lowered to a closed arrayop replacement."""
    rule = catalog["centered_4th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    # esd-t4h: 2nd-derivative weights now live inside the replacement `expr`.
    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"accuracy": "O(dx^4)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    # The 12*dx^2 denominator is the load-bearing coefficient.
    assert "12" in content
    assert '"dx"' in content
    assert '"stencil"' not in content
    assert '"order"' not in content


def test_centered_6th_deriv_uniform_well_formed(catalog):
    """centered_6th_deriv_uniform: O(dx^6) d2 lowered to a closed arrayop replacement."""
    rule = catalog["centered_6th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"accuracy": "O(dx^6)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    # The 180*dx^2 denominator is the load-bearing coefficient.
    assert "180" in content
    assert '"dx"' in content
    assert '"stencil"' not in content
    assert '"order"' not in content


def test_centered_8th_deriv_uniform_well_formed(catalog):
    """centered_8th_deriv_uniform: O(dx^8) d2 lowered to a closed arrayop replacement."""
    rule = catalog["centered_8th_deriv_uniform"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "d2"' in content
    assert '"cartesian"' in content
    assert '"accuracy": "O(dx^8)"' in content
    assert '"replacement"' in content
    assert '"arrayop"' in content
    assert '"index"' in content
    # The 5040*dx^2 denominator is the load-bearing coefficient.
    assert "5040" in content
    assert '"dx"' in content
    assert '"stencil"' not in content
    assert '"order"' not in content


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
