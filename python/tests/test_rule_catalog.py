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


def test_covariant_laplacian_cubed_sphere_well_formed(catalog):
    """covariant_laplacian_cubed_sphere uses cross_metric schema (§7.6, esd-p81)."""
    rule = catalog["covariant_laplacian_cubed_sphere"]
    assert rule.family == "finite_difference"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cubed_sphere"' in content
    assert '"op": "laplacian"' in content
    # Top-level rule uses cross_metric schema (§7.6) — no ad-hoc selectors.
    assert '"cross_metric"' in content
    assert '"terms"' in content
    # Per-axis stencils in the same file use stencil form.
    assert '"stencil"' in content
    assert '"axis": "xi"' in content
    assert '"axis": "eta"' in content
    # Per-axis stencils use the `selectors` array (one per axis).
    assert '"selectors"' in content
    # Cross-panel ghost handling lives in the grid accessor, NOT the selector:
    # neither the top-level rule nor its per-axis stencils carry a `panel` field.
    assert '"panel"' not in content
    # Metric-tensor components appear as cross_metric metric_component names.
    assert "ginv_xi_xi" in content
    assert "ginv_eta_eta" in content
    assert "ginv_xi_eta" in content
    assert "dJgxx_dxi" in content
    assert "dJgyy_deta" in content
    assert "dJgxe_dxi" in content
    assert "dJgxe_deta" in content


def test_ppm_edge_cubed_sphere_well_formed(catalog):
    rule = catalog["ppm_edge_cubed_sphere"]
    assert rule.family == "finite_volume"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"grid_family"' in content
    assert '"cubed_sphere"' in content
    assert '"stencil"' in content
    assert '"op": "reconstruct_panel_edge"' in content
    # FV3 two-sided edge formula on the gnomonic equidistant cubed sphere
    # collapses to a 4-cell linear stencil at offsets {-1, 0, 1, 2} along the
    # chosen axis ($x ∈ {xi, eta}); 1D-along-axis selectors composed via a
    # `selectors: [...]` array per SELECTOR_KINDS.md decision #13.
    assert '"selectors"' in content
    assert '"axis": "$x"' in content
    assert '"offset": -1' in content
    assert '"offset": 2' in content
    # Cross-panel ghost handling lives in the cubed_sphere grid accessor; the
    # rule must NOT carry a `panel` field.
    assert '"panel"' not in content
    # FV3 eq. 6.6 monotonicity clamp is documented separately as a non-linear
    # post-processing step; the rule keeps the linear stencil only.
    assert "monotonicity_constraint" in content


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

    # FV3 D-grid wind-field rules (dsc-247): vorticity (Stokes' theorem),
    # corner interpolation + Coriolis, kinetic energy, D→C interpolation,
    # and metric-corrected sub-grid sin(α) flux on the cubed sphere.
    fv3_rules = (
        "fv3_vorticity_cellmean",
        "fv3_vorticity_corner",
        "fv3_absolute_vorticity_cellmean",
        "fv3_kinetic_energy_cell",
        "fv3_d_to_c_xi",
        "fv3_d_to_c_eta",
        "fv3_sinsg_flux_xi",
        "fv3_sinsg_flux_eta",
    )
    for name in fv3_rules:
        assert name in catalog, f"missing rule {name!r}"
        assert catalog[name].family == "finite_volume"


def test_fv3_vorticity_cellmean_well_formed(catalog):
    rule = catalog["fv3_vorticity_cellmean"]
    assert rule.family == "finite_volume"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"applies_to"' in content
    assert '"op": "curl"' in content
    assert '"grid_family"' in content
    assert '"cubed_sphere"' in content
    assert '"stagger": "D"' in content
    assert '"emits_location": "cell_center"' in content
    assert '"requires_locations"' in content
    # D-grid stagger selectors at the four oriented edges of the primal cell
    # (decision #17 in SELECTOR_KINDS.md): u_edge / v_edge with axis xi/eta.
    assert '"stagger": "u_edge"' in content
    assert '"stagger": "v_edge"' in content
    assert '"axis": "xi"' in content
    assert '"axis": "eta"' in content
    assert '"offset": 0' in content
    assert '"offset": 1' in content
    # Per-cell metric bindings supplied by the cubed_sphere grid accessor:
    # physical edge lengths and primal-cell area.
    assert '"dx"' in content
    assert '"dy"' in content
    assert '"area"' in content
    # Cubed-sphere stagger enum carries the four staggering symbols.
    assert '"cubed_sphere_stagger"' in content


def test_fv3_vorticity_corner_well_formed(catalog):
    rule = catalog["fv3_vorticity_corner"]
    assert rule.family == "finite_volume"
    assert rule.path.is_file()

    content = _read(rule)
    assert '"emits_location": "corner"' in content
    # 4-point average over the four cells meeting at the corner: in-panel
    # offsets in {-1, 0} along both xi and eta axes (cubed_sphere selectors
    # follow decision #13 — composed via `selectors: [...]` array per entry).
    assert '"selectors"' in content
    assert '"offset": -1' in content
    assert '"offset": 0' in content
    # Cross-panel ghost lives in the accessor.
    assert '"panel"' not in content


def test_fv3_absolute_vorticity_cellmean_well_formed(catalog):
    rule = catalog["fv3_absolute_vorticity_cellmean"]
    assert rule.family == "finite_volume"

    content = _read(rule)
    assert '"op": "absolute_vorticity"' in content
    assert '"emits_location": "cell_center"' in content
    # ω_abs = ω_rel + 2·Omega_rot·sin(lat) — Coriolis bindings present.
    assert '"Omega_rot"' in content
    assert '"lat"' in content
    assert '"op": "sin"' in content


def test_fv3_kinetic_energy_cell_well_formed(catalog):
    rule = catalog["fv3_kinetic_energy_cell"]
    assert rule.family == "finite_volume"

    content = _read(rule)
    assert '"op": "kinetic_energy"' in content
    assert '"form": "closed_form_bilinear"' in content
    # Bilinear KE: closed-form `expression` over four named D-grid reads
    # (u_w, u_e, v_s, v_n) — not a sum of linear stencil rows.
    assert '"reads"' in content
    assert '"u_w"' in content
    assert '"u_e"' in content
    assert '"v_s"' in content
    assert '"v_n"' in content
    assert '"expression"' in content
    # Sub-grid metric at cell center (super-grid position 5).
    assert '"sin_a"' in content
    assert '"cos_a"' in content


@pytest.mark.parametrize(
    "rname,axis,stagger_out",
    [
        ("fv3_d_to_c_xi", "xi", "u_edge"),
        ("fv3_d_to_c_eta", "eta", "v_edge"),
    ],
)
def test_fv3_d_to_c_well_formed(catalog, rname, axis, stagger_out):
    rule = catalog[rname]
    assert rule.family == "finite_volume"

    content = _read(rule)
    assert '"op": "d_to_c"' in content
    assert f'"dim": "{axis}"' in content
    assert f'"emits_location": "{stagger_out}"' in content
    # Two-point centered average across the staggered interface.
    assert f'"axis": "{axis}"' in content
    assert '"offset": -1' in content
    assert '"offset": 0' in content


@pytest.mark.parametrize(
    "rname,axis,stagger,dl",
    [
        ("fv3_sinsg_flux_xi", "xi", "u_edge", "dx"),
        ("fv3_sinsg_flux_eta", "eta", "v_edge", "dy"),
    ],
)
def test_fv3_sinsg_flux_well_formed(catalog, rname, axis, stagger, dl):
    rule = catalog[rname]
    assert rule.family == "finite_volume"

    content = _read(rule)
    assert '"op": "metric_flux"' in content
    assert f'"dim": "{axis}"' in content
    assert '"form": "closed_form_upwind_blended"' in content
    assert f'"emits_location": "{stagger}"' in content
    # Upwind-blended closed form: sin_pos / sin_neg precomputed bindings,
    # |v| via the abs op, multiplied by the physical edge length and dt.
    assert '"sin_pos"' in content
    assert '"sin_neg"' in content
    assert f'"{dl}"' in content
    assert '"dt"' in content
    assert '"op": "abs"' in content


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
