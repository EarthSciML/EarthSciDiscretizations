"""Discretization-rule runtime for the Python binding.

This package mirrors ``src/rule_eval.jl`` (Julia binding): rule authors
express stencil coefficients, edge interpolants, and reconstructions as
JSON ExpressionNode ASTs (ESS spec §7 / §9.2), and the runtime exposes
a single AST entry point :func:`eval_coeff` plus a thin loader and the
applicators for the rule families currently in the catalog.

Cross-binding contract:

* :func:`eval_coeff` matches the Julia ``EarthSciDiscretizations.eval_coeff``
  passthrough (which delegates to ``EarthSciSerialization.evaluate``) for
  every op the rule catalog actually uses today: ``+ - * / ^``, the unary
  function set (``sin cos tan exp log sqrt abs``), and the constants
  ``pi`` / ``e``. Unbound variables raise :class:`UnboundVariableError`.
* :func:`load_rule` reads a rule JSON file and returns a :class:`Rule`.
  Single-stencil rules expose entries via :attr:`Rule.stencil`;
  multi-stencil rules (PPM-style; ESS §7.5) expose them via
  :attr:`Rule.sub_stencils`.
* :func:`apply_stencil_latlon` applies a ``selector.kind == "latlon"``
  stencil at a single cell with periodic-lon / bounded-lat indexing.
* :func:`apply_stencil_periodic_1d` applies a 1-D stencil (``cartesian``
  selectors) to a periodic sample vector with ``sub_stencil`` dispatch.
* :func:`parabola_reconstruct_periodic_1d` composes the named left- and
  right-edge sub-stencils of a multi-stencil rule into the
  Colella-Woodward (1984) PPM parabola and samples it at supplied
  sub-cell points -- the ``output_kind="reconstruction"`` /
  ``"parabola"`` regime exercised by ``ppm_reconstruction``.
"""

from .cartesian import (
    OUTPUT_KINDS,
    apply_stencil_periodic_1d,
    parabola_reconstruct_periodic_1d,
    reference_samples,
    resolve_sub_stencil,
)
from .evaluator import UnboundVariableError, eval_coeff
from .loader import Rule, StencilEntry, load_rule
from .stencil import apply_stencil_latlon

__all__ = [
    "OUTPUT_KINDS",
    "Rule",
    "StencilEntry",
    "UnboundVariableError",
    "apply_stencil_latlon",
    "apply_stencil_periodic_1d",
    "eval_coeff",
    "load_rule",
    "parabola_reconstruct_periodic_1d",
    "reference_samples",
    "resolve_sub_stencil",
]
