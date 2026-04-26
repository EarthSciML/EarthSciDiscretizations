"""Discretization-rule runtime for the Python binding.

This package mirrors ``src/rule_eval.jl`` (Julia binding): rule authors
express stencil coefficients, edge interpolants, and reconstructions as
JSON ExpressionNode ASTs (ESS spec §7 / §9.2), and the runtime exposes
a single AST entry point :func:`eval_coeff` plus a thin loader and a
stencil applicator for the ``latlon`` selector kind.

Cross-binding contract:

* :func:`eval_coeff` matches the Julia ``EarthSciDiscretizations.eval_coeff``
  passthrough (which delegates to ``EarthSciSerialization.evaluate``) for
  every op the rule catalog actually uses today: ``+ - * / ^``, the unary
  function set (``sin cos tan exp log sqrt abs``), and the constants
  ``pi`` / ``e``. Unbound variables raise :class:`UnboundVariableError`.
* :func:`load_rule` reads a rule JSON file and returns a
  :class:`Rule` with the ``applies_to``, ``grid_family``, ``combine``,
  ``accuracy``, and ``stencil`` blocks normalised.
* :func:`apply_stencil_latlon` applies a ``selector.kind == "latlon"``
  stencil to a ``(nlat, nlon)`` field with longitudinally-periodic
  neighbours and pole-aware lat-axis indexing, summing
  ``coeff * field[neighbor]`` per stencil entry under the rule's
  ``combine`` operator.

This is the seed of the cross-binding rule-evaluation contract for the
``centered_2nd_uniform_latlon`` pilot (bead ``dsc-ve1``) and is intended
to grow as additional rule families gain Python parity.
"""

from .evaluator import UnboundVariableError, eval_coeff
from .loader import Rule, StencilEntry, load_rule
from .stencil import apply_stencil_latlon

__all__ = [
    "Rule",
    "StencilEntry",
    "UnboundVariableError",
    "apply_stencil_latlon",
    "eval_coeff",
    "load_rule",
]
