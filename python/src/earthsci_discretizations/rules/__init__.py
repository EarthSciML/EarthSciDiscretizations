"""Discretization-rule runtime for the Python binding.

ESD does not carry a shadow evaluator: rule authors express stencil
coefficients, edge interpolants, and reconstructions as JSON
ExpressionNode ASTs (ESS spec §7 / §9.2), and the canonical evaluator
lives in ``earthsci_toolkit`` (ESS Python). This package exposes only:

* :func:`eval_coeff` — thin adapter over ``earthsci_toolkit.evaluate``
  for a single coefficient AST against a scalar binding map.
* :func:`load_rule` / :class:`Rule` / :class:`StencilEntry` — read a
  rule JSON file into a structured view (single-stencil or PPM-style
  multi-stencil; ESS §7.5).

Stencil application, ghost-cell synthesis, and parabola reconstruction
are owned by the canonical pipeline (``earthsci_toolkit``: ``discretize``
+ ``simulation`` / ``numpy_interpreter``); ESD does not duplicate them
here. The Julia binding's ``EarthSciDiscretizations.eval_coeff`` is the
mirror passthrough on the Julia side.
"""

from .evaluator import eval_coeff
from .loader import Rule, StencilEntry, load_rule

__all__ = [
    "Rule",
    "StencilEntry",
    "eval_coeff",
    "load_rule",
]
