"""Coefficient evaluator (thin adapter over ``earthsci_toolkit.evaluate``).

The discretization rule catalog stores coefficients as JSON-decoded
``ExpressionNode`` ASTs (``int`` / ``float`` / variable name / dict). ESS
``earthsci_toolkit.evaluate`` consumes the typed
:class:`earthsci_toolkit.ExprNode` form, so this module converts the dict
shape to ``ExprNode`` once and delegates — there is no parallel
operator-dispatch table here. The Julia binding does the same in
``src/rule_eval.jl``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from earthsci_toolkit.esm_types import ExprNode

from earthsci_toolkit import evaluate

__all__ = ["eval_coeff"]


_EXPR_NODE_FIELDS = {f for f in ExprNode.__dataclass_fields__ if f != "args"}


def _to_expr(node: Any) -> Any:
    if isinstance(node, bool):
        # ``bool`` subclasses ``int``; reject so boolean literals don't sneak
        # through as 0/1 (mirrors ESS parser behaviour).
        raise ValueError("Boolean literal is not a valid expression node")
    if isinstance(node, (int, float, str)):
        return node
    if isinstance(node, Mapping):
        op = node.get("op")
        if not isinstance(op, str):
            raise ValueError(f"AST dict node missing 'op' key: {node!r}")
        kwargs: dict[str, Any] = {"op": op}
        raw_args = node.get("args")
        if raw_args is not None:
            kwargs["args"] = [_to_expr(a) for a in raw_args]
        for k, v in node.items():
            if k in {"op", "args"}:
                continue
            if k in _EXPR_NODE_FIELDS:
                kwargs[k] = v
        return ExprNode(**kwargs)
    raise TypeError(
        f"Invalid AST node type {type(node).__name__}; expected number, string, or mapping"
    )


def eval_coeff(node: Any, bindings: Mapping[str, float]) -> float:
    """Evaluate a JSON-decoded ExpressionNode against scalar bindings.

    Returns the numeric value of the AST as a ``float``. Unbound variables
    raise :class:`ValueError` (propagated from
    :func:`earthsci_toolkit.evaluate`).
    """

    return float(evaluate(_to_expr(node), dict(bindings)))
