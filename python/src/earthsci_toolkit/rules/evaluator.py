"""ExpressionNode AST evaluator (Python port of the ESS evaluator surface).

Mirrors ``EarthSciSerialization.evaluate`` (see
``EarthSciSerialization.jl/src/expression.jl``) for the operator set used
by the discretization rule catalog. The bindings match
``src/rule_eval.jl`` in the Julia binding: ``eval_coeff(node, bindings)``
takes a JSON-decoded AST node (number / variable string / op-dict) and
returns a Python ``float``.

Op coverage (kept in sync with the Julia ``evaluate`` path):

* arithmetic: ``+ - * / ^``
* unary functions: ``sin cos tan exp log sqrt abs``
* constants (zero-arg ops): ``pi`` (alias ``π``), ``e``

Unbound variables raise :class:`UnboundVariableError`. Domain errors on
``log`` / ``sqrt`` raise :class:`ValueError` to match the ESS behaviour
where ``DomainError`` propagates out of ``evaluate``.

Ops that the catalog does not use yet (``fn``, ``const``, ``enum``, ``call``,
slope-ratio ``max``/``min`` until ESS lands them) are deliberately rejected
so that an unsupported op surfaces immediately rather than producing a
silent wrong answer.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

__all__ = ["UnboundVariableError", "eval_coeff"]


class UnboundVariableError(KeyError):
    """Raised when an AST references a variable not present in ``bindings``."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.variable_name = name

    def __str__(self) -> str:
        return f"UnboundVariableError: variable {self.variable_name!r} not found in bindings"


_UNARY_FUNCS: dict[str, Any] = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "exp": math.exp,
    "abs": abs,
}

_CONSTANTS: dict[str, float] = {
    "pi": math.pi,
    "π": math.pi,
    "e": math.e,
}


def eval_coeff(node: Any, bindings: Mapping[str, float]) -> float:
    """Evaluate a JSON-decoded ExpressionNode against scalar bindings.

    Parameters
    ----------
    node:
        ``int`` / ``float`` literal, ``str`` variable name, or ``dict``
        with ``"op"`` and ``"args"`` keys (the JSON form produced by
        ``json.load`` on a rule file's ``"coeff"`` field).
    bindings:
        Mapping of variable name to numeric value. Values are coerced to
        ``float``.

    Returns
    -------
    float
        The numeric value of the AST.
    """

    if isinstance(node, bool):
        # ``bool`` subclasses ``int`` in Python; reject explicitly to match
        # the ESS parser, which rejects boolean literals.
        raise ValueError("Boolean literal is not a valid expression node")
    if isinstance(node, (int, float)):
        return float(node)
    if isinstance(node, str):
        try:
            return float(bindings[node])
        except KeyError as exc:
            raise UnboundVariableError(node) from exc
    if isinstance(node, Mapping):
        op = node.get("op")
        if op is None:
            raise ValueError(f"AST dict node missing 'op' key: {node!r}")
        return _eval_op(str(op), node, bindings)
    raise TypeError(
        f"Invalid AST node type {type(node).__name__}; expected number, string, or mapping"
    )


def _eval_op(op: str, node: Mapping[str, Any], bindings: Mapping[str, float]) -> float:
    # Constants are the only zero-arg ops the catalog uses today; check
    # before walking ``args`` so a malformed args list doesn't shadow the
    # zero-arg branch.
    if op in _CONSTANTS:
        args = node.get("args") or []
        if len(args) != 0:
            raise ValueError(f"{op!r} constant takes no arguments, got {len(args)}")
        return _CONSTANTS[op]

    raw_args = node.get("args")
    if raw_args is None:
        raise ValueError(f"AST op {op!r} missing 'args'")
    args = [eval_coeff(a, bindings) for a in raw_args]
    n = len(args)

    if op == "+":
        if n == 1:
            return args[0]
        if n == 0:
            return 0.0
        return math.fsum(args)
    if op == "-":
        if n == 1:
            return -args[0]
        if n == 2:
            return args[0] - args[1]
        raise ValueError(f"Subtraction requires 1 or 2 arguments, got {n}")
    if op == "*":
        if n == 0:
            return 1.0
        if n == 1:
            return args[0]
        product = 1.0
        for v in args:
            product *= v
        return product
    if op == "/":
        if n != 2:
            raise ValueError(f"Division requires exactly 2 arguments, got {n}")
        if args[1] == 0.0:
            raise ZeroDivisionError("division by zero in AST evaluation")
        return args[0] / args[1]
    if op == "^":
        if n != 2:
            raise ValueError(f"Exponentiation requires exactly 2 arguments, got {n}")
        return args[0] ** args[1]
    if op in _UNARY_FUNCS:
        if n != 1:
            raise ValueError(f"{op!r} requires exactly 1 argument, got {n}")
        return _UNARY_FUNCS[op](args[0])
    if op == "log":
        if n != 1:
            raise ValueError(f"'log' requires exactly 1 argument, got {n}")
        if args[0] <= 0.0:
            raise ValueError(f"log argument must be positive, got {args[0]}")
        return math.log(args[0])
    if op == "sqrt":
        if n != 1:
            raise ValueError(f"'sqrt' requires exactly 1 argument, got {n}")
        if args[0] < 0.0:
            raise ValueError(f"sqrt argument must be non-negative, got {args[0]}")
        return math.sqrt(args[0])

    raise ValueError(f"Unsupported AST operator: {op!r}")
