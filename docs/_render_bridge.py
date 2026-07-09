#!/usr/bin/env python3
"""_render_bridge.py — line-oriented bridge to the official ESS display path.

Runs INSIDE the EarthSciAST earthsci_ast environment (the
pkg/earthsci-ast-py/.venv interpreter, with PYTHONPATH pointing at
pkg/earthsci-ast-py/src when the venv has no installed copy).
docs/generate_catalog.py spawns it when earthsci_ast is not importable
in-process, so the docs generator can use the official unicode renderer
without vendoring any display logic (AGENTS.md "single pathway": docs may
call official ESS display/pretty-print APIs; they never evaluate math).

Protocol: one JSON-encoded ESM expression per stdin line; one JSON object
{"ok": true, "text": ...} or {"ok": false, "error": ...} per stdout line.
"""

import json
import sys


def main() -> int:
    try:
        from earthsci_ast.display import to_unicode
        from earthsci_ast.parse import _parse_expression
    except Exception as exc:  # toolkit not importable in this interpreter
        sys.stdout.write(json.dumps({"ok": False, "error": f"import: {exc}"}) + "\n")
        sys.stdout.flush()
        return 1

    # Handshake so the parent knows the toolkit loaded.
    sys.stdout.write(json.dumps({"ok": True, "text": "ready"}) + "\n")
    sys.stdout.flush()

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            expr = json.loads(line)
            out = {"ok": True, "text": to_unicode(_parse_expression(expr))}
        except Exception as exc:  # any parse/display failure -> structural fallback
            out = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        sys.stdout.write(json.dumps(out) + "\n")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
