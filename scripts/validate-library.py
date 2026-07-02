#!/usr/bin/env python3
"""validate-library.py — schema validation + policy lint for the ESD standard library.

Checks every library file (grids/, regridding/, reprojection/, problems/) against:
  - the ESS JSON schema ($ESS_ROOT/esm-schema.json), and
  - the repo's policy lint (AGENTS.md "Metadata and tags contract"):

    L001  exactly one `esd:` kind tag (grid|stencil|rule|regrid|reproject|problem)
    L002  required prefix tags present for the kind
    L003  name agreement: filename stem == metadata.name == sole expression_templates
          key (esd:grid: file is grid.esm and metadata.name == parent directory name)
    L004  every `grid:` tag resolves to grids/<name>/grid.esm with matching metadata.name
    L005  stencils/rules live under the grid directory their `grid:` tag names
    L006  a rule's makearray regions tile the axes named by its `axes:` tag exactly
          (checked by folding metaparameter bounds at sampled sizes — no evaluation
          of any model math)
    L007  kind/payload agreement: grid|stencil|rule|regrid|reproject are pure
          template-library files; problem files carry models
    L008  esm version is 0.8.0

Usage:
  validate-library.py               validate + lint the library (exit 1 on findings)
  validate-library.py --invalid     check tests/invalid/lint/ fixtures against
                                    tests/invalid/lint/expected_errors.json

Requires: jsonschema (pip install jsonschema). CI installs it; locally the
EarthSciSerialization Python venv also carries it.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LIBRARY_DIRS = ("grids", "regridding", "reprojection", "problems")
KINDS = ("grid", "stencil", "rule", "regrid", "reproject", "problem")
REQUIRED_TAGS = {
    "grid": ["family:"],
    "stencil": ["family:", "grid:", "op:"],
    "rule": ["family:", "grid:", "op:", "order:", "bc:", "axes:"],
    "regrid": ["method:"],
    "reproject": ["crs:"],
    "problem": ["grid:"],
}
COMPONENT_KEYS = ("models", "reaction_systems", "data_loaders", "coupling", "domain")
# Sampled sizes for L006 metaparameter folding. Rules must tile for every size at or
# above the smallest stencil-supporting grid; one even and one odd sample catches
# parity mistakes in region arithmetic.
SAMPLE_SIZES = (16, 17)
MAX_COVERAGE_CELLS = 1_000_000


def ess_root() -> Path:
    out = subprocess.run(
        [str(REPO / "scripts" / "ess-locate.sh")], capture_output=True, text=True
    )
    if out.returncode != 0:
        sys.stderr.write(out.stderr)
        sys.exit(2)
    return Path(out.stdout.strip())


def load_json(path: Path):
    with open(path) as f:
        return json.load(f)


class Findings:
    def __init__(self) -> None:
        self.items: list[tuple[str, str, str]] = []  # (relpath, code, message)

    def add(self, path: Path, code: str, message: str) -> None:
        self.items.append((str(path.relative_to(REPO)), code, message))

    def report(self) -> int:
        for relpath, code, message in self.items:
            print(f"{relpath}: {code}: {message}")
        print(f"{len(self.items)} finding(s)")
        return 1 if self.items else 0


def prefix_tags(doc: dict) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for tag in doc.get("metadata", {}).get("tags", []):
        if ":" in tag:
            prefix, value = tag.split(":", 1)
            out.setdefault(prefix, []).append(value)
    return out


def fold(expr, bindings: dict[str, int]) -> int:
    """Fold a metaparameter expression to an int (spec §9.7.6 grammar)."""
    if isinstance(expr, bool) or not isinstance(expr, (int, str, dict)):
        raise ValueError(f"not a metaparameter expression: {expr!r}")
    if isinstance(expr, int):
        return expr
    if isinstance(expr, str):
        if expr not in bindings:
            raise ValueError(f"unbound metaparameter {expr!r}")
        return bindings[expr]
    op, args = expr.get("op"), [fold(a, bindings) for a in expr.get("args", [])]
    if op == "+":
        return sum(args)
    if op == "-":
        return -args[0] if len(args) == 1 else args[0] - sum(args[1:])
    if op == "*":
        r = 1
        for a in args:
            r *= a
        return r
    if op == "/":
        r = args[0]
        for a in args[1:]:
            if a == 0 or r % a != 0:
                raise ValueError(f"inexact division {r}/{a}")
            r //= a
        return r
    raise ValueError(f"bad metaparameter op {op!r}")


def find_makearray_bodies(doc: dict):
    """Yield (template_name, makearray_node) for every match rule whose body is a makearray."""
    for name, entry in (doc.get("expression_templates") or {}).items():
        if "match" in entry and isinstance(entry.get("body"), dict):
            body = entry["body"]
            if body.get("op") == "makearray":
                yield name, body


class Library:
    """The discovered library: files, parsed docs, and the grid registry."""

    def __init__(self) -> None:
        self.files: list[Path] = []
        for d in LIBRARY_DIRS:
            base = REPO / d
            if base.is_dir():
                self.files.extend(sorted(base.rglob("*.esm")))
        # Mesh-data files (grids/*/mesh/) are fixture data, not lintable library
        # entries; they are schema-validated only.
        self.entries = [p for p in self.files if "mesh" not in p.parts]
        self.docs: dict[Path, dict] = {}
        for p in self.files:
            try:
                self.docs[p] = load_json(p)
            except json.JSONDecodeError as e:
                self.docs[p] = {"__parse_error__": str(e)}
        self.grids: dict[str, Path] = {}
        for p in self.entries:
            doc = self.docs[p]
            if p.name == "grid.esm" and "__parse_error__" not in doc:
                name = doc.get("metadata", {}).get("name", "")
                self.grids[name] = p


def lint_file(path: Path, lib: Library, findings: Findings) -> None:
    doc = lib.docs[path]
    if "__parse_error__" in doc:
        findings.add(path, "PARSE", doc["__parse_error__"])
        return
    meta = doc.get("metadata", {})
    tags = prefix_tags(doc)
    stem = path.stem

    if doc.get("esm") != "0.8.0":
        findings.add(path, "L008", f"esm version {doc.get('esm')!r}, expected '0.8.0'")

    kinds = tags.get("esd", [])
    if len(kinds) != 1 or kinds[0] not in KINDS:
        findings.add(path, "L001", f"expected exactly one esd: kind tag, got {kinds}")
        return
    kind = kinds[0]

    for req in REQUIRED_TAGS[kind]:
        if req.rstrip(":") not in tags:
            findings.add(path, "L002", f"esd:{kind} requires a `{req}` tag")

    # L003 — name agreement.
    name = meta.get("name", "")
    templates = doc.get("expression_templates") or {}
    if kind == "grid":
        if path.name != "grid.esm":
            findings.add(path, "L003", "esd:grid files must be named grid.esm")
        if name != path.parent.name:
            findings.add(
                path, "L003",
                f"metadata.name {name!r} != grid directory {path.parent.name!r}",
            )
    elif kind == "problem":
        if name != stem:
            findings.add(path, "L003", f"metadata.name {name!r} != filename stem {stem!r}")
    else:
        if name != stem:
            findings.add(path, "L003", f"metadata.name {name!r} != filename stem {stem!r}")
        if list(templates) != [name]:
            findings.add(
                path, "L003",
                f"expected sole expression_templates key {name!r}, got {list(templates)}",
            )

    # L007 — kind/payload agreement.
    present = [k for k in COMPONENT_KEYS if k in doc]
    if kind == "problem":
        if "models" not in doc:
            findings.add(path, "L007", "esd:problem files must carry models")
    else:
        if present:
            findings.add(path, "L007", f"esd:{kind} must be a pure library file; found {present}")
        if not templates:
            findings.add(path, "L007", f"esd:{kind} must declare expression_templates")

    # L004 / L005 — grid resolution and placement.
    for grid_name in tags.get("grid", []):
        grid_path = lib.grids.get(grid_name)
        if grid_path is None:
            findings.add(path, "L004", f"grid:{grid_name} does not resolve to a grids/*/grid.esm")
        elif kind in ("stencil", "rule") and grid_path.parent not in path.parents:
            findings.add(
                path, "L005",
                f"esd:{kind} tagged grid:{grid_name} must live under {grid_path.parent.relative_to(REPO)}/",
            )

    if kind == "rule":
        lint_rule_regions(path, doc, tags, lib, findings)


def lint_rule_regions(path: Path, doc: dict, tags, lib: Library, findings: Findings) -> None:
    """L006 — makearray regions tile the axes named by the rule's axes: tag."""
    grid_names = tags.get("grid", [])
    axes_decl = tags.get("axes", [])
    if not grid_names or not axes_decl:
        return  # L002 already fired
    grid_path = lib.grids.get(grid_names[0])
    if grid_path is None:
        return  # L004 already fired
    grid = lib.docs[grid_path]
    axes = axes_decl[0].split(",") if axes_decl else []
    index_sets = grid.get("index_sets") or {}
    sizes_expr = {}
    for ax in axes:
        decl = index_sets.get(ax)
        if not isinstance(decl, dict) or decl.get("kind") != "interval":
            findings.add(path, "L006", f"axes:{ax} is not an interval index set of grid {grid_names[0]}")
            return
        sizes_expr[ax] = decl.get("size")
    metaparams = list((grid.get("metaparameters") or {}))

    for tmpl_name, makearray in find_makearray_bodies(doc):
        regions = makearray.get("regions", [])
        if any(len(region) != len(axes) for region in regions):
            findings.add(
                path, "L006",
                f"{tmpl_name}: region rank != len(axes:{','.join(axes)})",
            )
            continue
        for sample in SAMPLE_SIZES:
            bindings = {m: sample for m in metaparams}
            try:
                sizes = [fold(sizes_expr[ax], bindings) for ax in axes]
                boxes = [
                    [(fold(lo, bindings), fold(hi, bindings)) for lo, hi in region]
                    for region in regions
                ]
            except ValueError as e:
                findings.add(path, "L006", f"{tmpl_name}: fold failed at size {sample}: {e}")
                break
            cells = 1
            for s in sizes:
                cells *= s
            if cells > MAX_COVERAGE_CELLS:
                continue
            covered = bytearray(cells)
            ok = True
            for box in boxes:
                idxs = [range(lo, hi + 1) for lo, hi in box]
                if any((lo < 1 or hi > size or lo > hi) for (lo, hi), size in zip(box, sizes)):
                    findings.add(
                        path, "L006",
                        f"{tmpl_name}: region {box} exceeds axis extents {sizes} at size {sample}",
                    )
                    ok = False
                    break
                # Mark coverage.
                def mark(dim: int, flat: int) -> None:
                    if dim == len(idxs):
                        covered[flat] = 1
                        return
                    stride = 1
                    for s in sizes[dim + 1:]:
                        stride *= s
                    for i in idxs[dim]:
                        mark(dim + 1, flat + (i - 1) * stride)
                mark(0, 0)
            if ok and not all(covered):
                gaps = covered.count(0)
                findings.add(
                    path, "L006",
                    f"{tmpl_name}: regions leave {gaps} uncovered cell(s) at size {sample} "
                    f"(axes {axes}, sizes {sizes})",
                )
                break


def validate_schema(files, findings: Findings) -> None:
    try:
        import jsonschema
    except ImportError:
        sys.stderr.write(
            "error: jsonschema not importable; pip install jsonschema "
            "(or run with the EarthSciSerialization venv python)\n"
        )
        sys.exit(2)
    schema = load_json(ess_root() / "esm-schema.json")
    validator = jsonschema.Draft202012Validator(schema)
    for path in files:
        try:
            doc = load_json(path)
        except json.JSONDecodeError:
            continue  # PARSE finding handled by lint
        for err in validator.iter_errors(doc):
            findings.add(path, "SCHEMA", f"{'/'.join(map(str, err.path))}: {err.message[:200]}")


def run_validate() -> int:
    lib = Library()
    findings = Findings()
    validate_schema(lib.files, findings)
    for path in lib.entries:
        lint_file(path, lib, findings)
    print(f"validated {len(lib.files)} file(s) ({len(lib.entries)} library entries)")
    return findings.report()


def run_invalid() -> int:
    base = REPO / "tests" / "invalid" / "lint"
    expected = load_json(base / "expected_errors.json")
    lib = Library()
    failures = 0
    for rel, spec in sorted(expected.items()):
        path = base / rel
        findings = Findings()
        if not path.exists():
            print(f"MISSING fixture {rel}")
            failures += 1
            continue
        lib.docs[path] = load_json(path)
        lint_file(path, lib, findings)
        got = {code for _, code, _ in findings.items}
        want = set(spec["codes"])
        if not want <= got:
            print(f"FAIL {rel}: expected codes {sorted(want)}, got {sorted(got)}")
            failures += 1
        else:
            print(f"ok   {rel}: {sorted(want)}")
    print(f"{len(expected)} fixture(s), {failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(run_invalid() if "--invalid" in sys.argv else run_validate())
