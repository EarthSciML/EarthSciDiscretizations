"""ESM documents in this repository: which files are documents, and how to read them.

A library file (or a conformance fixture) may be written in any of the three
document formats EarthSciAST reads, and the format is chosen by the file's
extension exactly as EarthSciAST's loaders choose it for a root file or a
`ref` (esm-spec §1.2.1):

    .md   Markdown with LaTeX math
    .typ  Typst
    else  JSON (this repository's JSON documents are named `.esm`)

The format is never guessed from the text. The Markdown and Typst readers are
EarthSciAST's own public `parse_markdown` / `parse_typst`, which return the
same document (the JSON data model) that the JSON file would hold, so every
consumer after `read_document` is format-agnostic. JSON is read with the
standard library, which is what EarthSciAST's loader does for it, so a corpus
of `.esm` files needs no EarthSciAST import at all.

A name such as `grid` or a rule's stem may be written in one format only: two
documents `foo.esm` and `foo.typ` in one directory are ambiguous (see
`duplicate_stems`, linted as L009 by scripts/validate-library.py).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

#: Every extension a document may carry, JSON first.
DOCUMENT_SUFFIXES = (".esm", ".typ", ".md")


def is_document(path: Path) -> bool:
    return path.suffix in DOCUMENT_SUFFIXES and path.is_file()


def documents(directory: Path, *, recursive: bool = False) -> list[Path]:
    """The documents in `directory` (and below it when `recursive`), sorted."""
    if not directory.is_dir():
        return []
    found = directory.rglob("*") if recursive else directory.iterdir()
    return sorted(p for p in found if is_document(p))


def find_document(directory: Path, stem: str) -> Path | None:
    """The document named `stem` in `directory`, in whichever format it is
    written (the first by DOCUMENT_SUFFIXES order if, against L009, there are
    several), or None."""
    for suffix in DOCUMENT_SUFFIXES:
        p = directory / (stem + suffix)
        if p.is_file():
            return p
    return None


def duplicate_stems(paths: list[Path]) -> dict[Path, list[Path]]:
    """{directory/stem: [paths]} for every name written in more than one format."""
    by_stem: dict[Path, list[Path]] = {}
    for p in paths:
        by_stem.setdefault(p.with_suffix(""), []).append(p)
    return {k: sorted(v) for k, v in by_stem.items() if len(v) > 1}


def decode_document(text: str, path: Path | str) -> Any:
    """Decode a document's `text` in the format `path`'s extension names.

    Raises ValueError for text that is not a document in that format:
    `json.JSONDecodeError` for JSON, and EarthSciAST's `ParseError` (a
    ValueError, code `parse_error` or `removed_op`) for Markdown and Typst.
    """
    suffix = Path(path).suffix
    if suffix == ".md":
        from earthsci_ast import parse_markdown

        return parse_markdown(text)
    if suffix == ".typ":
        from earthsci_ast import parse_typst

        return parse_typst(text)
    return json.loads(text)


def read_document(path: Path) -> Any:
    """Read the document at `path` (see decode_document)."""
    return decode_document(Path(path).read_text(encoding="utf-8"), path)
