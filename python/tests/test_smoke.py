"""Smoke test: package imports and exposes the grids namespace."""

import earthsci_toolkit
from earthsci_toolkit import grids


def test_version():
    assert isinstance(earthsci_toolkit.__version__, str)
    assert earthsci_toolkit.__version__


def test_grids_namespace_importable():
    assert grids is earthsci_toolkit.grids
