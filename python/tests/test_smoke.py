"""Smoke test: package imports and exposes the grids namespace."""

import earthsci_discretizations
from earthsci_discretizations import grids


def test_version():
    assert isinstance(earthsci_discretizations.__version__, str)
    assert earthsci_discretizations.__version__


def test_grids_namespace_importable():
    assert grids is earthsci_discretizations.grids
