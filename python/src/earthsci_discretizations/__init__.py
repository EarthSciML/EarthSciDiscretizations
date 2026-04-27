"""earthsci_discretizations — Python binding for EarthSciDiscretizations.

The `grids` submodule houses per-family grid accessor runtimes conforming to
the cross-binding contract in `docs/GRIDS_API.md`.
"""

__version__ = "0.1.0"

from . import grids

__all__ = ["grids", "__version__"]
