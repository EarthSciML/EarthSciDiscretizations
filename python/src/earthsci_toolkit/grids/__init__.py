"""Grid accessor runtimes.

Per-family modules land here as Phase 3 implementation beads complete. Each
family exposes a generator function `earthsci_toolkit.grids.<family>(**opts)`
returning a `Grid` object per `docs/GRIDS_API.md` §2.4, §3.2.
"""

from .cubed_sphere import CubedSphereGrid, cubed_sphere

__all__ = ["cubed_sphere", "CubedSphereGrid"]
