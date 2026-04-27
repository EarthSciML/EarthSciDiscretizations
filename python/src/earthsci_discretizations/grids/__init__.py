"""Grid accessor runtimes.

Per-family modules land here as Phase 3 implementation beads complete. Each
family exposes a generator function `earthsci_discretizations.grids.<family>(**opts)`
returning a `Grid` object per `docs/GRIDS_API.md` §2.4, §3.2.
"""

from .arakawa import ArakawaGrid, BaseGrid, CartesianBase, arakawa
from .cartesian import CartesianGrid, cartesian
from .cubed_sphere import CubedSphereGrid, cubed_sphere
from .duo import DuoGrid, DuoLoader, duo
from .lat_lon import LatLonGrid, lat_lon
from .mpas import MpasGrid, MpasLoader, MpasMeshData, check_mesh, mpas, mpas_mesh_data
from .vertical import VerticalGrid, vertical

__all__ = [
    "arakawa",
    "ArakawaGrid",
    "BaseGrid",
    "CartesianBase",
    "cartesian",
    "CartesianGrid",
    "cubed_sphere",
    "CubedSphereGrid",
    "duo",
    "DuoGrid",
    "DuoLoader",
    "lat_lon",
    "LatLonGrid",
    "mpas",
    "MpasGrid",
    "MpasLoader",
    "MpasMeshData",
    "mpas_mesh_data",
    "check_mesh",
    "vertical",
    "VerticalGrid",
]
