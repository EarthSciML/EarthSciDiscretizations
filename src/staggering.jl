"""
C-grid staggering definitions.

The staggered-location vocabulary (`VarLocation`) lives here; it is the shared location
code used by the arakawa grid family and its declarative construction FAQ
(`src/arakawa_faq.jl`), as well as the DUO dual-geometry FAQ. The imperative
staggered-shape arithmetic (`grid_size` / `full_array_size` / `ghost_array_size`) that
used to live here was declarativized into the arakawa construction FAQ's
`_ark_location_shape` (S4) and deleted by S5 (esd-3we.5).
"""

@enum VarLocation CellCenter UEdge VEdge Corner
