r"""

**openpnm.topotools**

----

This module contains a selection of functions that deal specifically with
network topology.


"""

from .topotools import add_boundary_pores
from .topotools import bond_percolation
from .topotools import clone_pores
from .topotools import connect_pores
from .topotools import dimensionality
from .topotools import extend
from .topotools import find_path
from .topotools import find_surface_pores
from .topotools import find_neighbor_sites
from .topotools import find_neighbor_bonds
from .topotools import find_connected_sites
from .topotools import find_connecting_bonds
from .topotools import find_pore_to_pore_distance
from .topotools import find_clusters
from .topotools import find_complement
from .topotools import generate_base_points
from .topotools import iscoplanar
from .topotools import isoutside
from .topotools import issymmetric
from .topotools import ispercolating
from .topotools import istriu
from .topotools import istril
from .topotools import istriangular
from .topotools import label_faces
from .topotools import merge_networks
from .topotools import merge_pores
from .topotools import reduce_coordination
from .topotools import reflect_base_points
from .topotools import remove_isolated_clusters
from .topotools import rotate_coords
from .topotools import shear_coords
from .topotools import site_percolation
from .topotools import stitch
from .topotools import subdivide
from .topotools import template_cylinder_annulus
from .topotools import template_sphere_shell
from .topotools import trim
from .topotools import trim_occluded_throats
from .topotools import vor_to_am
from .topotools import tri_to_am
from .topotools import conns_to_am
from .plottools import plot_tutorial
from .plottools import plot_connections
from .plottools import plot_coordinates
from .plottools import plot_networkx
from .plottools import plot_vpython
