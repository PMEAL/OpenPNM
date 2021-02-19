r"""

**openpnm.topotools**

----

This module contains a selection of functions that deal specifically with
network topology.


"""

from .topotools import add_boundary_pores
from .topotools import clone_pores
from .topotools import connect_pores
from .topotools import dimensionality
from .topotools import extend
from .topotools import filter_pores_by_z
from .topotools import find_surface_pores
from .topotools import find_pore_to_pore_distance
from .topotools import generate_base_points
from .topotools import iscoplanar
from .topotools import isoutside
from .topotools import is_fully_connected
from .topotools import label_faces
from .topotools import merge_networks
from .topotools import merge_pores
from .topotools import reduce_coordination
from .topotools import reflect_base_points
from .topotools import rotate_coords
from .topotools import shear_coords
from .topotools import stitch
from .topotools import stitch_pores
from .topotools import subdivide
from .topotools import template_cylinder_annulus
from .topotools import template_sphere_shell
from .topotools import trim
from .topotools import trim_occluded_throats

from .perctools import ispercolating
from .perctools import remove_isolated_clusters
from .perctools import site_percolation
from .perctools import bond_percolation
from .perctools import find_clusters
from .perctools import find_path

from .graphtools import find_neighbor_sites
from .graphtools import find_neighbor_bonds
from .graphtools import find_connected_sites
from .graphtools import find_connecting_bonds
from .graphtools import find_complement
from .graphtools import istriu
from .graphtools import istril
from .graphtools import istriangular
from .graphtools import issymmetric
from .graphtools import _am_to_im
from .graphtools import _im_to_am
from .graphtools import tri_to_am
from .graphtools import vor_to_am
from .graphtools import conns_to_am

from .plottools import plot_tutorial
from .plottools import plot_connections
from .plottools import plot_coordinates
from .plottools import plot_networkx
from .plottools import plot_vpython
