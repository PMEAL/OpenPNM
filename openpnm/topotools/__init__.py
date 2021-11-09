r"""
====================================
TopoTools (:mod:`openpnm.topotools`)
====================================

This module contains a selection of functions that deal specifically with
network topology.

.. currentmodule:: openpnm.topotools

Topological manipulation
------------------------

.. autosummary::
   :template: mybase.rst
   :toctree: generated/
   :nosignatures:

   add_boundary_pores
   clone_pores
   connect_pores
   dimensionality
   extend
   filter_pores_by_z
   find_surface_pores
   find_pore_to_pore_distance
   from_cylindrical
   from_spherical
   generate_base_points
   get_shape
   get_spacing
   iscoplanar
   isoutside
   is_fully_connected
   label_faces
   merge_networks
   merge_pores
   reduce_coordination
   reflect_base_points
   rotate_coords
   shear_coords
   stitch
   stitch_pores
   subdivide
   template_cylinder_annulus
   template_sphere_shell
   to_cylindrical
   to_spherical
   trim
   trim_occluded_throats

Percolation tools
-----------------

.. autosummary::
   :template: mybase.rst
   :toctree: generated/
   :nosignatures:

   ispercolating
   remove_isolated_clusters
   site_percolation
   bond_percolation
   find_clusters
   find_path

Graph tools
-----------

.. autosummary::
   :template: mybase.rst
   :toctree: generated/
   :nosignatures:

   drop_sites
   find_neighbor_sites
   find_neighbor_bonds
   find_connected_sites
   find_connecting_bonds
   find_complement
   istriu
   istril
   istriangular
   issymmetric
   tri_to_am
   vor_to_am
   conns_to_am

Plot tools
----------

.. autosummary::
   :template: mybase.rst
   :toctree: generated/
   :nosignatures:

   plot_tutorial
   plot_connections
   plot_coordinates
   plot_networkx
   plot_network_jupyter
   generate_voxel_image

"""

from .topotools import *
from .perctools import *
from .graphtools import *
from .plottools import *
from . import generators
