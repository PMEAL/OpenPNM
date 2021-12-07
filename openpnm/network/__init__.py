r"""
Object model for storing topological information of the network
===============================================================

This module contains the ``GenericNetwork`` class, whose main purpose is to
manage the topological representation of the Network. It also houses a
collection of Network generators.

Available Network Generators
----------------------------

OpenPNM includes a variety of Network generators. The basically include two
families of topology: periodic lattices and tessellations of random points.

The GenericNetwork Class
------------------------

All of the above Network classes derive from the GenericNetwork class.  It is
a subclass of ``Base`` so contains methods for retrieving sets of pores based
on labels and so forth, but also contains the following additional methods
that are used soley for topological queries.

Pore networks require two essential pieces of information:

- the spatial location of pores
- the connectivity of which throats connect which pores

The ``GenericNetwork`` class and it's subclasses are responsible for storing,
managing, and utilizing this information.

Network topology is stored using `adjacency matrices
<https://en.wikipedia.org/wiki/Adjacency_matrix>`_. Moreover, this is stored
using a `sparse matrix format <https://en.wikipedia.org/wiki/Sparse_matrix>`_
known as COO. All netowrk objects store the COO matrix as ``'throat.conns'``.

The spatial location of each pore is stored in Cartesian coordinates [x, y, z],
under ``'pore.coords'``. All networks must be 3D, so even a 2D network must
have a z-component (but set to 0).

"""

from ._generic import GenericNetwork
from ._cubic import Cubic
from ._cubic_dual import CubicDual
from ._bravais import Bravais
from ._cubic_template import CubicTemplate
from ._delaunay_voronoi_dual import DelaunayVoronoiDual
from ._voronoi import Voronoi
from ._delaunay import Delaunay
from ._gabriel import Gabriel
