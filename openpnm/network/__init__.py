r"""

**openpnm.network**

----

This module contains the ``GenericNetwork`` class, whose main purpose is to
manage the topological representation of the Network.  It also houses a
collection of Network generators.

----

**Available Network Generators**

OpenPNM includes a variety of Network generators.  The basically include two
families of topology: periodic lattices and tessellations of random points.

+---------------------+-------------------------------------------------------+
| Generator Name      | Description                                           |
+=====================+=======================================================+
| Cubic               | Simple cubic lattice with connectivity from 6 to 26   |
+---------------------+-------------------------------------------------------+
| CubicDual           | Body centered cubic lattice plus face centered nodes  |
|                     | on the surfaces                                       |
+---------------------+-------------------------------------------------------+
| CubicTemplate       | Simple cubic lattice with arbitrary domain shape      |
|                     | specified by a template image                         |
+---------------------+-------------------------------------------------------+
| Bravais             | Crystal lattice types including fcc, bcc, sc, and hcp |
+---------------------+-------------------------------------------------------+
| Delaunay            | Random network formed by Delaunay tessellation of     |
|                     | arbitrary base points                                 |
+---------------------+-------------------------------------------------------+
| Voronoi             | Random network formed by Voronoi tessellation of      |
|                     | arbitrary base points                                 |
+---------------------+-------------------------------------------------------+
| Gabriel             | Random network formed by Gabriel tessellation of      |
|                     | arbitrary base points                                 |
+---------------------+-------------------------------------------------------+
| DelaunayVoronoiDual | Combined and interconnected Voronoi and Delaunay      |
|                     | tessellations                                         |
+---------------------+-------------------------------------------------------+

----

**The GenericNetwork Class**

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
<https://en.wikipedia.org/wiki/Adjacency_matrix>`_.  Moreover, this is stored
using a `sparse matrix format <https://en.wikipedia.org/wiki/Sparse_matrix>`_
known as COO.  All netowrk objects store the COO matrix as ``'throat.conns'``.

The spatial location of each pore is stored in Cartesian coordinates [x, y, z],
under ``'pore.coords'``.  All networks must be 3D, so even a 2D network must
have a z-component (but set to 0).

The following methods are implemented on ``GenericNetwork``, and look into
the ``'throat.conns'`` and ``'pore.coords'`` as needed.

+-------------------------+---------------------------------------------------+
| Method                  | Description                                       |
+=========================+===================================================+
| num_neighbors           | Counts the number of neighbors with a given label |
+-------------------------+---------------------------------------------------+
| find_neighbor_pores     | Gets indices of pores neighboring a given pore    |
+-------------------------+---------------------------------------------------+
| find_neighbor_throats   | Gets indices of neighbor throats to a given pore  |
+-------------------------+---------------------------------------------------+
| find_connected_pores    | Gets indices of pores connected by a given throat |
+-------------------------+---------------------------------------------------+
| find_connecting_throat  | Gets indices of the throat joining pairs of pores |
+-------------------------+---------------------------------------------------+
| find_nearby_pores       | Find all pores within given distance of given pore|
+-------------------------+---------------------------------------------------+
| create_adjacency_matrix | Generates a weighted adjacency matrix             |
+-------------------------+---------------------------------------------------+
| create_incidence_matrix | Creates a weighted incidence matrix               |
+-------------------------+---------------------------------------------------+
| get_adjacency_matrix    | Returns an adjacency matrix with default weights  |
+-------------------------+---------------------------------------------------+
| get_incidence_matrix    | Returns an incidence matrix with default weights  |
+-------------------------+---------------------------------------------------+
| check_network_health    | Check various aspects of topology for problems    |
+-------------------------+---------------------------------------------------+

"""

from .GenericNetwork import GenericNetwork
from .Cubic import Cubic
from .CubicDual import CubicDual
from .Bravais import Bravais
from .CubicTemplate import CubicTemplate
from .DelaunayVoronoiDual import DelaunayVoronoiDual
from .Voronoi import Voronoi
from .Delaunay import Delaunay
from .Gabriel import Gabriel
