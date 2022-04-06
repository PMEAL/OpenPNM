r"""

scikit-graph is a libary for working with networks or graphs. It is built
upon ``numpy`` and is compatible with ``scipy.sparse.csgraph``.  It offers
a protocol for working with more complex situations, specifically (1) when
nodes and edge have attributes which are used in computations that require
good performance, and (2) when the graph has a spatial component (i.e.
the nodes have coordinates).

The skgraph protocol is as follows:

    * All the information about a network or graph is stored in a standard
    python ``dict`` (e.g. ``net['node.coords']``)
    * The dict must contain ``edge.conns`` and ``node.coords``. Everything
    else is optional.
        * ``node.coords`` is an M-by-3 array containg the [x, y, z] locations
        of each node, where M is the number of nodes in the network.
        * ``edge.conns`` is an N-by-2 array of which nodes are connected to
        which, where N is the number of edges in the network.
    * node and edge indices are implied by their position the node or edge
    arrays. For e.g. g['node.coords'][0, :] contains the x,y,z coordinates
    of node 0.
    * The terms node and edge are used but all functions aime to be
    terminology agnostic, so that 'vertex.coords' can be used if really
    necessary.


"""


from . import generators
from . import io
from . import metrics
from . import operations
from . import queries
from . import simulations
from . import tools
from . import visualization
