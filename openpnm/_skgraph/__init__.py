r"""

scikit-graph is a libary for working with networks or graphs. It is built
upon ``numpy`` and is compatible with ``scipy.sparse.csgraph``.  It offers
a protocol for working with more complex situations, specifically (1) when
sites and bonds have attributes which are used in computations that require
good performance, and (2) when the graph has a spatial component (i.e.
the pores have coordinates).

The skgraph protocol is as follows:

    * All the information about a network or graph is stored in a standard
    python ``dict`` (e.g. ``net['coords']``)
    * The dict must contain ``conns`` and ``coords``. Everything else is
    optional.
        * ``coords`` is an M-by-3 array containg the [x, y, z] locations of
        each site, where M is the number of sites in the network.
        * ``conns`` is an N-by-2 array of which sites are connected to
        which, where N is the number of bonds in the network.


"""


from . import generators
from . import io
from . import metrics
from . import operations
from . import queries
from . import simulations
from . import tools
from . import visualization
