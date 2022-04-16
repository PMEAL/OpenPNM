r"""
------------
scikit-graph
------------

scikit-graph (or ``skgraph``) is a libary for working with graphs (or
networks). It is built on the following principles:

* Graphs are represented by Python dictionaries filled with numpy ``ndarrays``
* All functionality is in the form of a library of functions that operate on
  the data in the dictionaries with no assumptions about class attributes or
  methods. The data just needs to be formatted and named correctly.
* Nodes and edges have data attributes, such as volume or length, which are
  used in computations
* Networks can be very large (billions of nodes) so good numerical performance
  is critical, to the extent that numpy and related packages can provide

The ``skgraph`` protocol is as follows:

* All the information about a network or graph is stored in a standard
  python ``dict``.
* The ``keys`` must be named following the convention of
  ``'<node_prefix>.<attribute>'`` or ``'<edge_prefix>.<attribute>'``. (The
  prefix 'node' and 'edge' are customizable)
* The dict *must contain* ``'edge.conns'`` and ``'node.coords'``. Everything
  else is optional.
  * ``node.coords`` is an M-by-3 array containg the [x, y, z] locations
  of each node, where M is the number of nodes in the network.
  * ``edge.conns`` is an N-by-2 array of which nodes are connected to
  which, where N is the number of edges in the network.
* node and edge indices are implied by their position in the node or edge
  arrays. For e.g. g['node.coords'][0, :] contains the x,y,z coordinates
  of node 0.

"""
import numpy as _np
from openpnm.utils import PrintableDict as _pdict
from openpnm.utils import SettingsAttr as _Settings

settings = _Settings()
settings.node_prefix = 'node'
settings.edge_prefix = 'edge'
settings.missing_values = {'bool': False,
                           'int': _np.nan,
                           'float': _np.nan,
                           'object': None}


from . import generators
from . import io
from . import metrics
from . import operations
from . import queries
from . import simulations
from . import tools
from . import visualization


def info(g):
    d = _pdict(g)
    d._key = 'Attribute'
    d._value = 'Description'
    print(d)


from .tools import get_edge_prefix
from .tools import get_node_prefix












