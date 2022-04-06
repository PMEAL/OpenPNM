import logging
import numpy as np
import scipy.spatial as sptl
from openpnm.network import GenericNetwork
from openpnm._skgraph.generators import gabriel
from openpnm._skgraph.operations import trim_nodes


logger = logging.getLogger(__name__)
__all__ = ['Gabriel']


class Gabriel(GenericNetwork):
    r"""
    Random network formed by Gabriel tessellation of arbitrary base points

    This operates by performing a Deluanay tessellation, then removing
    connections that do not adhere to the definition of the `Gabriel graph
    <https://en.wikipedia.org/wiki/Gabriel_graph>`_

    This produces a network that has fewer throats than a Delaunay network.
    Since the longer-range throats tend to be removed this might be more
    realistic in some cases.

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        The size of the domain.  It's possible to create cubic, or 2D square
        domains by changing the domain ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, and
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y


    """

    def __init__(self, shape=[1, 1, 1], points=None, trim=True, **kwargs):
        super().__init__(**kwargs)
        net = gabriel(shape=shape, points=points,
                      node_prefix='pore', edge_prefix='throat')
        points = net['pore.coords']
        conns = net['throat.conns']
        # Find centroid of each pair of nodes
        c = points[conns]
        m = (c[:, 0, :] + c[:, 1, :])/2
        # Find radius of circle connecting each pair of nodes
        r = np.sqrt(np.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
        # Use KD-Tree to find distance to nearest neighbors
        tree = sptl.cKDTree(points)
        n = tree.query(x=m, k=1)[0]
        # Identify throats whose centroid is not near an unconnected node
        g = np.around(n, decimals=5) == np.around(r, decimals=5)
        net = trim_nodes(net, inds=~g)
        self.update(net)
