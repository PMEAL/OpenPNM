import numpy as np
from openpnm.network import GenericNetwork
from openpnm._skgraph.generators import delaunay, tools
from openpnm._skgraph.tools import isoutside
from openpnm._skgraph.operations import trim_nodes


__all__ = ['Delaunay']


class Delaunay(GenericNetwork):
    r"""
    Random network formed by Delaunay tessellation of arbitrary base points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        The size of the domain.  It's possible to create cubic as well as 2D
        square domains by changing the ``shape`` as follows:

            [x, y, z]
                will produce a normal cubic domain of dimension x, and and z
            [x, y, 0]
                will produce a 2D square domain of size x by y

    name : str
        An optional name for the object to help identify it.  If not given,
        one will be generated.

    See Also
    --------
    Gabriel
    Voronoi
    DelaunayVoronoiDual

    Notes
    -----
    This class always performs the tessellation on the full set of points, then
    trims any points that lie outside the given domain ``shape``.

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        points = tools.parse_points(shape=shape, points=points)
        net, tri = delaunay(points=points, shape=shape,
                            node_prefix='pore', edge_prefix='throat')
        Ps = isoutside(coords=net['pore.coords'], shape=shape)
        net = trim_nodes(g=net, inds=Ps, node_prefix='pore',
                         edge_prefix='throat')
        self.update(net)
        self['pore.all'] = np.ones(self['pore.coords'].shape[0], dtype=bool)
        self['throat.all'] = np.ones(self['throat.conns'].shape[0], dtype=bool)
