import numpy as np
from openpnm.network import Network
from openpnm.utils import Docorator
from openpnm._skgraph.generators import voronoi_delaunay_dual
from openpnm._skgraph.generators.tools import parse_points
from openpnm._skgraph.operations import trim_nodes


docstr = Docorator()
__all__ = ['Voronoi']


@docstr.dedent
class Voronoi(Network):
    r"""
    Random network formed by Voronoi tessellation of arbitrary base points

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be
        used, or a scalar value indicating the number of points to
        generate
    shape : array_like
        The size of the domain.  It's possible to create cubic as well as 2D
        square domains by changing the ``shape`` as follows:

            [x, y, z]
                will produce a normal cubic domain of dimension x, and z
            [x, y, 0]
                will produce a 2D square domain of size x by y

    %(Network.parameters)s

    Notes
    -----
    By definition these points will each lie in the center of a Voronoi cell,
    so they will not be the pore centers.  The number of pores in the
    returned network thus will differ from the number of points supplied

    """

    def __init__(self, shape, points, trim=True, reflect=True, **kwargs):
        # Clean-up input points
        super().__init__(**kwargs)
        points = parse_points(shape=shape, points=points, reflect=reflect)
        net, vor, tri = voronoi_delaunay_dual(points=points,
                                              shape=shape,
                                              trim=trim,
                                              node_prefix='pore',
                                              edge_prefix='throat')
        net = trim_nodes(network=net, inds=np.where(net['pore.delaunay'])[0])
        self.update(net)
