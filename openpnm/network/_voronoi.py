import numpy as np
from openpnm.network import GenericNetwork
from openpnm.utils import Docorator
from openpnm._skgraph.generators import voronoi
from openpnm._skgraph.generators.tools import parse_points


docstr = Docorator()
__all__ = ['Voronoi']


@docstr.dedent
class Voronoi(GenericNetwork):
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

    %(GenericNetwork.parameters)s

    Notes
    -----
    By definition these points will each lie in the center of a Voronoi cell,
    so they will not be the pore centers.  The number of pores in the
    returned network thus will differ from the number of points supplied

    """

    def __init__(self, shape, points, trim=True, **kwargs):
        # Clean-up input points
        points = parse_points(shape=shape, points=points)
        super().__init__(**kwargs)
        net, vor = voronoi(points=points, shape=shape, trim=trim,
                           node_prefix='pore', edge_prefix='throat')
        self['pore.all'] = np.ones(net['pore.coords'].shape[0], dtype=bool)
        self['throat.all'] = np.ones(net['throat.conns'].shape[0], dtype=bool)
        self.update(net)
