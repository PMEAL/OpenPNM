import numpy as np
from openpnm import topotools
from openpnm.utils import logging
from openpnm.network import GenericNetwork
from openpnm.network.generators import delaunay, tools
logger = logging.getLogger(__name__)


class Delaunay(GenericNetwork):
    r"""
    Random network formed by Delaunay tessellation of base points

    Parameters
    ----------
    points : array_like
        An array of coordinates indicating the [x, y, z] locations of each
        point to use in the tessellation.  Note that the points must be given
        in rectilinear coordinates regardless of which domain ``shape`` was
        specified.  To convert between coordinate systems see the
        ``convert_coords`` function in the ``openpnm.topotools`` module.
    shape : array_like
        The size of the domain.  It's possible to create cubic as well as 2D
        square domains by changing the ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, y
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.
    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

    See Also
    --------
    Gabriel
    Voronoi
    DelaunayVoronoiDual

    Notes
    -----
    This class always performs the tessellation on the full set of points,
    then trims any points that lie outside the given domain ``shape``.

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Clean-up input points
        super().__init__(**kwargs)
        points = tools.parse_points(points=points, shape=shape)
        net, tri = delaunay(points=points, shape=shape)
        net = tools.add_all_label(net)
        self.update(net)
        self._tri = tri
