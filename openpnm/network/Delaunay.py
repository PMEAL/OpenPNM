import numpy as np
from openpnm import topotools
from openpnm.utils import logging
from openpnm.network import GenericNetwork
from openpnm.network.generators import delaunay
logger = logging.getLogger(__name__)


class Delaunay(GenericNetwork):
    r"""
    Random network formed by Delaunay tessellation of arbitrary base points

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

        [x, y, z] - will produce a normal cubic domain of dimension x, and
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
    This class always performs the tessellation on the full set of points, then
    trims any points that lie outside the given domain ``shape``.

    Examples
    --------
    >>> import numpy as np
    >>> import scipy as sp
    >>> import openpnm as op
    >>> import matplotlib as mpl
    >>> mpl.use('Agg')

    Supplying custom specified points:

    >>> pts = np.random.rand(200, 3)
    >>> gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])
    >>> gn.Np
    200

    Which can be quickly visualized using:

    >>> fig = op.topotools.plot_connections(network=gn)

    .. image:: /../docs/_static/images/delaunay_network_given_points.png
        :align: center

    Upon visualization it can be seen that this network is not very cubic.
    There are a few ways to combat this, but none will make a truly square
    domain.  Points can be generated that lie outside the domain ``shape``
    and they will be automatically trimmed.

    >>> pts = np.random.rand(300, 3)*1.2 - 0.1  # Must have more points for same density
    >>> gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])
    >>> gn.Np < 300  # Confirm base points have been trimmed
    True

    And visualizing:

    >>> fig = op.topotools.plot_connections(network=gn)

    .. image:: /../docs/_static/images/delaunay_network_w_trimmed_points.png
        :align: center

    If a domain with random base points but flat faces is needed use
    ``Voronoi``.

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Clean-up input points
        super().__init__(**kwargs)
        self.update(delaunay(points=points, shape=shape))
        Np = self['pore.coords'][:, 0].shape[0]
        Nt = self['throat.conns'][:, 0].shape[0]
        self['pore.all'] = np.ones((Np, ), dtype=bool)
        self['throat.all'] = np.ones((Nt, ), dtype=bool)
