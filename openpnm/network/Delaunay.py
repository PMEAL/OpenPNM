from openpnm.network import DelaunayVoronoiDual
from openpnm import topotools
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Delaunay(DelaunayVoronoiDual):
    r"""
    Random network formed by Delaunay tessellation of arbitrary base points

    Parameters
    ----------
    num_points : scalar
        The number of points to place in the domain, which will become the
        pore centers after the tessellation is performed.  This value is
        ignored if ``points`` are given.

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
    >>> import openpnm as op
    >>> import scipy as sp
    >>> import matplotlib as mpl
    >>> mpl.use('Agg')

    Supplying custom specified points:

    >>> pts = sp.rand(200, 3)
    >>> gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])
    >>> gn.Np
    200

    Which can be quickly visualized using:

    >>> fig = op.topotools.plot_connections(network=gn)

    .. image:: /../docs/static/images/delaunay_network_given_points.png
        :align: center

    Upon visualization it can be seen that this network is not very cubic.
    There are a few ways to combat this, but none will make a truly square
    domain.  Points can be generated that lie outside the domain ``shape``
    and they will be automatically trimmed.

    >>> pts = sp.rand(300, 3)*1.2 - 0.1  # Must have more points for same density
    >>> gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])
    >>> gn.Np < 300  # Confirm base points have been trimmed
    True

    And visualizing:

    >>> fig = op.topotools.plot_connections(network=gn)

    .. image:: /../docs/static/images/delaunay_network_w_trimmed_points.png
        :align: center

    If a domain random base points, but truly flat faces is needed use
    ``Voronoi``.

    """

    def __init__(self, shape=None, num_points=None, **kwargs):
        # Clean-up input points
        points = kwargs.pop('points', None)
        points = self._parse_points(shape=shape,
                                    num_points=num_points,
                                    points=points)
        # Initialize network object
        super().__init__(shape=shape, points=points, **kwargs)
        topotools.trim(network=self, pores=self.pores(['voronoi']))
        pop = ['pore.voronoi', 'throat.voronoi', 'throat.interconnect',
               'pore.delaunay', 'throat.delaunay']
        for item in pop:
            del self[item]

        # Trim additional pores that are missed by the parent class's trimming
        Ps = topotools.isoutside(coords=self['pore.coords'], shape=shape)
        topotools.trim(network=self, pores=Ps)
