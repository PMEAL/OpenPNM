import scipy as sp
import scipy.spatial as sptl
from openpnm.network import Delaunay
from openpnm.topotools import trim
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Gabriel(Delaunay):
    r"""
    Creates a Gabriel network from a given set of base points.  This operates
    by performing a Deluanay tessellation, then removing connections that do
    not adhere to the definition of the `Gabriel graph
    <https://en.wikipedia.org/wiki/Gabriel_graph>`_

    This produces a network that has fewer throats than a Delaunay network.
    Since the longer-range throats tend to be removed this might be more
    realistic in some cases.

    Parameters
    ----------
    points : array_like
        An array of coordinates indicating the [x, y, z] locations of each
        point to use in the tessellation.  Note that the points must be given
        in rectilinear coordinates regardless of which domain ``shape`` was
        specified.  To convert between coordinate systems see the
        ``convert_coords`` function in the ``openpnm.topotools`` module.

    num_points : scalar
        The number of points to place in the domain, which will become the
        pore centers after the tessellation is performed.  This value is
        ignored if ``points`` are given.

    shape : array_like
        The size of the domain.  It's possible to create cubic, cylindrical,
        or spherical domains, as well as 2D square and circular by changing
        the domain ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, and
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y

        [r, z] - will produce a cylindrical domain with a radius of r and
        height of z

        [r, 0] - will produce a 2D circular domain with a radius of r

        [r] - will produce a spherical domain with a radius of r

    Examples
    --------
    >>> imoprt openpnm as op
    >>> import scipy as sp
    >>> pts = sp.rand(100, 3) * [1, 1, 0]  # Set z-axis to 0
    >>> gn = op.network.Gabriel(shape=[1, 1, 0], points=pts)
    >>> dn = op.network.Delaunay(shape=[1, 1, 0], points=pts)

    Now compare them side by side:

    >>> gn['pore.coords'] += [1, 0, 0]
    >>> topotools.merge_networks(dn, gn)
    >>> fig = topotools.plot_connections(dn)
    >>> fig = topotools.plot_coordinates(dn, c='r', s=100, fig=fig)

    .. image:: /../docs/static/images/gabriel_network.png
        :align: center

    """
    def __init__(self, **kwargs):
        # Generate Delaunay tessellation from super class, then trim
        super().__init__(**kwargs)
        points = self['pore.coords']
        conns = self['throat.conns']
        # Find centroid of each pair of nodes
        c = points[conns]
        m = (c[:, 0, :] + c[:, 1, :])/2
        # Find radius of circle connecting each pair of nodes
        r = sp.sqrt(sp.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
        # Use KD-Tree to find distance to nearest neighbors
        tree = sptl.cKDTree(points)
        n = tree.query(x=m, k=1)[0]
        # Identify throats whose centroid is not near an unconnected node
        g = sp.around(n, decimals=5) == sp.around(r, decimals=5)
        trim(self, throats=~g)
