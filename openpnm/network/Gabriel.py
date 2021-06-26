import numpy as np
import scipy.spatial as sptl
from openpnm.network import Delaunay
from openpnm.topotools import trim
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Gabriel(Delaunay):
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
        The size of the domain.  It's possible to create cubic, or 2D square
        domains by changing the domain ``shape`` as follows:

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

    Examples
    --------
    >>> import openpnm as op
    >>> import scipy as sp
    >>> import matplotlib.pyplot as plt
    >>> pts = np.random.rand(100, 3) * [1, 1, 0]  # Set z-axis to 0
    >>> gn = op.network.Gabriel(shape=[1, 1, 0], points=pts)
    >>> dn = op.network.Delaunay(shape=[1, 1, 0], points=pts)

    Now compare them side by side:

    >>> gn['pore.coords'] += [1, 0, 0]
    >>> op.topotools.merge_networks(dn, gn)
    >>> fig, ax = plt.subplots()
    >>> _ = op.topotools.plot_connections(dn, ax=ax)
    >>> _ = op.topotools.plot_coordinates(dn, c='r', s=100, ax=ax)

    .. image:: /../docs/_static/images/gabriel_network.png
        :align: center

    """

    def __init__(self, shape=[1, 1, 1], num_points=None, points=None, **kwargs):
        # Generate Delaunay tessellation from super class, then trim
        super().__init__(shape=shape, num_points=num_points, points=points, **kwargs)
        if 'pore.coords' in self.keys():
            points = self['pore.coords']
            conns = self['throat.conns']
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
            trim(self, throats=~g)
