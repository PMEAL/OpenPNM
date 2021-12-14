import numpy as np
import scipy.spatial as sptl
from openpnm.network import Delaunay
from openpnm.topotools import trim
from openpnm.utils import logging


logger = logging.getLogger(__name__)
__all__ = ['Gabriel']


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
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        The size of the domain.  It's possible to create cubic, or 2D square
        domains by changing the domain ``shape`` as follows:

        [x, y, z] - will produce a normal cubic domain of dimension x, and
        and z

        [x, y, 0] - will produce a 2D square domain of size x by y

    name : str
        An optional name for the object to help identify it. If not given,
        one will be generated.

    Examples
    --------
    .. plot::

       import openpnm as op
       import matplotlib.pyplot as plt

       pts = np.random.rand(100, 3) * [1, 1, 0]  # Set z-axis to 0
       gn = op.network.Gabriel(shape=[1, 1, 0], points=pts)
       dn = op.network.Delaunay(shape=[1, 1, 0], points=pts)

       # Now compare them side by side
       gn['pore.coords'] += [1, 0, 0]
       op.topotools.merge_networks(dn, gn)

       # Visualization
       fig, ax = plt.subplots(figsize=(6, 3))
       op.topotools.plot_connections(dn, ax=ax)
       op.topotools.plot_coordinates(dn, c='r', s=100, ax=ax)

       plt.axis("off")
       plt.show()

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Generate Delaunay tessellation from super class, then trim
        super().__init__(shape=shape, points=points, **kwargs)
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
