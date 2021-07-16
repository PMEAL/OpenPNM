import numpy as np
import scipy.spatial as sptl
from openpnm.network.generators import delaunay, gabriel
from openpnm.network import GenericNetwork
from openpnm.topotools import trim
from openpnm.utils import logging
logger = logging.getLogger(__name__)


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
    points : array_like, or scalar
        An array of coordinates indicating the [x, y, z] locations of each
        point to use in the tessellation.  Note that the points must be given
        in rectilinear coordinates regardless of which domain ``shape`` was
        specified.  To convert between coordinate systems see the
        ``convert_coords`` function in the ``openpnm.topotools`` module.

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

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Generate Delaunay tessellation from super class, then trim
        super().__init__(**kwargs)
        # Get delaunay tessellation
        dl = delaunay(shape=shape, points=points)
        # Pass delaunay to gabriel to trim offending throats
        gl = gabriel(dl)
        # Update self with network info
        self.update(gl)
        Np = self['pore.coords'][:, 0].shape[0]
        Nt = self['throat.conns'][:, 0].shape[0]
        self['pore.all'] = np.ones((Np, ), dtype=bool)
        self['throat.all'] = np.ones((Nt, ), dtype=bool)
