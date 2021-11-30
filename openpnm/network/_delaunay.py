from openpnm import topotools
from openpnm.utils import logging
from openpnm.network import DelaunayVoronoiDual
from auto_all import start_all, end_all


logger = logging.getLogger(__name__)
start_all()


class Delaunay(DelaunayVoronoiDual):
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

    Examples
    --------
    .. plot::

        import numpy as np
        import openpnm as op
        import matplotlib.pyplot as plt

        # Supplying custom specified points
        pts = np.random.rand(200, 3)
        gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])

        # Check the number of pores in 'gn'
        print(gn.Np)

        # Which can be quickly visualized using
        fig, ax = plt.subplots(figsize=(5, 5))
        op.topotools.plot_connections(network=gn, ax=ax)

        plt.show()

    Upon visualization it can be seen that this network is not very cubic.
    There are a few ways to combat this, but none will make a truly square
    domain. Points can be generated that lie outside the domain ``shape``
    and they will be automatically trimmed.

    .. plot::

        import numpy as np
        import openpnm as op
        import matplotlib.pyplot as plt

        # Must have more points for same density
        pts = np.random.rand(300, 3)*1.2 - 0.1
        gn = op.network.Delaunay(points=pts, shape=[1, 1, 1])

        # Confirm base points have been trimmed
        print(gn.Np < 300)

        # And visualizing
        fig, ax = plt.subplots(figsize=(5, 5))
        op.topotools.plot_connections(network=gn, ax=ax)

        plt.show()

    If a domain with random base points but flat faces is needed use
    ``Voronoi``.

    """

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        # Clean-up input points
        points = self._parse_points(shape=shape, points=points)
        super().__init__(shape=shape, points=points, **kwargs)

        # Initialize network object
        topotools.trim(network=self, pores=self.pores(['voronoi']))
        pop = ['pore.voronoi', 'throat.voronoi', 'throat.interconnect',
               'pore.delaunay', 'throat.delaunay']
        for item in pop:
            del self[item]

        # Trim additional pores that are missed by the parent class's trimming
        Ps = topotools.isoutside(coords=self['pore.coords'], shape=shape)
        topotools.trim(network=self, pores=Ps)

end_all()
