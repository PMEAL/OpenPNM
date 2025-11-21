import numpy as np
from openpnm.network import Network
from openpnm._skgraph.generators.tools import parse_points
from openpnm._skgraph.generators import voronoi_delaunay_dual
from openpnm.utils import Docorator


__all__ = ['DelaunayVoronoiDual']
docstr = Docorator()


@docstr.dedent
class DelaunayVoronoiDual(Network):
    r"""
    Combined and interconnected Voronoi and Delaunay tessellations

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        The size and shape of the domain:

        ========== ============================================================
        shape      result
        ========== ============================================================
        [x, y, z]  A 3D cubic domain of dimension x, y and z
        [x, y, 0]  A 2D square domain of size x by y
        ========== ============================================================

    trim : bool, optional
        If ``True`` then all vertices laying outside the domain will
        be removed. This is only useful if ``reflect=True``.
    reflect : bool, optional
        If ``True`` then the base points will be reflected across
        all the faces of the domain prior to performing the tessellation. This
        feature is best combined with ``trim=True``.
    f : float
        The fraction of points which should be reflected.  The default is 1 which
        reflects all the points in the domain, but this can lead to a lot of
        unnecessary points, so setting to 0.1 or 0.2 helps speed, but risks that
        the tessellation may not have smooth faces if not enough points are
        reflected.
    relaxation : int
        The number of time to iteratively relax the base points by moving them to
        the centroid of their respective Voronoi hulls. The default it 0.

    %(Network.parameters)s

    Notes
    -----
    A Delaunay tessellation is performed on a set of base points then the
    corresponding Voronoi diagram is generated.  Finally, each Delaunay node
    is connected to its neighboring Voronoi vertices to create interconnections
    between the two networks.

    All pores and throats are labelled according to their network (i.e.
    'pore.delaunay'), so they can be each assigned to a different Geometry.

    """

    def __init__(
        self,
        shape,
        points,
        trim=True,
        reflect=True,
        f=1,
        relaxation=0,
        **kwargs
    ):
        super().__init__(**kwargs)
        net, vor, tri = voronoi_delaunay_dual(
            shape=shape,
            points=points,
            trim=trim,
            reflect=reflect,
            relaxation=relaxation,
            f=f,
            node_prefix='pore',
            edge_prefix='throat',
        )
        self.update(net)
        self._post_init()
        self.vor = vor
        self.tri = tri

    def find_throat_facets(self, throats=None):
        r"""
        Finds the indicies of the Voronoi nodes that define the facet or
        ridge between the Delaunay nodes connected by the given throat.

        Parameters
        ----------
        throats : array_like
            The throats whose facets are sought.  The given throats should be
            from the 'delaunay' network. If no throats are specified, all
            'delaunay' throats are assumed.

        Notes
        -----
        The method is not well optimized as it scans through each given throat
        inside a for-loop, so it could be slow for large networks.

        """
        if throats is None:
            throats = self.throats('delaunay')
        temp = []
        tvals = self['throat.interconnect'].astype(int)
        am = self.create_adjacency_matrix(weights=tvals, fmt='lil',
                                          drop_zeros=True)
        for t in throats:
            P12 = self['throat.conns'][t]
            Ps = list(set(am.rows[P12][0]).intersection(am.rows[P12][1]))
            temp.append(Ps)
        return np.array(temp, dtype=object)

    def find_pore_hulls(self, pores=None):
        r"""
        Finds the indices of the Voronoi nodes that define the convex hull
        around the given Delaunay nodes.

        Parameters
        ----------
        pores : array_like
            The pores whose convex hull are sought.  The given pores should be
            from the 'delaunay' network.  If no pores are given, then the hull
            is found for all 'delaunay' pores.

        Notes
        -----
        This metod is not fully optimized as it scans through each pore in a
        for-loop, so could be slow for large networks.

        """
        if pores is None:
            pores = self.pores('delaunay')
        temp = []
        tvals = self['throat.interconnect'].astype(int)
        am = self.create_adjacency_matrix(weights=tvals, fmt='lil',
                                          drop_zeros=True)
        for p in pores:
            Ps = am.rows[p]
            temp.append(Ps)
        return np.array(temp, dtype=object)
