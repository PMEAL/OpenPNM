import numpy as np
import scipy.spatial as sptl
from openpnm import topotools
from openpnm.network import Network
from openpnm._skgraph.generators.tools import parse_points
from openpnm._skgraph.generators import voronoi_delaunay_dual


__all__ = ['DelaunayVoronoiDual']


class DelaunayVoronoiDual(Network):
    r"""
    Combined and interconnected Voronoi and Delaunay tessellations

    Parameters
    ----------
    points : array_like or int
        Can either be an N-by-3 array of point coordinates which will be used,
        or a scalar value indicating the number of points to generate
    shape : array_like
        The size and shape of the domain used for generating and trimming
        excess points. The coordinates are treated as the outer corner of a
        rectangle [x, y, z] whose opposite corner lies at [0, 0, 0].
        By default, a domain size of [1, 1, 1] is used. To create a 2D network
        set the Z-dimension to 0.

    Notes
    -----
    A Delaunay tessellation is performed on a set of base points then the
    corresponding Voronoi diagram is generated.  Finally, each Delaunay node
    is connected to it's neighboring Voronoi vertices to create interaction
    between the two networks.

    All pores and throats are labelled according to their network (i.e.
    'pore.delaunay'), so they can be each assigned to a different Geometry.

    The dual-nature of this network is meant for modeling transport in the void
    and solid space simultaneously by treating one network (i.e. Delaunay) as
    voids and the other (i.e. Voronoi) as solid.  Interaction such as heat
    transfer between the solid and void can be accomplished via the
    interconnections between the Delaunay and Voronoi nodes.

    """

    def __init__(self, shape, points, trim=True, reflect=True, **kwargs):
        super().__init__(**kwargs)
        points = parse_points(shape=shape, points=points, reflect=reflect)

        net, vor, tri = voronoi_delaunay_dual(shape=shape,
                                              points=points,
                                              trim=trim,
                                              node_prefix='pore',
                                              edge_prefix='throat')
        self.update(net)
        self._vor = vor
        self._tri = tri

    @property
    def tri(self):
        """A shortcut to get a handle to the Delanuay subnetwork"""
        if not hasattr(self, '_tri'):
            points = self._vor.points
            self._tri = sptl.Delaunay(points=points)
        return self._tri

    @property
    def vor(self):
        """A shortcut to get a handle to the Voronoi subnetwork"""
        return self._vor

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
