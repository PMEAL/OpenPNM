"""
===============================================================================
Voronoi: Generate a random network based on Voronoi tessellation of points
===============================================================================

"""

import scipy as sp
import scipy.spatial as sptl
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Voronoi(GenericNetwork):
    r"""

    """

    def __init__(self, num_cells, domain_size=[1, 1, 1], face_type='none',
                 base_points=None, **kwargs):
        r"""
        Generate a random network of pores connected via a Voronoi
        tessellation.

        Parameters
        ----------
        num_cells : int
            This is the number of basepoint around with the Voronoi
            tessellation is constructed.  The actual number of pores in the
            network is controlled by the number Voronoi vertices which cannot
            be known a priori (I think think)

        domain_size : array_like (3 x 1)
            The cubic domain where the base points lie.  Note that many of the
            Voronoi vertices will lie outside this domain, but these are
            trimmed based on the ``face_type`` argument.

        face_type : string
            Indicates how the faces of the network should be trimmed.  Options
            are:

                **'reflected'** : The base points are reflected across each
                face of the network resulting in Voronoi vertices that lie
                on the plane(s) of reflection and thus a flat surface.

                **'rough'** : Any Voronoi vertices that lie outside the
                domain are dropped from the network, leaving a rough exposed
                surface.

                **'intersected'** : All Voronoi edges that lie outside the
                domain are cut-off at the intersection with the domain face,
                and new Voronoi vertices are created at the domain face.  These
                new vertices are NOT connected to each other so no vertices
                will lie on the domain faces, thus the final result is not a
                rigorous Voronoi diagram, but the topology is a bit more
                natural.

                **'none'** : All vertices are kept as is.  This is useful when
                special ''base_points'' are supplied, and the final trimming,
                if any is performed by the user.

        """
        super().__init__(**kwargs)
        domain_size = sp.array(domain_size, ndmin=1)
        base_points = sp.rand(num_cells, 3)
        base_points = base_points*domain_size
        orig_points = base_points
        base_points = sp.vstack((base_points, [-1, 1, 1]*orig_points +
                                              [2, 0, 0]))
        base_points = sp.vstack((base_points, [1, -1, 1]*orig_points +
                                              [0, 2, 0]))
        base_points = sp.vstack((base_points, [1, 1, -1]*orig_points +
                                              [0, 0, 2]))
        base_points = sp.vstack((base_points, [-1, 1, 1]*orig_points))
        base_points = sp.vstack((base_points, [1, -1, 1]*orig_points))
        base_points = sp.vstack((base_points, [1, 1, -1]*orig_points))
        self.vor = sptl.Voronoi(points=base_points)
        # Extract pore coordinates from Voronoi vertices
        a = []
        [a.extend(self.vor.regions[i]) for i in self.vor.point_region[0:300]]
        inds = sp.unique(a)
        self.update({'pore.all': sp.ones_like(inds, dtype=bool)})
        self['pore.coords'] = self.vor.vertices[inds]
        # Extract throat connections from ridges
        am = sp.sparse.lil_matrix((self.Np, self.Np), dtype=int)
        region_mask = sp.zeros(len(self.vor.point_region) + 1, dtype=bool)
        inds = self.vor.point_region[0:300]
        region_mask[inds] = True
        i = 0
        for ridge in self.vor.ridge_vertices:
            if sp.all(region_mask[self.vor.ridge_points[i]]):
                ridge.append(ridge[0])
                conns = sp.vstack((ridge[0:-1], ridge[1:])).T
                am[conns[:, 0], conns[:, 1]] = 1























