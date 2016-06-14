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
            The cubic domain in which the base points lie.  Note that many of
            the Voronoi vertices will lie outside this domain, but these are
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

                **'intersected'** : (Not implemented yet!)
                All Voronoi edges that lie outside the
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

        # Reflect base points about domain to make flat faces up tessellation
        domain_size = sp.array(domain_size, ndmin=1)

        if face_type == 'rough':
            base_points = sp.rand(num_cells*(1.5**3), 3)
            base_points = base_points*(1.5*domain_size) - 0.25*domain_size
            inds = sp.all(base_points > [0, 0, 0], axis=1)
            inds = inds*sp.all(base_points < domain_size, axis=1)
            self.base_points = base_points[inds, :]
        elif face_type == 'reflected':
            base_points = sp.rand(num_cells, 3)
            base_points = base_points*domain_size
            self.base_points = base_points
            orig_points = base_points
            Nx, Ny, Nz = domain_size
            base_points = sp.vstack((base_points, [-1, 1, 1]*orig_points +
                                                  [2*Nx, 0, 0]))
            base_points = sp.vstack((base_points, [1, -1, 1]*orig_points +
                                                  [0, 2*Ny, 0]))
            base_points = sp.vstack((base_points, [1, 1, -1]*orig_points +
                                                  [0, 0, 2*Nz]))
            base_points = sp.vstack((base_points, [-1, 1, 1]*orig_points))
            base_points = sp.vstack((base_points, [1, -1, 1]*orig_points))
            base_points = sp.vstack((base_points, [1, 1, -1]*orig_points))
        else:
            raise Exception('Unrecognized face_type: ' + face_type)

        vor = sptl.Voronoi(points=base_points)
        internal_vertices = sp.zeros(vor.vertices.shape[0], dtype=bool)
        N = vor.vertices.shape[0]
        am = sp.sparse.lil_matrix((N, N), dtype=int)
        for item in vor.ridge_dict.keys():
            if sp.all(vor.vertices[vor.ridge_dict[item]] >= -0.1) and \
               sp.all(vor.vertices[vor.ridge_dict[item]] <= domain_size):
                internal_vertices[vor.ridge_dict[item]] = True
                vor.ridge_dict[item].append(vor.ridge_dict[item][0])
                hull = [vor.ridge_dict[item][0:-1], vor.ridge_dict[item][1:]]
                hull = sp.sort(sp.vstack(hull).T, axis=1)
                am[hull[:, 0], hull[:, 1]] = 1
        am = am[internal_vertices, :][:, internal_vertices]
        Nt = am.nnz
        Np = sp.sum(internal_vertices)
        am = am.tocoo()
        self.update({'pore.all': sp.ones((Np, ), dtype=bool)})
        self['pore.coords'] = vor.vertices[internal_vertices]
        self.update({'throat.all': sp.ones((Nt, ), dtype=bool)})
        self['throat.conns'] = sp.vstack([am.row, am.col]).T
