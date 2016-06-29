"""
===============================================================================
DelaunayVoronoiDual: Generate a random network with complementary Delaunay and
Voronoi networks, including connectings between them
===============================================================================

"""

import scipy as sp
import scipy.spatial as sptl
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class DelaunayVoronoiDual(GenericNetwork):
    r"""


    """

    def __init__(self, points, domain_size=[1, 1, 1], **kwargs):
        super().__init__(**kwargs)

        # Perform tessellation
        vor = sptl.Voronoi(points=points)

        # Combine points
        pts_vor = vor.vertices
        pts_all = sp.vstack((points, pts_vor))
        Npts = sp.size(points, 0)
        Nvor = sp.size(pts_vor, 0)
        Nall = Nvor + Npts

        # Create adjacency matrix in lil format for quick matrix construction
        am = sp.sparse.lil_matrix((Nall, Nall))
        for ridge in vor.ridge_dict.keys():
            # Make Delaunay-to-Delauny connections
            [am.rows[i].extend([ridge[0], ridge[1]]) for i in ridge]
            row = vor.ridge_dict[ridge]
            if -1 not in row:
                # Index Voronoi vertex numbers by Npts
                row = [i + Npts for i in row]
                # Make Voronoi-to-Delaunay connections
                [am.rows[i].extend(row) for i in ridge]
                # Make Voronoi-to-Voronoi connections
                row.append(row[0])
                [am.rows[row[i]].append(row[i+1]) for i in range(len(row)-1)]
                # Ensure connections are made symmetrically
                [am.rows[row[i+1]].append(row[i]) for i in range(len(row)-1)]
        # Finalize adjacency matrix by assigning data values to each location
        am.data = am.rows  # Values don't matter, only shape, so use 'rows'
        # Retrieve upper triangle and convert to csr to remove duplicates
        am = sp.sparse.triu(A=am, k=1, format='csr')
        # Convert to COO format for OpenPNM compatibility
        am = am.tocoo()

        # Translate adjacency matrix and points to OpenPNM format
        coords = pts_all
        conns = sp.vstack((am.row, am.col)).T
        Np = sp.size(coords, axis=0)
        Nt = sp.size(conns, axis=0)
        self.update({'pore.all': sp.ones((Np, ), dtype=bool)})
        self.update({'throat.all': sp.ones((Nt, ), dtype=bool)})

        # Apply labels to pores
        self['pore.coords'] = coords
        self['throat.conns'] = conns

        # Initialize pore and throat labels
        self['pore.delaunay'] = False
        self['throat.delaunay'] = False
        self['pore.voronoi'] = False
        self['throat.voronoi'] = False
        self['throat.interconnect'] = False
        self['pore.internal'] = False
        self['throat.internal'] = False
        self['pore.external'] = False
        self['throat.external'] = False
        self['pore.surface'] = False
        self['throat.surface'] = False
        self['pore.boundary'] = False
        self['throat.boundary'] = False

        # Label all pores and throats by type
        self['pore.delaunay'][0:Npts] = True
        self['pore.voronoi'][Npts:] = True
        # Label throats between Delaunay pores
        Ts = sp.all(self['throat.conns'] < Npts, axis=1)
        self['throat.delaunay'][Ts] = True
        # Label throats between Voronoi pores
        Ts = sp.all(self['throat.conns'] >= Npts, axis=1)
        self['throat.voronoi'][Ts] = True
        # Label throats connecting a Delaunay and a Voronoi pore
        Ts = self.throats(labels=['delaunay', 'voronoi'], mode='not')
        self['throat.interconnect'][Ts] = True

        # Label Delaunay pores as internal or external
        Ps = self.tomask(pores=self.pores('delaunay'))
        self['pore.external'][Ps*sp.any(self['pore.coords'] > domain_size, axis=1)] = True
        self['pore.external'][Ps*sp.any(self['pore.coords'] < [0, 0, 0], axis=1)] = True
        self['pore.internal'][~self['pore.external']*self['pore.delaunay']] = True
        # Label Delaunay boundary pores
        Ps = self.pores(['delaunay', 'internal'], mode='intersection')
        Ps = self.find_neighbor_pores(pores=Ps, mode='union')
        Ps = self.filter_by_label(pores=Ps, labels='delaunay')
        self['pore.boundary'][Ps] = True
        self['pore.external'][Ps] = False
        # Label internal Voronoi pores
        Ps = self.pores(['delaunay', 'internal'], mode='intersection')
        Ps = self.find_neighbor_pores(pores=Ps, mode='union')
        Ps = self.filter_by_label(pores=Ps, labels='voronoi')
        self['pore.internal'][Ps] = True
        self.trim(pores=self.pores(['internal', 'boundary'], mode='not'))

        # Label Voronoi surface pores and throats
        Ps = self.pores('boundary')
        Ps = self.find_neighbor_pores(pores=Ps)
        Ps = self.filter_by_label(pores=Ps, labels=['voronoi'])
        self['pore.surface'][Ps] = True
        Ts = self.find_neighbor_throats(pores=Ps, mode='intersection')
        self['throat.surface'][Ts] = True

        # Label Boundardy throats
        Ps = self.pores('internal')
        Ts = self.find_neighbor_throats(pores=Ps, mode='not_intersection')
        self['throat.boundary'][Ts] = True
