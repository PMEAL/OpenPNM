"""
===============================================================================
DelaunayVoronoiDual: Generate a random network with complementary Delaunay and
Voronoi networks, including connectings between them
===============================================================================

"""
import OpenPNM as op
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
        self['pore.coords'] = coords
        self['throat.conns'] = conns
        self['pore.coords'] = sp.around(self['pore.coords'], decimals=10)

        # Label all pores and throats by type
        self['pore.delaunay'] = False
        self['pore.delaunay'][0:Npts] = True
        self['pore.voronoi'] = False
        self['pore.voronoi'][Npts:] = True
        # Label throats between Delaunay pores
        self['throat.delaunay'] = False
        Ts = sp.all(self['throat.conns'] < Npts, axis=1)
        self['throat.delaunay'][Ts] = True
        # Label throats between Voronoi pores
        self['throat.voronoi'] = False
        Ts = sp.all(self['throat.conns'] >= Npts, axis=1)
        self['throat.voronoi'][Ts] = True
        # Label throats connecting a Delaunay and a Voronoi pore
        self['throat.interconnect'] = False
        Ts = self.throats(labels=['delaunay', 'voronoi'], mode='not')
        self['throat.interconnect'][Ts] = True

    def trim_domain(self, domain_size=None):
        r"""
        """
        # Note num neighbors for later
        self['pore._num_neighbors'] = self.num_neighbors(pores=self.Ps)
        if len(domain_size) == 1:  # Spherical
            # Trim external Delaunay pores
            r = sp.sqrt(sp.sum(self['pore.coords']**2, axis=1))
            Ps = (r > domain_size)*self['pore.delaunay']
            self.trim(pores=Ps)
            # Trim external Voronoi pores
            Ps = self.find_neighbor_pores(pores=self['pore.delaunay'])
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self.trim(pores=Ps)
        elif len(domain_size) == 2:  # Cylindrical
            # Trim external Delaunay pores outside radius
            r = sp.sqrt(sp.sum(self['pore.coords'][:, [0, 1]]**2, axis=1))
            Ps = (r > domain_size[0])*self['pore.delaunay']
            self.trim(pores=Ps)
            # Trim external Delaunay pores above and below cylinder
            Ps1 = self['pore.coords'][:, 2] > domain_size[1]
            Ps2 = self['pore.coords'][:, 2] < 0
            Ps = self['pore.delaunay']*(Ps1 + Ps2)
            self.trim(pores=Ps)
            # Trim external Voronoi pores
            Ps = self.find_neighbor_pores(pores=self['pore.delaunay'])
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self.trim(pores=Ps)
        elif len(domain_size) == 3:  # Rectilinear
            # Trim external Delaunay pores
            Ps1 = sp.any(self['pore.coords'] > domain_size, axis=1)
            Ps2 = sp.any(self['pore.coords'] < [0, 0, 0], axis=1)
            Ps = self['pore.delaunay']*(Ps1 + Ps2)
            self.trim(pores=Ps)
            # Trim external Voronoi pores
            Ps = self.find_neighbor_pores(pores=self['pore.delaunay'])
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self.trim(pores=Ps)

        # Label Voronoi surface pores and throats
        Nn = self['pore._num_neighbors'] != self.num_neighbors(pores=self.Ps)
        self['pore.surface'] = Nn*self['pore.voronoi']
        Ps = self.pores('surface')
        Ts = self.find_neighbor_throats(pores=Ps, mode='intersection')
        self['throat.surface'] = self.tomask(throats=Ts)

        # Find boundary Delaunay pores
        Ps = self.pores('surface')
        Ps = self.find_neighbor_pores(pores=Ps)
        self['pore.boundary'] = self.tomask(pores=Ps)

        # Label all remaining pores and throats as internal
        self['pore.internal'] = True
        self['throat.internal'] = True

        # Cleanup
        del self['pore._num_neighbors']
