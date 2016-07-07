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

    def __init__(self, points, domain_size, **kwargs):
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
        self['throat.conns'] = conns
        self['pore.coords'] = sp.around(coords, decimals=10)

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

        self._trim_domain(domain_size=domain_size)

    def _trim_domain(self, domain_size=None):
        r"""
        Trims pores that lie outside the specified domain.

        Parameters
        ----------
        domain_size : array_like
            The size and shape of the domain beyond which points should be
            trimmed. The argument is treated as follows:

            **sphere** : If a scalar or single element list is received, it's
            treated as the radius [r] of a sphere centered on [0, 0, 0].

            **cylinder** : If a two-element list is received it's treated as
            the radius and height of a cylinder [r, z] whose central axis
            starts at [0, 0, 0] and extends in the positive z-direction.

            **rectangle** : If a three element list is received, it's treated
            as the outer corner of rectangle [x, y, z] whose opposite corner
            lies at [0, 0, 0].

        Notes
        -----
        This function assumes that some Delaunay nodes exist outside the
        given ``domain_size``.  These points can either be the result of
        reflecting the base points or simply creating points beyond the
        domain.  Without these extra points the Voronoi network would contain
        points at inf.
        """
        # Label external pores for trimming below
        self['pore.external'] = False
        if len(domain_size) == 1:  # Spherical
            # Trim external Delaunay pores
            r = sp.sqrt(sp.sum(self['pore.coords']**2, axis=1))
            Ps = (r > domain_size)*self['pore.delaunay']
            self['pore.external'][Ps] = True
            # Trim external Voronoi pores
            Ps = ~self['pore.external']*self['pore.delaunay']
            Ps = self.find_neighbor_pores(pores=Ps)
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self['pore.external'][Ps] = True
        elif len(domain_size) == 2:  # Cylindrical
            # Trim external Delaunay pores outside radius
            r = sp.sqrt(sp.sum(self['pore.coords'][:, [0, 1]]**2, axis=1))
            Ps = (r > domain_size[0])*self['pore.delaunay']
            self['pore.external'][Ps] = True
            # Trim external Delaunay pores above and below cylinder
            Ps1 = self['pore.coords'][:, 2] > domain_size[1]
            Ps2 = self['pore.coords'][:, 2] < 0
            Ps = self['pore.delaunay']*(Ps1 + Ps2)
            self['pore.external'][Ps] = True
            # Trim external Voronoi pores
            Ps = ~self['pore.external']*self['pore.delaunay']
            Ps = self.find_neighbor_pores(pores=Ps)
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self['pore.external'][Ps] = True
        elif len(domain_size) == 3:  # Rectilinear
            # Trim external Delaunay pores
            Ps1 = sp.any(self['pore.coords'] > domain_size, axis=1)
            Ps2 = sp.any(self['pore.coords'] < [0, 0, 0], axis=1)
            Ps = self['pore.delaunay']*(Ps1 + Ps2)
            self['pore.external'][Ps] = True
            # Trim external Voronoi pores
            Ps = ~self['pore.external']*self['pore.delaunay']
            Ps = self.find_neighbor_pores(pores=Ps)
            Ps = ~self.tomask(pores=Ps)*self['pore.voronoi']
            self['pore.external'][Ps] = True

        # Label the internal pores & throats
        self['pore.internal'] = ~self['pore.external']
        Ps = self.pores('internal')
        Ts = self.find_neighbor_throats(pores=Ps, mode='intersection')
        self['throat.internal'] = self.tomask(throats=Ts)

        # Label boundary pores and throats
        Ps = self.find_neighbor_pores(pores=self.pores('internal'))
        Ps = self.filter_by_label(pores=Ps, labels='delaunay')
        self['pore.boundary'] = self.tomask(pores=Ps)
        self['pore.external'][Ps] = False
        Ts = self.find_neighbor_throats(pores=Ps)
        self['throat.boundary'] = self.tomask(throats=Ts)

        # Trim external pores
        Ps = self.pores('external')
        self.trim(pores=Ps)
        # Trim throats between boundary pores
        Ps = self.pores('boundary')
        Ts = self.find_neighbor_throats(pores=Ps, mode='intersection')
        self.trim(throats=Ts)

        # Move boundary pores to centroid of Voronoi facet
        for P in self.pores('boundary'):
            Ns = self.find_neighbor_pores(pores=P)
            Ns = self.filter_by_label(pores=Ns, labels='voronoi')
            coords = sp.mean(self['pore.coords'][Ns], axis=0)
            self['pore.coords'][P] = coords

        # Label Voronoi surface pores and throats
        Ps = self.pores('boundary')
        Ps = self.find_neighbor_pores(pores=Ps)
        Ps = self.filter_by_label(pores=Ps, labels='voronoi')
        self['pore.surface'] = self.tomask(pores=Ps)
        Ps = self.pores('surface')
        Ts = self.find_neighbor_throats(pores=Ps, mode='intersection')
        self['throat.surface'] = self.tomask(throats=Ts)

        # Clean-up
        del self['pore.external']


    def _find_facets(self):
        r"""
        """
        # Find boundary Delaunay pores using hidden info from trim_domain
        Ps = self.find_neighbor_pores(pores=self.pores('surface'))
        Ps = self.filter_by_label(pores=Ps, labels='delaunay')
        for pore in Ps:
            Ns = self.find_neighbor_pores(pores=pore)
            Ns = self.filter_by_label(pores=Ns, labels='surface')
            pts = self['pore.coords'][Ns]
            # Convert points to vectors
            vs = pts[0] - pts
            # Perform single value decomposition
            vs = vs - sp.mean(vs, axis=0)
            s = sp.linalg.svd(a=vs, overwrite_a=True, compute_uv=False)
            print(pore, s, sp.any(s < 1e-10))



