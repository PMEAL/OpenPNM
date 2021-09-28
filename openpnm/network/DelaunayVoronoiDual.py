import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
from openpnm import topotools
from openpnm.network.generators import voronoi_delaunay_dual, tools
from openpnm.utils import logging
logger = logging.getLogger(__name__)
from openpnm.network import GenericNetwork


class DelaunayVoronoiDual(GenericNetwork):
    r"""
    Combined and interconnected Voronoi and Delaunay tessellations

    Parameters
    ----------
    points : array_like (num_points x 3)
        A list of coordinates for pre-generated points, typically produced
        using ``generate_base_points`` in topotools.  Note that base points
        should extend beyond the domain so that degenerate Voronoi points
        can be trimmed.
    shape : array_like
        The size and shape of the domain used for generating and trimming
        excess points. The coordinates are treated as the outer corner of a
        rectangle [x, y, z] whose opposite corner lies at [0, 0, 0].

        By default, a domain size of [1, 1, 1] is used.  To create a 2D network
        set the Z-dimension to 0.

    name : string
        An optional name for the object to help identify it.  If not given,
        one will be generated.
    project : OpenPNM Project object, optional
        Each OpenPNM object must be part of a *Project*.  If none is supplied
        then one will be created and this Network will be automatically
        assigned to it.  To create a *Project* use ``openpnm.Project()``.

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

    def __init__(self, shape=[1, 1, 1], points=None, **kwargs):
        super().__init__(**kwargs)
        net, vor, tri = voronoi_delaunay_dual(points=points, shape=shape)
        net = tools.add_all_label(net)
        self.update(net)
        self._vor = vor
        self._tri = tri

        self._label_faces()
        # Label all pores and throats by type
        self['pore.delaunay'] = False
        self['pore.delaunay'][0:vor.npoints] = True
        self['pore.voronoi'] = False
        self['pore.voronoi'][vor.npoints:] = True
        # Label throats between Delaunay pores
        self['throat.delaunay'] = False
        Ts = np.all(self['throat.conns'] < vor.npoints, axis=1)
        self['throat.delaunay'][Ts] = True
        # Label throats between Voronoi pores
        self['throat.voronoi'] = False
        Ts = np.all(self['throat.conns'] >= vor.npoints, axis=1)
        self['throat.voronoi'][Ts] = True
        # Label throats connecting a Delaunay and a Voronoi pore
        self['throat.interconnect'] = False
        Ts = self.throats(labels=['delaunay', 'voronoi'], mode='not')
        self['throat.interconnect'][Ts] = True

        # Trim all pores that lie outside of the specified domain
        if self.settings['trim'] == False:
            pass
        else:
            self._trim_external_pores(shape=shape)

    def _trim_external_pores(self, shape):
        r'''
        '''
        # Find all pores within the domain
        Ps = topotools.isoutside(coords=self['pore.coords'], shape=shape)
        self['pore.external'] = False
        self['pore.external'][Ps] = True

        # Find which internal pores are delaunay
        Ps = (~self['pore.external'])*self['pore.delaunay']

        # Find all pores connected to an internal delaunay pore
        Ps = self.find_neighbor_pores(pores=Ps, include_input=True)

        # Mark them all as keepers
        self['pore.keep'] = False
        self['pore.keep'][Ps] = True

        # Trim all bad pores
        topotools.trim(network=self, pores=~self['pore.keep'])

        # Now label boundary pores
        self['pore.boundary'] = False
        self['pore.boundary'] = self['pore.delaunay']*self['pore.external']

        # Label Voronoi pores on boundary
        Ps = self.find_neighbor_pores(pores=self.pores('boundary'))
        Ps = self['pore.voronoi']*self.tomask(pores=Ps)
        self['pore.boundary'][Ps] = True

        # Label Voronoi and interconnect throats on boundary
        self['throat.boundary'] = False
        Ps = self.pores('boundary')
        Ts = self.find_neighbor_throats(pores=Ps, mode='xnor')
        self['throat.boundary'][Ts] = True

        # Trim throats between Delaunay boundary pores
        Ps = self.pores(labels=['boundary', 'delaunay'], mode='xnor')
        Ts = self.find_neighbor_throats(pores=Ps, mode='xnor')
        topotools.trim(network=self, throats=Ts)

        # Move Delaunay boundary pores to centroid of Voronoi facet
        Ps = self.pores(labels=['boundary', 'delaunay'], mode='xnor')
        for P in Ps:
            Ns = self.find_neighbor_pores(pores=P)
            Ns = Ps = self['pore.voronoi']*self.tomask(pores=Ns)
            coords = np.mean(self['pore.coords'][Ns], axis=0)
            self['pore.coords'][P] = coords

        self['pore.internal'] = ~self['pore.boundary']
        Ps = self.pores('internal')
        Ts = self.find_neighbor_throats(pores=Ps, mode='xnor')
        self['throat.internal'] = False
        self['throat.internal'][Ts] = True

        # Label surface pores and throats between boundary and internal
        Ts = self.throats(['boundary', 'internal'], mode='not')
        self['throat.surface'] = False
        self['throat.surface'][Ts] = True
        surf_pores = self['throat.conns'][Ts].flatten()
        surf_pores = np.unique(surf_pores[~self['pore.boundary'][surf_pores]])
        self['pore.surface'] = False
        self['pore.surface'][surf_pores] = True
        # Clean-up
        del self['pore.external']
        del self['pore.keep']

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


    def _label_faces(self):
        r'''
        Label the pores sitting on the faces of the domain in accordance with
        the conventions used for cubic etc.
        '''
        coords = np.around(self['pore.coords'], decimals=10)
        min_labels = ['front', 'left', 'bottom']
        max_labels = ['back', 'right', 'top']
        min_coords = np.amin(coords, axis=0)
        max_coords = np.amax(coords, axis=0)
        for ax in range(3):
            self['pore.' + min_labels[ax]] = coords[:, ax] == min_coords[ax]
            self['pore.' + max_labels[ax]] = coords[:, ax] == max_coords[ax]
