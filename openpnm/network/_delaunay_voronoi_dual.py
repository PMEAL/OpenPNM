import numpy as np
import scipy.spatial as sptl
from openpnm import topotools
from openpnm.network import GenericNetwork
from openpnm._skgraph.generators.tools import parse_points
from openpnm._skgraph.generators import voronoi_delaunay_dual


__all__ = ['DelaunayVoronoiDual']


class DelaunayVoronoiDual(GenericNetwork):
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

    def __init__(self, shape, points, trim=True, **kwargs):
        super().__init__(**kwargs)

        points = parse_points(shape=shape, points=points)

        # Deal with points that are only 2D...they break tessellations
        if points.shape[1] == 3 and len(np.unique(points[:, 2])) == 1:
            points = points[:, :2]

        net, vor, tri = voronoi_delaunay_dual(shape=shape, points=points)
        self.update(net)
        self._vor = vor
        self._tri = tri

        # Trim all pores that lie outside of the specified domain
        if trim:
            self._trim_external_pores(shape=shape)
            self._label_faces()

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

    def _trim_external_pores(self, shape):
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
        Ps = self['pore.voronoi']*self.to_mask(pores=Ps)
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
            Ns = Ps = self['pore.voronoi']*self.to_mask(pores=Ns)
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
        r"""
        Label the pores sitting on the faces of the domain in accordance
        with the conventions used for cubic etc.
        """
        coords = np.around(self['pore.coords'], decimals=10)
        min_labels = ['front', 'left', 'bottom']
        max_labels = ['back', 'right', 'top']
        min_coords = np.amin(coords, axis=0)
        max_coords = np.amax(coords, axis=0)
        for ax in range(3):
            self['pore.' + min_labels[ax]] = coords[:, ax] == min_coords[ax]
            self['pore.' + max_labels[ax]] = coords[:, ax] == max_coords[ax]

    def add_boundary_pores(self, labels=['top', 'bottom', 'front', 'back',
                                         'left', 'right'], offset=None):
        r"""
        Add boundary pores to the specified faces of the network.

        Pores are offset from the faces of the domain.

        Parameters
        ----------
        labels : str or list[str]
            The labels indicating the pores defining each face where
            boundary pores are to be added (e.g. 'left' or ['left', 'right'])
        offset : scalar or array_like
            The spacing of the network (e.g. [1, 1, 1]).  This must be
            given since it can be quite difficult to infer from the
            network, for instance if boundary pores have already added to
            other faces.

        """
        offset = np.array(offset)
        if offset.size == 1:
            offset = np.ones(3)*offset
        for item in labels:
            Ps = self.pores(item)
            coords = np.absolute(self['pore.coords'][Ps])
            axis = np.count_nonzero(np.diff(coords, axis=0), axis=0) == 0
            ax_off = np.array(axis, dtype=int)*offset
            if np.amin(coords) == np.amin(coords[:, np.where(axis)[0]]):
                ax_off = -1*ax_off
            topotools.add_boundary_pores(network=self, pores=Ps, offset=ax_off,
                                         apply_label=item + '_boundary')
