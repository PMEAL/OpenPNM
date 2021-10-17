import numpy as np
import scipy.sparse as sprs
from openpnm.network import GenericNetwork
from openpnm import topotools
from openpnm.network.generators import voronoi_delaunay_dual, tools
from openpnm.utils import logging
logger = logging.getLogger(__name__)


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

    def __init__(self, points, shape=[1, 1, 1], crop=True, **kwargs):
        super().__init__(**kwargs)
        net, vor, tri = voronoi_delaunay_dual(points=points, shape=shape, crop=False)

        # Label pores and throats
        Np = net['pore.coords'].shape[0]
        net['pore.delaunay'] = np.zeros(Np, dtype=bool)
        net['pore.delaunay'][0:vor.npoints] = True
        net['pore.voronoi'] = np.zeros(Np, dtype=bool)
        net['pore.voronoi'][vor.npoints:] = True
        # Label throats between Delaunay pores
        net['throat.delaunay'] = np.all(net['throat.conns'] < vor.npoints, axis=1)
        # Label throats between Voronoi pores
        net['throat.voronoi'] = np.all(net['throat.conns'] >= vor.npoints, axis=1)
        net['throat.interconnect'] = (~net['throat.delaunay']) * (~net['throat.voronoi'])
        net = tools.add_all_label(net)

        # Convert to an OpenPNM network object
        self.update(net)
        self._vor = vor
        self._tri = tri

        thresh =1e-6
        Ps = topotools.isoutside(coords=self.coords, shape=shape, thresh=thresh)
        self['pore.outside'] = Ps
        if crop:
            # Drop pores that are outside the domain, unless a neighbor is not
            Tdrop = np.all(Ps[self.conns], axis=1)
            topotools.trim(network=self, throats=Tdrop)
            Pdrop = np.where(self.num_neighbors(pores=self.Ps) == 0)[0]
            topotools.trim(network=self, pores=Pdrop)
            Pdrop = self.pores(['voronoi', 'outside'], mode='and')
            topotools.trim(network=self, pores=Pdrop)
            Pmove = self.pores(['delaunay', 'outside'], mode='and')
            for p in Pmove:
                Pn = self.find_neighbor_pores(p)
                Pn = self.filter_by_label(pores=Pn, labels=['voronoi'])
                if len(Pn) > 0:
                    self['pore.coords'][p] = self['pore.coords'][Pn].mean(axis=0)

        # Label faces, which will be nice and flat if crop was True
        topotools.label_faces(network=self)

        # Add some additional labels
        self['pore.boundary'] = np.copy(self['pore.outside'])
        Ps = self.find_neighbor_pores(pores=self.pores('boundary'))
        # Ps = self.filter_by_label(pores=Ps, labels=['delaunay'])
        self['pore.surface'] = self.tomask(pores=Ps)
        Ps = self.pores(['boundary', 'surface'], mode='not')
        self['pore.internal'] = self.tomask(pores=Ps)

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


if __name__ == "__main__":
    # np.random.seed(2)
    points = 1000
    shape = [1, 1, 0]
    points = topotools.generate_base_points(num_points=points, domain_size=shape, reflect=True)
    dvd = DelaunayVoronoiDual(points=points, shape=shape, crop=False)
    print(dvd)
    fig = topotools.plot_connections(dvd, throats=dvd.throats('voronoi'), c='g')
    fig = topotools.plot_connections(dvd, throats=dvd.throats('delaunay'), c='r', ax=fig)
    fig = topotools.plot_connections(dvd, throats=dvd.throats('interconnect'), c='b', ax=fig)
    # fig = topotools.plot_coordinates(dvd, pores=dvd.pores('internal'), c='k', ax=fig, markersize=100)
    fig = topotools.plot_coordinates(dvd, pores=dvd.pores('surface'), c='m', ax=fig, markersize=100)
    fig = topotools.plot_coordinates(dvd, pores=dvd.pores('boundary'), c='c', ax=fig, markersize=100)


