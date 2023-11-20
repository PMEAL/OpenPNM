import scipy.spatial as sptl
import scipy.sparse as sprs
import numpy as np
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes
from openpnm._skgraph.tools import isoutside, conns_to_am
from openpnm._skgraph.queries import find_neighbor_nodes


def voronoi_delaunay_dual(
    points,
    shape,
    trim=True,
    reflect=True,
    f=1,
    relaxation=0,
    node_prefix='node',
    edge_prefix='edge',
    return_tri=False,
):
    r"""
    Generate a dual Voronoi-Delaunay network from given base points

    Parameters
    ----------
    points : array_like or scalar
        The points to be tessellated.  If a scalar is given a set of points
        of that size is generated inside the given ``shape``.
    shape : array_like
        The size of the domain in which the points lie
    trim : bool, optional
        If ``True`` (default) then all points lying beyond the given domain
        shape will be removed
    reflect : bool, optionl
        If ``True`` (Default) then points are reflected across each face of the
        domain prior to performing the tessellation.
    f : float
        The fraction of points which should be reflected.  The default is 1 which
        reflects all the points in the domain, but this can lead to a lot of
        unnecessary points, so setting to 0.1 or 0.2 helps speed, but risks that
        the tessellation may not have smooth faces if not enough points are
        reflected.
    relaxation : int, optional (default = 0)
        The number of iterations to use for relaxing the base points. This is
        sometimes called `Lloyd's algorithm
        <https://en.wikipedia.org/wiki/Lloyd%27s_algorithm>`_. This function computes
        the new base points as the simple average of the Voronoi vertices instead
        of rigorously finding the center of mass, which is quite time consuming.
        To use the rigorous method, call the ``lloyd_relaxation`` function manually
        to obtain relaxed points, then pass the points directly to this funcion.
        The results are quite stable after only a few iterations.

    Returns
    -------
    network : dict
        A dictionary containing '<node_prefix>.coords' and '<edge_prefix>.conns'
    vor : Voronoi object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``
    tri : Delaunay object
        The Delaunay triangulation object produced by ``scipy.spatial.Delaunay``

    """
    # Generate a set of base points if scalar was given
    points = tools.parse_points(
        points=points,
        shape=shape,
        reflect=reflect,
        f=f,
    )

    # Generate mask to remove any dims with all 0's
    mask = ~np.all(points == 0, axis=0)

    # Perform tessellations
    vor = sptl.Voronoi(points=points[:, mask])
    for _ in range(relaxation):
        points = tools.lloyd_relaxation(vor, mode='fast')
        vor = sptl.Voronoi(points=points[:, mask])

    # Collect delaunay edges
    conns_del = vor.ridge_points
    # Deal with voronoi edges
    v = vor.ridge_vertices.copy()  # Assuming 'ridge' means facet between regions
    # Add row [0] to close the facet on itself, add -1 to break connection to
    # next facet in list as connections with -1 get deleted later
    _ = [row.extend([row[0], -1]) for row in v]
    v = np.hstack(v)
    conns_vor = np.vstack((v[:-1], v[1:])).T
    mask = np.any(conns_vor < 0, axis=1)
    conns_vor = conns_vor[~mask] + vor.npoints
    # Finally, get interconnecting edges
    idx = [vor.regions[vor.point_region[i]] for i in range(0, len(vor.regions)-1)]
    conns_inter = [([i]*len(idx[i]), idx[i]) for i in range(0, len(idx))]
    conns_inter = np.hstack(conns_inter).astype(int).T
    mask = np.any(conns_inter < 0, axis=1)
    conns_inter = conns_inter[~mask, :] + np.array([0, vor.npoints], dtype=int)
    conns = np.vstack((conns_del, conns_vor, conns_inter))

    # Tidy up
    am = conns_to_am(conns)
    conns = np.vstack((am.row, am.col)).T

    # Combine delaunay and voronoi points
    pts_all = np.vstack((vor.points, vor.vertices))
    # Rounding is crucial since some voronoi verts endup outside domain
    pts_all = np.around(pts_all, decimals=10)
    # Convert coords to 3D if necessary
    mask = ~np.all(points == 0, axis=0)
    if mask.sum() < 3:
        coords = np.zeros([pts_all.shape[0], 3], dtype=float)
        coords[:, mask] = pts_all
    else:
        coords = pts_all

    # Assign coords and conns to network dict
    network = {}
    network[node_prefix+'.coords'] = coords
    network[edge_prefix+'.conns'] = conns

    n_nodes = coords.shape[0]
    n_edges = conns.shape[0]

    # Label all pores and throats by type
    network[node_prefix+'.delaunay'] = np.zeros(n_nodes, dtype=bool)
    network[node_prefix+'.delaunay'][0:vor.npoints] = True
    network[node_prefix+'.voronoi'] = np.zeros(n_nodes, dtype=bool)
    network[node_prefix+'.voronoi'][vor.npoints:] = True
    # Label throats between Delaunay pores
    network[edge_prefix+'.delaunay'] = np.zeros(n_edges, dtype=bool)
    Ts = np.all(network[edge_prefix+'.conns'] < vor.npoints, axis=1)
    network[edge_prefix+'.delaunay'][Ts] = True
    # Label throats between Voronoi pores
    network[edge_prefix+'.voronoi'] = np.zeros(n_edges, dtype=bool)
    Ts = np.all(network[edge_prefix+'.conns'] >= vor.npoints, axis=1)
    network[edge_prefix+'.voronoi'][Ts] = True
    # Label throats connecting a Delaunay and a Voronoi pore
    Ts = np.sum(network[node_prefix+'.delaunay'][conns].astype(int), axis=1) == 1
    network[edge_prefix+'.interconnect'] = Ts

    if trim:
        # Find all delaunay nodes outside the domain
        Ps = isoutside(network=network, shape=shape)*network[node_prefix+'.delaunay']
        if np.any(Ps):  # only occurs if points were reflected
            # Find voronoi nodes connected to these and mark them as surface nodes
            inds = np.where(Ps)[0]
            Ns = find_neighbor_nodes(network=network, inds=inds)
            network[node_prefix+'.surface'] = np.zeros(n_nodes, dtype=bool)
            network[node_prefix+'.surface'][Ns] = True
            Ps = isoutside(network=network, shape=shape)
            inds = np.where(Ps)[0]
            network = trim_nodes(network=network, inds=inds)
        else:
            trim = isoutside(network=network, shape=shape)
            inds = np.where(trim)[0]
            network = trim_nodes(network=network, inds=inds)

    if return_tri:
        tri = sptl.Delaunay(points=points[:, mask])
    else:
        tri = None
    return network, vor, tri


if __name__ == "__main__":
    import openpnm as op
    pn, vor, tri = voronoi_delaunay_dual(
        points=500,
        shape=[1, 1, 0],
        trim=True,
        reflect=True,
        f=0.2,
        relaxation=0,
        node_prefix='pore',
        edge_prefix='throat',
    )
    net = op.network.Network()
    net.update(pn)
    h = None
    h = op.visualization.plot_connections(
        net, throats=pn['throat.delaunay'], c='b', ax=h)
    h = op.visualization.plot_connections(
        net, throats=pn['throat.voronoi'], c='g', ax=h)
    h = op.visualization.plot_connections(
        net, throats=pn['throat.interconnect'], c='r', ax=h)
    h = op.visualization.plot_coordinates(
        net, pores=pn['pore.voronoi'], c='g', markersize=150, ax=h)
    h = op.visualization.plot_coordinates(
        net, pores=pn['pore.delaunay'], c='b', markersize=150, ax=h)
