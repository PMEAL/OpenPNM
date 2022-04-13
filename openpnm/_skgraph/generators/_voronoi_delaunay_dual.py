import scipy.spatial as sptl
import scipy.sparse as sprs
import numpy as np
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes
from openpnm._skgraph.tools import isoutside, conns_to_am
from openpnm._skgraph import settings


def voronoi_delaunay_dual(points, shape, trim=True):
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

    Returns
    -------
    network : dict
        A dictionary containing 'node.coords' and 'edge.conns'
    vor : Voronoi object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``
    tri : Delaunay object
        The Delaunay triangulation object produced ``scipy.spatial.Delaunay``

    """
    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix
    # Generate a set of base points if scalar was given
    points = tools.parse_points(points=points, shape=shape)
    mask = ~np.all(points == 0, axis=0)

    # Perform tessellations
    vor = sptl.Voronoi(points=points[:, mask])
    tri = sptl.Delaunay(points=points[:, mask])

    # Combine points
    pts_all = np.vstack((vor.points, vor.vertices))
    Nall = np.shape(pts_all)[0]

    # Create adjacency matrix in lil format for quick construction
    am = sprs.lil_matrix((Nall, Nall))
    for ridge in vor.ridge_dict.keys():
        # Make Delaunay-to-Delaunay connections
        for i in ridge:
            am.rows[i].extend([ridge[0], ridge[1]])
        # Get Voronoi vertices for current ridge
        row = vor.ridge_dict[ridge].copy()
        # Index Voronoi vertex numbers by number of Delaunay points
        row = [i + vor.npoints for i in row if i > -1]
        # Make Voronoi-to-Delaunay connections
        for i in ridge:
            am.rows[i].extend(row)
        # Make Voronoi-to-Voronoi connections
        row.append(row[0])
        for i in range(len(row)-1):
            am.rows[row[i]].append(row[i+1])

    # Finalize adjacency matrix by assigning data values
    am.data = am.rows  # Values don't matter, only shape, so use 'rows'
    # Convert to COO format for direct acces to row and col
    am = am.tocoo()
    # Extract rows and cols
    conns = np.vstack((am.row, am.col)).T

    # Convert to sanitized adjacency matrix
    am = conns_to_am(conns)
    # Finally, retreive conns back from am
    conns = np.vstack((am.row, am.col)).T

    # Convert coords to 3D by adding col of 0's if necessary
    coords = np.around(pts_all, decimals=10)
    if np.any(mask == False):
        verts = np.zeros([np.shape(coords)[0], 3])
        for i, col in enumerate(np.where(mask)[0]):
            verts[:, col] = coords[:, i]
    else:
        verts = np.copy(coords)

    # Assign coords and conns to network dict
    network = {}
    network[node_prefix+'.coords'] = verts
    network[edge_prefix+'.conns'] = conns

    n_nodes = verts.shape[0]
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

    # Identify and trim pores outside the domain if requested
    if trim:
        Ps = isoutside(verts, shape=shape)
        network = trim_nodes(g=network, inds=np.where(Ps)[0])

    return network, vor, tri


if __name__ == "__main__":
    settings.node_prefix = 'node'
    settings.edge_prefix = 'edge'
    dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1])
    print(dvd.keys())
    print(dvd['node.coords'].shape)
    print(dvd['edge.conns'].shape)
    dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1], trim=True)
    print(dvd.keys())
    print(dvd['node.coords'].shape)
    print(dvd['edge.conns'].shape)
