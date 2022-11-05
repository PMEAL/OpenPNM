import scipy.spatial as sptl
import scipy.sparse as sprs
import numpy as np
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes
from openpnm._skgraph.tools import isoutside, conns_to_am
from openpnm._skgraph.queries import find_neighbor_nodes


def voronoi_delaunay_dual(points, shape, trim=True, reflect=True,
                          node_prefix='node', edge_prefix='edge'):
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
    # Generate a set of base points if scalar was given
    points = tools.parse_points(points=points, shape=shape, reflect=reflect)
    # Generate mask to remove any dims with all 0's
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

    # Convert coords to 3D if necessary
    # Rounding is crucial since some voronoi verts endup outside domain
    pts_all = np.around(pts_all, decimals=10)
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

    # Identify and trim nodes outside the domain if requested
    if trim:
        inside_all = ~isoutside(network, shape=shape)
        inside_delaunay = inside_all*network[node_prefix+'.delaunay']
        outside_delaunay = (~inside_all)*network[node_prefix+'.delaunay']
        neighbors = find_neighbor_nodes(network=network,
                                        inds=np.where(inside_delaunay)[0],
                                        include_input=True)
        trim = np.ones([network[node_prefix+'.coords'].shape[0], ], dtype=bool)
        trim[neighbors] = False  # Keep all neighbors to internal delaunay nodes
        trim[outside_delaunay] = True  # Re-add external delaunay nodes to trim
        network = trim_nodes(network=network, inds=np.where(trim)[0])

    return network, vor, tri


if __name__ == "__main__":
    from openpnm._skgraph.visualization import plot_edges

    # dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1])
    # print(dvd.keys())
    # print(dvd['node.coords'].shape)
    # print(dvd['edge.conns'].shape)
    # dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1], trim=True)
    # print(dvd.keys())
    # print(dvd['node.coords'].shape)
    # print(dvd['edge.conns'].shape)
    shape = [1]
    pts = tools.parse_points(points=1000, shape=shape, reflect=True)
    vn, vor, tri = voronoi_delaunay_dual(points=pts, shape=shape, trim=True)
    plot_edges(vn)
