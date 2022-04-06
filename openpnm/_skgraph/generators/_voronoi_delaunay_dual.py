import scipy.spatial as sptl
import scipy.sparse as sprs
import numpy as np
from openpnm._skgraph.generators import tools
from openpnm._skgraph.operations import trim_nodes
from openpnm.topotools import isoutside, conns_to_am


def voronoi_delaunay_dual(points, shape, crop=False, node_prefix='node', edge_prefix='edge'):
    r"""
    Generate a dual Voronoi-Delaunay network from given base points

    Parameters
    ----------
    points : array_like or scalar
        The points to be tessellated.  If a scalar is given a set of points
        of that size is generated inside the given ``shape``.
    shape : array_like
        The size of the domain in which the points lie
    crop : bool, optional (default is ``False``)
        If ``True`` then all points lying beyond the given domain shape will
        be removed
    node_prefix : str, optional
        If a custom prefix is used to indicate node arrays, such as ('site', or
        'vertex') it can be specified here.  The defaul it 'node'.
    edge_prefix : str, optional
        If a custom prefix is used to indicate site arrays, such as ('bond', or
        'link') it can be specified here.  The defaul it 'edge'.

    Returns
    -------
    network : dict
        A dictionary containing 'vert.coords' and 'edge.conns'
    vor : Voronoi object
        The Voronoi tessellation object produced by ``scipy.spatial.Voronoi``
    tri : Delaunay object
        The Delaunay triangulation object produced ``scipy.spatial.Delaunay``

    """
    # Generate a set of base points if number was given
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
    # Finally, retrieve conns back from am
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

    # Identify and trim pores outside the domain if requested
    if crop:
        Ps = isoutside(verts, shape=shape)
        network = trim_nodes(g=network, inds=np.where(Ps)[0])

    return network, vor, tri


if __name__ == "__main__":
    dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1])
    print(dvd.keys())
    print(dvd['node.coords'].shape)
    print(dvd['edge.conns'].shape)
    dvd, vor, tri = voronoi_delaunay_dual(points=50, shape=[1, 0, 1], crop=True)
    print(dvd.keys())
    print(dvd['node.coords'].shape)
    print(dvd['edge.conns'].shape)
