import scipy.spatial as sptl
import numpy as np


def gabriel(delaunay):
    # Find centroid or midpoint of each edge in conns
    c = delaunay['pore.coords'][delaunay['throat.conns']]
    m = (c[:, 0, :] + c[:, 1, :])/2
    # Find the radius sphere between each pair of nodes
    r = np.sqrt(np.sum((c[:, 0, :] - c[:, 1, :])**2, axis=1))/2
    # Use the kd-tree function in Scipy's spatial module
    tree = sptl.cKDTree(delaunay['pore.coords'])
    # Find the nearest point for each midpoint
    n = tree.query(x=m, k=1)[0]
    # If nearest point to m is at a distance r, then the edge is a Gabriel edge
    g = n >= r*(0.999)  # The factor is to avoid precision errors in the distances
    d = {}
    d.update(delaunay)
    # Reduce the connectivity to all True values found in g
    d['throat.conns'] = delaunay['throat.conns'][g]
    d['pore.coords'] = delaunay['pore.coords']
    return d
