import numpy as np
from openpnm._skgraph import settings


__all__ = [
    'join',
]


def join(g1, g2, L_max=0.99):
    r"""
    Joins two networks together topologically including new connections

    Parameters
    ----------
    g1 : dictionary
        A dictionary containing 'node.coords' and 'edge.conns'.
    g2 : dictionary
        A dictionary containing 'node.coords' and 'edge.conns'
    L_max : float
        The distance between nodes below which they are joined

    Returns
    -------
    network : dict
        A new dictionary containing coords and conns from both ``net1`` and
        ``net2``, plus new connections between the two input networks.

    Notes
    -----
    both ``g1`` and ``g2`` must use the same node and edge prefixes

    """
    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix

    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(g1[node_prefix+'.coords'])
    t2 = KDTree(g2[node_prefix+'.coords'])
    pairs = t1.query_ball_tree(t2, r=L_max)
    # Combine existing network data
    net3 = {}
    Np1 = g1[node_prefix+'.coords'].shape[0]
    Np2 = g2[node_prefix+'.coords'].shape[0]
    net3[node_prefix+'.coords'] = np.vstack((g1.pop(node_prefix+'.coords'),
                                             g2.pop(node_prefix+'.coords')))
    net3[edge_prefix+'.conns'] = np.vstack((g1.pop(edge_prefix+'.conns'),
                                            g2.pop(edge_prefix+'.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3[edge_prefix+'.conns'] = np.vstack((net3.pop(edge_prefix+'.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(g1.keys()).union(g2.keys())
    for item in keys:
        temp1 = g1.pop(item, np.zeros(Np1)*np.nan)
        temp2 = g2.pop(item, np.zeros(Np2)*np.nan)
        net3[item] = np.concatenate((temp1, temp2), axis=0)
    return net3
