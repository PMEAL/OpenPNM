import numpy as np
from openpnm._skgraph import settings


__all__ = [
    'join',
]


def join(net1, net2, L_max=0.99):
    r"""
    Joins two networks together topologically including new connections

    Parameters
    ----------
    net1 : dictionary
        A dictionary containing 'node.coords' and 'edge.conns'.
    net2 : dictionary
        A dictionary containing 'node.coords' and 'edge.conns'
    L_max : float
        The maximum distance between vertices below which they are called
        neighbors

    Returns
    -------
    network : dict
        A dictionary containing 'node.coords' vertices from both ``net1`` and
        ``net2``, and ``edge.conns`` with original connections plus new ones
        found during the join process.

    Notes
    -----
    This function uses ``scipy.spatial.KDTree``.

    """
    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix

    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1[node_prefix+'.coords'])
    t2 = KDTree(net2[node_prefix+'.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1[node_prefix+'.coords'].shape[0]
    Np2 = net2[node_prefix+'.coords'].shape[0]
    net3[node_prefix+'.coords'] = np.vstack((net1.pop(node_prefix+'.coords'),
                                             net2.pop(node_prefix+'.coords')))
    net3[edge_prefix+'.conns'] = np.vstack((net1.pop(edge_prefix+'.conns'),
                                            net2.pop(edge_prefix+'.conns') + Np1))
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
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, np.zeros(Np1)*np.nan)
        temp2 = net2.pop(item, np.zeros(Np2)*np.nan)
        net3[item] = np.concatenate((temp1, temp2), axis=0)
    return net3
