import numpy as np
from openpnm._skgraph import tools


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
    The returned network will use the same node and edge prefixes as ``g1``

    """
    node_prefix_1 = tools.get_node_prefix(g1)
    node_prefix_2 = tools.get_node_prefix(g2)
    if node_prefix_1 != node_prefix_2:
        for item in g2.keys():
            _, prop = item.split('.', 1)
            g2[node_prefix_1 + '.' + item] = g2.pop(item)
    node_prefix = node_prefix_1
    edge_prefix_1 = tools.get_edge_prefix(g1)
    edge_prefix_2 = tools.get_edge_prefix(g2)
    if edge_prefix_1 != edge_prefix_2:
        for item in g2.keys():
            _, prop = item.split('.', 1)
            g2[edge_prefix_1 + '.' + item] = g2.pop(item)
    edge_prefix = edge_prefix_1

    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(g1[node_prefix+'.coords'])
    t2 = KDTree(g2[node_prefix+'.coords'])
    pairs = t1.query_ball_tree(t2, r=L_max)
    # Combine existing network data
    net3 = {}
    Np1 = g1[node_prefix+'.coords'].shape[0]
    Np2 = g2[node_prefix+'.coords'].shape[0]
    Nt1 = g1[edge_prefix+'.conns'].shape[0]
    Nt2 = g2[edge_prefix+'.conns'].shape[0]
    net3[node_prefix+'.coords'] = \
        np.vstack((g1.pop(node_prefix+'.coords'),
                   g2.pop(node_prefix+'.coords')))
    net3[edge_prefix+'.conns'] = \
        np.vstack((g1.pop(edge_prefix+'.conns'),
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
    net3[edge_prefix+'.conns'] = \
        np.vstack((net3.pop(edge_prefix+'.conns'), conns))

    # Finally, expand any other data arrays on given networks
    items = list(g1.keys()) + list(g2.keys())
    N = {node_prefix: [Np1, Np2], edge_prefix: [Nt1, Nt2]}
    for item in items:
        prefix, attr = item.split('.', 1)
        if prefix == node_prefix:
            temp1 = g1.pop(item, None)
            temp2 = g2.pop(item, None)
        if prefix == edge_prefix:
            temp1 = g1.pop(item, None)
            temp2 = g2.pop(item, None)
        if temp1 is None:
            blank = False if temp2.dtype == bool else np.nan
            temp1 = np.zeros(N[prefix][0], dtype=type(blank))*blank
        elif temp2 is None:
            blank = False if temp1.dtype == bool else np.nan
            temp2 = np.zeros(N[prefix][1], dtype=type(blank))*blank
        net3[item] = np.concatenate((temp1, temp2), axis=0)
    return net3
