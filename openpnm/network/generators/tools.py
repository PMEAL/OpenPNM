import numpy as np


def trim(network, pores=None, throats=None):
    r"""
    Remove given pores or throats from a network

    Parameters
    ----------
    network : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns', describing
        the network
    pores : array_like
        The list of pores to trim.
    thraots : array_like
        The list of throats to trim

    Notes
    -----
    Only one of ``pores`` or ``throats`` can be given.  To trim both types,
    call the fuction twice.

    """
    if (pores is None) and (throats is None):
        raise Exception('Cannot trim pores and throats at the same time')
    if throats is not None:
        throats = np.atleast_1d(throats)
        keep = np.ones(network['throat.conns'].shape[0], dtype=bool)
        keep[throats] = False
        for item in network.keys():
            if item.startswith('throat'):
                network[item] = network[item][keep]
    elif pores is not None:
        pores = np.atleast_1d(pores)
        if pores.dtype == bool:
            pores = np.where(pores)[0]
        keep = np.ones(network['pore.coords'].shape[0], dtype=bool)
        keep[pores] = False
        for item in network.keys():
            if item.startswith('pore'):
                network[item] = network[item][keep]
        # Remove throats
        throats = np.any(np.isin(network['throat.conns'], pores), axis=1)
        network = trim(network, throats=throats)
        # Renumber throat conns
        remapping = np.cumsum(keep) - 1
        network['throat.conns'] = remapping[network['throat.conns']]
    return network


def join(net1, net2, L_max=0.99):
    r"""
    Joins two networks together topologically

    Parameters
    ----------
    net1 : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns'.
    net2 : dictionary
        A dictionary containing 'pore.coords' and 'throat.conns'
    L_max : float
        The maximum distance between pores below which they are called
        neighbors

    Returns
    -------
    net3 : dictionary
        A dictionary containing 'pore.coords' with pores from both ``net1`` and
        ``net2``, and ``throat.conns`` with original connections plus new ones
        found during the join process.

    Notes
    -----
    This function uses ``scipy.spatial.KDTree``.

    """
    # Perform neighbor query
    from scipy.spatial import KDTree
    t1 = KDTree(net1['pore.coords'])
    t2 = KDTree(net2['pore.coords'])
    pairs = t1.query_ball_tree(t2, r=0.99)
    # Combine existing network data
    net3 = {}
    Np1 = net1['pore.coords'].shape[0]
    Np2 = net2['pore.coords'].shape[0]
    net3['pore.coords'] = np.vstack((net1.pop('pore.coords'),
                                     net2.pop('pore.coords')))
    net3['throat.conns'] = np.vstack((net1.pop('throat.conns'),
                                      net2.pop('throat.conns') + Np1))
    # Convert kdtree result into new connections
    nnz = sum([len(row) for row in pairs])
    conns = np.zeros((nnz, 2), dtype=int)
    i = 0
    for j, row in enumerate(pairs):
        for col in row:
            conns[i, :] = j, col + Np1
            i += 1
    # Add new connections to network
    net3['throat.conns'] = np.vstack((net3.pop('throat.conns'), conns))
    # Finally, expand any other data arrays on given networks
    keys = set(net1.keys()).union(net2.keys())
    for item in keys:
        temp1 = net1.pop(item, None)
        temp2 = net2.pop(item, None)
        if temp1 is None:
            temp1 = np.zeros(Np1, dtype=temp2.dtype)
        else:
            temp2 = np.zeros(Np2, dtype=temp1.dtype)
        net3[item] = np.hstack((temp1, temp2)).T
    return net3
