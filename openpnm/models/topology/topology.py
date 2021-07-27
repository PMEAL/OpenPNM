r"""
Pore-scale models related to topology of the network

"""

import numpy as np


def coordination_number(target):
    r"""
    Compute the number of direct neighbors for each pore

    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals


def total_length(target):
    r"""
    Compute the pore-to-pore distance for each pair of connected pores
    """
    crds = target.network.coords
    conns = target.network.conns
    L = np.sqrt(np.sum(np.diff(crds[conns], axis=1)**2, axis=2)).flatten()
    return L
