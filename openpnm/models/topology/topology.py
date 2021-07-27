r"""
Pore-scale models related to topology of the network.
"""

import numpy as np

__all__ = ["coordination_number"]


def coordination_number(target):
    r"""
    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals

def total_length(target):
    r"""
    """
    crds = target.network.coords
    conns = target.network.conns
    L = np.sqrt(np.sum(np.diff(crds[conns], axis=1)**2, axis=2)).flatten()
    return L
