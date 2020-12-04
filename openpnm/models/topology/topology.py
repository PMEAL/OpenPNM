r"""
Pore-scale models related to topology of the network.
"""
__all__ = ["coordination_number"]


def coordination_number(target):
    r"""
    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals
