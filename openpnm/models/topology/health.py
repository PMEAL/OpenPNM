import numpy as np


def isolated_pores(target):
    r"""
    Identify pores that are not connected to any neighbors
    """
    network = target.network
    if 0:
        am = network.create_adjacency_matrix(fmt='lil')
        lens = [len(i) for i in am.rows]
        hits = np.where(lens == 0)
        check = np.zeros(network.Np, dtype=bool)
        check[hits] = True
    if 1:
        check = np.in1d(network.Ps, network.conns.flatten(), invert=True)
    return check


def oversize_throats(target,
                     pore_diameter='pore.diameter',
                     throat_diameter='throat.diameter'):
    r"""
    Identify throats that are larger than at least one connected pore
    """
    network = target.network
    Pdia12 = network[pore_diameter][network.conns]
    Tdia = network[throat_diameter]
    check = np.any(Pdia12 < np.atleast_2d(Tdia).T, axis=1)
    return check


def bidirectional_throats(target):
    r"""
    Identify pores that are connected by more thhan one throat
    """
    network = target.network
    check = network.conns[:, 1] > network.conns[:, 1]
    return check


def headless_throats(target):
    r"""
    """
    net = target.network
    check = np.zeros(net.Nt, dtype=bool)
    hits = np.where(net.conns[:, 0] > (net.Np - 1))
    check[hits] = True
    hits = np.where(net.conns[:, 1] > (net.Np - 1))
    check[hits] = True
    return check
