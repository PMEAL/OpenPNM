from numpy.linalg import norm as _norm
import numpy as _np


r"""
Pore-scale models related to topology of the network.
"""

__all__ = ["coordination_number",
           "pore_to_pore_distance"]


def coordination_number(target):
    r"""
    """
    network = target.network
    N = network.num_neighbors(pores=network.Ps, flatten=False)
    vals = N[network.pores(target.name)]
    return vals


def pore_to_pore_distance(target):
    r"""
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)
    return values

def nearest_neighbor_distance(target):
    r"""
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)

    data = values
    im = network.create_incidence_matrix()
    values = _np.ones((network.Np, ))*_np.inf
    _np.minimum.at(values, im.row, data[im.col])
    return _np.array(values)

def furthest_neighbor_distance(target):
    r"""
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    values = _norm(C1 - C2, axis=1)

    data = values
    im = network.create_incidence_matrix()
    values = _np.zeros((network.Np, ))
    _np.maximum.at(values, im.row, data[im.col])
    return _np.array(values)
