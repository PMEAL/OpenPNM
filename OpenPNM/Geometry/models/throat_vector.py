r"""
===============================================================================
Submodule -- throat_vector
===============================================================================

"""
import scipy as _sp


def pore_to_pore(geometry, network, **kwargs):
    r"""
    Calculates throat vector as straight path between connected pores.

    Notes
    -----
    There is an important impicit assumption here: the positive direction is
    taken as the direction from the pore with the lower index to the higher.
    This corresponds to the pores in the 1st and 2nd columns of the
    'throat.conns' array as stored on the etwork.
    """
    throats = network.throats(geometry.name)
    pores = network.find_connected_pores(throats, flatten=False)
    C0 = network['pore.coords'][pores, 0]
    C1 = network['pore.coords'][pores, 1]
    V = C1 - C0
    L = _sp.array(_sp.sqrt(_sp.sum(V[:, :]**2, axis=1)), ndmin=1)
    value = V/_sp.array(L, ndmin=2).T
    return value
