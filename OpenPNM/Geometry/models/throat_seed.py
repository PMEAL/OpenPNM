r"""
===============================================================================
Submodule -- throat_seeds
===============================================================================

"""
import scipy as _sp


def neighbor(network,throats,mode='min',**kwargs):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    pores = network.find_connected_pores(throats)
    pseeds = network['pore.seed'][pores]
    if mode == 'min':
        value = _sp.amin(pseeds,axis=1)
    if mode == 'max':
        value = _sp.amax(pseeds,axis=1)
    if mode == 'mean':
        value = _sp.mean(pseeds,axis=1)
    return value

