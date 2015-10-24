r"""
===============================================================================
throat_misc --  Miscillaneous and generic functions to apply to throats
===============================================================================

"""
import scipy as _sp


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to throats
    note: should this be called 'poisson'?
    """
    range_size = num_range[1] - num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(geometry.num_throats(),)
    value = value*range_size + range_min
    return value


def neighbor(geometry, network, pore_prop='pore.seed', mode='min', **kwargs):
    r"""
    Adopt a value based on the neighboring pores
    """
    throats = network.throats(geometry.name)
    P12 = network.find_connected_pores(throats)
    pvalues = network[pore_prop][P12]
    if mode == 'min':
        value = _sp.amin(pvalues, axis=1)
    if mode == 'max':
        value = _sp.amax(pvalues, axis=1)
    if mode == 'mean':
        value = _sp.mean(pvalues, axis=1)
    return value
