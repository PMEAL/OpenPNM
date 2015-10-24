r"""
===============================================================================
pore_misc -- miscillaneous and generic functions to apply to pores
===============================================================================

"""
import scipy as _sp


def constant(geometry, value, **kwargs):
    r"""
    Assign specified constant value.  This function is redundant and could be
    accomplished with geometry['pore.prop'] = value.
    """
    value = _sp.ones(geometry.num_pores(),)*value
    return value


def random(geometry, seed=None, num_range=[0, 1], **kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?
    """
    range_size = num_range[1]-num_range[0]
    range_min = num_range[0]
    _sp.random.seed(seed=seed)
    value = _sp.random.rand(geometry.num_pores(),)
    value = value*range_size + range_min
    return value


def neighbor(network, geometry, throat_prop='', mode='min', **kwargs):
    r"""
    Adopt the minimum seed value from the neighboring throats
    """
    Ps = geometry.pores()
    data = geometry[throat_prop]
    neighborTs = network.find_neighbor_throats(pores=Ps,
                                               flatten=False,
                                               mode='intersection')
    values = _sp.ones((_sp.shape(Ps)[0],))*_sp.nan
    if mode == 'min':
        for pore in Ps:
            values[pore] = _sp.amin(data[neighborTs[pore]])
    if mode == 'max':
        for pore in Ps:
            values[pore] = _sp.amax(data[neighborTs[pore]])
    if mode == 'mean':
        for pore in Ps:
            values[pore] = _sp.mean(data[neighborTs[pore]])
    return values
