r"""
===============================================================================
Submodule -- pore_misc
===============================================================================

"""
import scipy as _sp

def constant(pores,
             value,
             **kwargs):
    r"""
    Assign specified constant value
    """
    value = _sp.ones(_sp.shape(pores)[0])*value
    return value

def random(pores,
           seed=None,
           **kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    _sp.random.seed(seed)
    value = _sp.random.rand(_sp.shape(pores)[0])
    return value

def neighbor(network,
             pores,
             mode='min',
             pore_seed='pore.seed',
             **kwargs):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    pores = network.find_connected_pores(pores)
    tvalues = network[pore_seed][pores]
    if mode == 'min':
        value = _sp.amin(tvalues,axis=1)
    if mode == 'max':
        value = _sp.amax(tvalues,axis=1)
    if mode == 'mean':
        value = _sp.mean(tvalues,axis=1)
    return value
