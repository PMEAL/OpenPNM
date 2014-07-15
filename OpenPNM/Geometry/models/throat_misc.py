r"""
===============================================================================
Submodule -- pore_misc
===============================================================================

"""
import scipy as _sp

def constant(throats,value,**kwargs):
    r"""
    Assign specified constant value
    """
    value = _sp.ones(_sp.shape(throats)[0])*value
    return value

def random(throats,seed=None,**kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    _sp.random.seed(seed)
    value=_sp.random.rand(_sp.shape(throats)[0])
    return value

def neighbor(network,
             throats,
             pore_prop,
             mode='min',
             **kwargs):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    pores = network.find_connected_pores(throats)
    pvalues = network[pore_prop][pores]
    if mode == 'min':
        value = _sp.amin(pvalues,axis=1)
    if mode == 'max':
        value = _sp.amax(pvalues,axis=1)
    if mode == 'mean':
        value = _sp.mean(pvalues,axis=1)
    return value