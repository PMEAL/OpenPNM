r"""
===============================================================================
Submodule -- pore_misc
===============================================================================

"""
import scipy as _sp

def random(geometry,
           seed=None,
           **kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    _sp.random.seed(seed)
    value = _sp.random.rand(geometry.num_throats(),)
    return value

def neighbor(network,
             throats,
             pore_prop='pore.seed',
             mode='min',
             **kwargs):
    r"""
    Adopt a value based on the neighboring pores
    """
    P12 = network.find_connected_pores(throats)
    pvalues = network[pore_prop][P12]
    if mode == 'min':
        value = _sp.amin(pvalues,axis=1)
    if mode == 'max':
        value = _sp.amax(pvalues,axis=1)
    if mode == 'mean':
        value = _sp.mean(pvalues,axis=1)
    return value