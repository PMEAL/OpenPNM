r"""
===============================================================================
Submodule -- pore_seeds
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
    value=_sp.random.rand(_sp.shape(pores)[0])
    return value
    