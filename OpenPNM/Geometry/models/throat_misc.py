r"""
===============================================================================
Submodule -- pore_misc
===============================================================================

"""
import scipy as sp

def constant(throats,value,**kwargs):
    r"""
    Assign specified constant value
    """
    value = sp.ones(sp.shape(throats)[0])*value
    return value

def random(throats,seed=None,**kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    sp.random.seed(seed)
    value=sp.random.rand(sp.shape(throats)[0])
    return value
