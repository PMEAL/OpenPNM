r"""
===============================================================================
Submodule -- pore_misc
===============================================================================

"""
import scipy as sp

def constant(pores,value,**kwargs):
    r"""
    Assign specified constant value
    """
    value = sp.ones(sp.shape(pores)[0])*value
    return value

def random(pores,seed=None,**kwargs):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    sp.random.seed(seed)
    value=sp.random.rand(sp.shape(pores)[0])
    return value
