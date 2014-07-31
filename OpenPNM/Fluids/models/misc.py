r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic fluid properties

"""
import scipy as _sp

def constant(fluid,value,**kwargs):
    r"""
    Assigns specified constant value
    """
    temp = _sp.ones(_sp.shape(fluid.pores()))*value
    return temp

def random(fluid,seed=None,**kwargs):
    r"""
    Assigns specified constant value
    """
    _sp.random.seed(seed)
    value = _sp.random.rand(_sp.shape(fluid.pores())[0])
    return value
