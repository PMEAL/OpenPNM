r"""
===============================================================================
Submodule -- miscillaneous
===============================================================================

Models for applying basic phase properties

"""
import scipy as sp

def constant(phase,value,**kwargs):
    r"""
    Assigns specified constant value
    """
    temp = sp.ones(sp.shape(phase.pores()))*value
    return temp

def random(phase,seed=None,**kwargs):
    r"""
    Assigns specified constant value
    """
    sp.random.seed(seed)
    value = sp.random.rand(sp.shape(phase.pores())[0])
    return value
