
"""
module ThroatSeeds
===============================================================================

"""
import scipy as sp

def constant(value,**params):
    r"""
    Assign specified constant value
    """
    print('constant')
    return value

def random(**params):
    r"""
    Assign random number to throats
    note: should this be called 'poisson'?  
    """
    print('random: nothing yet')

def neighbor_min(**params):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    print('neighbor_min: nothing yet')

def neighbor_max(**params):
    r"""
    Adopt the maximum seed value from the neighboring pores
    """
    print('neighbor_max: nothing yet')  