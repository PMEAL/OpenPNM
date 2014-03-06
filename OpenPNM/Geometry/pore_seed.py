
"""
module pore_seeds
===============================================================================

"""
import scipy as sp

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assign specified constant value
    """
    network.set_pore_data(locations=geometry,prop=propname,data=value)

def random(geometry,
           network,
           propname,
           seed=None,
           **params):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    sp.random.seed(seed)
    Np = network.num_pores(geometry.name)
    value=sp.random.rand(Np)
    network.set_pore_data(locations=geometry,prop=propname,data=value)
    