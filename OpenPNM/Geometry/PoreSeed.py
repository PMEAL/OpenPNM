
"""
module PoreSeeds
===============================================================================

"""
import scipy as sp

def constant(network,value,**params):
    r"""
    Assign specified constant value
    """
    network.pore_properties['seed'] = value

def random(network,**params):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    Np = network.get_num_pores()
    network.pore_properties['seed'] = sp.random.rand(Np)
    