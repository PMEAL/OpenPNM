
"""
module pore_seeds
===============================================================================

"""
import scipy as sp

def constant(geometry,network,value,**params):
    r"""
    Assign specified constant value
    """
    network.set_pore_data(prop='seed',data= value)

def random(geometry,network,**params):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    Np = network.get_num_pores()
    network.set_pore_data(prop='seed',data=sp.random.rand(Np))
    