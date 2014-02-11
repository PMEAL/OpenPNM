
"""
module pore_seeds
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]
propname = propname.split('_')[1]

def constant(geometry,network,value,**params):
    r"""
    Assign specified constant value
    """
    network.set_pore_data(subdomain=geometry,prop=propname,data=value)

def random(geometry,network,**params):
    r"""
    Assign random number to pore bodies
    note: should this be called 'poisson'?  
    """
    Np = network.num_pores()
    value=sp.random.rand(Np)
    network.set_pore_data(subdomain=geometry,prop=propname,data=value)
    