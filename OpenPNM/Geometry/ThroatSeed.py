
"""
module ThroatSeeds
===============================================================================

"""
import scipy as sp

def constant(network, value,**params):
    r"""
    Assign specified constant value
    """
    network.throat_properties['seed'] = value

def random(**params):
    r"""
    Assign random number to throats
    note: should this be called 'poisson'?  
    """
    print('random: nothing yet')

def neighbor_min(network,**params):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    network.throat_properties['seed'] = sp.amin(network.pore_properties['seed'][network.throat_properties['connections']],1)

def neighbor_max(network,**params):
    r"""
    Adopt the maximum seed value from the neighboring pores
    """
    network.throat_properties['seed'] = sp.amax(network.pore_properties['seed'][network.throat_properties['connections']],1)
