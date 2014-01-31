
"""
module throat_seeds
===============================================================================

"""
import scipy as sp

def constant(geometry, network, value,**params):
    r"""
    Assign specified constant value
    """
    network.set_throat_data(prop='seed',data=value)

def random(geometry,network,**params):
    r"""
    Assign random number to throats
    note: should this be called 'poisson'?  
    """
    print('random: nothing yet')

def neighbor_min(geometry,network,**params):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    network.set_throat_data(prop='seed',data=sp.amin(network.get_pore_data(prop='seed')[network.get_throat_data(prop='connections')],1))

def neighbor_max(geometry,network,**params):
    r"""
    Adopt the maximum seed value from the neighboring pores
    """
    network.set_throat_data(prop='seed',data=sp.amax(network.get_pore_data(prop='seed')[network.get_throat_data(prop='connections')],1))
