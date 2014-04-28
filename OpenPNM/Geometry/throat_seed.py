
"""
module throat_seeds
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
    network.set_data(prop=propname,throats=geometry.throats,data=value)

def random(geometry,
           network,
           propname,
           seed=None,
           **params):
    r"""
    Assign random number to throats
    note: should this be called 'poisson'?  
    """
    print('random: nothing yet')

def neighbor_min(geometry,
                 network,
                 propname,
                 pore_seed='seed',
                 **params):
    r"""
    Adopt the minimum seed value from the neighboring pores
    """
    pseeds = network.get_data(prop=pore_seed,pores='all')
    conns = network.get_data(prop='connections',throats=geometry.throats)
    value = sp.amin(pseeds[conns],axis=1)
    network.set_data(prop=propname,throats=geometry.throats,data=value)

def neighbor_max(geometry,
                 network,
                 propname,
                 pore_seed='seed',
                 **params):
    r"""
    Adopt the maximum seed value from the neighboring pores
    """
    pseeds = network.get_data(prop=pore_seed,pores='all')
    conns = network.get_data(prop='connections',throats=geometry.throats)
    value = sp.amax(pseeds[conns],axis=1)
    network.set_data(prop=propname,throats=geometry.throats,data=value)
    
    
    