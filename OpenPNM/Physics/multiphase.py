r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as sp

def constant(physics,
             network,
             geometry,
             fluid,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def na(physics,
       network,
       geometry,
       fluid,
       propname,
       **params):
    r"""
    """
    value = -1
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

def conduit_conductance(physics,
                   network,
                   fluid,
                   geometry,
                   conductance,
                   shape = 'circular',
                   propname = 'conduit_conductance',
                   mode = 'strict',
                   factor = 1/1e3,
                   **params):
    r"""
    Calculate the diffusive conductance of conduits in network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties already be 
    calculated.

    """
    throat_value = fluid['throat.'+conductance]   
    
    throat_occupancy = fluid['throat.occupancy'] == 1
    pore_occupancy = fluid['pore.occupancy'] == 1
                     
    if (mode == 'loose'):
        closed_conduits = -throat_occupancy
    else:
        thoats_closed = -throat_occupancy
        connected_pores = network.find_connected_pores(geometry.throats())
        pores_1 = connected_pores[:,0]
        pores_2 = connected_pores[:,1]
        pores_1_closed = -pore_occupancy[pores_1]
        pores_2_closed = -pore_occupancy[pores_2]
        
        if(mode == 'medium'):
            closed_conduits = thoats_closed | (pores_1_closed & pores_2_closed)
            
        if(mode == 'strict'):
            closed_conduits = pores_1_closed | thoats_closed | pores_2_closed
    open_conduits = -closed_conduits
    print(closed_conduits.sum())
    value = throat_value*open_conduits + throat_value*closed_conduits/1.0e3
    
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

