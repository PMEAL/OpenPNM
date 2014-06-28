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
    value = fluid['throat.'+conductance]   
    
    throat_occupancy = list(fluid['throat.occupancy'] == 0)
    connected_pores = network.find_connected_pores(geometry.throats())
    
    if (mode == 'loose'):
        s = throat_occupancy
    else:
        pores_1 = connected_pores[:,0]
        pores_2 = connected_pores[:,1]
        pores_1_occupancy = list(fluid['pore.occupancy'][pores_1] == 0)
        pores_2_occupancy = list(fluid['pore.occupancy'][pores_2] == 0)
        
        if(mode == 'medium'):
            s = throat_occupancy or (pores_1_occupancy and pores_2_occupancy)
            
        if(mode == 'strict'):
            s = pores_1_occupancy or throat_occupancy or pores_2_occupancy
    s = sp.array(s)    
    value = value*(-s) + value*s/1.0e3
        
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

