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
                   occupied_condition = 'occupancy',
                   propname = 'conduit_conductance',
                   mode = 'strict',
                   factor = 1/1e3,
                   **params):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas. 
    
    This method "closes" conduits that are not sufficiently filled with the 
    specified fluid by multiplying the original conductance by a very small *factor*.

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
        The fluid of interest
    
    geometry : OpenPNM Geometry Object
        The geometry containing the conduits to be updated
        
    conductance : str
        The name of the conductance that should be updated
        
    occupied_condition : 'occupancy'
        The name of the pore and throat property that dictates whether conduit is "closed" or not
        
    propname : 'conduit_conductance'
        The name of the new throat property created by this method
        
    mode : 'strict' or 'medium' or 'loose'
        How agressive the method should be in "closing" conduits. 
        'strict' implies that if any pore or throat in the conduit is unoccupied by the given fluid, the conduit is closed.
        'medium' implies that if either the throat or both pores are unoccupied, the conduit is closed
        'loose' will only close the conduit if the throat is unoccupied.
        
    factor : 1/1e3
        The "closing factor" which becomes multiplied to the original conduit's conductance to severely limit transport.

    Notes
    -----
    This function requires that all the necessary fluid properties already be 
    calculated.

    """
    throats = geometry.throats()
    throat_value = fluid['throat.'+conductance][throats]   
    
    throat_occupancy = fluid['throat.'+occupied_condition][throats] == 1
    pore_occupancy = fluid['pore.'+occupied_condition] == 1
                     
    if (mode == 'loose'):
        closed_conduits = -throat_occupancy
    else:
        thoats_closed = -throat_occupancy[throats]
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
    value = throat_value*open_conduits + throat_value*closed_conduits*factor
    fluid.set_data(prop=propname,throats=geometry.throats(),data=value)

