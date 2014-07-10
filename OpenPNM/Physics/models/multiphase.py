r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as sp

def conduit_conductance(network,
                        fluid,
                        throats,
                        throat_conductance,
                        throat_occupancy,
                        pore_occupancy,
                        mode='strict',
                        factor=1e-5,
                        **kwargs):
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
    throat_value = fluid[throat_conductance][throats]   
    throat_occupancy = fluid[throat_occupancy][throats] == 1
    pores = network.find_connected_pores(throats,flatten=0)
    pore_occupancy = fluid[pore_occupancy][pores] == 1
                     
    if (mode == 'loose'):
        closed_conduits = -throat_occupancy
    else:
        thoats_closed = -throat_occupancy
        connected_pores = network.find_connected_pores(throats)
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
    return value

