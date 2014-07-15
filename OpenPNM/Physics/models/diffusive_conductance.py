r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as _sp

def bulk_diffusion(network,
                   fluid,
                   throats,
                   pore_molar_density='pore.molar_density',
                   pore_diffusivity='pore.diffusivity',
                   pore_area='pore.area',
                   pore_diameter='pore.diameter',
                   throat_area='throat.area',
                   throat_length='throat.length',
                   **kwargs):
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
    #Interpolate pore fluid property values to throats
    cp = fluid[pore_molar_density]
    ct = network.interpolate_data(data=cp)
    DABp = fluid[pore_diffusivity]
    DABt = network.interpolate_data(data=DABp)
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network.find_connected_pores(throats=throats)
    #Find g for half of pore 1
    parea = network[pore_area]
    pdia = network[pore_diameter]
    gp1 = ct*DABt*parea[Ps[:,0]]/(0.5*pdia[Ps[:,0]])
    gp1[~(gp1>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*parea[Ps[:,1]]/(0.5*pdia[Ps[:,1]])
    gp2[~(gp2>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    gt = ct*DABt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value

