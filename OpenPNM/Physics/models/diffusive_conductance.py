r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as _sp

def bulk_diffusion(network,fluid,throats,**kwargs):
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
    cp = fluid['pore.molar_density']
    ct = network.interpolate_data(data=cp)
    DABp = fluid['pore.diffusivity']
    DABt = network.interpolate_data(data=DABp)
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network.find_connected_pores(throats,flatten=0)
    #Find g for half of pore 1
    parea = network['pore.area']
    pdia = network['pore.diameter']
    gp1 = ct*DABt*parea[Ps[:,0]]/(0.5*pdia[Ps[:,0]])
    gp1[~(gp1>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*parea[Ps[:,1]]/(0.5*pdia[Ps[:,0]])
    gp2[~(gp2>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network['throat.area']
    tlen = network['throat.length']
    gt = ct*DABt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value

