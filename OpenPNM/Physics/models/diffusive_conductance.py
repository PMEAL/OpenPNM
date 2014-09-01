r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as _sp

def bulk_diffusion(physics,
                   phase,
                   network,
                   pore_molar_density='pore.molar_density',
                   pore_diffusivity='pore.diffusivity',
                   pore_area='pore.area',
                   pore_diameter='pore.diameter',
                   throat_area='throat.area',
                   throat_length='throat.length',
                   calc_pore_len=False,
                   **kwargs):
    r"""
    Calculate the diffusive conductance of conduits in network, where a 
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
        The phase of interest

    Notes
    -----
    This function requires that all the necessary phase properties already be 
    calculated.

    """    
    throats = phase.throats(physics.name)
    #Interpolate pore phase property values to throats
    cp = phase[pore_molar_density]
    ct = phase.interpolate_data(data=cp)
    DABp = phase[pore_diffusivity]
    DABt = phase.interpolate_data(data=DABp)
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network.find_connected_pores(throats=throats)
    parea = network[pore_area]
    pdia = network[pore_diameter]
    pcoords = network['pore.coords']
    if calc_pore_len:
        #Find half-lengths of each pore
        #   Find the pore-to-pore distance, minus the throat length
        lengths = _sp.sqrt(_sp.sum(_sp.square(pcoords[Ps[:,0]]-pcoords[Ps[:,1]]),1))-network[throat_length][throats]
        #   Calculate the fraction of that distance from the first pore    
        fractions = pdia[Ps[:,0]]/(pdia[Ps[:,0]]+pdia[Ps[:,1]])
        plen1 = lengths*fractions
        plen2 = lengths*(1-fractions)
        #remove any non-positive lengths
        plen1[plen1<=0]=1e-12
        plen2[plen2<=0]=1e-12
    else:        
        plen1 = (0.5*pdia[Ps[:,0]])
        plen2 = (0.5*pdia[Ps[:,1]])  
        #remove any non-positive lengths
        plen1[plen1<=0]=1e-12
        plen2[plen2<=0]=1e-12    
    #Find g for half of pore 1
    gp1 = ct*DABt*parea[Ps[:,0]]/plen1
    gp1[~(gp1>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = ct*DABt*parea[Ps[:,1]]/plen2
    gp2[~(gp2>0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    gt = ct*DABt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value


