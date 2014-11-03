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
                   throat_diameter='throat.diameter',
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
    (1) This function requires that all the necessary phase properties already 
    be calculated.
    
    (2) This function calculates the specified property for the *entire* 
    network then extracts the values for the appropriate throats at the end.
    
    """    
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    #Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    #Get the properties of every throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    #Interpolate pore phase property values to throats
    cp = phase[pore_molar_density]
    ct = phase.interpolate_data(data=cp)
    DABp = phase[pore_diffusivity]
    DABt = phase.interpolate_data(data=DABp)
    if calc_pore_len:
        #Find half-lengths of each pore
        pcoords = network['pore.coords']
        #   Find the pore-to-pore distance
        lengths = _sp.sqrt(_sp.sum(_sp.square(pcoords[Ps[:,0]]-pcoords[Ps[:,1]]),1))
        #update tlen to be the minimum of the original throat length, and the diameter ratio derived length
        tlen = _sp.minimum(lengths-2e-12,tlen)
        #   Calculate the fraction of the remaining distance for each pore
        len_rem = lengths - tlen
        sum_dia = pdia[Ps[:,0]]+pdia[Ps[:,1]]
        plen1 = len_rem*pdia[Ps[:,0]]/sum_dia
        plen2 = len_rem*pdia[Ps[:,1]]/sum_dia
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
    #remove any non-positive lengths
    tlen[tlen<=0] = 1e-12
    gt = ct*DABt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(physics.name)]
    return value


