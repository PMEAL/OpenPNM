r"""
===============================================================================
Submodule -- hydraulic_conductance
===============================================================================

"""

import scipy as _sp

def hagen_poiseuille(physics,
                     phase,
                     network,
                     pore_diameter='pore.diameter',
                     pore_viscosity='pore.viscosity',
                     throat_length='throat.length',
                     throat_diameter='throat.diameter',
                     calc_pore_len=False,
                     **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming cylindrical 
    geometry using the Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
    """
    throats = phase.throats(physics.name)
    mup = phase[pore_viscosity]
    mut = phase.interpolate_data(mup)
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network.find_connected_pores(throats=network.throats(),flatten=0)
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
    else:        
        plen1 = (0.5*pdia[Ps[:,0]])
        plen2 = (0.5*pdia[Ps[:,1]])   
    #Find g for half of pore 1
    gp1 = 2.28*plen1*(pdia[Ps[:,0]])**2/(8*mut)
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    #gp2 = 2.28*(pdia[pores[:,1]]/2)**4/(pdia[pores[:,1]]*mut)
    gp2 = _sp.pi*plen2*(pdia[Ps[:,1]])**2/(128*mut)
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    gt = 2.28*(tdia/2)**4/(2*tlen*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value


