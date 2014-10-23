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
    mup = phase[pore_viscosity]
    mut = phase.interpolate_data(mup)
    pdia = network[pore_diameter]
    #Get the properties of every throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
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
    gp1 = _sp.pi*(pdia[Ps[:,0]])**4/(128*plen1*mut)
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    #gp2 = 2.28*(pdia[pores[:,1]]/2)**4/(pdia[pores[:,1]]*mut)
    gp2 = _sp.pi*(pdia[Ps[:,1]])**4/(128*plen2*mut)
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    #remove any non-positive lengths
    tlen[tlen<=0] = 1e-12
    gt = _sp.pi*(tdia)**4/(128*tlen*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[phase.throats(physics.name)]
    return value


