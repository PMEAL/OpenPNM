r"""
===============================================================================
Submodule -- hydraulic_conductance
===============================================================================

"""

import scipy as _sp

def hagen_poiseuille(physics,
                     fluid,
                     network,
                     pore_diameter='pore.diameter',
                     pore_viscosity='pore.viscosity',
                     throat_length='throat.length',
                     throat_diameter='throat.diameter',
                     **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming square geometry 
    using a modified Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    throats = fluid.throats(physics.name)
    mup = fluid[pore_viscosity]
    mut = fluid.interpolate_data(mup)
    #Get Nt-by-2 list of pores connected to each throat
    Ps = network.find_connected_pores(throats=network.throats(),flatten=0)
    #Find g for half of pore 1
    pdia = network[pore_diameter]
    gp1 = 2.28*(pdia[Ps[:,0]])**3/(16*mut)
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    #gp2 = 2.28*(pdia[pores[:,1]]/2)**4/(pdia[pores[:,1]]*mut)
    gp2 = 2.28*(pdia[Ps[:,1]])**3/(16*mut)
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    gt = 2.28*(tdia/2)**4/(2*tlen*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value


