r"""
===============================================================================
Submodule -- hydraulic_conductance
===============================================================================

"""

import scipy as _sp

def hagen_poiseuille(network,fluid,throats,**kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming square geometry 
    using a modified Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    mut = fluid.get_data(prop='viscosity',throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    tind = network.throats()
    Ps = network.find_connected_pores(tind,flatten=0)
    #Find g for half of pore 1
    pdia = network['pore.diameter']
    gp1 = 2.28*(pdia[Ps[:,0]])**3/(16*mut)
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    #gp2 = 2.28*(pdia[pores[:,1]]/2)**4/(pdia[pores[:,1]]*mut)
    gp2 = 2.28*(pdia[Ps[:,1]])**3/(16*mut)
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network['throat.diameter']
    tlen = network['throat.length']
    gt = 2.28*(tdia/2)**4/(2*tlen*mut)
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value


