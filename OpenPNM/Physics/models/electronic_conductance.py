r"""
===============================================================================
Submodule -- electronic_conductance
===============================================================================

"""

import scipy as _sp

def series_resistors(network,fluid,throats,**kwargs):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    sigmat = fluid.get_data(prop='electrical_conductivity',throats='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    throats = network.get_throat_indices()
    pores = network.find_connected_pores(throats,flatten=0)
    #Find g for half of pore 1
    parea = network['pore.area']
    gp1 = sigmat*parea[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = sigmat*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network['throat.diameter']
    tlen = network['throat.length']
    gt = sigmat*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value

