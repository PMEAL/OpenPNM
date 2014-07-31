r"""
===============================================================================
Submodule -- electronic_conductance
===============================================================================

"""

import scipy as _sp

def series_resistors(network,
                     fluid,
                     throats,
                     pore_electrical_conductivity='pore.electrical_conductivity',
                     pore_area='pore.area',
                     pore_diameter='pore.diameter',
                     throat_area='throat.area',
                     throat_length='throat.length',
                     **kwargs):
    r"""
    Calculates the electronic conductance of throat assuming cylindrical geometry

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
    """
    sigmap = fluid[pore_electrical_conductivity]
    sigmat = fluid.interpolate_data(sigmap)
    #Get Nt-by-2 list of pores connected to each throat
    throats = network.throats()
    pores = network.find_connected_pores(throats,flatten=0)
    #Find g for half of pore 1
    parea = network[pore_area]
    pdia = network[pore_diameter]
    gp1 = sigmat*parea[pores[:,0]]/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = sigmat*parea[pores[:,1]]/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    gt = sigmat*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value

