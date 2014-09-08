r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""

import scipy as _sp

def series_resistors(physics,
                     phase,
                     network,
                     thermal_conductivity = 'pore.thermal_conductivity',
                     pore_diameter = 'pore.diameter',
                     pore_area = 'pore.area',
                     throat_area = 'throat.area',
                     throat_length = 'throat.length',
                     **kwargs):
    r"""
    Calculate the thermal conductance of void conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on size

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object
            The phase of interest

    """
    throats = phase.throats(physics.name)
    kp = phase[thermal_conductivity]
    kt = phase.interpolate_data(kp)
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.find_connected_pores(network.throats(),flatten=0)
    #Find g for half of pore 1
    pdia = network[pore_diameter]
    parea = network[pore_area]
    gp1 = kt*parea[pores[:,0]]/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = kt*parea[pores[:,1]]/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tarea = network[throat_area]
    tlen = network[throat_length]
    gt = kt*tarea/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value
    


