r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""

import scipy as _sp

def thermal_fluid(physics,
                  fluid,
                  network,
                  thermal_conductivity = 'pore.thermal_conductivity',
                  pore_diameter = 'pore.diameter',
                  throat_diameter = 'throat.diameter',
                  throat_length = 'throat.length',
                  **kwargs):
    r"""
    Calculate the thermal conductance of void conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on size

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
            The fluid of interest

    """
    throats = fluid.throats(physics.name)
    kt = fluid.get_data(prop=thermal_conductivity,pores='all',mode='interpolate')
    #Get Nt-by-2 list of pores connected to each throat
    pores = network.find_connected_pores(network.throats(),flatten=0)
    #Find g for half of pore 1
    pdia = network[pore_diameter]
    gp1 = kt*pdia[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = kt*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = _sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    gt = kt*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = value[throats]
    return value
    


