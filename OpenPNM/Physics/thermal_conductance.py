r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""

import scipy as sp

def constant(physics,
             network,
             geometry,
             fluid,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    fluid.set_throat_data(prop=propname,data=value,locations=geometry)

def na(physics,
       network,
       fluid,
       geometry,
       propname,
       **params):
    
    value = -1
    fluid.set_throat_data(prop=propname,data=value,locations=geometry)

def thermal_fluid(physics,
                  network,
                  fluid,
                  geometry,
                  propname,
                  thermal_conductivity='thermal_conductivity',
                  throat_diameter = 'diameter',
                  throat_length = 'length',
                  pore_diameter = 'diameter',
                  **params):
    r"""
    Calculate the thermal conductance of void conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on size

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : OpenPNM Fluid Object
            The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties have already been determined.

    """
    kp = fluid.get_pore_data(prop=thermal_conductivity)
    kt = network.interpolate_data(kp)

    #Get Nt-by-2 list of pores connected to each throat
    tind = network.get_throat_indices()
    pores = network.find_connected_pores(tind,flatten=0)
    #Find g for half of pore 1
    pdia = network.get_pore_data(prop=pore_diameter)
    gp1 = kt*pdia[pores[:,0]]**2/(0.5*pdia[pores[:,0]])
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = kt*pdia[pores[:,1]]**2/(0.5*pdia[pores[:,1]])
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    tdia = network.get_throat_data(prop=throat_diameter)
    tlen = network.get_throat_data(prop=throat_length)
    gt = kt*tdia**2/tlen
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    fluid.set_throat_data(prop=propname,data=value[geometry.throats()],locations=geometry)

def parallel_resistors(physics,
                       network,
                       geometry,
                       fluid,
                       propname,
                       thermal_conductivity='thermal_conductivity',
                       **params):
    r"""
    Calculate the thermal conductance of solid phase surrounding the void

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid : The fluid of interest

    Notes
    -----
    This function requires that all the necessary phase properties have already been determined.

    .. warning::
       This has not been fully implemented yet

    """
    kp = fluid.get_pore_data(prop=thermal_conductivity)
    kt = network.interpolate_data(kp)
    value = kt #A physical model of parallel resistors representing the solid phase surrouding each pore is required here
    mask = network.get_throat_indices(geometry)
    fluid.set_throat_data(prop=propname,data=value[mask],locations=geometry)


