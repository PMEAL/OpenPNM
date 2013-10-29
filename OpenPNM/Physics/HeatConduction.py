
"""
module HeatConduction
===============================================================================

"""

import OpenPNM
import scipy as sp


def ThermalConductance(network,fluid):
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
    try:
        kp = fluid.pore_conditions['thermal_conductivity']
    except:
        raise Exception('Thermal conductivity of the fluid has not been specified')
    kt = fluid.interpolate_throat_conditions(network,kp)

    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = kt*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = kt*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = kt*network.throat_properties['diameter']**2/(network.throat_properties['length'])
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    fluid.throat_conditions['thermal_conductance'] = g


def ThermalConductanceSolid(network,fluid):
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
    try:
        kp = fluid.pore_conditions['thermal_conductivity']
    except:
        raise Exception('Thermal conductivity of the '+fluid['name']+' has not been specified')
    kt = network.interpolate_throat_values(kp)

    g = kt #A physical model of parallel resistors representing the solid phase surrouding each pore is required here

    fluid.throat_conditions['thermal_conductance'] = g



