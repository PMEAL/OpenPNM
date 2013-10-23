
"""
module HeatConduction
===============================================================================

"""

import OpenPNM
import scipy as sp


def ThermalConductance(network,fluid_name):
    r"""
    Calculate the thermal conductance of void conduits in network ( 1/2 pore - full throat - 1/2 pore ) based on size

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid_name : string
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary fluid properties have already been determined.

    """
    try:
        kp = network.pore_conditions['thermal_conductivity'+'_'+fluid_name]
    except:
        raise Exception('Thermal conductivity of the '+fluid_name+' has not been specified')
    kt = network.interpolate_throat_values(kp)

    #Get Nt-by-2 list of pores connected to each throat
    pores = network.get_connected_pores(network.throat_properties['numbering'],flatten=0)
    #Find g for half of pore 1
    gp1 = kp*network.pore_properties['diameter'][pores[:,0]]**2/(network.pore_properties['diameter'][pores[:,0]]/2)
    gp1[~(gp1>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for half of pore 2
    gp2 = kp*network.pore_properties['diameter'][pores[:,1]]**2/(network.pore_properties['diameter'][pores[:,1]]/2)
    gp2[~(gp2>0)] = sp.inf #Set 0 conductance pores (boundaries) to inf
    #Find g for full throat
    gt = kt*network.throat_properties['diameter']**2/(network.throat_properties['length'])
    g = (1/gt + 1/gp1 + 1/gp2)**(-1)
    network.throat_conditions['thermal_conductance'+'_'+fluid_name] = g


def ThermalConductanceSolid(network,fluid_name):
    r"""
    Calculate the thermal conductance of solid phase surrounding the void

    Parameters
    ----------
    network : OpenPNM Network Object

    fluid_name : string
        The fluid of interest

    Notes
    -----
    This function requires that all the necessary phase properties have already been determined.

    .. warning::
       This has not been implemented yet

    """
    try:
        kp = network.pore_conditions['thermal_conductivity'+'_'+fluid_name]
    except:
        raise Exception('Thermal conductivity of the '+fluid_name+' has not been specified')
    kt = network.interpolate_throat_values(kp)

    g = kt #A physical model of parallel resistors representing the solid phase surrouding each pore is required here

    network.throat_conditions['thermal_conductance'+'_'+fluid_name] = g

def CarbonPaper_InPlane_ThermalConductivity_without_Teflon(network,density=355):

    k = (-4.91*-11*network.pore_properties['temperature']**3 + 1.42e-8*network.pore_properties['temperature']**2\
        -1.46e-6*network.pore_properties['temperature'] + 8.91e-5)*density\
        *OpenPNM.Physics.HeatConduction.CarbonPaper_HeatCapacity_with_Teflon(network,0)
    return k
def CarbonPaper_InPlane_ThermalConductivity_with_Teflon(network):

    k = -7.166e-6*network.pore_properties['temperature']**3 + 2.24e-3*network.pore_properties['temperature']**2\
        -0.237*network.pore_properties['temperature'] + 20.1
    return k

def CarbonPaper_HeatCapacity_with_Teflon(network,TeflonPercent):
    Cp_Carbon = 1.062e-6*network.pore_properties['temperature']**3 - 2.983e-3*network.pore_properties['temperature']**2\
                + 3.2*network.pore_properties['temperature'] + 639.66
    Cp_PTFE = 4*network.pore_properties['temperature'] + 1000
    Cp_total = TeflonPercent*Cp_PTFE + (1-TeflonPercent)*Cp_Carbon
    return Cp_total

