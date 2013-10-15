
"""
module HeatConduction
===============================================================================

"""

import OpenPNM
import scipy as sp


def FourierConductivity(network,k):
    r"""
    ---
    """
    gt = k*2*network.throat_properties['diameter']/network.lattice_spacing
    
    return gt
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
    
