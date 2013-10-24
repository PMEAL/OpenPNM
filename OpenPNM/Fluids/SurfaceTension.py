
"""
module SurfaceTension
===============================================================================

"""

import scipy as sp
import OpenPNM

def constant(fluid, value=0.072,**params):
    return 0.072

def na(fluid,**params):
    return 'n/a'

def Eotvos(fluid, k=2.1e-7, **params):
    r"""
    """
    Tc = fluid._fluid_recipe['Tc']
    T = fluid.pore_properties['temperature']
    Vm = 1/fluid.pore_properties['molar_volume']
    sigma = k*(Tc-T)/(Vm**(2/3))
    return sigma

def GuggenheimKatayama(fluid, K2=1, n=1.222, **params):
    r"""
    """
    T = fluid.pore_properties['temperature']
    Pc = fluid._fluid_recipe['Pc']
    Tc = fluid._fluid_recipe['Tc']
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    sigma = sigma_o*(1-T/Tc)**n
    return sigma

def BrockBird_scaling(fluid, sigma_o=0.072, To=298.,**params):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given reference temperature to temperature of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    sigma_o : float
        Surface tension at reference temperature (N/m)

    To : float
        Temperature at reference conditions (K)
    """
    Tc = fluid._fluid_dict['Tc']
    Ti = fluid.pore_conitions['temperature']
    Tro = To/Tc
    Tri = Ti/Tc
    sigma_i = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return sigma_i



