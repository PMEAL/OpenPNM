
"""
module SurfaceTension
===============================================================================

"""

import scipy as sp
import OpenPNM

def set_as(fluid1, fluid2, sigma=0.072):
    r"""
    Set the surface tension of fluid1 relative to fluid2

    Parameters
    ----------
    sigma : array_like
        The numerical value of surface tension

    fluid1 and fluid2 : The fluid pair of interest

    Notes
    -----
    This method sets the surface tension for both fluids simultaneously
    """
    sigma = sp.array(sigma,ndmin=1)
    if 'surface_tension' not in fluid1.keys():
        fluid1.update({'surface_tension': {}})
    if 'surface_tension' not in fluid2.keys():
        fluid2.update({'surface_tension': {}})
    fluid1['surface_tension'].update({fluid2['name']: sigma})
    fluid2['surface_tension'].update({fluid1['name']: sigma})

def Eotvos(fluid1, fluid2, T):
    r"""
    """
    k = 2.1e-7
    Tc1 = fluid1['critical_temperature']
    Vm = 1/fluid1['molar_volume']
    sigma = k*(Tc1-T)/(Vm**(2/3))
    OpenPNM.Fluids.SurfaceTension.set_as(fluid1,fluid2,sigma)

def GuggenheimKatayama(fluid1, fluid2, T, Tc1, Pc1, n=1.222):
    r"""
    """
    K2 = 1
    sigma_o = K2*Tc1**(1/3)*Pc1**(2/3)
    sigma = sigma_o*(1-T/Tc1)**n
    OpenPNM.Fluids.SurfaceTension.set_as(fluid1,fluid2,sigma)

def BrockBird_scaling(fluid,sigma_o=0.072,To=298.,Tc=647.3,**params):
    r"""
    Uses Brock_Bird model to adjust surface tension from reference temperature to temperature of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object
    
    Tc : float, array_like
        Critical temperature of the component (K)
    sigma_o : float, array_like
        Surface tension at reference temperature (N/m)
    To : float, array_like
        Temperature at reference conditions (K)
    """
    Ti = fluid.pore_conitions['temperature']   
    Tro = To/Tc
    Tri = Ti/Tc   
    sigma_i = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return sigma_i



