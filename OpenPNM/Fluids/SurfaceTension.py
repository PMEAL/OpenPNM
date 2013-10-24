
"""
module SurfaceTension
===============================================================================

"""

import scipy as sp
import OpenPNM

def constant(fluid,value,**params):
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
    return value

def na(fluid, **params):
    return

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







