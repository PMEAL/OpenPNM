r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

Some text here?

"""

import scipy as _sp

def empirical(fluid,a=[0],**kwargs):
    r"""
    Uses a polynomial fit of empirical data to calculate property
    """
    T = fluid['pore.temperature']
    value = _sp.zeros_like(T)
    for i in range(0,_sp.size(a)):
        value = value + a[i]*(T**i)
    return value

def eotvos(fluid,k,**kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    Tc = fluid['pore.Tc']
    T = fluid['pore.temperature']
    Vm = 1/fluid['pore.molar_density']
    value = k*(Tc-T)/(Vm**(2/3))
    return value

def guggenheim_katayama(fluid,K2,n,**kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    T = fluid['pore.temperature']
    Pc = fluid['pore.Pc']
    Tc = fluid.get_pore_data(prop='Tc')
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    return value

def brock_bird_scaling(fluid,sigma_o,To,**params):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given 
    reference temperature to temperature of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    sigma_o : float
        Surface tension at reference temperature (N/m)

    To : float
        Temperature at reference conditions (K)
    """
    Tc = fluid['pore.Tc']
    Ti = fluid['pore.temperature']
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return value



