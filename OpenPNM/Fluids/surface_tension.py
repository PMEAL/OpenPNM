r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

Some text here?

"""

import scipy as sp

def constant(fluid,network,propname,value,**params):
    r"""
    Assigns specified constant value
    """
    fluid.set_pore_data(prop=propname,data=value)

def na(fluid,network,propname,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  
    This ensurse stability of other methods 
    but introduces the possibility of being misused.
    """
    value = -1
    fluid.set_pore_data(prop=propname,data=value)

def empirical(fluid,network,propname,a=[0],**params):
    r"""
    Uses a polynomial fit of empirical data to calculate property
    """
    T = fluid.get_pore_data(prop='temperature')
    value = sp.zeros_like(T)
    for i in range(0,sp.size(a)):
        value = value + a[i]*(T**i)
    fluid.set_pore_data(prop=propname,data=value)

def Eotvos(fluid,network,propname,k,molar_density='molar_density', **params):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    Tc = fluid.get_pore_data(prop='Tc')
    T = fluid.get_pore_data(prop='temperature')
    Vm = 1/fluid.get_pore_data(prop=molar_density)
    value = k*(Tc-T)/(Vm**(2/3))
    fluid.set_pore_data(prop=propname,data=value)

def GuggenheimKatayama(fluid,network,propname,K2,n,**params):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    T = fluid.get_pore_data(prop='temperature')
    Pc = fluid.get_pore_data(prop='Pc')
    Tc = fluid.get_pore_data(prop='Tc')
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    fluid.set_pore_data(prop=propname,data=value)

def BrockBird_scaling(fluid,network,propname,sigma_o,To,**params):
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
    Tc = fluid.get_pore_data(prop='Tc')
    Ti = fluid.get_pore_data(prop='temperature')
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    fluid.set_pore_data(prop=propname,data=value)



