
"""
module thermal_conductance
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(fluid,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def na(fluid,network,**params):
    value = -1
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def Chung(fluid,network,Tc=132.64,Cv=1000,MW=0.0291,acentric=0.03,**params):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with low pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    Tc :  float, array_like
        Critical Temperature of the component (K)
    Cv :  float, array_like
        Heat capacity at constant volume (J/(mol.K))
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component

    """
    R = 8.314
    T = network.get_pore_data(phase=fluid,prop='temperature')
    mu = network.get_pore_data(phase=fluid,prop='viscosity')
    Tr = T/Tc
    z = 2.0 + 10.5*Tr**2
    beta = 0.7862 - 0.7109*acentric + 1.3168*acentric**2
    alpha = Cv/R -3/2
    s = 1 + alpha*((0.215+0.28288*alpha-1.061*beta+0.26665*z)/(0.6366+beta*z+1.061*alpha*beta))
    value = 3.75*s*(mu)*R/(MW)
    network.set_pore_data(phase=fluid,prop=propname,data=value)

def Sato(fluid,network,Tc=647.096,Tb=373.15,MW=0.0181,**params):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids from first principles at conditions of interest

    Parameters
    ----------
    Tc :  float, array_like
        Critical Temperature of interest (K)
    Tb :  float, array_like
        Boiling temperature of the component (K)
    MW : float, array_like
        Molecular weight of the component (kg/mol)

    """
    T = network.get_pore_data(phase=fluid,prop='temperature')
    Tbr = Tb/Tc
    Tr = T/Tc
    value = (1.11/((MW*1e3)**0.5))*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    network.set_pore_data(phase=fluid,prop=propname,data=value)
