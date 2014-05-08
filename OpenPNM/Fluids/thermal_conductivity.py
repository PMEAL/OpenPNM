r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

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

def Chung(fluid,network,propname,Cv,MW,acentric,viscosity='viscosity',**params):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with low pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    Cv :  float, array_like
        Heat capacity at constant volume (J/(mol.K))
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component

    """
    R = 8.314
    T = fluid.get_pore_data(prop='temperature')
    mu = fluid.get_pore_data(prop=viscosity)
    Tc = fluid.get_pore_data(prop='Tc')
    Tr = T/Tc
    z = 2.0 + 10.5*Tr**2
    beta = 0.7862 - 0.7109*acentric + 1.3168*acentric**2
    alpha = Cv/R -3/2
    s = 1 + alpha*((0.215+0.28288*alpha-1.061*beta+0.26665*z)/(0.6366+beta*z+1.061*alpha*beta))
    value = 3.75*s*(mu)*R/(MW)
    fluid.set_pore_data(prop=propname,data=value)

def Sato(fluid,network,propname,Tb,MW,**params):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids from first principles at conditions of interest

    Parameters
    ----------
    Tb :  float, array_like
        Boiling temperature of the component (K)
    MW : float, array_like
        Molecular weight of the component (kg/mol)

    """
    T = fluid.get_pore_data(prop='temperature')
    Tc = fluid.get_pore_data(prop='Tc')
    Tbr = Tb/Tc
    Tr = T/Tc
    value = (1.11/((MW*1e3)**0.5))*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    fluid.set_pore_data(prop=propname,data=value)
