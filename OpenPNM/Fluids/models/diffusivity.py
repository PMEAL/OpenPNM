r"""
===============================================================================
Submodule -- diffusivity
===============================================================================

"""
import scipy as _sp

def fuller(fluid,MA,MB,vA,vB,**kwargs): 
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first 
    principles at conditions of interest

    Parameters
    ----------
    MA : float, array_like
        Molecular weight of component A [kg/mol]
    MB : float, array_like
        Molecular weight of component B [kg/mol]
    vA:  float, array_like
        Sum of atomic diffusion volumes for component A
    vB:  float, array_like
        Sum of atomic diffusion volumes for component B
    """

    T = fluid['pore.temperature']
    P = fluid['pore.pressure']
    MAB = 2*(1/MA+1/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    return value

def fuller_scaling(fluid,DABo,To,Po,**kwargs):
    r"""
    Uses Fuller model to adjust a diffusion coefficient for gases from reference conditions to conditions of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To : float, array_like
        Pressure & temperature at reference conditions, respectively
    """
    Ti = fluid['pore.temperature']
    Pi = fluid['pore.pressure']
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    return value

def tyn_calus(fluid,VA,VB,sigma_A,sigma_B,**kwargs):
    r"""
    Uses Tyn_Calus model to estimate diffusion coefficient in a dilute liquid solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest (K)
    viscosity :  float, array_like
        Viscosity of solvent (Pa.s)
    VA : float, array_like
        Molar volume of component A at boiling temperature (m3/mol)
    VB : float, array_like
        Molar volume of component B at boiling temperature (m3/mol)
    sigmaA:  float, array_like
        Surface tension of component A at boiling temperature (N/m)
    sigmaB:  float, array_like
        Surface tension of component B at boiling temperature (N/m)

    """
    T = fluid['pore.temperature']
    mu = fluid['pore.viscosity']
    value = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T*(sigma_B/sigma_A)**0.15/(mu*1e3)
    return value

def tyn_calus_Scaling(fluid,DABo,To,mu_o,**kwargs):
    r"""
    Uses Tyn_Calus model to adjust a diffusion coeffciient for liquids from reference conditions to conditions of interest

    Parameters
    ----------
    fluid : OpenPNM Fluid Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    mu_o, To : float, array_like
        Viscosity & temperature at reference conditions, respectively
    """
    Ti = fluid['pore.temperature']
    mu_i = fluid['pore.viscosity']
    value = DABo*(Ti/To)*(mu_o/mu_i)
    return value
