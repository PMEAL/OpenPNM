
"""
module diffusivity
===============================================================================

"""
import scipy as sp
import os
propname = os.path.splitext(os.path.basename(__file__))[0]

def constant(fluid,network,value,**params):
    r"""
    Assigns specified constant value
    """
    network.pore_conditions[fluid.name+'_'+propname] = value

def na(fluid,network,**params):
    r"""
    Assigns nonsensical, but numerical value of -1.  This ensurse stability of other methods but introduces the possibility of being misused.
     """
    value = -1
    network.pore_conditions[fluid.name+'_'+propname] = value

def Fuller(fluid, network, MA, MB, vA, vB, **params): #MA=0.03199,MB=0.0291,vA=16.3,vB=19.7,**params):
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest [K]
    P :  float, array_like
        Pressure of interest [Pa]
    MA : float, array_like
        Molecular weight of component A [kg/mol]
    MB : float, array_like
        Molecular weight of component B [kg/mol]
    VA:  float, array_like
        Sum of atomic diffusion volumes for component A
    VB:  float, array_like
        Sum of atomic diffusion volumes for component B
    """
    T = network.pore_conditions[fluid.name+'_'+'temperature']
    P = network.pore_conditions[fluid.name+'_'+'pressure']
    MAB = 2*(1/MA+1/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    network.pore_conditions[fluid.name+'_'+propname] = value

def Fuller_scaling(fluid,network,DABo=2.09e-5,To=298.,Po=101325.,**params):
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
    Ti = network.pore_conditions[fluid.name+'_'+'temperature']
    Pi = network.pore_conditions[fluid.name+'_'+'pressure']
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    network.pore_conditions[fluid.name+'_'+propname] = value

def TynCalus(fluid,network,VA=0.018,VB=0.018,sigma_A=0.07197,sigma_B=0.07197,**params):
    r"""
    Uses Tyn_Calus model to estimate diffusion coefficient in a dilute liquid solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest (K)
    mu :  float, array_like
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
    T = network.pore_conditions[fluid.name+'_'+'temperature']
    mu = network.pore_conditions[fluid.name+'_'+'viscosity']
    value = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T*(sigma_B/sigma_A)**0.15/(mu*1e3)
    network.pore_conditions[fluid.name+'_'+propname] = value

def TynCalus_Scaling(fluid,network,DABo=2.09e-9,To=298.,mu_o=8.90e-4,**params):
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
    Ti = network.pore_conditions[fluid.name+'_'+'temperature']
    mu_i = network.pore_conditions[fluid.name+'_'+'viscosity']
    value = DABo*(Ti/To)*(mu_o/mu_i)
    network.pore_conditions[fluid.name+'_'+propname] = value