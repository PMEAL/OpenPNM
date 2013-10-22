
"""
module Diffusivity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(network,value,**params):
    return value

def na(network,**params):
    return 'n/a'

def Fuller(network,MA=31.99,MB=28.01,vA=16.6,vB=17.9,**params):
    r"""
    Uses the Fuller model to estimate diffusion coefficient for gases from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest (kelvin)
    P :  float, array_like
        Pressure of interest (bar)
    MA : float, array_like
        Molecular weight of component A (g/mole)
    MB : float, array_like
        Molecular weight of component B (g/mole)
    VA:  float, array_like
        Sum of the diffusion volumes for component A
    VB:  float, array_like
        Sum of the diffusion volumes for component B

    Notes
    -----
    The Fuller equation is:

    .. math::


    """
    T = network.pore_conditions['temperature']
    P = network.pore_conditions['pressure']
    MAB = 2*(1/MA+1/MB)**(-1)
    DAB = 0.00143*T**1.75/(P*MAB**0.5*(vA**(1/3)+vB**(1/3))**2)*1e-4
    return DAB


def Fuller_scaling(network,DABo=2.09e-5,To=298.,Po=101325.,Ti=298.,Pi=101325.,**params):
    r"""
    Uses the Fuller model to adjust a diffusion coeffciient for gases from reference conditions to conditions of interest

    Parameters
    ----------
    network : OpenPNM Network Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To, Pi & Ti : float, array_like
        Pressure & temperature at reference conditions, pressure & temperature at conditions of interest, respectively
    """
    DAB = DABo*(Ti/To)**1.75*(Po/Pi)
    return DAB

def TynCalus(network,T=298.,etta=0.890,VA=16.6,VB=17.9,sigmaB=71.97,sigmaA=22.27,**params):
    r"""
    Uses the Tyn Calus model to estimate diffusion coefficient in a dilute solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    T :  float, array_like
        Temperature of interest (kelvin)
    etta :  float, array_like
        Viscositty of solvent (cP)
    VA : float, array_like
        Molar volume of component A (g/s^2)
    VB : float, array_like
        Molar volume of component B at boiling temperature (g/s^2)
    sigmaA:  float, array_like
        Surface tension of component A at boiling temperature (mN/m)
    sigmaB:  float, array_like
        Surface tension of component B at boiling temperature (mN/m)
    Notes
    -----
    The Tyn Calus equation is:

    .. math::


    """
    T = network.pore_conditions['temperature']
    DAB = 8.93e-8*(VB**0.267/VA**0.433)*T*(sigmaB/sigmaA)**0.15/etta
    return DAB

def TynCalus_Scaling(network,DABo=2.09e-5,To=298.,ettao=0.890,Ti=298.,ettai=0.890,**params):
    r"""
    Uses the Tyn Calus model to adjust a diffusion coeffciient for liquids from reference conditions to conditions of interest

    Parameters
    ----------
    network : OpenPNM Network Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    ettao, To, ettai & Ti : float, array_like
        Viscosity & temperature at reference conditions, viscosity & temperature at conditions of interest, respectively
    """

    DAB = DABo*(Ti/To)*(ettao/ettai)
    return DAB