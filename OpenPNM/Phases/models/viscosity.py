# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- viscosity
===============================================================================

Models for predicting phase viscosity

"""
import scipy as sp

def poly(phase,a0,a1,a2,a3,**kwargs):
    r"""
    Uses a polynomial best fit correlation for the temperature dependance of 
    viscosity adapted from [1]_ with a range from 80 to 2000 K for air. The model
    can be used also with any gas or liquid if the correlation coefficients 
    are generated using any best fit software (e.g. Excel)
    
    Parameters
    ----------
    a0, a1, a2, a3 : float, array_like
            Coefficients of the viscosity polynomial model 
            (mu = a0 + a1*T + a2*T^2 + a3*T^3)
            where T is the temperature in Kelvin
            
    [1] Thermophysical Properties of Matter Vol. 11: 
    Viscosity,  Y.S. Touloukian, S.C. Saxena, and P. Hestermans
    IFI/Plenun, NY, 1970, ISBN 0-306067020-8 
    
    """
    T = phase['pore.temperature']
    value = a0+a1*T+a2*(T**2)+a3*(T**3)
    return value

def reynolds(phase,uo,b,**kwargs):
    r"""
    Uses exponential model by Reynolds [1]_ for the temperature dependance of 
    shear viscosity
    
    Parameters
    ----------
    u0, b : float, array_like
            Coefficients of the viscosity exponential model (mu = u0*Exp(-b*T)
            where T is the temperature in Kelvin
            
    [1] Reynolds O. (1886). Phil Trans Royal Soc London, v. 177, p.157.
    
    """
    T = phase['pore.temperature']
    value = uo*sp.exp(b*T)
    return value

def chung(phase,MW='molecular_weight',Tc='critical_temperature',Vc='critical_volume',**kwargs):
    r"""
    Uses Chung et al. [2]_ model to estimate viscosity for gases with low pressure 
    (much less than the critical pressure) at conditions of interest

    Parameters
    ----------
    Vc :  float, array_like
        Critical volume of the gas (m3/kmol)
    Tc :  float, array_like
        Critical temperature of the gas (K)
    MW : float, array_like
        Molecular weight of the gas (kg/kmol)
    
    [2] Chung, T.H., Lee, L.L., and Starling, K.E., Applications of Kinetic Gas 
        Theories and Multiparameter Correlation for Prediction of Dilute Gas 
        Viscosity and Thermal Conductivity”, Ind. Eng. Chem. Fundam.23:8, 1984.

    """
    T = phase['pore.temperature']
    MW = phase['pore.'+MW]
    Tc = phase['pore.'+Tc]
    Vc = phase['pore.'+Vc]
    Tr= T/Tc
    Tstar = 1.2593*Tr
    A = 1.161415; B = 0.14874; C = 0.52487; D = 0.77320; E = 2.16178; F = 2.43787
    omega = (A*(Tstar)**(-B)) + C*(sp.exp(-D*Tstar)) + E*(sp.exp(-F*Tstar))
    sigma=0.809*(Vc**(1/3))
    value = 26.69E-9*sp.sqrt(MW*T)/(omega*sigma**2)
    return value
