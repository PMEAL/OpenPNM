
"""
module Viscosity
===============================================================================

"""
import OpenPNM
import scipy as sp

def constant(fluid,value=0.001,**params):
    return value

def na(fluid,**params):
    return 'n/a'

def Reynolds(fluid,uo=0.001,b=0.1,**params):
    r"""
    Uses Reynolds method for the temperature dependance of shear viscosity 

    Parameters
    ----------
    u0, b :  float, array_like
            Coefficients of Reynolds method
        
    """     
    T = fluid.pore_conditions['temperature']
    mu = uo*sp.exp(-1*b*T)
    return mu
    
def Chung(fluid,Tc=132.65,Vc=92.35e-6,MW=0.0291,acentric=0,kappa=0,dipole=0,**params):
    r"""
    Uses Chung et al. model to estimate viscosity for gases with low pressure(not near the critical pressure) from first principles at conditions of interest 

    Parameters
    ----------
    Tc :  float, array_like
        Critical Temperature of the component (K)
    Vc :  float, array_like
        Critical volume of the component (m3/mol)
    MW : float, array_like
        Molecular weight of the component (kg/mol)
    acentric : float, array_like
        Acentric factor of the component 
    kappa : float, array_like
        Special correction for highly polar substances
    dipole :float, array_like
        Dipole moment (C.m)

    """                         
    T = fluid.pore_conditions['temperature']
    Tr= T/Tc
    Tstar = 1.2593*Tr
    A = 1.161415
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    sigma = (A*(Tstar)**(-B)) + C*(sp.exp(-D*Tstar)) + E*(sp.exp(-F*Tstar))
    dipole_r = 131.3*(dipole*2.997e29)/((Vc*1e6)*Tc)**0.5
    f = 1-0.2756*acentric + 0.059035*(dipole_r**4) + kappa
    mu = 40.785*f*(MW*1e3*T)**0.5/((Vc*1e6)**(2/3)*sigma)*1e-7
    
    return mu
    
   
    
    