r"""
===============================================================================
Submodule -- viscosity
===============================================================================

Models for predicting fluid viscosity

"""
import scipy as _sp

def reynolds(fluid,uo,b,**kwargs):
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
    T = fluid['pore.temperature']
    value = uo*_sp.exp(-1*b*T)
    return value

def chung(fluid,MW='pore.MW',Tc='pore.Tc',Vc='pore.Vc',**kwargs):
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
    T = fluid['pore.temperature']
    MW = fluid[MW]
    Tc = fluid[Tc]
    Vc = fluid[Vc]
    Tr= T/Tc
    Tstar = 1.2593*Tr
    A = 1.161415; B = 0.14874; C = 0.52487; D = 0.77320; E = 2.16178; F = 2.43787
    omega = (A*(Tstar)**(-B)) + C*(_sp.exp(-D*Tstar)) + E*(_sp.exp(-F*Tstar))
    value = 26.69E-9*sqrt(MW*T)/(omega*sigma**2)
    return value
