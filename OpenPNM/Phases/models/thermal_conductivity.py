# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- thermal_conductance
===============================================================================

"""
import scipy as sp

def water(phase,**kwargs):
    r"""
    Calculates thermal conductivity of pure water or seawater at atmospheric pressure
    using the correlation given in [1]_. Values at temperature higher 
    than the normal boiling temperature are calculated at the saturation pressure.

    Parameters
    ----------
    T, S: strings
        Property names where phase temperature and salinity are located.
    
    Returns
    -------
    k_sw, the thermal conductivity of water/seawater in [W/m.K]
    
    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 3 %
    
    References
    ----------
    [1] D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.

    """
    T = phase['pore.temperature']
    try:
        S = phase['pore.salinity']
    except:
        S = 0
    T68 = 1.00024*T;      # convert from T_90 to T_68" 
    SP = S/1.00472;       # convert from S to S_P"
    k_sw = 0.001*(10**(sp.log10(240+0.0002*SP)+0.434*(2.3-(343.5+0.037*SP)/T68)*((1-T68/(647.3+0.03*SP)))**(1/3)))
    value = k_sw
    return value

def air(phase,**kwargs):
    r"""
    Calculates thermal conductivity of air at atmospheric pressure
    using a polynominal correlation that fits the data given in [3]_.
    
    Parameters
    ----------
    T: strings
        Property names where phase temperature is located.  
            
    Returns
    -------
    k_air, the thermal conductivity of air in [W/m.K]
    
    Notes
    -----
    T must be in K. 
    VALIDITY: 273 < T < 1023 K
    ACCURACY: 0.2 %
    
    References
    ----------
    [3] E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity 
    Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of Thermophysics, 
    Vol. 25, No. 1, January 2004, pp. 21-69

    """
    T = phase['pore.temperature']
    a0=0.00422791; a1=0.0000789606; a2=-1.56383E-08
    k_air = a0 + a1*T + a2*(T**2)
    value = k_air
    return value

def chung(phase,
          Cv,
          MW,
          acentric,
          pore_viscosity='pore.viscosity',
          **kwargs):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with 
    low pressure(<10 bar) from first principles at conditions of interest

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
    T = phase['pore.temperature']
    mu = phase[pore_viscosity]
    Tc = phase['pore.Tc']
    Tr = T/Tc
    z = 2.0 + 10.5*Tr**2
    beta = 0.7862 - 0.7109*acentric + 1.3168*acentric**2
    alpha = Cv/R -3/2
    s = 1 + alpha*((0.215+0.28288*alpha-1.061*beta+0.26665*z)/(0.6366+beta*z+1.061*alpha*beta))
    value = 3.75*s*(mu)*R/(MW)
    return value

def sato(phase,Tb,MW,**params):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids 
    from first principles at conditions of interest

    Parameters
    ----------
    Tb :  float, array_like
        Boiling temperature of the component (K)
    MW : float, array_like
        Molecular weight of the component (kg/mol)

    """
    T = phase['pore.temperature']
    Tc = phase['pore.Tc']
    Tbr = Tb/Tc
    Tr = T/Tc
    value = (1.11/((MW*1e3)**0.5))*(3+20*(1-Tr)**(2/3))/(3+20*(1-Tbr)**(2/3))
    return value
