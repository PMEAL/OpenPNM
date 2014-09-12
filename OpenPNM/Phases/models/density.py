r"""
===============================================================================
Submodule -- density
===============================================================================

"""
import scipy as _sp

def ideal_gas(phase,**kwargs):
    r"""
    Uses ideal gas law to calculate the mass density of an ideal gas
 
    Parameters
    ----------
    P, T, MW: float, array_like
        P pressure of the gas in [Pa]
        T temperature of the gas in [K]
        MW molecular weight of the gas in [kg/kmole]
            
    Returns
    -------
    rho, the density in [mol/m3]
    
    """
   
    P = phase['pore.pressure']
    T = phase['pore.temperature']
    MW = phase['pore.molecular_weight']
    R = 8.31447
    value = P/(R*T)*MW
    return value

def water(phase,**kwargs):
    r"""
    Calculates density of pure water or seawater at atmospheric pressure
    using Eq. (8) given by Sharqawy et. al [1]_. Values at temperature higher 
    than the normal boiling temperature are calculated at the saturation pressure.

    Parameters
    ----------
    T, S: strings
        Property names where phase temperature and salinity are located.
    
    Returns
    -------
    rho_sw, the density of water/seawater in [kg/m3]
    
    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 0.1 %
    
    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water Treatment, 2010.

    """
    T = phase['pore.temperature']
    try:
        S = phase['pore.salinity']
    except:
        S = 0
    a1=9.9992293295E+02; a2=2.0341179217E-02; a3=-6.1624591598E-03; a4=2.2614664708E-05; a5=-4.6570659168E-08
    b1=8.0200240891E-01; b2=-2.0005183488E-03; b3=1.6771024982E-05; b4=-3.0600536746E-08; b5=-1.6132224742E-11
    TC = T-273.15
    rho_w = a1 + a2*TC + a3*TC**2 + a4*TC**3 + a5*TC**4;
    D_rho = b1*S + b2*S*TC + b3*S*(TC**2) + b4*S*(TC**3) + b5*(S**2)*(TC**2);
    rho_sw = rho_w + D_rho
    value = rho_sw
    return value
