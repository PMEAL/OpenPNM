r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

Some text here?

"""

import scipy as sp

def water(phase,**kwargs):
    r"""
    Calculates surface tension of pure water or seawater at atmospheric pressure
    using Eq. (28) given by Sharqawy et. al [1]_. Values at temperature higher 
    than the normal boiling temperature are calculated at the saturation pressure.

    Parameters
    ----------
    T, S: strings
        Property names where phase temperature and salinity are located.
    
    Returns
    -------
    sigma_sw, the surface tension of seawater in [N/m]
    
    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 313 K; 0 < S < 40 g/kg;
    ACCURACY: 0.2 %
    
    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water Treatment, 2010.

    """
    T = phase['pore.temperature']
    try:
        S = phase['pore.salinity']
    except:
        S = 0
    sigma_w = 0.2358*((1-(T/647.096))**1.256)*(1-0.625*(1-(T/647.096)));
    a1 = 2.2637334337E-04; a2 = 9.4579521377E-03; a3 = 3.3104954843E-02
    TC = T-273.15
    sigma_sw = sigma_w*(1+(a1*TC+a2)*sp.log(1+a3*S));    
    value = sigma_sw
    return value
    
def eotvos(phase,
           k,
           pore_molar_density='pore.molar_density',
           **kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    Tc = phase['pore.Tc']
    T = phase['pore.temperature']
    Vm = 1/phase[pore_molar_density]
    value = k*(Tc-T)/(Vm**(2/3))
    return value

def guggenheim_katayama(phase,
                        K2,
                        n,
                        **kwargs):
    r'''
    Documentation for this method is being updated, we are sorry for the inconvenience.
    '''
    T = phase['pore.temperature']
    Pc = phase['pore.Pc']
    Tc = phase['pore.Tc']
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    return value

def brock_bird_scaling(phase,
                       sigma_o,
                       To,
                       **params):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given 
    reference temperature to temperature of interest

    Parameters
    ----------
    phase : OpenPNM Phase Object

    sigma_o : float
        Surface tension at reference temperature (N/m)

    To : float
        Temperature at reference conditions (K)
    """
    Tc = phase['pore.Tc']
    Ti = phase['pore.temperature']
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return value



