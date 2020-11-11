r"""
"""
import numpy as np
import scipy as sp


def water(target, temperature='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates surface tension of pure water or seawater at atmospheric
    pressure using Eq. (28) given by Sharqawy et al. Values at
    temperature higher than the normal boiling temperature are calculated at
    the saturation pressure.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    temperature : string
        The dictionary key containing the temperature values.  Temperature must
        be in Kelvin for this emperical equation to work

    salinity : string
        The dictionary key containing the salinity values.  Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    value : NumPy ndarray
        Array containing surface tension of seawater in [N/m]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 313 K; 0 < S < 40 g/kg;
    ACCURACY: 0.2 %

    References
    ----------
    Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    sigma_w = 0.2358*((1-(T/647.096))**1.256)*(1-0.625*(1-(T/647.096)))
    a1 = 2.2637334337E-04
    a2 = 9.4579521377E-03
    a3 = 3.3104954843E-02
    TC = T-273.15
    sigma_sw = sigma_w*(1+(a1*TC+a2)*np.log(1+a3*S))
    value = sigma_sw
    return value


def eotvos(target, k, temperature='pore.temperature',
           critical_temperature='pore.critical_temperature',
           molar_density='pore.molar_density'):
    r"""
    Missing description

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    k : float
        Constant parameter specific to fluid

    temperature : string
        The dictionary key containing the temperature values (K)

    critical_temperature : string
        The dictionary key containing the critical temperature values (K)

    molar_density : string
        The dictionary key containing the molar density values (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing surface tension values [N/m]

    TODO: Needs description, and improve definition of k

    """
    Tc = target[critical_temperature]
    T = target[temperature]
    Vm = 1/target[molar_density]
    value = k*(Tc-T)/(Vm**(2/3))
    return value


def guggenheim_katayama(target, K2, n, temperature='pore.temperature',
                        critical_temperature='pore.critical_temperature',
                        critical_pressure='pore.critical_pressure'):
    r"""
    Missing description

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    K2 : scalar
        Fluid specific constant

    n : scalar
        Fluid specific constant

    temperature : string
        The dictionary key containing the temperature values (K)

    critical_temperature : string
        The dictionary key containing the critical temperature values (K)

    critical_pressure : string
        The dictionary key containing the critical pressure values (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing surface tension values [N/m]

    """
    T = target[temperature]
    Pc = target[critical_pressure]
    Tc = target[critical_temperature]
    sigma_o = K2*Tc**(1/3)*Pc**(2/3)
    value = sigma_o*(1-T/Tc)**n
    return value


def brock_bird_scaling(target, sigma_o, To, temperature='pore.temperature',
                       critical_temperature='pore.critical_temperature'):
    r"""
    Uses Brock_Bird model to adjust surface tension from it's value at a given
    reference temperature to temperature of interest

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    To : float
        Reference temperature (K)

    sigma_o : float
        Surface tension at reference temperature (N/m)

    temperature : string
        The dictionary key containing the temperature values (K)

    critical_temperature : string
        The dictionary key containing the critical temperature values (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing surface tension values scaled to the temperature [N/m]

    """
    Tc = target[critical_temperature]
    Ti = target[temperature]
    Tro = To/Tc
    Tri = Ti/Tc
    value = sigma_o*(1-Tri)**(11/9)/(1-Tro)**(11/9)
    return value
