import numpy as np
from openpnm.models.phase import _phasedocs


__all__ = [
    "water_correlation",
    "liquid_pure_antoine",
    "liquid_pure_lk",
]


@_phasedocs
def water_correlation(phase, T='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates vapor pressure of pure water or seawater given by [1] based on
    Raoult's law. The pure water vapor pressure is given by [2]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(salinity)s

    Returns
    -------


    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 473 K; 0 < S < 240 g/kg;
    ACCURACY: 0.5 percent

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.
    [2] ASHRAE handbook: Fundamentals, ASHRAE; 2005.

    """
    T = phase[T]
    if salinity in phase.keys():
        S = phase[salinity]
    else:
        S = 0
    a1 = -5.8002206E+03
    a2 = 1.3914993E+00
    a3 = -4.8640239E-02
    a4 = 4.1764768E-05
    a5 = -1.4452093E-08
    a6 = 6.5459673E+00
    Pv_w = np.exp((a1/T) + a2 + a3*T + a4*T**2 + a5*T**3 + a6*np.log(T))
    Pv_sw = Pv_w/(1+0.57357*(S/(1000-S)))
    value = Pv_sw
    return value


@_phasedocs
def liquid_pure_lk(
    phase,
    T='pore.temperature',
    Tc='param.critical_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    r"""
    Calculate the vapor pressure of a pure liquid using the correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s
    %(Pc)s
    %(omega)s

    Returns
    -------

    References
    ----------
    [1] Lee and Kesler
    """
    T = phase[T]
    Tc = phase[Tc]
    Tr = T/Tc
    Pc = phase[Pc]
    omega = phase[omega]
    f0 = 5.92714 - 6.09648/Tr - 1.28862*np.log(Tr) + 0.169347*(Tr**6)
    f1 = 15.2518 - 15.6875/Tr - 13.4721*np.log(Tr) + 0.43577*(Tr**6)
    Pr = np.exp(f0 + omega*f1)
    Pvap = Pr*Pc
    return Pvap


@_phasedocs
def liquid_pure_antoine(
    phase,
    T='pore.temperature',
):
    r"""
    Calculates the vapor pressure of a pure liquid using Antoine's equation

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s

    Returns
    -------

    References
    ----------
    [1] RPP

    Notes
    -----
    Coefficients are looked up from the ``chemicals`` package. If values for
    "extended" version are not found then the normal 3-coefficient values are
    sought.

    """
    # either antoine or extended antoine using constants from RPP
    CAS = phase.params['CAS']
    T = phase[T]
    from chemicals.vapor_pressure import Psat_data_AntoinePoling
    coeffs = Psat_data_AntoinePoling.loc[CAS]
    _, A, B, C, Tmin, Tmax = coeffs
    Pvap = 10**(A - B/(T + C))
    return Pvap
