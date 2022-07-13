import numpy as np
from openpnm.utils import Docorator
import chemicals as chem
from chemicals import numba_vectorized


docstr = Docorator()


__all__ = [
    "water_correlation",
    "liquid_pure_antoine",
    "liquid_pure_lee_kesler",
]


def water_correlation(target, temperature='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates vapor pressure of pure water or seawater given by [1] based on
    Raoult's law. The pure water vapor pressure is given by [2]

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : str
        The dictionary key containing the phase salinity values

    Returns
    -------
    %(models.phase.vapor_pressure.returns)s

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 473 K; 0 < S < 240 g/kg;
    ACCURACY: 0.5 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.
    [2] ASHRAE handbook: Fundamentals, ASHRAE; 2005.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
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


def liquid_pure_lee_kesler(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    critical_pressure='param.critical_pressure',
    acentric_factor='param.acentric_factor',
):
    r"""
    """
    T = target[temperature]
    Tc = target[critical_temperature]
    Tr = T/Tc
    Pc = target[critical_pressure]
    omega = target[acentric_factor]
    f0 = 5.92714 - 6.09648/Tr - 1.28862*np.log(Tr) + 0.169347*(Tr**6)
    f1 = 15.2518 - 15.6875/Tr - 13.4721*np.log(Tr) + 0.43577*(Tr**6)
    Pr = np.exp(f0 + omega*f1)
    Pvap = Pr*Pc
    return Pvap


def liquid_pure_antoine(
    target,
    temperature='pore.temperature'
):
    r"""
    """
    T = target[temperature]
    CAS = target.params['CAS']
    Tc = target['param.critical_temperature']
    try:
        coeffs = chem.vapor_pressure.Psat_data_AntoineExtended.loc[CAS]
        _, A, B, C, Tc, to, n, E, F, Tmin, Tmax = coeffs
        PV = numba_vectorized.TRC_Antoine_extended(T, A, B, C, n, E, F)
    except KeyError:
        coeffs = chem.vapor_pressure.Psat_data_AntoinePoling.loc[CAS]
        _, A, B, C, Tmin, Tmax = coeffs
        PV = 10**(A - B/(T + C))
    return PV
