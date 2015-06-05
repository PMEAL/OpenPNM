# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- vapor_pressure
===============================================================================

Methods for predicing the vapor pressure of pure species

"""
import scipy as sp


def antoine(phase, A, B, C, pore_temperature='pore.temperature', **kwargs):
    r"""
    Uses Antoine equation [1]_ to estimate vapor pressure of a pure component

    Parameters
    ----------
    A, B, C :  float, array_like
            Antoine vapor pressure coefficients for pure compounds. Note that
            these coefficients should be converted such that the Temperature
            is in [K] and the vapor pressure is in [Pa].

    pore_temperature : string
        The dictionary key containing the phase temperature values

    [1] Antoine, C. (1888), Vapor Pressure: a new relationship between pressure
        and temperature, Comptes Rendus des Séances de l'Académie des Sciences
        (in French) 107: 681–684, 778–780, 836–837

    """
    T = phase[pore_temperature]
    value = (10**(A-B/(C+T)))
    return value


def water(phase,
          pore_temperature='pore.temperature',
          pore_salinity = 'pore.salinity',
          **kwargs):
    r"""
    Calculates vapor pressure of pure water or seawater given by [1]_ based on
    Raoult's law. The pure water vapor pressure is given by [2]_

    Parameters
    ----------
    pore_temperature : strings
        The dictionary key containing the phase temperature values

    pore_salinity : strings
        The dictionary key containing the phase salinity values

    Returns
    -------
    The vapor pressure of water/seawater in [Pa]

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 473 K; 0 < S < 240 g/kg;
    ACCURACY: 0.5 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water
        Treatment, 2010.
    [2] ASHRAE handbook: Fundamentals, ASHRAE; 2005.

    """
    T = phase[pore_temperature]
    try:
        S = phase[pore_salinity]
    except:
        S = 0
    a1 = -5.8002206E+03
    a2 = 1.3914993E+00
    a3 = -4.8640239E-02
    a4 = 4.1764768E-05
    a5 = -1.4452093E-08
    a6 = 6.5459673E+00
    Pv_w = sp.exp((a1/T) + a2 + a3*T + a4*T**2 + a5*T**3 + a6*sp.log(T))
    Pv_sw = Pv_w/(1+0.57357*(S/(1000-S)))
    value = Pv_sw
    return value
