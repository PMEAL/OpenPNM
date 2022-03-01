import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.get_sections(base='models.phase.vapor_pressure',
                     sections=['Returns'])
@docstr.dedent
def antoine(target, A, B, C, temperature='pore.temperature'):
    r"""
    Uses Antoine equation to estimate vapor pressure of a pure component

    Parameters
    ----------
    %(models.target.parameters)s
    A, B, C :  scalars
        Antoine vapor pressure coefficients for pure compounds. Since virtually
        all Antoine coefficients are reported for units of mmHg and C for
        historical reasons, this method assumes these A, B and C values are for
        mmHg and C, but converts all properties internally to return Pascals.
    %(models.phase.T)s

    Returns
    -------
    value : ndarray
        Array containing vapor pressure values [Pa]

    References
    ----------
    ::

        Antoine, C. (1888), Vapor Pressure: a new relationship between pressure
        and temperature, Comptes Rendus des Séances de l'Académie des Sciences
        (in French) 107: 681–684, 778–780, 836–837

    """
    T = target[temperature] - 273.15
    value = (10**(A-B/(C+T)))/760*101325
    return value


def water(target, temperature='pore.temperature', salinity='pore.salinity'):
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
