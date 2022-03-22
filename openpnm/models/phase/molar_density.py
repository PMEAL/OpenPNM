import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
def standard(target, mol_weight='pore.molecular_weight',
             density='pore.density'):
    r"""
    Calculates the molar density from the molecular weight and mass density

    Parameters
    ----------
    %(models.target.parameters)s
    mol_weight : str
        The dictionary key containing the molecular weight in kg/mol
    density : str
        The dictionary key containing the density in kg/m3

    Returns
    -------
    value : ndarray
        A numpy ndrray containing molar density values [mol/m3]

    """
    MW = target[mol_weight]
    rho = target[density]
    value = rho/MW
    return value


@docstr.dedent
def ideal_gas(target, pressure='pore.pressure',
              temperature='pore.temperature'):
    r"""
    Uses ideal gas law to calculate the molar density of an ideal gas

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    value : ndarray
        A numpy ndarray containing molar density values [mol/m3]

    Notes
    -----
    This method uses the SI value for the ideal gas constant, hence the need to
    provide the temperature and pressure in SI.  In general, OpenPNM use SI
    throughout for consistency.

    """
    R = 8.31447
    P = target[pressure]
    T = target[temperature]
    value = P/(R*T)
    return value


@docstr.dedent
def vanderwaals(target, pressure='pore.pressure',
                temperature='pore.temperature',
                critical_pressure='pore.critical_pressure',
                critical_temperature='pore.critical_temperature'):
    r"""
    Uses Van der Waals equation of state to calculate the density of a real gas

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    %(models.phase.P)s
    critical_pressure : str
        The dictionary key containing the critical pressure values in Pascals
        (Pa)
    critical_temperature : str
        The dictionary key containing the critical temperature values in Kelvin
        (K)

    Returns
    -------
    value : ndarray
        A numpy ndarray containing molar density values [mol/m3]

    """

    P = target[pressure]/100000
    T = target[temperature]
    Pc = target[critical_pressure]/100000  # convert to bars
    Tc = target[critical_temperature]
    R = 83.1447
    a = 27*(R**2)*(Tc**2)/(64*Pc)
    b = R*Tc/(8*Pc)
    a1 = -1/b
    a2 = (R*T+b*P)/(a*b)
    a3 = -P/(a*b)
    a0 = np.ones(np.shape(a1))
    coeffs = np.vstack((a0, a1, a2, a3)).T
    density = np.array([np.roots(C) for C in coeffs])
    value = np.real(density[:, 2])*1e6  # Convert it to mol/m3
    return value
