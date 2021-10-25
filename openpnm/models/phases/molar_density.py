r"""
"""
import numpy as np


def standard(target, mol_weight='pore.molecular_weight',
             density='pore.density'):
    r"""
    Calculates the molar density from the molecular weight and mass density

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    mol_weight : string
        The dictionary key containing the molecular weight in kg/mol

    density : string
        The dictionary key containing the density in kg/m3

    Returns
    -------
    value : NumPy ndarray
        Array containing molar density values [mol/m3]

    """
    MW = target[mol_weight]
    rho = target[density]
    value = rho/MW
    return value


def ideal_gas(target, pressure='pore.pressure',
              temperature='pore.temperature'):
    r"""
    Uses ideal gas law to calculate the molar density of an ideal gas

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    temperature : string
        The dictionary key containing the density in kg/m3

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    Returns
    -------
    value : NumPy ndarray
        Array containing molar density values [mol/m3]

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


def vanderwaals(target, pressure='pore.pressure',
                temperature='pore.temperature',
                critical_pressure='pore.critical_pressure',
                critical_temperature='pore.critical_temperature'):
    r"""
    Uses Van der Waals equation of state to calculate the density of a real gas

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    critical_pressure : string
        The dictionary key containing the critical pressure values in Pascals
        (Pa)

    critical_temperature : string
        The dictionary key containing the critical temperature values in Kelvin
        (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing molar density values [mol/m3]

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
