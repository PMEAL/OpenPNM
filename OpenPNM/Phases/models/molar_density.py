r"""
===============================================================================
Submodule -- molar_density
===============================================================================

"""
import scipy as sp


def standard(phase,
             pore_MW='pore.molecular_weight',
             pore_density='pore.density',
             **kwargs):
    r"""
    Calculates the molar density from the molecular weight and mass density

    Parameters
    ----------
    phase : OpenPNM Phase Object
        The Phase object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides access
        to other necessary thermofluid properties.

    pore_MW : string
        The dictionary key containing the molecular weight in kg/mol

    pore_temperature : string
        The dictionary key containing the density in kg/m3
    """
    MW = phase[pore_MW]
    rho = phase[pore_density]
    value = rho/MW
    return value


def ideal_gas(phase,
              pore_pressure='pore.pressure',
              pore_temperature='pore.temperature',
              **kwargs):
    r"""
    Uses ideal gas law to calculate the molar density of an ideal gas

    Parameters
    ----------
    phase : OpenPNM Phase Object
        The Phase object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides access
        to other necessary thermofluid properties.

    pore_temperature : string
        The dictionary key containing the density in kg/m3

    pore_pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    Returns
    -------
    rho, the density in [mol/m3]

    Notes
    -----
    This method uses the SI value for the ideal gas constant, hence the need to
    provide the temperature and pressure in SI.  In general, OpenPNM use SI
    throughout for consistency.

    """
    R = 8.31447
    P = phase[pore_pressure]
    T = phase[pore_temperature]
    value = P/(R*T)
    return value


def vanderwaals(phase,
                pore_P='pore.pressure',
                pore_T='pore.temperature',
                pore_Pc='pore.critical_pressure',
                pore_Tc='pore.critical_temperature',
                **kwargs):
    r"""
    Uses Van der Waals equation of state to calculate the density of a real gas

    Parameters
    ----------
    pore_P : string
        The dictionary key containing the pressure values in Pascals (Pa)

    pore_T : string
        The dictionary key containing the temperature values in Kelvin (K)

    pore_Pc : string
        The dictionary key containing the critical pressure values in Pascals
        (Pa)

    pore_Tc : string
        The dictionary key containing the critical temperature values in Kelvin
        (K)

    Returns
    -------
    rho, the density in [mol/m3]

    Notes
    -----
    This equation and its constant coefficients are taken [1]_ which uses the
    cgs units system. All input parameters are expected in SI, then converted
    in the method.

    """

    P = phase[pore_P]/100000
    T = phase[pore_T]
    Pc = phase[pore_Pc]/100000  # convert to bars
    Tc = phase[pore_Tc]
    R = 83.1447
    a = 27*(R**2)*(Tc**2)/(64*Pc)
    b = R*Tc/(8*Pc)
    a1 = -1/b
    a2 = (R*T+b*P)/(a*b)
    a3 = -P/(a*b)
    a0 = sp.ones(sp.shape(a1))
    coeffs = sp.vstack((a0, a1, a2, a3)).T
    density = sp.array([sp.roots(C) for C in coeffs])
    value = sp.real(density[:, 2])*1e6  # Convert it to mol/m3
    return value
