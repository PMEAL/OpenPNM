r"""
"""


def fuller(target, MA, MB, vA, vB, temperature='pore.temperature',
           pressure='pore.pressure'):
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first
    principles at conditions of interest

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    MA : float, array_like
        Molecular weight of component A [kg/mol]

    MB : float, array_like
        Molecular weight of component B [kg/mol]

    vA:  float, array_like
        Sum of atomic diffusion volumes for component A

    vB:  float, array_like
        Sum of atomic diffusion volumes for component B

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing gas diffusion coefficient values [m2/s].

    """

    T = target[temperature]
    P = target[pressure]
    MAB = 2*(1.0/MA+1.0/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    return value


def fuller_scaling(target, DABo, To, Po, temperature='pore.temperature',
                   pressure='pore.pressure'):
    r"""
    Uses Fuller model to adjust a diffusion coefficient for gases from
    reference conditions to conditions of interest

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To : float, array_like
        Pressure & temperature at reference conditions, respectively

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing scaled gas diffusion coefficient values [m2/s].

    """
    Ti = target[temperature]
    Pi = target[pressure]
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    return value


def tyn_calus(target, VA, VB, sigma_A, sigma_B, temperature='pore.temperature',
              viscosity='pore.viscosity'):
    r"""
    Uses Tyn_Calus model to estimate diffusion coefficient in a dilute liquid
    solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    VA : float, array_like
        Molar volume of component A at boiling temperature (m3/mol)

    VB : float, array_like
        Molar volume of component B at boiling temperature (m3/mol)

    sigmaA:  float, array_like
        Surface tension of component A at boiling temperature (N/m)

    sigmaB:  float, array_like
        Surface tension of component B at boiling temperature (N/m)

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing liquid diffusion coefficient values [m2/s].

    """
    T = target[temperature]
    mu = target[viscosity]
    A = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T
    B = (sigma_B/sigma_A)**0.15/(mu*1e3)
    value = A*B
    return value


def tyn_calus_scaling(target, DABo, To, mu_o, viscosity='pore.viscosity',
                      temperature='pore.temperature'):
    r"""
    Uses Tyn_Calus model to adjust a diffusion coeffciient for liquids from
    reference conditions to conditions of interest

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    mu_o, To : float, array_like
        Viscosity & temperature at reference conditions, respectively

    pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    Returns
    -------
    value : NumPy ndarray
        Array containing scaled liquid diffusion coefficient values [m2/s].

    """
    Ti = target[temperature]
    mu_i = target[viscosity]
    value = DABo*(Ti/To)*(mu_o/mu_i)
    return value


def chapman_enskog(target, MA, MB, sigma_AB, omega_AB,
                   temperature='pore.temperature',
                   pressure='pore.pressure'):
    r"""
    Calculate gas phase diffusion coefficient using Chapman-Enskog equation.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    MA, MB : scalar
        The molecular mass of species A and B in units of kg/mol
    sigma_AB : scalar
        The collision diameter in units of Angstroms
    omega_AB : scalar
        The collision integral
    temperature : str
        Dictionary key to the phase temperature in K (default = 'pore.temperature')
    pressure : str
        Dictionary key to the phase pressure in Pa (default = 'pore.pressure')
    """
    T = target[temperature]
    P = target[pressure]
    DAB = 1.858e-3*(T**3*(0.001/MA + 0.001/MB))**0.5/(P*101325*sigma_AB**2*omega_AB)
    return DAB
