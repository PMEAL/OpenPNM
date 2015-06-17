r"""
===============================================================================
Submodule -- diffusivity
===============================================================================

"""
import scipy as _sp


def fuller(phase, MA, MB, vA, vB,
           pore_temperature='pore.temperature',
           pore_pressure='pore.pressure',
           **kwargs):
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first
    principles at conditions of interest

    Parameters
    ----------
    MA : float, array_like
        Molecular weight of component A [kg/mol]

    MB : float, array_like
        Molecular weight of component B [kg/mol]

    vA:  float, array_like
        Sum of atomic diffusion volumes for component A

    vB:  float, array_like
        Sum of atomic diffusion volumes for component B

    pore_pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    pore_temperature : string
        The dictionary key containing the temperature values in Kelvin (K)
    """

    T = phase[pore_temperature]
    P = phase[pore_pressure]
    MAB = 2*(1.0/MA+1.0/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    return value


def fuller_scaling(phase, DABo, To, Po,
                   pore_temperature='pore.temperature',
                   pore_pressure='pore.pressure',
                   **kwargs):
    r"""
    Uses Fuller model to adjust a diffusion coefficient for gases from
    reference conditions to conditions of interest

    Parameters
    ----------
    phase : OpenPNM Phase Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    Po, To : float, array_like
        Pressure & temperature at reference conditions, respectively

    pore_pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    pore_temperature : string
        The dictionary key containing the temperature values in Kelvin (K)
    """
    Ti = phase[pore_temperature]
    Pi = phase[pore_pressure]
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    return value


def tyn_calus(phase, VA, VB, sigma_A, sigma_B,
              pore_temperature='pore.temperature',
              pore_viscosity='pore.viscosity',
              **kwargs):
    r"""
    Uses Tyn_Calus model to estimate diffusion coefficient in a dilute liquid
    solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    VA : float, array_like
        Molar volume of component A at boiling temperature (m3/mol)

    VB : float, array_like
        Molar volume of component B at boiling temperature (m3/mol)

    sigmaA:  float, array_like
        Surface tension of component A at boiling temperature (N/m)

    sigmaB:  float, array_like
        Surface tension of component B at boiling temperature (N/m)

    pore_pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    pore_temperature : string
        The dictionary key containing the temperature values in Kelvin (K)

    """
    T = phase['pore.temperature']
    mu = phase['pore.viscosity']
    A = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T
    B = (sigma_B/sigma_A)**0.15/(mu*1e3)
    value = A*B
    return value


def tyn_calus_scaling(phase, DABo, To, mu_o,
                      pore_temperature='pore.temperature',
                      pore_viscosity='pore.viscosity',
                      **kwargs):
    r"""
    Uses Tyn_Calus model to adjust a diffusion coeffciient for liquids from
    reference conditions to conditions of interest

    Parameters
    ----------
    phase : OpenPNM Phase Object

    DABo : float, array_like
        Diffusion coefficient at reference conditions

    mu_o, To : float, array_like
        Viscosity & temperature at reference conditions, respectively

    pore_pressure : string
        The dictionary key containing the pressure values in Pascals (Pa)

    pore_temperature : string
        The dictionary key containing the temperature values in Kelvin (K)
    """
    Ti = phase[pore_temperature]
    mu_i = phase[pore_viscosity]
    value = DABo*(Ti/To)*(mu_o/mu_i)
    return value
