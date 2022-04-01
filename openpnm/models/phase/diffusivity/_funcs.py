from openpnm.utils import Docorator


docstr = Docorator()

__all__ = [
    "fuller",
    "fuller_scaling",
    "tyn_calus",
    "tyn_calus_scaling",
    "chapman_enskog"
]

@docstr.get_sections(base='models.phase.diffusivity', sections=['Returns'])
@docstr.dedent
def fuller(target, MA, MB, vA, vB, temperature='pore.temperature',
           pressure='pore.pressure'):
    r"""
    Uses Fuller model to estimate diffusion coefficient for gases from first
    principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    MA : float
        Molecular weight of component A [kg/mol]
    MB : float
        Molecular weight of component B [kg/mol]
    vA:  float
        Sum of atomic diffusion volumes for component A
    vB:  float, array_like
        Sum of atomic diffusion volumes for component B
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    diffusivities : ndarray
        A numpy ndarray containing diffusion coefficient values [m2/s].
    """
    T = target[temperature]
    P = target[pressure]
    MAB = 2*(1.0/MA+1.0/MB)**(-1)
    MAB = MAB*1e3
    P = P*1e-5
    value = 0.00143*T**1.75/(P*(MAB**0.5)*(vA**(1./3)+vB**(1./3))**2)*1e-4
    return value


@docstr.dedent
def fuller_scaling(target, DABo, To=298.0, Po=101325,
                   temperature='pore.temperature',
                   pressure='pore.pressure'):
    r"""
    Uses Fuller correlation to adjust a diffusion coefficient for gases from
    reference conditions to conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    DABo : float
        Diffusion coefficient at reference conditions
    Po : float
        Pressure at reference conditions
    To : float
        Temperature at reference conditions
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    Ti = target[temperature]
    Pi = target[pressure]
    value = DABo*(Ti/To)**1.75*(Po/Pi)
    return value


@docstr.dedent
def tyn_calus(target, VA, VB, sigma_A, sigma_B, temperature='pore.temperature',
              viscosity='pore.viscosity'):
    r"""
    Uses Tyn-Calus model to estimate diffusion coefficient in a dilute liquid
    solution of A in B from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    VA : float
        Molar volume of component A at boiling temperature (m3/mol)
    VB : float
        Molar volume of component B at boiling temperature (m3/mol)
    sigmaA : float
        Surface tension of component A at boiling temperature (N/m)
    sigmaB : float
        Surface tension of component B at boiling temperature (N/m)
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    T = target[temperature]
    mu = target[viscosity]
    A = 8.93e-8*(VB*1e6)**0.267/(VA*1e6)**0.433*T
    B = (sigma_B/sigma_A)**0.15/(mu*1e3)
    value = A*B
    return value


@docstr.dedent
def tyn_calus_scaling(target, DABo, To, mu_o,
                      viscosity='pore.viscosity',
                      temperature='pore.temperature'):
    r"""
    Uses Tyn-Calus model to adjust a diffusion coeffcient for liquids from
    reference conditions to conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    DABo : float
        Diffusion coefficient at reference conditions
    mu_o : float
        Viscosity at reference conditions
    To : float
        Temperature at reference conditions
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    Ti = target[temperature]
    mu_i = target[viscosity]
    value = DABo*(Ti/To)*(mu_o/mu_i)
    return value


@docstr.dedent
def chapman_enskog(target, MA, MB, sigma_AB, omega_AB,
                   temperature='pore.temperature',
                   pressure='pore.pressure'):
    r"""
    Calculate gas phase diffusion coefficient using Chapman-Enskog equation.

    Parameters
    ----------
    %(models.target.parameters)s
    MA, MB : scalar
        The molecular mass of species A and B in units of kg/mol
    sigma_AB : scalar
        The collision diameter in units of Angstroms
    omega_AB : scalar
        The collision integral
    %(models.phase.T)s
    %(models.phase.P)s

    Returns
    -------
    %(models.phase.diffusivity.returns)s
    """
    T = target[temperature]
    P = target[pressure]
    DAB = 1.858e-3*(T**3*(0.001/MA + 0.001/MB))**0.5/(P*101325*sigma_AB**2*omega_AB)
    return DAB
