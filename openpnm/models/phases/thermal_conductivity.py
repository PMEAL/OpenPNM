import numpy as np
from openpnm.utils import Docorator
from chemicals import numba_vectorized


docstr = Docorator()


@docstr.dedent
def water(target, temperature="pore.temperature", salinity="pore.salinity"):
    r"""
    Calculates thermal conductivity of pure water or seawater at atmospheric
    pressure using the correlation given by Jamieson and Tudhope. Values at
    temperature higher  the normal boiling temperature are calculated at the
    saturation pressure.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : string
        The dictionary key containing the salinity values.  Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity of water/seawater in [W/m.K]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
    thousand)

    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 3 %

    References
    ----------
    D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    T68 = 1.00024 * T  # convert from T_90 to T_68
    SP = S / 1.00472  # convert from S to S_P
    k_sw = 0.001 * (
        10
        ** (
            np.log10(240 + 0.0002 * SP)
            + 0.434
            * (2.3 - (343.5 + 0.037 * SP) / T68)
            * ((1 - T68 / (647.3 + 0.03 * SP))) ** (1 / 3)
        )
    )
    value = k_sw
    return value


@docstr.dedent
def chung(
    target,
    Cv="pore.heat_capacity",
    acentric_factor="pore.acentric_factor",
    mol_weight="pore.molecular_weight",
    viscosity="pore.viscosity",
    temperature="pore.temperature",
    critical_temperature="pore.critical_temperature",
):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with
    low pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    acentric_factor : string
        Dictionary key containing the acentric factor of the component
    Cv : string
        Dictionary key containing the heat capacity at constant volume
        (J/(mol.K))
    mol_weight : string
        Dictionary key containing the molecular weight of the component
        (kg/mol)
    viscosity : string
        The dictionary key containing the viscosity values (Pa.s)
    %(models.phase.T)s
    critical_temperatre: string
        The dictionary key containing the critical temperature values (K)

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity values in [W/m.K]

    """
    Cv = target[Cv]
    acentric = target[acentric_factor]
    MW = target[mol_weight]
    R = 8.314
    T = target[temperature]
    mu = target[viscosity]
    Tc = target[critical_temperature]
    Tr = T / Tc
    z = 2.0 + 10.5 * Tr ** 2
    beta = 0.7862 - 0.7109 * acentric + 1.3168 * acentric ** 2
    alpha = Cv / R - 3 / 2
    s = 1 + alpha * (
        (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
        / (0.6366 + beta * z + 1.061 * alpha * beta)
    )
    value = 3.75 * s * (mu) * R / (MW)
    return value


def sato(
    target,
    mol_weight="pore.molecular_weight",
    boiling_temperature="pore.boiling_point",
    temperature="pore.temperature",
    critical_temperature="pore.critical_temperature",
):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids
    from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    boiling_temperature :  string
        Dictionary key containing the toiling temperature of the component (K)
    mol_weight : string
        Dictionary key containing the molecular weight of the component
        (kg/mol)
    %(models.phase.T)s
    critical_temperature : string
        The dictionary key containing the critical temperature values (K)

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity values in [W/m.K]

    """
    T = target[temperature]
    Tc = target[critical_temperature]
    MW = target[mol_weight]
    Tbr = target[boiling_temperature] / Tc
    Tr = T / Tc
    value = (
        (1.11 / ((MW * 1e3) ** 0.5))
        * (3 + 20 * (1 - Tr) ** (2 / 3))
        / (3 + 20 * (1 - Tbr) ** (2 / 3))
    )
    return value


def liquid_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target['param.molecular_weight']
    Tb = target['param.boiling_temperature']
    Pc = target['param.critical_pressure']
    omega = target['param.acentric_factor']
    kL = numba_vectorized.Gharagheizi_liquid(T, MW, Tb, Pc, omega)
    return kL


def gas_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    MW = target['param.molecular_weight']
    Tb = target['param.boiling_temperature']
    Pc = target['param.critical_pressure']
    omega = target['param.acentric_factor']
    kG = numba_vectorized.Gharagheizi_gas(T, MW, Tb, Pc, omega)
    return kG
