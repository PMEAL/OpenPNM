import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    'water_correlation',
    'gas_pure',
    'liquid_pure',
    'gas_mixture',
    'liquid_mixture',
]


@docstr.dedent
def water_correlation(
    target,
    T="pore.temperature",
    salinity="pore.salinity"
):
    r"""
    Calculates thermal conductivity of pure water or seawater at atmospheric
    pressure using the correlation given by Jamieson and Tudhope. Values at
    temperature higher  the normal boiling temperature are calculated at the
    saturation pressure.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : str
        The dictionary key containing the salinity values.  Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity of water/seawater in [W/m.K]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
    thousand). The correlation is valid for 273 < T < 453 K and
    0 < S < 160 g/kg within 3% accuracy.

    References
    ----------
    D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.

    """
    T = target[T]
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


def gas_pure(
    target,
    T='pore.temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    # gharagheizi method
    T = target[T]
    MW = target[MW]
    Tb = target[Tb]
    # The following correction suggested by chemicals package author
    Pc = target[Pc]/10000
    omega = target[omega]
    B = (T + (2.0*omega + 2.0*T - 2.0*T*(2.0*omega + 3.2825)/Tb + 3.2825)
         / (2.0*omega + T - T*(2.0*omega + 3.2825)/Tb + 3.2825)
         - T*(2.0*omega + 3.2825)/Tb)
    A = (2*omega + T - (2*omega + 3.2825)*T/Tb + 3.2825)/(0.1*MW*Pc*T) \
        * (3.9752*omega + 0.1*Pc + 1.9876*B + 6.5243)**2
    k = 7.9505e-4 + 3.989e-5*T - 5.419e-5*MW + 3.989e-5*A
    return k


def liquid_pure(
    target,
    T='pore.temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    # gharagheizi metho
    T = target[T]
    MW = target[MW]
    Tb = target[Tb]
    Pc = target[Pc]/100000
    omega = target[omega]
    B = 16.0407*MW + 2.0*Tb - 27.9074
    A = 3.8588*(MW**8)*(1.0045*B + 6.5152*MW - 8.9756)
    k = (1e-4)*(10*omega + 2*Pc - 2*T + 4 + 1.908*(Tb + 1.009*(B**2)/(MW**2))
                + 3.9287*(MW**4)/(B**4) + A/(B**8))
    return k


def liquid_pure_sato_riedel(
    target,
    T="pore.temperature",
    Tc='param.critical_temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
):
    r"""
    Uses Sato et al. model to estimate thermal conductivity for pure liquids
    from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity values in [W/m.K]

    """
    T = target[T]
    Tc = target[Tc]
    MW = target[MW]
    Tbr = target[Tb]/Tc
    Tr = T / Tc
    value = ((1.1053 / ((MW) ** 0.5)) * (3 + 20 * (1 - Tr) ** (2 / 3))
             / (3 + 20 * (1 - Tbr) ** (2 / 3)))
    return value


def liquid_mixture_DIPPR9I(
    target,
    rhos='pore.density.*',
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    raise NotImplementedError("This function is not ready yet")
    from chemicals import rho_to_Vm
    xs = target['pore.mole_fraction']
    kLs = target.get_comp_vals(ks)
    # kL = numba_vectorized.DIPPR9I(xs, kLs)  # Another one that doesn't work
    Vm = [rho_to_Vm(c[rhos], c[MWs])
          for c in target.components.keys()]
    denom = np.sum([xs[i]*Vm[i] for i in range(len(xs))], axis=0)
    phis = np.array([xs[i]*Vm[i] for i in range(len(xs))])/denom
    kij = 2/np.sum([1/kLs[i] for i in range(len(xs))], axis=0)
    kmix = np.zeros_like(xs[0])
    N = len(xs)
    for i in range(N):
        for j in range(N):
            kmix += phis[i]*phis[j]*kij
    return kmix


def liquid_mixture(
    target,
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    # DIPPR9H
    xs = target['pore.mole_fraction']
    MW = target.get_comp_vals(MWs)
    ks = target.get_comp_vals(ks)
    num = np.vstack([xs[k]*MW[k]/(ks[k]**2) for k in xs.keys()]).sum(axis=0)
    denom = np.vstack([xs[k]*MW[k] for k in xs.keys()]).sum(axis=0)
    temp = num/denom
    kmix = (1/temp)**0.5
    return kmix


def gas_mixture(
    target,
    T='pore.temperature',
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    # Wassiljew_Herning_Zipperer
    T = target[T]
    ys = target['pore.mole_fraction']
    kGs = target.get_comp_vals(ks)
    MWs = target.get_comp_vals(MWs)
    kmix = np.zeros_like(T)
    for i, ki in enumerate(ys.keys()):
        num = ys[ki]*kGs[ki]
        denom = 0.0
        for j, kj in enumerate(ys.keys()):
            A = np.sqrt(MWs[kj]/MWs[ki])
            denom += ys[ki]*A
        kmix += num/denom
    return kmix
