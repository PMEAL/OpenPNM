import numpy as np
from openpnm.models.phase import _phasedocs


__all__ = [
    'water_correlation',
    'gas_pure_gismr',
    'liquid_pure_gismr',
    'liquid_pure_sr',
    'gas_mixture_whz',
    'liquid_mixture_DIPPR9H',
]


@_phasedocs
def water_correlation(
    phase,
    T="pore.temperature",
    salinity="pore.salinity",
):
    r"""
    Calculates thermal conductivity of pure water or seawater at atmospheric
    pressure using the correlation given [1].

    Values at temperature higher the normal boiling temperature are calculated
    at the saturation pressure.

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(salinity)s

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity of water/seawater in [W/m.K]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
    thousand). The correlation is valid for 273 < T < 453 K and
    0 < S < 160 g/kg within 3 percent accuracy.

    References
    ----------
    [1] D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.

    """
    T = phase[T]
    if salinity in phase.keys():
        S = phase[salinity]
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


@_phasedocs
def gas_pure_gismr(
    phase,
    T='pore.temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    r"""
    Calculate the thermal conductivty of a pure gas using the correlation in
    [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(MW)s
    %(Tb)s
    %(Pc)s
    %(omega)s

    Returns
    -------

    References
    ----------
    [1] gharagheizi method: doi:10.1002/aic.13938

    """
    T = phase[T]
    MW = phase[MW]
    Tb = phase[Tb]
    # The following correction suggested by chemicals package author
    Pc = phase[Pc]/10000
    omega = phase[omega]
    B = (T + (2.0*omega + 2.0*T - 2.0*T*(2.0*omega + 3.2825)/Tb + 3.2825)
         / (2.0*omega + T - T*(2.0*omega + 3.2825)/Tb + 3.2825)
         - T*(2.0*omega + 3.2825)/Tb)
    A = (2*omega + T - (2*omega + 3.2825)*T/Tb + 3.2825)/(0.1*MW*Pc*T) \
        * (3.9752*omega + 0.1*Pc + 1.9876*B + 6.5243)**2
    k = 7.9505e-4 + 3.989e-5*T - 5.419e-5*MW + 3.989e-5*A
    return k


@_phasedocs
def liquid_pure_gismr(
    phase,
    T='pore.temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    r"""
    Calculates the thermal conductivity of a pure liquid using the correlation
    in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(MW)s
    %(Tb)s
    %(Pc)s
    %(omega)s

    Returns
    -------

    References
    ----------
    [1] gharagheizi method: doi:10.1002/aic.13938

    """
    T = phase[T]
    MW = phase[MW]
    Tb = phase[Tb]
    Pc = phase[Pc]/100000
    omega = phase[omega]
    B = 16.0407*MW + 2.0*Tb - 27.9074
    A = 3.8588*(MW**8)*(1.0045*B + 6.5152*MW - 8.9756)
    k = (1e-4)*(10*omega + 2*Pc - 2*T + 4 + 1.908*(Tb + 1.009*(B**2)/(MW**2))
                + 3.9287*(MW**4)/(B**4) + A/(B**8))
    return k


@_phasedocs
def liquid_pure_sr(
    phase,
    T="pore.temperature",
    Tc='param.critical_temperature',
    MW='param.molecular_weight',
    Tb='param.boiling_temperature',
):
    r"""
    Calculates the thermal conductivity of a pure liquid using the correlation
    in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s
    %(MW)s
    %(Tb)s

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity values in [W/m.K]

    References
    ----------
    [1] Sato et al

    """
    T = phase[T]
    Tc = phase[Tc]
    MW = phase[MW]
    Tbr = phase[Tb]/Tc
    Tr = T / Tc
    value = ((1.1053 / ((MW) ** 0.5)) * (3 + 20 * (1 - Tr) ** (2 / 3))
             / (3 + 20 * (1 - Tbr) ** (2 / 3)))
    return value


@_phasedocs
def liquid_mixture_DIPPR9I(
    phase,
    rhos='pore.density.*',
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    r"""
    Calculates the thermal conductivity of a liquid mixture using the
    correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(rhos)s
    %(ks)s
    %(MWs)s

    Returns
    -------

    References
    ----------
    [1] DIPPR9I

    """
    raise NotImplementedError("This function is not ready yet")
    from chemicals import rho_to_Vm
    xs = phase['pore.mole_fraction']
    kLs = phase.get_comp_vals(ks)
    # kL = numba_vectorized.DIPPR9I(xs, kLs)  # Another one that doesn't work
    Vm = [rho_to_Vm(c.get_comp_vals(rhos), c.get_comp_vals(MWs))
          for c in phase.components.keys()]
    denom = np.sum([xs[i]*Vm[i] for i in range(len(xs))], axis=0)
    phis = np.array([xs[i]*Vm[i] for i in range(len(xs))])/denom
    kij = 2/np.sum([1/kLs[i] for i in range(len(xs))], axis=0)
    kmix = np.zeros_like(xs[0])
    N = len(xs)
    for i in range(N):
        for j in range(N):
            kmix += phis[i]*phis[j]*kij
    return kmix


@_phasedocs
def liquid_mixture_DIPPR9H(
    phase,
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    r"""
    Calculates the thermal conductivity of a liquid mixture using the
    correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(ks)s
    %(MWs)s

    Returns
    -------

    References
    ----------
    [1] DIPPR9H

    """
    xs = phase['pore.mole_fraction']
    MW = phase.get_comp_vals(MWs)
    ks = phase.get_comp_vals(ks)
    num = np.vstack([xs[k]*MW[k]/(ks[k]**2) for k in xs.keys()]).sum(axis=0)
    denom = np.vstack([xs[k]*MW[k] for k in xs.keys()]).sum(axis=0)
    temp = num/denom
    kmix = (1/temp)**0.5
    return kmix


@_phasedocs
def gas_mixture_whz(
    phase,
    T='pore.temperature',
    ks='pore.thermal_conductivity.*',
    MWs='param.molecular_weight.*',
):
    r"""
    Calculates the viscosity of a gas mixture using the correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(ks)s
    %(MWs)s

    Returns
    -------

    References
    ----------
    [1] Wassiljew, Herning & Zipperer
    """
    T = phase[T]
    ys = phase['pore.mole_fraction']
    kGs = phase.get_comp_vals(ks)
    MWs = phase.get_comp_vals(MWs)
    kmix = np.zeros_like(T)
    for i, ki in enumerate(ys.keys()):
        num = ys[ki]*kGs[ki]
        denom = 0.0
        for j, kj in enumerate(ys.keys()):
            A = np.sqrt(MWs[kj]/MWs[ki])
            denom += ys[ki]*A
        kmix += num/denom
    return kmix
