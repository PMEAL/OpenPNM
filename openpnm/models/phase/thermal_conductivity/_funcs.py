import numpy as np
from openpnm.utils import Docorator
from chemicals import rho_to_Vm


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
    temperature="pore.temperature",
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


def gas_pure(
    target,
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
    boiling_temperature='param.boiling_temperature',
    critical_pressure='param.critical_pressure',
    acentric_factor='param.acentric_factor',
):
    # gharagheizi method
    T = target[temperature]
    MW = target[molecular_weight]
    Tb = target[boiling_temperature]
    # The following correction suggested by chemicals package author
    Pc = target[critical_pressure]/10000
    omega = target[acentric_factor]
    B = (T + (2.0*omega + 2.0*T - 2.0*T*(2.0*omega + 3.2825)/Tb + 3.2825)
         / (2.0*omega + T - T*(2.0*omega + 3.2825)/Tb + 3.2825)
         - T*(2.0*omega + 3.2825)/Tb)
    A = (2*omega + T - (2*omega + 3.2825)*T/Tb + 3.2825)/(0.1*MW*Pc*T) \
        * (3.9752*omega + 0.1*Pc + 1.9876*B + 6.5243)**2
    k = 7.9505e-4 + 3.989e-5*T - 5.419e-5*MW + 3.989e-5*A
    return k


def liquid_pure(
    target,
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
    boiling_temperature='param.boiling_temperature',
    critical_pressure='param.critical_pressure',
    acentric_factor='param.acentric_factor',
):
    # gharagheizi metho
    T = target[temperature]
    MW = target[molecular_weight]
    Tb = target[boiling_temperature]
    Pc = target[critical_pressure]/100000
    omega = target[acentric_factor]
    B = 16.0407*MW + 2.0*Tb - 27.9074
    A = 3.8588*(MW**8)*(1.0045*B + 6.5152*MW - 8.9756)
    k = (1e-4)*(10*omega + 2*Pc - 2*T + 4 + 1.908*(Tb + 1.009*(B**2)/(MW**2))
                + 3.9287*(MW**4)/(B**4) + A/(B**8))
    return k


@docstr.dedent
def gas_pure_chung(
    target,
    temperature="pore.temperature",
    viscosity="pore.viscosity",
    Cv="pore.constant_volume_molar_heat_capacity",
    acentric_factor='param.acentric_factor',
    molecular_weight='param.molecular_weight',
    critical_temperature='param.critical_temperature',
):
    r"""
    Uses Chung et al. model to estimate thermal conductivity for gases with
    low pressure(<10 bar) from first principles at conditions of interest

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    viscosity : str
        The dictionary key containing the viscosity values (Pa.s)
    Cv : str
        Dictionary key containing the heat capacity at constant volume
        (J/(mol.K))

    Returns
    -------
    value : ndarray
        A numpy ndarray containing thermal conductivity values in [W/m.K]

    """
    omega = target[acentric_factor]
    MW = target[molecular_weight]
    Tc = target[critical_temperature]
    Cv = target[Cv]
    T = target[temperature]
    mu = target[viscosity]
    R = 8.314
    Tr = T / Tc
    z = 2.0 + 10.5 * Tr ** 2
    beta = 0.7862 - 0.7109 * omega + 1.3168 * omega ** 2
    alpha = Cv / R - 3 / 2
    s = 1 + alpha * (
        (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
        / (0.6366 + beta * z + 1.061 * alpha * beta)
    )
    value = 3.75 * s * (mu) * R / (MW)
    return value


def liquid_pure_sato_riedel(
    target,
    temperature="pore.temperature",
    critical_temperature='param.critical_temperature',
    molecular_weight='param.molecular_weight',
    boiling_temperature='param.boiling_temperature',
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
    T = target[temperature]
    Tc = target[critical_temperature]
    MW = target[molecular_weight]
    Tbr = target[boiling_temperature] / Tc
    Tr = T / Tc
    value = ((1.1053 / ((MW) ** 0.5)) * (3 + 20 * (1 - Tr) ** (2 / 3))
             / (3 + 20 * (1 - Tbr) ** (2 / 3)))
    return value


def liquid_mixture_DIPPR9I(
    target,
    density='pore.density',
    thermal_conductivity='pore.thermal_conductivity',
    molecular_weight='param.molecular_weight',
):
    xs = target['pore.mole_fraction']
    kLs = [c[thermal_conductivity] for c in target.components.values()]
    # kL = numba_vectorized.DIPPR9I(xs, kLs)  # Another one that doesn't work
    Vm = [rho_to_Vm(c[density], c[molecular_weight])
          for c in target.components.values()]
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
    thermal_conductivity='pore.thermal_conductivity',
    molecular_weight='param.molecular_weight',
):
    # DIPPR9H
    xs = target['pore.mole_fraction']
    MW = target.get_comp_vals('param.molecular_weight')
    ks = target.get_comp_vals('pore.thermal_conductivity')
    num = np.vstack([xs[k]*MW[k]/(ks[k]**2) for k in xs.keys()]).sum(axis=0)
    denom = np.vstack([xs[k]*MW[k] for k in xs.keys()]).sum(axis=0)
    temp = num/denom
    kmix = (1/temp)**0.5
    return kmix


def gas_mixture(
    target,
    temperature='pore.temperature',
    thermal_conductivity='pore.thermal_conductivity',
    molecular_weight='param.molecular_weight',
):
    # Wassiljew_Herning_Zipperer
    T = target[temperature]
    ys = target['pore.mole_fraction']
    kGs = target.get_comp_vals(thermal_conductivity)
    MWs = target.get_comp_vals(molecular_weight)
    kmix = np.zeros_like(T)
    for i, ki in enumerate(ys.keys()):
        num = ys[ki]*kGs[ki]
        denom = 0.0
        for j, kj in enumerate(ys.keys()):
            A = np.sqrt(MWs[kj]/MWs[ki])
            denom += ys[ki]*A
        kmix += num/denom
    return kmix


if __name__ == "__main__":
    import chemicals as chem
    import openpnm as op
    from numpy.testing import assert_allclose

    pn = op.network.Demo()

    chem.thermal_conductivity.Gharagheizi_gas

    ch4 = op.phase.Species(network=pn, species='methane')
    ch4.add_model(propname='pore.thermal_conductivity',
                  model=gas_pure)
    k_calc = ch4['pore.thermal_conductivity'][0]
    k_ref = chem.thermal_conductivity.Gharagheizi_gas(
        T=ch4['pore.temperature'][0],
        MW=ch4['param.molecular_weight'],
        Tb=ch4['param.boiling_temperature'],
        Pc=ch4['param.critical_pressure'],
        omega=ch4['param.acentric_factor'],
    )
    assert_allclose(k_ref, k_calc, rtol=1e-10, atol=0)

    h2o = op.phase.Species(network=pn, species='water')
    h2o.add_model(propname='pore.thermal_conductivity',
                  model=liquid_pure)
    k_calc = h2o['pore.thermal_conductivity'][0]
    k_ref = chem.thermal_conductivity.Gharagheizi_liquid(
        T=h2o['pore.temperature'][0],
        MW=h2o['param.molecular_weight'],
        Tb=h2o['param.boiling_temperature'],
        Pc=h2o['param.critical_pressure'],
        omega=h2o['param.acentric_factor'],
    )
    assert_allclose(k_ref, k_calc, rtol=1e-10, atol=0)

    h2o = op.phase.Species(network=pn, species='water')
    h2o.add_model(propname='pore.thermal_conductivity',
                  model=liquid_pure)
    etoh = op.phase.Species(network=pn, species='ethanol')
    etoh.add_model(propname='pore.thermal_conductivity',
                   model=liquid_pure)

    vodka = op.phase.LiquidMixture(network=pn, components=[h2o, etoh])
    vodka.x(h2o.name, 0.5)
    vodka.x(etoh.name, 0.5)
    vodka.add_model(propname='pore.thermal_conductivity',
                    model=liquid_mixture_DIPPR9H)
    k_ref = chem.thermal_conductivity.DIPPR9H(
        ws=np.vstack(list(vodka['pore.mole_fraction'].values()))[:, 0],
        ks=np.vstack(list(vodka.get_comp_vals('pore.thermal_conductivity').values()))[:, 0],
    )
    k_calc = vodka['pore.thermal_conductivity'][0]
    assert_allclose(k_ref, k_calc, rtol=1e-10, atol=0)






























