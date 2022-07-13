import numpy as np
from openpnm.models.phase.mixtures import mixing_rule


__all__ = [
    'water_correlation',
    'liquid_mixture',
    'liquid_pure',
    'gas_mixture',
    'gas_pure',
]


def water_correlation(
    target,
    temperature='pore.temperature',
    salinity='pore.salinity'
):
    r"""
    Calculates viscosity of pure water or seawater at atmospheric pressure
    using Eq. (22) given by Sharqawy et. al [1]. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation
    pressure.

    Parameters
    ----------
    target : Phase
        The object for which these values are being calculated.
    temperature : str
        The dictionary key containing the temperature values. Temperature must
        be in Kelvin for this emperical equation to work.
    salinity : str
        The dictionary key containing the salinity values. Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    mu : ndarray
        Array containing viscosity of water or seawater in [kg/m.s] or [Pa.s]

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 150 g/kg;
    ACCURACY: 1.5 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
        Water Treatment, 2010.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    TC = T-273.15
    S = S/1000
    a1 = 1.5700386464E-01
    a2 = 6.4992620050E+01
    a3 = -9.1296496657E+01
    a4 = 4.2844324477E-05
    mu_w = a4 + 1/(a1*(TC+a2)**2+a3)
    a5 = 1.5409136040E+00
    a6 = 1.9981117208E-02
    a7 = -9.5203865864E-05
    a8 = 7.9739318223E+00
    a9 = -7.5614568881E-02
    a10 = 4.7237011074E-04
    A = a5 + a6*T + a7*T**2
    B = a8 + a9*T + a10*T**2
    mu_sw = mu_w*(1 + A*S + B*S**2)
    value = mu_sw
    return value


def air_correlation(
    target,
    temperature='pore.temperature',
    molar_density='pore.molar_density',
):
    r"""
    """
    raise Exception("This function does not work yet")
    # Get props from object
    T = target[temperature]
    rho = target[molar_density]
    # Declare given constants
    MW = 28.9586  # g/mol
    Tc = 132.6312  # K
    rho_c = 10447.7  # mol/m3
    sigma = 0.36  # nm
    e_k = 103.3  # K
    # Compute a few definitions
    delta = rho/rho_c
    T_star = T/e_k
    tau = np.atleast_2d(Tc/T)
    # Declare the summation variables
    ind = np.atleast_2d([0, 1, 2, 3, 4]).T
    N_i = np.atleast_2d([10.72, 1.122, 0.002019, -8.878, -0.02916]).T
    t_i = np.atleast_2d([0.2, 0.05, 2.4, 0.6, 3.6]).T
    d_i = np.atleast_2d([1, 4, 9, 1, 8]).T
    b_i = np.atleast_2d([0.431, -0.4623, 0.08406, 0.005341, -0.00331]).T
    gamma_i = np.atleast_2d([0, 0, 0, 1, 1]).T
    # Start crunching numbers
    omega = np.exp(np.sum(b_i*(np.atleast_2d(np.log(T_star))**(ind)), axis=0))
    mu_o = 0.9266958*(MW*T)**0.5/((sigma**2) * omega)
    temp = N_i * (tau**t_i) * (delta**d_i) * np.exp(-gamma_i*(delta**gamma_i))
    mu_t = np.sum(temp, axis=0)
    mu = mu_o + mu_t
    return mu


def gas_pure_stiel_thodos(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    critical_pressure='param.critical_pressure',
    molecular_weight='param.molecular_weight',
):
    T = target[temperature]
    Tc = target[critical_temperature]
    Pc = target[critical_pressure]/101325
    MW = target[molecular_weight]

    mu = np.zeros_like(T)
    Tr = T/Tc
    zeta = Tc**(1/6)/((MW**0.5)*(Pc**(2/3)))
    mu_hi = (17.78e-5*(4.58*Tr - 1.67)**0.625)/zeta
    mu_lo = (34e-5*Tr**0.94)/zeta
    mask = Tr > 1.5
    mu[mask] = mu_hi[mask]
    mu[~mask] = mu_lo[~mask]
    return mu/1000


def gas_pure(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    critical_pressure='param.critical_pressure',
    molecular_weight='param.molecular_weight',
):
    r"""
    """
    # Gharagheizi method from chemicals
    T = target[temperature]
    MW = target[molecular_weight]
    Pc = target[critical_pressure]
    Tc = target[critical_temperature]
    Tr = T/Tc
    mu = (1e-5)*Pc*Tr + (0.091 - 0.477/MW)*T + \
        MW*((1e-5)*Pc - 8.0*(MW**2)/(T**2))*(10.7639/Tc - 4.1929/T)
    mu = mu*1e-7
    mu = np.clip(mu, a_min=1e-7, a_max=np.inf)
    return mu


def gas_pure_chung(
    target,
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
    critical_temperature='param.critical_temperature',
    critical_volume='param.critical_volume'
):
    r"""
    Uses Chung et al. [1] model to estimate viscosity for gases at low
    pressure (much less than the critical pressure) at conditions of interest.

    Parameters
    ----------
    target : Phase
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.
    temperatre: str
        The dictionary key containing the temperature values (K)
    critical_temperature : str
        The dictionary key containing the temperature values (K)
    molecular_weight: str
        The dictionary key containing the molecular weight values (kg/mol)
    critical_volume : str
        The dictionary key containing the critical volume values (m3/kmol)

    Returns
    -------
    value : ndarray
        Array containing viscosity values based on Chung model in [Pa.s].

    References
    ----------
    [1] Chung, T.H., Lee, L.L., and Starling, K.E., Applications of Kinetic Gas
        Theories and Multiparameter Correlation for Prediction of Dilute Gas
        Viscosity and Thermal Conductivity”, Ind. Eng. Chem. Fundam.23:8, 1984.

    """
    T = target[temperature]
    MW = target[molecular_weight]
    Tc = target[critical_temperature]
    Vc = target[critical_volume]
    Tr = T / Tc
    Tstar = 1.2593*Tr
    A = 1.161415
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    omega = (A*(Tstar)**(-B)) + C*(np.exp(-D*Tstar)) + E*(np.exp(-F*Tstar))
    sigma = 0.809*(Vc**(1/3))
    value = 26.69E-9*np.sqrt(MW*T)/(omega*sigma**2)
    return value


def gas_mixture(
    target,
    viscosity='pore.viscosity',
    molecular_weight='param.molecular_weight',
):
    r"""
    Computes the viscosity of a gas mixture using a custom written version
    of ``chemicals.viscosity.Herning_Zipperer``

    Parameters
    ----------
    target : dict
        The openpnm object to which this model applies
    molecular_weight : str
        The location of the molecular weights of each species. The default is
        ``'param.molecular_weight'``
    viscosity : str
        The locaiton of the viscosities for each individual species.  The
        default is ``'pore.viscosity'``.

    Returns
    -------
    mu : ndarray
        An ndarray of Np length containing the viscosity of the mixture in
        each pore
    """
    MWs = target.get_comp_vals(molecular_weight)
    mus = target.get_comp_vals(viscosity)
    xs = target['pore.mole_fraction']
    num = 0.0
    denom = 0.0
    for k in xs.keys():
        num += xs[k]*mus[k]*MWs[k]**0.5
        denom += xs[k]*MWs[k]**0.5
    mu = num/denom
    return mu


def liquid_pure(
    target,
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
    critical_temperature='param.critical_temperature',
    critical_pressure='param.critical_pressure',
    acentric_factor='param.acentric_factor',
):
    r"""
    """
    # Letsou_Stiel
    T = target[temperature]
    MW = target[molecular_weight]
    Tc = target[critical_temperature]
    Pc = target[critical_pressure]
    omega = target[acentric_factor]
    Tr = T/Tc
    zeta = 2173.424 * (Tc**(1/6))/((MW**(0.5))*(Pc**(2/3)))
    zeta0 = (1.5174 - 2.135*Tr + 0.75*(Tr**2))*1e-5
    zeta1 = (4.2552 - 7.674*Tr + 3.4*(Tr**2))*1e-5
    mu = (zeta0 + omega*zeta1)/zeta
    return mu


def liquid_mixture(
    target,
    prop='pore.viscosity',
    mode='logarithmic',
    power=1,
):
    return mixing_rule(target=target, prop=prop, mode=mode, power=power)


liquid_mixture.__doc__ = mixing_rule.__doc__


if __name__ == "__main__":

    import chemicals as chem
    import openpnm as op
    from numpy.testing import assert_allclose

    pn = op.network.Demo()

    ch4 = op.phase.Species(network=pn, species='methane')
    ch4.add_model(propname='pore.viscosity',
                  model=gas_pure_stiel_thodos)
    mu_ref = chem.viscosity.Stiel_Thodos(
        T=ch4['pore.temperature'][0],
        Tc=ch4['param.critical_temperature'],
        Pc=ch4['param.critical_pressure'],
        MW=ch4['param.molecular_weight'],
    )
    mu_calc = ch4['pore.viscosity'][0]
    assert_allclose(mu_ref, mu_calc, rtol=1e-10, atol=0)

    co2 = op.phase.Species(network=pn, species='co2')
    co2.add_model(propname='pore.viscosity',
                  model=gas_pure)
    mu_ref = chem.viscosity.viscosity_gas_Gharagheizi(
        T=co2['pore.temperature'][0],
        Tc=co2['param.critical_temperature'],
        Pc=co2['param.critical_pressure'],
        MW=co2['param.molecular_weight'],
    )
    mu_calc = co2['pore.viscosity'][0]
    assert_allclose(mu_ref, mu_calc, rtol=1e-10, atol=0)

    fluegas = op.phase.GasMixture(network=pn, components=[co2, ch4])
    fluegas.y(ch4, 0.1)
    fluegas.y(co2, 0.9)
    fluegas.add_model(propname='pore.viscosity',
                      model=gas_mixture)
    mu_ref = chem.viscosity.Herning_Zipperer(
        zs=np.vstack(list(fluegas['pore.mole_fraction'].values()))[:, 0],
        mus=np.vstack(list(fluegas.get_comp_vals('pore.viscosity').values()))[:, 0],
        MWs=np.vstack(list(fluegas.get_comp_vals('param.molecular_weight').values()))[:, 0],
    )
    mu_calc = fluegas['pore.viscosity'][0]
    assert_allclose(mu_ref, mu_calc, rtol=1e-10, atol=0)

    cbz = op.phase.Species(network=pn, species='chlorobenzene')
    cbz.add_model(propname='pore.viscosity',
                  model=liquid_pure)
    mu_calc = cbz['pore.viscosity'][0]
    mu_ref = chem.viscosity.Letsou_Stiel(
        T=cbz['pore.temperature'][0],
        Tc=cbz['param.critical_temperature'],
        Pc=cbz['param.critical_pressure'],
        MW=cbz['param.molecular_weight'],
        omega=cbz['param.acentric_factor'],
    )
    assert_allclose(mu_ref, mu_calc, rtol=1e-10, atol=0)























