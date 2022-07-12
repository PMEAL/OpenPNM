import numpy as np
from openpnm.models.phase.mixtures import mixing_rule


__all__ = [
    'water_correlation',
    'liquid_chung',
    'liquid_mixture',
    'gas_mixture',
    'gas_gharagheizi',
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


def liquid_chung(
    target,
    temperature='pore.temperature',
    mol_weight='param.molecular_weight',
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
    mol_weight: str
        The dictionary key containing the molecular weight values (kg/mol)
    critical_volume : str
        The dictionary key containing the critical volume values (m3/kmol)

    Returns
    -------
    value : ndarray
        Array containing viscosity values based on Chung model [kg/m.s].

    References
    ----------
    [1] Chung, T.H., Lee, L.L., and Starling, K.E., Applications of Kinetic Gas
        Theories and Multiparameter Correlation for Prediction of Dilute Gas
        Viscosity and Thermal Conductivity”, Ind. Eng. Chem. Fundam.23:8, 1984.

    """
    T = target[temperature]
    MW = target[mol_weight]
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


def gas_gharagheizi(
    target,
):
    pass


def gas_mixture(
    target,
    MWs='param.molecular_weight',
    mus='pore.viscosity',
):
    r"""
    Computes the viscosity of a gas mixture using a custom written version
    of ``chemicals.viscosity.HerningZipperer``

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
    MWs = target.get_comp_vals(MWs)
    mus = target.get_comp_vals(mus)
    xs = target['pore.mole_fraction']
    num = 0.0
    denom = 0.0
    for k in xs.keys():
        num += xs[k]*mus[k]*MWs[k]**0.5
        denom += xs[k]*MWs[k]**0.5
    mu = num/denom
    return mu


def liquid_mixture(
    target,
    prop='pore.viscosity',
    mode='logarithmic',
    power=1,
):
    return mixing_rule(target=target, prop=prop, mode=mode, power=power)


liquid_mixture.__doc__ = mixing_rule.__doc__
