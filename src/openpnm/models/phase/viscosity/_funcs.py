import numpy as np
from openpnm.models.phase.mixtures import mixing_rule
from openpnm.models.phase import _phasedocs


__all__ = [
    'water_correlation',
    'gas_pure_st',
    'gas_pure_gesmr',
    'gas_mixture_hz',
    'liquid_pure_ls',
    'liquid_mixture_xweighted',
]


@_phasedocs
def water_correlation(
    phase,
    T='pore.temperature',
    salinity='pore.salinity'
):
    r"""
    Calculates viscosity of pure water or seawater at atmospheric pressure.

    This correlation uses Eq. (22) given by Sharqawy et. al [1]. Values at
    temperature higher than the normal boiling temperature are calculated
    at the saturation pressure.

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(salinity)s

    Returns
    -------
    mu : ndarray
        Array containing viscosity of water or seawater in [kg/m.s] or [Pa.s]

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 150 g/kg;
    ACCURACY: 1.5 percent

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
        Water Treatment, 2010.

    """
    T = phase[T]
    if salinity in phase.keys():
        S = phase[salinity]
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


@_phasedocs
def air_correlation(
    phase,
    T='pore.temperature',
    n_V='pore.molar_density',
):
    r"""
    Calculates the viscosity of air at given conditions

    This model uses the correlation from [1].

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(n_V)s

    Returns
    -------

    References
    ----------
    [1] I forget

    """
    raise Exception("This function does not work yet")
    # Get props from object
    T = phase[T]
    rho = phase[n_V]
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


@_phasedocs
def gas_pure_st(
    phase,
    T='pore.temperature',
    Tc='param.critical_temperature',
    Pc='param.critical_pressure',
    MW='param.molecular_weight',
):
    r"""
    Calculates the viscosity of a pure gas using the correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s
    %(Pc)s
    %(MW)s

    Returns
    -------

    References
    ----------
    [1] Stiel and Thodos

    """
    #  stiel-thodos
    T = phase[T]
    Tc = phase[Tc]
    Pc = phase[Pc]/101325
    MW = phase[MW]

    mu = np.zeros_like(T)
    Tr = T/Tc
    zeta = Tc**(1/6)/((MW**0.5)*(Pc**(2/3)))
    mu_hi = (17.78e-5*(4.58*Tr - 1.67)**0.625)/zeta
    mu_lo = (34e-5*Tr**0.94)/zeta
    mask = Tr > 1.5
    mu[mask] = mu_hi[mask]
    mu[~mask] = mu_lo[~mask]
    return mu/1000


@_phasedocs
def gas_pure_gesmr(
    phase,
    T='pore.temperature',
    Tc='param.critical_temperature',
    Pc='param.critical_pressure',
    MW='param.molecular_weight',
):
    r"""
    Calculates the viscosity of a pure gas using the correlation  [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(Tc)s
    %(Pc)s
    %(MW)s

    Returns
    -------

    References
    ----------
    [1] Gharagheizi et al
    """
    T = phase[T]
    MW = phase[MW]
    Pc = phase[Pc]
    Tc = phase[Tc]
    Tr = T/Tc
    mu = (1e-5)*Pc*Tr + (0.091 - 0.477/MW)*T + \
        MW*((1e-5)*Pc - 8.0*(MW**2)/(T**2))*(10.7639/Tc - 4.1929/T)
    mu = mu*1e-7
    mu = np.clip(mu, a_min=1e-7, a_max=np.inf)
    return mu


@_phasedocs
def gas_mixture_hz(
    phase,
    mus='pore.viscosity.*',
    MWs='param.molecular_weight.*',
):
    r"""
    Computes the viscosity of a gas mixture using the correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(mus)s
    %(MWs)s

    Returns
    -------
    mu : ndarray
        An ndarray of Np length containing the viscosity of the mixture in
        each pore

    References
    ----------
    [1] Herning and Zipperer
    """
    MWs = phase.get_comp_vals(MWs)
    mus = phase.get_comp_vals(mus)
    xs = phase['pore.mole_fraction']
    num = 0.0
    denom = 0.0
    for k in xs.keys():
        num += xs[k]*mus[k]*MWs[k]**0.5
        denom += xs[k]*MWs[k]**0.5
    mu = num/denom
    return mu


@_phasedocs
def liquid_pure_ls(
    phase,
    T='pore.temperature',
    MW='param.molecular_weight',
    Tc='param.critical_temperature',
    Pc='param.critical_pressure',
    omega='param.acentric_factor',
):
    r"""
    Computes the viscosity of a pure liquid using the correlation in [1]

    Parameters
    ----------
    %(phase)s
    %(T)s
    %(MW)s
    %(Tc)s
    %(Pc)s
    %(omega)s

    Returns
    -------
    mu : ndarray
        An ndarray of Np length containing the viscosity of the mixture in
        each pore

    References
    ----------
    [1] Letsou and Stiel

    """
    T = phase[T]
    MW = phase[MW]
    Tc = phase[Tc]
    Pc = phase[Pc]
    omega = phase[omega]
    Tr = T/Tc
    zeta = 2173.424 * (Tc**(1/6))/((MW**(0.5))*(Pc**(2/3)))
    zeta0 = (1.5174 - 2.135*Tr + 0.75*(Tr**2))*1e-5
    zeta1 = (4.2552 - 7.674*Tr + 3.4*(Tr**2))*1e-5
    mu = (zeta0 + omega*zeta1)/zeta
    return mu


def liquid_mixture_xweighted(
    phase,
    prop='pore.viscosity.*',
    mode='logarithmic',
    power=1,
):
    return mixing_rule(phase=phase, prop=prop, mode=mode, power=power)


liquid_mixture_xweighted.__doc__ = mixing_rule.__doc__
