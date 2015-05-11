# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- viscosity
===============================================================================

Models for predicting phase viscosity

"""
import scipy as sp


def water(phase, **kwargs):
    r"""
    Calculates viscosity of pure water or seawater at atmospheric pressure
    using Eq. (22) given by Sharqawy et. al [1]_. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation pressure.

    Parameters
    ----------
    T, S: strings
        Property names where phase temperature and salinity are located.

    Returns
    -------
    mu_sw, the viscosity of water/seawater in [kg/m.s]

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 150 g/kg;
    ACCURACY: 1.5 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and Water
        Treatment, 2010.

    """
    T = phase['pore.temperature']
    try:
        S = phase['pore.salinity']
    except:
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


def reynolds(phase, uo, b, **kwargs):
    r"""
    Uses exponential model by Reynolds [1]_ for the temperature dependance of
    shear viscosity

    Parameters
    ----------
    u0, b : float, array_like
            Coefficients of the viscosity exponential model (mu = u0*Exp(-b*T)
            where T is the temperature in Kelvin

    [1] Reynolds O. (1886). Phil Trans Royal Soc London, v. 177, p.157.

    """
    T = phase['pore.temperature']
    value = uo*sp.exp(b*T)
    return value


def chung(phase, MW='molecular_weight', Tc='critical_temperature',
          Vc='critical_volume', **kwargs):
    r"""
    Uses Chung et al. [2]_ model to estimate viscosity for gases with low pressure
    (much less than the critical pressure) at conditions of interest

    Parameters
    ----------
    Vc :  float, array_like
        Critical volume of the gas (m3/kmol)
    Tc :  float, array_like
        Critical temperature of the gas (K)
    MW : float, array_like
        Molecular weight of the gas (kg/kmol)

    [2] Chung, T.H., Lee, L.L., and Starling, K.E., Applications of Kinetic Gas
        Theories and Multiparameter Correlation for Prediction of Dilute Gas
        Viscosity and Thermal Conductivity”, Ind. Eng. Chem. Fundam.23:8, 1984.

    """
    T = phase['pore.temperature']
    MW = phase['pore.'+MW]
    Tc = phase['pore.'+Tc]
    Vc = phase['pore.'+Vc]
    Tr = T / Tc
    Tstar = 1.2593*Tr
    A = 1.161415
    B = 0.14874
    C = 0.52487
    D = 0.77320
    E = 2.16178
    F = 2.43787
    omega = (A*(Tstar)**(-B)) + C*(sp.exp(-D*Tstar)) + E*(sp.exp(-F*Tstar))
    sigma = 0.809*(Vc**(1/3))
    value = 26.69E-9*sp.sqrt(MW*T)/(omega*sigma**2)
    return value
