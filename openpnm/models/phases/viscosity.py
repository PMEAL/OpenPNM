r"""
"""
import numpy as np


def water(target, temperature='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates viscosity of pure water or seawater at atmospheric pressure
    using Eq. (22) given by Sharqawy et. al [1]. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation
    pressure.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated. This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    temperature : string
        The dictionary key containing the temperature values. Temperature must
        be in Kelvin for this emperical equation to work. Can be either a pore
        or throat array.

    salinity : string
        The dictionary key containing the salinity values. Salinity must be
        expressed in g of salt per kg of solution (ppt). Can be either a
        pore or throat array, but must be consistent with ``temperature``.

    Returns
    -------
    value : NumPy ndarray
        Array containing viscosity of water/seawater in [kg/m.s]

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


def reynolds(target, u0, b, temperature='pore.temperature'):
    r"""
    Uses exponential model by Reynolds [1] for the temperature dependance of
    shear viscosity

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    u0, b : float, array_like
            Coefficients of the viscosity exponential model (mu = u0*Exp(-b*T)
            where T is the temperature in Kelvin

    temperature : string
        The dictionary key containing the temperature values (K).  Can be
        either a pore or throat array.

    Returns
    -------
    value : NumPy ndarray
        Array containing viscosity values based on Reynolds model.

    [1] Reynolds O. (1886). Phil Trans Royal Soc London, v. 177, p.157.

    """
    value = u0*np.exp(b*target[temperature])
    return value


def chung(target, temperature='pore.temperature',
          mol_weight='pore.molecular_weight',
          critical_temperature='pore.critical_temperature',
          critical_volume='pore.critical_volume'):
    r"""
    Uses Chung et al. [1] model to estimate viscosity for gases at low
    pressure (much less than the critical pressure) at conditions of interest.

    Parameters
    ----------
    target : OpenPNM Object
        The object for which these values are being calculated.  This
        controls the length of the calculated array, and also provides
        access to other necessary thermofluid properties.

    temperatre: string
        The dictionary key containing the temperature values (K)

    critical_temperature : string
        The dictionary key containing the temperature values (K)

    mol_weight: string
        The dictionary key containing the molecular weight values (kg/mol)

    critical_volume : string
        The dictionary key containing the critical volume values (m3/kmol)

    Returns
    -------
    value : NumPy ndarray
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
