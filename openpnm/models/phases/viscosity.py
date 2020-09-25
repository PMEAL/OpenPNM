r"""

.. autofunction:: openpnm.models.phases.viscosity.water
.. autofunction:: openpnm.models.phases.viscosity.reynolds
.. autofunction:: openpnm.models.phases.viscosity.chung

"""
import numpy as np


def powerlaw_fluid(target,
                   pore_area="pore.area",
                   throat_area="throat.area",
                   pore_viscosity_min="pore.viscosity_min",
                   throat_viscosity_min="throat.viscosity_min",
                   pore_viscosity_max="pore.viscosity_max",
                   throat_viscosity_max="throat.viscosity_max",
                   conduit_lengths="throat.conduit_lengths",
                   conduit_shape_factors="throat.flow_shape_factors",
                   pore_consistency="pore.consistency",
                   throat_consistency="throat.consistency",
                   pore_flow_index="pore.flow_index",
                   throat_flow_index="throat.flow_index",
                   pore_pressure="pore.pressure",
                   shape="cylinder"):
    r"""

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_area : string
        Dictionary key of the pore area values
    throat_area : string
        Dictionary key of the throat area values
    pore_viscosity_min : string
        Dictionary key of the pore minimum viscosity values
    throat_viscosity_min : string
        Dictionary key of the throat minimum viscosity values
    pore_viscosity_max : string
        Dictionary key of the pore maximum viscosity values
    throat_viscosity_max : string
        Dictionary key of the throat maximum viscosity values
    conduit_lengths : string
        Dictionary key of the conduit length values
    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values
    pore_consistency : string
        Dictionary key of the pore fluid consistency values
    throat_consistency : string
        Dictionary key of the throat fluid consistency values
    pore_flow_index : string
        Dictionary key of the pore fluid flow index values
    throat_flow_index : string
        Dictionary key of the throat fluid flow index values
    pore_pressure : string
        Dictionary key of the pore pressure values

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + ".pore1"][throats]
    Lt = network[conduit_lengths + ".throat"][throats]
    L2 = network[conduit_lengths + ".pore2"][throats]
    # # Preallocating g
    # g1, g2, gt = np.zeros((3, len(Lt)))
    # # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    # g1[~m1] = g2[~m2] = gt[~mt] = np.inf

    # Check if pressure field exists
    try:
        phase[pore_pressure]
    except KeyError:
        phase[pore_pressure] = 0
    P = phase[pore_pressure]

    # Interpolate pore phase property values to throats
    try:
        mu_mint = phase[throat_viscosity_min][throats]
    except KeyError:
        mu_mint = phase.interpolate_data(propname=pore_viscosity_min)[throats]
    try:
        mu_maxt = phase[throat_viscosity_max][throats]
    except KeyError:
        mu_maxt = phase.interpolate_data(propname=pore_viscosity_max)[throats]
    try:
        Ct = phase[throat_consistency][throats]
    except KeyError:
        Ct = phase.interpolate_data(propname=pore_consistency)[throats]
    try:
        nt = phase[throat_flow_index][throats]
    except KeyError:
        nt = phase.interpolate_data(propname=pore_flow_index)[throats]
    # Interpolate throat phase property values to pores
    try:
        mu_min1 = phase[pore_viscosity_min][cn[:, 0]]
        mu_min2 = phase[pore_viscosity_min][cn[:, 1]]
    except KeyError:
        mu_min1 = phase.interpolate_data(propname=throat_viscosity_min)[cn[:, 0]]
        mu_min2 = phase.interpolate_data(propname=throat_viscosity_min)[cn[:, 1]]
    try:
        mu_max1 = phase[pore_viscosity_max][cn[:, 0]]
        mu_max2 = phase[pore_viscosity_max][cn[:, 1]]
    except KeyError:
        mu_max1 = phase.interpolate_data(propname=throat_viscosity_max)[cn[:, 0]]
        mu_max2 = phase.interpolate_data(propname=throat_viscosity_max)[cn[:, 1]]
    try:
        C1 = phase[pore_consistency][cn[:, 0]]
        C2 = phase[pore_consistency][cn[:, 1]]
    except KeyError:
        C1 = phase.interpolate_data(propname=throat_consistency)[cn[:, 0]]
        C2 = phase.interpolate_data(propname=throat_consistency)[cn[:, 1]]
    try:
        n1 = phase[pore_flow_index][cn[:, 0]]
        n2 = phase[pore_flow_index][cn[:, 1]]
    except KeyError:
        n1 = phase.interpolate_data(propname=throat_flow_index)[cn[:, 0]]
        n2 = phase.interpolate_data(propname=throat_flow_index)[cn[:, 1]]
    # Interpolate pore pressure values to throats
    Pt = phase.interpolate_data(propname=pore_pressure)[throats]

    # Pressure differences dP
    dP1 = np.absolute(P[cn[:, 0]] - Pt)
    dP2 = np.absolute(P[cn[:, 1]] - Pt)
    dPt = np.absolute(np.diff(P[cn], axis=1).squeeze())

    dP1 = dP1.clip(min=1e-20)
    dP2 = dP2.clip(min=1e-20)
    dPt = dPt.clip(min=1e-20)

    # Apparent viscosities
    mu1, mu2, mut = np.zeros((3, len(Lt)))

    mu1[m1] = mu_min1+(dP1 ** (1 - 1 / n1) * C1 ** (1 / n1))[m1] / (
        (4 * n1 / (3 * n1 + 1))[m1]
        * (2 * L1[m1] / ((A1[m1] / np.pi) ** 0.5)) ** (1 - 1 / n1[m1]))

    mu2[m2] = mu_min2+(dP2 ** (1 - 1 / n2) * C2 ** (1 / n2))[m2] / (
        (4 * n2 / (3 * n2 + 1))[m2]
        * (2 * L2[m2] / ((A2[m2] / np.pi) ** 0.5)) ** (1 - 1 / n2[m2]))

    mut[mt] = mu_mint+(dPt ** (1 - 1 / nt) * Ct ** (1 / nt))[mt] / (
        (4 * nt / (3 * nt + 1))[mt]
        * (2 * Lt[mt] / ((At[mt] / np.pi) ** 0.5)) ** (1 - 1 / nt[mt]))

    # Bound the apparent viscosity
    mu1[m1] = np.minimum(np.maximum(mu1[m1], mu_min1[m1]), mu_max1[m1])
    mu2[m2] = np.minimum(np.maximum(mu2[m2], mu_min2[m2]), mu_max2[m2])
    mut[mt] = np.minimum(np.maximum(mut[mt], mu_mint[mt]), mu_maxt[mt])

    mu = {'pore1': mu1, 'pore2': mu2, 'throat': mut}

    return mu


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
