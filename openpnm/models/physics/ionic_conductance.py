r"""
Pore-scale models for calculating ionic conductance of conduits.
"""
import numpy as _np
from openpnm.utils import logging
from openpnm.models.physics.utils import generic_transport_conductance
logger = logging.getLogger(__name__)

__all__ = ["generic_ionic_poisson_laplace", "poisson", "electroneutrality"]


def generic_ionic_poisson_laplace(target,
                                  pore_conductivity='pore.permittivity',
                                  throat_conductivity='throat.permittivity',
                                  size_factors='throat.diffusive_size_factors'):
    r"""
    Calculate the ionic conductance of conduits in network (using the Poisson
    equation for charge conservation), where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    conduit_lengths : string
        Dictionary key of the conduit length values

    size_factors : string
        Dictionary key of the conduit DIFFUSION size shape factor values

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties already
    be calculated.

    """
    epsilon0 = 8.854187817e-12
    gic = epsilon0*generic_transport_conductance(target=target,
                                                 pore_conductivity=pore_conductivity,
                                                 throat_conductivity=throat_conductivity,
                                                 size_factors=size_factors)
    return gic


def generic_ionic_electroneutrality(target,
                                    pore_area='pore.area',
                                    throat_area='throat.area',
                                    pore_diffusivity='pore.diffusivity',
                                    throat_diffusivity='throat.diffusivity',
                                    conduit_lengths='throat.conduit_lengths',
                                    size_factors="throat.diffusive_size_factors",
                                    pore_volume='pore.volume',
                                    pore_temperature='pore.temperature',
                                    throat_temperature='throat.temperature',
                                    pore_valence='pore.valence',
                                    throat_valence='throat.valence',
                                    pore_concentration='pore.concentration',
                                    ions=[]):
    r"""
    Calculate the ionic conductance of conduits in network (assuming
    electroneutrality for charge conservation), where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

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

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    size_factors : string
        Dictionary key of the conduit DIFFUSION size factor values

    pore_volume : string
        Dictionary key of the pore volume values

    pore_temperature : string
        Dictionary key of the pore temperature values

    throat_temperature : string
        Dictionary key of the throat temperature values

    pore_valence : string
       Dictionary key of the pore ionic species valence values

    throat_valence : string
       Dictionary key of the throat ionic species valence values

    pore_concentration : string
       Dictionary key of the pore ionic species concentration values

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F_integral = network[size_factors]
    # Electroneutrality
    F = 96485.3329
    R = 8.3145
    # Getting pores volumes
    Vol1 = network[pore_volume][cn[:, 0]]
    Vol2 = network[pore_volume][cn[:, 1]]
    # Preallocating g
    g1, g2, gt = _np.zeros((3, len(cn)))
    # Interpolate pore phase property values to throats
    Tt = phase[throat_temperature][throats]
    try:
        T1 = phase[pore_temperature][cn[:, 0]]
        T2 = phase[pore_temperature][cn[:, 1]]
    except KeyError:
        T1 = phase.interpolate_data(propname=throat_temperature)[cn[:, 0]]
        T2 = phase.interpolate_data(propname=throat_temperature)[cn[:, 1]]
    # Iterate over all ions present in the solution
    if ions == []:
        logger.error('List of ions must be provided')
    for i in ions:
        i = '.'+i
        # Check if a concetration field is defined
        try:
            c1 = phase[pore_concentration+i][cn[:, 0]]
            c2 = phase[pore_concentration+i][cn[:, 1]]
        except KeyError:
            c1 = _np.zeros((network.Np))[cn[:, 0]]
            c2 = _np.zeros((network.Np))[cn[:, 1]]
        ct = (c1*Vol1 + c2*Vol2)/(Vol1 + Vol2)
        Dt = phase[throat_diffusivity+i][throats]
        Vt = phase[throat_valence+i][throats]

        try:
            D1 = phase[pore_diffusivity+i][cn[:, 0]]
            D2 = phase[pore_diffusivity+i][cn[:, 1]]
            V1 = phase[pore_valence+i][cn[:, 0]]
            V2 = phase[pore_valence+i][cn[:, 1]]
        except KeyError:
            D1 = phase.interpolate_data(
                propname=throat_diffusivity+i)[cn[:, 0]]
            D2 = phase.interpolate_data(
                propname=throat_diffusivity+i)[cn[:, 1]]
            V1 = phase.interpolate_data(
                propname=throat_valence+i)[cn[:, 0]]
            V2 = phase.interpolate_data(
                propname=throat_valence+i)[cn[:, 1]]
        g1[throats] += F**2 * V1**2 * (D1*c1) / (R * T1)
        g2[throats] += F**2 * V2**2 * (D2*c2) / (R * T2)
        gt[throats] += F**2 * Vt**2 * (Dt*ct) / (R * Tt)
    # Preallocating g_inv
    g_inv1, g_inv2, g_invt = _np.zeros((3, len(cn)))
    f1, f2, ft = [gi != 0 for gi in [g1, g2, gt]]
    g_inv1[~f1] = g_inv2[~f2] = g_invt[~ft] = _np.inf
    g_inv1[f1] = 1/g1[f1]
    g_inv2[f2] = 1/g2[f2]
    g_invt[ft] = 1/gt[ft]
    # Apply shape factors and calculate the final conductance
    g = g_inv1/F_integral[f"{size_factors}.pore1"][throats] + g_inv2/F_integral[f"{size_factors}.pore2"][throats] + g_invt/F_integral[f"{size_factors}.throat"][throats]
    g[g != 0] = g[g != 0]**(-1)
    return g


def poisson(target,
            pore_area='pore.area',
            throat_area='throat.area',
            conduit_lengths='throat.conduit_lengths',
            conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    Calculate the ionic conductance of conduits in network (using the Poisson
    equation for charge conservation), where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

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

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _np.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _np.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Poisson or Laplace
    epsilon0 = 8.854187817e-12
    epsilonr = phase['pore.permittivity'][0]
    g1[m1] = epsilon0 * epsilonr * (A1)[m1] / L1[m1]
    g2[m2] = epsilon0 * epsilonr * (A2)[m2] / L2[m2]
    gt[mt] = epsilon0 * epsilonr * (At)[mt] / Lt[mt]
    # Preallocating g_inv
    g_inv1, g_inv2, g_invt = _np.zeros((3, len(Lt)))
    f1, f2, ft = [gi != 0 for gi in [g1, g2, gt]]
    g_inv1[~f1] = g_inv2[~f2] = g_invt[~ft] = _np.inf
    g_inv1[f1] = 1/g1[f1]
    g_inv2[f2] = 1/g2[f2]
    g_invt[ft] = 1/gt[ft]
    # Apply shape factors and calculate the final conductance
    g = g_inv1/SF1 + g_inv2/SF2 + g_invt/SFt
    g[g != 0] = g[g != 0]**(-1)
    return g


def electroneutrality(target,
                      pore_area='pore.area',
                      throat_area='throat.area',
                      pore_diffusivity='pore.diffusivity',
                      throat_diffusivity='throat.diffusivity',
                      conduit_lengths='throat.conduit_lengths',
                      conduit_shape_factors='throat.poisson_shape_factors',
                      pore_volume='pore.volume',
                      pore_temperature='pore.temperature',
                      throat_temperature='throat.temperature',
                      pore_valence='pore.valence',
                      throat_valence='throat.valence',
                      pore_concentration='pore.concentration',
                      ions=[]):
    r"""
    Calculate the ionic conductance of conduits in network (assuming
    electroneutrality for charge conservation), where a conduit is
    ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

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

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    pore_volume : string
        Dictionary key of the pore volume values

    pore_temperature : string
        Dictionary key of the pore temperature values

    throat_temperature : string
        Dictionary key of the throat temperature values

    pore_valence : string
       Dictionary key of the pore ionic species valence values

    throat_valence : string
       Dictionary key of the throat ionic species valence values

    pore_concentration : string
       Dictionary key of the pore ionic species concentration values

    Returns
    -------
    g : ndarray
        Array containing ionic conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Getting conduit lengths
    L1 = network[conduit_lengths + '.pore1'][throats]
    Lt = network[conduit_lengths + '.throat'][throats]
    L2 = network[conduit_lengths + '.pore2'][throats]
    # Preallocating g
    g1, g2, gt = _np.zeros((3, len(Lt)))
    # Setting g to inf when Li = 0 (ex. boundary pores)
    # INFO: This is needed since area could also be zero, which confuses NumPy
    m1, m2, mt = [Li != 0 for Li in [L1, L2, Lt]]
    g1[~m1] = g2[~m2] = gt[~mt] = _np.inf
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors+'.pore1'][throats]
        SFt = phase[conduit_shape_factors+'.throat'][throats]
        SF2 = phase[conduit_shape_factors+'.pore2'][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Electroneutrality
    F = 96485.3329
    R = 8.3145
    # Getting pores volumes
    Vol1 = network[pore_volume][cn[:, 0]]
    Vol2 = network[pore_volume][cn[:, 1]]
    Tt = phase[throat_temperature][throats]
    try:
        T1 = phase[pore_temperature][cn[:, 0]]
        T2 = phase[pore_temperature][cn[:, 1]]
    except KeyError:
        T1 = phase.interpolate_data(propname=throat_temperature)[cn[:, 0]]
        T2 = phase.interpolate_data(propname=throat_temperature)[cn[:, 1]]
    # Iterate over all ions present in the solution
    if ions == []:
        logger.error('List of ions must be provided')
    for i in ions:
        i = '.'+i
        # Check if a concetration field is defined
        try:
            c1 = phase[pore_concentration+i][cn[:, 0]]
            c2 = phase[pore_concentration+i][cn[:, 1]]
        except KeyError:
            c1 = _np.zeros((network.Np))[cn[:, 0]]
            c2 = _np.zeros((network.Np))[cn[:, 1]]
        ct = (c1*Vol1 + c2*Vol2)/(Vol1 + Vol2)
        Dt = phase[throat_diffusivity+i][throats]
        Vt = phase[throat_valence+i][throats]
        try:
            D1 = phase[pore_diffusivity+i][cn[:, 0]]
            D2 = phase[pore_diffusivity+i][cn[:, 1]]
            V1 = phase[pore_valence+i][cn[:, 0]]
            V2 = phase[pore_valence+i][cn[:, 1]]
        except KeyError:
            D1 = phase.interpolate_data(
                propname=throat_diffusivity+i)[cn[:, 0]]
            D2 = phase.interpolate_data(
                propname=throat_diffusivity+i)[cn[:, 1]]
            V1 = phase.interpolate_data(
                propname=throat_valence+i)[cn[:, 0]]
            V2 = phase.interpolate_data(
                propname=throat_valence+i)[cn[:, 1]]

        g1[m1] += F**2 * V1**2 * (D1*A1*c1)[m1] / (R * T1 * L1[m1])
        g2[m2] += F**2 * V2**2 * (D2*A2*c2)[m1] / (R * T2 * L2[m2])
        gt[mt] += F**2 * Vt**2 * (Dt*At*ct)[mt] / (R * Tt * Lt[mt])
    # Preallocating g_inv
    g_inv1, g_inv2, g_invt = _np.zeros((3, len(Lt)))
    f1, f2, ft = [gi != 0 for gi in [g1, g2, gt]]
    g_inv1[~f1] = g_inv2[~f2] = g_invt[~ft] = _np.inf
    g_inv1[f1] = 1/g1[f1]
    g_inv2[f2] = 1/g2[f2]
    g_invt[ft] = 1/gt[ft]
    # Apply shape factors and calculate the final conductance
    g = g_inv1/SF1 + g_inv2/SF2 + g_invt/SFt
    g[g != 0] = g[g != 0]**(-1)
    return g
