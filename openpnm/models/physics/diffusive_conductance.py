r"""
Pore-scale models for calculating the diffusive conductance of conduits.
"""
import numpy as _np
import scipy.constants as _const
from openpnm.models.physics.utils import generic_transport_conductance

__all__ = [
    "generic_diffusive",
    "ordinary_diffusion",
    "mixed_diffusion",
    "taylor_aris_diffusion",
    "classic_ordinary_diffusion",
    "multiphase_diffusion"
]


def generic_diffusive(target,
                      pore_diffusivity="pore.diffusivity",
                      throat_diffusivity="throat.diffusivity",
                      size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    size_factors: str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    """
    return generic_transport_conductance(target=target,
                                         pore_conductivity=pore_diffusivity,
                                         throat_conductivity=throat_diffusivity,
                                         size_factors=size_factors)


def ordinary_diffusion(target,
                       pore_diffusivity="pore.diffusivity",
                       throat_diffusivity="throat.diffusivity",
                       size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ).

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    size_factors: str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    """
    return generic_transport_conductance(target=target,
                                         pore_conductivity=pore_diffusivity,
                                         throat_conductivity=throat_diffusivity,
                                         size_factors=size_factors)


def mixed_diffusion(target,
                    pore_diameter="pore.diameter",
                    throat_diameter="throat.diameter",
                    pore_diffusivity="pore.diffusivity",
                    throat_diffusivity="throat.diffusivity",
                    pore_temperature="pore.temperature",
                    throat_temperature="throat.temperature",
                    pore_molecular_weight="pore.molecular_weight",
                    throat_molecular_weight="throat.molecular_weight",
                    size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network with
    Knudsen correction.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore ). See
    Notes section for the limitations of this method.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diameter : str
        Dictionary key of the pore diameter values.
    throat_diameter : str
        Dictionary key of the throat diameter values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    pore_temperature : str
        Dictionary key of the pore temperature values.
    throat_temperature : str
        Dictionary key of the throat temperature values.
    pore_molecular_weigth : str
        Dictionary key of the pore molecular weight values.
    throat_molecular_weigth : str
        Dictionary key of the throat molecular weight values.
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function is ONLY suitable for dilute mixtures and NOT those with
    concentrated species.

    This function requires that all the necessary phase properties already
    be calculated.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F = network[size_factors]
    Dt = phase[throat_diffusivity][throats]
    D1, D2 = phase[pore_diffusivity][cn].T
    # Calculating Knudsen diffusivity
    d1, d2 = network[pore_diameter][cn].T
    dt = network[throat_diameter][throats]
    MW1, MW2 = phase[pore_molecular_weight][cn].T
    MWt = phase[throat_molecular_weight][throats]
    T1, T2 = phase[pore_temperature][cn].T
    Tt = phase[throat_temperature][throats]
    DK1 = d1 / 3 * (8 * _const.R * T1 / _const.pi / MW1) ** 0.5
    DK2 = d2 / 3 * (8 * _const.R * T2 / _const.pi / MW2) ** 0.5
    DKt = dt / 3 * (8 * _const.R * Tt / _const.pi / MWt) ** 0.5
    # Calculate mixed diffusivity
    D1e = (1 / DK1 + 1 / D1) ** (-1)
    D2e = (1 / DK2 + 1 / D2) ** (-1)
    Dte = (1 / DKt + 1 / Dt) ** (-1)
    if isinstance(F, dict):
        g1 = D1e * F[f"{size_factors}.pore1"][throats]
        gt = Dte * F[f"{size_factors}.throat"][throats]
        g2 = D2e * F[f"{size_factors}.pore2"][throats]
        gtotal = 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        gtotal = Dte * F[throats]
    return gtotal


def taylor_aris_diffusion(target,
                          pore_area="pore.area",
                          throat_area="throat.area",
                          pore_diffusivity="pore.diffusivity",
                          throat_diffusivity="throat.diffusivity",
                          pore_pressure="pore.pressure",
                          throat_hydraulic_conductance="throat.hydraulic_conductance",
                          size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network
    considering the Taylor-Aris effect (effect of flow on diffusion).

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore )

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_area : str
        Dictionary key of the pore area values.
    throat_area : str
        Dictionary key of the throat area values.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    pore_pressure : str
        Dictionary key of the pore pressure values.
    throat_hydraulic_conductance : str
        Dictionary key of the throat hydraulic_conductance values.
    size_factors: str
        Dictionary key of the conduit size factors' values.
    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties are
    already calculated.
    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F = network[size_factors]
    # Getting equivalent areas
    A1 = network[pore_area][cn[:, 0]]
    At = network[throat_area][throats]
    A2 = network[pore_area][cn[:, 1]]
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity][throats]
    # Fetch properties for calculating Peclet
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]
    Qt = -gh * _np.diff(P[cn], axis=1).squeeze()
    # Find fluid velocity in elements
    u1 = Qt / A1
    u2 = Qt / A2
    ut = Qt / At
    # Find Peclet number in elements
    Pe1 = u1 * ((4 * A1 / _np.pi) ** 0.5) / D1
    Pe2 = u2 * ((4 * A2 / _np.pi) ** 0.5) / D2
    Pet = ut * ((4 * At / _np.pi) ** 0.5) / Dt
    if isinstance(F, dict):
        g1 = D1 * (1 + (Pe1 ** 2) / 192) * F[f"{size_factors}.pore1"][throats]
        gt = Dt * (1 + (Pet ** 2) / 192) * F[f"{size_factors}.throat"][throats]
        g2 = D2 * (1 + (Pe2 ** 2) / 192) * F[f"{size_factors}.pore2"][throats]
        gtotal = 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        gtotal = Dt * (1 + (Pet ** 2) / 192) * F[throats]
    return gtotal


def classic_ordinary_diffusion(target,
                               pore_molar_density="pore.molar_density",
                               throat_molar_density="throat.molar_density",
                               pore_diffusivity="pore.diffusivity",
                               throat_diffusivity="throat.diffusivity",
                               pore_diameter="pore.diameter",
                               pore_area="pore.area",
                               size_factors="throat.diffusive_size_factors"):
    r"""
    Legacy method that calculates the diffusive conductance of conduits.

    A conduit is defined as ( 1/2 pore - full throat - 1/2 pore )

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_molar_density : str
        Dictionary key pointing to pore molar density values.
    throat_molar_density : str
        Dictionary key pointing to throat molar density values.
    pore_diffusivity : str
        Dictionary key pointing to pore diffusivity values.
    throat_diffusivity : str
        Dictionary key pointing to throat diffusivity values.
    pore_diameter : str
        Dictionary key pointing to pore diameter values.
    pore_area : str
        Dictionary key pointing to pore area values.
    size_factors: str
        Dictionary key of the conduit diffusive size factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This method requires that all the necessary phase properties already
    be calculated.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F = network[size_factors]
    # Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Get the properties of every throat
    DABt = phase[throat_diffusivity][throats]
    ct = phase[throat_molar_density][throats]
    # Get pore lengths
    plen1 = 0.5 * pdia[cn[:, 0]]
    plen2 = 0.5 * pdia[cn[:, 1]]
    # Remove any non-positive lengths
    plen1[plen1 <= 1e-12] = 1e-12
    plen2[plen2 <= 1e-12] = 1e-12
    # Find g for half of pore 1
    gp1 = ct * DABt * parea[cn[:, 0]] / plen1
    gp1[_np.isnan(gp1)] = _np.inf
    gp1[~(gp1 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = ct * DABt * parea[cn[:, 1]] / plen2
    gp2[_np.isnan(gp2)] = _np.inf
    gp2[~(gp2 > 0)] = _np.inf  # Set 0 conductance pores (boundaries) to inf
    if isinstance(F, dict):
        gt = ct * DABt * F[f"{size_factors}.throat"][throats]
        # Set 0 conductance pores (boundaries) to inf
        gt[~(gt > 0)] = _np.inf
        gtotal = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
    else:
        gtotal = ct * DABt * F[throats]
    return gtotal


def multiphase_diffusion(target,
                         pore_diffusivity="pore.diffusivity",
                         throat_diffusivity="throat.diffusivity",
                         size_factors="throat.diffusive_size_factors",
                         partition_coef_global="throat.partition_coef.all"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    Parameters
    ----------
    target : GenericPhysics
        Physics object with which this model is associated.
    pore_diffusivity : str
        Dictionary key of the pore diffusivity values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    ndarray
        Array (Nt by 2) containing diffusive conductance values for
        conduits in the geometry attached to the given physics object.

    Notes
    -----
    This method assumes that phase["partition_coef"] contains information on
    binary phase partitioning. See ``MultiPhase`` class documentation for more
    information.

    """
    network = target.network
    throats = target.throats(target=network)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]
    F = network[size_factors]
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase[throat_diffusivity][throats]

    # Find g for half of pore 1, throat, and half of pore 2 + apply size factors
    g1 = D1 * F[f"{size_factors}.pore1"][throats]
    g2 = D2 * F[f"{size_factors}.pore2"][throats]
    gt = Dt * F[f"{size_factors}.throat"][throats]

    # Apply Henry's partitioning coefficient
    # Note 1: m12 = G21*c1 - G12*c2
    K12 = phase[partition_coef_global][throats]
    G21 = (1/g1 + 0.5/gt + K12 * (1/g2 + 0.5/gt)) ** (-1)
    G12 = K12 * G21

    return _np.vstack((G12, G21)).T
