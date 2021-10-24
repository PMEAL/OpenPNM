r"""
Pore-scale models for calculating the diffusive conductance of conduits.
"""
import numpy as _np
import scipy.constants as _const
from openpnm.models.physics.utils import _poisson_conductance

__all__ = [
    "generic_diffusive",
    "ordinary_diffusion",
    "mixed_diffusion",
    "taylor_aris_diffusion",
    "multiphase_diffusion"
]


def generic_diffusive(target,
                      pore_diffusivity="pore.diffusivity",
                      throat_diffusivity="throat.diffusivity",
                      size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

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
    return _poisson_conductance(target=target,
                                pore_conductivity=pore_diffusivity,
                                throat_conductivity=throat_diffusivity,
                                size_factors=size_factors)


def ordinary_diffusion(target,
                       pore_diffusivity="pore.diffusivity",
                       throat_diffusivity="throat.diffusivity",
                       size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore).

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
    return _poisson_conductance(target=target,
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
                    molecular_weight="pore.molecular_weight",
                    size_factors="throat.diffusive_size_factors"):
    r"""
    Calculates the diffusive conductance of conduits in network with
    Knudsen correction.

    A conduit is defined as (1/2 pore - full throat - 1/2 pore). See
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
    molecular_weigth : str
        Dictionary key of the pore molecular weight values.
    size_factors: str
        Dictionary key of the conduit diffusive shape factors' values.

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
    # Fetch openpnm objects
    network = target.network
    phase = target.project.find_phase(target)

    # Fetch model parameters
    Dp = phase[pore_diffusivity]
    Dt = phase[throat_diffusivity]
    dp = network[pore_diameter]
    dt = network[throat_diameter]
    MWp = phase[molecular_weight]
    MWt = phase.interpolate_data(propname=molecular_weight)
    Tp = phase[pore_temperature]
    Tt = phase[throat_temperature]

    # Calculate Knudsen contribution
    DKp = dp/3 * (8*_const.R*Tp / (_const.pi*MWp)) ** 0.5
    DKt = dt/3 * (8*_const.R*Tt / (_const.pi*MWt)) ** 0.5

    # Calculate mixed diffusivity
    Dp_eff = (1/DKp + 1/Dp) ** -1
    Dt_eff = (1/DKt + 1/Dt) ** -1

    return _poisson_conductance(target=target,
                                pore_conductivity=Dp_eff,
                                throat_conductivity=Dt_eff,
                                size_factors=size_factors)


def taylor_aris_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diffusivity="pore.diffusivity",
    pore_pressure="pore.pressure",
    throat_hydraulic_conductance="throat.hydraulic_conductance",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
):
    r"""
    Calculates the diffusive conductance of conduits in network
    considering the Taylor-Aris effect (effect of flow on diffusion).

    A conduit is defined as (1/2 pore - full throat - 1/2 pore)

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
    pore_pressure : str
        Dictionary key of the pore pressure values.
    throat_hydraulic_conductance : str
        Dictionary key of the throat hydraulic_conductance values.
    throat_diffusivity : str
        Dictionary key of the throat diffusivity values.
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

    Returns
    -------
    ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties are
    already calculated.

    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

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
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
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
    # Find g for half of pore 1, throat, and half of pore 2
    g1 = D1 * (1 + (Pe1 ** 2) / 192) * A1 / L1
    g2 = D2 * (1 + (Pe2 ** 2) / 192) * A2 / L2
    gt = Dt * (1 + (Pet ** 2) / 192) * At / Lt
    # Ensure infinite conductance for elements with zero length
    g1[L1 == 0] = _np.inf
    g2[L2 == 0] = _np.inf
    gt[Lt == 0] = _np.inf
    # Apply shape factors and calculate the final conductance
    return (1 / gt / SFt + 1 / g1 / SF1 + 1 / g2 / SF2) ** (-1)


def multiphase_diffusion(
    target,
    pore_area="pore.area",
    throat_area="throat.area",
    pore_diffusivity="pore.diffusivity",
    throat_diffusivity="throat.diffusivity",
    conduit_lengths="throat.conduit_lengths",
    conduit_shape_factors="throat.poisson_shape_factors",
    partition_coef_global="throat.partition_coef.all"
):
    r"""
    Calculates the diffusive conductance of conduits in network.

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
    conduit_lengths : str
        Dictionary key of the conduit length values.
    conduit_shape_factors : str
        Dictionary key of the conduit diffusive shape factors' values.

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
    # Getting equivalent areas
    A1, A2 = network[pore_area][cn].T
    At = network[throat_area][throats]
    # Getting conduit lengths
    L1 = network[f"{conduit_lengths}.pore1"][throats]
    Lt = network[f"{conduit_lengths}.throat"][throats]
    L2 = network[f"{conduit_lengths}.pore2"][throats]
    # Getting shape factors
    try:
        SF1 = phase[conduit_shape_factors + ".pore1"][throats]
        SFt = phase[conduit_shape_factors + ".throat"][throats]
        SF2 = phase[conduit_shape_factors + ".pore2"][throats]
    except KeyError:
        SF1 = SF2 = SFt = 1.0
    # Interpolate pore phase property values to throats
    D1, D2 = phase[pore_diffusivity][cn].T
    Dt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    # Find g for half of pore 1, throat, and half of pore 2 + apply shape factors
    g1 = (D1 * A1) / L1 * SF1
    g2 = (D2 * A2) / L2 * SF2
    gt = (Dt * At) / Lt * SFt
    # Apply Henry's partitioning coefficient
    # Note 1: m12 = G21*c1 - G12*c2
    K12 = phase[partition_coef_global][throats]
    G21 = (1/g1 + 0.5/gt + K12 * (1/g2 + 0.5/gt)) ** (-1)
    G12 = K12 * G21
    return _np.vstack((G12, G21)).T
