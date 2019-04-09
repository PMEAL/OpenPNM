r"""

.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion

"""

import scipy as _sp
from .misc import generic_conductance


def ordinary_diffusion(target,
                       pore_area='pore.area',
                       throat_area='throat.area',
                       pore_diffusivity='pore.diffusivity',
                       throat_diffusivity='throat.diffusivity',
                       conduit_lengths='throat.conduit_lengths',
                       conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    Returns
    -------
    g : ndarray
        Array containing diffusive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

    """
    return generic_conductance(target=target, transport_type='diffusion',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_diffusivity,
                               throat_diffusivity=throat_diffusivity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)


def taylor_aris_diffusion(
        target,
        pore_area='pore.area',
        throat_area='throat.area',
        pore_diffusivity='pore.diffusivity',
        pore_pressure='pore.pressure',
        throat_hydraulic_conductance='throat.hydraulic_conductance',
        throat_diffusivity='throat.diffusivity',
        conduit_lengths='throat.conduit_lengths',
        conduit_shape_factors='throat.poisson_shape_factors'):
    r"""
    """
    return generic_conductance(
        target=target,
        transport_type='taylor_aris_diffusion',
        pore_area=pore_area,
        throat_area=throat_area,
        pore_diffusivity=pore_diffusivity,
        throat_diffusivity=throat_diffusivity,
        conduit_lengths=conduit_lengths,
        conduit_shape_factors=conduit_shape_factors,
        pore_pressure=pore_pressure,
        throat_hydraulic_conductance=throat_hydraulic_conductance)


def classic_ordinary_diffusion(target,
                               pore_molar_density='pore.molar_density',
                               pore_diffusivity='pore.diffusivity',
                               pore_area='pore.area',
                               pore_diameter='pore.diameter',
                               throat_area='throat.area',
                               throat_length='throat.length',
                               throat_diameter='throat.diameter',
                               shape_factor='throat.shape_factor',
                               **kwargs):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas
    Parameters
    ----------
    network : OpenPNM Network Object
    phase : OpenPNM Phase Object
        The phase of interest
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
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    parea = network[pore_area]
    pdia = network[pore_diameter]
    # Get the properties of every throat
    tdia = network[throat_diameter]
    tarea = _sp.pi * (tdia / 2) ** 2
    tlen = network[throat_length]
    # Interpolate pore phase property values to throats
    DABt = phase.interpolate_data(propname=pore_diffusivity)[throats]
    ct = phase.interpolate_data(propname=pore_molar_density)[throats]
    # Get pore lengths
    plen1 = (0.5 * pdia[Ps[:, 0]])
    plen2 = (0.5 * pdia[Ps[:, 1]])
    # Remove any non-positive lengths
    plen1[plen1 <= 0] = 1e-12
    plen2[plen2 <= 0] = 1e-12
    # Find g for half of pore 1
    gp1 = ct * DABt * parea[Ps[:, 0]] / plen1
    gp1[_sp.isnan(gp1)] = _sp.inf
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for half of pore 2
    gp2 = ct * DABt * parea[Ps[:, 1]] / plen2
    gp2[_sp.isnan(gp2)] = _sp.inf
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat, remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except:
        sf = _sp.ones(network.num_throats())
    sf[_sp.isnan(sf)] = 1.0
    gt = (1 / sf) * ct * DABt * tarea / tlen
    # Set 0 conductance pores (boundaries) to inf
    gt[~(gt > 0)] = _sp.inf
    value = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
    return value
