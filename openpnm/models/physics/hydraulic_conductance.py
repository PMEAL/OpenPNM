r"""

.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille

"""

import scipy as _sp
from .misc import generic_conductance


def hagen_poiseuille(target,
                     pore_area='pore.area',
                     throat_area='throat.area',
                     pore_viscosity='pore.viscosity',
                     throat_viscosity='throat.viscosity',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.flow_shape_factors'):
    r"""
    Calculate the hydraulic conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_viscosity : string
        Dictionary key of the pore viscosity values

    throat_viscosity : string
        Dictionary key of the throat viscosity values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_shape_factors : string
        Dictionary key of the conduit FLOW shape factor values

    Returns
    -------
    g : ndarray
        Array containing hydraulic conductance values for conduits in the
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
    return generic_conductance(target=target, transport_type='flow',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)


def hagen_poiseuille_power_law(
        target,
        pore_area='pore.area',
        throat_area='throat.area',
        pore_viscosity='pore.viscosity',
        throat_viscosity='throat.viscosity',
        conduit_lengths='throat.conduit_lengths',
        conduit_shape_factors='throat.flow_shape_factors',
        pore_consistency='pore.consistency',
        throat_consistency='throat.consistency',
        pore_flow_index='pore.flow_index',
        throat_flow_index='throat.flow_index',
        pore_pressure='pore.pressure'):
    r"""
    """
    return generic_conductance(target=target, transport_type='flow_power_law',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity=pore_viscosity,
                               throat_diffusivity=throat_viscosity,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors,
                               pore_consistency=pore_consistency,
                               throat_consistency=throat_consistency,
                               pore_flow_index=pore_flow_index,
                               throat_flow_index=throat_flow_index,
                               pore_pressure=pore_pressure)


def classic_hagen_poiseuille(target,
                             pore_diameter='pore.diameter',
                             pore_viscosity='pore.viscosity',
                             throat_length='throat.length',
                             throat_diameter='throat.diameter',
                             shape_factor='throat.shape_factor',
                             **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming cylindrical
    geometry using the Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : OpenPNM Phase Object

    Notes
    -----
    This function calculates the specified property for the *entire* network
    then extracts the values for the appropriate throats at the end.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    # Get Nt-by-2 list of pores connected to each throat
    Ps = network['throat.conns']
    # Get properties in every pore in the network
    phase = target.project.find_phase(target)
    mup = phase[pore_viscosity]
    mut = phase.interpolate_data(propname=pore_viscosity)[throats]
    pdia = network[pore_diameter]
    # Get pore lengths
    plen1 = (0.5*pdia[Ps[:, 0]])
    plen2 = (0.5*pdia[Ps[:, 1]])
    # Remove any non-positive lengths
    plen1[plen1 <= 0] = 1e-12
    plen2[plen2 <= 0] = 1e-12
    # Find g for half of pore 1
    gp1 = _sp.pi*(pdia[Ps[:, 0]])**4/(128*plen1*mut)
    gp1[_sp.isnan(gp1)] = _sp.inf
    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf

    # Find g for half of pore 2
    gp2 = _sp.pi*(pdia[Ps[:, 1]])**4/(128*plen2*mut)
    gp2[_sp.isnan(gp2)] = _sp.inf
    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    # Find g for full throat
    tdia = network[throat_diameter]
    tlen = network[throat_length]
    # Remove any non-positive lengths
    tlen[tlen <= 0] = 1e-12
    # Get shape factor
    try:
        sf = network[shape_factor]
    except:
        sf = _sp.ones(network.num_throats())
    sf[_sp.isnan(sf)] = 1.0
    gt = (1/sf)*_sp.pi*(tdia)**4/(128*tlen*mut)
    gt[~(gt > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return value
