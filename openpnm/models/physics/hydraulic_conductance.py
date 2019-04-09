r"""

.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille

"""

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
