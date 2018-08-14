from .misc import generic_conductance


def series_resistors(target,
                     pore_conductivity='pore.electrical_conductivity',
                     throat_conductivity='throat.electrical_conductivity',
                     pore_area='pore.area',
                     throat_area='throat.area',
                     conduit_lengths='throat.conduit_lengths',
                     conduit_shape_factors='throat.diffusion_shape_factor'):
    r"""
    Calculate the electrical conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_thermal_conductivity : string
        Dictionary key of the pore thermal conductivity values

    throat_thermal_conductivity : string
        Dictionary key of the throat thermal conductivity values

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

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
    return generic_conductance(target=target, mechanism='diffusion',
                               pore_diffusivity=pore_conductivity,
                               throat_diffusivity=throat_conductivity,
                               pore_area=pore_area,
                               throat_area=throat_area,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)
