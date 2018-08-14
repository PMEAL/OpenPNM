from .misc import generic_conductance


def ordinary_diffusion(target,
                       pore_diffusivity='pore.diffusivity',
                       throat_diffusivity='throat.diffusivity',
                       pore_area='pore.area',
                       throat_area='throat.area',
                       conduit_lengths='throat.conduit_lengths',
                       conduit_shape_factors='throat.diffusion_shape_factor'):
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
                               pore_diffusivity=pore_diffusivity,
                               throat_diffusivity=throat_diffusivity,
                               pore_area=pore_area,
                               throat_area=throat_area,
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors=conduit_shape_factors)
