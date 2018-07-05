from .misc import poisson_conductance


def ordinary_diffusion(target,
                       pore_diffusivity='pore.diffusivity',
                       throat_diffusivity='throat.diffusivity',
                       throat_equivalent_area='throat.equivalent_area',
                       throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the diffusive conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

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

    throat_equivalent_area : string
        Dictionary key of the throat equivalent area values

    throat_conduit_lengths : string
        Dictionary key of the throat conduit lengths

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    """
    return poisson_conductance(target=target,
                               pore_diffusivity=pore_diffusivity,
                               throat_diffusivity=throat_diffusivity,
                               throat_equivalent_area=throat_equivalent_area,
                               throat_conduit_lengths=throat_conduit_lengths)
