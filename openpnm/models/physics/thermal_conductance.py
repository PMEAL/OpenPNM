from .misc import poisson_conductance


def series_resistors(target,
                     pore_thermal_conductivity='pore.thermal_conductivity',
                     throat_thermal_conductivity='throat.thermal_conductivity',
                     throat_equivalent_area='throat.equivalent_area',
                     throat_conduit_lengths='throat.conduit_lengths'):
    r"""
    Calculate the thermal conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas

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
                               pore_diffusivity=pore_thermal_conductivity,
                               throat_diffusivity=throat_thermal_conductivity,
                               throat_equivalent_area=throat_equivalent_area,
                               throat_conduit_lengths=throat_conduit_lengths)
