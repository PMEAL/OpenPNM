r"""
Pore-scale models for calculating the conductance of conduits.
"""

def _poisson_conductance(target,
                         pore_conductivity=None,
                         throat_conductivity=None,
                         size_factors=None):
    r"""
    Calculates the conductance of the conduits in the network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    target : GenericPhysics
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_conductivity : str
        Dictionary key of the pore conductivity values
    throat_conductivity : str
        Dictionary key of the throat conductivity values
    size_factors: str
        Dictionary key of the conduit size factors' values.

    Returns
    -------
    g : ndarray
        Array containing conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    This function requires that all the necessary phase properties
    already be calculated.

    """
    network = target.network
    domain = target._domain
    throats = domain.throats(target.name)
    phase = target.project.find_phase(target)
    cn = network.conns[throats]
    # TODO: Uncomment the following 2 lines once #2087 is merged
    # Dt = phase[throat_conductivity][throats]
    # D1, D2 = phase[pore_conductivity][cn].T
    # TODO: Delete the following 2 lines once #2087 is merged
    Dt = _parse_input(phase, throat_conductivity)[throats]
    D1, D2 = _parse_input(phase, pore_conductivity)[cn].T
    # If individual size factors for conduit constiuents are known
    try:
        F1, Ft, F2 = network[size_factors].T
        g1 = D1 * F1[throats]
        gt = Dt * Ft[throats]
        g2 = D2 * F2[throats]
        return 1 / (1 / g1 + 1 / gt + 1 / g2)
    except ValueError:
        # Otherwise, i.e., the size factor for the entire conduit is only known
        F = network[size_factors]
        return Dt * F[throats]


def _parse_input(obj, arg):
    """
    Returns obj[arg] if arg is string, otherwise returns arg.
    """
    return obj[arg] if isinstance(arg, str) else arg


def _get_key_props(phase=None,
                   diameter="throat.diameter",
                   surface_tension="pore.surface_tension",
                   contact_angle="pore.contact_angle"):
    """
    Returns desired properties in the correct format! See Notes.

    Notes
    -----
    Many of the methods are generic to pores and throats. Some information may
    be stored on either the pore or throat and needs to be interpolated.
    This is a helper method to return the properties in the correct format.

    TODO: Check for method to convert throat to pore data

    """
    element = diameter.split(".")[0]
    if element == "pore":
        if "throat" in surface_tension:
            sigma = phase.interpolate_data(propname=surface_tension)
        else:
            sigma = phase[surface_tension]
        if "throat" in contact_angle:
            theta = phase.interpolate_data(propname=contact_angle)
        else:
            theta = phase[contact_angle]
    if element == "throat":
        if "pore" in surface_tension:
            sigma = phase.interpolate_data(propname=surface_tension)
        else:
            sigma = phase[surface_tension]
        if "pore" in contact_angle:
            theta = phase.interpolate_data(propname=contact_angle)
        else:
            theta = phase[contact_angle]

    return element, sigma, theta
