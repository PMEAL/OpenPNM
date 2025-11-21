r"""
Pore-scale models for calculating the conductance of conduits.
"""
from numpy import vstack


__all__ = [
    '_poisson_conductance',
    '_get_key_props',
]


def _poisson_conductance(phase,
                         pore_conductivity=None,
                         throat_conductivity=None,
                         size_factors=None):
    r"""
    Calculates the conductance of the conduits in the network, where a
    conduit is (1/2 pore - full throat - 1/2 pore). See the notes section.

    Parameters
    ----------
    phase : OpenPNM Phase
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
    network = phase.network
    cn = network.conns
    Dt = phase[throat_conductivity]
    D1, D2 = phase[pore_conductivity][cn].T
    # If individual size factors for conduit constiuents are known
    SF = network[size_factors]
    if SF.ndim == 2:
        F1, Ft, F2 = SF.T
        g1 = D1 * F1
        gt = Dt * Ft
        g2 = D2 * F2
        return 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        # Otherwise, i.e., the size factor for the entire conduit is only known
        F = network[size_factors]
        return Dt * F


def _get_key_props(phase=None,
                   diameter="throat.diameter",
                   surface_tension="throat.surface_tension",
                   contact_angle="throat.contact_angle"):
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
