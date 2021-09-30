r"""
Pore-scale models for calculating the conductance of conduits.
"""
__all__ = ["generic_transport_conductance"]


def generic_transport_conductance(target,
                                  pore_conductivity=None,
                                  throat_conductivity=None,
                                  size_factors=None):
    r"""
    Calculate the conductance of conduits in network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_conductivity : string
        Dictionary key of the pore conductivity values

    throat_conductivity : string
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
    This function requires that all the necessary phase properties already
    be calculated.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network['throat.conns'][throats]
    F = network[size_factors]
    Dt = phase[throat_conductivity][throats]
    try:
        D1 = phase[pore_conductivity][cn[:, 0]]
        D2 = phase[pore_conductivity][cn[:, 1]]
    except KeyError:
        D1 = phase.interpolate_data(propname=throat_conductivity)[cn[:, 0]]
        D2 = phase.interpolate_data(propname=throat_conductivity)[cn[:, 1]]
    if isinstance(F, dict):
        g1 = D1 * F[f"{size_factors}.pore1"][throats]
        gt = Dt * F[f"{size_factors}.throat"][throats]
        g2 = D2 * F[f"{size_factors}.pore2"][throats]
        gtotal = 1 / (1 / g1 + 1 / gt + 1 / g2)
    else:
        gtotal = Dt * F[throats]
    return gtotal
