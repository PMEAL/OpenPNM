r"""
Pore-scale models for calculating the conductance of conduits.
"""
__all__ = ["_poisson_conductance"]


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
    F = network[size_factors]
    # TODO: Uncomment the following 2 lines once #2087 is merged
    # Dt = phase[throat_conductivity][throats]
    # D1, D2 = phase[pore_conductivity][cn].T
    # TODO: Delete the following 2 lines once #2087 is merged
    Dt = _parse_input(phase, throat_conductivity)[throats]
    D1, D2 = _parse_input(phase, pore_conductivity)[cn].T
    # If individual size factors for conduit constiuents are known
    if isinstance(F, dict):
        g1 = D1 * F[f"{size_factors}.pore1"][throats]
        gt = Dt * F[f"{size_factors}.throat"][throats]
        g2 = D2 * F[f"{size_factors}.pore2"][throats]
        return 1 / (1 / g1 + 1 / gt + 1 / g2)
    # Otherwise, i.e., the size factor for the entire conduit is only known
    return Dt * F[throats]


def _parse_input(obj, arg):
    return obj[arg] if isinstance(arg, str) else arg
