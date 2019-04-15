import numpy as _np


def mix_and_match(target, prop, phases, occupancy):
    r"""
    Return the given property by looking it up from a list of given phases
    based on occupancy.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    prop : string
        The dictionary key to the array containing the pore/throat property to
        be used in the calculation.

    phases : list
        List of OpenPNM phase objects over which the given `prop` is to be
        averaged out.

    occupancy : string
        The dictionary key to the array containing the occupancy associated
        with each of the given ``phases``.

    Returns
    -------
    weighted_average : ndarray
        Weighted average of the given `prop` averaged over `phases`.

    """
    # Hack for ModelsMixin to not complain (cyclic dep)
    prop = prop.strip("_")

    values = _np.zeros_like(phases[0][prop])

    for phase in phases:
        mask = phase[occupancy]
        values[mask] = phase[prop][mask]

    return values
