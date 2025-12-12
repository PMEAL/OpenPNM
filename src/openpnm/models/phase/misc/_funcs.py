import numpy as np
from openpnm.models.phase import _phasedocs


__all__ = [
    "mix_and_match",
    "mole_to_mass_fraction",
]


@_phasedocs
def mix_and_match(
    phase,
    prop,
    phases,
    occupancy,
):
    r"""
    Return the given property by looking it up from a list of given phases
    based on occupancy.

    Parameters
    ----------
    %(phase)s
    prop : str
        The dictionary key to the array containing the pore/throat property to
        be used in the calculation.
    phases : list
        List of Phases over which the given `prop` is to be
        averaged out.
    occupancy : str
        The dictionary key to the array containing the occupancy associated
        with each of the given ``phases``.

    Returns
    -------
    weighted_average : ndarray
        Weighted average of the given `prop` averaged over `phases`.

    """
    # Hack for ModelsMixin to not complain (cyclic dep)
    prop = prop.strip("_")
    values = np.zeros_like(phases[0][prop])

    for phase in phases:
        mask = phase[occupancy]
        values[mask] = phase[prop][mask]
    return values


@_phasedocs
def mole_to_mass_fraction(
    phase,
    MWs='param.molecular_weight.*',
):
    r"""
    Convert mole fraction to mass fraction

    Parameters
    ----------
    %(phase)s
    %(MWs)s

    Returns
    -------

    """
    MWs = phase.get_comp_vals(MWs)
    xs = phase['pore.mole_fraction']
    ms = {}
    # Find the actual masses in each pore
    for c in xs.keys():
        ms[c] = xs[c]*MWs[c]
    # Normalize component mass by total mass in each pore
    denom = np.vstack(list(ms.values())).sum(axis=0)
    for c in xs.keys():
        ms[c] = ms[c]/denom
    return ms
