import numpy as np
from openpnm.utils import HealthDict


__all__ = [
    'check_data_health',
    'check_network_health',
]


def check_data_health(obj):
    r"""
    Checks the health of pore and throat data arrays.

    Parameters
    ----------
    obj : Base
        A handle of the object to be checked

    Returns
    -------
    health : dict
        Returns a HealthDict object which is a basic dictionary with an
        added ``health`` attribute that is ``True`` is all entries in the
        dict are deemed healthy (empty lists), or ``False`` otherwise.

    """
    health = HealthDict()
    for item in obj.props():
        health[item] = []
        if obj[item].dtype == 'O':
            health[item] = 'No checks on object'
        elif np.sum(np.isnan(obj[item])) > 0:
            health[item] = 'Has NaNs'
        elif np.shape(obj[item])[0] != obj._count(item.split('.', 1)[0]):
            health[item] = 'Wrong Length'
    return health


def check_network_health(network):
    r"""
    This method checks the topological health of the network

    The following aspects are checked for:

        (1) Isolated pores
        (2) Disconnected clusters of pores
        (3) Duplicate throats
        (4) Headless throats
        (5) Bidirectional throats

    Returns
    -------
    health : dict
        A dictionary containing the offending pores or throat numbers
        under each named key.

    Notes
    -----
    It also returns a list of which pores and throats should be trimmed
    from the network to restore health.  This list is a suggestion only,
    and is based on keeping the largest cluster and trimming the others.

    - Does not yet check for duplicate pores
    - Does not yet suggest which throats to remove
    - This is just a 'check' and does not 'fix' the problems it finds

    """
    import openpnm.models.network as mods

    health = HealthDict()

    # Check for headless throats
    headless = mods.headless_throats(network)
    health['headless_throats'] = np.where(headless)[0].tolist()

    # Check for throats that loop back onto the same pore
    looped = mods.looped_throats(network)
    health['looped_throats'] = np.where(looped)[0].tolist()

    # Check for individual isolated pores
    isolated = mods.isolated_pores(network)
    health['isolated_pores'] = np.where(isolated)[0].tolist()

    # Check for separated clusters of pores
    size = mods.cluster_size(network)
    mx = np.max(size)
    health['disconnected_pores'] = np.where(size < mx)[0].tolist()

    # Check for duplicate throats
    dupes = mods.duplicate_throats(network)
    health['duplicate_throats'] = np.where(dupes)[0].tolist()

    # Check for bidirectional throats
    bidir = mods.bidirectional_throats(network)
    health['bidirectional_throats'] = np.where(bidir)[0].tolist()

    return health
