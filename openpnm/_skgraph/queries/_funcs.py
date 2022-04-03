import numpy as np


__all__ = [
    'find_interface_bonds',
    'filter_by_z',
]


def find_interface_bonds(conns, P1, P2):
    """
    Finds all bonds between two sets pores with the given labels.

    Parameters
    ----------
    network : GenericNetwork
        The network on which the query is to be performed
    P1 : array_like
        The first set of pores
    P2 : array_like
        The second set of pores

    Returns
    -------
    ndarray
        List of interface throats between the two given sets of pores

    Notes
    -----
    This method requires that the two labels do not share any pores.

    """
    if np.intersect1d(P1, P2).size != 0:
        raise Exception("P1 and P2 must not share any pores.")
    Ts1 = find_neighbor_bonds(P1, mode="xor")
    Ts2 = find_neighbor_bonds(P2, mode="xor")
    return np.intersect1d(Ts1, Ts2)


def filter_by_z(conns, sites, z=1):
    r"""
    Find sites with a given number of neighbors

    Parameters
    ----------
    conns : ndarray
        The network on which the query is to be performed
    sites : array_like
        The sites to be filtered
    z : int
        The coordination number to filter by

    Returns
    -------
    pores : array_like
        A list of pores which satisfy the criteria

    """
    pores = network._parse_indices(pores)
    Nz = network.num_neighbors(pores=pores)
    orphans = np.where(Nz == z)[0]
    hits = pores[orphans]
    return hits
