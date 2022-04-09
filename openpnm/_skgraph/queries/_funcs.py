import numpy as np


__all__ = [
    'find_interface_edges',
    'filter_by_z',
]


def find_interface_edges(conns, P1, P2):
    """
    Finds all bonds between two sets of nodes

    Parameters
    ----------
    conns : ndarray
        The connections of the sparse adjacency matrix in COO format
    P1 : array_like
        The first set of nodes
    P2 : array_like
        The second set of nodes

    Returns
    -------
    ndarray
        List of interface edges between the two given sets of nodes

    """
    if np.intersect1d(P1, P2).size != 0:
        raise Exception("P1 and P2 must not share any pores.")
    Ts1 = find_neighbor_bonds(P1, mode="xor")
    Ts2 = find_neighbor_bonds(P2, mode="xor")
    return np.intersect1d(Ts1, Ts2)


def filter_by_z(conns, nodes, z=1):
    r"""
    Find nodes with a given number of neighbors

    Parameters
    ----------
    conns : ndarray
        The connections of the sparse adjacency matrix in COO format
    nodes : array_like
        The nodes to be filtered
    z : int
        The coordination number to filter by

    Returns
    -------
    pores : array_like
        A list of pores which satisfy the criteria

    """
    Nz = num_neighbors(nodes=nodes)
    orphans = np.where(Nz == z)[0]
    hits = nodes[orphans]
    return hits
