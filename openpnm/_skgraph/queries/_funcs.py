import numpy as np


__all__ = [
    'find_interface_edges',
    'filter_by_z',
    'find_connecting_edges',
    'find_neighbor_nodes',
]


def find_neighbor_nodes(nodes, am, flatten=True, include_input=False,
                        logic='or'):
    r"""
    Given a symmetric adjacency matrix, finds all nodes that are connected
    to the input nodes.

    Parameters
    ----------
    nodes : array_like
        The list of node indices for which neighbors should be found
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        nodes *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).
    flatten : bool
        If ``True`` (default) the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.
    include_input : bool
        If ``False`` (default) the input nodes will be removed from the result.

    logic : str
        Specifies logic to filter the resulting list.  Options are:

        ======= ===============================================================
        logic   Description
        ======= ===============================================================
        'or'    (default) All neighbors of the input nodes.  This is also
                known as the 'union' in set theory or 'any' in boolean logic.
                Both keywords are accepted and treated as 'or'.
        'xor'   Only neighbors of one and only one input nodes.  This is
                useful for finding the nodes that are not *shared* by any of
                the input nodes.  'exclusive_or' is also accepted.
        'xnor'  Neighbors that are shared by two or more input nodes. This
                is equivalent to finding all neighbors with 'or', minus those
                found with 'xor', and is useful for finding neighbors that the
                inputs have in common.  'nxor' is also accepted.
        'and'   Only neighbors shared by all input nodes.  This is also
                known as 'intersection' in set theory and (somtimes) as 'all'
                in boolean logic.  Both keywords are accepted and treated as
                'and'.
        ======= ===============================================================

    Returns
    -------
    nodes : ndarray
        An array containing the neighboring nodes filtered by the given logic.  If
        ``flatten`` is ``False`` then the result is a list of lists containing the
        neighbors of each input site.

    Notes
    -----
    The ``logic`` options are applied to neighboring nodes only, thus it is not
    possible to obtain nodes that are part of the global set but not neighbors.
    This is because (a) the list global nodes might be very large, and (b) it
    is not possible to return a list of neighbors for each input site if global
    nodes are considered.

    """
    if am.format != 'lil':
        am = am.tolil(copy=False)
    nodes = np.array(nodes, ndmin=1)
    if len(nodes) == 0:
        return []
    am_coo = am.tocoo()
    n_nodes = am.shape[0]
    rows = am.rows[nodes].tolist()
    if len(rows) == 0:
        return []
    neighbors = am_coo.col[np.in1d(am_coo.row, nodes)]
    if logic in ['or', 'union', 'any']:
        neighbors = np.unique(neighbors)
    elif logic in ['xor', 'exclusive_or']:
        neighbors = np.unique(np.where(np.bincount(neighbors) == 1)[0])
    elif logic in ['xnor', 'nxor']:
        neighbors = np.unique(np.where(np.bincount(neighbors) > 1)[0])
    elif logic in ['and', 'all', 'intersection']:
        neighbors = set(neighbors)
        [neighbors.intersection_update(i) for i in rows]
        neighbors = np.array(list(neighbors), dtype=np.int64, ndmin=1)
    else:
        raise Exception('Specified logic is not implemented')
    # Deal with removing inputs or not
    mask = np.zeros(shape=n_nodes, dtype=bool)
    mask[neighbors] = True
    if not include_input:
        mask[nodes] = False
    # Finally flatten or not
    if flatten:
        neighbors = np.where(mask)[0]
    else:
        if neighbors.size > 0:
            for i in range(len(rows)):
                vals = np.array(rows[i], dtype=np.int64)
                rows[i] = vals[mask[vals]]
            neighbors = rows
        else:
            neighbors = [np.array([], dtype=int) for i in range(len(nodes))]
    return neighbors


def find_connecting_edges(nodes, am):
    r"""
    Finds the edges that connects each pair of given nodes

    Parameters
    ----------
    nodes : array_like
        A 2-column vector containing pairs of node indices on each row
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        nodes *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    Returns
    -------
    edges : ndarray
        An ndarry the same length as P1 (and P2) with each element
        containing the edge number that connects the corresponding nodes,
        or `nan`` if nodes are not connected.

    Notes
    -----
    The adjacency matrix is converted to the ``DOK`` format internally if
    needed, so if this format is already available it should be provided to
    save time.

    """
    if am.format != 'dok':
        am = am.todok(copy=True)
    nodes = np.array(nodes, ndmin=2)
    if nodes.size == 0:
        return []
    z = tuple(zip(nodes[:, 0], nodes[:, 1]))
    neighbors = np.array([am.get(z[i], np.nan) for i in range(len(z))])
    return neighbors


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
