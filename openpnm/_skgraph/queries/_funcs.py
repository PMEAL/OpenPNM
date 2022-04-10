import numpy as np
import scipy.sparse as sprs
from openpnm._skgraph.tools import istriu


__all__ = [
    'find_common_edges',
    'filter_by_z',
    'find_connecting_edges',
    'find_neighbor_nodes',
    'find_neighbor_edges',
    'find_connected_nodes',
    'find_complementary_nodes',
    'find_complementary_edges',
]


def find_complementary_edges(inds, am, asmask=False):
    r"""
    Finds the complementary edges to a given set of inputs

    Parameters
    ----------
    bonds : array_like
        A list of edge indices for which the complement is sought
    am : scipy.sparse matrix
        The adjacency matrix of the network.
    asmask : bool
        If set to ``True`` the result is returned as a boolean mask of the
        correct length with ``True`` values indicate the complements.  The
        default is ``False`` which returns a list of indices instead.

    Returns
    -------
    An array containing indices of the edges that are not part of the input
    list

    """
    inds = np.unique(inds)
    N = int(am.nnz / 2)
    mask = np.ones(shape=N, dtype=bool)
    mask[inds] = False
    if asmask:
        return mask
    else:
        return np.arange(N)[mask]


def find_complementary_nodes(inds, am, asmask=False):
    r"""
    Finds the complementary nodes to a given set of inputs

    Parameters
    ----------
   inds : array_like (optional)
        A list of indices for which the complement is sought
    am : scipy.sparse matrix
        The adjacency matrix of the network
    asmask : bool
        If set to ``True`` the result is returned as a boolean mask of the
        correct length with ``True`` values indicate the complements. The
        default is ``False`` which returns a list of indices instead.

    Returns
    -------
    An array containing indices of the nodes that are not part of the input
    list

    """
    inds = np.unique(inds)
    N = am.shape[0]
    mask = np.ones(shape=N, dtype=bool)
    mask[inds] = False
    if asmask:
        return mask
    else:
        return np.arange(N)[mask]


def find_connected_nodes(inds, am, flatten=True, logic='or'):
    r"""
    Finds which nodes are connected to the edges

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        nodes *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).
    flatten : bool (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input edge.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.
    logic : str
        Specifies logic to filter the resulting list.  Options are:

        ======= ===============================================================
        logic   Description
        ======= ===============================================================
        'or'    (default) All neighbors of the inputs.  This is also known as
                the 'union' in set theory or 'any' in boolean logic. Both
                keywords are accepted and treated as 'or'.
        'xor'   Only neighbors of one and only one inputs.  This is useful for
                finding neighbors that are not *shared* by any of the input
                nodes. 'exclusive_or' is also accepted.
        'xnor'  Neighbors that are shared by two or more inputs . This is
                equivalent to finding all neighbors with 'or', minus those
                found with 'xor', and is useful for finding neighbors that the
                inputs have in common.  'nxor' is also accepted.
        'and'   Only neighbors shared by all inputs. This is also known as
                'intersection' in set theory and (somtimes) as 'all' in
                boolean logic.  Both keywords are accepted and treated as
                'and'.
        ======= ===============================================================

    Returns
    -------
    An array containing the connected sites, filtered by the given logic.  If
    ``flatten`` is ``False`` then the result is a list of lists containing the
    neighbors of each given input edge.  In this latter case, nodes that
    have been removed by the given logic are indicated by ``nans``, thus the
    array is of type ``float`` and is not suitable for indexing.

    """
    if am.format != 'coo':
        raise Exception('Adjacency matrix must be in COO format')
    edges = np.array(inds, ndmin=1)
    if len(edges) == 0:
        return []
    # This function only uses the upper triangular portion, so make sure it
    # is sorted properly first
    if not istriu(am):
        am = sprs.triu(am, k=1)
    neighbors = np.hstack((am.row[edges], am.col[edges])).astype(np.int64)
    if neighbors.size:
        n_sites = np.amax(neighbors)
    if logic in ['or', 'union', 'any']:
        neighbors = np.unique(neighbors)
    elif logic in ['xor', 'exclusive_or']:
        neighbors = np.unique(np.where(np.bincount(neighbors) == 1)[0])
    elif logic in ['xnor']:
        neighbors = np.unique(np.where(np.bincount(neighbors) > 1)[0])
    elif logic in ['and', 'all', 'intersection']:
        temp = np.vstack((am.row[edges], am.col[edges])).T.tolist()
        temp = [set(pair) for pair in temp]
        neighbors = temp[0]
        [neighbors.intersection_update(pair) for pair in temp[1:]]
        neighbors = np.array(list(neighbors), dtype=np.int64, ndmin=1)
    else:
        raise Exception('Specified logic is not implemented')
    if flatten is False:
        if neighbors.size:
            mask = np.zeros(shape=n_sites + 1, dtype=bool)
            mask[neighbors] = True
            temp = np.hstack((am.row[edges], am.col[edges])).astype(np.int64)
            temp[~mask[temp]] = -1
            inds = np.where(temp == -1)[0]
            if len(inds):
                temp = temp.astype(float)
                temp[inds] = np.nan
            temp = np.reshape(a=temp, newshape=[len(edges), 2], order='F')
            neighbors = temp
        else:
            neighbors = [np.array([], dtype=np.int64) for i in range(len(edges))]
    return neighbors


def find_neighbor_edges(inds, im=None, am=None, flatten=True, logic='or'):
    r"""
    Finds all edges that are connected to the given input nodes

    Parameters
    ----------
    im : scipy.sparse matrix
        The incidence matrix of the network.  Must be shaped as (N-nodes,
        N-edges), with non-zeros indicating which nodes are connected. Either
        ``am`` or ``im`` must be given.  Passing in ``im`` is slower, but more
        powerful as it allows for an unflattened list of neighbors.
    am : scipy.sparse matrix (optional)
        The adjacency matrix of the network. Either ``am`` or ``im`` must be
        given.  Passing in ``am`` is faster, but does not allow for an
        unflattened list.
    flatten : bool (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input node.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.
    logic : str
        Specifies logic to filter the resulting list.  Options are:

        ======= ===============================================================
        logic   Description
        ======= ===============================================================
        'or'    (default) All neighbors of the inputs.  This is also known as
                the 'union' in set theory or 'any' in boolean logic. Both
                keywords are accepted and treated as 'or'.
        'xor'   Only neighbors of one and only one inputs.  This is useful for
                finding neighbors that are not *shared* by any of the input
                nodes. 'exclusive_or' is also accepted.
        'xnor'  Neighbors that are shared by two or more inputs . This is
                equivalent to finding all neighbors with 'or', minus those
                found with 'xor', and is useful for finding neighbors that the
                inputs have in common.  'nxor' is also accepted.
        'and'   Only neighbors shared by all inputs. This is also known as
                'intersection' in set theory and (somtimes) as 'all' in
                boolean logic.  Both keywords are accepted and treated as
                'and'.
        ======= ===============================================================

    Returns
    -------
    An array containing the neighboring edges filtered by the given logic. If
    ``flatten`` is ``False`` then the result is a list of lists containing the
    neighbors of each given input node.

    Notes
    -----
    The ``logic`` options are applied to neighboring edges only, thus it is not
    possible to obtain edges that are part of the global set but not neighbors.
    This is because (a) the list of global edges might be very large, and
    (b) it is not possible to return a list of neighbors for each input site
    if global sites are considered.

    """
    if im is not None:
        if im.format != 'lil':
            im = im.tolil(copy=False)
        rows = [im.rows[i] for i in np.array(inds, ndmin=1, dtype=np.int64)]
        if len(rows) == 0:
            return []
        neighbors = np.hstack(rows).astype(np.int64)
        n_bonds = int(im.nnz / 2)
        if logic in ['or', 'union', 'any']:
            neighbors = np.unique(neighbors)
        elif logic in ['xor', 'exclusive_or']:
            neighbors = np.unique(np.where(np.bincount(neighbors) == 1)[0])
        elif logic in ['xnor', 'shared']:
            neighbors = np.unique(np.where(np.bincount(neighbors) > 1)[0])
        elif logic in ['and', 'all', 'intersection']:
            neighbors = set(neighbors)
            [neighbors.intersection_update(i) for i in rows]
            neighbors = np.array(list(neighbors), dtype=int, ndmin=1)
        else:
            raise Exception('Specified logic is not implemented')
        if (flatten is False):
            if (neighbors.size > 0):
                mask = np.zeros(shape=n_bonds, dtype=bool)
                mask[neighbors] = True
                for i in range(len(rows)):
                    vals = np.array(rows[i], dtype=np.int64)
                    rows[i] = vals[mask[vals]]
                neighbors = rows
            else:
                neighbors = [np.array([], dtype=np.int64) for i in range(len(inds))]
        return neighbors
    elif am is not None:
        if am.format != 'coo':
            am = am.tocoo(copy=False)
        if not istriu(am):
            am = sprs.triu(am, k=1)
        if flatten is False:
            raise Exception('flatten cannot be used with an adjacency matrix')
        Ps = np.zeros(am.shape[0], dtype=bool)
        Ps[inds] = True
        conns = np.vstack((am.row, am.col)).T
        if logic in ['or', 'union', 'any']:
            neighbors = np.any(Ps[conns], axis=1)
        elif logic in ['xor', 'exclusive_or']:
            neighbors = np.sum(Ps[conns], axis=1) == 1
        elif logic in ['xnor', 'shared']:
            neighbors = np.all(Ps[conns], axis=1)
        elif logic in ['and', 'all', 'intersection']:
            raise Exception('Specified logic is not implemented')
        else:
            raise Exception('Specified logic is not implemented')
        neighbors = np.where(neighbors)[0]
        return neighbors
    else:
        raise Exception('Either the incidence or the adjacency matrix must be specified')


def find_neighbor_nodes(inds, am, flatten=True, include_input=False,
                        logic='or'):
    r"""
    Given a symmetric adjacency matrix, finds all nodes that are connected
    to the input nodes.

    Parameters
    ----------
    inds : array_like
        A list of node indices for which neighbors should be found
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
        'or'    (default) All neighbors of the inputs.  This is also known as
                the 'union' in set theory or 'any' in boolean logic. Both
                keywords are accepted and treated as 'or'.
        'xor'   Only neighbors of one and only one inputs.  This is useful for
                finding neighbors that are not *shared* by any of the input
                nodes. 'exclusive_or' is also accepted.
        'xnor'  Neighbors that are shared by two or more inputs . This is
                equivalent to finding all neighbors with 'or', minus those
                found with 'xor', and is useful for finding neighbors that the
                inputs have in common.  'nxor' is also accepted.
        'and'   Only neighbors shared by all inputs. This is also known as
                'intersection' in set theory and (somtimes) as 'all' in
                boolean logic.  Both keywords are accepted and treated as
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
    nodes = np.array(inds, ndmin=1)
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


def find_connecting_edges(inds, am):
    r"""
    Finds the edges that connects each pair of given nodes

    Parameters
    ----------
    inds : array_like
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
    nodes = np.array(inds, ndmin=2)
    if nodes.size == 0:
        return []
    z = tuple(zip(nodes[:, 0], nodes[:, 1]))
    neighbors = np.array([am.get(z[i], np.nan) for i in range(len(z))])
    return neighbors


def find_common_edges(conns, inds_1, inds_2):
    """
    Finds all bonds between two sets of nodes

    Parameters
    ----------
    conns : ndarray
        The connections of the sparse adjacency matrix in COO format
    inds2 : array_like
        A list of indices defining the first set of nodes
    P2 : array_like
        A list of indices defining the second set of nodes

    Returns
    -------
    ndarray
        List of interface edges between the two given sets of nodes

    """
    if np.intersect1d(inds_1, inds_2).size != 0:
        raise Exception("inds_1 and inds_2 must not share any pores.")
    edges_1 = find_neighbor_edges(inds_1, mode="xor")
    edges_2 = find_neighbor_edges(inds_2, mode="xor")
    return np.intersect1d(edges_1, edges_2)


def filter_by_z(conns, nodes, z=1):
    r"""
    Find nodes with a given number of neighbors

    Parameters
    ----------
    conns : ndarray
        The connections of the sparse adjacency matrix in COO format
    inds : array_like
        A list containing the indices of the nodes to be filtered
    z : int
        The coordination number by which to filter

    Returns
    -------
    pores : array_like
        A list of pores which satisfy the criteria

    """
    Nz = num_neighbors(inds=inds)
    orphans = np.where(Nz == z)[0]
    hits = nodes[orphans]
    return hits
