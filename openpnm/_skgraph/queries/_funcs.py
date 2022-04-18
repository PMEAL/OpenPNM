import numpy as np
import scipy.sparse as sprs
from scipy.sparse import csgraph
from openpnm._skgraph.tools import istriu, conns_to_am, dict_to_am, dict_to_im
from openpnm._skgraph.tools import get_node_prefix, get_edge_prefix


__all__ = [
    'find_common_edges',
    'filter_by_z',
    'find_connecting_edges',
    'find_neighbor_nodes',
    'find_neighbor_edges',
    'find_connected_nodes',
    'find_complementary_nodes',
    'find_complementary_edges',
    'find_path',
    'find_coordination',
]


def find_complementary_edges(inds, g, asmask=False):
    r"""
    Finds the complementary edges to a given set of inputs

    Parameters
    ----------
    inds : array_like
        A list of edge indices for which the complement is sought
    g : dict
        Graph dictionary
    asmask : bool
        If set to ``True`` the result is returned as a boolean mask of the
        correct length with ``True`` values indicate the complements.  The
        default is ``False`` which returns a list of indices instead.

    Returns
    -------
    An array containing indices of the edges that are not part of the input
    list

    """
    edge_prefix = get_edge_prefix(g)
    inds = np.unique(inds)
    N = g[edge_prefix+'.conns'].shape[0]
    mask = np.ones(shape=N, dtype=bool)
    mask[inds] = False
    if asmask:
        return mask
    else:
        return np.arange(N)[mask]


def find_complementary_nodes(inds, g, asmask=False):
    r"""
    Finds the complementary nodes to a given set of inputs

    Parameters
    ----------
    inds : array_like (optional)
        A list of indices for which the complement is sought
    g : dict
        Graph dictionary
    asmask : bool
        If set to ``True`` the result is returned as a boolean mask of the
        correct length with ``True`` values indicate the complements. The
        default is ``False`` which returns a list of indices instead.

    Returns
    -------
    An array containing indices of the nodes that are not part of the input
    list

    """
    node_prefix = get_node_prefix(g)
    inds = np.unique(inds)
    N = g[node_prefix+'.coords'].shape[0]
    mask = np.ones(shape=N, dtype=bool)
    mask[inds] = False
    if asmask:
        return mask
    else:
        return np.arange(N)[mask]


def find_connected_nodes(inds, g=None, am=None, flatten=True, logic='or'):
    r"""
    Finds which nodes are connected to a given set of edges

    Parameters
    ----------
    inds : array_like
        A list of edges indices whose connected nodes are sought
    g : dict, optional
        The graph dictionary. Either this or ``am`` must be given.
    am : scipy.sparse matrix, optional.
        The adjacency matrix of the network. Either this or ``g`` must be given.
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
                'intersection' in set theory and (sometimes) as 'all' in
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
    edges = np.array(inds, ndmin=1)
    if len(edges) == 0:  # Short-circuit this fuction if edges is empty
        return []
    if g is not None:
        am = dict_to_am(g)
    elif am is not None:
        if am.format != 'coo':
            am = am.tocoo()
    else:
        raise Exception('either g or am must be given')
    # This function only uses the upper triangular portion, so make sure it
    # is sorted properly first
    if not istriu(am):
        am = sprs.triu(am, k=1)
    neighbors = np.hstack((am.row[edges], am.col[edges]))
    if neighbors.size:
        n_sites = np.amax(neighbors)
    if logic in ['or', 'union', 'any']:
        neighbors = np.unique(neighbors)
    elif logic in ['xor', 'exclusive_or']:
        neighbors = np.unique(np.where(np.bincount(neighbors) == 1)[0])
    elif logic in ['xnor', 'nxor']:
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


def find_neighbor_edges(inds, g, flatten=True, logic='or'):
    r"""
    Finds all edges that are connected to the given input nodes

    Parameters
    ----------
    inds : array_like (optional)
        A list of node indices whose neighbor edges are sought
    g : dict
        The graph dictionary
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
    if flatten == False:
        im = dict_to_im(g)
        am = None
    else:
        am = dict_to_am(g)
        im = None
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


def find_neighbor_nodes(inds, g=None, am=None, flatten=True, include_input=False,
                        logic='or'):
    r"""
    Ginds all nodes that are directly connected to the input nodes

    Parameters
    ----------
    inds : array_like
        A list of node indices whose neighbors are sought
    g : dict, optional
        The graph dictionary. Either this or ``am`` must be provided.
    am : scipy.sparse matrix, optional
        The adjacency matrix of the network. Must be symmetrical such that if
        nodes *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i). Either this or ``g`` must be provided.
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
    nodes = np.array(inds, ndmin=1)
    # Short-circuit the function if the input list is already empty
    if len(nodes) == 0:
        return []
    if g is not None:
        am = dict_to_am(g)
    am_coo = am.tocoo()
    if am.format != 'lil':
        am = am.tolil(copy=False)
    rows = am.rows[nodes].tolist()
    if len(rows) == 0:
        return []
    n_nodes = am.shape[0]
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


def find_connecting_edges(inds, g=None, am=None):
    r"""
    Finds the edge that connects each pair of given nodes

    Parameters
    ----------
    inds : array_like
        A 2-column vector containing pairs of node indices
    g : dict, optional
        The graph dictionary.  Either this or ``am`` must be provided
    am : scipy.sparse matrix, optional
        The adjacency matrix of the network. Must be symmetrical such that if
        nodes *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i). Either this or ``g`` must be provided.

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
    nodes = np.array(inds, ndmin=2)
    # Short-circuit function if nodes is an empty list
    if nodes.size == 0:
        return []
    if g is not None:
        am = dict_to_am(g)
    elif am is not None:
        pass
    else:
        raise Exception('Either g or am must be provided')
    if am.format != 'dok':
        am = am.todok(copy=True)
    z = tuple(zip(nodes[:, 0], nodes[:, 1]))
    neighbors = np.array([am.get(z[i], np.nan) for i in range(len(z))])
    return neighbors


def find_common_edges(g, inds_1, inds_2):
    """
    Finds edges shared between two sets of nodes

    Parameters
    ----------
    g : dict
        The graph dictionary
    inds_1 : array_like
        A list of indices defining the first set of nodes
    inds_2 : array_like
        A list of indices defining the second set of nodes

    Returns
    -------
    edges : ndarray
        List of edge indices connecting the two given sets of nodes

    """
    if np.intersect1d(inds_1, inds_2).size != 0:
        raise Exception("inds_1 and inds_2 must not share any nodes")
    edges_1 = find_neighbor_edges(inds=inds_1, g=g, logic="xor")
    edges_2 = find_neighbor_edges(inds=inds_2, g=g, logic="xor")
    return np.intersect1d(edges_1, edges_2)


def filter_by_z(g, inds, z=1):
    r"""
    Filters a list of nodes to those with a given number of neighbors

    Parameters
    ----------
    g : dict
        The graph dictionary
    inds : array_like
        A list containing the indices of the nodes to be filtered
    z : int
        The coordination number by which to filter

    Returns
    -------
    inds : array_like
        A list of node indices which satisfy the criteria

    """
    inds = np.array(inds)
    coordination = find_coordination(g)
    hits = coordination == z
    inds = inds[hits[inds]]
    return inds


def find_coordination(g, nodes=None):
    r"""
    Find the coordination number of nodes

    Parameters
    ----------
    g : dict
        The graph dictionary
    nodes : array_like, optional
        The nodes for which coordination is sought. If not provided then
        coordination for *all* nodes is returned

    Returns
    -------
    z : ndarray
        An array containing the number of neighbors for each given node

    Notes
    -----
    This gives the same result for both directed and undirected graphs since
    'coordination' is direction independent. It also gives the same result
    for multigraphs since it counts the number of neighboring nodes, not
    the number of incoming/outgoing edges.
    """
    am = dict_to_am(g, directed=False)
    z = am.getnnz(axis=1)
    if nodes is None:
        return z
    else:
        return z[np.array(nodes)]


def find_path(g, pairs, weights=None):
    r"""
    Find the shortest path between pairs of nodes

    Parameters
    ----------
    g : dict
        The graph dictionary
    pairs : array_like
        An N x 2 array containing N pairs of nodes between which the shortest
        path is sought
    weights : ndarray, optional
        The edge weights to use when traversing the path. If not provided
        then 1's will be used.

    Returns
    -------
    A dictionary containing both the nodes and edges that define the
    shortest path connecting each pair of input nodes.

    Notes
    -----
    The shortest path is found using Dijkstra's algorithm included in the
    ``scipy.sparse.csgraph`` module

    TODO: The returned throat path contains the correct values, but not
    necessarily in the true order

    """
    am = dict_to_am(g)
    if weights is not None:
        am.data = np.ones_like(am.row, dtype=int)
    pairs = np.array(pairs, ndmin=2)
    paths = csgraph.dijkstra(csgraph=am, indices=pairs[:, 0],
                             return_predecessors=True, min_only=False)[1]
    am.data = np.hstack(2*[np.arange(am.data.size/2)]).astype(int)
    dok = am.todok()
    nodes = []
    edges = []
    for row in range(0, np.shape(pairs)[0]):
        j = pairs[row][1]
        ans = []
        while paths[row][j] > -9999:
            ans.append(j)
            j = paths[row][j]
        ans.append(pairs[row][0])
        ans.reverse()
        nodes.append(np.array(ans, dtype=int))

        edges.append(np.array(Ts, dtype=int))
    return {'nodes': nodes, 'edges': edges}
