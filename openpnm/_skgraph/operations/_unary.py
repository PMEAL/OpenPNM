import numpy as np
import scipy.sparse as sprs
from openpnm._skgraph import settings, tools


__all__ = [
    'trim_nodes',
    'trim_edges',
    'drop_nodes_from_am',
    'add_nodes',
    'add_edges',
    'split_edges',
]


def add_nodes(network, new_coords):
    r"""
    Given a list of node coordinates, add them to the network

    Parameters
    ----------
    network : dict
        A dictionary containing the node and edge attributes as ndarrays
    new_coords : ndarray
        The N-by-3 array of coordinates of the new nodes

    Returns
    -------
    network : dict
        The network dictionary with the new nodes added to the end. Note that
        any existing node attributes are also extended and filled with default
        values specified in ``settings.default_values``

    """

    # TODO: At present this does not work with empty networks.  It's a bit
    # challenging to determine what size each 'empty' array should become.
    # It's fine for 1D arrays like ``np.array([])``, but for N-dimensional
    # stuff it's trickier. For instance ``np.array([[], []])`` has a shape
    # of (2, 0) so the empty dimension is the wrong one since all the array
    # extending in this function occurs on the 1st axis.
    g = network
    node_prefix = tools.get_node_prefix(g)
    coords = np.atleast_2d(new_coords)
    Nnew = coords.shape[0]
    for k, v in g.items():
        if k.startswith(node_prefix):
            dval = None
            for t in settings.missing_values.keys():
                if v.dtype == t:
                    dval = settings.missing_values[t]
            blank = np.repeat(v[:1, ...], Nnew, axis=0)*dval
            blank.fill(dval)
            g[k] = np.concatenate((v, blank), axis=0).astype(v.dtype)
    # Lastly, overwrite the -Nnew elements of coords with the given values
    g[node_prefix+'.coords'][-Nnew:] = np.array(coords)
    return g


def add_edges(network, new_conns):
    r"""
    Given a list of edge connections, add them to the network

    Parameters
    ----------
    network : dict
        A dictionary containing the node and edge attributes as ndarrays
    new_conns : ndarray
        The N-by-2 array of connections betwween existing nodes

    Returns
    -------
    network : dict
        The network dictionary with the new edges added to the end. Note that
        any existing edge attributes are also extended and filled with default
        values specified in ``settings.default_values``.

    """

    # TODO: At present this does not work with empty networks.  It's a bit
    # challenging to determine what size each 'empty' array should become.
    # It's fine for 1D arrays like ``np.array([])``, but for N-dimensional
    # stuff it's trickier. For instance ``np.array([[], []])`` has a shape
    # of (2, 0) so the empty dimension is the wrong one since all the array
    # extending in this function occurs on the 1st axis.
    g = network
    edge_prefix = tools.get_edge_prefix(g)
    conns = np.atleast_2d(new_conns)
    Nnew = conns.shape[0]
    for k, v in g.items():
        if k.startswith(edge_prefix):
            dval = None
            for t in settings.missing_values.keys():
                if v.dtype == t:
                    dval = settings.missing_values[t]
            blank = np.repeat(v[:1, ...], Nnew, axis=0)*dval
            blank.fill(dval)
            g[k] = np.concatenate((v, blank), axis=0).astype(v.dtype)
    # Lastly, overwrite the -Nnew elements of coords with the given values
    g[edge_prefix+'.conns'][-Nnew:] = np.array(conns)
    return g


def trim_edges(network, inds):
    r"""
    Removes given edges from a graph or network

    Parameters
    ----------
    network : dictionary
        A dictionary containing coords, conns and other attributes
    inds : array_like
        The edge indices to be trimmed in the form of a 1D list or boolean
        mask with ``True`` values indicating indices to trim.

    Returns
    -------
    network : dict
        The dictionary with all edge arrays trimmed accordingly

    """
    g = network
    edge_prefix = tools.get_edge_prefix(g)
    N_bonds = g[edge_prefix+'.conns'].shape[0]
    inds = np.atleast_1d(inds)
    keep = np.ones(N_bonds, dtype=bool)
    keep[inds] = False
    for k, v in g.items():
        if k.startswith(edge_prefix):
            g[k] = v[keep]
    return g


def trim_nodes(network, inds):
    r"""
    Removes given nodes and any connected edges from a graph or network

    Parameters
    ----------
    network : dict
        A dictionary containing coords, conns and other attributes
    inds : array_like
        The node indices to be trimmed in the form of a 1D list or boolean
        mask with ``True`` values indicating indices to trim.

    Returns
    -------
    network : dict
        The dictionary with all nodes arrays trimmed accordingly, all edges
        trimmed that were connected to the trimmed nodes, and the 'edge.conns'
        array renumbered so edges point to the updated node indices.

    """
    node_prefix = tools.get_node_prefix(network)
    edge_prefix = tools.get_edge_prefix(network)
    N_sites = network[node_prefix+'.coords'].shape[0]
    inds = np.atleast_1d(inds)
    if inds.dtype == bool:
        inds = np.where(inds)[0]
    keep = np.ones(N_sites, dtype=bool)
    keep[inds] = False
    for k, v in network.items():
        if k.startswith(node_prefix):
            network[k] = v[keep]
    # Remove edges
    edges = np.any(np.isin(network[edge_prefix+'.conns'], inds), axis=1)
    network = trim_edges(network, inds=edges)
    # Renumber conns
    remapping = np.cumsum(keep) - 1
    network[edge_prefix+'.conns'] = remapping[network[edge_prefix+'.conns']]
    return network


def drop_nodes_from_am(am, inds):
    r"""
    Update adjacency matrix after dropping nodes

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network in COO format.
    inds : array_like
        A list of which nodes indices to drop.  Can either be integer indices
        or a boolean mask with ``True`` indicating which locations to drop

    Returns
    -------
    am : ndarray
        An updated adjacency matrix with nodes and headless edges removed,
        and node indices updated accordingly
    dropped_edges : ndarray
        A boolean array with ``True`` values indicating which edges
        were rendered headless. This can be used to drop invalid edges
        from other arrays (i.e. array = array[~dropped_edges]).

    """
    nodes = np.array(inds)
    if nodes.dtype != bool:
        inds = np.copy(nodes)
        nodes = np.zeros(am.shape[0], dtype=bool)
        nodes[inds] = True
    node_mask = ~nodes
    conns = np.vstack((am.row, am.col)).T
    node_id = np.cumsum(node_mask) - 1
    edge_mask = ~np.all(node_mask[conns], axis=1)
    conns = node_id[conns[~edge_mask]]
    am = sprs.coo_matrix((am.data[~edge_mask], (conns[:, 0], conns[:, 1])))
    return am, edge_mask


def split_edges(network):
    r"""
    Inserts an new node between each existing node and joins with new edges

    Parameters
    ----------
    network : dict
        The dictionary containing the network connections and coordinates

    Returns
    -------
    result : tuple
        A tuple containing ``new_conns`` and optionally ``new_coords`` if
        ``coords`` was provided.

        ============== ========================================================
        Value          Description
        ============== ========================================================
        ``new_conns``  A new adjacency matrix in COO format with new nodes
                       added between each original node. If edge 1 connected
                       nodes 1 and 2, then row 1 of the new sparse adjacency
                       matrix will be [1, Nt + 1], and row Nt + 1 will be
                       [Nt + 1, 2].
        ``new_coords`` A and updated list of node coordinates with the new
                       nodes appended to the end.  The coordinates of the new
                       nodes are taken as the average of the two nodes between
                       which they were inserted.
        ============== ========================================================

    """
    g = network
    node_prefix = tools.get_node_prefix(g)
    edge_prefix = tools.get_edge_prefix(g)
    conns = g[edge_prefix + '.conns']
    coords = g[node_prefix + '.coords']
    Nt = conns.shape[0]
    Np = conns.max() + 1
    new_conns = np.zeros([2*Nt, 2], dtype=int)
    new_conns[:Nt, :] = np.vstack((conns[:, 0], np.arange(Np, Np+Nt))).T
    new_conns[Nt:, :] = np.vstack((np.arange(Np, Np+Nt), conns[:, 1])).T
    result = (new_conns, )
    if coords is not None:
        Np = coords.shape[0]
        new_coords = np.zeros([Np + Nt, 3], dtype=float)
        new_coords[:Np, :] = coords
        new_coords[Np:, :] = np.mean(coords[conns], axis=1)
        result = (new_conns, new_coords)
    return result
