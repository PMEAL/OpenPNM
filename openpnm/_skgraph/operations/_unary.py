import numpy as np
import scipy.sparse as sprs
from openpnm._skgraph import settings


__all__ = [
    'trim_nodes',
    'trim_edges',
    'drop_nodes_from_am',
    'add_nodes',
    'add_edges',
]


def add_nodes(g, coords):
    r"""
    Given a list of node coordinates, add them to the network

    Parameters
    ----------
    g : dict
        A dictionary containing the node and edge attributes as ndarrays
    coords : ndarray
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
    # extending in this functionoccurs on the 1st axis.

    node_prefix = settings.node_prefix
    Nnew = coords.shape[0]
    for k, v in g.items():
        if k.startswith(node_prefix):
            blank = np.repeat(v[:1, ...], Nnew, axis=0)
            dval = None
            for t in settings.missing_values.keys():
                if v.dtype == t:
                    dval = settings.missing_values[t]
            blank.fill(dval)
            g[k] = np.concatenate((v, blank), axis=0)
    # Lastly, overwrite the -Nnew elements of coords with the given values
    g[node_prefix+'.coords'][-Nnew:] = np.array(coords)
    return g


def add_edges(g, conns):
    r"""
    Given a list of edge connections, add them to the network

    Parameters
    ----------
    network : dict
        A dictionary containing the node and edge attributes as ndarrays
    conns : ndarray
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
    # extending in this functionoccurs on the 1st axis.

    edge_prefix = settings.edge_prefix
    Nnew = np.array(conns).shape[0]
    for k, v in g.items():
        if k.startswith(edge_prefix):
            blank = np.repeat(v[:1, ...], Nnew, axis=0)
            dval = None
            for t in settings.missing_values.keys():
                if v.dtype == t:
                    dval = settings.missing_values[t]
            blank.fill(dval)
            g[k] = np.concatenate((v, blank), axis=0)
    # Lastly, overwrite the -Nnew elements of coords with the given values
    g[edge_prefix+'.conns'][-Nnew:] = np.array(conns)
    return g


def trim_edges(g, inds):
    r"""
    Removes given edge from a graph or network

    Parameters
    ----------
    g : dictionary
        A dictionary containing 'edge.conns' and other edge attributes in the
        form of 'edge.<attribute>'.
    inds : array_like
        The edge indices to be trimmed in the form of a 1D list or boolean
        mask with ``True`` values indicating indices to trim.

    Returns
    -------
    g : dict
        The dictionary with all edge arrays trimmed accordingly

    """
    edge_prefix = settings.edge_prefix
    N_bonds = g[edge_prefix+'.conns'].shape[0]
    inds = np.atleast_1d(inds)
    keep = np.ones(N_bonds, dtype=bool)
    keep[inds] = False
    for k, v in g.items():
        if k.startswith(edge_prefix):
            g[k] = v[keep]
    return g


def trim_nodes(g, inds):
    r"""
    Removes given nodes and any connected edges from a graph or network

    Parameters
    ----------
    g : dictionary
        A dictionary containing 'node.coods' and other node attributes in the
        form of 'node.<attribute>'.
    inds: array_like
        The node indices to be trimmed in the form of a 1D list or boolean
        mask with ``True`` values indicating indices to trim.

    Returns
    -------
    g : dict
        The dictionary with all nodes arrays trimmed accordingly, all edges
        trimmed that were connected to the trimmed nodes, and the 'edge.conns'
        array renumbered so edges point to the updated node indices.

    """
    node_prefix = settings.node_prefix
    edge_prefix = settings.edge_prefix
    N_sites = g[node_prefix+'.coords'].shape[0]
    inds = np.atleast_1d(inds)
    keep = np.ones(N_sites, dtype=bool)
    keep[inds] = False
    for k, v in g.items():
        if k.startswith(node_prefix):
            g[k] = v[keep]
    # Remove edges
    edges = np.any(np.isin(g[edge_prefix+'.conns'], inds), axis=1)
    g = trim_edges(g, inds=edges)
    # Renumber throat conns
    remapping = np.cumsum(keep) - 1
    g[edge_prefix+'.conns'] = remapping[g[edge_prefix+'.conns']]
    return g


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
