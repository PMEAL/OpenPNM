import numpy as np
import scipy.sparse as sprs
from openpnm._skgraph import settings


__all__ = [
    'trim_nodes',
    'trim_edges',
    'drop_nodes_from_am',
]


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
    for item in g.keys():
        if item.startswith(edge_prefix):
            g[item] = g[item][keep]
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
    if inds.dtype == bool:
        inds = np.where(inds)[0]
    keep = np.ones(N_sites, dtype=bool)
    keep[inds] = False
    for item in g.keys():
        if item.startswith(node_prefix):
            g[item] = g[item][keep]
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
