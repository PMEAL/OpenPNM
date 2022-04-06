import numpy as np


__all__ = [
    'trim_nodes',
    'trim_edges',
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
    N_bonds = g['edge.conns'].shape[0]
    inds = np.atleast_1d(inds)
    keep = np.ones(N_bonds, dtype=bool)
    keep[inds] = False
    for item in g.keys():
        if item.startswith('edge'):
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
    N_sites = g['node.coords'].shape[0]
    inds = np.atleast_1d(inds)
    if inds.dtype == bool:
        inds = np.where(inds)[0]
    keep = np.ones(N_sites, dtype=bool)
    keep[inds] = False
    for item in g.keys():
        if item.startswith('node'):
            g[item] = g[item][keep]
    # Remove edges
    edges = np.any(np.isin(g['edge.conns'], inds), axis=1)
    g = trim_edges(g, inds=edges)
    # Renumber throat conns
    remapping = np.cumsum(keep) - 1
    g['edge.conns'] = remapping[g['edge.conns']]
    return g
