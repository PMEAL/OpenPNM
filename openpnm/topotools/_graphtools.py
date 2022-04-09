import logging
import warnings
import numpy as np
import scipy.sparse as sprs
from openpnm.utils import Workspace
from openpnm._skgraph import tools


logger = logging.getLogger(__name__)
ws = Workspace()
__all__ = [
    'find_neighbor_sites',
    'find_neighbor_bonds',
    'find_connected_sites',
    'find_connecting_bonds',
    'find_complement',
    'istriu',
    'istril',
    'istriangular',
    'issymmetric',
    'tri_to_am',
    'vor_to_am',
    'conns_to_am',
    'drop_sites',
    '_am_to_im',
    '_im_to_am',
]


def find_neighbor_sites(sites, am, flatten=True, include_input=False,
                        logic='or'):
    r"""
    Given a symmetric adjacency matrix, finds all sites that are connected
    to the input sites.

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    flatten : bool
        If ``True`` (default) the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    include_input : bool
        If ``False`` (default) the input sites will be removed from the result.

    logic : str
        Specifies logic to filter the resulting list.  Options are:

        **'or'** : (default) All neighbors of the input sites.  This is also
        known as the 'union' in set theory or 'any' in boolean logic.  Both
        keywords are accepted and treated as 'or'.

        **'xor'** : Only neighbors of one and only one input site.  This is
        useful for finding the sites that are not shared by any of the input
        sites.  'exclusive_or' is also accepted.

        **'xnor'** : Neighbors that are shared by two or more input sites. This
        is equivalent to finding all neighbors with 'or', minus those found
        with 'xor', and is useful for finding neighbors that the inputs have
        in common.  'nxor' is also accepted.

        **'and'** : Only neighbors shared by all input sites.  This is also
        known as 'intersection' in set theory and (somtimes) as 'all' in
        boolean logic.  Both keywords are accepted and treated as 'and'.

    Returns
    -------
    An array containing the neighboring sites filtered by the given logic.  If
    ``flatten`` is ``False`` then the result is a list of lists containing the
    neighbors of each input site.

    See Also
    --------
    find_complement

    Notes
    -----
    The ``logic`` options are applied to neighboring sites only, thus it is not
    possible to obtain sites that are part of the global set but not neighbors.
    This is because (a) the list global sites might be very large, and (b) it
    is not possible to return a list of neighbors for each input site if global
    sites are considered.

    """
    if am.format != 'lil':
        am = am.tolil(copy=False)
    sites = np.array(sites, ndmin=1)
    if len(sites) == 0:
        return []
    am_coo = am.tocoo()
    n_sites = am.shape[0]
    rows = am.rows[sites].tolist()
    if len(rows) == 0:
        return []
    neighbors = am_coo.col[np.in1d(am_coo.row, sites)]
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
    mask = np.zeros(shape=n_sites, dtype=bool)
    mask[neighbors] = True
    if not include_input:
        mask[sites] = False
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
            neighbors = [np.array([], dtype=int) for i in range(len(sites))]
    return neighbors


def find_neighbor_bonds(sites, im=None, am=None, flatten=True, logic='or'):
    r"""
    Given an incidence matrix, finds all sites that are connected to the
    input sites.

    This function accepts either an incidence or adjacency matrix.

    Parameters
    ----------
    im : scipy.sparse matrix
        The incidence matrix of the network.  Must be shaped as (N-sites,
        N-bonds), with non-zeros indicating which sites are connected. Either
        ``am`` or ``im`` must be given.  Passing in ``im`` is slower, but more
        powerful as it allows for an unflattened list of neighbors.

    am : scipy.sparse matrix (optional)
        The adjacency matrix of the network. Either ``am`` or ``im`` must be
        given.  Passing in ``am`` is faster, but does not allow for an
        unflattened list.

    flatten : bool (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    logic : str
        Specifies logic to filter the resulting list.  Options are:

        **'or'** : (default) All neighbors of the input sites.  This is also
        known as the 'union' in set theory or 'any' in boolean logic.  Both
        keywords are accepted and treated as 'or'.

        **'xor'** : Only neighbors of one and only one input site.  This is
        useful for finding the sites that are not shared by any of the input
        sites.  'exclusive_or' is also accepted'.

        **'xnor'** : Neighbors that are shared by two or more input sites. This
        is equivalent to finding all neighbors with 'or', minus those found
        with 'xor', and is useful for finding neighbors that the inputs have
        in common.  'nxor' is also accepted.

        **'and'** : Only neighbors shared by all input sites.  This is also
        known as 'intersection' in set theory and (somtimes) as 'all' in
        boolean logic.  Both keywords are accepted and treated as 'and'.

    Returns
    -------
    An array containing the neighboring bonds filtered by the given logic.  If
    ``flatten`` is ``False`` then the result is a list of lists containing the
    neighbors of each given input site.

    See Also
    --------
    find_complement

    Notes
    -----
    The ``logic`` options are applied to neighboring bonds only, thus it is not
    possible to obtain bonds that are part of the global set but not neighbors.
    This is because (a) the list of global bonds might be very large, and
    (b) it is not possible to return a list of neighbors for each input site
    if global sites are considered.

    """
    if im is not None:
        if im.format != 'lil':
            im = im.tolil(copy=False)
        rows = [im.rows[i] for i in np.array(sites, ndmin=1, dtype=np.int64)]
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
                neighbors = [np.array([], dtype=np.int64) for i in range(len(sites))]
        return neighbors
    elif am is not None:
        if am.format != 'coo':
            am = am.tocoo(copy=False)
        if not istriu(am):
            am = sprs.triu(am, k=1)
        if flatten is False:
            raise Exception('flatten cannot be used with an adjacency matrix')
        Ps = np.zeros(am.shape[0], dtype=bool)
        Ps[sites] = True
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


def find_connected_sites(bonds, am, flatten=True, logic='or'):
    r"""
    Given an adjacency matrix, finds which sites are connected to the input
    bonds.

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    flatten : bool (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    logic : str
        Specifies logic to filter the resulting list.  Options are:

        **'or'** : (default) All neighbors of the input bonds.  This is also
        known as the 'union' in set theory or (sometimes) 'any' in boolean
        logic.  Both keywords are accepted and treated as 'or'.

        **'xor'** : Only neighbors of one and only one input bond.  This is
        useful for finding the sites that are not shared by any of the input
        bonds.  'exclusive_or' is also accepted.

        **'xnor'** : Neighbors that are shared by two or more input bonds. This
        is equivalent to finding all neighbors with 'or', minus those found
        with 'xor', and is useful for finding neighbors that the inputs have
        in common.  'nxor' is also accepted.

        **'and'** : Only neighbors shared by all input bonds.  This is also
        known as 'intersection' in set theory and (somtimes) as 'all' in
        boolean logic.  Both keywords are accepted and treated as 'and'.

    Returns
    -------
    An array containing the connected sites, filtered by the given logic.  If
    ``flatten`` is ``False`` then the result is a list of lists containing the
    neighbors of each given input bond.  In this latter case, sites that
    have been removed by the given logic are indicated by ``nans``, thus the
    array is of type ``float`` and is not suitable for indexing.

    See Also
    --------
    find_complement

    """
    if am.format != 'coo':
        raise Exception('Adjacency matrix must be in COO format')
    bonds = np.array(bonds, ndmin=1)
    if len(bonds) == 0:
        return []
    # This function only uses the upper triangular portion, so make sure it
    # is sorted properly first
    if not istriu(am):
        am = sprs.triu(am, k=1)
    neighbors = np.hstack((am.row[bonds], am.col[bonds])).astype(np.int64)
    if neighbors.size:
        n_sites = np.amax(neighbors)
    if logic in ['or', 'union', 'any']:
        neighbors = np.unique(neighbors)
    elif logic in ['xor', 'exclusive_or']:
        neighbors = np.unique(np.where(np.bincount(neighbors) == 1)[0])
    elif logic in ['xnor']:
        neighbors = np.unique(np.where(np.bincount(neighbors) > 1)[0])
    elif logic in ['and', 'all', 'intersection']:
        temp = np.vstack((am.row[bonds], am.col[bonds])).T.tolist()
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
            temp = np.hstack((am.row[bonds], am.col[bonds])).astype(np.int64)
            temp[~mask[temp]] = -1
            inds = np.where(temp == -1)[0]
            if len(inds):
                temp = temp.astype(float)
                temp[inds] = np.nan
            temp = np.reshape(a=temp, newshape=[len(bonds), 2], order='F')
            neighbors = temp
        else:
            neighbors = [np.array([], dtype=np.int64) for i in range(len(bonds))]
    return neighbors


def find_connecting_bonds(sites, am):
    r"""
    Given pairs of sites, finds the bonds which connects each pair.

    Parameters
    ----------
    sites : array_like
        A 2-column vector containing pairs of site indices on each row.

    am : scipy.sparse matrix
        The adjacency matrix of the network.  Must be symmetrical such that if
        sites *i* and *j* are connected, the matrix contains non-zero values
        at locations (i, j) and (j, i).

    Returns
    -------
    Returns a list the same length as P1 (and P2) with each element
    containing the throat number that connects the corresponding pores,
    or `None`` if pores are not connected.

    Notes
    -----
    The returned list can be converted to an ndarray, which will convert
    the ``None`` values to ``nan``.  These can then be found using
    ``numpy.isnan``.

    """
    if am.format != 'dok':
        am = am.todok(copy=False)
    sites = np.array(sites, ndmin=2)
    if sites.size == 0:
        return []
    z = tuple(zip(sites[:, 0], sites[:, 1]))
    neighbors = [am.get(z[i], None) for i in range(len(z))]
    return neighbors


def find_complement(am, sites=None, bonds=None, asmask=False):
    r"""
    Finds the complementary sites (or bonds) to a given set of inputs

    Parameters
    ----------
    am : scipy.sparse matrix
        The adjacency matrix of the network.
    sites : array_like (optional)
        The set of sites for which the complement is sought
    bonds : array_like (optional)
        The set of bonds for which the complement is sought
    asmask : bool
        If set to ``True`` the result is returned as a boolean mask of the
        correct length with ``True`` values indicate the complements.  The
        default is ``False`` which returns a list of indices instead.

    Returns
    -------
    An array containing indices of the sites (or bonds) that are not part of
    the input list.

    Notes
    -----
    Either ``sites`` or ``bonds`` must be specified

    """
    if (sites is not None) and (bonds is None):
        inds = np.unique(sites)
        N = am.shape[0]
    elif (bonds is not None) and (sites is None):
        inds = np.unique(bonds)
        N = int(am.nnz / 2)
    elif (bonds is not None) and (sites is not None):
        raise Exception('Only one of sites or bonds can be specified')
    else:
        raise Exception('Either sites or bonds must be specified')
    mask = np.ones(shape=N, dtype=bool)
    mask[inds] = False
    if asmask:
        return mask
    else:
        return np.arange(N)[mask]


def istriu(am):
    return tools.istriu(am)


istriu.__doc__ = tools.istriu.__doc__


def istril(am):
    return tools.istril(am)


istril.__doc__ = tools.istril.__doc__


def istriangular(am):
    return tools.istriangular(am)


istriangular.__doc__ = tools.istriangular.__doc__


def issymmetric(am):
    return tools.issymmetric(am)


issymmetric.__doc__ = tools.issymmetric.__doc__


def _am_to_im(am):
    return tools.am_to_im(am)


_am_to_im.__doc__ = tools.am_to_im.__doc__


def _im_to_am(im):
    return tools.im_to_am(im)


_im_to_am.__doc__ = tools.im_to_am.__doc__


def tri_to_am(tri):
    return tools.tri_to_am(tri=tri)


tri_to_am.__doc__ = tools.tri_to_am.__doc__


def vor_to_am(vor):
    return tools.vor_to_am(vor=vor)


vor_to_am.__doc__ = tools.vor_to_am.__doc__


def conns_to_am(*args, **kwargs):
    return tools.conns_to_am(*args, **kwargs)


conns_to_am.__doc__ = tools.conns_to_am.__doc__


def drop_sites(am, sites):
    return tools.drop_nodes_from_am(am=am, sites=sites)


drop_sites.__doc__ = tools.drop_nodes_from_am.__doc__
