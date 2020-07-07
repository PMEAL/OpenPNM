import warnings
import numpy as np
import scipy as sp
import scipy.sparse as sprs
import scipy.ndimage as spim
from scipy.sparse import csgraph
from scipy.spatial import ConvexHull
from openpnm.utils import PrintableDict, logging, Workspace
logger = logging.getLogger(__name__)
ws = Workspace()


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

    flatten : boolean
        If ``True`` (default) the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    include_input : boolean
        If ``False`` (default) the input sites will be removed from the result.

    logic : string
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
        neighbors = np.array(list(neighbors), dtype=sp.int64, ndmin=1)
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
                vals = np.array(rows[i], dtype=sp.int64)
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

    flatten : boolean (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    logic : string
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
        rows = [im.rows[i] for i in np.array(sites, ndmin=1, dtype=sp.int64)]
        if len(rows) == 0:
            return []
        neighbors = np.hstack(rows).astype(sp.int64)
        n_bonds = int(im.nnz/2)
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
                    vals = np.array(rows[i], dtype=sp.int64)
                    rows[i] = vals[mask[vals]]
                neighbors = rows
            else:
                neighbors = [np.array([], dtype=sp.int64) for i in range(len(sites))]
        return neighbors
    elif am is not None:
        if am.format != 'coo':
            am = am.tocoo(copy=False)
        if not istriu(am):
            am = sp.sparse.triu(am, k=1)
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
        raise Exception('Either the incidence or the adjacency matrix must '
                        + 'must be specified')


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

    flatten : boolean (default is ``True``)
        Indicates whether the returned result is a compressed array of all
        neighbors, or a list of lists with each sub-list containing the
        neighbors for each input site.  Note that an *unflattened* list might
        be slow to generate since it is a Python ``list`` rather than a Numpy
        array.

    logic : string
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
        am = sp.sparse.triu(am, k=1)
    neighbors = np.hstack((am.row[bonds], am.col[bonds])).astype(sp.int64)
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
        neighbors = np.array(list(neighbors), dtype=sp.int64, ndmin=1)
    else:
        raise Exception('Specified logic is not implemented')
    if flatten is False:
        if neighbors.size:
            mask = np.zeros(shape=n_sites+1, dtype=bool)
            mask[neighbors] = True
            temp = np.hstack((am.row[bonds], am.col[bonds])).astype(sp.int64)
            temp[~mask[temp]] = -1
            inds = np.where(temp == -1)[0]
            if len(inds):
                temp = temp.astype(float)
                temp[inds] = sp.nan
            temp = np.reshape(a=temp, newshape=[len(bonds), 2], order='F')
            neighbors = temp
        else:
            neighbors = [np.array([], dtype=sp.int64) for i in range(len(bonds))]
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
    The returned list can be converted to an ND-array, which will convert
    the ``None`` values to ``nan``.  These can then be found using
    ``scipy.isnan``.

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

    asmask : boolean
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
        N = int(am.nnz/2)
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
    r"""
    Returns ``True`` is the sparse adjacency matrix is upper triangular
    """
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return np.all(am.row <= am.col)


def istril(am):
    r"""
    Returns ``True`` is the sparse adjacency matrix is lower triangular
    """
    if am.shape[0] != am.shape[1]:
        print('Matrix is not square, triangularity is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return np.all(am.row >= am.col)


def istriangular(am):
    r"""
    Returns ``True`` is the sparse adjacency matrix is either upper or lower
    triangular
    """
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    return istril(am) or istriu(am)


def issymmetric(am):
    r"""
    A method to check if a square matrix is symmetric
    Returns ``True`` if the sparse adjacency matrix is symmetric
    """
    if am.shape[0] != am.shape[1]:
        logger.warning('Matrix is not square, symmetrical is irrelevant')
        return False
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    if istril(am) or istriu(am):
        return False
    # Compare am with its transpose, element wise
    sym = ((am != am.T).size) == 0
    return sym


def am_to_im(am):
    r"""
    Convert an adjacency matrix into an incidence matrix
    """
    if am.shape[0] != am.shape[1]:
        raise Exception('Adjacency matrices must be square')
    if am.format != 'coo':
        am = am.tocoo(copy=False)
    conn = np.vstack((am.row, am.col)).T
    row = conn[:, 0]
    data = am.data
    col = np.arange(np.size(am.data))
    if istriangular(am):
        row = np.append(row, conn[:, 1])
        data = np.append(data, data)
        col = np.append(col, col)
    im = sprs.coo.coo_matrix((data, (row, col)), (row.max()+1, col.max()+1))
    return im


def im_to_am(im):
    r"""
    Convert an incidence matrix into an adjacency matrix
    """
    if im.shape[0] == im.shape[1]:
        print('Warning: Received matrix is square which is unlikely')
    if im.shape[0] > im.shape[1]:
        print('Warning: Received matrix has more sites than bonds')
    if im.format != 'coo':
        im = im.tocoo(copy=False)


def tri_to_am(tri):
    r"""
    Given a Delaunay Triangulation object from Scipy's ``spatial`` module,
    converts to a sparse adjacency matrix network representation.

    Parameters
    ----------
    tri : Delaunay Triangulation Object
        This object is produced by ``scipy.spatial.Delaunay``

    Returns
    -------
    A sparse adjacency matrix in COO format.  The network is undirected
    and unweighted, so the adjacency matrix is upper-triangular and all the
    weights are set to 1.

    """
    # Create an empty list-of-list matrix
    lil = sprs.lil_matrix((tri.npoints, tri.npoints))
    # Scan through Delaunay triangulation to retrieve pairs
    indices, indptr = tri.vertex_neighbor_vertices
    for k in range(tri.npoints):
        lil.rows[k] = indptr[indices[k]:indices[k+1]]
    # Convert to coo format
    lil.data = lil.rows  # Just a dummy array to make things work properly
    coo = lil.tocoo()
    # Set weights to 1's
    coo.data = np.ones_like(coo.data)
    # Remove diagonal, and convert to csr remove duplicates
    am = sp.sparse.triu(A=coo, k=1, format='csr')
    # The convert back to COO and return
    am = am.tocoo()
    return am


def vor_to_am(vor):
    r"""
    Given a Voronoi tessellation object from Scipy's ``spatial`` module,
    converts to a sparse adjacency matrix network representation in COO format.

    Parameters
    ----------
    vor : Voronoi Tessellation object
        This object is produced by ``scipy.spatial.Voronoi``

    Returns
    -------
    A sparse adjacency matrix in COO format.  The network is undirected
    and unweighted, so the adjacency matrix is upper-triangular and all the
    weights are set to 1.

    """
    # Create adjacency matrix in lil format for quick matrix construction
    N = vor.vertices.shape[0]
    rc = [[], []]
    for ij in vor.ridge_dict.keys():
        row = vor.ridge_dict[ij].copy()
        # Make sure voronoi cell closes upon itself
        row.append(row[0])
        # Add connections to rc list
        rc[0].extend(row[:-1])
        rc[1].extend(row[1:])
    rc = np.vstack(rc).T
    # Make adj mat upper triangular
    rc = np.sort(rc, axis=1)
    # Remove any pairs with ends at infinity (-1)
    keep = ~np.any(rc == -1, axis=1)
    rc = rc[keep]
    data = np.ones_like(rc[:, 0])
    # Build adj mat in COO format
    M = N = np.amax(rc) + 1
    am = sprs.coo_matrix((data, (rc[:, 0], rc[:, 1])), shape=(M, N))
    # Remove diagonal, and convert to csr remove duplicates
    am = sp.sparse.triu(A=am, k=1, format='csr')
    # The convert back to COO and return
    am = am.tocoo()
    return am


def conns_to_am(conns, shape=None, force_triu=True, drop_diag=True,
                drop_dupes=True, drop_negs=True):
    r"""
    Converts a list of connections into a Scipy sparse adjacency matrix

    Parameters
    ----------
    conns : array_like, N x 2
        The list of site-to-site connections

    shape : list, optional
        The shape of the array.  If none is given then it is inferred from the
        maximum value in ``conns`` array.

    force_triu : boolean
        If True (default), then all connections are assumed undirected, and
        moved to the upper triangular portion of the array

    drop_diag : boolean
        If True (default), then connections from a site and itself are removed.

    drop_dupes : boolean
        If True (default), then all pairs of sites sharing multiple connections
        are reduced to a single connection.

    drop_negs : boolean
        If True (default), then all connections with one or both ends pointing
        to a negative number are removed.

    """
    if force_triu:  # Sort connections to [low, high]
        conns = np.sort(conns, axis=1)
    if drop_negs:  # Remove connections to -1
        keep = ~np.any(conns < 0, axis=1)
        conns = conns[keep]
    if drop_diag:  # Remove connections of [self, self]
        keep = np.where(conns[:, 0] != conns[:, 1])[0]
        conns = conns[keep]
    # Now convert to actual sparse array in COO format
    data = np.ones_like(conns[:, 0], dtype=int)
    if shape is None:
        N = conns.max() + 1
        shape = (N, N)
    am = sprs.coo_matrix((data, (conns[:, 0], conns[:, 1])), shape=shape)
    if drop_dupes:  # Convert to csr and back too coo
        am = am.tocsr()
        am = am.tocoo()
    # Perform one last check on adjacency matrix
    missing = np.where(np.bincount(conns.flatten()) == 0)[0]
    if np.size(missing) or np.any(am.col.max() < (shape[0] - 1)):
        warnings.warn('Some nodes are not connected to any bonds')
    return am


def isoutside(coords, shape):
    r"""
    Identifies points that lie outside the specified region.

    Parameters
    ----------
    domain_size : array_like
        The size and shape of the domain beyond which points should be
        trimmed. The argument is treated as follows:

        **sphere** : If a scalar or single element list is received, it's
        treated as the radius [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as
        the radius and height of a cylinder [r, z] whose central axis
        starts at [0, 0, 0] and extends in the positive z-direction.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner
        lies at [0, 0, 0].

    Returns
    -------
    An Np-long mask of True values indicating pores that lie outside the
    domain.

    """
    # Label external pores for trimming below
    if len(shape) == 1:  # Spherical
        # Find external points
        r = np.sqrt(np.sum(coords**2, axis=1))
        Ps = r > shape[0]
    elif len(shape) == 2:  # Cylindrical
        # Find external pores outside radius
        r = np.sqrt(np.sum(coords[:, [0, 1]]**2, axis=1))
        Ps = r > shape[0]
        # Find external pores above and below cylinder
        if shape[1] > 0:
            Ps = Ps + (coords[:, 2] > shape[1])
            Ps = Ps + (coords[:, 2] < 0)
        else:
            pass
    elif len(shape) == 3:  # Rectilinear
        shape = np.array(shape, dtype=float)
        try:
            lo_lim = shape[:, 0]
            hi_lim = shape[:, 1]
        except IndexError:
            lo_lim = np.array([0, 0, 0])
            hi_lim = shape
        Ps1 = np.any(coords > hi_lim, axis=1)
        Ps2 = np.any(coords < lo_lim, axis=1)
        Ps = Ps1 + Ps2
    return Ps


def ispercolating(am, inlets, outlets, mode='site'):
    r"""
    Determines if a percolating clusters exists in the network spanning
    the given inlet and outlet sites

    Parameters
    ----------
    am : adjacency_matrix
        The adjacency matrix with the ``data`` attribute indicating
        if a bond is occupied or not

    inlets : array_like
        An array of indices indicating which sites are part of the inlets

    outlets : array_like
        An array of indices indicating which sites are part of the outlets

    mode : string
        Indicates which type of percolation to apply, either `'site'` or
        `'bond'`

    """
    if am.format != 'coo':
        am = am.to_coo()
    ij = np.vstack((am.col, am.row)).T
    if mode.startswith('site'):
        occupied_sites = np.zeros(shape=am.shape[0], dtype=bool)
        occupied_sites[ij[am.data].flatten()] = True
        clusters = site_percolation(ij, occupied_sites)
    elif mode.startswith('bond'):
        occupied_bonds = am.data
        clusters = bond_percolation(ij, occupied_bonds)
    ins = np.unique(clusters.sites[inlets])
    if ins[0] == -1:
        ins = ins[1:]
    outs = np.unique(clusters.sites[outlets])
    if outs[0] == -1:
        outs = outs[1:]
    hits = np.in1d(ins, outs)
    return np.any(hits)


def remove_isolated_clusters(labels, inlets):
    r"""
    Finds cluster labels not attached to the inlets, and sets them to
    unoccupied (-1)

    Parameters
    ----------
    labels : tuple of site and bond labels
        This information is provided by the ``site_percolation`` or
        ``bond_percolation`` functions

    inlets : array_like
        A list of which sites are inlets.  Can be a boolean mask or an
        array of indices.

    Returns
    -------
    A tuple containing a list of site and bond labels, with all clusters
    not connected to the inlet sites set to not occupied.

    """
    # Identify clusters of invasion sites
    inv_clusters = np.unique(labels.sites[inlets])
    # Remove cluster numbers == -1, if any
    inv_clusters = inv_clusters[inv_clusters >= 0]
    # Find all pores in invading clusters
    p_invading = np.in1d(labels.sites, inv_clusters)
    labels.sites[~p_invading] = -1
    t_invading = np.in1d(labels.bonds, inv_clusters)
    labels.bonds[~t_invading] = -1
    return labels


def site_percolation(ij, occupied_sites):
    r"""
    Calculates the site and bond occupancy status for a site percolation
    process given a list of occupied sites.

    Parameters
    ----------
    ij : array_like
        An N x 2 array of [site_A, site_B] connections.  If two connected
        sites are both occupied they are part of the same cluster, as it
        the bond connecting them.

    occupied_sites : boolean
        A list indicating whether sites are occupied or not

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates unoccupied.

    Notes
    -----
    The ``connected_components`` function of scipy.csgraph will give ALL
    sites a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    from collections import namedtuple
    import scipy.stats as spst

    Np = np.size(occupied_sites)
    occupied_bonds = np.all(occupied_sites[ij], axis=1)
    adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    clusters[~occupied_sites] = -1
    s_labels = spst.rankdata(clusters + 1, method="dense") - 1
    if np.any(~occupied_sites):
        s_labels -= 1
    b_labels = np.amin(s_labels[ij], axis=1)
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def bond_percolation(ij, occupied_bonds):
    r"""
    Calculates the site and bond occupancy status for a bond percolation
    process given a list of occupied bonds.

    Parameters
    ----------
    ij : array_like
        An N x 2 array of [site_A, site_B] connections.  A site is
        considered occupied if any of it's connecting bonds are occupied.

    occupied_bonds: boolean
        A list indicating whether a bond is occupied or not

    Returns
    -------
    A tuple containing a list of site and bond labels, indicating which
    cluster each belongs to.  A value of -1 indicates uninvaded.

    Notes
    -----
    The ``connected_components`` function of scipy.csgraph will give ALL
    sites a cluster number whether they are occupied or not, so this
    function essentially adjusts the cluster numbers to represent a
    percolation process.

    """
    from collections import namedtuple
    Np = np.amax(ij) + 1
    # Find occupied sites based on occupied bonds
    # (the following 2 lines are not needed but worth keeping for future ref)
    # occupied_sites = np.zeros([Np, ], dtype=bool)
    # np.add.at(occupied_sites, ij[occupied_bonds].flatten(), True)
    adj_mat = sprs.csr_matrix((occupied_bonds, (ij[:, 0], ij[:, 1])),
                              shape=(Np, Np))
    adj_mat.eliminate_zeros()
    clusters = csgraph.connected_components(csgraph=adj_mat, directed=False)[1]
    # Clusters of size 1 only occur if all a site's bonds are uninvaded
    valid_clusters = np.bincount(clusters) > 1
    mapping = -np.ones(shape=(clusters.max()+1, ), dtype=int)
    mapping[valid_clusters] = np.arange(0, valid_clusters.sum())
    s_labels = mapping[clusters]
    # Bond inherit the cluster number of its connected sites
    b_labels = np.amin(s_labels[ij], axis=1)
    # Set bond cluster to -1 if not actually occupied
    b_labels[~occupied_bonds] = -1
    tup = namedtuple('cluster_labels', ('sites', 'bonds'))
    return tup(s_labels, b_labels)


def rotate_coords(network, a=0, b=0, c=0, R=None):
    r"""
    Rotates coordinates a given amount about each axis

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose pore coordinates should be transformed
    a : scalar
        The amount in degrees to rotate about the x-axis
    b : scalar
        The amount in degrees to rotate about the y-axis
    c : scalar
        The amount in degrees to rotate about the z-axis
    R : array_like
        Rotation matrix.  Must be a 3-by-3 matrix since pore coordinates are
        always in 3D.  If this is given then the other individual arguments
        are ignored.

    See Also
    --------
    rotate_coords

    Notes
    -----
    It is possible to rotate about anyof the three axes by specifying ``a``,
    ``b``, and/or ``c``.  In this case each rotation is applied in sequence.

    """
    if R is None:
        if a:
            R = np.array([[1, 0, 0],
                          [0, np.cos(np.deg2rad(a)), -np.sin(np.deg2rad(a))],
                          [0, np.sin(np.deg2rad(a)), np.cos(np.deg2rad(a))]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
        if b:
            R = np.array([[np.cos(np.deg2rad(b)), 0, -np.sin(np.deg2rad(b))],
                          [0, 1, 0],
                          [np.sin(np.deg2rad(b)), 0, np.cos(np.deg2rad(b))]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
        if c:
            R = np.array([[np.cos(np.deg2rad(c)), -np.sin(np.deg2rad(c)), 0],
                          [np.sin(np.deg2rad(c)), np.cos(np.deg2rad(c)), 0],
                          [0, 0, 1]])
            network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                                  axes=(1, 1))
    else:
        network['pore.coords'] = np.tensordot(network['pore.coords'], R,
                                              axes=(1, 1))


def shear_coords(network, ay=0, az=0, bx=0, bz=0, cx=0, cy=0, S=None):
    r"""
    Shears the coordinates a given amount about along axis

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose pore coordinates should be transformed
    ay : scalar
        The factor by which to shear along the x-axis as a function of y
    az : scalar
        The factor by which to shear along the x-axis as a function of z
    bx : scalar
        The factor by which to shear along the y-axis as a function of x
    bz : scalar
        The factor by which to shear along the y-axis as a function of z
    cx : scalar
        The factor by which to shear along the z-axis  as a function of x
    cy : scalar
        The factor by which to shear along the z-axis as a function of y
    S : array_like
        The shear matrix.  Must be a 3-by-3 matrix since pore coordinates are
        always in 3D.  If this is given then the other individual arguments
        are ingnored.

    See Also
    --------
    rotate_coords

    Notes
    -----
    The shear along the i *th* -axis is given as i\* = i + aj.  This means
    the new i coordinate is the old one plus some linear factor *a* in the
    j *th* direction.

    The values of ``a``, ``b``, and ``c`` are essentially the inverse of the
    slope to be formed by the neighboring layers of sheared pores.  A value of
    0 means no shear, and neighboring points are stacked directly on top of
    each other; a value of 1 means they form a 45 degree diagonal, and so on.

    If ``S`` is given, then is should be of the form:

    ::

        S = [[1 , ay, az],
             [bx, 1 , bz],
             [cx, cy, 1 ]]

        where any of the off-diagonal components can be 0

    """
    if S is None:
        S = np.array([[1, ay, az],
                      [bx, 1, bz],
                      [cx, cy, 1]])
    network['pore.coords'] = (S@network['pore.coords'].T).T


def trim(network, pores=[], throats=[]):
    '''
    Remove pores or throats from the network.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network from which pores or throats should be removed

    pores (or throats) : array_like
        The indices of the of the pores or throats to be removed from the
        network.

    Notes
    -----
    This is an in-place operation, meaning the received Network object will
    be altered directly.


    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn.Np
    125
    >>> pn.Nt
    300
    >>> op.topotools.trim(network=pn, pores=[1])
    >>> pn.Np
    124
    >>> pn.Nt
    296

    '''
    pores = network._parse_indices(pores)
    throats = network._parse_indices(throats)
    Pkeep = np.copy(network['pore.all'])
    Tkeep = np.copy(network['throat.all'])
    if np.size(pores) > 0:
        Pkeep[pores] = False
        if not np.any(Pkeep):
            raise Exception('Cannot delete ALL pores')
        # # Performing customized find_neighbor_throats which is much faster, but
        # # not general for other types of queries
        # temp = np.in1d(network['throat.conns'].flatten(), pores)
        # temp = np.reshape(temp, (network.Nt, 2))
        # Ts = np.any(temp, axis=1)
        # Ts = network.Ts[Ts]
        Ts = network.find_neighbor_throats(pores=~Pkeep, mode='union')
        if len(Ts) > 0:
            Tkeep[Ts] = False
    if np.size(throats) > 0:
        Tkeep[throats] = False
        # The following IF catches the special case of deleting ALL throats
        # It removes all throat props, adds 'all', and skips rest of function
        if not np.any(Tkeep):
            logger.info('Removing ALL throats from network')
            for item in network.keys():
                if item.split('.')[0] == 'throat':
                    del network[item]
            network['throat.all'] = np.array([], ndmin=1)
            return

    # Temporarily store throat conns and pore map for processing later
    Np_old = network.Np
    Nt_old = network.Nt
    Pkeep_inds = np.where(Pkeep)[0]
    Tkeep_inds = np.where(Tkeep)[0]
    Pmap = np.ones((network.Np,), dtype=int)*-1
    tpore1 = network['throat.conns'][:, 0]
    tpore2 = network['throat.conns'][:, 1]

    # Delete specified pores and throats from all objects
    for obj in network.project[::-1]:
        if (obj.Np == Np_old) and (obj.Nt == Nt_old):
            Ps = Pkeep_inds
            Ts = Tkeep_inds
        else:
            Ps = obj.map_pores(pores=Pkeep, origin=network)
            Ts = obj.map_throats(throats=Tkeep, origin=network)
        for key in list(obj.keys()):
            temp = obj.pop(key)
            if key.split('.')[0] == 'throat':
                obj.update({key: temp[Ts]})
            if key.split('.')[0] == 'pore':
                obj.update({key: temp[Ps]})

    # Remap throat connections
    Pmap[Pkeep] = np.arange(0, np.sum(Pkeep))
    Tnew1 = Pmap[tpore1[Tkeep]]
    Tnew2 = Pmap[tpore2[Tkeep]]
    network.update({'throat.conns': np.vstack((Tnew1, Tnew2)).T})

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def extend(network, pore_coords=[], throat_conns=[], labels=[]):
    r'''
    Add individual pores and/or throats to the network from a list of coords
    or conns.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network to which pores or throats should be added

    pore_coords : array_like
        The coordinates of the pores to add

    throat_conns : array_like
        The throat connections to add

    labels : string, or list of strings, optional
        A list of labels to apply to the new pores and throats

    Notes
    -----
    This needs to be enhanced so that it increases the size of all pore
    and throat props and labels on ALL associated Phase objects.  At the
    moment it throws an error is there are any associated Phases.

    This is an in-place operation, meaning the received Network object will
    be altered directly.

    '''
    if len(network.project.phases()) > 0:
        raise Exception('Project has active Phases, copy network to a new '
                        + 'project and try again')

    Np_old = network.num_pores()
    Nt_old = network.num_throats()
    Np = Np_old + int(np.size(pore_coords)/3)
    Nt = Nt_old + int(np.size(throat_conns)/2)
    # Adjust 'all' labels
    del network['pore.all'], network['throat.all']
    network['pore.all'] = np.ones((Np,), dtype=bool)
    network['throat.all'] = np.ones((Nt,), dtype=bool)
    # Add coords and conns
    if np.size(pore_coords) > 0:
        coords = np.vstack((network['pore.coords'], pore_coords))
        network['pore.coords'] = coords
    if np.size(throat_conns) > 0:
        conns = np.vstack((network['throat.conns'], throat_conns))
        network['throat.conns'] = conns

    # Increase size of any prop or label arrays aready on Network
    for item in list(network.keys()):
        if item.split('.')[1] not in ['coords', 'conns', 'all', '_id']:
            N = network._count(element=item.split('.')[0])
            arr = network.pop(item)
            s = arr.shape
            if arr.dtype == bool:
                network[item] = np.zeros(shape=(N, *s[1:]), dtype=bool)
            else:
                network[item] = np.ones(shape=(N, *s[1:]), dtype=float)*sp.nan
            # This is a temporary work-around until I learn to handle 2+ dims
            network[item][:arr.shape[0]] = arr

    # Apply labels, if supplied
    if labels != []:
        # Convert labels to list if necessary
        if type(labels) is str:
            labels = [labels]
        for label in labels:
            # Remove pore or throat from label, if present
            label = label.split('.')[-1]
            if np.size(pore_coords) > 0:
                Ps = sp.r_[Np_old:Np]
                if 'pore.'+label not in network.labels():
                    network['pore.'+label] = False
                network['pore.'+label][Ps] = True
            if np.size(throat_conns) > 0:
                Ts = sp.r_[Nt_old:Nt]
                if 'throat.'+label not in network.labels():
                    network['throat.'+label] = False
                network['throat.'+label][Ts] = True

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def reduce_coordination(network, z):
    r"""
    """
    # Find minimum spanning tree using random weights
    am = network.create_adjacency_matrix(weights=np.random.rand(network.Nt),
                                         triu=False)
    mst = csgraph.minimum_spanning_tree(am, overwrite=True)
    mst = mst.tocoo()

    # Label throats on spanning tree to avoid deleting them
    Ts = network.find_connecting_throat(mst.row, mst.col)
    Ts = np.hstack(Ts)
    network['throat.mst'] = False
    network['throat.mst'][Ts] = True

    # Trim throats not on the spanning tree to acheive desired coordination
    Ts = sp.random.permutation(network.throats('mst', mode='nor'))
    Ts = Ts[:int(network.Nt - network.Np*(z/2))]
    trim(network=network, throats=Ts)


def label_faces(network, tol=0.0, label='surface'):
    r"""
    Finds pores on the surface of the network and labels them according to
    whether they are on the *top*, *bottom*, etc.  This function assumes the
    network is cubic in shape (i.e. with six flat sides)

    Parameters
    ----------
    network : OpenPNM Network object
        The network to apply the labels

    tol : scalar
        The tolerance for defining what counts as a surface pore, which is
        specifically meant for random networks.  All pores with ``tol`` of
        the maximum or minimum along each axis are counts as pores.  The
        default is 0.

    label : string
        An identifying label to isolate the pores on the faces of the network.
        The default is 'surface'.  Surface pores can be found using
        ``find_surface_pores``.

    """
    label = label.split('.', 1)[-1]
    if 'pore.'+label not in network.labels():
        find_surface_pores(network, label=label)
    Psurf = network['pore.'+label]
    crds = network['pore.coords']
    xmin, xmax = np.amin(crds[:, 0]), np.amax(crds[:, 0])
    xspan = xmax - xmin
    ymin, ymax = np.amin(crds[:, 1]), np.amax(crds[:, 1])
    yspan = ymax - ymin
    zmin, zmax = np.amin(crds[:, 2]), np.amax(crds[:, 2])
    zspan = zmax - zmin
    dims = dimensionality(network)
    if dims[0]:
        network['pore.front'] = (crds[:, 0] <= (xmin + tol*xspan)) * Psurf
        network['pore.back'] = (crds[:, 0] >= (xmax - tol*xspan)) * Psurf
    if dims[1]:
        network['pore.left'] = (crds[:, 1] <= (ymin + tol*yspan)) * Psurf
        network['pore.right'] = (crds[:, 1] >= (ymax - tol*yspan)) * Psurf
    if dims[2]:
        network['pore.top'] = (crds[:, 2] >= (zmax - tol*zspan)) * Psurf
        network['pore.bottom'] = (crds[:, 2] <= (zmin + tol*zspan)) * Psurf


def find_surface_pores(network, markers=None, label='surface'):
    r"""
    Find the pores on the surface of the domain by performing a Delaunay
    triangulation between the network pores and some external ``markers``. All
    pores connected to these external marker points are considered surface
    pores.

    Parameters
    ----------
    network: OpenPNM Network Object
        The network for which the surface pores are to be found

    markers: array_like
        3 x N array of the marker coordinates to use in the triangulation.  The
        labeling is performed in one step, so all points are added, and then
        any pores connected to at least one marker is given the provided label.
        By default, this function will automatically generate 6 points outside
        each axis of the network domain.

        Users may wish to specify a single external marker point and provide an
        appropriate label in order to identify specific faces.  For instance,
        the marker may be *above* the domain, and the label might be
        'top_surface'.

    label : string
        The label to apply to the pores.  The default is 'surface'.

    Notes
    -----
    This function does not check whether the given markers actually lie outside
    the domain, allowing the labeling of *internal* sufaces.

    If this method fails to mark some surface pores, consider sending more
    markers on each face.

    Examples
    --------
    >>> import openpnm as op
    >>> net = op.network.Cubic(shape=[5, 5, 5])
    >>> op.topotools.find_surface_pores(network=net)
    >>> net.num_pores('surface')
    98

    When cubic networks are created, the surfaces are already labeled:

    >>> net.num_pores(['top','bottom', 'left', 'right', 'front','back'])
    98

    """
    import scipy.spatial as sptl
    dims = dimensionality(network)
    coords = network['pore.coords'][:, dims]
    if markers is None:
        # normalize coords to a 1 unit cube centered on origin
        coords -= np.amin(coords, axis=0)
        coords /= np.amax(coords, axis=0)
        coords -= 0.5
        npts = max((network.Np/10, 100))
        if sum(dims) == 1:
            network['pore.'+label] = True
            return
        if sum(dims) == 2:
            r = 0.75
            theta = sp.linspace(0, 2*sp.pi, int(npts), dtype=float)
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            markers = np.vstack((x, y)).T
        if sum(dims) == 3:
            r = 1.00
            indices = np.arange(0, int(npts), dtype=float) + 0.5
            phi = np.arccos(1 - 2*indices/npts)
            theta = sp.pi * (1 + 5**0.5) * indices
            x = r*np.cos(theta) * np.sin(phi)
            y = r*np.sin(theta) * np.sin(phi)
            z = r*np.cos(phi)
            markers = np.vstack((x, y, z)).T
    else:
        if sum(dims) == 1:
            pass
        if sum(dims) == 2:
            markers = np.atleast_2d(markers)
            if markers.shape[1] != 2:
                raise Exception('Network appears planar, so markers must be 2D')
        if sum(dims) == 3:
            markers = np.atleast_2d(markers)
            if markers.shape[1] != 3:
                raise Exception('Markers must be 3D for this network')
    pts = np.vstack((coords, markers))
    tri = sptl.Delaunay(pts, incremental=False)
    (indices, indptr) = tri.vertex_neighbor_vertices
    for k in range(network.Np, tri.npoints):
        neighbors = indptr[indices[k]:indices[k+1]]
        inds = np.where(neighbors < network.Np)
        neighbors = neighbors[inds]
        if 'pore.'+label not in network.keys():
            network['pore.'+label] = False
        network['pore.'+label][neighbors] = True


def dimensionality(network):
    r"""
    Checks the dimensionality of the network

    Parameters
    ----------
    network : OpenPNM Network object
        The network whose dimensionality is to be checked

    Returns
    -------
    Returns an 3-by-1 array containing ``True`` for each axis that contains
    multiple values, indicating that the pores are spatially distributed
    in that direction.
    """
    xyz = network["pore.coords"]
    xyz_unique = [np.unique(xyz[:, i]) for i in range(3)]
    return np.array([elem.size != 1 for elem in xyz_unique])


def clone_pores(network, pores, labels=['clone'], mode='parents'):
    r"""
    Clones the specified pores and adds them to the network

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object to which the new pores are to be added

    pores : array_like
        List of pores to clone

    labels : string, or list of strings
        The labels to apply to the clones, default is 'clone'

    mode : string
        Controls the connections between parents and clones.  Options are:

        - 'parents': (Default) Each clone is connected only to its parent
        - 'siblings': Clones are only connected to each other in the same
                      manner as parents were connected
        - 'isolated': No connections between parents or siblings
    """
    if len(network.project.geometries()) > 0:
        logger.warning('Network has active Geometries, new pores must be '
                       'assigned a Geometry')
    if len(network.project.phases()) > 0:
        raise Exception('Network has active Phases, cannot proceed')

    if type(labels) == str:
        labels = [labels]
    network._parse_indices(pores)
    Np = network.Np
    Nt = network.Nt
    # Clone pores
    parents = np.array(pores, ndmin=1)
    pcurrent = network['pore.coords']
    pclone = pcurrent[pores, :]
    pnew = np.concatenate((pcurrent, pclone), axis=0)
    Npnew = np.shape(pnew)[0]
    clones = np.arange(Np, Npnew)
    # Add clone labels to network
    for item in labels:
        network['pore.'+item] = False
        network['throat.'+item] = False
    # Add connections between parents and clones
    if mode == 'parents':
        tclone = np.vstack((parents, clones)).T
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'siblings':
        ts = network.find_neighbor_throats(pores=pores, mode='xnor')
        tclone = network['throat.conns'][ts] + network.num_pores()
        extend(network=network, pore_coords=pclone, throat_conns=tclone)
    if mode == 'isolated':
        extend(network=network, pore_coords=pclone)
    # Apply provided labels to cloned pores
    for item in labels:
        network['pore.'+item][network.pores('all') >= Np] = True
        network['throat.'+item][network.throats('all') >= Nt] = True

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def merge_networks(network, donor=[]):
    r"""
    Combine multiple networks into one without doing any topological
    manipulations (such as stiching nearby pores to each other).

    Parameters
    ----------
    network : OpenPNM Network Object
        The network to which all the other networks should be added.

    donor : OpenPNM Network Object or list of Objects
        The network object(s) to add to the given network

    Notes
    -----
    This methods does *not* attempt to stitch the networks topologically.

    See Also
    --------
    extend
    trim
    stitch

    """
    if type(donor) == list:
        donors = donor
    else:
        donors = [donor]

    for donor in donors:
        network['pore.coords'] = np.vstack((network['pore.coords'],
                                            donor['pore.coords']))
        network['throat.conns'] = np.vstack((network['throat.conns'],
                                             donor['throat.conns']
                                             + network.Np))
        p_all = np.ones((np.shape(network['pore.coords'])[0],), dtype=bool)
        t_all = np.ones((np.shape(network['throat.conns'])[0],), dtype=bool)
        network.update({'pore.all': p_all})
        network.update({'throat.all': t_all})
        for key in set(network.keys()).union(set(donor.keys())):
            if key.split('.')[1] not in ['conns', 'coords', '_id', 'all']:
                if key in network.keys():
                    pop_flag = False
                    # If key not on donor add it first with dummy values to
                    # simplify merging later
                    if key not in donor.keys():
                        logger.debug('Adding ' + key + ' to donor')
                        if network[key].dtype == bool:  # Deal with labels
                            donor[key] = False
                        else:  # Deal with numerical data
                            element = key.split('.')[0]
                            shape = list(network[key].shape)
                            N = donor._count(element)
                            shape[0] = N
                            donor[key] = sp.empty(shape=shape)*sp.nan
                        pop_flag = True
                    # Then merge it with existing array on network
                    if len(network[key].shape) == 1:
                        temp = np.hstack((network[key], donor[key]))
                    else:
                        temp = np.vstack((network[key], donor[key]))
                    network[key] = temp
                    if pop_flag:
                        donor.pop(key, None)
                else:
                    # If key not on network add it first
                    logger.debug('Adding ' + key + ' to network')
                    if donor[key].dtype == bool:
                        network[key] = False
                    else:
                        data_shape = list(donor[key].shape)
                        pore_prop = True if key.split(".")[0] == "pore" else False
                        data_shape[0] = network.Np if pore_prop else network.Nt
                        network[key] = sp.empty(data_shape) * sp.nan
                    # Then append donor values to network
                    s = np.shape(donor[key])[0]
                    network[key][-s:] = donor[key]

    # Clear adjacency and incidence matrices which will be out of date now
    network._am.clear()
    network._im.clear()


def stitch(network, donor, P_network, P_donor, method='nearest',
           len_max=sp.inf, label_suffix='', label_stitches='stitched'):
    r'''
    Stitches a second a network to the current network.

    Parameters
    ----------
    networK : OpenPNM Network Object
        The Network to which to donor Network will be attached

    donor : OpenPNM Network Object
        The Network to stitch on to the current Network

    P_network : array_like
        The pores on the current Network

    P_donor : array_like
        The pores on the donor Network

    label_suffix : str or None
        Some text to append to each label in the donor Network before
        inserting them into the recipient.  The default is to append no
        text, but a common option would be to append the donor Network's
        name. To insert none of the donor labels, use ``None``.

    label_stitches : str or list of strings
        The label to apply to the newly created 'stitch' throats.  The
        defaul is 'stitched'.  If performing multiple stitches in a row it
        might be helpful to the throats created during each step uniquely
        for later identification.

    len_max : float
        Set a length limit on length of new throats

    method : string (default = 'nearest')
        The method to use when making pore to pore connections. Options are:

        - 'radius' : Connects each pore on the recipient network to the
                     nearest pores on the donor network, within ``len_max``
        - 'nearest' : Connects each pore on the recipienet network to the
                      nearest pore on the donor network.

    Notes
    -----
    Before stitching it is necessary to translate the pore coordinates of
    one of the Networks so that it is positioned correctly relative to the
    other.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn2 = op.network.Cubic(shape=[5, 5, 5])
    >>> [pn.Np, pn.Nt]
    [125, 300]
    >>> [pn2.Np, pn2.Nt]
    [125, 300]
    >>> pn2['pore.coords'][:, 2] += 5.0
    >>> op.topotools.stitch(network=pn, donor=pn2, P_network=pn.pores('top'),
    ...                     P_donor=pn2.pores('bottom'), method='radius',
    ...                     len_max=1.0)
    >>> [pn.Np, pn.Nt]
    [250, 625]

    '''
    # Ensure Networks have no associated objects yet
    if (len(network.project) > 1) or (len(donor.project) > 1):
        raise Exception('Cannot stitch a Network with active objects')
    if isinstance(label_stitches, str):
        label_stitches = [label_stitches]
    for s in label_stitches:
        if s not in network.keys():
            network['throat.' + s] = False
    # Get the initial number of pores and throats
    N_init = {}
    N_init['pore'] = network.Np
    N_init['throat'] = network.Nt
    if method == 'nearest':
        P1 = P_network
        P2 = P_donor + N_init['pore']  # Increment pores on donor
        C1 = network['pore.coords'][P_network]
        C2 = donor['pore.coords'][P_donor]
        D = sp.spatial.distance.cdist(C1, C2)
        [P1_ind, P2_ind] = np.where(D == D.min(axis=0))
        conns = np.vstack((P1[P1_ind], P2[P2_ind])).T
    elif method == 'radius':
        P1 = P_network
        P2 = P_donor + N_init['pore']  # Increment pores on donor
        C1 = network['pore.coords'][P_network]
        C2 = donor['pore.coords'][P_donor]
        D = sp.spatial.distance.cdist(C1, C2)
        [P1_ind, P2_ind] = np.where(D <= len_max)
        conns = np.vstack((P1[P1_ind], P2[P2_ind])).T
    else:
        raise Exception('<{}> method not supported'.format(method))

    # Enter donor's pores into the Network
    extend(network=network, pore_coords=donor['pore.coords'])

    # Enter donor's throats into the Network
    extend(network=network, throat_conns=donor['throat.conns'] + N_init['pore'])

    # Trim throats that are longer then given len_max
    # C1 = network['pore.coords'][conns[:, 0]]
    # C2 = network['pore.coords'][conns[:, 1]]
    # L = np.sum((C1 - C2)**2, axis=1)**0.5
    # conns = conns[L <= len_max]

    # Add donor labels to recipient network
    if label_suffix is not None:
        if label_suffix != '':
            label_suffix = '_'+label_suffix
        for label in donor.labels():
            element = label.split('.')[0]
            locations = np.where(network._get_indices(element)
                                 >= N_init[element])[0]
            if label + label_suffix not in network.keys():
                network[label + label_suffix] = False
            network[label+label_suffix][locations] = donor[label]

    # Add the new stitch throats to the Network
    extend(network=network, throat_conns=conns, labels=label_stitches)

    # Remove donor from Workspace, if present
    # This check allows for the reuse of a donor Network multiple times
    for sim in list(ws.values()):
        if donor in sim:
            del ws[sim.name]


def connect_pores(network, pores1, pores2, labels=[], add_conns=True):
    r'''
    Returns the possible connections between two groups of pores, and optionally
    makes the connections.

    See ``Notes`` for advanced usage.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores1 : array_like
        The first group of pores on the network

    pores2 : array_like
        The second group of pores on the network

    labels : list of strings
        The labels to apply to the new throats.  This argument is only needed
        if ``add_conns`` is True.

    add_conns : bool
        Indicates whether the connections should be added to the supplied
        network (default is True).  Otherwise, the connections are returned
        as an Nt x 2 array that can be passed directly to ``extend``.

    Notes
    -----
    (1) The method also works if ``pores1`` and ``pores2`` are list of lists,
    in which case it consecutively connects corresponding members of the two
    lists in a 1-to-1 fashion. Example: pores1 = [[0, 1], [2, 3]] and
    pores2 = [[5], [7, 9]] leads to creation of the following connections:

    ::

        0 --> 5     2 --> 7     3 --> 7
        1 --> 5     2 --> 9     3 --> 9

    (2) If you want to use the batch functionality, make sure that each element
    within ``pores1`` and ``pores2`` are of type list or ndarray.

    (3) It creates the connections in a format which is acceptable by
    the default OpenPNM connections ('throat.conns') and either adds them to
    the network or returns them.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> pn.Nt
    300
    >>> op.topotools.connect_pores(network=pn, pores1=[22, 32],
    ...                            pores2=[16, 80, 68])
    >>> pn.Nt
    306
    >>> pn['throat.conns'][300:306]
    array([[16, 22],
           [22, 80],
           [22, 68],
           [16, 32],
           [32, 80],
           [32, 68]])

    '''
    # Assert that `pores1` and `pores2` are list of lists
    try:
        len(pores1[0])
    except (TypeError, IndexError):
        pores1 = [pores1]
    try:
        len(pores2[0])
    except (TypeError, IndexError):
        pores2 = [pores2]

    if len(pores1) != len(pores2):
        raise Exception('Running in batch mode! pores1 and pores2 must be'
                        + ' of the same length.')

    arr1, arr2 = [], []
    for ps1, ps2 in zip(pores1, pores2):
        size1 = np.size(ps1)
        size2 = np.size(ps2)
        arr1.append(sp.repeat(ps1, size2))
        arr2.append(np.tile(ps2, size1))
    conns = np.vstack([np.concatenate(arr1), np.concatenate(arr2)]).T
    if add_conns:
        extend(network=network, throat_conns=conns, labels=labels)
    else:
        return conns


def find_pore_to_pore_distance(network, pores1=None, pores2=None):
    r'''
    Find the distance between all pores on set 1 to each pore in set 2

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object containing the pore coordinates

    pores1 : array_like
        The pore indices of the first set

    pores2 : array_Like
        The pore indices of the second set.  It's OK if these indices are
        partially or completely duplicating ``pores``.

    Returns
    -------
    A distance matrix with ``len(pores1)`` rows and ``len(pores2)`` columns.
    The distance between pore *i* in ``pores1`` and *j* in ``pores2`` is
    located at *(i, j)* and *(j, i)* in the distance matrix.

    Notes
    -----
    This function computes and returns a distance matrix, which is a dense
    matrix of size Np_1 by Np_2, so can get large.  For distances between
    larger sets a KD-tree approach would be better, which is available in
    ``scipy.spatial``.

    '''
    from scipy.spatial.distance import cdist
    p1 = np.array(pores1, ndmin=1)
    p2 = np.array(pores2, ndmin=1)
    coords = network['pore.coords']
    return cdist(coords[p1], coords[p2])


def subdivide(network, pores, shape, labels=[]):
    r'''
    It trim the pores and replace them by cubic networks with the sent shape.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores : array_like
        The first group of pores to be replaced

    shape : array_like
        The shape of cubic networks in the target locations

    Notes
    -----
    - It works only for cubic networks.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 6, 5], spacing=0.001)
    >>> pn.Np
    150
    >>> nano_pores = [2, 13, 14, 15]
    >>> op.topotools.subdivide(network=pn, pores=nano_pores, shape=[4, 7, 3],
    ...                        labels='nano')
    >>> pn.Np
    482

    '''
    mro = network._mro()
    if 'Cubic' not in mro:
        raise Exception('Subdivide is only supported for Cubic Networks')
    from openpnm.network import Cubic
    pores = network._parse_indices(pores)

    # Checks to find boundary pores in the selected pores
    if 'pore.boundary' in network.labels():
        if (np.in1d(pores, network.pores('boundary'))).any():
            raise Exception('boundary pores cannot be subdivided!')
    if not hasattr(network, '_subdivide_flag'):
        network._subdivide_flag = True
    else:
        raise Exception('The network has subdivided pores, so the method \
                         does not support another subdivision')
    # Assigning right shape and division
    if np.size(shape) != 2 and np.size(shape) != 3:
        raise Exception('Subdivide not implemented for Networks other than 2D \
                         and 3D')
    elif np.size(shape) == 3 and 1 not in shape:
        div = np.array(shape, ndmin=1)
        single_dim = None
    else:
        single_dim = np.where(np.array(network.shape) == 1)[0]
        if np.size(single_dim) == 0:
            single_dim = None
        if np.size(shape) == 3:
            div = np.array(shape, ndmin=1)
        else:
            div = np.zeros(3, dtype=sp.int32)
            if single_dim is None:
                dim = 2
            else:
                dim = single_dim
            div[dim] = 1
            div[-np.array(div, ndmin=1, dtype=bool)] = np.array(shape, ndmin=1)

    # Creating small network and handling labels
    networkspacing = network.spacing
    new_netspacing = networkspacing/div
    new_net = Cubic(shape=div, spacing=new_netspacing)
    main_labels = ['left', 'right', 'front', 'back', 'top', 'bottom']
    if single_dim is not None:
        label_groups = np.array([['front', 'back'],
                                 ['left', 'right'],
                                 ['top', 'bottom']])
        non_single_labels = label_groups[np.array([0, 1, 2]) != single_dim]
    for l in main_labels:
        new_net['pore.surface_' + l] = False
        network['pore.surface_' + l] = False
        if single_dim is None:
            new_net['pore.surface_' + l][new_net.pores(labels=l)] = True
        else:
            for ind in [0, 1]:
                loc = (non_single_labels[ind] == l)
                temp_pores = new_net.pores(non_single_labels[ind][loc])
                new_net['pore.surface_' + l][temp_pores] = True

    old_coords = np.copy(new_net['pore.coords'])
    if labels == []:
        labels = ['pore.subdivided_' + new_net.name]
    for P in pores:
        # Shifting the new network to the right location and attaching it to
        # the main network
        shift = network['pore.coords'][P] - networkspacing/2
        new_net['pore.coords'] += shift
        Pn = network.find_neighbor_pores(pores=P)
        try:
            Pn_new_net = network.pores(labels)
        except KeyError:
            Pn_new_net = []
        Pn_old_net = Pn[~np.in1d(Pn, Pn_new_net)]
        Np1 = network.Np
        extend(pore_coords=new_net['pore.coords'],
               throat_conns=new_net['throat.conns'] + Np1,
               labels=labels, network=network)

        # Moving the temporary labels to the big network
        for l in main_labels:
            network['pore.surface_'+l][Np1:] = new_net['pore.surface_'+l]

        # Stitching the old pores of the main network to the new extended pores
        surf_pores = network.pores('surface_*')
        surf_coord = network['pore.coords'][surf_pores]
        for neighbor in Pn:
            neighbor_coord = network['pore.coords'][neighbor]
            dist = [round(sp.inner(neighbor_coord-x, neighbor_coord-x),
                          20) for x in surf_coord]
            nearest_neighbor = surf_pores[dist == np.amin(dist)]
            if neighbor in Pn_old_net:
                coplanar_labels = network.labels(pores=nearest_neighbor)
                new_neighbors = network.pores(coplanar_labels,
                                              mode='and')
                # This might happen to the edge of the small network
                if np.size(new_neighbors) == 0:
                    labels = network.labels(pores=nearest_neighbor,
                                            mode='and')
                    common_label = [l for l in labels if 'surface_' in l]
                    new_neighbors = network.pores(common_label)
            elif neighbor in Pn_new_net:
                new_neighbors = nearest_neighbor
            connect_pores(network=network, pores1=neighbor,
                          pores2=new_neighbors, labels=labels)

        # Removing temporary labels
        for l in main_labels:
            network['pore.surface_' + l] = False
        new_net['pore.coords'] = np.copy(old_coords)

    label_faces(network=network)
    for l in main_labels:
        del network['pore.surface_'+l]
    trim(network=network, pores=pores)
    ws = network.project.workspace
    ws.close_project(new_net.project)


def trim_occluded_throats(network, mask='all'):
    r"""
    Remove throats with zero area from the network and also remove
    pores that are isolated (as a result or otherwise)

    Parameters
    ----------
    network : OpenPNM Network Object

    mask : string
        Applies routine only to pores and throats with this label
    """
    occluded_ts = network['throat.area'] == 0
    if np.sum(occluded_ts) > 0:
        occluded_ts *= network["throat."+mask]
        trim(network=network, throats=occluded_ts)


def merge_pores(network, pores, labels=['merged']):
    r"""
    Combines a selection of pores into a new single pore located at the
    centroid of the selected pores and connected to all of their neighbors.

    Parameters
    ----------
    network : OpenPNM Network Object

    pores : array_like
        The list of pores which are to be combined into a new single pore

    labels : string or list of strings
        The labels to apply to the new pore and new throat connections

    Notes
    -----
    (1) The method also works if a list of lists is passed, in which case
    it consecutively merges the given selections of pores.

    (2) The selection of pores should be chosen carefully, preferrable so that
    they all form a continuous cluster.  For instance, it is recommended
    to use the ``find_nearby_pores`` method to find all pores within a
    certain distance of a given pore, and these can then be merged without
    causing any abnormal connections.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[20, 20, 1])
    >>> Ps = pn.find_nearby_pores(pores=111, r=5, flatten=True)
    >>> op.topotools.merge_pores(network=pn, pores=Ps, labels=['merged'])
    >>> print(pn.Np)
    321
    >>> pn.pores('merged')
    array([320])
    >>> pn.num_throats('merged')
    32

    """
    # Assert that `pores` is list of lists
    try:
        len(pores[0])
    except (TypeError, IndexError):
        pores = [pores]

    N = len(pores)
    NBs, XYZs = [], []

    for Ps in pores:
        temp = network.find_neighbor_pores(pores=Ps, mode='union', flatten=True,
                                           include_input=False)
        NBs.append(temp)
        points = np.concatenate((temp, Ps))
        XYZs.append(hull_centroid(network["pore.coords"][points]))

    extend(network, pore_coords=XYZs, labels=labels)
    Pnew = network.Ps[-N::]

    # Possible throats between new pores: This only happens when running in
    # batch mode, i.e. multiple groups of pores are to be merged. In case
    # some of these groups share elements, possible throats between the
    # intersecting elements is not captured and must be added manually.
    pores_set = [set(items) for items in pores]
    NBs_set = [set(items) for items in NBs]
    ps1, ps2 = [], []
    from itertools import combinations
    for i, j in combinations(range(N), 2):
        if not NBs_set[i].isdisjoint(pores_set[j]):
            ps1.append([network.Ps[-N+i]])
            ps2.append([network.Ps[-N+j]])

    # Add (possible) connections between the new pores
    connect_pores(network, pores1=ps1, pores2=ps2, labels=labels)
    # Add connections between the new pores and the rest of the network
    connect_pores(network, pores2=np.split(Pnew, N), pores1=NBs, labels=labels)
    # Trim merged pores from the network
    trim(network=network, pores=np.concatenate(pores))


def hull_centroid(points):
    r"""
    Computes centroid of the convex hull enclosing the given coordinates.

    Parameters
    ----------
    points : Np by 3 ndarray
        Coordinates (xyz)

    Returns
    -------
    A 3 by 1 Numpy array containing coordinates of the centroid.

    """
    dim = [np.unique(points[:, i]).size != 1 for i in range(3)]
    hull = ConvexHull(points[:, dim])
    centroid = points.mean(axis=0)
    centroid[dim] = hull.points[hull.vertices].mean(axis=0)

    return centroid


def _template_sphere_disc(dim, outer_radius, inner_radius):
    r"""
    This private method generates an image array of a sphere/shell-disc/ring.
    It is useful for passing to Cubic networks as a ``template`` to make
    networks with desired shapes.

    Parameters
    ----------
    dim : int
        Network dimension

    outer_radius : int
        Number of the nodes in the outer radius of the network

    inner_radius : int
        Number of the nodes in the inner radius of the network

    Returns
    -------
    A Numpy array containing 1's to demarcate the desired shape, and 0's
    elsewhere.

    """
    rmax = np.array(outer_radius, ndmin=1)
    rmin = np.array(inner_radius, ndmin=1)
    ind = 2 * rmax - 1
    coord = sp.indices((ind * np.ones(dim, dtype=int)))
    coord = coord - (ind - 1)/2
    x = coord[0, :]
    y = coord[1, :]
    if dim == 2:
        img = (x ** 2 + y ** 2) < rmax ** 2
    elif dim == 3:
        z = coord[2, :]
        img = (x ** 2 + y ** 2 + z ** 2) < rmax ** 2
    if rmin[0] != 0:
        if dim == 2:
            img_min = (x ** 2 + y ** 2) > rmin ** 2
        elif dim == 3:
            img_min = (x ** 2 + y ** 2 + z ** 2) > rmin ** 2
        img = img * img_min
    return img


def template_sphere_shell(outer_radius, inner_radius=0, dim=3):
    r"""
    This method generates an image array of a sphere-shell. It is useful for
    passing to Cubic networks as a ``template`` to make spherical shaped
    networks.

    Parameters
    ----------
    outer_radius : int
        Number of nodes in the outer radius of the sphere.

    inner_radius : int
        Number of nodes in the inner radius of the shell.  a value of 0 will
        result in a solid sphere.

    dim : scalar
        Controls the number of dimensions of the result.  3 returns a sphere,
        while 2 returns a disk.

    Returns
    -------
    A Numpy array containing 1's to demarcate the sphere-shell, and 0's
    elsewhere.

    """
    img = _template_sphere_disc(dim=dim, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    return img


def template_cylinder_annulus(height, outer_radius, inner_radius=0):
    r"""
    This method generates an image array of a disc-ring.  It is useful for
    passing to Cubic networks as a ``template`` to make circular-shaped 2D
    networks.

    Parameters
    ----------
    height : int
        The height of the cylinder

    outer_radius : int
        Number of nodes in the outer radius of the cylinder

    inner_radius : int
        Number of the nodes in the inner radius of the annulus.  A value of 0
        will result in a solid cylinder.

    Returns
    -------
    A Numpy array containing 1's to demarcate the disc-ring, and 0's
    elsewhere.

    """

    img = _template_sphere_disc(dim=2, outer_radius=outer_radius,
                                inner_radius=inner_radius)
    img = np.tile(np.atleast_3d(img), reps=height)
    return img


def generate_base_points(num_points, domain_size, density_map=None,
                         reflect=True):
    r"""
    Generates a set of base points for passing into the Tessellation-based
    Network classes.  The points can be distributed in spherical, cylindrical,
    or rectilinear patterns, as well as 2D and 3D (disks and squares).

    Parameters
    ----------
    num_points : scalar
        The number of base points that lie within the domain.  Note that the
        actual number of points returned will be larger, with the extra points
        lying outside the domain.

    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.

    density_map : array, optional
        A an array that contains fractional values (0 < i < 1) indicating the
        liklihood that a point in that region should be kept.  The size of this
        array can be anything, but the shape must match the ``domain_size``;
        that is for a 3D network the shape of the ``density_map`` can be
        [10, 10, 10] or [50, 50, 50], depending on how important the resolution
        of the density distribution is.  For a 2D network the ``density_map``
        should be [10, 10].

        When specifying a custom probabiliy map is it recommended to also set
        values outside the given domain to zero.  If not, then the correct
        shape will still be returned, but with too few points in it.

    reflect : boolean
        If True, the the base points are generated as specified, the reflected
        about each face of the domain.  This essentially tricks the
        tessellation functions into creating smooth flat faces at the
        boundaries once these excess pores are trimmed.

    Notes
    -----
    The reflection approach tends to create larger pores near the surfaces, so
    it might be necessary to use the ``density_map`` argument to specify a
    slightly higher density of points near the surfaces.

    The ``Voronoi``, ``Delaunay``, ``Gabriel``, and ``DelunayVoronoiDual``
    classes can *techncially* handle base points with spherical or cylindrical
    domains, but the reflection across round surfaces does not create perfect
    Voronoi cells so the surfaces will not be smooth.


    Examples
    --------
    The following generates a spherical array with higher values near the core.
    It uses a distance transform to create a sphere of radius 10, then a
    second distance transform to create larger values in the center away from
    the sphere surface.  These distance values could be further skewed by
    applying a power, with values higher than 1 resulting in higher values in
    the core, and fractional values smoothinging them out a bit.

    >>> import openpnm as op
    >>> import scipy as sp
    >>> import scipy.ndimage as spim
    >>> im = np.ones([21, 21, 21], dtype=int)
    >>> im[10, 10, 10] = 0
    >>> im = spim.distance_transform_edt(im) <= 20  # Create sphere of 1's
    >>> prob = spim.distance_transform_edt(im)
    >>> prob = prob / np.amax(prob)  # Normalize between 0 and 1
    >>> pts = op.topotools.generate_base_points(num_points=50,
    ...                                         domain_size=[1, 1, 1],
    ...                                         density_map=prob)
    >>> net = op.network.DelaunayVoronoiDual(points=pts, shape=[1, 1, 1])

    """
    def _try_points(num_points, prob):
        prob = np.atleast_3d(prob)
        prob = np.array(prob)/np.amax(prob)  # Ensure prob is normalized
        base_pts = []
        N = 0
        while N < num_points:
            pt = sp.random.rand(3)  # Generate a point
            # Test whether to keep it or not
            [indx, indy, indz] = np.floor(pt*np.shape(prob)).astype(int)
            if sp.random.rand(1) <= prob[indx][indy][indz]:
                base_pts.append(pt)
                N += 1
        base_pts = np.array(base_pts)
        return base_pts

    if len(domain_size) == 1:  # Spherical
        domain_size = np.array(domain_size)
        r = domain_size[0]
        if density_map is None:
            # Make an image of a sphere filled with ones and use _try_points
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, 20] = 0
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to spherical coordinates
        [X, Y, Z] = np.array(base_pts - [0.5, 0.5, 0.5]).T
        r = 2*np.sqrt(X**2 + Y**2 + Z**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        phi = 2*np.arctan(np.sqrt(X**2 + Y**2)/Z)
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, phi] = [r[inds], theta[inds], phi[inds]]
        # Reflect base points across perimeter
        if reflect:
            r, theta, phi = reflect_base_points(np.vstack((r, theta, phi)),
                                                domain_size)
        # Convert to Cartesean coordinates
        X = r*np.cos(theta)*np.sin(phi)
        Y = r*np.sin(theta)*np.sin(phi)
        Z = r*np.cos(phi)
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 2:  # Cylindrical or Disk
        domain_size = np.array(domain_size)
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            density_map[20, 20, :] = 0
            if domain_size[1] == 0:  # Disk
                density_map = density_map[:, :, 0]
            density_map = spim.distance_transform_edt(density_map) < 20
        base_pts = _try_points(num_points, density_map)
        # Convert to cylindrical coordinates
        [X, Y, Z] = np.array(base_pts - [0.5, 0.5, 0]).T  # Center on z-axis
        r = 2*np.sqrt(X**2 + Y**2)*domain_size[0]
        theta = 2*np.arctan(Y/X)
        z = Z*domain_size[1]
        # Trim points outside the domain (from improper prob images)
        inds = r <= domain_size[0]
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        inds = ~((z > domain_size[1]) + (z < 0))
        [r, theta, z] = [r[inds], theta[inds], z[inds]]
        if reflect:
            r, theta, z = reflect_base_points(np.vstack([r, theta, z]),
                                              domain_size)
        # Convert to Cartesean coordinates
        X = r*np.cos(theta)
        Y = r*np.sin(theta)
        Z = z
        base_pts = np.vstack([X, Y, Z]).T

    elif len(domain_size) == 3:  # Cube or square
        if density_map is None:
            density_map = np.ones([41, 41, 41])
            if domain_size[2] == 0:
                density_map = density_map[:, :, 0]
        base_pts = _try_points(num_points, density_map)
        base_pts = base_pts*domain_size
        if reflect:
            base_pts = reflect_base_points(base_pts, domain_size)

    return base_pts


def reflect_base_points(base_pts, domain_size):
    r'''
    Helper function for relecting a set of points about the faces of a
    given domain.

    Parameters
    ----------
    base_pts : array_like
        The coordinates of the base_pts to be reflected in the coordinate
        system corresponding to the the domain as follows:

        **spherical** : [r, theta, phi]
        **cylindrical** or **circular** : [r, theta, z]
        **rectangular** or **square** : [x, y, z]

    domain_size : list or array
        Controls the size and shape of the domain, as follows:

        **sphere** : If a single value is received, its treated as the radius
        [r] of a sphere centered on [0, 0, 0].

        **cylinder** : If a two-element list is received it's treated as the
        radius and height of a cylinder [r, z] positioned at [0, 0, 0] and
        extending in the positive z-direction.  If the z dimension is 0, a
        disk of radius r is created.

        **rectangle** : If a three element list is received, it's treated
        as the outer corner of rectangle [x, y, z] whose opposite corner lies
        at [0, 0, 0].  If the z dimension is 0, a rectangle of size X-by-Y is
        created.

    '''
    domain_size = np.array(domain_size)
    if len(domain_size) == 1:
        r, theta, phi = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        phi = np.hstack([phi, phi])
        base_pts = np.vstack((r, theta, phi))
    if len(domain_size) == 2:
        r, theta, z = base_pts
        new_r = 2*domain_size[0] - r
        r = np.hstack([r, new_r])
        theta = np.hstack([theta, theta])
        z = np.hstack([z, z])
        if domain_size[1] != 0:  # If not a disk
            r = np.hstack([r, r, r])
            theta = np.hstack([theta, theta, theta])
            z = np.hstack([z, -z, 2-z])
        base_pts = np.vstack((r, theta, z))
    elif len(domain_size) == 3:
        Nx, Ny, Nz = domain_size
        # Reflect base points about all 6 faces
        orig_pts = base_pts
        base_pts = np.vstack((base_pts,
                              [-1, 1, 1] * orig_pts + [2.0 * Nx, 0, 0]))
        base_pts = np.vstack((base_pts, [-1, 1, 1] * orig_pts))
        base_pts = np.vstack((base_pts,
                              [1, -1, 1] * orig_pts + [0, 2.0 * Ny, 0]))
        base_pts = np.vstack((base_pts, [1, -1, 1] * orig_pts))
        if domain_size[2] != 0:
            base_pts = np.vstack((base_pts,
                                  [1, 1, -1] * orig_pts + [0, 0, 2.0 * Nz]))
            base_pts = np.vstack((base_pts, [1, 1, -1] * orig_pts))
    return base_pts


def find_clusters(network, mask=[], t_labels=False):
    r"""
    Identify connected clusters of pores in the network.  This method can
    also return a list of throat cluster numbers, which correspond to the
    cluster numbers of the pores to which the throat is connected.  Either
    site and bond percolation can be considered, see description of input
    arguments for details.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network

    mask : array_like, boolean
        A list of active bonds or sites (throats or pores).  If the mask is
        Np long, then the method will perform a site percolation, and if
        the mask is Nt long bond percolation will be performed.

    Returns
    -------
    A tuple containing an Np long list of pore cluster labels, and an Nt-long
    list of throat cluster labels.  The label numbers correspond such that
    pores and throats with the same label are part of the same cluster.

    Examples
    --------
    >>> import openpnm as op
    >>> from scipy import rand
    >>> pn = op.network.Cubic(shape=[25, 25, 1])
    >>> pn['pore.seed'] = rand(pn.Np)
    >>> pn['throat.seed'] = rand(pn.Nt)


    """
    # Parse the input arguments
    mask = np.array(mask, ndmin=1)
    if mask.dtype != bool:
        raise Exception('Mask must be a boolean array of Np or Nt length')

    # If pore mask was given perform site percolation
    if np.size(mask) == network.Np:
        (p_clusters, t_clusters) = _site_percolation(network, mask)
    # If pore mask was given perform bond percolation
    elif np.size(mask) == network.Nt:
        (p_clusters, t_clusters) = _bond_percolation(network, mask)
    else:
        raise Exception('Mask received was neither Nt nor Np long')

    return (p_clusters, t_clusters)


def _site_percolation(network, pmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Find throats that produce site percolation
    conns = np.copy(network['throat.conns'])
    conns[:, 0] = pmask[conns[:, 0]]
    conns[:, 1] = pmask[conns[:, 1]]
    # Only if both pores are True is the throat set to True
    tmask = np.all(conns, axis=1)

    # Perform the clustering using scipy.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Adjust cluster numbers such that non-invaded pores are labelled -1
    # Note: The following line also takes care of assigning cluster numbers
    # to single isolated invaded pores
    p_clusters = (clusters + 1)*(pmask) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']]
    ind = (t_clusters[:, 0] == t_clusters[:, 1])
    t_clusters = t_clusters[:, 0]
    # Label non-invaded throats with -1
    t_clusters[~ind] = -1

    return (p_clusters, t_clusters)


def _bond_percolation(network, tmask):
    r"""
    This private method is called by 'find_clusters'
    """
    # Perform the clustering using scipy.csgraph
    csr = network.create_adjacency_matrix(weights=tmask, fmt='csr',
                                          drop_zeros=True)
    clusters = sprs.csgraph.connected_components(csgraph=csr,
                                                 directed=False)[1]

    # Convert clusters to a more usable output:
    # Find pores attached to each invaded throats
    Ps = network.find_connected_pores(throats=tmask, flatten=True)
    # Adjust cluster numbers such that non-invaded pores are labelled -1
    p_clusters = (clusters + 1)*(network.tomask(pores=Ps).astype(int)) - 1
    # Label invaded throats with their neighboring pore's label
    t_clusters = clusters[network['throat.conns']][:, 0]
    # Label non-invaded throats with -1
    t_clusters[~tmask] = -1

    return (p_clusters, t_clusters)


def add_boundary_pores(network, pores, offset=None, move_to=None,
                       apply_label='boundary'):
    r"""
    This method uses ``clone_pores`` to clone the input pores, then shifts
    them the specified amount and direction, then applies the given label.

    Parameters
    ----------
    pores : array_like
        List of pores to offset.  If no pores are specified, then it
        assumes that all surface pores are to be cloned.

    offset : 3 x 1 array
        The distance in vector form which the cloned boundary pores should
        be offset.  Either this, or ``move_to`` must be specified.

    move_to : 3 x 1 array
        The location to move the boundary pores to.  A value of ``None``
        indicates that no translation should be applied in that axis.  For
        instance, ``[None, None, 0]`` indicates that the boundary pores should
        moved along the z-axis to the specified location.  Either this or
        ``offset`` must be specified.

    apply_label : string
        This label is applied to the boundary pores.  Default is
        'boundary'.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[5, 5, 5])
    >>> print(pn.Np)  # Confirm initial Network size
    125
    >>> Ps = pn.pores('top')  # Select pores on top face
    >>> pn.add_boundary_pores(labels=['top'])
    >>> print(pn.Np)  # Confirm addition of 25 new pores
    150

    """
    # Parse the input pores
    Ps = np.array(pores, ndmin=1)
    if Ps.dtype is bool:
        Ps = network.toindices(Ps)
    if np.size(pores) == 0:  # Handle an empty array if given
        return np.array([], dtype=sp.int64)
    # Clone the specifed pores
    clone_pores(network=network, pores=Ps)
    newPs = network.pores('pore.clone')
    del network['pore.clone']
    newTs = network.throats('clone')
    del network['throat.clone']
    if offset is not None:  # Offset the cloned pores
        network['pore.coords'][newPs] += offset
    if move_to is not None:  # Move the cloned pores
        for i, d in enumerate(move_to):
            if d is not None:
                temp = network['pore.coords'][newPs]
                temp[:, i] = d
                network['pore.coords'][newPs] = temp
    # Apply labels to boundary pores (trim leading 'pores' if present)
    label = apply_label.split('.')[-1]
    plabel = 'pore.' + label
    tlabel = 'throat.' + label
    network[plabel] = False
    network[plabel][newPs] = True
    network[tlabel] = False
    network[tlabel][newTs] = True


def find_path(network, pore_pairs, weights=None):
    r"""
    Find the shortest path between pairs of pores.

    Parameters
    ----------
    network : OpenPNM Network Object
        The Network object on which the search should be performed

    pore_pairs : array_like
        An N x 2 array containing N pairs of pores for which the shortest
        path is sought.

    weights : array_like, optional
        An Nt-long list of throat weights for the search.  Typically this
        would be the throat lengths, but could also be used to represent
        the phase configuration.  If no weights are given then the
        standard topological connections of the Network are used.

    Returns
    -------
    A dictionary containing both the pores and throats that define the
    shortest path connecting each pair of input pores.

    Notes
    -----
    The shortest path is found using Dijkstra's algorithm included in the
    scipy.sparse.csgraph module

    TODO: The returned throat path contains the correct values, but not
    necessarily in the true order

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[3, 3, 3])
    >>> a = op.topotools.find_path(network=pn, pore_pairs=[[0, 4], [0, 10]])
    >>> a['pores']
    [array([0, 1, 4]), array([ 0,  1, 10])]
    >>> a['throats']
    [array([ 0, 19]), array([ 0, 37])]
    """
    Ps = np.array(pore_pairs, ndmin=2)
    if weights is None:
        weights = np.ones_like(network.Ts)
    graph = network.create_adjacency_matrix(weights=weights, fmt='csr',
                                            drop_zeros=False)
    paths = csgraph.dijkstra(csgraph=graph, indices=Ps[:, 0],
                             return_predecessors=True)[1]
    pores = []
    throats = []
    for row in range(0, np.shape(Ps)[0]):
        j = Ps[row][1]
        ans = []
        while paths[row][j] > -9999:
            ans.append(j)
            j = paths[row][j]
        ans.append(Ps[row][0])
        ans.reverse()
        pores.append(np.array(ans, dtype=int))
        Ts = network.find_neighbor_throats(pores=ans, mode='xnor')
        throats.append(np.array(Ts, dtype=int))
    pdict = PrintableDict
    dict_ = pdict(**{'pores': pores, 'throats': throats})
    return dict_


def iscoplanar(coords):
    r'''
    Determines if given pores are coplanar with each other

    Parameters
    ----------
    coords : array_like
        List of pore coords to check for coplanarity.  At least 3 pores are
        required.

    Returns
    -------
    A boolean value of whether given points are coplanar (True) or not (False)
    '''
    coords = np.array(coords, ndmin=1)
    if np.shape(coords)[0] < 3:
        raise Exception('At least 3 input pores are required')

    Px = coords[:, 0]
    Py = coords[:, 1]
    Pz = coords[:, 2]

    # Do easy check first, for common coordinate
    if np.shape(np.unique(Px))[0] == 1:
        return True
    if np.shape(np.unique(Py))[0] == 1:
        return True
    if np.shape(np.unique(Pz))[0] == 1:
        return True

    # Perform rigorous check using vector algebra
    # Grab first basis vector from list of coords
    n1 = np.array((Px[1] - Px[0], Py[1] - Py[0], Pz[1] - Pz[0])).T
    n = np.array([0.0, 0.0, 0.0])
    i = 1
    while n.sum() == 0:
        if i >= (np.size(Px) - 1):
            logger.warning('No valid basis vectors found')
            return False
        # Chose a secon basis vector
        n2 = np.array((Px[i+1] - Px[i], Py[i+1] - Py[i], Pz[i+1] - Pz[i])).T
        # Find their cross product
        n = np.cross(n1, n2)
        i += 1
    # Create vectors between all other pairs of points
    r = np.array((Px[1:-1] - Px[0], Py[1:-1] - Py[0], Pz[1:-1] - Pz[0]))
    # Ensure they all lie on the same plane
    n_dot = np.dot(n, r)

    if np.sum(np.absolute(n_dot)) == 0:
        return True
    else:
        return False
